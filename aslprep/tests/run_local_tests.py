#!/usr/bin/env python3
"""Run tests locally in the ASLPrep test image, the way CircleCI does.

Build the image from this checkout first::

    docker build --target test -t pennlinc/aslprep:test .

The checkout is mounted read-only at CI's path and used as the working directory, so code
edits need no rebuild (only Dockerfile or pixi.lock changes do). Test data, including the
aslscan fixtures, are read from ``--data-dir`` (default ``aslprep/tests/test_data``);
generate the fixtures on the host first (see ``python -m aslprep.tests.aslscan_fixtures -h``).
Outputs go to ``aslprep/tests/pytests/``.
"""

import argparse
import os
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CI_SOURCE = '/tmp/src/aslprep'  # noqa: S108 - the path inside the container, as in CI
FORWARDED_ENV = ('ASLPREP_REQUIRE_FIXTURES',)


def _get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__.split('\n\n')[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-k', dest='test_regex', metavar='PATTERN', help='Test pattern.')
    parser.add_argument('-m', dest='test_mark', metavar='LABEL', help='Test mark label.')
    parser.add_argument(
        '--data-dir',
        type=Path,
        default=REPO / 'aslprep' / 'tests' / 'test_data',
        help='Test data directory, mounted read-only at /data.',
    )
    parser.add_argument('--image', default='pennlinc/aslprep:test', help='Docker image.')
    parser.add_argument('--cpus', type=int, default=4, help='CIRCLE_CPUS for the run.')
    return parser


def run_tests(
    test_regex=None, test_mark=None, data_dir=None, image='pennlinc/aslprep:test', cpus=4
):
    """Run pytest in the test image with CI's mounts and options."""
    data_dir = Path(data_dir or REPO / 'aslprep' / 'tests' / 'test_data').resolve()
    pytests = REPO / 'aslprep' / 'tests' / 'pytests'
    mounts = {
        REPO: f'{CI_SOURCE}:ro',
        data_dir: '/data:ro',
        pytests / 'out': '/out',
        pytests / 'work': '/work',
        pytests / 'test-results': '/test-results',
    }
    docker = ['docker', 'run', '--rm', '-e', f'CIRCLE_CPUS={cpus}']
    for name in FORWARDED_ENV:
        if name in os.environ:
            docker += ['-e', f'{name}={os.environ[name]}']
    for host, container in mounts.items():
        Path(host).mkdir(parents=True, exist_ok=True)
        docker += ['-v', f'{host}:{container}']
    docker += ['-w', CI_SOURCE]

    # Show which aslprep the container imports (it must be the mounted checkout).
    subprocess.run(
        [
            *docker,
            '--entrypoint',
            'python',
            image,
            '-c',
            'import aslprep; print(aslprep.__file__)',
        ],
        check=True,
    )
    cmd = [
        *docker,
        '--entrypoint',
        'pytest',
        image,
        '--strict-markers',
        '--strict-config',
        '-rP',
        '-o',
        'log_cli=true',
        '--junitxml=/test-results/local.xml',
        '--data_dir=/data',
        '--output_dir=/out',
        '--working_dir=/work',
    ]
    if test_regex:
        cmd += ['-k', test_regex]
    if test_mark:
        cmd += ['-m', test_mark]
    cmd.append('aslprep/tests')
    subprocess.run(cmd, check=True)


def _main(argv=None):
    options = _get_parser().parse_args(argv)
    run_tests(**vars(options))


if __name__ == '__main__':
    _main()
