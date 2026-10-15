#!/usr/bin/env python3
"""Run a harness candidate on simulated fixtures and score each against its truth.

For every fixture directory under ``<data dir>/aslscan`` (or the ``--names`` given), run
``harness.py`` in the ASLPrep image, then ``score.py`` here, and write
``<out>/<label>/<fixture>/score.json``. The resampled series are deleted after scoring, and
with ``--delete-fixtures`` so is the fixture itself (disk space is limited; fixtures can be
regenerated from their recipe and seed).

Usage::

    python scripts/hmc_probe/batch_fixtures.py --data-dir <dir> --src <code dir> \\
        --out <dir> --label c0 [--candidate c0] [--names ...] [--jobs 6] [--delete-fixtures]
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_panel import IMAGE, SITE  # noqa: E402

#: niimath (v1.0.20260924, official release binary), mounted for the ``niimath`` engines; the
#: ASLPrep image does not include it.
NIIMATH_DIR = Path('/home/tsalo/hmc-probe/tools/niimath-1.0.20260924')


def run_one(fixture, args):
    out = args.out / args.label / fixture.name
    if (out / 'score.json').exists() and not args.rerun:
        return fixture.name, 'done'
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    asl = next((fixture / 'sub-01' / 'perf').glob('*_asl.nii.gz'))
    cmd = ['docker', 'run', '--rm', '-u', f'{os.getuid()}:{os.getgid()}', '-e', 'HOME=/tmp']
    cmd += ['-e', 'FSLOUTPUTTYPE=NIFTI_GZ']
    for host, container, mode in (
        (HERE, '/probe', 'ro'),
        (args.src / 'aslprep', SITE, 'ro'),
        (fixture, '/fx', 'ro'),
        (out, '/out', 'rw'),
    ):
        cmd += ['-v', f'{host}:{container}:{mode}']
    if (NIIMATH_DIR / 'niimath').exists():
        cmd += ['-v', f'{NIIMATH_DIR}:/opt/niimath:ro', '-e', 'NIIMATH=/opt/niimath/niimath']
    cmd += ['--entrypoint', '/app/.pixi/envs/aslprep/bin/python', IMAGE, '/probe/harness.py']
    cmd += ['--asl', f'/fx/sub-01/perf/{asl.name}', '--candidate', args.candidate, '--out', '/out']
    cmd += args.harness_args
    with (out / 'harness.log').open('w') as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)
    if proc.returncode:
        return fixture.name, f'harness FAILED ({proc.returncode})'
    from score import score

    try:
        result = score(fixture, out)
    except Exception as exc:  # noqa: BLE001 - record and continue the batch
        (out / 'score_error.txt').write_text(repr(exc))
        return fixture.name, f'score FAILED: {exc!r}'[:200]
    (out / 'score.json').write_text(json.dumps(result, indent=2) + '\n')
    (out / 'corrected.nii.gz').unlink(missing_ok=True)
    if args.delete_fixtures:
        shutil.rmtree(fixture, ignore_errors=True)
    return fixture.name, 'ok'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--src', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--candidate', default='c0')
    parser.add_argument('--names', nargs='*')
    parser.add_argument('--jobs', type=int, default=6)
    parser.add_argument('--rerun', action='store_true')
    parser.add_argument('--delete-fixtures', action='store_true')
    parser.add_argument(
        '--harness-args', nargs=argparse.REMAINDER, default=[], help='Passed to harness.py.'
    )
    args = parser.parse_args(argv)
    args.out, args.src = args.out.resolve(), args.src.resolve()
    root = args.data_dir.resolve() / 'aslscan'
    fixtures = sorted(
        p
        for p in root.iterdir()
        if (p / 'manifest.json').exists() and (not args.names or p.name in args.names)
    )
    print(f'{len(fixtures)} fixtures, {args.jobs} at a time', flush=True)
    with ThreadPoolExecutor(args.jobs) as pool:
        for name, status in pool.map(lambda f: run_one(f, args), fixtures):
            print(f'{name}: {status}', flush=True)


if __name__ == '__main__':
    main()
