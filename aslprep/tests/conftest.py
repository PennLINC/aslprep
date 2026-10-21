"""Fixtures for the CircleCI tests."""

import base64
import os

import pytest


def pytest_addoption(parser):
    """Collect pytest parameters for running tests."""
    import tempfile

    parser.addoption('--working_dir', action='store', default=tempfile.mkdtemp())
    parser.addoption('--data_dir', action='store')
    parser.addoption('--output_dir', action='store')


# Set up the commandline options as fixtures
@pytest.fixture(scope='session')
def data_dir(request):
    """Grab data directory."""
    return request.config.getoption('--data_dir')


@pytest.fixture(scope='session')
def working_dir(request):
    """Grab working directory."""
    workdir = request.config.getoption('--working_dir')
    os.makedirs(workdir, exist_ok=True)
    return workdir


@pytest.fixture(scope='session')
def output_dir(request):
    """Grab output directory."""
    outdir = request.config.getoption('--output_dir')
    os.makedirs(outdir, exist_ok=True)
    return outdir


@pytest.fixture(scope='session')
def datasets(data_dir):
    """Locate downloaded datasets."""
    return {
        'qtab': os.path.join(data_dir, 'qtab'),
        'nonqtab': os.path.join(data_dir, 'nonqtab'),
    }


@pytest.fixture(scope='session', autouse=True)
def _fslicense(working_dir):
    """Set the FreeSurfer license as an environment variable."""
    FS_LICENSE = os.path.join(working_dir, 'license.txt')
    os.environ['FS_LICENSE'] = FS_LICENSE
    LICENSE_CODE = (
        'bWF0dGhldy5jaWVzbGFrQHBzeWNoLnVjc2IuZWR1CjIwNzA2CipDZmVWZEg1VVQ4clkKRlNCWVouVWtlVElDdwo='
    )
    with open(FS_LICENSE, 'w') as f:
        f.write(base64.b64decode(LICENSE_CODE).decode())


# TEMPORARY (PR 767 CI diagnosis; to be reverted): after each test, can `python` on PATH
# still import surfa, as mri_synthstrip needs?
_SURFA_PROBE = []


def _surfa_probe():
    import shutil
    import subprocess

    proc = subprocess.run(
        ['python', '-c', 'import sys, surfa; print(sys.executable, surfa.__file__)'],
        capture_output=True,
        text=True,
        check=False,
    )
    tail = (proc.stdout or proc.stderr).strip().splitlines()[-1:]
    return proc.returncode, shutil.which('python'), tail


@pytest.fixture(autouse=True)
def _surfa_probe_fixture(request):
    yield
    code, which, tail = _surfa_probe()
    if code or not _SURFA_PROBE:
        _SURFA_PROBE.append((request.node.nodeid, code, which, tail, os.environ.get('PATH')))


def pytest_terminal_summary(terminalreporter):
    for nodeid, code, which, tail, path in _SURFA_PROBE[:4]:
        terminalreporter.write_line(
            f'SURFA-PROBE after {nodeid}: exit {code}, python {which}, {tail}, PATH {path}'
        )
