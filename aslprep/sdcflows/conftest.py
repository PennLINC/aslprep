"""py.test configuration"""
import os
from pathlib import Path
import numpy
import pytest
from bids.layout import BIDSLayout

test_data_env = os.getenv('TEST_DATA_HOME', str(Path.home() / 'sdcflows-tests'))
test_output_dir = os.getenv('TEST_OUTPUT_DIR')
test_workdir = os.getenv('TEST_WORK_DIR')

layouts = {p.name: BIDSLayout(str(p), validate=False, derivatives=True)
           for p in Path(test_data_env).glob('*') if p.is_dir()}


def pytest_report_header(config):
    msg = "Datasets found: %s" % ', '.join([v.root for v in layouts.values()])
    if test_output_dir is not None:
        msg += '\nOutput folder: %s' % Path(test_output_dir).resolve()
    return msg


@pytest.fixture(autouse=True)
def add_np(doctest_namespace):
    doctest_namespace['np'] = numpy
    doctest_namespace['os'] = os
    doctest_namespace['Path'] = Path
    for key, val in list(layouts.items()):
        doctest_namespace[key] = Path(val.root)


@pytest.fixture
def workdir():
    return None if test_workdir is None else Path(test_workdir)


@pytest.fixture
def output_path():
    return None if test_output_dir is None else Path(test_output_dir)


@pytest.fixture
def bids_layouts():
    return layouts
