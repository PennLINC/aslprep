"""Fixtures for the CircleCI tests."""
import os

import pytest


def pytest_addoption(parser):
    """Collect pytest parameters for running tests."""
    parser.addoption("--working_dir", action="store", default="/tmp")
    parser.addoption("--data_dir", action="store")
    parser.addoption("--output_dir", action="store")


# Set up the commandline options as fixtures
@pytest.fixture
def data_dir(request):
    """Grab data directory."""
    return request.config.getoption("--data_dir")


@pytest.fixture
def working_dir(request):
    """Grab working directory."""
    return request.config.getoption("--working_dir")


@pytest.fixture
def output_dir(request):
    """Grab output directory."""
    return request.config.getoption("--output_dir")


@pytest.fixture
def datasets(data_dir):
    """Locate downloaded datasets."""
    dsets = {}
    dsets["dset"] = os.path.join(data_dir, "dset")
    dsets["smriprep"] = os.path.join(dsets["dset"], "derivatives/smriprep")
    return dsets
