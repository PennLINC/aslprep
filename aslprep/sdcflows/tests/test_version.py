"""Test _version.py."""
import sys
from collections import namedtuple
from pkg_resources import DistributionNotFound
from importlib import reload
import sdcflows


def test_version_scm0(monkeypatch):
    """Retrieve the version via setuptools_scm."""

    class _version:
        __version__ = "10.0.0"

    monkeypatch.setitem(sys.modules, "sdcflows._version", _version)
    reload(sdcflows)
    assert sdcflows.__version__ == "10.0.0"


def test_version_scm1(monkeypatch):
    """Retrieve the version via pkg_resources."""
    monkeypatch.setitem(sys.modules, "sdcflows._version", None)

    def _dist(name):
        Distribution = namedtuple("Distribution", ["name", "version"])
        return Distribution(name, "success")

    monkeypatch.setattr("pkg_resources.get_distribution", _dist)
    reload(sdcflows)
    assert sdcflows.__version__ == "success"


def test_version_scm2(monkeypatch):
    """Check version could not be interpolated."""
    monkeypatch.setitem(sys.modules, "sdcflows._version", None)

    def _raise(name):
        raise DistributionNotFound("No get_distribution mock")

    monkeypatch.setattr("pkg_resources.get_distribution", _raise)
    reload(sdcflows)
    assert sdcflows.__version__ == "unknown"
