"""SDCflows - :abbr:`SDC (susceptibility distortion correction)` by DUMMIES, for dummies.

The copy of SDCflows used by ASLPrep has been modified by the PennLINC developers.
Please see the ASLPrep license file for SDCflows's licensing.
"""
__packagename__ = "sdcflows"
__copyright__ = "2020, The NiPreps developers"
try:
    from ._version import __version__
except ModuleNotFoundError:
    from pkg_resources import get_distribution, DistributionNotFound
    try:
        __version__ = get_distribution(__packagename__).version
    except DistributionNotFound:
        __version__ = "unknown"
    del get_distribution
    del DistributionNotFound
