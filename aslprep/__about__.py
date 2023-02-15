# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Base module variables."""
from aslprep._version import get_versions

__version__ = get_versions()["version"]
del get_versions

__packagename__ = "aslprep"
__copyright__ = "Copyright 2020, PENN LINC"
__credits__ = (
    "Contributors: please check the ``.zenodo.json`` file at the top-level folder"
    "of the repository"
)
__url__ = "https://github.com/pennlinc/aslprep"

DOWNLOAD_URL = f"https://github.com/pennlinc/{__packagename__}/archive/{__version__}.tar.gz"
