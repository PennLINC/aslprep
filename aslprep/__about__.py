# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Base module variables."""

try:
    from aslprep._version import __version__
except ImportError:
    __version__ = '0+unknown'

__packagename__ = 'aslprep'
__copyright__ = 'Copyright 2023, The ASLPrep Developers'
__credits__ = (
    'Contributors: please check the ``.zenodo.json`` file at the top-level folder '
    'of the repository.'
)
__url__ = 'https://github.com/PennLINC/aslprep'

DOWNLOAD_URL = f'https://github.com/PennLINC/{__packagename__}/archive/{__version__}.tar.gz'
