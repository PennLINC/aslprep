#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
This pipeline is developed by the Poldrack lab at Stanford University
(https://poldracklab.stanford.edu/) for use at
the Center for Reproducible Neuroscience (http://reproducibility.stanford.edu/),
as well as for open-source software distribution.
"""

from .__about__ import (  # noqa
    __version__,
    __copyright__,
    __credits__,
    __packagename__,
)

import warnings

# cmp is not used by aslprep, so ignore nipype-generated warnings
warnings.filterwarnings('ignore', r'cmp not installed')
warnings.filterwarnings('ignore', r'This has not been fully tested. Please report any failures.')
warnings.filterwarnings('ignore', r"can't resolve package from __spec__ or __package__")
warnings.simplefilter('ignore', DeprecationWarning)
warnings.simplefilter('ignore', ResourceWarning)
