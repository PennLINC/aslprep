# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Nipype interfaces for aslprep."""

from aslprep.interfaces import (
    ants,
    bids,
    cbf,
    confounds,
    parcellation,
    plotting,
    reference,
    reports,
    utility,
)

__all__ = [
    'ants',
    'bids',
    'cbf',
    'confounds',
    'parcellation',
    'plotting',
    'reference',
    'reports',
    'utility',
]
