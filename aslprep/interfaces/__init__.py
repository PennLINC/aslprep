# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

# Load modules for compatibility
from ..niworkflows.interfaces import (
    bids, cifti, freesurfer, images, itk, surf, utils)

from .reports import SubjectSummary, FunctionalSummary, AboutSummary
from .confounds import GatherConfounds, ASLSummary
from .multiecho import T2SMap
from ..niworkflows.interfaces.plotting import (CBFSummary, CBFtsSummary)


class DerivativesDataSink(bids.DerivativesDataSink):
    out_path_base = 'aslprep'


__all__ = [
    'bids',
    'cifti',
    'freesurfer',
    'images',
    'itk',
    'surf',
    'utils',
    'SubjectSummary',
    'FunctionalSummary',
    'AboutSummary',
    'GatherConfounds',
    'ASLSummary',
    'CBFSummary',
    'CBFtsSummary',
    'T2SMap',
    'DerivativesDataSink',
]
