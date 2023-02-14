# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

# Load modules for compatibility
from aslprep.niworkflows.interfaces import bids, cifti, freesurfer, images, itk, surf, utils
from aslprep.niworkflows.interfaces.plotting import CBFSummary, CBFtsSummary
from aslprep.interfaces.confounds import ASLSummary, GatherConfounds
from aslprep.interfaces.multiecho import T2SMap
from aslprep.interfaces.reports import AboutSummary, FunctionalSummary, SubjectSummary


class DerivativesDataSink(bids.DerivativesDataSink):
    out_path_base = "aslprep"


__all__ = [
    "bids",
    "cifti",
    "freesurfer",
    "images",
    "itk",
    "surf",
    "utils",
    "SubjectSummary",
    "FunctionalSummary",
    "AboutSummary",
    "GatherConfounds",
    "ASLSummary",
    "CBFSummary",
    "CBFtsSummary",
    "T2SMap",
    "DerivativesDataSink",
]
