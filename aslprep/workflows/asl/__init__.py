# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for the ASL processing portions of ASLPrep."""
from aslprep.workflows.asl import base, cbf, confounds, hmc, registration, resampling

__all__ = [
    "base",
    "cbf",
    "confounds",
    "hmc",
    "registration",
    "resampling",
]
