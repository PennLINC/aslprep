# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Interfaces for calculating and collecting confounds."""
import os

import numpy as np
from nipype import logging
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.utils.misc import normalize_mc_params

from aslprep.utils.confounds import _gather_confounds

LOGGER = logging.getLogger("nipype.interface")


class _GatherConfoundsInputSpec(BaseInterfaceInputSpec):
    signals = File(exists=True, desc="input signals")
    dvars = File(exists=True, desc="file containing DVARS")
    rmsd = File(exists=True, desc="input RMS framewise displacement")
    std_dvars = File(exists=True, desc="file containing standardized DVARS")
    fd = File(exists=True, desc="input framewise displacement")
    motion = File(exists=True, desc="input motion parameters")


class _GatherConfoundsOutputSpec(TraitedSpec):
    confounds_file = File(exists=True, desc="output confounds file")
    confounds_list = traits.List(traits.Str, desc="list of headers")


class GatherConfounds(SimpleInterface):
    """Combine various sources of confounds in one TSV file."""

    input_spec = _GatherConfoundsInputSpec
    output_spec = _GatherConfoundsOutputSpec

    def _run_interface(self, runtime):
        combined_out, confounds_list = _gather_confounds(
            signals=self.inputs.signals,
            dvars=self.inputs.dvars,
            rmsd=self.inputs.rmsd,
            std_dvars=self.inputs.std_dvars,
            fdisp=self.inputs.fd,
            motion=self.inputs.motion,
            newpath=runtime.cwd,
        )
        self._results["confounds_file"] = combined_out
        self._results["confounds_list"] = confounds_list
        return runtime


class _GatherCBFConfoundsInputSpec(BaseInterfaceInputSpec):
    signals = File(exists=True, desc="input signals")
    score = File(exists=True, desc="SCORE outlier index")


class _GatherCBFConfoundsOutputSpec(TraitedSpec):
    confounds_file = File(exists=True, desc="output confounds file")
    confounds_list = traits.List(traits.Str, desc="list of headers")


class GatherCBFConfounds(SimpleInterface):
    """Combine various sources of confounds in one TSV file."""

    input_spec = _GatherCBFConfoundsInputSpec
    output_spec = _GatherCBFConfoundsOutputSpec

    def _run_interface(self, runtime):
        combined_out, confounds_list = _gather_confounds(
            signals=self.inputs.signals,
            dvars=None,
            rmsd=None,
            std_dvars=None,
            fdisp=None,
            motion=None,
            score=self.inputs.score,
            newpath=runtime.cwd,
        )
        self._results["confounds_file"] = combined_out
        self._results["confounds_list"] = confounds_list
        return runtime


class _NormalizeMotionParamsInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="the input parameters file")
    format = traits.Enum("FSL", "AFNI", "FSFAST", "NIPY", usedefault=True, desc="output format")


class _NormalizeMotionParamsOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="written file path")


class NormalizeMotionParams(SimpleInterface):
    """Convert input motion parameters into the designated convention.

    NOTE: Modified by Taylor Salo to support 1D arrays.
    """

    input_spec = _NormalizeMotionParamsInputSpec
    output_spec = _NormalizeMotionParamsOutputSpec

    def _run_interface(self, runtime):
        mpars = np.loadtxt(self.inputs.in_file)  # mpars is N_t x 6

        # Support single-volume motion parameters
        if mpars.ndim == 1:
            mpars = mpars[None, :]

        mpars = np.apply_along_axis(
            func1d=normalize_mc_params, axis=1, arr=mpars, source=self.inputs.format
        )
        self._results["out_file"] = os.path.join(runtime.cwd, "motion_params.txt")
        np.savetxt(self._results["out_file"], mpars)
        return runtime
