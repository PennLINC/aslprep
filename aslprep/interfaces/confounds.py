# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Interfaces for calculating and collecting confounds."""
from nipype import logging
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)

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
