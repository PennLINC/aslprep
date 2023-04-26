# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Interfaces for calculating and collecting confounds."""
import numpy as np
import pandas as pd
from nipype import logging
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.utils.filemanip import fname_presuffix

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


class _ZigZagCorrectionInputSpec(BaseInterfaceInputSpec):
    movpar_file = File(exists=True, mandatory=True, desc="HMC motion parameters TSV.")
    aslcontext = File(exists=True, mandatory=True, desc="ASL context TSV.")


class _ZigZagCorrectionOutputSpec(TraitedSpec):
    movpar_file = File(exists=True, desc="Corrected HMC motion parameters TSV.")
    rmsd_file = File(exists=True, desc="RMSD based on corrected HMC motion parameters.")


class ZigZagCorrection(SimpleInterface):
    """Correct motion parameters for control-label intensity differences using zig-zag method."""

    input_spec = _ZigZagCorrectionInputSpec
    output_spec = _ZigZagCorrectionOutputSpec

    def _run_interface(self, runtime):
        movpar_file = self.inputs.movpar_file
        movpars_df = pd.read_table(movpar_file)
        aslcontext = pd.read_table(self.inputs.aslcontext)
        if aslcontext.shape[0] != movpars_df.shape[0]:
            raise ValueError(
                f"Number of volumes in aslcontext ({aslcontext.shape[0]}) != "
                f"number of volumes in motion parameters ({movpars_df.shape[0]})"
            )

        if not all(aslcontext["volume_type"].isin(["control", "label"])):
            raise ValueError("aslcontext contains volumes other than 'control' and 'label'.")

        voltypes = aslcontext["volume_type"].tolist()
        regressor = np.array([-1 if i == "label" else 1 for i in voltypes])

        # Mean-center the motion parameters before performing the regression.
        movpars = movpars_df.to_numpy()
        movpars_mean = np.mean(movpars, axis=0)
        movpars_mc = movpars - movpars_mean

        # Regress the mean-centered motion parameters on the zig-zag regressor.
        betas = np.linalg.lstsq(regressor[:, None], movpars_mc)[0]
        predicted_movpars = np.dot(regressor[:, None], betas)

        # Remove the predicted motion parameters from the actual motion parameters.
        corrected_movpars = movpars - predicted_movpars

        # Save out the results
        self._results["movpar_file"] = fname_presuffix(
            self.inputs.movpar_file,
            suffix="_zigzag",
            newpath=runtime.cwd,
        )
        corrected_movpars_df = pd.DataFrame(columns=movpars_df.columns, data=corrected_movpars)
        corrected_movpars_df.to_csv(self._results["movpar_file"], sep="\t", index=False)

        self._results["rmsd_file"] = fname_presuffix(
            self.inputs.movpar_file,
            suffix="_rmsd",
            newpath=runtime.cwd,
        )
        rmsd_df = pd.DataFrame(columns=movpars_df.columns, data=corrected_movpars)
        rmsd_df.to_csv(self._results["rmsd_file"], sep="\t", index=False)

        return runtime
