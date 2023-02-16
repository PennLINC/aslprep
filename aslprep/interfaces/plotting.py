"""Plotting interfaces."""
import numpy as np
import pandas as pd
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    isdefined,
    traits,
)
from nipype.utils.filemanip import fname_presuffix

from aslprep.utils.plotting import ASLPlot, CBFPlot, CBFtsPlot


class _ASLSummaryInputSpec(BaseInterfaceInputSpec):
    in_func = File(exists=True, mandatory=True, desc="")
    in_mask = File(exists=True, desc="")
    in_segm = File(exists=True, desc="")
    in_spikes_bg = File(exists=True, mandatory=True, desc="")
    fd = File(exists=True, mandatory=True, desc="")
    fd_thres = traits.Float(0.2, usedefault=True, desc="")
    dvars = File(exists=True, mandatory=True, desc="")
    outliers = File(exists=True, mandatory=True, desc="")
    tr = traits.Either(None, traits.Float, usedefault=True, desc="the TR")


class _ASLSummaryOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="written file path")


class ASLSummary(SimpleInterface):
    """Prepare an fMRI summary plot for the report."""

    input_spec = _ASLSummaryInputSpec
    output_spec = _ASLSummaryOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = fname_presuffix(
            self.inputs.in_func,
            suffix="_aslplot.svg",
            use_ext=False,
            newpath=runtime.cwd,
        )

        dataframe = pd.DataFrame(
            {
                "outliers": np.loadtxt(self.inputs.outliers, usecols=[0]).tolist(),
                # Pick non-standardize dvars (col 1)
                # First timepoint is NaN (difference)
                "DVARS": [np.nan]
                + np.loadtxt(self.inputs.dvars, skiprows=1, usecols=[1]).tolist(),
                # First timepoint is zero (reference volume)
                "FD": [0.0] + np.loadtxt(self.inputs.fd, skiprows=1, usecols=[0]).tolist(),
            }
        )

        fig = ASLPlot(
            self.inputs.in_func,
            mask_file=self.inputs.in_mask if isdefined(self.inputs.in_mask) else None,
            seg_file=self.inputs.in_segm if isdefined(self.inputs.seg_file) else None,
            spikes_files=[self.inputs.in_spikes_bg],
            tr=self.inputs.tr,
            data=dataframe[["outliers", "DVARS", "FD"]],
            units={"outliers": "%", "FD": "mm"},
            vlines={"FD": [self.inputs.fd_thres]},
        ).plot()
        fig.savefig(self._results["out_file"], bbox_inches="tight")
        return runtime


class _CBFSummaryInputSpec(BaseInterfaceInputSpec):
    cbf = File(exists=True, mandatory=True, desc="")
    label = traits.Str(exists=True, mandatory=True, desc="label")
    vmax = traits.Int(exists=True, default_value=90, mandatory=True, desc="max value of asl")
    ref_vol = File(exists=True, mandatory=True, desc="")


class _CBFSummaryOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="written file path")


class CBFSummary(SimpleInterface):
    """Prepare an CBF summary plot for the report."""

    input_spec = _CBFSummaryInputSpec
    output_spec = _CBFSummaryOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = fname_presuffix(
            self.inputs.cbf, suffix="_cbfplot.svg", use_ext=False, newpath=runtime.cwd
        )
        CBFPlot(
            cbf=self.inputs.cbf,
            label=self.inputs.label,
            ref_vol=self.inputs.ref_vol,
            vmax=self.inputs.vmax,
            outfile=self._results["out_file"],
        ).plot()
        # fig.savefig(self._results['out_file'], bbox_inches='tight')
        return runtime


class _CBFtsSummaryInputSpec(BaseInterfaceInputSpec):
    cbf_ts = File(exists=True, mandatory=True, desc=" cbf time series")
    conf_file = File(exists=True, mandatory=False, desc="confound file ")
    score_file = File(exists=True, mandatory=False, desc="scorexindex file ")
    seg_file = File(exists=True, mandatory=True, desc="seg_file")
    tr = traits.Float(desc="TR", mandatory=True)


class _CBFtsSummaryOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="written file path")


class CBFtsSummary(SimpleInterface):
    """Prepare an CBF summary plot for the report."""

    input_spec = _CBFtsSummaryInputSpec
    output_spec = _CBFtsSummaryOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = fname_presuffix(
            self.inputs.cbf_ts, suffix="_cbfcarpetplot.svg", use_ext=False, newpath=runtime.cwd
        )
        fig = CBFtsPlot(
            cbf_file=self.inputs.cbf_ts,
            seg_file=self.inputs.seg_file,
            scoreindex=self.inputs.score_file,
            tr=self.inputs.tr,
        ).plot()
        fig.savefig(self._results["out_file"], bbox_inches="tight")
        return runtime
