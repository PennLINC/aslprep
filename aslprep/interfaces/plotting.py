"""Visualization tools."""
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.utils.filemanip import fname_presuffix

from aslprep.utils.viz import CBFPlot, CBFtsPlot


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
