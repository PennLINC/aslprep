"""Plotting interfaces."""
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
    in_func = File(exists=True, mandatory=True, desc="input ASL time-series (4D file)")
    in_mask = File(exists=True, desc="3D brain mask")
    in_segm = File(exists=True, desc="resampled segmentation")
    confounds_file = File(exists=True, desc="BIDS' _confounds.tsv file")

    str_or_tuple = traits.Either(
        traits.Str,
        traits.Tuple(traits.Str, traits.Either(None, traits.Str)),
        traits.Tuple(traits.Str, traits.Either(None, traits.Str), traits.Either(None, traits.Str)),
    )
    confounds_list = traits.List(
        str_or_tuple,
        minlen=1,
        desc="list of headers to extract from the confounds_file",
    )
    tr = traits.Either(None, traits.Float, usedefault=True, desc="the repetition time")


class _ASLSummaryOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="written file path")


class ASLSummary(SimpleInterface):
    """Copy the x-form matrices from `hdr_file` to `out_file`.

    Clearly that's wrong.
    """

    input_spec = _ASLSummaryInputSpec
    output_spec = _ASLSummaryOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = fname_presuffix(
            self.inputs.in_func,
            suffix="_aslplot.svg",
            use_ext=False,
            newpath=runtime.cwd,
        )

        dataframe = pd.read_csv(
            self.inputs.confounds_file,
            sep="\t",
            index_col=None,
            dtype="float32",
            na_filter=True,
            na_values="n/a",
        )

        headers = []
        units = {}
        names = {}

        for conf_el in self.inputs.confounds_list:
            if isinstance(conf_el, (list, tuple)):
                headers.append(conf_el[0])
                if conf_el[1] is not None:
                    units[conf_el[0]] = conf_el[1]

                if len(conf_el) > 2 and conf_el[2] is not None:
                    names[conf_el[0]] = conf_el[2]
            else:
                headers.append(conf_el)

        if not headers:
            data = None
            units = None
        else:
            data = dataframe[headers]

        colnames = data.columns.ravel().tolist()

        for name, newname in list(names.items()):
            colnames[colnames.index(name)] = newname

        data.columns = colnames

        fig = ASLPlot(
            self.inputs.in_func,
            mask_file=self.inputs.in_mask if isdefined(self.inputs.in_mask) else None,
            seg_file=(self.inputs.in_segm if isdefined(self.inputs.in_segm) else None),
            tr=self.inputs.tr,
            data=data,
            units=units,
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
            self.inputs.cbf,
            suffix="_cbfplot.svg",
            use_ext=False,
            newpath=runtime.cwd,
        )
        CBFPlot(
            cbf=self.inputs.cbf,
            label=self.inputs.label,
            ref_vol=self.inputs.ref_vol,
            vmax=self.inputs.vmax,
            outfile=self._results["out_file"],
        ).plot()
        return runtime


class _CBFtsSummaryInputSpec(BaseInterfaceInputSpec):
    cbf_ts = File(exists=True, mandatory=True, desc=" cbf time series")
    confounds_file = File(exists=True, mandatory=False, desc="confound file ")
    score_outlier_index = File(exists=True, mandatory=False, desc="scorexindex file ")
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
            self.inputs.cbf_ts,
            suffix="_cbfcarpetplot.svg",
            use_ext=False,
            newpath=runtime.cwd,
        )
        fig = CBFtsPlot(
            cbf_file=self.inputs.cbf_ts,
            seg_file=self.inputs.seg_file,
            score_outlier_index=self.inputs.score_outlier_index,
            tr=self.inputs.tr,
        ).plot()
        fig.savefig(self._results["out_file"], bbox_inches="tight")
        return runtime


class _CBFByTissueTypePlotInputSpec(BaseInterfaceInputSpec):
    cbf = File(exists=True, mandatory=True, desc="")
    seg_file = File(exists=True, mandatory=True, desc="Segmentation file")


class _CBFByTissueTypePlotOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="written file path")


class CBFByTissueTypePlot(SimpleInterface):
    """Prepare an CBF summary plot for the report."""

    input_spec = _CBFByTissueTypePlotInputSpec
    output_spec = _CBFByTissueTypePlotOutputSpec

    def _run_interface(self, runtime):
        import matplotlib.pyplot as plt
        import seaborn as sns
        from nilearn import image, masking

        self._results["out_file"] = fname_presuffix(
            self.inputs.cbf,
            suffix="_cbfplot.svg",
            use_ext=False,
            newpath=runtime.cwd,
        )

        dfs = []
        for i_tissue_type, tissue_type in enumerate(["GM", "WM", "CSF"]):
            tissue_type_val = i_tissue_type + 1
            mask_img = image.math_img(
                f"(img == {tissue_type_val}).astype(int)",
                img=self.inputs.seg_file,
            )
            tissue_type_vals = masking.apply_mask(self.inputs.cbf, mask_img)
            df = pd.DataFrame(
                columns=["CBF\n(mL/100 g/min)", "Tissue Type"],
                data=list(
                    map(list, zip(*[tissue_type_vals, [tissue_type] * tissue_type_vals.size]))
                ),
            )
            dfs.append(df)

        df = pd.concat(dfs, axis=0)

        # Create the plot
        with sns.axes_style("whitegrid"), sns.plotting_context(font_scale=3):
            fig, ax = plt.subplots(figsize=(16, 8))
            sns.despine(ax=ax, bottom=True, left=True)
            sns.boxenplot(
                x="Tissue Type",
                y="CBF\n(mL/100 g/min)",
                data=df,
                width=0.6,
                showfliers=True,
                palette={"GM": "#1b60a5", "WM": "#2da467", "CSF": "#9d8f25"},
                ax=ax,
            )
            fig.tight_layout()
            fig.savefig(self._results["out_file"])

        return runtime
