"""Plotting interfaces."""
import nibabel as nb
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
from niworkflows.utils.timeseries import _cifti_timeseries, _nifti_timeseries

from aslprep.utils.plotting import CBFPlot, fMRIPlot


class _ASLCarpetPlotInputSpec(BaseInterfaceInputSpec):
    in_nifti = File(exists=True, mandatory=True, desc="input BOLD (4D NIfTI file)")
    in_cifti = File(exists=True, desc="input BOLD (CIFTI dense timeseries)")
    in_segm = File(exists=True, desc="volumetric segmentation corresponding to in_nifti")
    confounds_file = File(exists=True, desc="BIDS' _confounds.tsv file")

    str_or_tuple = traits.Either(
        traits.Str,
        traits.Tuple(traits.Str, traits.Either(None, traits.Str)),
        traits.Tuple(traits.Str, traits.Either(None, traits.Str), traits.Either(None, traits.Str)),
    )
    confounds_list = traits.Either(
        traits.List(str_or_tuple, minlen=1),
        None,
        desc="list of headers to extract from the confounds_file",
    )
    tr = traits.Either(None, traits.Float, usedefault=True, desc="the repetition time")
    drop_trs = traits.Int(0, usedefault=True, desc="dummy scans")


class _ASLCarpetPlotOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="written file path")


class ASLCarpetPlot(SimpleInterface):
    """Create a combined carpet/line plot for ASL time series data.

    Copy the x-form matrices from `hdr_file` to `out_file`.
    """

    input_spec = _ASLCarpetPlotInputSpec
    output_spec = _ASLCarpetPlotOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = fname_presuffix(
            self.inputs.in_nifti, suffix="_fmriplot.svg", use_ext=False, newpath=runtime.cwd
        )

        has_cifti = isdefined(self.inputs.in_cifti)

        # Read input object and create timeseries + segments object
        seg_file = self.inputs.in_segm if isdefined(self.inputs.in_segm) else None
        dataset, segments = _nifti_timeseries(
            nb.load(self.inputs.in_nifti),
            nb.load(seg_file),
            remap_rois=False,
            labels=(
                ("WM+CSF", "Edge")
                if has_cifti
                else ("Ctx GM", "dGM", "sWM+sCSF", "dWM+dCSF", "Cb", "Edge")
            ),
        )

        # Process CIFTI
        if has_cifti:
            cifti_data, cifti_segments = _cifti_timeseries(nb.load(self.inputs.in_cifti))

            if seg_file is not None:
                # Append WM+CSF and Edge masks
                cifti_length = cifti_data.shape[0]
                dataset = np.vstack((cifti_data, dataset))
                segments = {k: np.array(v) + cifti_length for k, v in segments.items()}
                cifti_segments.update(segments)
                segments = cifti_segments
            else:
                dataset, segments = cifti_data, cifti_segments

        dataframe = pd.read_table(
            self.inputs.confounds_file,
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
            data = data.rename(columns=names)

        fig = fMRIPlot(
            dataset,
            segments=segments,
            tr=self.inputs.tr,
            confounds=data,
            units=units,
            nskip=self.inputs.drop_trs,
            paired_carpet=has_cifti,
        ).plot()
        fig.savefig(self._results["out_file"], bbox_inches="tight")
        fig.clf()
        return runtime


class _CBFSummaryPlotInputSpec(BaseInterfaceInputSpec):
    cbf = File(exists=True, mandatory=True, desc="")
    label = traits.Str(mandatory=True, desc="label")
    vmax = traits.Int(default_value=90, mandatory=True, desc="max value of asl")
    ref_vol = File(exists=True, mandatory=True, desc="")


class _CBFSummaryPlotOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="written file path")


class CBFSummaryPlot(SimpleInterface):
    """Prepare an CBF summary plot for the report.

    This plot restricts CBF values to -20 (if there are negative values) or 0 (if not) to 100.
    """

    input_spec = _CBFSummaryPlotInputSpec
    output_spec = _CBFSummaryPlotOutputSpec

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
                y="CBF\n(mL/100 g/min)",
                data=df,
                width=0.6,
                showfliers=True,
                palette={"GM": "#1b60a5", "WM": "#2da467", "CSF": "#9d8f25"},
                hue="Tissue Type",
                legend=False,
                ax=ax,
            )
            fig.tight_layout()
            fig.savefig(self._results["out_file"])
            fig.clf()

        return runtime
