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
from niworkflows.viz.plots import fMRIPlot


class _ASLSummaryInputSpec(BaseInterfaceInputSpec):
    in_nifti = File(exists=True, mandatory=True, desc="input BOLD (4D NIfTI file)")
    in_cifti = File(exists=True, desc="input BOLD (CIFTI dense timeseries)")
    in_segm = File(exists=True, desc="volumetric segmentation corresponding to in_nifti")
    confounds_file = File(exists=True, desc="BIDS' _confounds.tsv file")

    str_or_tuple = traits.Either(
        traits.Str,
        traits.Tuple(traits.Str, traits.Either(None, traits.Str)),
        traits.Tuple(traits.Str, traits.Either(None, traits.Str), traits.Either(None, traits.Str)),
    )
    confounds_list = traits.List(
        str_or_tuple, minlen=1, desc="list of headers to extract from the confounds_file"
    )
    tr = traits.Either(None, traits.Float, usedefault=True, desc="the repetition time")
    drop_trs = traits.Int(0, usedefault=True, desc="dummy scans")


class _ASLSummaryOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="written file path")


class ASLSummary(SimpleInterface):
    """Copy the x-form matrices from `hdr_file` to `out_file`."""

    input_spec = _ASLSummaryInputSpec
    output_spec = _ASLSummaryOutputSpec

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

        data = data.rename(columns=names)

        fig = fMRIPlot(
            dataset,
            segments=segments,
            tr=self.inputs.tr,
            confounds=data,
            units=units,
            nskip=self.inputs.drop_trs,
            paired_carpet=has_cifti,
            # The main change from fMRIPrep's usage is that detrend is False for ASL.
            detrend=False,
        ).plot()
        fig.savefig(self._results["out_file"], bbox_inches="tight")
        return runtime
