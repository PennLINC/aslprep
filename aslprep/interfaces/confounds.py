# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Interfaces for calculating and collecting confounds."""
import pandas as pd
from nipype import logging
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    isdefined,
    traits,
)
from nipype.utils.filemanip import fname_presuffix

from aslprep.niworkflows.viz.plots import ASLPlot
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
        str_or_tuple, minlen=1, desc="list of headers to extract from the confounds_file"
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
            self.inputs.in_func, suffix="_aslplot.svg", use_ext=False, newpath=runtime.cwd
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
