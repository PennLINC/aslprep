"""Utility interfaces for ASLPrep."""
import pandas as pd
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.utils.filemanip import fname_presuffix


class _CombineMotionParametersInputSpec(BaseInterfaceInputSpec):
    m0type = traits.Str()
    processing_target = traits.Str()
    aslcontext = File(exists=True)
    control_file_mat_file = traits.Either(File(exists=True), None)
    control_file_par_file = traits.Either(File(exists=True), None)
    control_file_rms_file = traits.Either(File(exists=True), None)
    label_file_mat_file = traits.Either(File(exists=True), None)
    label_file_par_file = traits.Either(File(exists=True), None)
    label_file_rms_file = traits.Either(File(exists=True), None)
    deltam_file_mat_file = traits.Either(File(exists=True), None)
    deltam_file_par_file = traits.Either(File(exists=True), None)
    deltam_file_rms_file = traits.Either(File(exists=True), None)
    cbf_file_mat_file = traits.Either(File(exists=True), None)
    cbf_file_par_file = traits.Either(File(exists=True), None)
    cbf_file_rms_file = traits.Either(File(exists=True), None)
    m0scan_file_mat_file = traits.Either(File(exists=True), None)
    m0scan_file_par_file = traits.Either(File(exists=True), None)
    m0scan_file_rms_file = traits.Either(File(exists=True), None)


class _CombineMotionParametersOutputSpec(TraitedSpec):
    combined_mat_file = File(exists=True)
    combined_par_file = File(exists=True)
    combined_rms_file = File(exists=True)


class CombineMotionParameters(SimpleInterface):
    """Combine motion parameter files from MCFLIRT across image types."""

    input_spec = _CombineMotionParametersInputSpec
    output_spec = _CombineMotionParametersOutputSpec

    def _run_interface(self, runtime):
        aslcontext = pd.read_table(self.inputs.aslcontext)
        files_to_combine = sorted(list(set(aslcontext["volume_type"].tolist())))

        out_par = [None] * aslcontext.shape[0]
        out_mat = [None] * aslcontext.shape[0]
        out_rms = [None] * aslcontext.shape[0]
        for file_to_combine in files_to_combine:
            mat_file = getattr(self.inputs, f"{file_to_combine}_mat_file")
            par_file = getattr(self.inputs, f"{file_to_combine}_par_file")
            rms_file = getattr(self.inputs, f"{file_to_combine}_rms_file")
            idx = aslcontext.loc[
                aslcontext["volume_type"] == file_to_combine.replace("_file", "")
            ].index.values

            with open(par_file, "r") as fo:
                par = fo.readlines()

            with open(mat_file, "r") as fo:
                mat = fo.readlines()

            with open(rms_file, "r") as fo:
                rms = fo.readlines()

            for i_vol, vol_idx in enumerate(idx):
                out_par[vol_idx] = par[i_vol]
                out_mat[vol_idx] = mat[i_vol]
                out_rms[vol_idx] = rms[i_vol]

        self._results["combined_par_file"] = fname_presuffix(
            par_file,
            suffix="_combined",
            newpath=runtime.cwd,
            use_ext=True,
        )
        with open(self._results["combined_par_file"], "w") as fo:
            fo.write(out_par)

        self._results["combined_mat_file"] = fname_presuffix(
            mat_file,
            suffix="_combined",
            newpath=runtime.cwd,
            use_ext=True,
        )
        with open(self._results["combined_mat_file"], "w") as fo:
            fo.write(out_mat)

        self._results["combined_rms_file"] = fname_presuffix(
            rms_file,
            suffix="_combined",
            newpath=runtime.cwd,
            use_ext=True,
        )
        with open(self._results["combined_rms_file"], "w") as fo:
            fo.write(out_rms)

        return runtime
