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
    control_mat_file = traits.Either(traits.List(File(exists=True)), None)
    control_par_file = traits.Either(File(exists=True), None)
    control_rms_file = traits.Either(File(exists=True), None)
    label_mat_file = traits.Either(traits.List(File(exists=True)), None)
    label_par_file = traits.Either(File(exists=True), None)
    label_rms_file = traits.Either(File(exists=True), None)
    deltam_mat_file = traits.Either(traits.List(File(exists=True)), None)
    deltam_par_file = traits.Either(File(exists=True), None)
    deltam_rms_file = traits.Either(File(exists=True), None)
    cbf_mat_file = traits.Either(traits.List(File(exists=True)), None)
    cbf_par_file = traits.Either(File(exists=True), None)
    cbf_rms_file = traits.Either(File(exists=True), None)
    m0scan_mat_file = traits.Either(traits.List(File(exists=True)), None)
    m0scan_par_file = traits.Either(File(exists=True), None)
    m0scan_rms_file = traits.Either(File(exists=True), None)


class _CombineMotionParametersOutputSpec(TraitedSpec):
    combined_mat_file = traits.List(File(exists=True))
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
        out_mat_files = [None] * aslcontext.shape[0]
        out_rms = [None] * aslcontext.shape[0]
        for file_to_combine in files_to_combine:
            mat_files = getattr(self.inputs, f"{file_to_combine}_mat_file")
            par_file = getattr(self.inputs, f"{file_to_combine}_par_file")
            rms_file = getattr(self.inputs, f"{file_to_combine}_rms_file")
            idx = aslcontext.loc[aslcontext["volume_type"] == file_to_combine].index.values

            with open(par_file, "r") as fo:
                par = fo.readlines()

            with open(rms_file, "r") as fo:
                rms = fo.readlines()

            for i_vol, vol_idx in enumerate(idx):
                out_par[vol_idx] = par[i_vol]
                out_mat_files[vol_idx] = mat_files[i_vol]
                out_rms[vol_idx] = rms[i_vol]

        self._results["combined_par_file"] = fname_presuffix(
            par_file,
            suffix="_combined",
            newpath=runtime.cwd,
            use_ext=True,
        )
        with open(self._results["combined_par_file"], "w") as fo:
            fo.write("".join(out_par))

        self._results["combined_mat_file"] = out_mat_files

        self._results["combined_rms_file"] = fname_presuffix(
            rms_file,
            suffix="_combined",
            newpath=runtime.cwd,
            use_ext=True,
        )
        with open(self._results["combined_rms_file"], "w") as fo:
            fo.write("".join(out_rms))

        return runtime