"""Utility interfaces for ASLPrep."""
import os

import nibabel as nb
import pandas as pd
from nilearn import image
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec
from nipype.utils.filemanip import fname_presuffix, load_json, save_json

from aslprep import config
from aslprep.utils.misc import reduce_metadata_lists


class _ReduceASLFilesInputSpec(BaseInterfaceInputSpec):
    asl_file = File(exists=True, mandatory=True, desc="ASL file to split.")
    aslcontext = File(exists=True, mandatory=True, desc="aslcontext TSV.")
    processing_target = traits.Str()
    metadata = traits.Dict()


class _ReduceASLFilesOutputSpec(TraitedSpec):
    asl_file = File(exists=True, desc="Modified ASL file.")
    aslcontext = File(exists=True, desc="Modified aslcontext file.")
    metadata = traits.Dict()


class ReduceASLFiles(SimpleInterface):
    """Split ASL data into different files."""

    input_spec = _ReduceASLFilesInputSpec
    output_spec = _ReduceASLFilesOutputSpec

    def _run_interface(self, runtime):
        aslcontext = pd.read_table(self.inputs.aslcontext)
        asl_img = nb.load(self.inputs.asl_file)
        assert asl_img.shape[3] == aslcontext.shape[0]

        if self.inputs.processing_target == "controllabel":
            files_to_keep = ["control", "label", "m0scan"]
        elif self.inputs.processing_target == "deltam":
            files_to_keep = ["deltam", "m0scan"]
        else:
            files_to_keep = ["cbf", "m0scan"]

        asl_idx = aslcontext.loc[aslcontext["volume_type"].isin(files_to_keep)].index.values
        asl_idx = asl_idx.astype(int)
        self._results["metadata"] = reduce_metadata_lists(self.inputs.metadata, asl_idx)

        asl_img = image.index_img(asl_img, asl_idx)

        self._results["asl_file"] = fname_presuffix(
            self.inputs.asl_file,
            suffix="_reduced",
            newpath=runtime.cwd,
            use_ext=True,
        )
        asl_img.to_filename(self._results["asl_file"])

        aslcontext = aslcontext.loc[asl_idx]
        self._results["aslcontext"] = fname_presuffix(
            self.inputs.aslcontext,
            suffix="_reduced",
            newpath=runtime.cwd,
            use_ext=True,
        )
        aslcontext.to_csv(self._results["aslcontext"], sep="\t", index=False)

        return runtime


class _RMSDiffInputSpec(FSLCommandInputSpec):
    matrixfile1 = File(
        exists=True,
        position=0,
        argstr="%s",
        desc="First matrix file.",
        mandatory=True,
    )
    matrixfile2 = File(
        exists=True,
        position=1,
        argstr="%s",
        desc="Second matrix file.",
        mandatory=True,
    )
    ref_vol = File(
        exists=True,
        position=2,
        argstr="%s",
        desc="Reference volume.",
        mandatory=True,
    )


class _RMSDiffOutputSpec(TraitedSpec):
    rmsd = traits.Float()


class RMSDiff(FSLCommand):
    """Run rmsdiff."""

    _cmd = "rmsdiff"
    input_spec = _RMSDiffInputSpec
    output_spec = _RMSDiffOutputSpec

    def aggregate_outputs(self, runtime=None, needed_outputs=None):  # noqa: U100
        """Taken from nipype.interfaces.afni.preprocess.ClipLevel."""
        outputs = self._outputs()

        outfile = os.path.join(os.getcwd(), "stat_result.json")

        if runtime is None:
            try:
                rmsd = load_json(outfile)["stat"]
            except IOError:
                return self.run().outputs
        else:
            rmsd = []
            for line in runtime.stdout.split("\n"):
                if line:
                    values = line.split()
                    if len(values) > 1:
                        rmsd.append([float(val) for val in values])
                    else:
                        rmsd.extend([float(val) for val in values])

            if len(rmsd) == 1:
                rmsd = rmsd[0]

            save_json(outfile, dict(stat=rmsd))

        outputs.rmsd = rmsd

        return outputs


class _PairwiseRMSDiffInputSpec(BaseInterfaceInputSpec):
    in_files = traits.List(
        File(exists=True),
        desc="Matrix files to compare with each other.",
        mandatory=True,
    )
    ref_file = File(
        exists=True,
        desc="Reference volume.",
        mandatory=True,
    )


class _PairwiseRMSDiffOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="Output txt file.")


class PairwiseRMSDiff(SimpleInterface):
    """Run rmsdiff on each contiguous pair of transform files to build a txt file of rmsd values.

    This interface uses :class:`~aslprep.interfaces.utility.RMSDiff` internally, which may not be
    a proper nipype pattern.
    """

    input_spec = _PairwiseRMSDiffInputSpec
    output_spec = _PairwiseRMSDiffOutputSpec

    def _run_interface(self, runtime):
        rmsd = []
        for i_file in range(len(self.inputs.in_files) - 1):
            j_file = i_file + 1
            file1 = self.inputs.in_files[i_file]
            file2 = self.inputs.in_files[j_file]
            # run rmsdiff
            rmsdiff = RMSDiff(matrixfile1=file1, matrixfile2=file2, ref_vol=self.inputs.ref_file)
            res = rmsdiff.run()
            assert isinstance(res.outputs.rmsd, float)
            rmsd.append(str(res.outputs.rmsd))

        self._results["out_file"] = fname_presuffix(
            self.inputs.ref_file,
            suffix="_rmsd.txt",
            newpath=runtime.cwd,
            use_ext=False,
        )
        with open(self._results["out_file"], "w") as fo:
            fo.write("\n".join(rmsd))

        return runtime


class _CombineMotionParametersInputSpec(BaseInterfaceInputSpec):
    m0type = traits.Str()
    processing_target = traits.Str()
    aslcontext = File(exists=True)
    control_mat_files = traits.Either(traits.List(File(exists=True)), None)
    control_par_file = traits.Either(File(exists=True), None)
    label_mat_files = traits.Either(traits.List(File(exists=True)), None)
    label_par_file = traits.Either(File(exists=True), None)
    deltam_mat_files = traits.Either(traits.List(File(exists=True)), None)
    deltam_par_file = traits.Either(File(exists=True), None)
    cbf_mat_files = traits.Either(traits.List(File(exists=True)), None)
    cbf_par_file = traits.Either(File(exists=True), None)
    m0scan_mat_files = traits.Either(traits.List(File(exists=True)), None)
    m0scan_par_file = traits.Either(File(exists=True), None)


class _CombineMotionParametersOutputSpec(TraitedSpec):
    mat_file_list = traits.List(File(exists=True))
    combined_par_file = File(exists=True)


class CombineMotionParameters(SimpleInterface):
    """Combine motion parameter files from MCFLIRT across image types."""

    input_spec = _CombineMotionParametersInputSpec
    output_spec = _CombineMotionParametersOutputSpec

    def _run_interface(self, runtime):
        aslcontext = pd.read_table(self.inputs.aslcontext)
        files_to_combine = sorted(list(set(aslcontext["volume_type"].tolist())))

        out_par = [None] * aslcontext.shape[0]
        out_mat_files = [None] * aslcontext.shape[0]
        for file_to_combine in files_to_combine:
            mat_files = getattr(self.inputs, f"{file_to_combine}_mat_files")
            par_file = getattr(self.inputs, f"{file_to_combine}_par_file")
            idx = aslcontext.loc[aslcontext["volume_type"] == file_to_combine].index.values

            with open(par_file, "r") as fo:
                par = fo.readlines()

            for i_vol, vol_idx in enumerate(idx):
                out_par[vol_idx] = par[i_vol]
                out_mat_files[vol_idx] = mat_files[i_vol]

        self._results["combined_par_file"] = fname_presuffix(
            par_file,
            suffix="_combined",
            newpath=runtime.cwd,
            use_ext=True,
        )
        with open(self._results["combined_par_file"], "w") as fo:
            fo.write("".join(out_par))

        self._results["mat_file_list"] = out_mat_files

        return runtime


class _SplitOutVolumeTypeInputSpec(BaseInterfaceInputSpec):
    volumetype = traits.Str()
    aslcontext = File(exists=True)
    asl_file = File(exists=True)


class _SplitOutVolumeTypeOutputSpec(TraitedSpec):
    out_file = File(exists=True)


class SplitOutVolumeType(SimpleInterface):
    """Split out a specific volume type from the ASL file."""

    input_spec = _SplitOutVolumeTypeInputSpec
    output_spec = _SplitOutVolumeTypeOutputSpec

    def _run_interface(self, runtime):
        aslcontext = pd.read_table(self.inputs.aslcontext)
        volumetype_df = aslcontext.loc[aslcontext["volume_type"] == self.inputs.volumetype]
        volumetype_idx = volumetype_df.index.tolist()
        if len(volumetype_idx) == 0:
            raise ValueError(f"No volumes found for {self.inputs.volumetype}")

        out_img = image.index_img(self.inputs.asl_file, volumetype_idx)
        self._results["out_file"] = fname_presuffix(
            self.inputs.asl_file,
            suffix=f"_{self.inputs.volumetype}",
            newpath=runtime.cwd,
            use_ext=True,
        )
        out_img.to_filename(self._results["out_file"])

        return runtime


class _SplitReferenceTargetInputSpec(BaseInterfaceInputSpec):
    aslcontext = File(exists=True, required=True)
    asl_file = File(exists=True, required=True)


class _SplitReferenceTargetOutputSpec(TraitedSpec):
    out_file = File(exists=True)


class SplitReferenceTarget(SimpleInterface):
    """Split out just the optimal volume type for reference files from the overall ASL file.

    This means grabbing M0 volumes if they're available, or control volumes if *they're* available,
    or deltams, or cbfs.
    """

    input_spec = _SplitReferenceTargetInputSpec
    output_spec = _SplitReferenceTargetOutputSpec

    def _run_interface(self, runtime):
        aslcontext = pd.read_table(self.inputs.aslcontext)
        volume_types = aslcontext["volume_type"].values
        if "m0scan" in volume_types:
            ref_target = "m0scan"
        elif "control" in volume_types:
            ref_target = "control"
        elif "deltam" in volume_types:
            ref_target = "deltam"
        elif "cbf" in volume_types:
            ref_target = "cbf"
        else:
            raise ValueError(volume_types)

        config.loggers.interface.warning(f"Selected {ref_target} for reference.")

        volumetype_idx = aslcontext["volume_type"].loc[volume_types == ref_target].index.tolist()
        out_img = image.index_img(self.inputs.asl_file, volumetype_idx)
        self._results["out_file"] = fname_presuffix(
            self.inputs.asl_file,
            suffix=f"_{ref_target}",
            newpath=runtime.cwd,
            use_ext=True,
        )
        out_img.to_filename(self._results["out_file"])

        return runtime