"""GE-specific interfaces."""
import os

import nibabel as nb
import numpy as np
import pandas as pd
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.utils.filemanip import fname_presuffix

from aslprep.utils.misc import gen_reference


class _GeReferenceFileInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="asl_file")
    in_metadata = traits.Dict(exists=True, mandatory=True, desc="metadata for asl or deltam ")
    bids_dir = traits.Str(exits=True, mandatory=True, desc=" bids directory")
    fwhm = traits.Float(
        exits=False, mandatory=False, default_value=5, desc="smoothing kernel for M0"
    )


class _GeReferenceFileOutputSpec(TraitedSpec):
    ref_file = File(exists=True, mandatory=True, desc="ref file")
    m0_file = File(exists=True, mandatory=True, desc="m0 file")


class GeReferenceFile(SimpleInterface):
    """Generate a reference grid to resample image to new space, but with original resolution."""

    input_spec = _GeReferenceFileInputSpec
    output_spec = _GeReferenceFileOutputSpec

    def _run_interface(self, runtime):
        filex = os.path.abspath(self.inputs.in_file)
        aslcontext1 = filex.replace("_asl.nii.gz", "_aslcontext.tsv")
        aslcontext = pd.read_csv(aslcontext1)
        vol_types = aslcontext["volume_type"].tolist()
        control_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "control"]
        m0_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "m0scan"]
        deltam_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "deltam"]
        cbf_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "CBF"]

        allasl = nb.load(self.inputs.in_file)
        dataasl = allasl.get_fdata()

        if self.inputs.in_metadata["M0Type"] == "Separate":
            m0file = self.inputs.in_file.replace("asl.nii.gz", "m0scan.nii.gz")
            reffile = gen_reference(m0file, fwhm=self.inputs.fwhm, newpath=runtime.cwd)
            m0file = reffile

        elif self.inputs.in_metadata["M0Type"] == "Included":
            modata2 = dataasl[:, :, :, m0_volume_idx]
            m0filename = fname_presuffix(
                self.inputs.in_file, suffix="_mofile", newpath=os.getcwd()
            )
            m0obj = nb.Nifti1Image(modata2, allasl.affine, allasl.header)
            m0obj.to_filename(m0filename)
            reffile = gen_reference(m0filename, fwhm=self.inputs.fwhm, newpath=runtime.cwd)
            m0file = reffile

        elif self.inputs.in_metadata["M0Type"] == "Estimate":
            m0num = float(self.inputs.in_metadata["M0Estimate"])
            if len(deltam_volume_idx) > 0:
                modata2 = dataasl[:, :, :, deltam_volume_idx]
            elif len(control_volume_idx) > 0:
                modata2 = dataasl[:, :, :, control_volume_idx]
            elif len(cbf_volume_idx) > 0:
                modata2 = dataasl[:, :, :, cbf_volume_idx]

            m0filename = fname_presuffix(
                self.inputs.in_file, suffix="_m0file", newpath=os.getcwd()
            )
            if len(modata2.shape) > 3:
                mdata = np.mean(modata2, axis=3)
            else:
                mdata = modata2

            m0file_data = m0num * np.ones_like(mdata)

            m0obj = nb.Nifti1Image(m0file_data, allasl.affine, allasl.header)
            m0obj.to_filename(m0filename)
            m0file = m0filename

            reffilename = fname_presuffix(
                self.inputs.in_file, suffix="_refile", newpath=os.getcwd()
            )
            refobj = nb.Nifti1Image(mdata, allasl.affine, allasl.header)

            refobj.to_filename(reffilename)
            reffile = gen_reference(reffilename, fwhm=0, newpath=runtime.cwd)

        elif self.inputs.in_metadata["M0Type"] == "Absent":
            if len(control_volume_idx) > 0:
                modata2 = dataasl[:, :, :, control_volume_idx]
            m0filename = fname_presuffix(
                self.inputs.in_file, suffix="_m0file", newpath=os.getcwd()
            )
            if len(modata2.shape) > 3:
                mdata = np.mean(modata2, axis=3)
            else:
                mdata = modata2

            m0obj = nb.Nifti1Image(mdata, allasl.affine, allasl.header)
            m0obj.to_filename(m0filename)
            m0file = m0filename

            reffilename = fname_presuffix(
                self.inputs.in_file, suffix="_refile", newpath=os.getcwd()
            )
            refobj = nb.Nifti1Image(mdata, allasl.affine, allasl.header)

            refobj.to_filename(reffilename)
            reffile = gen_reference(
                reffilename,
                fwhm=self.inputs.fwhm,
                newpath=runtime.cwd,
            )

        else:
            raise RuntimeError("no path way to obtain real m0scan")

        self._results["ref_file"] = reffile
        self._results["m0_file"] = m0file

        return runtime