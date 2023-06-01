"""GE-specific interfaces."""
import os

import nibabel as nb
import numpy as np
import pandas as pd
from nibabel.processing import smooth_image
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.utils.filemanip import fname_presuffix


class _GeReferenceFileInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="asl_file")
    metadata = traits.Dict(mandatory=True, desc="metadata for asl or deltam ")
    fwhm = traits.Float(
        mandatory=False,
        default_value=5,
        desc="smoothing kernel for M0",
    )
    m0scan = traits.Either(
        File(exists=True),
        None,
        mandatory=True,
        desc="M0 file.",
    )
    m0scan_metadata = traits.Either(
        traits.Dict,
        None,
        mandatory=True,
        desc="metadata for M0 scan. Only defined if M0Type is 'Separate'.",
    )
    aslcontext = File(exists=True, mandatory=True, desc="aslcontext file.")


class _GeReferenceFileOutputSpec(TraitedSpec):
    ref_file = File(exists=True, mandatory=True, desc="ref file")
    m0_file = File(exists=True, mandatory=True, desc="Averaged and smoothed m0 file")
    m0tr = traits.Either(
        traits.Float,
        None,
        desc="RepetitionTimePreparation for M0 scans.",
    )


class GeReferenceFile(SimpleInterface):
    """Generate a reference grid to resample image to new space, but with original resolution.

    TODO: Extract M0 in a separate step, to be shared between the GE and non-GE pipelines.
    """

    input_spec = _GeReferenceFileInputSpec
    output_spec = _GeReferenceFileOutputSpec

    def _run_interface(self, runtime):
        metadata = self.inputs.metadata
        aslcontext_df = pd.read_table(self.inputs.aslcontext)
        vol_types = aslcontext_df["volume_type"].tolist()
        control_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "control"]
        deltam_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "deltam"]
        cbf_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "cbf"]

        asl_img = nb.load(self.inputs.in_file)
        asl_data = asl_img.get_fdata()

        m0_file = fname_presuffix(self.inputs.in_file, suffix="_m0file", newpath=os.getcwd())
        ref_file = fname_presuffix(self.inputs.in_file, suffix="_ref", newpath=os.getcwd())

        if metadata["M0Type"] == "Separate":
            # Average and smooth the M0 data
            m0_img = nb.load(self.inputs.m0scan)
            m0_data = m0_img.get_fdata()
            if m0_data.ndim > 3:
                m0_data = np.mean(m0_data, axis=3)

            mean_img = nb.Nifti1Image(m0_data, asl_img.affine, asl_img.header)
            smoothed_img = smooth_image(mean_img, fwhm=self.inputs.fwhm)
            smoothed_img.to_filename(m0_file)

            self._results["m0_file"] = m0_file
            self._results["ref_file"] = m0_file  # The reference file is the averaged, smoothed M0

            m0tr = self.inputs.m0scan_metadata["RepetitionTimePreparation"]
            if np.array(m0tr).size > 1 and np.std(m0tr) > 0:
                raise ValueError("M0 scans have variable TR. ASLPrep does not support this.")

        elif metadata["M0Type"] == "Included":
            m0_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "m0scan"]
            m0_data = asl_data[:, :, :, m0_volume_idx]

            # Average and smooth the M0 data
            mean_m0_data = np.mean(m0_data, axis=3)
            mean_img = nb.Nifti1Image(mean_m0_data, asl_img.affine, asl_img.header)
            smoothed_img = smooth_image(mean_img, fwhm=self.inputs.fwhm)
            smoothed_img.to_filename(m0_file)

            self._results["m0_file"] = m0_file
            self._results["ref_file"] = m0_file  # The reference file is the averaged, smoothed M0

            if np.array(metadata["RepetitionTimePreparation"]).size > 1:
                m0tr = np.array(metadata["RepetitionTimePreparation"])[m0_volume_idx]
            else:
                m0tr = metadata["RepetitionTimePreparation"]

            if np.array(m0tr).size > 1 and np.std(m0tr) > 0:
                raise ValueError("M0 scans have variable TR. ASLPrep does not support this.")

        elif metadata["M0Type"] == "Estimate":
            if deltam_volume_idx:
                idx = deltam_volume_idx
            elif control_volume_idx:
                idx = control_volume_idx
            elif cbf_volume_idx:
                idx = cbf_volume_idx

            selected_data = asl_data[:, :, :, idx]
            mean_data = np.mean(selected_data, axis=3)

            # XXX: Why not use an existing brain mask here?
            m0_data = metadata["M0Estimate"] * np.ones(mean_data.shape)
            m0_img = nb.Nifti1Image(m0_data, asl_img.affine, asl_img.header)
            m0_img.to_filename(m0_file)

            # The reference image is the mean deltaM, or control, or CBF data.
            ref_img = nb.Nifti1Image(mean_data, asl_img.affine, asl_img.header)
            ref_img.to_filename(ref_file)

            self._results["m0_file"] = m0_file
            self._results["ref_file"] = ref_file

            m0tr = None

        elif metadata["M0Type"] == "Absent":
            if not control_volume_idx:
                raise RuntimeError("M0 could not be estimated from control volumes.")

            control_data = asl_data[:, :, :, control_volume_idx]
            m0_data = np.mean(control_data, axis=3)

            mean_img = nb.Nifti1Image(m0_data, asl_img.affine, asl_img.header)
            mean_img.to_filename(m0_file)  # XXX: Why not smooth the "M0" image?
            smoothed_img = smooth_image(mean_img, fwhm=self.inputs.fwhm)
            smoothed_img.to_filename(ref_file)

            self._results["m0_file"] = m0_file
            self._results["ref_file"] = ref_file

            m0tr = None

        else:
            raise RuntimeError("no path way to obtain real m0scan")

        self._results["m0tr"] = m0tr

        return runtime
