"""Interfaces originally contained in the aslprep niworkflows copy.

These interfaces must be removed and replaced with imports of current niworkflows
interfaces ASAP.
"""
import os

import nibabel as nb
import numpy as np
from nipype.interfaces.base import (
    traits,
    isdefined,
    TraitedSpec,
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    InputMultiObject,
)
from nipype.interfaces import fsl, afni


class _EstimateReferenceImageInputSpec(BaseInterfaceInputSpec):
    in_file = InputMultiObject(
        File(exists=True),
        mandatory=True,
        desc=(
            "4D EPI file. If multiple files "
            "are provided, they are assumed "
            "to represent multiple echoes "
            "from the same run."
        ),
    )
    sbref_file = InputMultiObject(
        File(exists=True),
        desc=(
            "Single band reference image. "
            "If multiple files are provided, "
            "they are assumed to represent "
            "multiple echoes."
        ),
    )
    mc_method = traits.Enum(
        "AFNI",
        "FSL",
        usedefault=True,
        desc="Which software to use to perform motion correction",
    )
    multiecho = traits.Bool(
        False,
        usedefault=True,
        desc=(
            "If multiecho data was supplied, data from "
            "the first echo will be selected."
        ),
    )


class _EstimateReferenceImageOutputSpec(TraitedSpec):
    ref_image = File(exists=True, desc="3D reference image")
    n_volumes_to_discard = traits.Int(
        desc="Number of detected non-steady "
        "state volumes in the beginning of "
        "the input file"
    )


class EstimateReferenceImage(SimpleInterface):
    """
    Generate a reference 3D map from BOLD and SBRef EPI images for BOLD datasets.

    Given a 4D BOLD file or one or more 3/4D SBRefs, estimate a reference
    image for subsequent motion estimation and coregistration steps.
    For the case of BOLD datasets, it estimates a number of T1w saturated volumes
    (non-steady state at the beginning of the scan) and calculates the median
    across them.
    Otherwise (SBRefs or detected zero non-steady state frames), a median of
    of a subset of motion corrected volumes is used.
    If the input reference (BOLD or SBRef) is 3D already, it just returns a
    copy of the image with the NIfTI header extensions removed.

    LIMITATION: If one wants to extract the reference from several SBRefs
    with several echoes each, the first echo should be selected elsewhere
    and run this interface in ``multiecho = False`` mode.
    """

    input_spec = _EstimateReferenceImageInputSpec
    output_spec = _EstimateReferenceImageOutputSpec

    def _run_interface(self, runtime):
        is_sbref = isdefined(self.inputs.sbref_file)
        ref_input = self.inputs.sbref_file if is_sbref else self.inputs.in_file

        if self.inputs.multiecho:
            if len(ref_input) < 2:
                input_name = "sbref_file" if is_sbref else "in_file"
                raise ValueError("Argument 'multiecho' is True but "
                                 f"'{input_name}' has only one element.")
            else:
                # Select only the first echo (see LIMITATION above for SBRefs)
                ref_input = ref_input[:1]
        elif not is_sbref and len(ref_input) > 1:
            raise ValueError("Input 'in_file' cannot point to more than one file "
                             "for single-echo BOLD datasets.")

        # Build the nibabel spatial image we will work with
        ref_im = []
        for im_i in ref_input:
            nib_i = nb.squeeze_image(nb.load(im_i))
            if nib_i.dataobj.ndim == 3:
                ref_im.append(nib_i)
            elif nib_i.dataobj.ndim == 4:
                ref_im += nb.four_to_three(nib_i)
        ref_im = nb.squeeze_image(nb.concat_images(ref_im))

        # Volumes to discard only makes sense with BOLD inputs.
        if not is_sbref:
            n_volumes_to_discard = _get_vols_to_discard(ref_im)
            out_ref_fname = os.path.join(runtime.cwd, "ref_bold.nii.gz")
        else:
            n_volumes_to_discard = 0
            out_ref_fname = os.path.join(runtime.cwd, "ref_sbref.nii.gz")

        # Set interface outputs
        self._results["n_volumes_to_discard"] = n_volumes_to_discard
        self._results["ref_image"] = out_ref_fname

        # Slicing may induce inconsistencies with shape-dependent values in extensions.
        # For now, remove all. If this turns out to be a mistake, we can select extensions
        # that don't break pipeline stages.
        ref_im.header.extensions.clear()

        # If reference is only 1 volume, return it directly
        if ref_im.dataobj.ndim == 3:
            ref_im.to_filename(out_ref_fname)
            return runtime

        if n_volumes_to_discard == 0:
            if ref_im.shape[-1] > 40:
                ref_im = nb.Nifti1Image(
                    ref_im.dataobj[:, :, :, 20:40], ref_im.affine, ref_im.header
                )

            ref_name = os.path.join(runtime.cwd, "slice.nii.gz")
            ref_im.to_filename(ref_name)
            if self.inputs.mc_method == "AFNI":
                res = afni.Volreg(
                    in_file=ref_name,
                    args="-Fourier -twopass",
                    zpad=4,
                    outputtype="NIFTI_GZ",
                ).run()
            elif self.inputs.mc_method == "FSL":
                res = fsl.MCFLIRT(
                    in_file=ref_name, ref_vol=0, interpolation="sinc"
                ).run()
            mc_slice_nii = nb.load(res.outputs.out_file)

            median_image_data = np.median(mc_slice_nii.get_fdata(), axis=3)
        else:
            median_image_data = np.median(
                ref_im.dataobj[:, :, :, :n_volumes_to_discard], axis=3
            )

        nb.Nifti1Image(median_image_data, ref_im.affine, ref_im.header).to_filename(
            out_ref_fname
        )
        return runtime


def _get_vols_to_discard(img):
    from nipype.algorithms.confounds import is_outlier

    data_slice = img.dataobj[:, :, :, :50]
    global_signal = data_slice.mean(axis=0).mean(axis=0).mean(axis=0)
    return is_outlier(global_signal)
