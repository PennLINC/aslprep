"""ANTS interfaces."""

import nibabel as nb
import numpy as np
from nipype.interfaces.base import traits
from nipype.utils.filemanip import fname_presuffix
from niworkflows.interfaces.fixes import (
    FixHeaderApplyTransforms,
    _FixTraitApplyTransformsInputSpec,
)


class _ApplyTransformsInputSpec(_FixTraitApplyTransformsInputSpec):
    # Nipype's version doesn't have GenericLabel
    interpolation = traits.Enum(
        'Linear',
        'NearestNeighbor',
        'CosineWindowedSinc',
        'WelchWindowedSinc',
        'HammingWindowedSinc',
        'LanczosWindowedSinc',
        'MultiLabel',
        'Gaussian',
        'BSpline',
        'GenericLabel',
        argstr='%s',
        usedefault=True,
    )


class ApplyTransforms(FixHeaderApplyTransforms):
    """A modified version of FixHeaderApplyTransforms from niworkflows.

    The niworkflows version of ApplyTransforms "fixes the resampled image header
    to match the xform of the reference image".
    This modification overrides the allowed interpolation values,
    since FixHeaderApplyTransforms doesn't support GenericLabel,
    which is preferred over MultiLabel.
    """

    input_spec = _ApplyTransformsInputSpec

    def _run_interface(self, runtime):
        # Run normally
        self.inputs.output_image = fname_presuffix(
            self.inputs.input_image,
            suffix='_trans.nii.gz',
            newpath=runtime.cwd,
            use_ext=False,
        )
        runtime = super()._run_interface(runtime)
        return runtime


class RegridImage(ApplyTransforms):
    """Resample an image onto the reference grid, unless it is already on that grid.

    ANTs' ``Gaussian`` interpolation smooths the image even under an identity transform,
    so an image that already has the reference's shape and affine is passed through unchanged.
    """

    def _run_interface(self, runtime):
        self._on_grid = _same_grid(self.inputs.input_image, self.inputs.reference_image)
        if self._on_grid:
            return runtime

        return super()._run_interface(runtime)

    def _list_outputs(self):
        if getattr(self, '_on_grid', False):
            return {'output_image': self.inputs.input_image}

        return super()._list_outputs()


def _same_grid(in_file, ref_file, atol=1e-4):
    """Check that two images share their spatial shape and affine.

    Parameters
    ----------
    in_file, ref_file : :obj:`str`
        Paths to the images to compare. Either may be 4D.
    atol : :obj:`float`
        Absolute tolerance on the elements of the affines.

    Returns
    -------
    :obj:`bool`
    """
    img, ref = nb.load(in_file), nb.load(ref_file)
    return img.shape[:3] == ref.shape[:3] and np.allclose(
        img.affine, ref.affine, rtol=0, atol=atol
    )
