"""Interfaces for gradient nonlinearity correction.

The displacement field is built with TORTOISE V4's ``CreateNonlinearityDisplacementMap``.
These wrappers are adapted from QSIPrep (https://github.com/PennLINC/qsiprep),
which documents the hazards in the upstream binary that shape them:

1. The coefficient file is the *first* positional argument
   (``mk_displacement(argv[1], img, is_GE)`` in ``src/tools/gradnonlin/mk_displacementMaps.cxx``).
2. A fourth argument is read with ``(bool)atoi(argv[4])``, so it is appended only for GE data.
3. The output is the field TORTOISE calls ``gradwarp_field_inv``:
   the one TORTOISE itself resamples with.
   It maps a point in the corrected image to where that point was recorded in the raw image,
   so it is used as-is and must **not** be inverted.
"""

import os

import nibabel as nb
import numpy as np
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    CommandLine,
    CommandLineInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    isdefined,
    traits,
)


class _CreateNonlinearityDisplacementMapInputSpec(CommandLineInputSpec):
    coeff_file = File(
        exists=True,
        mandatory=True,
        argstr='%s',
        position=0,
        desc='Scanner gradient coefficient file (.grad, .dat, or .gc)',
    )
    ref_image = File(
        exists=True,
        mandatory=True,
        argstr='%s',
        position=1,
        desc='3D image defining the grid the field is generated on',
    )
    out_field = traits.Str(
        'gradwarp_field.nii',
        usedefault=True,
        argstr='%s',
        position=2,
        desc='Output displacement field name. Must end in .nii.',
    )
    # No argstr: appended in _parse_inputs only when True. See module docstring.
    is_ge = traits.Bool(False, usedefault=True, desc='Whether the scanner is a GE scanner')


class _CreateNonlinearityDisplacementMapOutputSpec(TraitedSpec):
    out_field = File(exists=True, desc='ITK displacement field (LPS, mm) on the reference grid')


class CreateNonlinearityDisplacementMap(CommandLine):
    """Expand gradient coefficients into a displacement field with TORTOISE.

    The output is TORTOISE's ``gradwarp_field_inv``. Do not invert it.
    """

    input_spec = _CreateNonlinearityDisplacementMapInputSpec
    output_spec = _CreateNonlinearityDisplacementMapOutputSpec
    _cmd = 'CreateNonlinearityDisplacementMap'

    def _parse_inputs(self, skip=None):
        parsed = super()._parse_inputs(skip=skip)
        if self.inputs.is_ge:
            parsed.append('1')
        return parsed

    def _list_outputs(self):
        return {'out_field': os.path.abspath(self.inputs.out_field)}


class _MaskWarpDimensionsInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='ITK displacement field')
    ref_image = File(
        exists=True,
        desc=(
            'Image whose slice axis defines the through-plane direction. '
            'Required unless warp_dim is "3D".'
        ),
    )
    slice_axis = traits.Enum(
        'k',
        'i',
        'j',
        usedefault=True,
        desc='Voxel axis of ref_image along which slices are stacked',
    )
    warp_dim = traits.Enum(
        '3D',
        '2D',
        '1D',
        usedefault=True,
        desc=(
            'Which displacement components to keep. '
            '"3D" keeps all, "2D" zeroes the through-plane component, '
            'and "1D" keeps only the through-plane component.'
        ),
    )


class _MaskWarpDimensionsOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='Displacement field with components zeroed')


class MaskWarpDimensions(SimpleInterface):
    """Remove the displacement components a scanner has already corrected.

    The through-plane direction is the physical slice normal of ``ref_image``.
    TORTOISE zeroes world (LPS) components instead,
    which is only equivalent for axial acquisitions.
    """

    input_spec = _MaskWarpDimensionsInputSpec
    output_spec = _MaskWarpDimensionsOutputSpec

    def _run_interface(self, runtime):
        img = nb.load(self.inputs.in_file)
        data = img.get_fdata(dtype='float32')
        if self.inputs.warp_dim != '3D':
            if not isdefined(self.inputs.ref_image):
                raise ValueError(f'ref_image is required for warp_dim={self.inputs.warp_dim}.')
            normal = slice_normal_lps(
                nb.load(self.inputs.ref_image).affine, self.inputs.slice_axis
            )
            # ITK vector fields are (X, Y, Z, 1, 3): the last axis holds the LPS components.
            through_plane = (data @ normal)[..., np.newaxis] * normal
            data = through_plane if self.inputs.warp_dim == '1D' else data - through_plane

        out_img = nb.Nifti1Image(data, img.affine, img.header)
        # Without the vector intent, ANTs reads a 5D image as an all-zero field.
        out_img.header.set_intent('vector')
        out_file = os.path.join(runtime.cwd, 'gradwarp_field_masked.nii.gz')
        out_img.to_filename(out_file)
        self._results['out_file'] = out_file
        return runtime


def slice_normal_lps(affine, slice_axis='k'):
    """Return the unit slice normal of an image, in LPS (ITK) world coordinates."""
    column = np.asarray(affine, dtype='float64')[:3, 'ijk'.index(slice_axis)]
    return (column / np.linalg.norm(column) * [-1.0, -1.0, 1.0]).astype('float32')
