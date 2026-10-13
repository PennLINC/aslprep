"""Tests for the aslprep.interfaces.ants submodule."""

import shutil

import nibabel as nb
import numpy as np
import pytest
from nipype.pipeline import engine as pe

from aslprep.interfaces import ants


def _write(path, shape, affine, seed=0):
    data = np.random.default_rng(seed).normal(1000, 100, shape).astype('f4')
    nb.Nifti1Image(data, affine).to_filename(path)
    return str(path)


def test_regridimage_passes_through_on_grid(tmp_path):
    """An image already on the reference grid is returned untouched, without calling ANTs."""
    affine = np.diag([3.0, 3.0, 4.0, 1.0])
    ref = _write(tmp_path / 'ref.nii.gz', (8, 8, 6), affine)

    # A 4D M0 scan on the grid of a 3D reference still counts as on-grid.
    for name, shape in (('m0_3d.nii.gz', (8, 8, 6)), ('m0_4d.nii.gz', (8, 8, 6, 2))):
        m0 = _write(tmp_path / name, shape, affine, seed=1)
        regrid = pe.Node(
            ants.RegridImage(
                input_image=m0,
                reference_image=ref,
                interpolation='Gaussian',
                transforms=['identity'],
            ),
            name=f'regrid_{name.split(".")[0]}',
            base_dir=str(tmp_path),
        )
        # ANTs is not needed (and not called) in this case.
        results = regrid.run()
        assert results.outputs.output_image == m0


@pytest.mark.parametrize(
    ('shape', 'affine'),
    [
        ((16, 16, 12), np.diag([1.5, 1.5, 2.0, 1.0])),  # finer grid
        ((8, 8, 6), np.array([[3, 0, 0, 0.5], [0, 3, 0, 0], [0, 0, 4, 0], [0, 0, 0, 1.0]])),
        ((8, 8, 5), np.diag([3.0, 3.0, 4.0, 1.0])),  # same affine, different shape
    ],
)
def test_same_grid_detects_differences(tmp_path, shape, affine):
    """A different shape, resolution or origin is not the same grid."""
    ref_affine = np.diag([3.0, 3.0, 4.0, 1.0])
    ref = _write(tmp_path / 'ref.nii.gz', (8, 8, 6), ref_affine)
    assert ants._same_grid(ref, ref)

    other = _write(tmp_path / 'other.nii.gz', shape, affine)
    assert not ants._same_grid(other, ref)


@pytest.mark.skipif(shutil.which('antsApplyTransforms') is None, reason='ANTs is not installed')
def test_regridimage_resamples_off_grid(tmp_path):
    """An image on another grid is resampled onto the reference grid."""
    ref = _write(tmp_path / 'ref.nii.gz', (8, 8, 6), np.diag([3.0, 3.0, 4.0, 1.0]))
    m0 = _write(tmp_path / 'm0.nii.gz', (16, 16, 12), np.diag([1.5, 1.5, 2.0, 1.0]), seed=1)

    regrid = ants.RegridImage(
        input_image=m0,
        reference_image=ref,
        interpolation='Gaussian',
        transforms=['identity'],
    )
    results = regrid.run(cwd=tmp_path)
    assert results.outputs.output_image != m0
    assert ants._same_grid(results.outputs.output_image, ref)
