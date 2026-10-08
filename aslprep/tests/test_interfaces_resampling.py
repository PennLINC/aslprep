"""Tests for aslprep.interfaces.resampling."""

import nibabel as nb
import numpy as np
import pytest
from fmriprep.interfaces.resampling import ResampleSeries as FMRIPrepResampleSeries

from aslprep.interfaces.resampling import GradwarpField, ResampleSeries

SHAPE = (12, 14, 10)
ZOOMS = (2.0, 2.5, 3.0)


def _affine():
    affine = np.diag(ZOOMS + (1.0,))
    affine[:3, 3] = [-11.0, -16.0, -13.0]
    return affine


def _write_series(path, nvols=3, seed=0):
    rng = np.random.default_rng(seed)
    grid = np.indices(SHAPE).astype('float32')
    base = np.sin(grid[0] / 3) + np.cos(grid[1] / 4) + 0.5 * np.sin(grid[2] / 2) + 3
    data = np.stack([base + 0.1 * rng.standard_normal(SHAPE) for _ in range(nvols)], axis=-1)
    nb.Nifti1Image(data.astype('float32'), _affine()).to_filename(path)
    return path


def _write_itk_field(path, deltas_lps):
    """Write an (X, Y, Z, 1, 3) ITK displacement field from (X, Y, Z, 3) LPS deltas."""
    data = deltas_lps.reshape(SHAPE + (1, 3)).astype('float32')
    img = nb.Nifti1Image(data, _affine())
    img.header.set_intent('vector')
    img.to_filename(path)
    return path


def _run(interface, tmp_path, name):
    workdir = tmp_path / name
    workdir.mkdir()
    result = interface.run(cwd=str(workdir))
    return np.asarray(nb.load(result.outputs.out_file).dataobj)


@pytest.mark.parametrize('jacobian', [True, False])
def test_zero_field_matches_fmriprep(tmp_path, jacobian):
    """A zero field must reproduce fMRIPrep's ResampleSeries exactly."""
    series = _write_series(tmp_path / 'series.nii.gz')
    field = _write_itk_field(tmp_path / 'field.nii.gz', np.zeros(SHAPE + (3,)))

    expected = _run(
        FMRIPrepResampleSeries(in_file=series, ref_file=series, jacobian=jacobian),
        tmp_path,
        'fmriprep',
    )
    result = _run(
        ResampleSeries(
            in_file=series,
            ref_file=series,
            jacobian=jacobian,
            gradwarp_field=field,
        ),
        tmp_path,
        'aslprep',
    )
    np.testing.assert_allclose(result, expected, rtol=1e-5, atol=1e-5)


def test_no_field_is_fmriprep(tmp_path):
    """Without a field, the subclass defers to fMRIPrep."""
    series = _write_series(tmp_path / 'series.nii.gz')
    expected = _run(
        FMRIPrepResampleSeries(in_file=series, ref_file=series, jacobian=False),
        tmp_path,
        'fmriprep',
    )
    result = _run(
        ResampleSeries(in_file=series, ref_file=series, jacobian=False),
        tmp_path,
        'aslprep',
    )
    np.testing.assert_array_equal(result, expected)


def test_constant_field_shifts(tmp_path):
    """A constant one-voxel displacement in z samples the next slice."""
    series = _write_series(tmp_path / 'series.nii.gz')
    deltas = np.zeros(SHAPE + (3,))
    deltas[..., 2] = ZOOMS[2]  # LPS z == RAS z
    field = _write_itk_field(tmp_path / 'field.nii.gz', deltas)

    result = _run(
        ResampleSeries(in_file=series, ref_file=series, jacobian=False, gradwarp_field=field),
        tmp_path,
        'aslprep',
    )
    source = np.asarray(nb.load(series).dataobj)
    np.testing.assert_allclose(result[:, :, :-1], source[:, :, 1:], rtol=1e-4, atol=1e-4)
    # Beyond the edge of the source, the default mode fills with zeros
    assert np.allclose(result[:, :, -1], 0)


def test_lps_x_is_flipped(tmp_path):
    """A positive ITK (LPS) x displacement moves toward the right (negative RAS x)."""
    deltas = np.zeros(SHAPE + (3,))
    deltas[..., 0] = 1.5
    field = GradwarpField(nb.load(_write_itk_field(tmp_path / 'field.nii.gz', deltas)))
    points = np.array([[0.0, 0.0, 0.0], [3.0, -2.0, 1.0]])
    np.testing.assert_allclose(field.map(points), points + [-1.5, 0, 0])


def test_jacobian_of_linear_field(tmp_path):
    """u(x) = a * x along RAS x has a Jacobian determinant of 1 + a everywhere."""
    a = 0.05
    affine = _affine()
    ijk = np.indices(SHAPE).reshape(3, -1)
    xyz = (affine[:3, :3] @ ijk + affine[:3, 3:4]).T.reshape(SHAPE + (3,))
    deltas = np.zeros(SHAPE + (3,))
    deltas[..., 0] = -a * xyz[..., 0]  # RAS x displacement of a * x, stored in LPS
    field = GradwarpField(nb.load(_write_itk_field(tmp_path / 'field.nii.gz', deltas)))

    np.testing.assert_allclose(field.jacobian(), 1 + a, rtol=1e-5)


def test_jacobian_modulation(tmp_path):
    """Intensities are scaled by the field's Jacobian unless disabled."""
    a = 0.05
    affine = _affine()
    ijk = np.indices(SHAPE).reshape(3, -1)
    xyz = (affine[:3, :3] @ ijk + affine[:3, 3:4]).T.reshape(SHAPE + (3,))
    deltas = np.zeros(SHAPE + (3,))
    deltas[..., 0] = -a * xyz[..., 0]
    field = _write_itk_field(tmp_path / 'field.nii.gz', deltas)
    series = _write_series(tmp_path / 'series.nii.gz')

    kwargs = {'in_file': series, 'ref_file': series, 'jacobian': False, 'gradwarp_field': field}
    modulated = _run(ResampleSeries(**kwargs), tmp_path, 'on')
    unmodulated = _run(ResampleSeries(gradwarp_jacobian=False, **kwargs), tmp_path, 'off')
    np.testing.assert_allclose(modulated, unmodulated * (1 + a), rtol=1e-4, atol=1e-5)


def test_bad_field_shape(tmp_path):
    img = nb.Nifti1Image(np.zeros(SHAPE, dtype='float32'), _affine())
    with pytest.raises(ValueError, match='ITK displacement field'):
        GradwarpField(img)


@pytest.mark.parametrize('gradwarp_jacobian', [False, True])
def test_field_is_evaluated_after_head_motion(tmp_path, gradwarp_jacobian):
    """The scanner-fixed field is evaluated in each volume's (moved) frame: G(H(x)).

    The image is linear in world x, the field displaces by a * x along x, and the second
    volume is translated by t along x. Sampling world point H(x) + u(H(x)) gives
    (1 + a) * (x + t), whereas evaluating the field before motion would give (1 + a) * x + t.
    """
    import nitransforms as nt

    from aslprep.interfaces.resampling import resample_image

    a, t = 0.05, 4.0
    affine = _affine()
    ijk = np.indices(SHAPE).reshape(3, -1)
    xyz = (affine[:3, :3] @ ijk + affine[:3, 3:4]).T.reshape(SHAPE + (3,))
    data = np.stack([xyz[..., 0]] * 2, axis=-1).astype('float32')
    source = nb.Nifti1Image(data, affine)

    deltas = np.zeros(SHAPE + (3,))
    deltas[..., 0] = -a * xyz[..., 0]  # RAS x displacement of a * x, stored in LPS
    field = GradwarpField(nb.load(_write_itk_field(tmp_path / 'field.nii.gz', deltas)))

    shift = np.eye(4)
    shift[0, 3] = t
    hmc = nt.linear.LinearTransformsMapping([np.eye(4), shift])

    resampled = np.asarray(
        resample_image(
            source=source,
            target=source,
            transforms=nt.TransformChain([hmc]),
            fieldmap=None,
            pe_info=None,
            gradwarp=field,
            jacobian=False,
            gradwarp_jacobian=gradwarp_jacobian,
            # Linear interpolation is exact for linear data
            order=1,
        ).dataobj
    )

    scale = (1 + a) if gradwarp_jacobian else 1.0
    x = xyz[..., 0]
    # Stay away from the edges, where sampling falls outside the source
    inner = (slice(4, -4), slice(1, -1), slice(1, -1))
    np.testing.assert_allclose(resampled[..., 0][inner], scale * (1 + a) * x[inner], atol=1e-3)
    np.testing.assert_allclose(
        resampled[..., 1][inner], scale * (1 + a) * (x[inner] + t), atol=1e-3
    )
