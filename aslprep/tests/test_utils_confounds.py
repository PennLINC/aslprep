"""Tests for the aslprep.utils.confounds module."""

import numpy as np
import pytest

from aslprep.utils.confounds import dispersion_index


def _make_masks(shape, gm_n=2, wm_n=2, csf_n=2):
    """Create simple non-overlapping boolean masks with the requested voxel counts."""
    gm_mask = np.zeros(shape, dtype=bool)
    wm_mask = np.zeros(shape, dtype=bool)
    csf_mask = np.zeros(shape, dtype=bool)

    flat_gm = gm_mask.ravel()
    flat_wm = wm_mask.ravel()
    flat_csf = csf_mask.ravel()

    flat_gm[:gm_n] = True
    flat_wm[gm_n : gm_n + wm_n] = True
    flat_csf[gm_n + wm_n : gm_n + wm_n + csf_n] = True

    return (
        gm_mask.reshape(shape),
        wm_mask.reshape(shape),
        csf_mask.reshape(shape),
    )


def test_dispersion_index_normal():
    """dispersion_index returns a finite, non-negative float for valid inputs."""
    rng = np.random.default_rng(0)
    shape = (5, 5, 5)
    cbf_image = rng.random(shape).astype(np.float32) + 0.1  # ensure positive, non-zero values

    gm_mask, wm_mask, csf_mask = _make_masks(shape)

    result = dispersion_index(wm_mask, gm_mask, csf_mask, cbf_image)
    assert np.isfinite(result)
    assert result >= 0


@pytest.mark.parametrize(
    ('desc', 'gm_n', 'wm_n', 'csf_n', 'match'),
    [
        ('GM has 1 voxel', 1, 2, 2, 'GM mask contains 1 voxel'),
        ('WM has 1 voxel', 2, 1, 2, 'WM mask contains 1 voxel'),
        ('CSF has 1 voxel', 2, 2, 1, 'CSF mask contains 1 voxel'),
        ('GM has 0 voxels', 0, 2, 2, 'GM mask contains 0 voxel'),
    ],
)
def test_dispersion_index_insufficient_tissue_voxels(desc, gm_n, wm_n, csf_n, match):
    """dispersion_index raises ValueError when any tissue mask has fewer than 2 voxels."""
    rng = np.random.default_rng(0)
    shape = (5, 5, 5)
    cbf_image = rng.random(shape).astype(np.float32)

    gm_mask, wm_mask, csf_mask = _make_masks(shape, gm_n=gm_n, wm_n=wm_n, csf_n=csf_n)

    with pytest.raises(ValueError, match=match):
        dispersion_index(wm_mask, gm_mask, csf_mask, cbf_image)


def test_dispersion_index_zero_mean_gm_cbf():
    """dispersion_index raises ValueError when mean GM CBF is zero."""
    shape = (5, 5, 5)
    cbf_image = np.zeros(shape, dtype=np.float32)

    gm_mask, wm_mask, csf_mask = _make_masks(shape)

    with pytest.raises(ValueError, match='Mean GM CBF is zero'):
        dispersion_index(wm_mask, gm_mask, csf_mask, cbf_image)
