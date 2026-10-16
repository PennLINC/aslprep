"""Tests for aslprep.utils.pvc (design spec Section 7.1)."""

import numpy as np
import pytest
from scipy import ndimage

from aslprep.utils import pvc


def _rotation(rx, ry, rz):
    cx, sx = np.cos(rx), np.sin(rx)
    cy, sy = np.cos(ry), np.sin(ry)
    cz, sz = np.cos(rz), np.sin(rz)
    r_x = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    r_y = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    r_z = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return r_z @ r_y @ r_x


def _sphere_tpm(shape=(48, 48, 48), centre=(24.3, 23.7, 24.1), radius=12.0, width=2.0):
    """A smooth 'tissue' probability: 1 inside a sphere, falling to 0 over ``width`` mm."""
    grid = np.indices(shape).astype(float)
    dist = np.sqrt(sum((g - c) ** 2 for g, c in zip(grid, centre, strict=True)))
    return np.clip((radius + width / 2 - dist) / width, 0, 1)


# --- Effective resolution (test 7) ---


@pytest.mark.parametrize(
    ('metadata', 'source'),
    [
        ({'MRAcquisitionType': '2D', 'PulseSequenceType': 'EPI'}, 'auto:2d'),
        ({'MRAcquisitionType': '3D', 'PulseSequenceType': '3D_GRASE'}, 'auto:grase-single-shot'),
        (
            {'MRAcquisitionType': '3D', 'PulseSequenceType': '3Dgrase', 'NumberShots': 1},
            'auto:grase-single-shot',
        ),
        (
            {'MRAcquisitionType': '3D', 'PulseSequenceType': '3Dgrase', 'NumberShots': 4},
            'nominal-fallback',
        ),
        ({'MRAcquisitionType': '3D', 'PulseSequenceType': 'spiral'}, 'auto:spiral-petr2018'),
        ({'MRAcquisitionType': '3D', 'PulseSequenceType': 'EPI'}, 'nominal-fallback'),
        ({}, 'nominal-fallback'),
    ],
)
def test_resolve_effective_resolution_auto(metadata, source):
    """Metadata cases map to the stated source labels and factors."""
    out = pvc.resolve_effective_resolution('auto', metadata, (3.0, 3.0, 4.0), slice_axis=2)
    assert out['source'] == source
    if source in ('auto:2d', 'nominal-fallback'):
        assert out['fwhm'] == (3.0, 3.0, 4.0)
    else:
        factors = pvc.EFFECTIVE_RESOLUTION_FACTORS[source]
        np.testing.assert_allclose(out['fwhm'], [3 * factors[0], 3 * factors[1], 4 * factors[2]])
    if source in ('nominal-fallback', 'auto:spiral-petr2018'):
        assert out['warnings']


def test_resolve_effective_resolution_spiral_factors():
    """The spiral factors come from Petr et al. (2018)'s published FWHM."""
    factors = pvc.EFFECTIVE_RESOLUTION_FACTORS['auto:spiral-petr2018']
    assert factors[2] == pytest.approx(11.8 / 4)
    assert factors[0] == factors[1] == pytest.approx((4.6 + 4.3) / 2 / 3.8)


def test_resolve_effective_resolution_slice_axis():
    """The slice factor follows the slice-encoding axis."""
    meta = {'MRAcquisitionType': '3D', 'PulseSequenceType': 'GRASE'}
    out = pvc.resolve_effective_resolution('auto', meta, (4.0, 3.0, 3.0), slice_axis=0)
    np.testing.assert_allclose(out['fwhm'], [4 * 1.38, 3 * 1.1, 3 * 1.1])

    out = pvc.resolve_effective_resolution([3.5, 3.6, 9.0], {}, (3.0, 4.0, 3.0), slice_axis=1)
    assert out['source'] == 'user'
    assert out['fwhm'] == (3.5, 9.0, 3.6)


def test_resolve_effective_resolution_rejects_bad_values():
    """User FWHMs below the voxel size, or non-positive, are rejected."""
    with pytest.raises(ValueError, match='finer than the voxel size'):
        pvc.resolve_effective_resolution([2.0, 3.0, 4.0], {}, (3.0, 3.0, 4.0), slice_axis=2)
    for bad in (-3.0, np.nan, np.inf):
        with pytest.raises(ValueError, match='finite and positive'):
            pvc.resolve_effective_resolution([3.0, bad, 4.0], {}, (3.0, 3.0, 4.0), slice_axis=2)
    with pytest.raises(ValueError, match='three'):
        pvc.resolve_effective_resolution([3.0, 4.0], {}, (3.0, 3.0, 4.0), slice_axis=2)


def test_resolve_effective_resolution_nominal():
    """Nominal is a box of one voxel."""
    out = pvc.resolve_effective_resolution('nominal', {}, (3.0, 3.0, 4.0), slice_axis=2)
    assert out['source'] == 'nominal'
    assert out['fwhm'] == (3.0, 3.0, 4.0)


def test_acquisition_voxel_size_warning():
    """A reconstruction that interpolated is flagged."""
    meta = {'AcquisitionVoxelSize': [4.0, 4.0, 8.0]}
    out = pvc.resolve_effective_resolution('nominal', meta, (2.0, 2.0, 4.0), slice_axis=2)
    assert any('AcquisitionVoxelSize' in w for w in out['warnings'])


@pytest.mark.parametrize('relative_fwhm', [1.1, 1.21, 1.38, 2.95])
def test_box_gaussian_fwhm(relative_fwhm):
    """The continuous box-Gaussian kernel has the requested FWHM (within 0.1%)."""
    width = 3.7
    sigma = pvc.solve_box_gaussian_sigma(relative_fwhm * width, width)
    assert sigma > 0
    fwhm = pvc.box_gaussian_fwhm(width, sigma)
    assert fwhm == pytest.approx(relative_fwhm * width, rel=1e-3)


def test_box_gaussian_sigma_box_only():
    """A FWHM of one voxel is the box alone, and narrower requests are rejected."""
    assert pvc.solve_box_gaussian_sigma(4.0, 4.0) == 0
    with pytest.raises(ValueError, match='at least the box width'):
        pvc.solve_box_gaussian_sigma(3.0, 4.0)


@pytest.mark.parametrize('relative_fwhm', [1.0, 1.1, 1.38, 2.95])
def test_discrete_kernel_matches_continuous(relative_fwhm):
    """The bin-integrated kernel matches the continuous kernel averaged over each bin."""
    width = 4.0
    subdiv = 9
    step = width / subdiv
    sigma = pvc.solve_box_gaussian_sigma(relative_fwhm * width, width)
    kernel = pvc.discrete_box_gaussian_kernel(width, sigma, step)
    assert kernel.size % 2 == 1
    np.testing.assert_allclose(kernel, kernel[::-1])
    assert kernel.sum() == pytest.approx(1.0)

    # Independent reference: fine quadrature of the continuous profile within each bin
    n_half = (kernel.size - 1) // 2
    fine = 201
    reference = []
    for j in range(-n_half, n_half + 1):
        x = j * step + (np.arange(fine) + 0.5) / fine * step - step / 2
        reference.append(pvc.box_gaussian_profile(x, width, sigma).mean())
    reference = np.array(reference)
    reference /= reference.sum()
    assert np.abs(kernel - reference).sum() < 1e-3

    # Truncation at width / 2 + 4 sigma loses almost no mass
    if sigma > 0:
        x = np.linspace(-width / 2 - 10 * sigma, width / 2 + 10 * sigma, 20001)
        prof = pvc.box_gaussian_profile(x, width, sigma)
        inside = np.abs(x) <= width / 2 + 4 * sigma
        assert prof[~inside].sum() / prof.sum() < 1e-4


def test_choose_subdivision():
    """The subdivision is the smallest odd value with a small enough subvoxel."""
    assert pvc.choose_subdivision(4.0, 0.5) == 9
    assert pvc.choose_subdivision(3.0, 1.0) == 3
    assert pvc.choose_subdivision(3.0, 0.5) == 7


# --- PV-map construction (test 8) ---


def _reference_box_average(tpm, tpm_affine, ref_shape, ref_affine, ref2anat, fine=15):
    """Average the linearly interpolated TPM over each reference voxel by fine quadrature."""
    offsets = (np.arange(fine) - (fine - 1) / 2) / fine
    ref2vox = np.linalg.inv(tpm_affine) @ ref2anat @ ref_affine
    out = np.zeros(ref_shape)
    sub = np.stack(np.meshgrid(offsets, offsets, offsets, indexing='ij'), -1).reshape(-1, 3)
    for idx in np.ndindex(*ref_shape):
        pts = np.asarray(idx)[None, :] + sub
        pts_h = np.hstack([pts, np.ones((pts.shape[0], 1))])
        vox = (ref2vox @ pts_h.T)[:3]
        out[idx] = ndimage.map_coordinates(tpm, vox, order=1, mode='constant', cval=0.0).mean()
    return out


@pytest.mark.parametrize(
    'ref2anat',
    [
        np.eye(4),
        np.array([[1, 0, 0, 0.37], [0, 1, 0, -0.81], [0, 0, 1, 0.52], [0, 0, 0, 1.0]]),
        np.block(
            [
                [_rotation(0.12, -0.08, 0.2), np.array([[0.4], [-0.3], [0.7]])],
                [np.zeros((1, 3)), np.ones((1, 1))],
            ]
        ),
    ],
    ids=['identity', 'translated', 'oblique'],
)
def test_build_pv_maps_nominal_matches_reference(ref2anat):
    """Nominal PV maps equal the voxel average of the linearly interpolated TPM."""
    tpm = _sphere_tpm()
    tpm_affine = np.eye(4)
    ref_affine = np.diag([4.0, 4.0, 4.0, 1.0])
    ref_affine[:3, 3] = [2.0, 1.5, 2.5]
    ref_shape = (11, 11, 11)

    pv_maps, info = pvc.build_pv_maps(
        [tpm], tpm_affine, ref_shape, ref_affine, ref2anat, fwhm=(4.0, 4.0, 4.0)
    )
    reference = _reference_box_average(tpm, tpm_affine, ref_shape, ref_affine, ref2anat)
    assert np.abs(pv_maps[0] - reference).max() < 0.01
    assert info['sigma_mm'] == [0.0, 0.0, 0.0]


def test_build_pv_maps_convergence():
    """Tripling the subdivision changes the PV maps by less than 0.01."""
    tpm = _sphere_tpm()
    ref_affine = np.diag([4.0, 4.0, 5.0, 1.0])
    ref2anat = np.block(
        [
            [_rotation(0.05, 0.1, -0.15), np.array([[0.2], [0.6], [-0.4]])],
            [np.zeros((1, 3)), np.ones((1, 1))],
        ]
    )
    kwargs = {
        'tpms': [tpm],
        'tpm_affine': np.eye(4),
        'ref_shape': (11, 11, 9),
        'ref_affine': ref_affine,
        'ref2anat': ref2anat,
        'fwhm': (4.4, 4.4, 6.9),
    }
    coarse, _ = pvc.build_pv_maps(max_subvoxel=1.0, **kwargs)
    fine, _ = pvc.build_pv_maps(max_subvoxel=1.0 / 3, **kwargs)
    assert np.abs(coarse[0] - fine[0]).max() < 0.01


@pytest.mark.parametrize('slice_axis', [0, 1, 2])
def test_build_pv_maps_slabs_match_whole_volume(slice_axis):
    """Slab processing gives the same result as processing the whole volume at once."""
    tpm = _sphere_tpm()
    kwargs = {
        'tpms': [tpm, 1 - tpm],
        'tpm_affine': np.eye(4),
        'ref_shape': (11, 10, 9),
        'ref_affine': np.diag([4.0, 4.0, 4.0, 1.0]),
        'ref2anat': np.eye(4),
        'fwhm': (4.4, 5.0, 8.0),
        'max_subvoxel': 1.0,
        'slice_axis': slice_axis,
    }
    whole, _ = pvc.build_pv_maps(max_points=10**9, **kwargs)
    slabbed, _ = pvc.build_pv_maps(max_points=20_000, **kwargs)
    for a, b in zip(whole, slabbed, strict=True):
        np.testing.assert_allclose(a, b, atol=1e-6)


def test_build_pv_maps_zero_outside_fov():
    """Tissue outside the anatomical field of view counts as zero."""
    tpm = np.ones((20, 20, 20))
    pv_maps, _ = pvc.build_pv_maps(
        [tpm],
        np.eye(4),
        (6, 6, 6),
        np.diag([4.0, 4.0, 4.0, 1.0]),
        np.eye(4),
        fwhm=(4.0, 4.0, 4.0),
        max_subvoxel=0.5,
    )
    # Voxel 5 spans 18-22 mm along each axis; the field is linearly interpolated to zero
    # beyond index 19.
    assert pv_maps[0][2, 2, 2] == pytest.approx(1.0)
    assert pv_maps[0][5, 2, 2] < 0.5


# --- Local regression (tests 1-5) ---


def _two_tissue_phantom(shape=(20, 20, 3), seed=0):
    """Smoothly varying GM/WM fractions with constant tissue amplitudes."""
    rng = np.random.default_rng(seed)
    x = np.linspace(0, 1, shape[0])[:, None, None]
    y = np.linspace(0, 1, shape[1])[None, :, None]
    pv_gm = 0.2 + 0.6 * (0.5 + 0.5 * np.sin(6 * x + 3 * y)) * np.ones(shape)
    pv_wm = (1 - pv_gm) * (0.6 + 0.3 * rng.random(shape))
    return pv_gm, pv_wm


def test_lr_exact_recovery():
    """Noiseless data under the model are recovered exactly in two-tissue voxels."""
    pv_gm, pv_wm = _two_tissue_phantom()
    y = 54.0 * pv_gm + 16.0 * pv_wm
    mask = np.ones(y.shape, bool)
    fit = pvc.local_regression_pvc(y, pv_gm, pv_wm, mask, mask, kernel_size=(5, 5, 1))
    two = fit['status'] == 1
    assert two.mean() > 0.5
    np.testing.assert_allclose(fit['beta_gm'][two], 54.0, rtol=5e-3)
    np.testing.assert_allclose(fit['beta_wm'][two], 16.0, rtol=5e-3)


def test_lr_through_plane_kernel_3d():
    """Without slice timing, a through-plane kernel recovers constant amplitudes too."""
    pv_gm, pv_wm = _two_tissue_phantom(shape=(16, 16, 8))
    y = 54.0 * pv_gm + 16.0 * pv_wm
    mask = np.ones(y.shape, bool)
    fit = pvc.local_regression_pvc(y, pv_gm, pv_wm, mask, mask, kernel_size=(5, 5, 3))
    two = fit['status'] == 1
    np.testing.assert_allclose(fit['beta_gm'][two], 54.0, rtol=5e-3)


def test_lr_psf_mismatch_biases_gm():
    """Regressing on PV maps that ignore the readout blur biases GM by more than 2%."""
    rng = np.random.default_rng(3)
    hi_gm = ndimage.gaussian_filter((rng.random((80, 80, 8)) > 0.6).astype(float), 1.0)
    hi_wm = 1 - hi_gm
    blur = (2.0, 2.0, 0)
    pv_gm_true = ndimage.gaussian_filter(hi_gm, blur)[::4, ::4]
    pv_wm_true = ndimage.gaussian_filter(hi_wm, blur)[::4, ::4]
    pv_gm_nominal = hi_gm[::4, ::4]
    pv_wm_nominal = hi_wm[::4, ::4]
    y = 60 * pv_gm_true + 20 * pv_wm_true
    mask = np.ones(y.shape, bool)

    exact = pvc.local_regression_pvc(y, pv_gm_true, pv_wm_true, mask, mask)
    nominal = pvc.local_regression_pvc(y, pv_gm_nominal, pv_wm_nominal, mask, mask)
    both = (exact['status'] == 1) & (nominal['status'] == 1)
    assert both.sum() > 20
    np.testing.assert_allclose(exact['beta_gm'][both], 60, rtol=5e-3)
    bias = np.median(nominal['beta_gm'][both] / 60 - 1)
    assert abs(bias) > 0.02


def test_lr_identifiability():
    """Collinear fractions are not identifiable; pure tissue falls back to one tissue."""
    shape = (9, 9, 1)
    mask = np.ones(shape, bool)

    # 50:50 everywhere: the GM-only fit would report 80 for GM 60 and WM 20.
    half = np.full(shape, 0.5)
    fit = pvc.local_regression_pvc(60 * half + 20 * half, half, half, mask, mask)
    assert np.all(fit['status'] == 4)
    assert np.all(fit['beta_gm'] == 0)

    # Pure GM everywhere
    fit = pvc.local_regression_pvc(
        np.full(shape, 60.0), np.ones(shape), np.zeros(shape), mask, mask
    )
    assert np.all(fit['status'] == 2)
    np.testing.assert_allclose(fit['beta_gm'], 60)

    # No tissue support
    fit = pvc.local_regression_pvc(np.zeros(shape), np.zeros(shape), np.zeros(shape), mask, mask)
    assert np.all(fit['status'] == 4)

    # Outside the mask
    out_mask = mask.copy()
    out_mask[0, 0, 0] = False
    fit = pvc.local_regression_pvc(
        np.full(shape, 60.0), np.ones(shape), np.zeros(shape), out_mask, out_mask
    )
    assert fit['status'][0, 0, 0] == 0


def test_lr_invalid_centre():
    """An invalid centre in a well-conditioned neighbourhood is not fitted."""
    pv_gm, pv_wm = _two_tissue_phantom()
    y = 54.0 * pv_gm + 16.0 * pv_wm
    mask = np.ones(y.shape, bool)
    valid = mask.copy()
    valid[10, 10, 1] = False
    fit = pvc.local_regression_pvc(y, pv_gm, pv_wm, mask, valid)
    assert fit['status'][10, 10, 1] == 4
    assert fit['beta_gm'][10, 10, 1] == 0
    # Its neighbours are still fitted, without it.
    assert fit['status'][11, 10, 1] == 1


def test_lr_matches_brute_force():
    """The box-filter implementation matches a per-voxel least-squares fit."""
    rng = np.random.default_rng(42)
    shape = (12, 11, 3)
    pv_gm = rng.random(shape)
    pv_wm = rng.random(shape) * (1 - pv_gm)
    y = 1e6 * (0.5 + rng.random(shape))
    y[rng.random(shape) < 0.1] = 0.0  # genuine zeros
    y[2, 3, 1] = np.nan
    mask = rng.random(shape) > 0.1
    valid = mask & (rng.random(shape) > 0.05)
    kernel = (5, 3, 1)

    fit = pvc.local_regression_pvc(y, pv_gm, pv_wm, mask, valid, kernel_size=kernel, kappa_max=1e4)
    ok = valid & mask & np.isfinite(y) & ((pv_gm + pv_wm) > 0)
    half = [k // 2 for k in kernel]
    checked = 0
    for idx in zip(*np.nonzero(fit['status'] == 1), strict=True):
        sl = tuple(
            slice(max(i - h, 0), min(i + h + 1, n))
            for i, h, n in zip(idx, half, shape, strict=True)
        )
        sel = ok[sl]
        design = np.column_stack([pv_gm[sl][sel], pv_wm[sl][sel]])
        target = y[sl][sel]
        beta, *_ = np.linalg.lstsq(design, target, rcond=None)
        rss = np.sum((target - design @ beta) ** 2)
        np.testing.assert_allclose([fit['beta_gm'][idx], fit['beta_wm'][idx]], beta, rtol=1e-10)
        np.testing.assert_allclose(
            fit['rmse'][idx] ** 2 * (sel.sum() - 2), rss, rtol=1e-8, atol=1e-10 * np.sum(target**2)
        )
        checked += 1
    assert checked > 50


def test_lr_keeps_small_residuals():
    """A near-exact fit keeps its small, positive residual."""
    shape = (5, 5, 1)
    pv_gm = np.ones(shape)
    pv_wm = np.zeros(shape)
    y = 60 + 0.001 * np.where(np.indices(shape).sum(0) % 2, 1, -1)
    mask = np.ones(shape, bool)
    fit = pvc.local_regression_pvc(y, pv_gm, pv_wm, mask, mask)
    centre = (2, 2, 0)
    assert fit['status'][centre] == 2
    expected = np.sum((y - y.mean()) ** 2)
    assert fit['rmse'][centre] ** 2 * (25 - 1) == pytest.approx(expected, rel=1e-6)
    assert fit['rmse'][centre] > 0


def test_lr_noise_amplification():
    """The empirical spread of the GM estimate matches sigma * noise amplification."""
    pv_gm, pv_wm = _two_tissue_phantom(shape=(12, 12, 1))
    mask = np.ones(pv_gm.shape, bool)
    rng = np.random.default_rng(7)
    sigma = 2.0
    estimates = []
    for _ in range(500):
        y = 54 * pv_gm + 16 * pv_wm + rng.normal(0, sigma, pv_gm.shape)
        fit = pvc.local_regression_pvc(y, pv_gm, pv_wm, mask, mask)
        estimates.append(fit['beta_gm'])
    estimates = np.array(estimates)
    centre = (6, 6, 0)
    assert fit['status'][centre] == 1
    assert estimates[:, 6, 6, 0].std() == pytest.approx(sigma * fit['noise_gm'][centre], rel=0.1)


def test_lr_kernel_validation():
    """Even or non-positive kernel sizes are rejected."""
    arr = np.zeros((5, 5, 1))
    with pytest.raises(ValueError, match='odd'):
        pvc.local_regression_pvc(arr, arr, arr, arr > -1, arr > -1, kernel_size=(4, 5, 1))


def test_calibrate_single_tissue():
    """A single-tissue kernel gives the uncorrected CBF, y / M0, and failures are excluded."""
    shape = (7, 7, 1)
    m0 = np.full(shape, 0.9)
    y = np.full(shape, 60.0)
    mask = np.ones(shape, bool)
    fit = pvc.local_regression_pvc(y, np.ones(shape), np.zeros(shape), mask, mask)
    out = pvc.calibrate_pvc(fit, m0, np.ones(shape), np.zeros(shape))
    np.testing.assert_allclose(out['cbf_gm'], 60 / 0.9)
    assert out['valid_gm'].all()
    assert not out['valid_wm'].any()
    assert np.all(out['cbf_wm'] == 0)


def test_kernel_size_to_axes_and_slice_timing():
    """Kernel sizes map onto array axes, and slice timing variation is detected."""
    assert pvc.kernel_size_to_axes((5, 5, 1), 2) == (5, 5, 1)
    assert pvc.kernel_size_to_axes((5, 3, 1), 0) == (1, 5, 3)
    assert pvc.get_slice_axis({'SliceEncodingDirection': 'j-'}) == 1
    assert pvc.get_slice_axis({}) == 2
    assert pvc.has_varying_slice_timing({'SliceTiming': [0, 0.05, 0.1]})
    assert not pvc.has_varying_slice_timing({'SliceTiming': [0, 0, 0]})
    assert not pvc.has_varying_slice_timing({})


def _three_tissue_phantom(shape=(20, 20, 3), seed=0):
    """Random but well-mixed GM, WM and CSF fractions that sum to one."""
    rng = np.random.default_rng(seed)
    raw = rng.random((3, *shape)) + 0.2
    raw /= raw.sum(axis=0)
    return raw[0], raw[1], raw[2]


@pytest.mark.parametrize('csf_amplitude', [0.0, 12.0])
def test_lr_csf_exact_recovery(csf_amplitude):
    """The three-tissue fit recovers the CSF amplitude, whether or not it is zero."""
    pv_gm, pv_wm, pv_csf = _three_tissue_phantom()
    y = 54.0 * pv_gm + 16.0 * pv_wm + csf_amplitude * pv_csf
    mask = np.ones(y.shape, bool)
    fit = pvc.local_regression_csf(y, pv_gm, pv_wm, pv_csf, mask, mask)
    full = fit['status'] == 1
    assert full.mean() > 0.5
    np.testing.assert_allclose(fit['beta_csf'][full], csf_amplitude, atol=1e-6 * 54)
    np.testing.assert_allclose(fit['rmse'][full], 0, atol=1e-8 * 54)


def test_lr_csf_matches_brute_force():
    """The three-tissue fit matches a per-voxel least-squares fit."""
    rng = np.random.default_rng(5)
    pv_gm, pv_wm, pv_csf = _three_tissue_phantom(shape=(12, 11, 2), seed=5)
    y = 50 * pv_gm + 20 * pv_wm + rng.normal(0, 3, pv_gm.shape)
    y[3, 3, 0] = np.nan
    mask = rng.random(y.shape) > 0.05
    valid = mask & (rng.random(y.shape) > 0.05)
    kernel = (5, 5, 1)
    fit = pvc.local_regression_csf(y, pv_gm, pv_wm, pv_csf, mask, valid, kernel_size=kernel)
    ok = valid & np.isfinite(y)
    checked = 0
    for idx in zip(*np.nonzero(fit['status'] == 1), strict=True):
        sl = tuple(
            slice(max(i - 2 * (k > 1), 0), min(i + 2 * (k > 1) + 1, n))
            for i, k, n in zip(idx, kernel, y.shape, strict=True)
        )
        sel = ok[sl]
        design = np.column_stack([pv_gm[sl][sel], pv_wm[sl][sel], pv_csf[sl][sel]])
        beta, *_ = np.linalg.lstsq(design, y[sl][sel], rcond=None)
        rss = np.sum((y[sl][sel] - design @ beta) ** 2)
        assert fit['beta_csf'][idx] == pytest.approx(beta[2], rel=1e-8, abs=1e-8)
        assert fit['rmse'][idx] ** 2 * (sel.sum() - 3) == pytest.approx(rss, rel=1e-8)
        g_inv = np.linalg.inv(design.T @ design)
        assert fit['noise_csf'][idx] == pytest.approx(np.sqrt(g_inv[2, 2]), rel=1e-8)
        checked += 1
    assert checked > 50


def test_lr_csf_reduced_fits_and_identifiability():
    """Absent tissues are dropped only when they cannot leak into the CSF estimate."""
    shape = (9, 9, 1)
    mask = np.ones(shape, bool)
    zeros = np.zeros(shape)

    # Pure CSF: a CSF-only fit
    fit = pvc.local_regression_csf(np.full(shape, 5.0), zeros, zeros, np.ones(shape), mask, mask)
    assert np.all(fit['status'] == 5)
    np.testing.assert_allclose(fit['beta_csf'], 5.0)

    # GM and CSF mixed, no WM: a GM + CSF fit
    rng = np.random.default_rng(2)
    pv_gm = rng.uniform(0.2, 0.8, shape)
    pv_csf = 1 - pv_gm
    fit = pvc.local_regression_csf(60 * pv_gm, pv_gm, zeros, pv_csf, mask, mask)
    assert np.all(fit['status'] == 2)
    np.testing.assert_allclose(fit['beta_csf'], 0, atol=1e-8)

    # CSF everywhere at a constant 50:50 with GM: not identifiable, and no leaked value.
    half = np.full(shape, 0.5)
    fit = pvc.local_regression_csf(30 + 0 * half, half, zeros, half, mask, mask)
    assert np.all(fit['status'] == 4)
    assert np.all(fit['beta_csf'] == 0)

    # No CSF in the kernel: not identifiable
    fit = pvc.local_regression_csf(60 * np.ones(shape), np.ones(shape), zeros, zeros, mask, mask)
    assert np.all(fit['status'] == 4)


def test_calibrate_csf_and_summary():
    """CSF CBF is calibrated by the centre M0 and summarized where CSF dominates."""
    shape = (7, 7, 1)
    mask = np.ones(shape, bool)
    zeros = np.zeros(shape)
    fit = pvc.local_regression_csf(np.full(shape, 4.0), zeros, zeros, np.ones(shape), mask, mask)
    out = pvc.calibrate_csf(fit, np.full(shape, 0.8), np.ones(shape))
    np.testing.assert_allclose(out['cbf_csf'], 5.0)
    assert out['valid_csf'].all()
    summary = pvc.summarize_csf_qc(out['cbf_csf'], out['valid_csf'], np.ones(shape))
    assert summary == {'median': 5.0, 'iqr': 0.0, 'n': 49}
    empty = pvc.summarize_csf_qc(out['cbf_csf'], out['valid_csf'], np.zeros(shape))
    assert empty['n'] == 0
    assert np.isnan(empty['median'])


def test_normalized_resample_ratio():
    """Normalized resampling keeps a constant next to invalid zeros."""
    # A 1D 'warp' half-way between voxels: linear weights 0.5/0.5, then 0.75/0.25.
    value = np.array([60.0, 60.0, 0.0])
    valid = np.array([1.0, 1.0, 0.0])
    weights = np.array([[0.5, 0.5, 0], [0, 0.75, 0.25], [0, 0.25, 0.75]])
    warped_data = weights @ (value * valid)
    warped_valid = weights @ valid
    out = pvc.normalized_resample_ratio(warped_data, warped_valid)
    np.testing.assert_allclose(out, [60, 60, 0])
    # Plain interpolation would have given 45 for the second sample.
    assert (weights @ value)[1] == pytest.approx(45)
