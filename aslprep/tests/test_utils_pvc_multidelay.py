"""Tests for multi-delay local-regression PVC (design spec Section 5.7.7)."""

import numpy as np
import pytest

from aslprep.utils import pvc

KINETICS = {'labeleff': 0.85, 't1blood': 1.65, 'is_casl': True, 'tau': 1.8}
#: An identifiable schedule: delays that bracket both tissues' arrival times.
PLDS = np.arange(0.25, 3.01, 0.25)


def _mixture(shape=(16, 16, 1), seed=0):
    rng = np.random.default_rng(seed)
    x = np.linspace(0, 1, shape[0])[:, None, None]
    y = np.linspace(0, 1, shape[1])[None, :, None]
    pv_gm = 0.2 + 0.6 * (0.5 + 0.5 * np.sin(5 * x + 3 * y)) * np.ones(shape)
    pv_wm = (1 - pv_gm) * (0.6 + 0.3 * rng.random(shape))
    return pv_gm, pv_wm


def _tissue_curve(cbf, att, m0, plds=PLDS, abv=0.0, abat=0.5, **kwargs):
    kinetics = {**KINETICS, **kwargs}
    return pvc.gkm_predict(
        np.array([cbf]),
        np.array([att]),
        np.array([abat]),
        np.array([abv]),
        np.array([m0]),
        plds[None, :],
        **kinetics,
    )[0]


def test_series_regression_matches_per_volume():
    """A series shares one design; each volume's amplitudes match a 3D fit of that volume."""
    pv_gm, pv_wm = _mixture()
    rng = np.random.default_rng(1)
    amplitudes = rng.uniform(5, 50, (2, 4))
    y = pv_gm[..., None] * amplitudes[0] + pv_wm[..., None] * amplitudes[1]
    y += rng.normal(0, 0.5, y.shape)
    mask = np.ones(pv_gm.shape, bool)
    series = pvc.local_regression_pvc(y, pv_gm, pv_wm, mask, mask)
    assert series['beta_gm'].shape == (*pv_gm.shape, 4)
    rss = 0
    for i in range(4):
        single = pvc.local_regression_pvc(y[..., i], pv_gm, pv_wm, mask, mask)
        np.testing.assert_array_equal(series['status'], single['status'])
        np.testing.assert_allclose(series['beta_gm'][..., i], single['beta_gm'], rtol=1e-12)
        np.testing.assert_allclose(series['noise_gm'], single['noise_gm'])
        dof = single['n'] - np.where(single['status'] == 1, 2, 1)
        rss = rss + single['rmse'] ** 2 * dof
    fitted = (series['status'] > 0) & (series['status'] < 4)
    dof = series['n'] - np.where(series['status'] == 1, 2, 1)
    np.testing.assert_allclose(
        series['rmse'][fitted] ** 2 * 4 * dof[fitted], rss[fitted], rtol=1e-10
    )


def test_series_average_identity():
    """The mean of per-volume amplitudes equals the amplitude of the mean (linearity)."""
    pv_gm, pv_wm = _mixture()
    rng = np.random.default_rng(2)
    y = rng.normal(30, 5, (*pv_gm.shape, 6))
    mask = np.ones(pv_gm.shape, bool)
    series = pvc.local_regression_pvc(y, pv_gm, pv_wm, mask, mask)
    mean = pvc.local_regression_pvc(y.mean(axis=-1), pv_gm, pv_wm, mask, mask)
    two = series['status'] == 1
    np.testing.assert_allclose(series['beta_gm'][two].mean(axis=-1), mean['beta_gm'][two])


def test_series_validity_requires_every_volume():
    """A non-finite sample in any volume excludes that voxel from every volume's fit."""
    pv_gm, pv_wm = _mixture()
    y = np.ones((*pv_gm.shape, 3)) * 10
    y[8, 8, 0, 2] = np.nan
    mask = np.ones(pv_gm.shape, bool)
    fit = pvc.local_regression_pvc(y, pv_gm, pv_wm, mask, mask)
    assert fit['status'][8, 8, 0] == 4


@pytest.mark.parametrize('arterial', [False, True])
def test_multidelay_exact_recovery(arterial):
    """Unmixing each delay, then fitting each tissue, recovers the calibrated tissue amplitude.

    The estimand is f_t * M0_t / m_bar_c (centre-M0 calibration), and ATT_t.
    """
    pv_gm, pv_wm = _mixture()
    m0_gm, m0_wm = 1000.0, 800.0
    abv = 0.01 if arterial else 0.0
    curve_gm = _tissue_curve(60.0, 1.0, m0_gm, abv=abv)
    curve_wm = _tissue_curve(20.0, 1.6, m0_wm, abv=abv)
    y = pv_gm[..., None] * curve_gm + pv_wm[..., None] * curve_wm
    m0_mixed = pv_gm * m0_gm + pv_wm * m0_wm + (1 - pv_gm - pv_wm) * 1500.0
    mask = np.ones(pv_gm.shape, bool)

    fit = pvc.local_regression_pvc(y, pv_gm, pv_wm, mask, mask)
    two = fit['status'] == 1
    assert two.sum() > 50
    np.testing.assert_allclose(
        fit['beta_gm'][two],
        np.broadcast_to(curve_gm, (two.sum(), PLDS.size)),
        rtol=1e-6,
        atol=1e-9,
    )

    for tissue, cbf_true, att_true, m0_true in (
        ('gm', 60.0, 1.0, m0_gm),
        ('wm', 20.0, 1.6, m0_wm),
    ):
        curves = fit[f'beta_{tissue}'][two]
        m0_c = m0_mixed[two]
        result = pvc.fit_tissue_curves(
            curves, m0_c, np.broadcast_to(PLDS, curves.shape), **KINETICS
        )
        assert np.all(result['status'] == 1)
        np.testing.assert_allclose(result['cbf'], cbf_true * m0_true / m0_c, rtol=0.01)
        np.testing.assert_allclose(result['att'], att_true, atol=0.01)
        assert np.all(result['att_upper'] - result['att_lower'] < 0.02)
        assert np.all(result['kinres'] < 1e-3)


def test_multidelay_slice_timing_single_slice_kernel():
    """With slice timing, a single-slice kernel shares its delays, so recovery is still exact."""
    pv_gm, pv_wm = _mixture(shape=(16, 16, 3))
    slice_times = np.array([0.0, 0.1, 0.2])
    plds = PLDS[None, None, None, :] + slice_times[None, None, :, None]
    m0_gm, m0_wm = 1000.0, 800.0
    y = np.zeros((*pv_gm.shape, PLDS.size))
    for k in range(3):
        y[:, :, k] = pv_gm[:, :, k, None] * _tissue_curve(
            60.0, 1.0, m0_gm, plds=plds[0, 0, k]
        ) + pv_wm[:, :, k, None] * _tissue_curve(20.0, 1.6, m0_wm, plds=plds[0, 0, k])
    mask = np.ones(pv_gm.shape, bool)
    fit = pvc.local_regression_pvc(y, pv_gm, pv_wm, mask, mask, kernel_size=(5, 5, 1))
    two = fit['status'] == 1
    result = pvc.fit_tissue_curves(
        fit['beta_gm'][two],
        np.full(two.sum(), m0_gm),
        np.broadcast_to(plds, (*pv_gm.shape, PLDS.size))[two],
        **KINETICS,
    )
    np.testing.assert_allclose(result['cbf'], 60.0, rtol=0.01)
    np.testing.assert_allclose(result['att'], 1.0, atol=0.01)


def test_fit_tissue_curves_edge_cases():
    """Negative curves, curves above the voxel bound, unidentifiable ATT and bad inputs."""
    m0 = 1000.0
    curves = np.stack(
        [
            _tissue_curve(60.0, 1.0, m0),
            -_tissue_curve(60.0, 1.0, m0),  # no positive tissue signal
            _tissue_curve(360.0, 1.0, m0),  # above the voxel fit's 300 bound
            _tissue_curve(2000.0, 1.0, m0),  # above the LR bound
            np.full(PLDS.size, np.nan),
        ]
    )
    result = pvc.fit_tissue_curves(
        curves, np.full(5, m0), np.broadcast_to(PLDS, curves.shape), **KINETICS
    )
    np.testing.assert_array_equal(result['status'], [1, 2, 1, 3, 4])
    assert result['cbf'][1] == 0
    # Not clipped at 300: the raised bound keeps the kinetics intact.
    assert result['cbf'][2] == pytest.approx(360.0, rel=0.01)
    assert result['att'][2] == pytest.approx(1.0, abs=0.01)
    assert result['cbf'][4] == 0


def test_fit_tissue_curves_unidentifiable_att():
    """With every delay longer than ATT, the ATT bounds are wide (ATT is not claimed)."""
    plds = np.arange(1.5, 3.01, 0.25)
    curve = _tissue_curve(60.0, 0.5, 1000.0, plds=plds)
    result = pvc.fit_tissue_curves(curve[None, :], np.array([1000.0]), plds[None, :], **KINETICS)
    assert result['status'][0] == 1
    assert result['att_upper'][0] - result['att_lower'][0] > 0.5
    assert result['cbf'][0] == pytest.approx(60.0, rel=0.01)


def test_csf_kinetic_projection_is_signed():
    """The CSF diagnostic keeps its sign, in units of the voxel's uncorrected CBF."""
    template = np.stack([_tissue_curve(50.0, 1.2, 1000.0), _tissue_curve(40.0, 1.2, 1000.0)])
    curves = np.stack([0.1 * template[0], -0.05 * template[1]])
    cbf, kinres, defined = pvc.csf_kinetic_projection(curves, template, np.array([50.0, 40.0]))
    np.testing.assert_allclose(cbf, [5.0, -2.0])
    np.testing.assert_allclose(kinres, 0, atol=1e-12)
    assert defined.all()
    cbf, _, defined = pvc.csf_kinetic_projection(curves, np.zeros_like(template), np.zeros(2))
    assert not defined.any()
    assert np.all(cbf == 0)


def test_csf_series_regression():
    """The three-tissue QC fit works on a series and recovers a CSF curve."""
    rng = np.random.default_rng(4)
    raw = rng.random((3, 16, 16, 1)) + 0.2
    raw /= raw.sum(axis=0)
    pv_gm, pv_wm, pv_csf = raw
    csf_curve = np.linspace(-2, 3, 5)
    y = (
        pv_gm[..., None] * np.linspace(40, 60, 5)
        + pv_wm[..., None] * np.linspace(10, 20, 5)
        + pv_csf[..., None] * csf_curve
    )
    mask = np.ones(pv_gm.shape, bool)
    fit = pvc.local_regression_csf(y, pv_gm, pv_wm, pv_csf, mask, mask)
    full = fit['status'] == 1
    assert full.sum() > 50
    np.testing.assert_allclose(
        fit['beta_csf'][full], np.broadcast_to(csf_curve, (full.sum(), 5)), atol=1e-6
    )
