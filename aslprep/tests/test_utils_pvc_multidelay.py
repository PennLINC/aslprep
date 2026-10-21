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


@pytest.mark.parametrize('partition_coefficient', [0.9, 1.0])
def test_fit_tissue_curves_partition_coefficient(partition_coefficient):
    """The fit and its kinetic residual use the same partition coefficient.

    It is 1 for an M0Estimate, which is already the M0 of blood.
    """
    kinetics = {**KINETICS, 'partition_coefficient': partition_coefficient}
    curve = _tissue_curve(60.0, 1.0, 1000.0, **kinetics)
    result = pvc.fit_tissue_curves(curve[None, :], np.array([1000.0]), PLDS[None, :], **kinetics)
    assert result['status'][0] == 1
    np.testing.assert_allclose(result['cbf'], 60.0, rtol=0.01)
    assert result['kinres'][0] < 1e-3


@pytest.mark.parametrize('abv', [0.01, 0.05, 0.1])
def test_fit_tissue_curves_pure_arterial(abv):
    """A curve with only arterial signal has no tissue CBF, so its ATT is not claimed.

    The solver returns ~1e-13 rather than exactly 0 here, which must count as the lower bound.
    """
    curve = _tissue_curve(0.0, 1.0, 1000.0, abv=abv, abat=0.4)
    result = pvc.fit_tissue_curves(curve[None, :], np.array([1000.0]), PLDS[None, :], **KINETICS)
    assert result['status'][0] == 2
    assert result['cbf'][0] == 0


def test_fit_tissue_curves_at_ceiling():
    """A curve at the CBF ceiling is flagged, although the solver stops just below it."""
    curve = _tissue_curve(pvc.LR_CBF_MAX, 1.0, 1000.0)
    result = pvc.fit_tissue_curves(curve[None, :], np.array([1000.0]), PLDS[None, :], **KINETICS)
    assert result['status'][0] == 3


def test_voxel_delays_refuses_non_finite_timing():
    """Non-finite delays or slice times are refused, not turned into zeros by the masker."""
    import nibabel as nb
    from nilearn import maskers

    from aslprep.utils.cbf import voxel_delays

    img = nb.Nifti1Image(np.zeros((2, 2, 3, 3)), np.eye(4))
    masker = maskers.NiftiMasker(mask_img=nb.Nifti1Image(np.ones((2, 2, 3)), np.eye(4))).fit()
    with pytest.raises(ValueError, match='PostLabelingDelay'):
        voxel_delays({'PostLabelingDelay': [0.5, np.nan, 1.5]}, img, masker, 12, True)
    with pytest.raises(ValueError, match='SliceTiming'):
        voxel_delays(
            {'PostLabelingDelay': [0.5, 1.0, 1.5], 'SliceTiming': [0, np.nan, 0.2]},
            img,
            masker,
            12,
            True,
        )
    plds, _ = voxel_delays(
        {'PostLabelingDelay': [0.5, 1.0, 1.5], 'SliceTiming': [0, 0.1, 0.2]}, img, masker, 12, True
    )
    assert plds.shape == (12, 3)
    assert np.all(np.isfinite(plds))


@pytest.mark.parametrize('bad', [np.nan, np.inf, -1.0, 0.0])
def test_multipld_bolus_refuses_impossible_bolus(bad):
    """A labeling duration or TI1 that is not finite and positive is refused."""
    from aslprep.utils.cbf import multipld_bolus

    with pytest.raises(ValueError, match='LabelingDuration'):
        multipld_bolus({'LabelingDuration': bad}, is_casl=True)
    with pytest.raises(ValueError, match='LabelingDuration'):
        multipld_bolus({'LabelingDuration': [1.8, bad, 1.8]}, is_casl=True)
    with pytest.raises(ValueError, match='BolusCutOffDelayTime'):
        multipld_bolus(
            {'BolusCutOffTechnique': 'QUIPSSII', 'BolusCutOffDelayTime': bad}, is_casl=False
        )
    assert multipld_bolus({'LabelingDuration': [1.8, 1.5]}, is_casl=True) == ([1.8, 1.5], None)
    assert multipld_bolus(
        {'BolusCutOffTechnique': 'QUIPSSII', 'BolusCutOffDelayTime': 0.7}, is_casl=False
    ) == (None, 0.7)


def test_estimate_lr_pvc_mem_gb():
    """The memory estimate scales with the series size."""
    small = pvc.estimate_lr_pvc_mem_gb((16, 16, 8, 84))
    large = pvc.estimate_lr_pvc_mem_gb((64, 64, 30, 84), num_threads=3)
    assert small < 1.2
    assert large > 2.0  # the 1.33 GB of regression arrays measured in review, plus overhead


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
    """The CSF diagnostic is the signed CBF of the curve against a unit-CBF basis."""
    basis = np.stack([_tissue_curve(1.0, 1.2, 1000.0), _tissue_curve(1.0, 0.8, 1200.0)])
    curves = np.stack([5.0 * basis[0], -2.0 * basis[1]])
    cbf, kinres, defined = pvc.csf_kinetic_projection(curves, basis)
    np.testing.assert_allclose(cbf, [5.0, -2.0])
    np.testing.assert_allclose(kinres, 0, atol=1e-12)
    assert defined.all()
    cbf, _, defined = pvc.csf_kinetic_projection(curves, np.zeros_like(basis))
    assert not defined.any()
    assert np.all(cbf == 0)
    curves[1, 0] = np.nan
    _, _, defined = pvc.csf_kinetic_projection(curves, basis)
    np.testing.assert_array_equal(defined, [True, False])


def test_fit_tissue_curves_large_abv():
    """A tissue curve with an arterial volume above the voxel fit's 0.1 bound is recovered.

    With the voxel bound, the excess arterial signal was absorbed into tissue CBF and ATT.
    """
    curve = _tissue_curve(60.0, 1.0, 1000.0, abv=0.12, abat=0.5)
    result = pvc.fit_tissue_curves(curve[None, :], np.array([1000.0]), PLDS[None, :], **KINETICS)
    assert result['status'][0] == 1
    assert result['cbf'][0] == pytest.approx(60.0, rel=0.01)
    assert result['att'][0] == pytest.approx(1.0, abs=0.01)


def test_fit_tissue_curves_abv_at_bound():
    """A fit whose aBV is at its bound is flagged (status 5), so its CBF is not claimed."""
    curve = _tissue_curve(60.0, 1.0, 1000.0, abv=0.12, abat=0.5)
    result = pvc.fit_tissue_curves(
        curve[None, :], np.array([1000.0]), PLDS[None, :], abv_max=0.05, **KINETICS
    )
    assert result['status'][0] == 5


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
