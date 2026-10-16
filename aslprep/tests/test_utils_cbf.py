"""Tests for the aslprep.utils.cbf module."""

import numpy as np
import pytest

from aslprep.utils import cbf


def test_estimate_att_pcasl_3d():
    """Smoke test aslprep.utils.cbf.estimate_att_pcasl with 3D data (no slice timing)."""
    n_voxels, n_unique_plds = 1000, 6
    base_plds = np.linspace(0.25, 1.5, n_unique_plds)[None, :]
    plds = np.repeat(base_plds, n_voxels, axis=0)
    deltam_arr = np.random.random((n_voxels, n_unique_plds))
    tau = np.full(n_unique_plds, 1.4)
    att = cbf.estimate_att_pcasl(
        deltam_arr=deltam_arr,
        plds=plds,
        lds=tau,
        t1blood=1.6,
        t1tissue=1.3,
    )
    assert att.shape == (n_voxels,)


def test_estimate_att_pcasl_2d():
    """Smoke test aslprep.utils.cbf.estimate_att_pcasl with slice-shifted PLDs."""
    n_voxels, n_unique_plds, n_slice_times = 1000, 6, 50
    base_plds = np.linspace(0.25, 1.5, n_unique_plds)
    slice_times = np.linspace(0, 3, n_slice_times)
    slice_shifted_plds = slice_times[:, None] + base_plds[None, :]
    plds = np.repeat(slice_shifted_plds, n_voxels // n_slice_times, axis=0)
    deltam_arr = np.random.random((n_voxels, n_unique_plds))
    tau = np.full(n_unique_plds, 1.4)
    att = cbf.estimate_att_pcasl(
        deltam_arr=deltam_arr,
        plds=plds,
        lds=tau,
        t1blood=1.6,
        t1tissue=1.3,
    )
    assert att.shape == (n_voxels,)


def test_estimate_cbf_pcasl_multipld():
    """Smoke test aslprep.utils.cbf.estimate_cbf_pcasl_multipld with slice-shifted PLDs."""
    n_voxels, n_unique_plds, n_slice_times = 1000, 6, 50
    base_plds = np.linspace(0.25, 1.5, n_unique_plds)
    slice_times = np.linspace(0, 3, n_slice_times)
    slice_shifted_plds = slice_times[:, None] + base_plds[None, :]
    plds = np.repeat(slice_shifted_plds, n_voxels // n_slice_times, axis=0)
    deltam_arr = np.random.random((n_voxels, n_unique_plds))
    scaled_m0data = np.random.random(n_voxels)

    cbf.estimate_cbf_pcasl_multipld(
        deltam_arr=deltam_arr,
        scaled_m0data=scaled_m0data,
        plds=plds,
        tau=np.array(1.4),
        labeleff=0.7,
        t1blood=1.6,
        t1tissue=1.3,
        unit_conversion=6000,
        partition_coefficient=0.9,
    )


def _score_phantom(n_vol=40, seed=0, soft=False):
    """Build a 16 x 16 x 4 phantom with GM (CBF 60), WM (20), and a CSF strip (0)."""
    rng = np.random.default_rng(seed)
    gm = np.zeros((16, 16, 4))
    wm = np.zeros_like(gm)
    csf = np.zeros_like(gm)
    gm[:7], wm[9:], csf[7:9] = 1, 1, 1
    truth = 60 * gm + 20 * wm
    if soft:
        # Probability maps that never reach exactly 1
        gm, wm, csf = gm * 0.98, wm * 0.98, csf * 0.98

    cbf_ts = truth[..., None] + 15 * rng.standard_normal((16, 16, 4, n_vol))
    return cbf_ts, gm, wm, csf, np.ones_like(gm)


def _pooled_variance(cbf_ts, keep, gm, wm, csf, thresh=0.7):
    mean_cbf = cbf_ts[..., keep].mean(axis=-1)
    return sum(
        ((tpm >= thresh).sum() - 1) * np.var(mean_cbf[tpm >= thresh]) for tpm in (gm, wm, csf)
    )


def test_getcbfscore_soft_tpms():
    """SCORE's second pass must use the thresholded masks, not voxels where the TPM equals 1."""
    cbf_ts, gm, wm, csf, mask = _score_phantom(seed=7)
    bad = cbf_ts[..., 0] + 200 * gm
    cbf_ts = np.concatenate([cbf_ts[..., 1:], bad[..., None]], axis=-1)
    _, index_binary = cbf._getcbfscore(cbf_ts, wm, gm, csf, mask)

    soft = [tpm * 0.98 for tpm in (gm, wm, csf)]
    _, index_soft = cbf._getcbfscore(cbf_ts, soft[1], soft[0], soft[2], mask)
    np.testing.assert_array_equal(index_binary, index_soft)
    assert index_soft[-1] == 1


def test_getcbfscore_order_invariant():
    """SCORE's retained set should not depend on the order of the volumes."""
    cbf_ts, gm, wm, csf, mask = _score_phantom(seed=7)
    bad = cbf_ts[..., 0] + 200 * gm
    rest = [cbf_ts[..., i] for i in range(1, cbf_ts.shape[3])]

    _, index_first = cbf._getcbfscore(np.stack([bad] + rest, axis=-1), wm, gm, csf, mask)
    _, index_last = cbf._getcbfscore(np.stack(rest + [bad], axis=-1), wm, gm, csf, mask)

    # The outlier is rejected by the first pass, regardless of position
    assert index_first[0] == 1
    assert index_last[-1] == 1
    # The remaining volumes get the same labels
    np.testing.assert_array_equal(index_first[1:], index_last[:-1])

    rng = np.random.default_rng(0)
    order = rng.permutation(cbf_ts.shape[3])
    _, index = cbf._getcbfscore(cbf_ts, wm, gm, csf, mask)
    _, index_perm = cbf._getcbfscore(cbf_ts[..., order], wm, gm, csf, mask)
    np.testing.assert_array_equal(index[order], index_perm)


def test_getcbfscore_never_increases_variance():
    """SCORE's second pass must not keep a removal that raised the pooled variance."""
    for seed in range(10):
        cbf_ts, gm, wm, csf, mask = _score_phantom(seed=seed)
        _, index = cbf._getcbfscore(cbf_ts, wm, gm, csf, mask)
        v_final = _pooled_variance(cbf_ts, index == 0, gm, wm, csf)
        v_start = _pooled_variance(cbf_ts, index != 1, gm, wm, csf)
        assert v_final <= v_start


def _score_phantom_with_artifact(seed=3):
    """Add a zero-mean checkerboard to 12 volumes, which SCORE's second pass partly removes."""
    cbf_ts, gm, wm, csf, mask = _score_phantom(seed=seed)
    checker = np.indices(gm.shape).sum(axis=0) % 2 * 2 - 1.0
    cbf_ts[..., :12] += 60 * checker[..., None]
    return cbf_ts, gm, wm, csf, mask


def test_getcbfscore_empty_tissue_mask():
    """An empty tissue mask must not disable SCORE's structural-correlation pass."""
    cbf_ts, gm, wm, _, mask = _score_phantom_with_artifact()
    # A CSF mask with no voxels and one with a single voxel inside GM both carry no
    # within-tissue variance, so they must give the same labels
    csf_empty = np.zeros_like(gm)
    csf_one = np.zeros_like(gm)
    csf_one[0, 0, 0] = 1
    _, index_empty = cbf._getcbfscore(cbf_ts, wm, gm, csf_empty, mask)
    _, index_one = cbf._getcbfscore(cbf_ts, wm, gm, csf_one, mask)
    np.testing.assert_array_equal(index_empty, index_one)
    assert np.any(index_empty == 2)


def test_getcbfscore_constant_volume():
    """A volume with an undefined correlation must never be picked as the most correlated."""
    cbf_ts, gm, wm, csf, mask = _score_phantom_with_artifact()
    # Its GM mean is near the median, so the first pass keeps it. Rounding keeps it exactly
    # constant, so its correlation is NaN rather than roundoff noise.
    flat = np.full(gm.shape, np.round(np.median(cbf_ts[gm >= 0.7].mean(axis=0))))
    cbf_ts = np.concatenate([flat[..., None], cbf_ts], axis=-1)
    _, index = cbf._getcbfscore(cbf_ts, wm, gm, csf, mask)
    assert index[0] == 0
    assert np.any(index[1:] == 2)


def _simulate_multipld(is_casl, plds, bolus, params, m0, labeleff=0.7, t1blood=1.65):
    """Simulate delta-M for each voxel with the scalar reference models."""
    n_voxels, n_volumes = plds.shape
    model = cbf.calculate_deltam_pcasl if is_casl else cbf.calculate_deltam_pasl
    deltam = np.empty((n_voxels, n_volumes))
    for i_voxel in range(n_voxels):
        X = np.zeros((n_volumes, 6))
        X[0, :4] = labeleff, t1blood, m0[i_voxel] / 0.9, m0[i_voxel] / 0.9
        X[:, 4] = plds[i_voxel]
        if is_casl:
            X[:, 5] = bolus
        else:
            X[0, 5] = bolus

        deltam[i_voxel] = model(X, *params[i_voxel])

    return deltam


def _random_params(rng, n_voxels):
    return np.column_stack(
        [
            rng.uniform(10, 120, n_voxels),
            rng.uniform(0.3, 2.5, n_voxels),
            rng.uniform(0.2, 2.0, n_voxels),
            rng.uniform(0, 0.05, n_voxels),
        ]
    )


def _flatten_fit(fit):
    """Stack every array returned by fit_deltam_multipld, including the bounds."""
    *estimates, bounds = fit
    return np.vstack([*estimates, *(bounds[key] for key in sorted(bounds))])


def test_gkm_bases_match_reference_models():
    """The vectorized tissue and arterial terms reproduce the scalar reference models."""
    rng = np.random.default_rng(0)
    n_voxels = 200
    params = _random_params(rng, n_voxels)
    m0 = rng.uniform(500, 1500, n_voxels)
    for is_casl, delays, bolus in (
        (True, np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0]), np.array([1.8, 1.8, 1.8, 1.5, 1.5, 1.5])),
        (False, np.array([0.6, 1.0, 1.4, 1.8, 2.2, 2.6]), 0.7),
    ):
        plds = np.tile(delays, (n_voxels, 1)) + rng.uniform(0, 0.5, (n_voxels, 1))
        expected = _simulate_multipld(is_casl, plds, bolus, params, m0)
        kwargs = {'is_casl': is_casl, 'labeleff': 0.7, 't1blood': 1.65, 'xp': np}
        tissue = cbf._gkm_tissue_basis(params[:, 1:2], plds, bolus, **kwargs)
        arterial = cbf._gkm_arterial_basis(params[:, 2:3], plds, bolus, **kwargs)
        predicted = (m0 / 0.9)[:, None] * (params[:, 0:1] * tissue + params[:, 3:4] * arterial)
        np.testing.assert_allclose(predicted, expected, rtol=1e-12, atol=1e-12)


def test_solve_box_lsq_matches_lsq_linear():
    """The closed-form bounded solve agrees with SciPy's bounded linear least squares."""
    from scipy.optimize import lsq_linear

    rng = np.random.default_rng(1)
    (c_lo, c_hi), _, _, (v_lo, v_hi) = cbf.MULTIPLD_BOUNDS
    for _ in range(500):
        design = np.abs(rng.normal(size=(6, 2))) * [1e-4, 1.0]
        design[rng.random(6) < 0.3, 1] = 0
        y = design @ [rng.uniform(-50, 400), rng.uniform(-0.05, 0.2)] + rng.normal(0, 3e-3, 6)
        t, a = design.T
        c, v, loss = cbf._solve_box_lsq(t @ y, a @ y, t @ t, a @ a, t @ a, 0.0, np)
        reference = lsq_linear(design, y, bounds=([c_lo, v_lo], [c_hi, v_hi]), tol=1e-14)
        reference_loss = np.sum((design @ reference.x - y) ** 2) - y @ y
        assert c_lo <= c <= c_hi
        assert v_lo <= v <= v_hi
        assert loss <= reference_loss + 1e-10 * (y @ y)


def test_fit_deltam_multipld_recovers_noiseless_data():
    """Noiseless data are fit exactly in every voxel.

    The fitted parameters are not compared to the truth, because the arterial term can absorb
    part of the tissue signal, so different parameters can fit the data equally well.
    """
    rng = np.random.default_rng(2)
    n_voxels = 300
    params = _random_params(rng, n_voxels)
    m0 = rng.uniform(500, 1500, n_voxels)
    for is_casl, delays, bolus in (
        (True, np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0]), 1.8),
        (False, np.array([0.6, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0]), 0.7),
    ):
        # Slice-timing shifts give each voxel its own delays
        plds = np.tile(delays, (n_voxels, 1)) + rng.choice(np.linspace(0, 0.5, 10), (n_voxels, 1))
        deltam = _simulate_multipld(is_casl, plds, bolus, params, m0)
        fit_cbf, fit_att, fit_abat, fit_abv, failures, bounds = cbf.fit_deltam_multipld(
            deltam_arr=deltam,
            scaled_m0data=m0,
            plds=plds,
            labeleff=0.7,
            t1blood=1.65,
            partition_coefficient=0.9,
            is_casl=is_casl,
            tau=bolus if is_casl else None,
            ti1=None if is_casl else bolus,
            chunk_size=128,
        )
        assert not failures.any()
        fitted = np.column_stack([fit_cbf, fit_att, fit_abat, fit_abv])
        predicted = _simulate_multipld(is_casl, plds, bolus, fitted, m0)
        np.testing.assert_allclose(predicted, deltam, atol=1e-5 * np.abs(deltam).max())
        # The estimates lie within their bounds, and the aBAT bounds are model edges
        assert np.all(bounds['att_lower'] <= fit_att)
        assert np.all(fit_att <= bounds['att_upper'])
        assert np.all(bounds['abat_lower'] <= fit_abat)
        assert np.all(fit_abat <= bounds['abat_upper'])
        edges = cbf._delay_edges(plds, bolus, is_casl, np)
        assert np.all(np.isclose(bounds['abat_lower'][:, None], edges).any(axis=1))


def test_fit_deltam_multipld_repeated_volumes():
    """Repeating every volume leaves the quality of the fit unchanged.

    Parameters are not compared directly, because rounding can break ties between equally good
    fits differently.
    """
    rng = np.random.default_rng(3)
    n_voxels = 50
    plds = np.tile([0.5, 1.0, 1.5, 2.0, 2.5], (n_voxels, 1))
    m0 = rng.uniform(500, 1500, n_voxels)
    deltam = _simulate_multipld(True, plds, 1.8, _random_params(rng, n_voxels), m0)
    deltam += rng.normal(0, 0.5, deltam.shape)
    kwargs = {
        'scaled_m0data': m0,
        'labeleff': 0.7,
        't1blood': 1.65,
        'partition_coefficient': 0.9,
        'is_casl': True,
        'tau': 1.8,
    }
    single = cbf.fit_deltam_multipld(deltam_arr=deltam, plds=plds, **kwargs)
    repeated = cbf.fit_deltam_multipld(
        deltam_arr=np.hstack([deltam, deltam]),
        plds=np.hstack([plds, plds]),
        **kwargs,
    )
    sse = [
        np.sum((_simulate_multipld(True, plds, 1.8, np.column_stack(fit[:4]), m0) - deltam) ** 2)
        for fit in (single, repeated)
    ]
    np.testing.assert_allclose(sse[1], sse[0], rtol=1e-8)


def test_fit_deltam_multipld_invalid_voxels():
    """Voxels with non-positive M0 or non-finite data are flagged as failures."""
    plds = np.tile([0.5, 1.0, 1.5, 2.0], (4, 1))
    deltam = np.ones((4, 4))
    deltam[2, 1] = np.nan
    m0 = np.array([1000.0, 0.0, 1000.0, -5.0])
    estimates = cbf.fit_deltam_multipld(
        deltam_arr=deltam,
        scaled_m0data=m0,
        plds=plds,
        labeleff=0.7,
        t1blood=1.65,
        partition_coefficient=0.9,
        is_casl=True,
        tau=[1.8, 1.8, 1.8, 1.8],
    )
    np.testing.assert_array_equal(estimates[4], [False, True, True, True])
    for estimate in [*estimates[:4], *estimates[5].values()]:
        assert np.isfinite(estimate[0])
        assert np.isnan(estimate[1:]).all()


def test_fit_deltam_multipld_threads():
    """Fitting chunks in several threads gives exactly the same estimates as one thread."""
    rng = np.random.default_rng(4)
    n_voxels = 200
    plds = np.tile([0.5, 1.0, 1.5, 2.0, 2.5], (n_voxels, 1))
    m0 = rng.uniform(500, 1500, n_voxels)
    m0[::17] = 0  # failed voxels must not shift the others
    deltam = _simulate_multipld(True, plds, 1.8, _random_params(rng, n_voxels), m0)
    deltam += rng.normal(0, 0.5, deltam.shape)
    kwargs = {
        'deltam_arr': deltam,
        'scaled_m0data': m0,
        'plds': plds,
        'labeleff': 0.7,
        't1blood': 1.65,
        'partition_coefficient': 0.9,
        'is_casl': True,
        'tau': 1.8,
        'chunk_size': 16,
    }
    serial = cbf.fit_deltam_multipld(**kwargs)
    threaded = cbf.fit_deltam_multipld(**kwargs, n_jobs=4)
    np.testing.assert_array_equal(_flatten_fit(threaded), _flatten_fit(serial))


def test_fit_deltam_multipld_unidentifiable_att():
    """An ATT shorter than every delay is reported as the midpoint of the range it could take.

    Once all of the labeled blood has arrived by the first delay, the signal no longer depends
    on ATT, so any ATT from zero to that delay fits equally well.
    """
    n_voxels = 20
    rng = np.random.default_rng(5)
    first_delays = rng.uniform(0.5, 1.0, (n_voxels, 1))
    plds = first_delays + np.array([0.0, 0.5, 1.0, 1.5, 2.0])
    params = np.column_stack(
        [
            rng.uniform(20, 90, n_voxels),
            rng.uniform(0.05, 0.45, n_voxels),
            np.full(n_voxels, 1.0),
            np.zeros(n_voxels),
        ]
    )
    m0 = rng.uniform(500, 1500, n_voxels)
    deltam = _simulate_multipld(True, plds, 1.8, params, m0)
    fit_cbf, fit_att, _, fit_abv, failures, bounds = cbf.fit_deltam_multipld(
        deltam_arr=deltam,
        scaled_m0data=m0,
        plds=plds,
        labeleff=0.7,
        t1blood=1.65,
        partition_coefficient=0.9,
        is_casl=True,
        tau=1.8,
    )
    assert not failures.any()
    np.testing.assert_allclose(fit_att, first_delays[:, 0] / 2, atol=1e-9)
    np.testing.assert_allclose(fit_cbf, params[:, 0], rtol=1e-6)
    np.testing.assert_allclose(fit_abv, 0, atol=1e-9)
    np.testing.assert_allclose(bounds['att_lower'], 0, atol=1e-9)
    # Just past the first delay, the fit is worse by less than the tolerance that defines a tie
    np.testing.assert_allclose(bounds['att_upper'], first_delays[:, 0], atol=1e-5)
    # Without arterial signal, aBAT is not constrained at all
    np.testing.assert_allclose(bounds['abat_lower'], 0, atol=1e-9)
    np.testing.assert_allclose(bounds['abat_upper'], cbf.MULTIPLD_BOUNDS[2][1], atol=1e-9)


def test_fit_deltam_multipld_att_bounds():
    """The ATT bounds contain the true ATT, and collapse onto it where the data locate it.

    The arterial term can absorb the part of the tissue signal that ATT changes,
    so some voxels have a range of equally good ATTs even though ATT lies between the delays.
    """
    n_voxels = 20
    rng = np.random.default_rng(7)
    plds = np.tile([0.5, 1.0, 1.5, 2.0, 2.5], (n_voxels, 1))
    params = np.column_stack(
        [
            rng.uniform(20, 90, n_voxels),
            rng.uniform(1.1, 1.9, n_voxels),
            np.full(n_voxels, 1.0),
            np.zeros(n_voxels),
        ]
    )
    m0 = rng.uniform(500, 1500, n_voxels)
    deltam = _simulate_multipld(True, plds, 1.8, params, m0)
    _, fit_att, _, _, _, bounds = cbf.fit_deltam_multipld(
        deltam_arr=deltam,
        scaled_m0data=m0,
        plds=plds,
        labeleff=0.7,
        t1blood=1.65,
        partition_coefficient=0.9,
        is_casl=True,
        tau=1.8,
    )
    assert np.all(bounds['att_lower'] <= params[:, 1] + 1e-5)
    assert np.all(params[:, 1] <= bounds['att_upper'] + 1e-5)
    assert np.all(bounds['att_lower'] <= fit_att)
    assert np.all(fit_att <= bounds['att_upper'])
    located = bounds['att_upper'] - bounds['att_lower'] < 1e-4
    assert located.any()
    assert not located.all()
    np.testing.assert_allclose(fit_att[located], params[located, 1], atol=1e-4)


def test_fit_deltam_multipld_no_signal():
    """Without any signal, neither ATT nor aBAT is constrained."""
    plds = np.tile([0.5, 1.0, 1.5, 2.0], (3, 1))
    for is_casl in (True, False):
        fit_cbf, _, _, fit_abv, _, bounds = cbf.fit_deltam_multipld(
            deltam_arr=np.zeros((3, 4)),
            scaled_m0data=np.full(3, 1000.0),
            plds=plds,
            labeleff=0.7,
            t1blood=1.65,
            partition_coefficient=0.9,
            is_casl=is_casl,
            tau=1.8 if is_casl else None,
            ti1=None if is_casl else 0.7,
        )
        np.testing.assert_array_equal(fit_cbf, 0)
        np.testing.assert_array_equal(fit_abv, 0)
        for i_param, name in ((1, 'att'), (2, 'abat')):
            lower, upper = cbf.MULTIPLD_BOUNDS[i_param]
            np.testing.assert_allclose(bounds[f'{name}_lower'], lower, atol=1e-9)
            np.testing.assert_allclose(bounds[f'{name}_upper'], upper, atol=1e-9)


def test_fit_deltam_multipld_abv_at_bound():
    """A (P)CASL arterial signal that needs the largest aBV restricts aBAT within its interval.

    Only aBV * exp(-aBAT / T1b) is identifiable, so an aBAT later in the interval needs a larger
    aBV for the same signal. Once that aBV would exceed its bound, the later aBATs are ruled out.
    """
    t1blood = 1.65
    abv_max = cbf.MULTIPLD_BOUNDS[3][1]
    plds = np.tile([0.5, 1.0, 1.5, 2.0, 2.5], (3, 1))
    # The aBAT interval is 0.5-1.0 s. The first voxel can take any aBAT in it,
    # the second only those up to 0.6 s, and the third only those just above the lower edge.
    true_abat = np.array([0.6, 0.6, 0.501])
    true_abv = np.array([0.02, abv_max, abv_max])
    params = np.column_stack([np.full(3, 60.0), np.full(3, 1.7), true_abat, true_abv])
    m0 = np.full(3, 1000.0)
    deltam = _simulate_multipld(True, plds, 1.8, params, m0)
    fit_cbf, fit_att, fit_abat, fit_abv, _, bounds = cbf.fit_deltam_multipld(
        deltam_arr=deltam,
        scaled_m0data=m0,
        plds=plds,
        labeleff=0.7,
        t1blood=t1blood,
        partition_coefficient=0.9,
        is_casl=True,
        tau=1.8,
    )
    fitted = np.column_stack([fit_cbf, fit_att, fit_abat, fit_abv])
    predicted = _simulate_multipld(True, plds, 1.8, fitted, m0)
    np.testing.assert_allclose(predicted, deltam, atol=1e-7 * np.abs(deltam).max())
    assert np.all(fit_abv <= abv_max)
    np.testing.assert_allclose(bounds['abat_lower'], 0.5, atol=1e-9)
    np.testing.assert_allclose(bounds['abat_upper'], [1.0, 0.6, 0.501], atol=1e-6)
    np.testing.assert_allclose(fit_abat, [0.75, 0.55, 0.5005], atol=1e-6)
    # The identifiable product is recovered in every voxel
    np.testing.assert_allclose(
        fit_abv * np.exp(-fit_abat / t1blood),
        true_abv * np.exp(-true_abat / t1blood),
        rtol=1e-6,
    )


def test_fit_deltam_multipld_abv_beyond_bound():
    """An arterial signal too large for the aBV bound puts aBAT at the bottom of its interval.

    The arterial signal per unit aBV is largest there, so the bounded fit is best there.
    The interval does not include its lower edge, so aBAT is reported just above it.
    """
    plds = np.tile([0.5, 1.0, 1.5, 2.0, 2.5], (1, 1))
    abv_max = cbf.MULTIPLD_BOUNDS[3][1]
    m0 = np.full(1, 1000.0)
    deltam = _simulate_multipld(True, plds, 1.8, np.array([[60.0, 1.7, 0.6, 2 * abv_max]]), m0)
    fit_cbf, fit_att, fit_abat, fit_abv, _, bounds = cbf.fit_deltam_multipld(
        deltam_arr=deltam,
        scaled_m0data=m0,
        plds=plds,
        labeleff=0.7,
        t1blood=1.65,
        partition_coefficient=0.9,
        is_casl=True,
        tau=1.8,
    )
    assert fit_abat[0] > 0.5
    np.testing.assert_allclose(fit_abat, 0.5, atol=1e-8)
    np.testing.assert_allclose(fit_abv, abv_max)
    np.testing.assert_allclose(bounds['abat_lower'], 0.5)
    np.testing.assert_allclose(bounds['abat_upper'], fit_abat)

    def sse(abat):
        fitted = np.column_stack([fit_cbf, fit_att, [abat], fit_abv])
        return np.sum((_simulate_multipld(True, plds, 1.8, fitted, m0) - deltam) ** 2)

    # The reported aBAT is in the interval, and fits better than its midpoint does
    assert sse(fit_abat[0]) < sse(0.75)
    assert sse(fit_abat[0]) < sse(0.5)


def test_fit_gkm_chunk_abat_bounds_stop_at_abv_bound():
    """The aBAT bounds do not cross the part of an interval that the aBV bound rules out.

    With two delays, several aBAT intervals fit these data exactly. In the reported interval,
    the longest aBATs would need more than the largest aBV, so the range ends inside it,
    even though the next interval fits as well.
    """
    kwargs = {'is_casl': True, 'labeleff': 0.7, 't1blood': 1.65, 'xp': np}
    y = np.array([[0.054, 0.052]])
    plds = np.array([[0.5, 1.0]])
    bolus = np.array([[1.8, 1.8]])
    cbf_, att, abat, abv, _, _, abat_lower, abat_upper = cbf._fit_gkm_chunk(
        y, plds, bolus, np.ones((1, 2)), **kwargs
    )
    assert abat_lower[0] <= abat[0] <= abat_upper[0] < 2.3
    tissue = cbf._gkm_tissue_basis(att[:, None], plds, bolus, **kwargs)

    def sse(abat, abv):
        arterial = cbf._gkm_arterial_basis(np.array([[abat]]), plds, bolus, **kwargs)
        return np.sum((cbf_[0] * tissue + abv * arterial - y) ** 2)

    # Every aBAT in the bounds fits as well as the estimate, with aBV within its bound
    best = sse(abat[0], abv[0])
    for fraction in (0.01, 0.5, 0.99):
        other = abat_lower[0] + fraction * (abat_upper[0] - abat_lower[0])
        other_abv = abv[0] * np.exp((other - abat[0]) / 1.65)
        assert other_abv <= cbf.MULTIPLD_BOUNDS[3][1]
        np.testing.assert_allclose(sse(other, other_abv), best, atol=1e-20)


def test_get_array_module(monkeypatch, tmp_path):
    """The GPU array library is returned only when it is importable and finds a device."""
    import os
    import sys

    from aslprep import config
    from aslprep.tests.tests import stand_in_cupy

    monkeypatch.delenv('CUPY_CACHE_DIR', raising=False)
    monkeypatch.setattr(config.execution, 'work_dir', tmp_path)
    assert cbf.get_array_module() is np
    assert cbf.get_array_module(use_gpu=False) is np
    assert 'CUPY_CACHE_DIR' not in os.environ

    fake = stand_in_cupy()
    monkeypatch.setitem(sys.modules, 'cupy', fake)
    assert cbf.get_array_module(use_gpu=True) is fake
    # Compiled kernels are cached in the working directory, unless the user chose a location
    assert os.environ['CUPY_CACHE_DIR'] == str(tmp_path / 'cupy_cache')
    monkeypatch.setenv('CUPY_CACHE_DIR', '/somewhere/else')
    cbf.get_array_module(use_gpu=True)
    assert os.environ['CUPY_CACHE_DIR'] == '/somewhere/else'

    monkeypatch.setitem(sys.modules, 'cupy', stand_in_cupy(n_devices=0))
    with pytest.raises(RuntimeError, match='did not find a CUDA device'):
        cbf.get_array_module(use_gpu=True)

    broken = stand_in_cupy()

    def _raise():
        raise OSError('no driver')

    broken.cuda.runtime.getDeviceCount = _raise
    monkeypatch.setitem(sys.modules, 'cupy', broken)
    with pytest.raises(RuntimeError, match='could not query CUDA devices: no driver'):
        cbf.get_array_module(use_gpu=True)

    # None in sys.modules makes the import fail, as if CuPy were not installed
    monkeypatch.setitem(sys.modules, 'cupy', None)
    with pytest.raises(RuntimeError, match='requires CuPy'):
        cbf.get_array_module(use_gpu=True)


def test_fit_deltam_multipld_other_array_module(monkeypatch):
    """Another array library gives the same estimates, with chunks sized by memory."""
    from aslprep.tests.tests import stand_in_cupy

    rng = np.random.default_rng(6)
    n_voxels = 60
    plds = np.tile([0.5, 1.0, 1.5, 2.0, 2.5], (n_voxels, 1))
    m0 = rng.uniform(500, 1500, n_voxels)
    deltam = _simulate_multipld(True, plds, 1.8, _random_params(rng, n_voxels), m0)
    deltam += rng.normal(0, 0.5, deltam.shape)
    kwargs = {
        'deltam_arr': deltam,
        'scaled_m0data': m0,
        'plds': plds,
        'labeleff': 0.7,
        't1blood': 1.65,
        'partition_coefficient': 0.9,
        'is_casl': True,
        'tau': 1.8,
    }
    expected = cbf.fit_deltam_multipld(**kwargs)

    chunk_sizes = []
    fit_chunk = cbf._fit_gkm_chunk

    def _record(y, **chunk_kwargs):
        chunk_sizes.append(y.shape[0])
        return fit_chunk(y=y, **chunk_kwargs)

    monkeypatch.setattr(cbf, '_fit_gkm_chunk', _record)
    # 5 delays give 11 intervals, so each voxel needs 11 * 11 * 11 elements: 25 voxels fit
    monkeypatch.setattr(cbf, 'MULTIPLD_GPU_CHUNK_ELEMENTS', 25 * 11**3)
    estimates = cbf.fit_deltam_multipld(**kwargs, xp=stand_in_cupy())
    assert chunk_sizes == [25, 25, 10]
    np.testing.assert_array_equal(_flatten_fit(estimates), _flatten_fit(expected))


def test_zero_low_m0():
    """M0 far below the brain's level is zeroed; everything else is untouched."""
    m0 = np.full(1000, 100.0)
    m0[:300] = np.linspace(-5, 9.9, 300)  # air: noise around zero, below 10 % of 100
    m0[300:320] = np.linspace(10, 60, 20)  # partial volume at the edge: kept
    m0[320] = np.nan

    out = cbf.zero_low_m0(m0)
    assert np.all(out[:300] == 0)
    np.testing.assert_array_equal(out[300:320], m0[300:320])
    assert np.isnan(out[320])
    assert np.all(out[321:] == 100)
    # The input is not modified.
    assert m0[0] == -5

    # A higher threshold removes more of the edge.
    assert np.count_nonzero(cbf.zero_low_m0(m0, threshold=0.5) == 0) > 300

    # A constant M0 (an M0Estimate) is never zeroed, and non-finite data pass through.
    np.testing.assert_array_equal(cbf.zero_low_m0(np.full(10, 3.0)), np.full(10, 3.0))
    assert np.isnan(cbf.zero_low_m0(np.full(10, np.nan))).all()


@pytest.mark.parametrize(('haematocrit', 'expected'), [(0.47, 1.70), (0.20, 2.09)])
def test_t1blood_from_haematocrit(haematocrit, expected):
    """The model reproduces the arterial examples in Hales et al. (2016), at 3 T."""
    assert cbf.t1blood_from_haematocrit(haematocrit, 3.0) == pytest.approx(expected, abs=0.005)


@pytest.mark.parametrize(('field_strength', 'default'), [(1.5, 1.35), (3, 1.65), (7, 2.087)])
def test_estimate_t1_haematocrit(field_strength, default):
    """A measured haematocrit scales the default T1blood, which a typical one leaves unchanged."""
    t1blood, t1tissue = cbf.estimate_t1({'MagneticFieldStrength': field_strength})
    assert t1blood == default

    typical = cbf.estimate_t1({'MagneticFieldStrength': field_strength, 'Haematocrit': 0.42})
    assert typical == pytest.approx((default, t1tissue))

    # Fewer red blood cells lengthen T1blood, more shorten it; T1tissue is unaffected.
    low = cbf.estimate_t1({'MagneticFieldStrength': field_strength, 'Haematocrit': 0.20})
    high = cbf.estimate_t1({'MagneticFieldStrength': field_strength, 'Haematocrit': 0.47})
    assert low[0] > default > high[0]
    assert low[1] == high[1] == t1tissue
    expected_ratio = cbf.t1blood_from_haematocrit(0.20, field_strength) / (
        cbf.t1blood_from_haematocrit(0.42, field_strength)
    )
    assert low[0] == pytest.approx(default * expected_ratio)


def test_estimate_t1_haematocrit_value_at_3t():
    """The documented 3 T value for a haematocrit of 0.20."""
    t1blood, _ = cbf.estimate_t1({'MagneticFieldStrength': 3, 'Haematocrit': 0.20})
    assert t1blood == pytest.approx(1.948, abs=0.0005)


@pytest.mark.parametrize(
    ('metadata', 'expected'),
    [
        ({'MagneticFieldStrength': 3}, False),
        ({'MagneticFieldStrength': 3, 'Haematocrit': 0.4}, True),
        ({'MagneticFieldStrength': 1.5, 'Haematocrit': 0.4}, True),
        ({'MagneticFieldStrength': 7, 'Haematocrit': 0.4}, True),
        # Outside the field strengths the model was fitted on
        ({'MagneticFieldStrength': 1.0, 'Haematocrit': 0.4}, False),
        ({'MagneticFieldStrength': 9.4, 'Haematocrit': 0.4}, False),
    ],
)
def test_haematocrit_applies(metadata, expected):
    """The haematocrit is only used within the model's 1.5-7 T range."""
    assert cbf.haematocrit_applies(metadata) is expected
