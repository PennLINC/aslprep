"""Tests for the aslprep.utils.cbf module."""

import numpy as np

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
        fit_cbf, fit_att, fit_abat, fit_abv, failures = cbf.fit_deltam_multipld(
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
    for estimate in estimates[:4]:
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
    for serial_estimate, threaded_estimate in zip(serial, threaded, strict=True):
        np.testing.assert_array_equal(threaded_estimate, serial_estimate)


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
    fit_cbf, fit_att, _, fit_abv, failures = cbf.fit_deltam_multipld(
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
