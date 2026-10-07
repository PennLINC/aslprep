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
