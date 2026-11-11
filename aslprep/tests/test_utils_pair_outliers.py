"""Tests for the aslprep.utils.pair_outliers module."""

import numpy as np
import pytest

from aslprep.tests.test_utils_cbf import _random_params, _simulate_multipld
from aslprep.utils import pair_outliers

FIT_KWARGS = {
    'labeleff': 0.7,
    't1blood': 1.65,
    'partition_coefficient': 0.9,
}


def test_predict_deltam_multipld_matches_reference_models():
    """Predictions reproduce the scalar reference models, with slice-shifted delays."""
    rng = np.random.default_rng(0)
    n_voxels = 100
    params = _random_params(rng, n_voxels)
    m0 = rng.uniform(500, 1500, n_voxels)
    for is_casl, delays, bolus in (
        (True, np.array([0.5, 1.0, 1.0, 2.0, 2.5]), 1.8),
        (False, np.array([0.6, 1.0, 1.4, 1.4, 2.6]), 0.7),
    ):
        plds = np.tile(delays, (n_voxels, 1)) + rng.uniform(0, 0.5, (n_voxels, 1))
        expected = _simulate_multipld(is_casl, plds, bolus, params, m0)
        predicted = pair_outliers.predict_deltam_multipld(
            *params.T,
            scaled_m0data=m0,
            plds=plds,
            is_casl=is_casl,
            tau=bolus if is_casl else None,
            ti1=None if is_casl else bolus,
            **FIT_KWARGS,
        )
        np.testing.assert_allclose(predicted, expected, rtol=1e-10, atol=1e-10)


def test_delay_groups_follow_fit_grouping():
    """Pairs group by delay and labeling duration, not by slice shift."""
    plds = np.array([[1.0, 1.5, 1.0, 1.5, 1.0], [1.2, 1.7, 1.2, 1.7, 1.2]])
    groups = pair_outliers._delay_groups(plds, np.array([1.8, 1.8, 1.8, 1.8, 1.5]), True)
    assert groups[0] == groups[2]
    assert groups[1] == groups[3]
    assert len(set(groups)) == 3


def _noisy_pairs(rng, n_voxels=400, repeats=4, noise=0.5):
    delays = np.repeat([0.5, 1.0, 1.5, 2.0, 2.5], repeats)
    plds = np.tile(delays, (n_voxels, 1))
    params = _random_params(rng, n_voxels)
    m0 = rng.uniform(800, 1200, n_voxels)
    clean = _simulate_multipld(True, plds, 1.8, params, m0)
    deltam = clean + rng.normal(0, noise, clean.shape)
    groups = pair_outliers._delay_groups(plds, np.full(delays.size, 1.8), True)
    return deltam, clean, plds, m0, groups


@pytest.mark.parametrize('shrinkage', [0.0, 2.0, np.inf])
def test_score_pairs_ranks_corrupted_pairs(shrinkage):
    """Corrupted pairs have the largest energies; a pair without signal loses its amplitude."""
    rng = np.random.default_rng(1)
    deltam, clean, _, _, groups = _noisy_pairs(rng)
    deltam[:, 6] += 3.0
    deltam[:, 13] -= clean[:, 13]
    keep = np.ones(deltam.shape[1], dtype=bool)
    scores = pair_outliers.score_pairs(deltam, clean, groups, keep, shrinkage=shrinkage)
    assert set(np.argsort(scores['energy'])[-2:]) == {6, 13}
    clean_pairs = np.setdiff1d(np.arange(deltam.shape[1]), [6, 13])
    # The median-based noise estimate is a little high with few pairs, so clean pairs score
    # a little under one.
    assert 0.5 < np.median(scores['energy'][clean_pairs]) < 1.2
    assert scores['energy'][6] > 10 * np.max(scores['energy'][clean_pairs])
    if np.isinf(shrinkage):
        # Against the model alone, the pair without labeling has lost all of its signal.
        np.testing.assert_allclose(scores['amplitude'][13], -1, atol=0.05)
    assert np.argmin(scores['amplitude_t']) == 13


def test_score_pairs_same_delay_correction_needs_three_others():
    """With fewer than three other pairs at a delay, pairs are predicted by the model alone."""
    rng = np.random.default_rng(5)
    deltam, clean, _, _, groups = _noisy_pairs(rng, repeats=3)
    keep = np.ones(deltam.shape[1], dtype=bool)
    model = pair_outliers.score_pairs(deltam, clean, groups, keep, shrinkage=np.inf)
    empirical = pair_outliers.score_pairs(deltam, clean, groups, keep, shrinkage=0.0)
    np.testing.assert_allclose(empirical['energy'], model['energy'])

    deltam, clean, _, _, groups = _noisy_pairs(rng, repeats=4)
    keep = np.ones(deltam.shape[1], dtype=bool)
    model = pair_outliers.score_pairs(deltam, clean, groups, keep, shrinkage=np.inf)
    empirical = pair_outliers.score_pairs(deltam, clean, groups, keep, shrinkage=0.0)
    assert not np.allclose(empirical['energy'], model['energy'])


def test_score_pairs_smoothing_raises_coherent_residuals():
    """Smoothing raises the energy of a spatially coherent residual relative to noise."""
    rng = np.random.default_rng(2)
    deltam, clean, _, _, groups = _noisy_pairs(rng)
    deltam[:, 6] += 0.5
    keep = np.ones(deltam.shape[1], dtype=bool)

    def smooth(values):
        # Voxels are in a line here; average each with its neighbors.
        kernel = np.ones(5) / 5
        return np.column_stack(
            [np.convolve(col, kernel, mode='same') for col in np.asarray(values).T]
        )

    scores = pair_outliers.score_pairs(
        deltam, clean, groups, keep, shrinkage=np.inf, smooth=smooth
    )
    others = np.delete(np.arange(deltam.shape[1]), 6)
    raw_ratio = scores['energy'][6] / np.median(scores['energy'][others])
    smoothed_ratio = scores['energy_smooth'][6] / np.median(scores['energy_smooth'][others])
    assert smoothed_ratio > 2 * raw_ratio


def test_flag_pairs_effect_size_floor():
    """Tiny differences between nearly identical pairs are not flagged."""
    energy = np.ones(12)
    energy[3] = 1.01
    scores = {'energy': energy, 'amplitude': np.zeros(12), 'amplitude_t': np.zeros(12)}
    flags, _ = pair_outliers.flag_pairs(scores, np.zeros(12, dtype=int))
    assert not flags.any()

    energy[3] = 3.0
    flags, _ = pair_outliers.flag_pairs(scores, np.zeros(12, dtype=int))
    assert np.flatnonzero(flags).tolist() == [3]


def test_flag_pairs_keeps_every_delay():
    """The least severe flagged pairs are kept, so that no delay loses all of its pairs."""
    energy = np.ones(12)
    energy[[0, 1]] = [10.0, 5.0]
    groups = np.array([0, 0] + [1] * 10)
    scores = {'energy': energy, 'amplitude': np.zeros(12), 'amplitude_t': np.zeros(12)}
    flags, _ = pair_outliers.flag_pairs(scores, groups, min_pairs_per_delay=1)
    assert np.flatnonzero(flags).tolist() == [0]
    flags, _ = pair_outliers.flag_pairs(scores, groups, min_pairs_per_delay=2)
    assert not flags.any()


@pytest.mark.parametrize('shrinkage', [2.0, np.inf])
def test_detect_outlier_pairs_multipld(shrinkage):
    """Corrupted pairs are excluded from the fit, and clean data lose no pairs."""
    rng = np.random.default_rng(3)
    deltam, clean, plds, m0, _ = _noisy_pairs(rng, n_voxels=300)
    kwargs = {
        'scaled_m0data': m0,
        'plds': plds,
        'is_casl': True,
        'tau': 1.8,
        'shrinkage': shrinkage,
        **FIT_KWARGS,
    }
    result = pair_outliers.detect_outlier_pairs_multipld(deltam_arr=deltam, **kwargs)
    assert result['keep'].all()

    corrupted = deltam.copy()
    corrupted[:, 6] += rng.normal(0, 3.0, deltam.shape[0])
    corrupted[:, 13] -= clean[:, 13]
    result = pair_outliers.detect_outlier_pairs_multipld(deltam_arr=corrupted, **kwargs)
    assert np.flatnonzero(~result['keep']).tolist() == [6, 13]
    assert result['fit'][0].shape == (deltam.shape[0],)
    assert result['n_fits'] >= 2


def test_detect_outlier_pairs_multipld_pasl():
    """PASL data (inversion times and TI1) are handled, and a corrupted pair is excluded."""
    rng = np.random.default_rng(6)
    n_voxels = 300
    delays = np.repeat([0.8, 1.2, 1.6, 2.0, 2.4], 4)
    plds = np.tile(delays, (n_voxels, 1))
    params = _random_params(rng, n_voxels)
    m0 = rng.uniform(800, 1200, n_voxels)
    clean = _simulate_multipld(False, plds, 0.7, params, m0)
    deltam = clean + rng.normal(0, 0.5, clean.shape)
    deltam[:, 9] += rng.normal(0, 3.0, n_voxels)
    result = pair_outliers.detect_outlier_pairs_multipld(
        deltam_arr=deltam, scaled_m0data=m0, plds=plds, is_casl=False, ti1=0.7, **FIT_KWARGS
    )
    assert np.flatnonzero(~result['keep']).tolist() == [9]


def test_detect_outlier_pairs_multipld_per_volume_tau():
    """Pairs with the same delay but different labeling durations are separate groups."""
    rng = np.random.default_rng(7)
    deltam, _, plds, m0, _ = _noisy_pairs(rng, n_voxels=200)
    tau = np.full(deltam.shape[1], 1.8)
    tau[::4] = 1.5  # One pair at each delay has a shorter labeling duration
    deltam[:, 5] += 3.0
    result = pair_outliers.detect_outlier_pairs_multipld(
        deltam_arr=deltam, scaled_m0data=m0, plds=plds, is_casl=True, tau=tau, **FIT_KWARGS
    )
    # The pairs with tau = 1.5 are each alone in their group, so they are always kept.
    assert result['keep'][::4].all()
    assert not result['keep'][5]


def test_detect_outlier_pairs_multipld_input_checks():
    """max_iter below one and an empty scoring mask are rejected."""
    rng = np.random.default_rng(8)
    deltam, _, plds, m0, _ = _noisy_pairs(rng, n_voxels=50)
    kwargs = {'scaled_m0data': m0, 'plds': plds, 'is_casl': True, 'tau': 1.8, **FIT_KWARGS}
    with pytest.raises(ValueError, match='max_iter'):
        pair_outliers.detect_outlier_pairs_multipld(deltam_arr=deltam, max_iter=0, **kwargs)
    with pytest.raises(ValueError, match='No voxels'):
        pair_outliers.detect_outlier_pairs_multipld(
            deltam_arr=deltam, mask=np.zeros(deltam.shape[0], dtype=bool), **kwargs
        )


def test_detect_outlier_pairs_multipld_one_pair_per_delay():
    """With one pair per delay, no pair can be spared, so nothing is scored or excluded."""
    rng = np.random.default_rng(9)
    deltam, _, plds, m0, _ = _noisy_pairs(rng, n_voxels=50, repeats=1)
    deltam[:, 2] += 5.0
    result = pair_outliers.detect_outlier_pairs_multipld(
        deltam_arr=deltam, scaled_m0data=m0, plds=plds, is_casl=True, tau=1.8, **FIT_KWARGS
    )
    assert result['keep'].all()
    assert result['scores'] is None
    assert result['n_fits'] == 1


def test_detect_outlier_pairs_multipld_too_few_pairs():
    """With fewer than MIN_PAIRS pairs, every pair is kept."""
    rng = np.random.default_rng(4)
    deltam, _, plds, m0, _ = _noisy_pairs(rng, n_voxels=20, repeats=1)
    result = pair_outliers.detect_outlier_pairs_multipld(
        deltam_arr=deltam[:, :3],
        scaled_m0data=m0,
        plds=plds[:, :3],
        is_casl=True,
        tau=1.8,
        **FIT_KWARGS,
    )
    assert result['keep'].all()
    assert result['scores'] is None
