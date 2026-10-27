"""Tests for the tissue-T1 kinetic model in aslprep.utils.cbf."""

import numpy as np
import pytest

from aslprep.tests.truth_models import tissue_t1_signal
from aslprep.utils import cbf

T1BLOOD = 1.65
LABELEFF = 0.7
M0 = 1000.0
#: Apparent T1s of gray and white matter at 3 T, with the default nominal CBFs
T1PRIME = np.array(
    [
        cbf.apparent_tissue_t1(cbf.T1TISSUE_BY_FIELD_STRENGTH[3], cbf.NOMINAL_CBF['gm']),
        cbf.apparent_tissue_t1(cbf.T1WM_BY_FIELD_STRENGTH[3], cbf.NOMINAL_CBF['wm']),
    ]
)
PCASL = (True, 1.8, np.arange(0.5, 4.01, 0.5))
PASL = (False, 0.8, np.arange(0.6, 3.61, 0.3))


def _simulate(is_casl, bolus, delays, cbf_values, atts, t1primes, weights):
    """Simulate delta-M as a sum of tissues, each with its own CBF, ATT, T1' and weight.

    The signal comes from numerical integration (:func:`tissue_t1_signal`), not from the
    closed-form basis under test.
    """
    times = delays + bolus if is_casl else delays
    signal = np.zeros(delays.size)
    for tissue_cbf, att, t1prime, weight in zip(cbf_values, atts, t1primes, weights, strict=True):
        tissue = [tissue_t1_signal(t, att, bolus, T1BLOOD, t1prime, is_casl) for t in times]
        signal += weight * tissue_cbf * np.array(tissue)
    return 2 * LABELEFF * (M0 / 0.9) / 6000 * signal


def _fit(deltam, delays, is_casl, bolus, **kwargs):
    deltam = np.atleast_2d(deltam)
    return cbf.fit_deltam_multipld(
        deltam_arr=deltam,
        scaled_m0data=np.full(deltam.shape[0], M0),
        plds=np.broadcast_to(delays, deltam.shape).copy(),
        labeleff=LABELEFF,
        t1blood=T1BLOOD,
        partition_coefficient=0.9,
        is_casl=is_casl,
        tau=bolus if is_casl else None,
        ti1=None if is_casl else bolus,
        **kwargs,
    )


@pytest.mark.parametrize('is_casl', [True, False])
def test_relaxed_basis_matches_numerical_integration(is_casl):
    """The closed-form tissue terms match direct integration of the delivered label."""
    rng = np.random.default_rng(0)
    bolus = 1.8 if is_casl else 0.8
    delays = np.linspace(0.3, 4.0, 15)
    times = delays + bolus if is_casl else delays
    for att in rng.uniform(0.1, 3.0, 20):
        for t1prime in (*T1PRIME, T1BLOOD, 2.5):
            basis = cbf._relaxed_tissue_basis(
                np.array([[att]]),
                delays[None],
                bolus,
                is_casl,
                T1BLOOD,
                np.array([[t1prime]]),
                np,
            )[0]
            expected = [tissue_t1_signal(t, att, bolus, T1BLOOD, t1prime, is_casl) for t in times]
            np.testing.assert_allclose(basis, expected, rtol=1e-10, atol=1e-14)


@pytest.mark.parametrize('is_casl', [True, False])
@pytest.mark.parametrize('relative_change', [0.0, 1e-12, -1e-12, 1e-9])
def test_relaxed_basis_reduces_to_consensus(is_casl, relative_change):
    """With T1' at (or within rounding of) blood T1, the tissue-T1 basis is the consensus one."""
    rng = np.random.default_rng(1)
    bolus = 1.8 if is_casl else 0.8
    att = rng.uniform(0, 3, (50, 1))
    delays = np.tile(np.linspace(0.3, 4.0, 12), (50, 1))
    kwargs = {'is_casl': is_casl, 'labeleff': LABELEFF, 't1blood': T1BLOOD, 'xp': np}
    consensus = cbf._gkm_tissue_basis(att, delays, bolus, **kwargs)
    relaxed = cbf._gkm_tissue_basis(
        att,
        delays,
        bolus,
        **kwargs,
        t1prime=np.full((50, 1), T1BLOOD * (1 + relative_change)),
        tissue_weights=np.ones((50, 1)),
    )
    scale = np.abs(consensus).max()
    np.testing.assert_allclose(relaxed, consensus, rtol=1e-8, atol=1e-12 * scale)


@pytest.mark.parametrize('is_casl', [True, False])
def test_relaxed_basis_is_continuous_at_edges(is_casl):
    """The tissue-T1 basis is continuous where ATT crosses each delay edge."""
    bolus = 1.8 if is_casl else 0.8
    delays = np.array([[0.8, 1.6, 2.4]])
    edges = np.concatenate([delays, delays + bolus] if is_casl else [delays - bolus, delays])
    for edge in edges.ravel():
        below, above = (
            cbf._relaxed_tissue_basis(
                np.array([[att]]),
                delays,
                bolus,
                is_casl,
                T1BLOOD,
                np.array([[T1PRIME[1]]]),
                np,
            )
            for att in (edge - 1e-9, edge + 1e-9)
        )
        np.testing.assert_allclose(below, above, atol=1e-8)


def test_mixed_basis_is_weighted_sum():
    """A mixed voxel's basis is the weighted sum of its tissues' bases."""
    att = np.array([[0.7], [1.1], [1.6]])
    delays = np.tile(np.arange(0.5, 4.01, 0.5), (3, 1))
    weights = np.array([[1.0, 0.0], [0.5, 0.5], [0.2, 0.8]])
    kwargs = {'is_casl': True, 'labeleff': LABELEFF, 't1blood': T1BLOOD, 'xp': np}
    mixed = cbf._gkm_tissue_basis(
        att,
        delays,
        1.8,
        **kwargs,
        t1prime=np.tile(T1PRIME, (3, 1)),
        tissue_weights=weights,
    )
    pure = [
        cbf._gkm_tissue_basis(
            att,
            delays,
            1.8,
            **kwargs,
            t1prime=np.full((3, 1), t1prime),
            tissue_weights=np.ones((3, 1)),
        )
        for t1prime in T1PRIME
    ]
    np.testing.assert_allclose(mixed, weights[:, :1] * pure[0] + weights[:, 1:] * pure[1])


def test_apparent_tissue_t1():
    """Flow shortens the apparent T1, with CBF in mL/100 g/min."""
    assert cbf.apparent_tissue_t1(1.607, 0.0) == pytest.approx(1.607)
    assert cbf.apparent_tissue_t1(1.607, 60.0) == pytest.approx(1 / (1 / 1.607 + 0.01 / 0.9))
    assert cbf.apparent_tissue_t1(1.607, 60.0, 1.0) == pytest.approx(1 / (1 / 1.607 + 0.01))


def test_tissue_mixing_weights():
    """Weights are partial volume times nominal CBF, normalized, without CSF."""
    gm = np.array([1.0, 0.0, 0.5, 0.3, 0.0, 0.005, np.nan])
    wm = np.array([0.0, 1.0, 0.5, 0.3, 0.0, 0.004, 0.5])
    weights = cbf.tissue_mixing_weights(gm, wm, 60.0, 20.0)
    np.testing.assert_allclose(
        weights,
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [0.75, 0.25],
            # CSF fills the rest of the voxel, but does not change the weights
            [0.75, 0.25],
            # Almost no gray or white matter: the gray matter basis
            [1.0, 0.0],
            [1.0, 0.0],
            # A missing partial volume is treated as zero
            [0.0, 1.0],
        ],
    )


@pytest.mark.parametrize(('is_casl', 'bolus', 'delays'), [PCASL, PASL])
def test_fit_tissue_t1_recovers_truth(is_casl, bolus, delays):
    """Noise-free data from the tissue-T1 model give back its CBF and ATT.

    The consensus model, fit to the same data, underestimates CBF and ATT in white matter,
    as on aslscan's simulations (gh-744).
    """
    cases = (
        ([1.0, 0.0], 60.0, 1.0),  # gray matter
        ([0.0, 1.0], 20.0, 1.3),  # white matter
        ([0.6, 0.4], 45.0, 1.1),  # a mixed voxel
    )
    deltam = np.vstack(
        [
            _simulate(is_casl, bolus, delays, [true_cbf] * 2, [att] * 2, T1PRIME, weights)
            for weights, true_cbf, att in cases
        ]
    )
    weights = np.array([case[0] for case in cases])
    true_cbf = np.array([case[1] for case in cases])
    true_att = np.array([case[2] for case in cases])
    fit_cbf, fit_att, _, fit_abv, failures, _ = _fit(
        deltam,
        delays,
        is_casl,
        bolus,
        t1prime=T1PRIME,
        tissue_weights=weights,
    )
    assert not failures.any()
    np.testing.assert_allclose(fit_cbf, true_cbf, rtol=1e-3)
    np.testing.assert_allclose(fit_att, true_att, atol=1e-3)
    np.testing.assert_allclose(fit_abv, 0, atol=1e-6)

    consensus_cbf, consensus_att, *_ = _fit(deltam[1:2], delays, is_casl, bolus)
    assert consensus_cbf[0] < 0.7 * true_cbf[1]
    assert consensus_att[0] < true_att[1] - 0.05


def test_fit_tissue_t1_unequal_tissue_atts():
    """With one ATT per voxel, a mix of tissues with different ATTs is fit approximately.

    The model gives gray and white matter one ATT, so this records how far that approximation is
    from a voxel whose tissues have aslscan's default ATTs (0.8 s and 1.2 s).
    """
    is_casl, bolus, delays = PCASL
    deltam = _simulate(is_casl, bolus, delays, [60.0, 20.0], [0.8, 1.2], T1PRIME, [0.6, 0.4])
    true_cbf = 0.6 * 60 + 0.4 * 20
    weights = cbf.tissue_mixing_weights(np.array([0.6]), np.array([0.4]), 60.0, 20.0)
    fit_cbf, fit_att, *_ = _fit(
        deltam,
        delays,
        is_casl,
        bolus,
        t1prime=T1PRIME,
        tissue_weights=weights,
    )
    assert fit_cbf[0] == pytest.approx(true_cbf, rel=0.02)
    assert 0.8 < fit_att[0] < 1.2


def test_fit_tissue_t1_mixed_basis_beats_average_t1():
    """A voxel-averaged T1 misfits a mixed voxel that the sum of tissue bases fits exactly."""
    is_casl, bolus, delays = PCASL
    weights = np.array([[0.75, 0.25]])
    deltam = _simulate(is_casl, bolus, delays, [44.0] * 2, [1.0] * 2, T1PRIME, weights[0])
    mixed_cbf, mixed_att, *_ = _fit(
        deltam,
        delays,
        is_casl,
        bolus,
        t1prime=T1PRIME,
        tissue_weights=weights,
    )
    average_t1 = 1 / np.sum(weights[0] / T1PRIME)
    average_cbf, *_ = _fit(deltam, delays, is_casl, bolus, t1prime=average_t1)
    assert mixed_cbf[0] == pytest.approx(44.0, rel=1e-3)
    assert abs(mixed_att[0] - 1.0) < 1e-3
    # About 5% for this voxel
    assert abs(average_cbf[0] / 44.0 - 1) > 0.03


def test_fit_tissue_t1_sensitivity_to_t1():
    """An error in tissue T1 biases CBF, in opposite directions for shorter and longer T1s.

    In gray matter, a 10% error in T1 moves CBF by about 10% (gh-744).
    """
    is_casl, bolus, delays = PCASL
    deltam = _simulate(is_casl, bolus, delays, [60.0], [1.0], T1PRIME[:1], [1.0])
    ratios = []
    for scale in (0.9, 1.0, 1.1):
        t1prime = cbf.apparent_tissue_t1(cbf.T1TISSUE_BY_FIELD_STRENGTH[3] * scale, 60.0)
        fit_cbf, *_ = _fit(deltam, delays, is_casl, bolus, t1prime=t1prime)
        ratios.append(fit_cbf[0] / 60.0)

    short, exact, long = ratios
    assert exact == pytest.approx(1.0, rel=1e-3)
    assert 1.03 < short < 1.2
    assert 0.8 < long < 0.97


@pytest.mark.parametrize('is_casl', [True, False])
def test_fit_tissue_t1_chunks_and_threads(is_casl):
    """The tissue-T1 fit, including its ATT interval, does not depend on chunks or threads."""
    rng = np.random.default_rng(3)
    n_voxels = 40
    bolus = 1.8 if is_casl else 0.8
    delays = np.tile(np.arange(0.5, 3.6, 0.5), (n_voxels, 1)) + rng.uniform(0, 0.4, (n_voxels, 1))
    weights = cbf.tissue_mixing_weights(
        rng.uniform(0, 1, n_voxels),
        rng.uniform(0, 1, n_voxels),
        60.0,
        20.0,
    )
    kwargs = {'is_casl': is_casl, 'labeleff': LABELEFF, 't1blood': T1BLOOD, 'xp': np}
    tissue = cbf._gkm_tissue_basis(
        rng.uniform(0.3, 2.0, (n_voxels, 1)),
        delays,
        bolus,
        **kwargs,
        t1prime=np.tile(T1PRIME, (n_voxels, 1)),
        tissue_weights=weights,
    )
    deltam = (M0 / 0.9) * rng.uniform(20, 80, (n_voxels, 1)) * tissue
    deltam += rng.normal(0, 0.05 * np.abs(deltam).max(), deltam.shape)

    def fit(**kwargs):
        cbf_, att, abat, abv, _, bounds = cbf.fit_deltam_multipld(
            deltam_arr=deltam,
            scaled_m0data=np.full(n_voxels, M0),
            plds=delays,
            labeleff=LABELEFF,
            t1blood=T1BLOOD,
            partition_coefficient=0.9,
            is_casl=is_casl,
            tau=bolus if is_casl else None,
            ti1=None if is_casl else bolus,
            att_interval=True,
            t1prime=T1PRIME,
            tissue_weights=weights,
            **kwargs,
        )
        return np.vstack([cbf_, att, abat, abv, *(bounds[key] for key in sorted(bounds))])

    reference = fit(chunk_size=n_voxels)
    np.testing.assert_array_equal(fit(chunk_size=7), reference)
    np.testing.assert_array_equal(fit(chunk_size=7, n_jobs=3), reference)


def test_fit_tissue_t1_other_array_module():
    """Another array library gives the same tissue-T1 estimates, including the ATT interval."""
    from aslprep.tests.tests import stand_in_cupy

    is_casl, bolus, delays = PCASL
    deltam = np.vstack(
        [
            _simulate(is_casl, bolus, delays, [60.0] * 2, [att] * 2, T1PRIME, [0.7, 0.3])
            for att in (0.6, 1.0, 1.4)
        ]
    )
    deltam += np.random.default_rng(5).normal(0, 0.02 * np.abs(deltam).max(), deltam.shape)
    kwargs = {
        't1prime': T1PRIME,
        'tissue_weights': np.array([0.7, 0.3]),
        'att_interval': True,
    }
    expected = _fit(deltam, delays, is_casl, bolus, **kwargs)
    estimates = _fit(deltam, delays, is_casl, bolus, xp=stand_in_cupy(), **kwargs)
    for value, expected_value in zip(estimates[:5], expected[:5], strict=True):
        np.testing.assert_array_equal(value, expected_value)

    for key, bound in expected[5].items():
        np.testing.assert_array_equal(estimates[5][key], bound)


def test_fit_tissue_t1_gpu_matches_cpu():
    """The tissue-T1 fit gives the same results on a GPU."""
    cupy = pytest.importorskip('cupy')
    try:
        if cupy.cuda.runtime.getDeviceCount() < 1:
            pytest.skip('No CUDA device')
    except Exception:  # noqa: BLE001
        pytest.skip('No CUDA device')

    is_casl, bolus, delays = PCASL
    deltam = np.vstack(
        [
            _simulate(is_casl, bolus, delays, [60.0] * 2, [att] * 2, T1PRIME, [0.7, 0.3])
            for att in (0.6, 1.0, 1.4)
        ]
    )
    kwargs = {'t1prime': T1PRIME, 'tissue_weights': np.array([0.7, 0.3])}
    cpu = _fit(deltam, delays, is_casl, bolus, **kwargs)
    gpu = _fit(deltam, delays, is_casl, bolus, xp=cupy, **kwargs)
    for cpu_value, gpu_value in zip(cpu[:4], gpu[:4], strict=True):
        np.testing.assert_allclose(gpu_value, cpu_value, rtol=1e-9, atol=1e-12)


def test_fit_tissue_t1_scalar_arguments():
    """A single tissue can be given as scalars, to the fit and to gkm_predict."""
    from aslprep.utils.pvc import gkm_predict

    is_casl, bolus, delays = PCASL
    deltam = _simulate(is_casl, bolus, delays, [20.0], [1.3], T1PRIME[1:], [1.0])
    fit_cbf, fit_att, *_ = _fit(
        deltam,
        delays,
        is_casl,
        bolus,
        t1prime=T1PRIME[1],
        tissue_weights=1.0,
    )
    assert fit_cbf[0] == pytest.approx(20.0, rel=1e-3)
    assert fit_att[0] == pytest.approx(1.3, abs=1e-3)
    predicted = gkm_predict(
        np.array([20.0]),
        np.array([1.3]),
        np.array([0.0]),
        np.array([0.0]),
        np.array([M0]),
        delays[None],
        labeleff=LABELEFF,
        t1blood=T1BLOOD,
        is_casl=is_casl,
        tau=bolus,
        t1prime=T1PRIME[1],
        tissue_weights=1.0,
    )
    np.testing.assert_allclose(predicted[0], deltam, rtol=1e-10)


def test_fit_tissue_t1_att_interval():
    """The noise-aware ATT interval uses the tissue-T1 model.

    In white matter at long delays, the tissue-T1 interval contains the true ATT in most voxels.
    The consensus model misfits these data, which inflates its pooled noise estimate,
    so its intervals are much wider.
    """
    is_casl, bolus, delays = PCASL
    n_voxels = 200
    clean = _simulate(is_casl, bolus, delays, [20.0], [1.3], T1PRIME[1:], [1.0])
    rng = np.random.default_rng(7)
    deltam = clean + rng.normal(0, 0.05 * np.abs(clean).max(), (n_voxels, delays.size))

    def interval(**kwargs):
        *_, bounds = _fit(deltam, delays, is_casl, bolus, att_interval=True, **kwargs)
        lower, upper = bounds['att_lower_ci'], bounds['att_upper_ci']
        return np.mean((lower <= 1.3) & (1.3 <= upper)), np.median(upper - lower)

    tissue_coverage, tissue_width = interval(t1prime=T1PRIME[1])
    _, consensus_width = interval()
    assert tissue_coverage > 0.8
    assert tissue_width < 0.75 * consensus_width


def test_fit_tissue_t1_att_before_every_delay():
    """With ATT shorter than every delay, a single tissue's ATT trades off against CBF.

    After arrival, ATT scales the signal by exp(ATT (1 / T1' - 1 / T1b)) at every delay, so a
    range of ATTs fits equally well, as in the consensus model. Unlike the consensus model,
    CBF then depends on the ATT reported from that range (its midpoint).
    """
    is_casl, bolus = True, 1.8
    delays = np.arange(1.0, 4.01, 0.5)
    deltam = _simulate(is_casl, bolus, delays, [40.0], [0.1], T1PRIME[1:], [1.0])
    fit_cbf, fit_att, _, _, _, bounds = _fit(deltam, delays, is_casl, bolus, t1prime=T1PRIME[1])
    # The bisection locates the ends of the range to a few microseconds
    assert bounds['att_lower'][0] == pytest.approx(0.0, abs=1e-5)
    assert bounds['att_upper'][0] == pytest.approx(delays[0], abs=1e-5)
    assert fit_att[0] == pytest.approx(delays[0] / 2, abs=1e-5)
    rate = 1 / T1PRIME[1] - 1 / T1BLOOD
    assert fit_cbf[0] == pytest.approx(40.0 * np.exp(rate * (0.1 - fit_att[0])), rel=1e-3)


def test_fit_tissue_t1_abv_beyond_bound():
    """The aBV bound and its flag work with the tissue-T1 basis.

    An arterial signal too large for the aBV bound puts aBAT at the bottom of its interval,
    as with the consensus model, and the voxel is flagged.
    """
    is_casl, bolus = True, 1.8
    delays = np.array([0.5, 1.0, 1.5, 2.0, 2.5])
    abv_max = cbf.MULTIPLD_BOUNDS[3][1]
    tissue = _simulate(is_casl, bolus, delays, [20.0], [1.7], T1PRIME[1:], [1.0])
    arterial = cbf._gkm_arterial_basis(
        np.array([[0.6]]),
        delays[None],
        bolus,
        is_casl,
        LABELEFF,
        T1BLOOD,
        np,
    )[0]
    deltam = tissue + (M0 / 0.9) * 2 * abv_max * arterial
    _, _, fit_abat, fit_abv, failures, bounds = _fit(
        deltam,
        delays,
        is_casl,
        bolus,
        t1prime=T1PRIME[1],
    )
    assert fit_abv[0] == pytest.approx(abv_max)
    assert fit_abat[0] > 0.5
    assert fit_abat[0] == pytest.approx(0.5, abs=1e-8)
    assert bounds['abat_lower'][0] == pytest.approx(0.5)
    assert cbf.abv_at_bound(fit_abv, failures)[0]


def test_fit_tissue_t1_rejects_bad_t1prime():
    """Non-positive apparent T1s are refused."""
    is_casl, bolus, delays = PCASL
    deltam = np.ones((1, delays.size))
    with pytest.raises(ValueError, match='t1prime must be finite and positive'):
        _fit(deltam, delays, is_casl, bolus, t1prime=np.array([1.0, 0.0]))


def test_gkm_predict_tissue_t1():
    """gkm_predict uses the same tissue-T1 basis as the fit."""
    from aslprep.utils.pvc import gkm_predict

    is_casl, bolus, delays = PCASL
    deltam = _simulate(is_casl, bolus, delays, [50.0], [1.2], T1PRIME[1:], [1.0])
    predicted = gkm_predict(
        np.array([50.0]),
        np.array([1.2]),
        np.array([0.0]),
        np.array([0.0]),
        np.array([M0]),
        delays[None],
        labeleff=LABELEFF,
        t1blood=T1BLOOD,
        is_casl=is_casl,
        tau=bolus,
        t1prime=T1PRIME[1],
    )
    np.testing.assert_allclose(predicted[0], deltam, rtol=1e-10)


def test_fit_tissue_curves_tissue_t1():
    """Local-regression tissue curves are fit with their own tissue's T1'."""
    from aslprep.utils.pvc import fit_tissue_curves

    is_casl, bolus, delays = PCASL
    curve = _simulate(is_casl, bolus, delays, [20.0], [1.3], T1PRIME[1:], [1.0])
    kwargs = {
        'labeleff': LABELEFF,
        't1blood': T1BLOOD,
        'is_casl': is_casl,
        'tau': bolus,
        'partition_coefficient': 0.9,
    }
    relaxed = fit_tissue_curves(
        curve[None],
        np.array([M0]),
        delays[None],
        t1prime=T1PRIME[1],
        **kwargs,
    )
    assert relaxed['cbf'][0] == pytest.approx(20.0, rel=1e-3)
    assert relaxed['att'][0] == pytest.approx(1.3, abs=1e-3)
    assert relaxed['kinres'][0] < 1e-6
    consensus = fit_tissue_curves(curve[None], np.array([M0]), delays[None], **kwargs)
    assert consensus['kinres'][0] > 10 * relaxed['kinres'][0]
