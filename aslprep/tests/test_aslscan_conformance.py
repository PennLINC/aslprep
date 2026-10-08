"""Fast conformance tier: ASLPrep's CBF interfaces on raw simulated data (spec Section 4.7).

``ExtractCBF`` and ``ComputeCBF`` run directly on each fast fixture, with no preprocessing,
and their CBF is compared with an independent implementation of the documented model
(:mod:`aslprep.tests.truth_models`). Mutated expectations must fail the same bounds, which
shows the comparison is sensitive. Agreement with the physical truth is recorded but not
asserted here (Tier B), except for the known background-suppression gap.
"""

import json

import nibabel as nb
import numpy as np
import pytest

from aslprep.tests import aslscan_fixtures as af
from aslprep.tests import truth_models as tm

#: Tier A bounds for the exact (unpreprocessed) comparison (spec Section 4.7).
MEDIAN_TOL = 1e-3
P99_TOL = 1e-2
MIN_VOXELS = 500

SINGLE_DELAY = sorted(n for n in af.RECIPES if n.startswith('fast_') and 'multipld' not in n)
MULTI_DELAY = sorted(n for n in af.RECIPES if n.startswith('fast_') and 'multipld' in n)
#: Multi-delay Tier A bounds (spec Section 4.7) and the size of the refitted voxel subset.
MULTI_CBF_TOL = 0.02
MULTI_ATT_TOL = 0.05  # s
MULTI_FINITE = 0.98
MULTI_SUBSET = 300
MUTATIONS = ('swap', 'reverse_slices', 'pld_shift', 'm0_scale', 'efficiency_ignored')

#: Known disagreements with the documented model, as strict expected failures (spec Section 8).
KNOWN_TIER_A = {
    'fast_pasl_q2tips': (
        'ASLPrep quantifies single-delay Q2TIPS with exp(TI2 / T1b) and no slice-time shift '
        '(aslprep/interfaces/cbf.py, ComputeCBF); Alsop 2015 eq. 2 uses the inversion time TI. '
        'Here CBF is 0.81 of the white-paper value, and the white paper matches the truth as '
        'well as QUIPSS II does. Reported; pending a decision.'
    ),
}


def _known(names):
    return [
        pytest.param(n, marks=pytest.mark.xfail(strict=True, reason=KNOWN_TIER_A[n]))
        if n in KNOWN_TIER_A
        else n
        for n in names
    ]


def _fixture(name, data_dir):
    try:
        return af.fixture_dir(name, data_dir)
    except af.FixtureUnavailable as exc:
        pytest.skip(str(exc))


class Run:
    """One fast fixture's raw inputs, ASLPrep's CBF from them, and the truth."""

    def __init__(self, name, fixture, workdir, fwhm=0.0):
        from aslprep.interfaces.cbf import ComputeCBF, ExtractCBF

        self.recipe = af.load_recipe(name)
        perf = fixture / 'sub-01' / 'perf'
        stem = f'sub-01_acq-{self.recipe.acq}'
        self.asl_file = perf / f'{stem}_asl.nii.gz'
        self.metadata = json.loads((perf / f'{stem}_asl.json').read_text())
        self.context = (perf / f'{stem}_aslcontext.tsv').read_text().split()[1:]
        m0_file = perf / f'{stem}_m0scan.nii.gz'
        self.m0_file = m0_file if m0_file.exists() else None
        self.m0_metadata = (
            json.loads((perf / f'{stem}_m0scan.json').read_text()) if self.m0_file else None
        )
        self.m0_scale = float(self.recipe.m0_divisor)
        self.fwhm = fwhm

        extract = ExtractCBF(
            name_source=str(self.asl_file),
            asl_file=str(self.asl_file),
            metadata=self.metadata,
            aslcontext=str(perf / f'{stem}_aslcontext.tsv'),
            m0scan=str(self.m0_file) if self.m0_file else None,
            m0scan_metadata=self.m0_metadata,
            fwhm=fwhm,
        ).run(cwd=str(workdir))
        compute = ComputeCBF(
            deltam=extract.outputs.out_file,
            metadata=extract.outputs.metadata,
            m0_scale=self.m0_scale,
            m0_file=extract.outputs.m0_file,
            cbf_only=False,
        )
        if extract.outputs.m0tr is not None:
            compute.inputs.m0tr = extract.outputs.m0tr
        result = compute.run(cwd=str(workdir))
        self.cbf = nb.load(result.outputs.mean_cbf).get_fdata()
        att = result.outputs.att
        self.att = nb.load(att).get_fdata() if att else None

        gt = perf / 'ground-truth'
        self.simulation = json.loads((gt / 'simulation.json').read_text())
        self.perfusion = nb.load(gt / 'sub-01_desc-perfusion_gt.nii.gz').get_fdata()
        self.pv = {
            t: nb.load(gt / f'sub-01_desc-pv{t}_gt.nii.gz').get_fdata()
            for t in ('GM', 'WM', 'CSF')
        }
        img = nb.load(self.asl_file)
        self.asl, self.affine = img.get_fdata(), img.affine
        self.m0scan = nb.load(self.m0_file).get_fdata() if self.m0_file else None

    def expected(self, mutation=None, alpha=None):
        return tm.expected_cbf(
            self.asl,
            self.context,
            self.metadata,
            self.affine,
            m0scan=self.m0scan,
            m0_metadata=self.m0_metadata,
            m0_scale=self.m0_scale,
            fwhm=self.fwhm,
            alpha=alpha,
            mutation=mutation,
        )

    def valid(self, reference):
        pv = self.pv
        brain = pv['GM'] + pv['WM'] + pv['CSF'] >= 0.5
        mask = brain & (pv['CSF'] < 0.5) & np.isfinite(reference) & (reference > 5)
        assert mask.sum() >= MIN_VOXELS, f'only {mask.sum()} valid voxels'
        return mask

    def deviation(self, reference, mask=None):
        """(|median ratio - 1|, p99 |ratio - 1|) of ASLPrep's CBF against ``reference``."""
        mask = self.valid(reference) if mask is None else mask
        ratio = self.cbf[mask] / reference[mask]
        return abs(np.median(ratio) - 1), np.percentile(np.abs(ratio - 1), 99)

    def tier_b(self):
        """Median of ASLPrep CBF / true perfusion in GM- and WM-dominant voxels."""
        out = {}
        for tissue in ('GM', 'WM'):
            mask = self.valid(self.perfusion) & (self.pv[tissue] >= 0.7)
            # nanmedian: a voxel whose multi-delay fit failed is NaN in ASLPrep's output
            out[tissue] = float(np.nanmedian(self.cbf[mask] / self.perfusion[mask]))
        return out


@pytest.fixture(scope='module')
def runs(data_dir, tmp_path_factory):
    cache = {}

    def get(name, fwhm=0.0):
        key = (name, fwhm)
        if key not in cache:
            fixture = _fixture(name, data_dir)
            cache[key] = Run(name, fixture, tmp_path_factory.mktemp(f'{name}-{fwhm}'), fwhm)
        return cache[key]

    return get


@pytest.mark.parametrize('name', _known(SINGLE_DELAY))
def test_tier_a_exact(name, runs, record_property):
    run = runs(name)
    median, p99 = run.deviation(run.expected())
    record_property('tier_a_median_deviation', median)
    record_property('tier_b', run.tier_b())
    assert median <= MEDIAN_TOL, f'|median(cbf / cbf_A) - 1| = {median:.4g} > {MEDIAN_TOL}'
    assert p99 <= P99_TOL, f'p99 |cbf / cbf_A - 1| = {p99:.4g} > {P99_TOL}'


@pytest.mark.parametrize('mutation', MUTATIONS)
@pytest.mark.parametrize('name', SINGLE_DELAY)
def test_tier_a_detects_mutations(name, mutation, runs):
    """Each deliberate error in the expectation must break the Tier A bound."""
    run = runs(name)
    if mutation == 'reverse_slices' and len(set(run.metadata.get('SliceTiming', [0]))) < 2:
        pytest.skip('no slice timing to reverse')
    if mutation == 'swap' and 'label' not in run.context:
        pytest.skip('delta-M volumes have no label-control order to swap')
    # score on the unmutated expectation's mask (a swapped sign would otherwise empty it)
    mask = run.valid(run.expected())
    median, p99 = run.deviation(run.expected(mutation=mutation), mask)
    assert median > MEDIAN_TOL or p99 > P99_TOL, f'{mutation} went undetected'


@pytest.mark.parametrize('name', ['fast_pcasl_seq', 'fast_m0_included', 'fast_m0_absent'])
def test_tier_a_with_m0_smoothing(name, runs):
    """ASLPrep's default 5 mm M0 smoothing, reproduced with the same nibabel call."""
    run = runs(name, fwhm=5.0)
    median, p99 = run.deviation(run.expected())
    assert median <= MEDIAN_TOL
    assert p99 <= P99_TOL


def test_background_suppression_gap(runs):
    """Known gap (spec Section 8): ASLPrep assumes 0.95 per suppression pulse.

    The simulator attenuates the label by |1 - 2 x 0.95| = 0.9 per pulse, so with two pulses
    ASLPrep's CBF is 0.9^2 / 0.95^2 = 0.8975 of the physically calibrated value. The size is
    asserted first, so an unexpected change fails instead of passing as an expected failure.
    """
    run = runs('fast_bs_le_absent')
    physical = run.expected(alpha=tm.labeling_efficiency_physical(run.simulation))
    mask = run.valid(physical)
    measured = float(np.median(run.cbf[mask] / physical[mask]))
    predicted = (0.9 / 0.95) ** 2
    assert measured == pytest.approx(predicted, rel=0.02), (
        f'the suppression gap changed: measured {measured:.4f}, predicted {predicted:.4f}'
    )
    pytest.xfail(
        f'ASLPrep assumes 0.95 per suppression pulse; the simulator applies 0.9 '
        f'(CBF ratio {measured:.4f}, predicted {predicted:.4f})'
    )


def test_q2tips_convention_is_as_diagnosed(runs):
    """Pin the cause of the known Q2TIPS disagreement, so a partial change cannot hide.

    ASLPrep's current single-delay Q2TIPS CBF equals the white-paper form with ``TI2`` in the
    exponent instead of the slice-shifted inversion time.
    """
    run = runs('fast_pasl_q2tips')
    md = run.metadata
    ti1, ti2 = md['BolusCutOffDelayTime']
    t1b = tm.T1_BLOOD[md['MagneticFieldStrength']]
    deltam = tm.deltam_mean(run.asl, run.context)
    convention = (
        tm.UNIT_CONVERSION
        * tm.LAMBDA
        * deltam
        * np.exp(ti2 / t1b)
        / (2 * tm.labeling_efficiency(md) * ti1 * run.m0scan)
    )
    median, p99 = run.deviation(convention)
    assert median <= MEDIAN_TOL
    assert p99 <= P99_TOL


def _multi_subset(run):
    """A deterministic subset of truth-valid voxels for the reference fits."""
    pv = run.pv
    mask = (pv['GM'] + pv['WM'] + pv['CSF'] >= 0.5) & (pv['CSF'] < 0.5) & (run.perfusion > 5)
    voxels = np.argwhere(mask)
    step = max(1, len(voxels) // MULTI_SUBSET)
    return voxels[::step][:MULTI_SUBSET], mask


def _multi_deviation(run, expected, voxels, reference):
    """(|median CBF ratio - 1|, median |ATT difference|) over the subset.

    Voxels are chosen from the *unmutated* ``reference`` fit (finite, CBF > 5) where ASLPrep's
    own fit succeeded; its failures are counted separately by the finite-fraction check. A
    mutated expectation that fails or collapses there counts as a disagreement (infinite).
    """
    i, j, k = voxels.T
    cbf, att = run.cbf[i, j, k], run.att[i, j, k]
    ok = np.isfinite(reference[:, 0]) & (reference[:, 0] > 5) & np.isfinite(cbf)
    assert ok.sum() >= MULTI_SUBSET // 2, f'only {ok.sum()} usable voxels in the subset'
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = cbf[ok] / expected[ok, 0]
    ratio[~np.isfinite(ratio)] = np.inf
    att_diff = np.abs(att[ok] - expected[ok, 1])
    att_diff[~np.isfinite(att_diff)] = np.inf
    return abs(np.median(ratio) - 1), np.median(att_diff)


def _multi_expected(run, voxels, mutation=None):
    return tm.expected_multi_delay(
        run.asl,
        run.context,
        run.metadata,
        run.affine,
        voxels,
        m0scan=run.m0scan,
        m0_metadata=run.m0_metadata,
        m0_scale=run.m0_scale,
        mutation=mutation,
    )


_REFERENCE_FITS = {}


def _reference_fit(name, run, voxels):
    if name not in _REFERENCE_FITS:
        _REFERENCE_FITS[name] = _multi_expected(run, voxels)
    return _REFERENCE_FITS[name]


@pytest.mark.parametrize('name', MULTI_DELAY)
def test_tier_a_multi_delay(name, runs, record_property):
    run = runs(name)
    voxels, mask = _multi_subset(run)
    finite = np.isfinite(run.cbf[mask]).mean()
    reference = _reference_fit(name, run, voxels)
    cbf_dev, att_dev = _multi_deviation(run, reference, voxels, reference)
    record_property('tier_a', {'cbf': cbf_dev, 'att': att_dev, 'finite': finite})
    record_property('tier_b', run.tier_b())
    assert finite >= MULTI_FINITE, f'finite fraction {finite:.3f} < {MULTI_FINITE}'
    assert cbf_dev <= MULTI_CBF_TOL, f'|median(cbf / cbf_A) - 1| = {cbf_dev:.4g}'
    assert att_dev <= MULTI_ATT_TOL, f'median |att - att_A| = {att_dev:.4g} s'


@pytest.mark.parametrize('mutation', ['swap', 'pld_shift', 'm0_scale'])
@pytest.mark.parametrize('name', MULTI_DELAY)
def test_tier_a_multi_delay_detects_mutations(name, mutation, runs):
    run = runs(name)
    voxels, _ = _multi_subset(run)
    reference = _reference_fit(name, run, voxels)
    mutated = _multi_expected(run, voxels, mutation)
    cbf_dev, att_dev = _multi_deviation(run, mutated, voxels, reference)
    assert cbf_dev > MULTI_CBF_TOL or att_dev > MULTI_ATT_TOL, f'{mutation} went undetected'
