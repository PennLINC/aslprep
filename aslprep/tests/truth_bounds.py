"""Acceptance ceilings for the truth-scored integration tests (spec Section 4.7).

Ceilings are fixed a priori from what each method should achieve; they are never fitted to
ASLPrep's output. A ceiling may only be loosened in a commit that explains why the method
cannot meet it. Every assertion goes through :func:`check`, which fails when the metric or
its ceiling is missing, so a test cannot pass silently.

Regression bands (narrower, measured) are a separate, optional layer: :data:`BANDS`, filled
only after the ceilings pass over at least five runs (see the developer documentation).
``python -m aslprep.tests.truth_bounds propose truth_score.json ...`` prints proposals.
"""

import math
import sys

from aslprep.tests.truth_scoring import lookup

#: name -> (score path, lower, upper, unit, rationale). ``None`` leaves a side open.
#: Paths may contain ``{t}`` (tissue) or ``{space}``, filled in by the caller.
CEILINGS = {
    'mask_coverage': (
        ('native', 'coverage', 'mask_coverage'),
        0.95,
        None,
        'fraction',
        "ASLPrep's brain mask must cover the brain",
    ),
    'finite': (
        ('native', 'coverage', 'finite'),
        0.999,
        None,
        'fraction',
        'single-delay CBF is closed-form; every brain voxel is finite',
    ),
    'finite_multi': (
        ('native', 'coverage', 'finite'),
        0.98,
        None,
        'fraction',
        'the multi-delay fit may fail in a few voxels',
    ),
    'n_dominant': (
        ('native', 'coverage', 'n_{t}'),
        300,
        None,
        'voxels',
        'enough tissue-dominant voxels for stable medians',
    ),
    # End to end (raw data through the documented model): the median is unbiased by
    # preprocessing. Per-voxel tails are report-only: motion correction resamples label and
    # control volumes separately, and delta-M (about 1 % of the signal) amplifies interpolation
    # differences at tissue boundaries (15-20 % at the 95th percentile on motion-free data).
    'tier_a_median': (
        ('native', 'tier_a', 'median_dev'),
        None,
        0.05,
        'ratio',
        'preprocessing (HMC, M0 registration) leaves the median within 5 % of the model',
    ),
    # Quantification given ASLPrep's own preprocessed series, in voxels 10 mm inside its mask.
    'tier_a_quant_median': (
        ('native', 'tier_a_quant', 'median_dev'),
        None,
        0.01,
        'ratio',
        'the CBF model and calibration match the documented model given the same series',
    ),
    'tier_a_quant_p95': (
        ('native', 'tier_a_quant', 'p95_abs_dev'),
        None,
        0.05,
        'ratio',
        'only M0 registration and interpolation separate the two',
    ),
    'tier_a_att': (
        ('native', 'tier_a_quant', 'att_median_abs_diff'),
        None,
        0.1,
        's',
        'ATT from the same model fitted to the same preprocessed series',
    ),
    'aslref_pose': (
        ('frames', 'aslref_pose_rms_mm'),
        None,
        0.5,
        'mm',
        'without motion, native outputs compare voxel to voxel only if the aslref is within a '
        'seventh of a voxel of the static frame',
    ),
    'tier_b': (
        ('native', 'tier_b', '{t}', 'median'),
        None,
        0.05,
        'ratio',
        'within the Tier A ceiling of the independent expectation '
        "(checked against ('native', 'expected_ratio', t))",
    ),
    # Coregistration accuracy scales with resolution, so it is bounded in units of the coarsest
    # acquisition voxel (1 mm / 1 degree absolute ceilings, as in qsiprep, assumed 2-3 mm data;
    # the 8 mm slices of a GE 3D spiral cannot be held to them).
    'coreg_rot': (
        ('frames', 'coreg', 'rot_arc_voxels'),
        None,
        0.25,
        'voxels',
        'rigid coregistration under --sloppy: the rotation error, as arc length at 70 mm, '
        'within a quarter of the coarsest voxel',
    ),
    'coreg_rms': (
        ('frames', 'coreg', 'rms_voxels'),
        None,
        0.25,
        'voxels',
        'rigid coregistration under --sloppy: RMS displacement over brain voxels within a '
        'quarter of the coarsest voxel',
    ),
    'motion_median': (
        ('frames', 'motion', 'rms_error_median_mm'),
        None,
        0.5,
        'mm',
        'MCFLIRT recovers rigid motion of a few mm to sub-voxel accuracy',
    ),
    'motion_max': (
        ('frames', 'motion', 'rms_error_max_mm'),
        None,
        1.5,
        'mm',
        'worst volume, including spikes',
    ),
    'space': (
        ('{space}', '{t}', 'median'),
        None,
        0.15,
        'ratio',
        'GM only: two resamplings at ASL resolution blur tissue boundaries; within 0.15 of the '
        'native expectation. WM is not bounded: GM signal spilling into it dominates its ratio '
        '(F1: 0.56 against a native 0.39), so alignment is checked by space_alignment instead',
    ),
    'space_finite': (
        ('{space}', '{t}', 'finite'),
        0.99,
        None,
        'fraction',
        'standard-space CBF covers the truth voxels',
    ),
    'space_alignment': (
        ('{space}', 'r_native'),
        0.8,
        None,
        'r',
        "correlation with ASLPrep's native CBF at the same anatomical points; a 4 mm "
        'misregistration drops it to about 0.55 on these phantoms (r_native_shifted_4mm)',
    ),
    # Pairwise ceilings: used with check_pair, which supplies both paths.
    'scorescrub': (
        None,
        None,
        0.05,
        'ratio',
        "without outliers, SCORE and SCRUB stay within 5 % of the mean CBF's accuracy",
    ),
}

#: (recipe, name) -> replacement (lower, upper, rationale). Empty: no recipe is exempt.
OVERRIDES = {}

#: (recipe, name) -> (lower, upper, comment). Measured regression bands; see module docstring.
BANDS = {}


def _limits(name, recipe):
    if name not in CEILINGS:
        raise AssertionError(f'no ceiling named {name!r}: add one to truth_bounds.CEILINGS')
    path, lo, hi, unit, why = CEILINGS[name]
    if (recipe, name) in OVERRIDES:
        lo, hi, why = OVERRIDES[(recipe, name)]
    return path, lo, hi, unit, why


def _fail_message(label, value, lo, hi, unit, why):
    bounds = ' and '.join(
        part for part in (lo is not None and f'>= {lo}', hi is not None and f'<= {hi}') if part
    )
    return f'{label} = {value:.4g} {unit}, expected {bounds} ({why})'


def _check_value(label, value, lo, hi, unit, why):
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        raise AssertionError(f'{label} is missing or not finite ({value!r})')
    if (lo is not None and value < lo) or (hi is not None and value > hi):
        raise AssertionError(_fail_message(label, value, lo, hi, unit, why))


def check(score, name, recipe, reference=None, **fields):
    """Assert one ceiling. ``fields`` fill the path's ``{t}``/``{space}`` placeholders.

    With ``reference`` (a score path), the metric is compared as ``|value - reference|``.
    """
    path, lo, hi, unit, why = _limits(name, recipe)
    if path is None:
        raise AssertionError(f'{name!r} is a pairwise ceiling; use check_pair')
    path = tuple(p.format(**fields) for p in path)
    value = lookup(score, path)
    label = '.'.join(path)
    if reference is not None:
        ref_path = tuple(p.format(**fields) for p in reference)
        ref = lookup(score, ref_path)
        if ref is None or value is None:
            raise AssertionError(f'{label} or {".".join(ref_path)} is missing')
        value, label = abs(value - ref), f'|{label} - {".".join(ref_path)}|'
    _check_value(label, value, lo, hi, unit, why)
    band = BANDS.get((recipe, name))
    if band is not None:
        _check_value(f'{label} (regression band)', value, band[0], band[1], unit, band[2])


def check_pair(score, name, recipe, path, reference, **fields):
    """Assert ``|score[path] - score[reference]|`` against the ceiling ``name``."""
    _, lo, hi, unit, why = _limits(name, recipe)
    path = tuple(p.format(**fields) for p in path)
    reference = tuple(p.format(**fields) for p in reference)
    value, ref = lookup(score, path), lookup(score, reference)
    label = f'|{".".join(path)} - {".".join(reference)}|'
    if value is None or ref is None:
        raise AssertionError(f'{label}: a value is missing')
    _check_value(label, abs(value - ref), lo, hi, unit, why)


def propose(paths):
    """Print each numeric metric's values across score files (a developer aid)."""
    import json

    scores = [json.loads(open(p).read()) for p in paths]

    def walk(node, prefix=()):
        if isinstance(node, dict):
            for key, value in node.items():
                yield from walk(value, (*prefix, key))
        elif isinstance(node, int | float) and not isinstance(node, bool):
            yield prefix, node

    values = {}
    for score in scores:
        for path, value in walk(score):
            values.setdefault(path, []).append(value)
    for path, vals in sorted(values.items()):
        spread = max(vals) - min(vals)
        print(
            f'{".".join(map(str, path))}: n={len(vals)} min={min(vals):.4g} '
            f'max={max(vals):.4g} spread={spread:.3g}'
        )


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == 'propose':
        propose(sys.argv[2:])
    else:
        print('usage: python -m aslprep.tests.truth_bounds propose truth_score.json ...')
