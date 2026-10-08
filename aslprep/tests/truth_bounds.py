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
    'tier_a_median': (
        ('native', 'tier_a', 'median_abs_dev'),
        None,
        0.05,
        'ratio',
        'preprocessing (HMC, M0 registration) perturbs the exact model by < 5 %',
    ),
    'tier_a_p95': (
        ('native', 'tier_a', 'p95_abs_dev'),
        None,
        0.15,
        'ratio',
        'edge voxels move under resampling; 95 % stay within 15 %',
    ),
    'tier_a_att': (
        ('native', 'tier_a', 'att_median_abs_diff'),
        None,
        0.15,
        's',
        'ATT from the same fit on preprocessed data',
    ),
    'tier_b': (
        ('native', 'tier_b', '{t}', 'median'),
        None,
        0.05,
        'ratio',
        'within the Tier A ceiling of the independent expectation '
        "(checked against ('native', 'expected_ratio', t))",
    ),
    'coreg_rot': (
        ('frames', 'coreg', 'rot_deg'),
        None,
        1.0,
        'deg',
        'rigid coregistration under --sloppy (as qsiprep)',
    ),
    'coreg_rms': (
        ('frames', 'coreg', 'rms_mm'),
        None,
        1.0,
        'mm',
        'rigid coregistration under --sloppy, RMS over brain voxels',
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
        'two resamplings at ASL resolution blur tissue boundaries; within 0.15 of '
        'the native expectation',
    ),
    'space_finite': (
        ('{space}', '{t}', 'finite'),
        0.99,
        None,
        'fraction',
        'standard-space CBF covers the truth voxels',
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
