#!/usr/bin/env python3
"""Gate 1 for the coregistration probe's scorer: known answers must score as expected.

Plan: Task 2 ("Validation") and "Decision gates" 1. Three parts:

- ``synthetic``: the formulas of ``scoring`` on random oblique rigid transforms. The truth scores
  zero, a known error scores its size, a transform read in the inverse direction is caught, and a
  wrong composition order scores far from zero.
- ``images``: the image path on one real run's work directory (``run_real``). A direct
  registration of a reference to itself is the identity; one to a copy with a known moved header
  recovers the move, and every arm's inconsistency over that pair equals the direct registration's
  own error; the semi-synthetic runner (moved T1w header, production's chain and the ANTs arms
  rerun) returns each arm's own solution composed with the move; the carried-mask Dice drops when
  the solution is shifted or inverted.
- ``fixtures``: the gross-failure rules and the direct registration against aslscan truth. All
  development fixtures share one T1w, its tissue maps and one anatomical offset, so two fixtures
  form a pair of runs whose true direct registration is ``T_a T_b^-1`` (``T``: truth pulls).

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/validate_scoring.py --out <dir> synthetic
    python scripts/coreg_mc/validate_scoring.py --out <dir> images --work <run work dir> \\
        --license <license.txt>
    python scripts/coreg_mc/validate_scoring.py --out <dir> fixtures --fixtures <aslscan> \\
        --runs <fixture runs> --work <fixture work dirs> --names rc_pcasl ...
"""

import argparse
import itertools
import json
import shutil
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import run_ants  # noqa: E402
import run_real  # noqa: E402
import scoring  # noqa: E402
from score_full import _transform_error  # noqa: E402

SEED = 20261022
#: The known moves of the image checks: a few mm and degrees about the brain centre.
ASL_MOVE = {'translation': (3.0, -2.0, 1.5), 'rotation_deg': (3.0, -2.0, 4.0)}
T1W_MOVE = {'translation': (-4.0, 3.0, 2.0), 'rotation_deg': (-4.0, 3.0, 3.0)}
#: Sizes (RMS mm) of the deliberately wrong arms of the fixture checks, either side of the 5 mm
#: gross-failure limit.
PLANTED_MM = (2.0, 4.0, 6.0, 8.0, 12.0)
#: How close each direct-registration engine must come to a known answer on the real run. FLIRT
#: lands 0.2 mm from the identity even for an image registered to itself (ANTs 0.009 mm).
DIRECT_TOLERANCE_MM = {'ANTs': 0.05, 'FLIRT': 0.5}


class Checks:
    """Named pass/fail checks with their values, printed and written as JSON."""

    def __init__(self):
        self.items = []

    def add(self, name, ok, **values):
        self.items.append({'check': name, 'ok': bool(ok), **values})
        parts = [f'{k}={v:.4g}' if isinstance(v, float) else f'{k}={v}' for k, v in values.items()]
        shown = ', '.join(parts)
        print(f'{"PASS" if ok else "FAIL"} {name}: {shown}', flush=True)

    @property
    def passed(self):
        return all(i['ok'] for i in self.items)


def random_rigid(rng, mm, deg, centre=(0, 0, 0)):
    return scoring.rigid(rng.uniform(-mm, mm, 3), rng.uniform(-deg, deg, 3), centre)


def synthetic(checks):
    rng = np.random.default_rng(SEED)
    points = rng.normal(0, 40, (2000, 3))
    for i in range(5):
        truth = random_rigid(rng, 20, 20)
        checks.add(f'truth scores zero ({i})', scoring.distance(truth, truth, points) < 1e-9)
        shift = np.array([3.0, -4.0, 0.0])
        err = scoring.distance(truth, truth @ scoring.rigid(shift), points)
        checks.add(f'a 5 mm shift scores 5 mm ({i})', abs(err - 5.0) < 1e-9, rms_mm=err)
        inv = _transform_error(np.linalg.inv(truth), truth, points)
        checks.add(
            f'an inverted reading is caught ({i})',
            inv['rms_mm_if_inverted'] < 1e-9 and inv['rms_mm'] > 1.0,
            rms_mm=inv['rms_mm'],
        )

        c_b = random_rigid(rng, 20, 20)
        d_ab = random_rigid(rng, 5, 5)
        c_a = d_ab @ c_b
        checks.add(
            f'consistent pair scores zero ({i})',
            scoring.inconsistency(d_ab, c_a, c_b, points) < 1e-9,
        )
        p = random_rigid(rng, 3, 3)
        got = scoring.inconsistency(d_ab, d_ab @ p @ c_b, c_b, points)
        want = scoring.distance(np.eye(4), p, points)
        checks.add(f'known inconsistency recovered ({i})', abs(got - want) < 1e-9, rms_mm=got)
        wrong = scoring.distance(d_ab, c_b @ np.linalg.inv(c_a), points)
        checks.add(f'wrong composition order is far from zero ({i})', wrong > 1.0, rms_mm=wrong)
        floor = scoring.noise_floor(d_ab, np.linalg.inv(d_ab), points)
        checks.add(f'exact forward and reverse give a zero floor ({i})', floor < 1e-9)

        c = random_rigid(rng, 20, 20)
        m = random_rigid(rng, 10, 10)
        moved = points @ m[:3, :3].T + m[:3, 3]
        exact = scoring.semi_synthetic_error(c @ np.linalg.inv(m), c, m, moved)
        checks.add(f'semi-synthetic truth scores zero ({i})', exact < 1e-9)
        wrong = scoring.semi_synthetic_error(c @ m, c, m, moved)
        checks.add(
            f'semi-synthetic wrong convention is far from zero ({i})',
            wrong > 1.0,
            rms_mm=wrong,
        )
        q = random_rigid(rng, 3, 3)
        got = scoring.semi_synthetic_error(c @ np.linalg.inv(m) @ q, c, m, moved)
        want = scoring.distance(np.eye(4), np.linalg.inv(q), moved)
        checks.add(
            f'semi-synthetic known error recovered ({i})',
            abs(got - want) < 1e-9,
            rms_mm=got,
        )

    # The adopted gross-failure rule, branch by branch: (partners, mask flag) -> verdict.
    cases = [
        ([6.0, 7.0], False, 'fail'),
        ([6.0, 1.0], True, 'pass'),
        ([6.0], False, 'ambiguous'),
        ([1.0], True, 'pass'),
        ([], True, 'fail'),
        ([], False, 'pass'),
    ]
    got = [scoring.gross_verdict(values, mask) for values, mask, _ in cases]
    checks.add(
        'adopted gross rule: every partner, single partner ambiguous, mask without partners',
        got == [want for *_, want in cases],
        verdicts=got,
    )


def obliquity_deg(affine):
    """The largest angle between a voxel axis and its nearest world axis."""
    axes = affine[:3, :3] / np.sqrt((affine[:3, :3] ** 2).sum(0))
    return float(np.degrees(np.arccos(np.abs(axes).max(0))).max())


def images(checks, args):
    work = args.work
    stored = json.loads((work / 'transforms.json').read_text())
    arms = {k: np.asarray(v) for k, v in stored.items() if np.shape(v) == (4, 4)}
    asl = nb.load(work / 'conventional.nii.gz')
    t1 = nb.load(work / 't1w_brain.nii.gz')
    points = run_real.brain_points(work)
    checks.add(
        'headers',
        True,
        asl_obliquity_deg=obliquity_deg(asl.affine),
        t1w_obliquity_deg=obliquity_deg(t1.affine),
    )
    # The semi-synthetic check below reruns the harness's reproduction of production, so the
    # reproduction must first equal ASLPrep's stored transform on this run.
    r0_vs = scoring.distance(arms['R0 production'], arms['FSL conventional default'], points)
    checks.add(
        "the harness reproduces ASLPrep's stored transform",
        r0_vs < run_real.R0_REPRODUCTION_MM,
        rms_mm=r0_vs,
    )
    out = args.out / 'images'
    shutil.rmtree(out, ignore_errors=True)

    to_self = run_real.direct(work, work, out / 'self', 'conventional')
    for engine, d_self in zip(('ANTs', 'FLIRT'), to_self, strict=True):
        err = scoring.distance(d_self, np.eye(4), points)
        name = f'direct registration to itself is the identity ({engine})'
        checks.add(name, err < DIRECT_TOLERANCE_MM[engine], rms_mm=err)

    # Run b: run a's ASL images with the header M A (world moved by M); the pull b -> a is M^-1.
    centre = points.mean(0)
    m_asl = scoring.rigid(centre=centre, **ASL_MOVE)
    moved = out / 'moved_asl'
    moved.mkdir(parents=True)
    for name in run_real.ASL_SIDE:
        scoring.move_header(work / f'{name}.nii.gz', moved / f'{name}.nii.gz', m_asl)

    points_b = points @ m_asl[:3, :3].T + m_asl[:3, 3]
    start = scoring.distance(np.eye(4), np.linalg.inv(m_asl), points_b)
    for frame, image in run_real.FRAMES.items():
        d_ab, d_ab_flirt = run_real.direct(work, moved, out / f'{frame}_ab', image)
        d_ba, _ = run_real.direct(moved, work, out / f'{frame}_ba', image)
        err = scoring.distance(d_ab, np.linalg.inv(m_asl), points_b)
        err_flirt = scoring.distance(d_ab_flirt, np.linalg.inv(m_asl), points_b)
        floor = scoring.noise_floor(d_ab, d_ba, points)
        checks.add(
            f'direct registration recovers a known move ({frame})',
            err < DIRECT_TOLERANCE_MM['ANTs'] and err_flirt < DIRECT_TOLERANCE_MM['FLIRT'],
            ants_mm=err,
            flirt_mm=err_flirt,
            start_mm=start,
            closure_mm=floor,
        )
        spread = [scoring.inconsistency(d_ab, c, m_asl @ c, points_b) - err for c in arms.values()]
        checks.add(
            f'inconsistency equals the direct error for every arm ({frame})',
            np.abs(spread).max() < 1e-6,
            max_diff=float(np.abs(spread).max()),
        )

    # Semi-synthetic: the T1w side moved by M; every arm rerun from its own initialization.
    t1_points = scoring.mask_points(work / 't1w_mask.nii.gz', every=7)
    m_t1 = scoring.rigid(centre=t1_points.mean(0), **T1W_MOVE)
    t1_moved = t1_points @ m_t1[:3, :3].T + m_t1[:3, 3]
    found, fallbacks = run_real.semi_synthetic(
        work,
        m_t1,
        out / 'semi',
        list(run_ants.ARMS),
        args.license,
    )
    move_mm = scoring.distance(np.eye(4), m_t1, t1_points)
    rows = []
    for arm, c_moved in found.items():
        err = scoring.semi_synthetic_error(c_moved, arms[arm], m_t1, t1_moved)
        wrong = scoring.distance(c_moved, arms[arm] @ m_t1, t1_moved)
        unmoved = scoring.distance(c_moved, arms[arm], t1_moved)
        rows.append(
            {
                'arm': arm,
                'error_mm': err,
                'wrong_convention_mm': wrong,
                'as_if_unmoved_mm': unmoved,
                'fallback': arm in fallbacks,
            },
        )
    semi = pd.DataFrame(rows)
    semi.to_csv(args.out / 'images_semi_synthetic.tsv', sep='\t', index=False)
    print(semi.round(3).to_string(index=False))
    production = semi.set_index('arm').loc['FSL conventional default']
    checks.add(
        'semi-synthetic: production recovers the move',
        production['error_mm'] < 0.5 and production['wrong_convention_mm'] > move_mm,
        error_mm=float(production['error_mm']),
        move_mm=move_mm,
    )
    no_fallback = semi[~semi['fallback']]
    checks.add(
        'semi-synthetic: every arm is far closer to the right convention than the wrong one',
        bool((no_fallback['error_mm'] < no_fallback['wrong_convention_mm'] / 4).all()),
        worst_error_mm=float(no_fallback['error_mm'].max()),
    )

    carried = scoring.CarriedMask(work / 'mask.nii.gz', work / 't1w_mask.nii.gz')
    c = arms['R0 production']
    dice = carried.dice(c)
    shifted = carried.dice(c @ scoring.rigid(translation=(5, 0, 0)))
    inverted = carried.dice(np.linalg.inv(c))
    checks.add(
        'carried mask: a 5 mm shift and an inverted reading lower the Dice',
        shifted < dice and inverted < dice,
        dice=dice,
        shifted=shifted,
        inverted=inverted,
        calibration=carried.calibration(c),
    )


def planted(truth, points, rng):
    """Deliberately wrong arms: the truth composed with rigid errors of known RMS size.

    Each error is a random rotation axis and translation direction, scaled so the RMS
    displacement over ``points`` is the arm's size (``PLANTED_MM``); one draw per fixture.
    """
    centre = points.mean(0)
    out = {}
    for size in PLANTED_MM:
        error = scoring.scaled_rigid(rng.normal(size=6), size, points, centre)
        out[f'planted {size:g} mm'] = truth @ error
    return out


def fixtures(checks, args):
    rng = np.random.default_rng(SEED)
    data = {}
    for name in args.names:
        work = args.work / name
        fixture_dir = args.fixtures / name
        truth, points, production, voxel = run_ants.fixture_truth(fixture_dir, args.runs / name)
        arms, _fallbacks = run_real.collect(work, list(run_ants.ARMS))
        arms = {'R0 production': production} | arms
        missing = sorted(set(run_real.expected_arms(list(run_ants.ARMS))) - set(arms))
        checks.add(f'every arm has a transform ({name})', not missing, missing=missing)
        arms |= planted(truth, points, rng)
        errors = {arm: scoring.distance(c, truth, points) for arm, c in arms.items()}
        sidecar = next((fixture_dir / 'sub-01' / 'perf').glob('*_asl.json'))
        meta = json.loads(sidecar.read_text())
        data[name] = {
            'work': work,
            'truth': truth,
            'points': points,
            'arms': arms,
            'errors': errors,
            'voxel': voxel,
            'meta': meta,
        }

    t1 = {name: (args.fixtures / name / 'sub-01' / 'anat') for name in args.names}
    digests = {name: next(d.glob('*_T1w.nii.gz')).read_bytes() for name, d in t1.items()}
    distinct = len({hash(v) for v in digests.values()})
    checks.add('fixtures share one T1w', distinct == 1, n=len(digests))

    # Pairs of runs with known truth, among fixtures of the same contrast (2D readouts).
    group = [n for n in args.names if data[n]['meta'].get('MRAcquisitionType') == '2D']
    pair_rows = []
    cons_rows = []
    directs = {}
    for a, b in itertools.permutations(group, 2):
        da, db = data[a], data[b]
        d_true = da['truth'] @ np.linalg.inv(db['truth'])
        points_a = run_real.brain_points(da['work'])
        points_b = run_real.brain_points(db['work'])
        forward = args.out / 'fixtures' / f'{a}__to__{b}'
        reverse = args.out / 'fixtures' / f'{b}__to__{a}'
        d_ab, d_ab_flirt = run_real.direct(da['work'], db['work'], forward, 'conventional')
        d_ba, _ = run_real.direct(db['work'], da['work'], reverse, 'conventional')
        directs[(a, b)] = d_ab
        closure = scoring.noise_floor(d_ab, d_ba, points_a)
        engines = scoring.distance(d_ab, d_ab_flirt, points_b)
        pair_rows.append(
            {
                'run_a': a,
                'run_b': b,
                'direct_error_mm': scoring.distance(d_ab, d_true, points_b),
                'flirt_error_mm': scoring.distance(d_ab_flirt, d_true, points_b),
                'floor_closure_mm': closure,
                'floor_engines_mm': engines,
                'floor_mm': max(closure, engines),
            },
        )
        cons_rows += [
            {
                'dataset': 'fixtures',
                'run_a': a,
                'run_b': b,
                'frame': 'conventional',
                'arm': arm,
                'inconsistency_mm': scoring.inconsistency(d_ab, c_a, c_b, points_b),
                'with_true_direct_mm': scoring.inconsistency(d_true, c_a, c_b, points_b),
                'truth_error_a_mm': da['errors'][arm],
                'truth_error_b_mm': db['errors'][arm],
            }
            for arm, c_a in sorted(da['arms'].items())
            if (c_b := db['arms'].get(arm)) is not None
        ]
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(args.out / 'fixtures_direct.tsv', sep='\t', index=False)
    pd.DataFrame(cons_rows).to_csv(args.out / 'fixtures_consistency.tsv', sep='\t', index=False)
    print(pairs.round(3).to_string(index=False))
    # The closure alone is judged here: the combined floor also carries FLIRT's own imprecision,
    # which would make the check easy to pass.
    checks.add(
        'ANTs direct registration error is within twice its closure plus 0.2 mm',
        bool((pairs['direct_error_mm'] <= 2 * pairs['floor_closure_mm'] + 0.2).all()),
        worst_error_mm=float(pairs['direct_error_mm'].max()),
        median_error_mm=float(pairs['direct_error_mm'].median()),
        median_closure_mm=float(pairs['floor_closure_mm'].median()),
    )
    # The second engine must itself be sane (a mirrored FLIRT solution is 100+ mm off); how often
    # the combined floor covers the true error is reported, since a floor is an estimate, not a
    # bound.
    checks.add(
        'FLIRT direct registration is within 1 mm of the truth on every pair',
        bool((pairs['flirt_error_mm'] <= 1.0).all()),
        median_flirt_error_mm=float(pairs['flirt_error_mm'].median()),
        worst_flirt_error_mm=float(pairs['flirt_error_mm'].max()),
        median_floor_mm=float(pairs['floor_mm'].median()),
        median_engines_mm=float(pairs['floor_engines_mm'].median()),
        floor_covers_error=float((pairs['direct_error_mm'] <= pairs['floor_mm']).mean()),
    )

    rule_rows = mask_rule_rows(data)
    rule_rows += one_run_rows(data, group, directs, rng)
    rule_rows += every_run_rows(data, cons_rows)
    rule_rows += adopted_rows(data, cons_rows, rule_rows)
    rules = pd.DataFrame(rule_rows)
    rules['truth_gross'] = rules['truth_error_mm'] > scoring.GROSS_INCONSISTENCY_MM
    rules.to_csv(args.out / 'fixtures_gross.tsv', sep='\t', index=False)
    summary = rules.groupby(['scenario', 'rule']).apply(confusion, include_groups=False)
    summary.to_csv(args.out / 'fixtures_rules.tsv', sep='\t')
    print(summary.to_string())

    # How well the rules detect is for Task 3 to judge; these are the mechanics any version of
    # them must keep, so a regression in the rules fails gate 1.
    def flagged(scenario, rule, arm_pattern, run=None):
        sel = rules[(rules['scenario'] == scenario) & (rules['rule'] == rule)]
        sel = sel[sel['arm'].str.contains(arm_pattern, regex=False)]
        if run is not None:
            sel = sel[sel['fixture'] == run]
        return sel['flag']

    for kind in ('axes', 'rigid'):
        caught = flagged('all arms', f'mask ({kind})', 'planted 12 mm')
        name = f'mask rule ({kind}) flags every 12 mm error'
        checks.add(name, bool(caught.all()), n=len(caught))
        clean = flagged('all arms', f'mask ({kind})', 'planted 2 mm')
        checks.add(f'mask rule ({kind}) never flags a 2 mm error', not clean.any(), n=len(clean))

    one = rules[rules['scenario'] == 'one run wrong']
    one = one[one['rule'] == 'inconsistency (every partner)']
    target = one.apply(lambda r: r['arm'].startswith(f'{r["fixture"]} '), axis=1)
    size = one['arm'].str.extract(r'planted (\d+) mm', expand=False).astype(float)
    big = size > scoring.GROSS_INCONSISTENCY_MM
    checks.add(
        'every-partner rule: a lone run wrong by 6 mm or more is flagged',
        bool(one.loc[target & big, 'flag'].all()),
        n=int((target & big).sum()),
    )
    checks.add(
        'every-partner rule: good runs beside a wrong one are never flagged',
        not one.loc[~target, 'flag'].any(),
        n=int((~target).sum()),
    )


def confusion(df):
    """Flags against the truth (error above the gross limit), by run."""
    return pd.Series(
        {
            'caught': int((df['flag'] & df['truth_gross']).sum()),
            'missed': int((~df['flag'] & df['truth_gross']).sum()),
            'false_alarms': int((df['flag'] & ~df['truth_gross']).sum()),
            'ambiguous': int((df.get('verdict') == 'ambiguous').sum()),
            'clean': int((~df['flag'] & ~df['truth_gross']).sum()),
            'missed_arms': ', '.join(sorted(df.loc[~df['flag'] & df['truth_gross'], 'arm'])),
        },
    )


def mask_rule_rows(data):
    """The carried-mask rule on every fixture and arm, with both calibrations."""
    rows = []
    for name, d in data.items():
        carried = scoring.CarriedMask(d['work'] / 'mask.nii.gz', d['work'] / 't1w_mask.nii.gz')
        for kind in ('axes', 'rigid'):
            rows += [
                {
                    'scenario': 'all arms',
                    'rule': f'mask ({kind})',
                    'fixture': name,
                    'arm': arm,
                    'truth_error_mm': d['errors'][arm],
                    'flag': result['gross_mask'],
                    'dice': result['dice'],
                    'best_dice': result['best_dice'],
                    'limit': result['limit'],
                }
                for arm, result in scoring.gross_by_mask(carried, d['arms'], kind=kind).items()
            ]
    return rows


def flags(values, limit=scoring.GROSS_INCONSISTENCY_MM):
    return {
        'inconsistency (any partner)': any(v > limit for v in values),
        'inconsistency (every partner)': bool(values) and all(v > limit for v in values),
    }


def one_run_rows(data, group, directs, rng, base='ANTs Q1'):
    """One run wrong at a time: the target carries the truth with a planted error.

    The other runs carry their own ``base`` solution; the target and its partners are scored.
    """
    asl_points = {n: run_real.brain_points(data[n]['work']) for n in group}

    def pair_value(pulls, run, partner):
        """The larger inconsistency of the two ordered pairs."""
        both = [
            scoring.inconsistency(directs[(a, b)], pulls[a], pulls[b], asl_points[b])
            for a, b in ((run, partner), (partner, run))
        ]
        return max(both)

    rows = []
    for target in group:
        d = data[target]
        centre = d['points'].mean(0)
        for size in PLANTED_MM:
            error = scoring.scaled_rigid(rng.normal(size=6), size, d['points'], centre)
            pulls = {n: data[n]['arms'][base] for n in group} | {target: d['truth'] @ error}
            truth_err = {n: data[n]['errors'][base] for n in group}
            truth_err[target] = scoring.distance(pulls[target], d['truth'], d['points'])
            for run in group:
                values = [pair_value(pulls, run, partner) for partner in group if partner != run]
                base_row = {
                    'scenario': 'one run wrong',
                    'fixture': run,
                    'arm': f'{target} planted {size:g} mm',
                    'truth_error_mm': truth_err[run],
                }
                rows += [base_row | {'rule': r, 'flag': f} for r, f in flags(values).items()]
    return rows


def every_run_rows(data, cons_rows):
    """Every run of an arm wrong at once (the planted arms on all fixtures, and the real arms)."""
    rows = []
    for name in sorted({r['run_a'] for r in cons_rows}):
        for arm in data[name]['arms']:
            partners = run_real.partner_inconsistency(cons_rows, 'fixtures', name, arm)
            base_row = {
                'scenario': 'every run wrong' if arm.startswith('planted') else 'real arms',
                'fixture': name,
                'arm': arm,
                'truth_error_mm': data[name]['errors'][arm],
            }
            rule_flags = flags(list(partners.values()))
            rows += [base_row | {'rule': r, 'flag': f} for r, f in rule_flags.items()]
    return rows


def adopted_rows(data, cons_rows, rule_rows):
    """The adopted rule (``scoring.gross_verdict``) on every fixture and arm, every run at once.

    Fixtures without a partner fall back on the rigid-calibration mask flags in ``rule_rows``.
    """
    mask = {(r['fixture'], r['arm']): r['flag'] for r in rule_rows if r['rule'] == 'mask (rigid)'}
    rows = []
    for name, d in data.items():
        for arm in d['arms']:
            partners = run_real.partner_inconsistency(cons_rows, 'fixtures', name, arm)
            verdict = scoring.gross_verdict(partners.values(), mask[(name, arm)])
            rows.append(
                {
                    'scenario': 'every run wrong' if arm.startswith('planted') else 'real arms',
                    'rule': 'adopted',
                    'fixture': name,
                    'arm': arm,
                    'truth_error_mm': d['errors'][arm],
                    'flag': verdict != 'pass',
                    'verdict': verdict,
                },
            )
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--out', type=Path, required=True)
    sub = p.add_subparsers(dest='part', required=True)
    sub.add_parser('synthetic')
    im = sub.add_parser('images')
    im.add_argument('--work', type=Path, required=True)
    im.add_argument('--license', type=Path, required=True)
    fx = sub.add_parser('fixtures')
    fx.add_argument('--fixtures', type=Path, required=True)
    fx.add_argument('--runs', type=Path, required=True)
    fx.add_argument('--work', type=Path, required=True)
    fx.add_argument('--names', nargs='+', required=True)
    args = p.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    checks = Checks()
    if args.part == 'synthetic':
        synthetic(checks)
    elif args.part == 'images':
        images(checks, args)
    else:
        fixtures(checks, args)
    report = json.dumps(checks.items, indent=1, default=float)
    (args.out / f'{args.part}_checks.json').write_text(report + '\n')
    print('ALL PASS' if checks.passed else 'SOME FAILED')
    return 0 if checks.passed else 1


if __name__ == '__main__':
    sys.exit(main())
