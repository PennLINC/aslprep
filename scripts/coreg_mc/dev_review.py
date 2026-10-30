#!/usr/bin/env python3
"""The development review batch of Task 3 (plan: "Blinded review figures", Task 3).

Two comparisons on development real runs, outside the runs every arm fails:

- ``R1 vs R2``: R1 (the conventional reference, BBR sign by the frozen contrast rule) against R2
  (P-mean, default sign), on runs where R2 is eligible;
- ``sign``: the conventional reference with the default sign against the opposite sign, which
  checks R1's rule on real data.

A pair is a candidate when its two transforms differ by more than the run's separation: the
larger of ``MIN_SEPARATION_MM`` (the smallest known-answer size) and the run's median noise floor
over its partners (its family's median when it has none). The tie distance itself is frozen after
this batch, from the reviewer's known-answer sensitivity. Pairs are drawn by seed, round-robin over
families, at most one per subject per comparison. Known-answer pairs are R1's solution against
itself moved by 1, 2 and 4 mm, on runs drawn the same way.

Writes ``<out>/batch.json`` and ``<out>/figures`` (what the review page loads), and
``<out>/key.json`` and ``<out>/selection.tsv`` (never published with the page).

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/dev_review.py --freeze <freeze dir> --dev-real <dir> --out <dir>
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'hmc_probe'))

import review  # noqa: E402
import run_real  # noqa: E402
import scoring  # noqa: E402

SEED = 20261030
PER_COMPARISON = 12
KNOWN_PER_SIZE = 5
MIN_SEPARATION_MM = min(review.KNOWN_MM)
R2_ARM = 'FSL P-mean default'


def r1_arm(contrast):
    """R1's frozen rule: the opposite sign when white matter is brighter than grey."""
    return f'FSL conventional {"opposite" if contrast > 0 else "default"}'


def floors(dev_real):
    """Each run's median noise floor over its partners (conventional frame)."""
    cons = run_real_tables(dev_real, 'consistency')
    cons = cons[(cons['frame'] == 'conventional') & (cons['arm'] == 'R0 production')]
    keep = ['dataset', 'run', 'floor_mm']
    both = pd.concat(
        [cons.rename(columns={side: 'run'})[keep] for side in ('run_a', 'run_b')],
        ignore_index=True,
    )
    return both.groupby(['dataset', 'run'])['floor_mm'].median()


def run_real_tables(folder, suffix):
    tables = []
    for f in sorted(folder.glob(f'*_{suffix}.tsv')):
        try:
            tables.append(pd.read_table(f))
        except pd.errors.EmptyDataError:
            continue

    return pd.concat(tables, ignore_index=True)


def candidates(real, dev_real):
    """One row per run and comparison, with the two transforms' distance and the separation."""
    runs = real.drop_duplicates(['dataset', 'run'])
    runs = runs[~runs['shared_failure']]
    floor = floors(dev_real)
    keys = zip(runs['dataset'], runs['run'], strict=True)
    runs = runs.assign(floor_mm=[floor.get(k, np.nan) for k in keys])
    family_floor = runs.groupby('family')['floor_mm'].median()
    rows = []
    for row in runs.itertuples(index=False):
        work = dev_real / 'work' / f'{row.dataset}__{row.run}'
        arms, _ = run_real.collect(work, [])
        points = run_real.brain_points(work)
        run_floor = row.floor_mm if np.isfinite(row.floor_mm) else family_floor.get(row.family)
        separation = max(MIN_SEPARATION_MM, run_floor if np.isfinite(run_floor) else 0.0)
        pairs = [('sign', 'FSL conventional default', 'FSL conventional opposite')]
        if str(row.elig_eligible) == 'True':
            pairs.append(('R1 vs R2', r1_arm(row.contrast), R2_ARM))

        for comparison, a, b in pairs:
            if a not in arms or b not in arms:
                continue

            rows.append(
                {
                    'comparison': comparison,
                    'dataset': row.dataset,
                    'subject': row.subject,
                    'run': row.run,
                    'family': row.family,
                    'arm_a': a,
                    'arm_b': b,
                    'distance_mm': scoring.distance(arms[a], arms[b], points),
                    'separation_mm': separation,
                    'work': str(work),
                },
            )
    out = pd.DataFrame(rows)
    out['candidate'] = out['distance_mm'] > out['separation_mm']
    return out


def draw(pool, n, rng):
    """Up to ``n`` rows, round-robin over families, at most one per subject, random by ``rng``."""
    pool = pool.iloc[rng.permutation(len(pool))]
    by_family = {f: list(g.index) for f, g in pool.groupby('family', sort=True)}
    chosen, subjects = [], set()
    while len(chosen) < n and any(by_family.values()):
        for family in sorted(by_family):
            queue = by_family[family]
            while queue:
                i = queue.pop(0)
                subject = (pool.at[i, 'dataset'], pool.at[i, 'subject'])
                if subject not in subjects:
                    subjects.add(subject)
                    chosen.append(i)
                    break

            if len(chosen) == n:
                break
    return pool.loc[chosen]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--freeze', type=Path, required=True)
    p.add_argument('--dev-real', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--seed', type=int, default=SEED)
    args = p.parse_args(argv)
    real = pd.read_table(args.freeze / 'real.tsv')
    fsl = real[real['arm'].str.startswith('FSL')]
    every = fsl.groupby(['dataset', 'run'])['verdict'].apply(lambda s: (s == 'fail').all())
    shared = {k for k, v in every.items() if v}
    real['shared_failure'] = [k in shared for k in zip(real['dataset'], real['run'], strict=True)]
    manifest = pd.read_table(HERE / 'manifest.tsv')[['dataset', 'run', 'subject']]
    real = real.merge(manifest, on=['dataset', 'run'], how='left')

    cand = candidates(real, args.dev_real)
    rng = np.random.default_rng(args.seed)
    chosen = pd.concat(
        [
            draw(cand[cand['candidate'] & (cand['comparison'] == c)], PER_COMPARISON, rng)
            for c in ('R1 vs R2', 'sign')
        ],
    )
    cand['chosen'] = cand.index.isin(chosen.index)
    print(cand.groupby(['comparison', 'family'])[['candidate', 'chosen']].sum().to_string())

    pairs = []
    for row in chosen.itertuples(index=False):
        arms, _ = run_real.collect(Path(row.work), [])
        label = f'{row.dataset}/{row.run}'
        pairs.append(
            {
                'run': label,
                'work': Path(row.work),
                'arm_a': row.arm_a,
                'arm_b': row.arm_b,
                'pull_a': arms[row.arm_a],
                'pull_b': arms[row.arm_b],
            },
        )

    runs = real.drop_duplicates(['dataset', 'run'])
    runs = runs[~runs['shared_failure']].reset_index(drop=True)
    runs['work'] = [
        str(args.dev_real / 'work' / f'{d}__{r}')
        for d, r in zip(runs['dataset'], runs['run'], strict=True)
    ]
    known_runs = draw(runs, KNOWN_PER_SIZE * len(review.KNOWN_MM), rng)
    sizes = np.repeat(review.KNOWN_MM, KNOWN_PER_SIZE)
    known = []
    for row, size in zip(known_runs.itertuples(index=False), sizes, strict=True):
        arms, _ = run_real.collect(Path(row.work), [])
        arm = r1_arm(row.contrast)
        known.append(
            {
                'run': f'{row.dataset}/{row.run}',
                'work': Path(row.work),
                'arm': arm,
                'pull': arms[arm],
                'size_mm': float(size),
            },
        )

    batch, key = review.build_batch(pairs, known, args.out, args.seed)
    cand.drop(columns='work').to_csv(args.out / 'selection.tsv', sep='\t', index=False)
    words = {'FSL', 'ANTs', 'BBR', 'conventional', 'opposite', 'default', 'P-mean', 'moved'}
    words |= {w for k in key.values() for w in (k['run'], *k['run'].split('/'))}
    found = review.leaks(args.out, words)
    kinds = pd.Series([k['kind'] for k in key.values()]).value_counts().to_dict()
    print(json.dumps({'figures': len(batch), 'kinds': kinds}))
    print('LEAKS: ' + ', '.join(found) if found else 'NO LEAKS')
    return 1 if found else 0


if __name__ == '__main__':
    sys.exit(main())
