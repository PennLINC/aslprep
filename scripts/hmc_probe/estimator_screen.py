#!/usr/bin/env python3
"""Task 5 estimator screen: the best variant of each estimator family (plan: "Bounding the
combinations", stage 1).

Every arm in ``<scores>/<label>/`` is assigned to a family by its label prefix (``c1``, ``c2``,
``c4``, ``c8``; ``c0`` and the #513 comparators are reported but are not families that go
forward). Each arm is summarized over the development fixtures it shares with C0 by the median
of the two screening metrics fixed in the plan, relative pose (``rel_max``, metric 1) and pair
error (``pair_max``, metric 3), and ranked within its family by the sum of its ranks on the
two. The edge delta-M error is reported alongside but does not rank. Arms missing fixtures are
flagged.

Usage::

    python estimator_screen.py --scores <fixture-scores> [--out screen.tsv]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from report import load_scores  # noqa: E402

FAMILIES = ('c1', 'c2', 'c4', 'c8')
REFERENCE = ('c0', 'c3', 'fixvol', 'perfref')
METRICS = ('rel_max', 'pair_max')
#: Reported but never short-listed (plan: C8's diagnostic-only runs; the leave-one-out
#: template is a null test of C1, not a variant).
DIAGNOSTIC = ('c8-noreg', 'c8-prod', 'c8-prod-huber', 'c8-spm', 'c1-loo')


def family(label):
    for f in FAMILIES:
        if label == f or label.startswith(f'{f}-'):
            return f
    return label if label in REFERENCE else None


def screen(df):
    # Downstream-policy runs (``-d0``/``-d1``) belong to the policy screen, not to this one.
    df = df[~df.label.str.endswith(('-d0', '-d1'))]
    df = df[df.label.map(family).notna()].copy()
    df['family'] = df.label.map(family)
    fixtures = sorted(set(df.query("label == 'c0'").fixture))
    df = df[df.fixture.isin(fixtures)]
    summary = (
        df.groupby(['family', 'label'])
        .agg(
            n=('fixture', 'nunique'),
            rel_max=('rel_max', 'median'),
            pair_max=('pair_max', lambda s: float(np.nanmedian(s)) if s.notna().any() else np.nan),
            dm_edge=('dm_edge', 'median'),
        )
        .reset_index()
    )
    summary['complete'] = summary.n == len(fixtures)
    summary['diagnostic'] = summary.label.isin(DIAGNOSTIC)
    summary['rank'] = np.nan
    eligible = summary[summary.family.isin(FAMILIES) & ~summary.diagnostic]
    for _, g in eligible.groupby('family'):
        ranks = sum(g[m].rank(method='min') for m in METRICS)
        # Ties on the rank sum are broken by the pair error, which drives subtraction.
        order = pd.DataFrame({'r': ranks, 'p': g.pair_max}).sort_values(['r', 'p'])
        summary.loc[order.index, 'rank'] = range(1, len(order) + 1)
    return summary, len(fixtures)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--scores', type=Path, required=True)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args(argv)
    summary, n = screen(load_scores(args.scores.parent))
    if args.out:
        summary.to_csv(args.out, sep='\t', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(f'{n} development fixtures')
        print(summary.sort_values(['family', 'rank', 'label']).round(3).to_string(index=False))
    best = summary[summary.complete & summary['rank'].eq(1)].groupby('family').label.first()
    print('short list:', ['c0', *best.tolist()])


if __name__ == '__main__':
    main()
