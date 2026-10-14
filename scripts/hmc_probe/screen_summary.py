#!/usr/bin/env python3
"""Summarize the motion screen (``batch_harness.py`` output) by family and motion stratum.

Reads ``<runs>/<label>/<dataset>/<run>/motion.json`` and prints, per acquisition family, the
number of runs and the distribution of mean framewise displacement, with each run assigned to
a stratum (low / medium / high) by the tertiles of mean framewise displacement within its
family. Writes the table to ``--tsv``.

These are the baseline's own estimates, so the strata describe what production sees, not
the true motion (the plan's caveat on the motion screen).

Usage::

    python scripts/hmc_probe/screen_summary.py --runs <dir> --label c0 [--tsv out.tsv]
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--runs', type=Path, required=True)
    parser.add_argument('--label', default='c0')
    parser.add_argument('--tsv', type=Path)
    args = parser.parse_args(argv)
    rows = [
        json.loads(f.read_text()) for f in sorted((args.runs / args.label).glob('*/*/motion.json'))
    ]
    df = pd.DataFrame(rows)
    df['stratum'] = ''
    for _, sub in df.groupby('family'):
        if len(sub) >= 3:
            q1, q2 = np.quantile(sub['fd_mean_mm'], [1 / 3, 2 / 3])
            df.loc[sub.index, 'stratum'] = np.where(
                sub['fd_mean_mm'] <= q1, 'low', np.where(sub['fd_mean_mm'] <= q2, 'medium', 'high')
            )
        else:
            df.loc[sub.index, 'stratum'] = 'n<3'
    print(f'{len(df)} runs')
    summary = df.groupby('family').agg(
        runs=('run', 'size'),
        datasets=('dataset', 'nunique'),
        fd_median=('fd_mean_mm', 'median'),
        fd_p90=('fd_mean_mm', lambda x: np.quantile(x, 0.9)),
        fd_max=('fd_max_mm', 'max'),
        pair_median=('pair_median_mm', 'median'),
        offset_median=('offset_median_mm', 'median'),
    )
    with pd.option_context('display.width', 200, 'display.max_columns', 20):
        print(summary.round(3))
        print('\nStillest runs per family (candidates for semi-synthetic motion):')
        print(
            df.sort_values('fd_mean_mm')
            .groupby('family')
            .head(3)[['family', 'dataset', 'run', 'n_volumes', 'fd_mean_mm', 'fd_max_mm']]
            .round(3)
            .to_string(index=False)
        )
        print('\nMost moving runs:')
        print(
            df.sort_values('fd_mean_mm', ascending=False)
            .head(12)[['family', 'dataset', 'run', 'fd_mean_mm', 'fd_max_mm', 'pair_max_mm']]
            .round(3)
            .to_string(index=False)
        )
    if args.tsv:
        df.to_csv(args.tsv, sep='\t', index=False)


if __name__ == '__main__':
    main()
