#!/usr/bin/env python3
"""Tabulate fixture scores across candidates (``<scores>/<label>/<fixture>/score.json``).

Noise realizations (``<name>__s<seed>``) are pooled: the table gives their mean and, with
more than one, the spread (max - min). Columns follow the plan's metrics: relative pose
(1), reference consistency (2), pair (3), delta-M against the mean-frame and reference-frame
oracles in the interior and the outer 0-3 mm (4), and framewise displacement error (7).

Usage::

    python scripts/hmc_probe/compare.py --scores <dir> [--labels c0 c1 ...] [--tsv out.tsv]
"""

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

METRICS = {
    'rel_max': ('relative_pose', 'to_mean', 'max_mm'),
    'ref_med': ('reference_consistency', 'median_mm'),
    'pair_max': ('pair', 'max_mm'),
    'fd_mae': ('framewise', 'mean_abs_error_mm'),
    'dm_int': ('deltam', 'oracle_mean_frame', 'interior', 'rmse_rel_gm'),
    'dm_edge': ('deltam', 'oracle_mean_frame', 'shell_0_3mm', 'rmse_rel_gm'),
    'dmref_int': ('deltam', 'oracle_reference_frame', 'interior', 'rmse_rel_gm'),
    'dmref_edge': ('deltam', 'oracle_reference_frame', 'shell_0_3mm', 'rmse_rel_gm'),
}


def _get(d, path):
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return np.nan
        d = d[k]
    return np.nan if d is None else float(d)


def collect(scores, labels):
    rows = defaultdict(lambda: defaultdict(list))
    for label in labels:
        for f in sorted((scores / label).glob('*/score.json')):
            fixture = re.sub(r'__s\d+$', '', f.parent.name)
            s = json.loads(f.read_text())
            for m, path in METRICS.items():
                rows[(fixture, label)][m].append(_get(s, path))
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--scores', type=Path, required=True)
    parser.add_argument('--labels', nargs='*')
    parser.add_argument('--tsv', type=Path)
    args = parser.parse_args(argv)
    labels = args.labels or sorted(p.name for p in args.scores.iterdir() if p.is_dir())
    rows = collect(args.scores, labels)
    header = ['fixture', 'label', 'n', *METRICS]
    lines = ['\t'.join(header)]
    print(f'{"fixture":28s} {"label":8s} n ' + ' '.join(f'{m:>10s}' for m in METRICS))
    for (fixture, label), vals in sorted(rows.items()):
        n = len(next(iter(vals.values())))
        cells = []
        for m in METRICS:
            v = np.asarray(vals[m], float)
            mean = np.nanmean(v) if np.isfinite(v).any() else np.nan
            spread = np.nanmax(v) - np.nanmin(v) if np.isfinite(v).sum() > 1 else np.nan
            cells.append((mean, spread))
        print(
            f'{fixture:28s} {label:8s} {n} '
            + ' '.join(
                f'{m:6.3f}' + (f'±{s:<3.2f}' if np.isfinite(s) else '    ') for m, s in cells
            )
        )
        lines.append('\t'.join([fixture, label, str(n), *(f'{m:.4f}' for m, _ in cells)]))
    if args.tsv:
        args.tsv.write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
