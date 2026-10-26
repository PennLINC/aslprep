#!/usr/bin/env python3
"""Do truth-free control-label difference metrics rank HMC candidates as the truth does?

For every scored harness run in ``<scores>/<label>/<fixture>/`` (its ``deltam.nii.gz``: one
control-minus-label volume per pair, on the run's reference grid), with a brain mask estimated
from the run's own reference (nilearn's EPI mask, as on real data):

- ``ssd_brain``: mean over brain voxels of the squared mean delta-M;
- ``ssd_pairs``: mean over brain voxels and pairs of each pair's squared delta-M;
- ``var_pairs``: mean over brain voxels of the pair-to-pair variance of delta-M (blind to an
  offset shared by every pair);
- ``ssd_edge``: as ``ssd_pairs``, inside the brain within 3 voxels of the mask edge;
- ``ssd_outer``: as ``ssd_pairs``, in a 3-voxel shell outside the mask, where true perfusion is
  near zero.

Each is compared, within each fixture across candidates, with the truth-based edge delta-M
error (``dm_edge`` in ``score.json``) by Spearman correlation; the median over fixtures is the
summary. Fixtures without pairs are skipped.

Usage::

    python ssd_check.py --scores <fixture-scores> --out <tsv>
"""

import argparse
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
from scipy.ndimage import binary_dilation, binary_erosion
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from perfusion import epi_mask  # noqa: E402
from report import load_scores  # noqa: E402

METRICS = ('ssd_brain', 'ssd_brain_shared', 'ssd_pairs', 'var_pairs', 'ssd_edge', 'ssd_outer')


def metrics(run_dir):
    dm = np.asanyarray(nb.load(run_dir / 'deltam.nii.gz').dataobj, np.float32)
    dm = dm[..., None] if dm.ndim == 3 else dm
    brain = epi_mask(run_dir / 'aslref.nii.gz')
    edge = brain & ~binary_erosion(brain, iterations=3)
    outer = binary_dilation(brain, iterations=3) & ~brain
    sq = dm**2
    return {
        'ssd_brain': float((dm.mean(-1)[brain] ** 2).mean()),
        'ssd_pairs': float(sq[brain].mean()),
        'var_pairs': float(dm[brain].var(-1).mean()) if dm.shape[-1] > 1 else np.nan,
        'ssd_edge': float(sq[edge].mean()),
        'ssd_outer': float(sq[outer].mean()),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--scores', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument(
        '--exclude', nargs='*', default=['-d0', '-d1'], help='Label suffixes to skip.'
    )
    args = parser.parse_args(argv)
    scores = load_scores(args.scores.parent)
    scores = scores[~scores.label.str.endswith(tuple(args.exclude))]
    rows = []
    for r in scores.itertuples():
        if not (r.path / 'deltam.nii.gz').exists() or not np.isfinite(r.dm_edge):
            continue
        rows.append(
            {
                'label': r.label,
                'fixture': r.fixture,
                'dm_edge': r.dm_edge,
                'pair_max': r.pair_max,
                **metrics(r.path),
            }
        )
    df = pd.DataFrame(rows)
    # The same metric inside one mask per fixture: the union of every arm's estimated mask.
    shared = {}
    for fx, g in scores[scores.fixture.isin(set(df.fixture))].groupby('fixture'):
        masks = [epi_mask(p / 'aslref.nii.gz') for p in g.path if (p / 'deltam.nii.gz').exists()]
        shapes = {m.shape for m in masks}
        shared[fx] = np.logical_or.reduce(masks) if len(shapes) == 1 else None
    vals = []
    for r in df.itertuples():
        mask = shared.get(r.fixture)
        path = scores[(scores.label == r.label) & (scores.fixture == r.fixture)].path.iloc[0]
        dm = np.asanyarray(nb.load(path / 'deltam.nii.gz').dataobj, np.float32)
        dm = dm[..., None] if dm.ndim == 3 else dm
        ok = mask is not None and mask.shape == dm.shape[:3]
        vals.append(float((dm.mean(-1)[mask] ** 2).mean()) if ok else np.nan)
    df['ssd_brain_shared'] = vals
    df.to_csv(args.out, sep='\t', index=False)
    corr = []
    for fx, g in df.groupby('fixture'):
        if len(g) < 5:
            continue
        corr.append(
            {
                'fixture': fx,
                'n_arms': len(g),
                **{m: spearmanr(g[m], g.dm_edge, nan_policy='omit')[0] for m in METRICS},
            }
        )
    c = pd.DataFrame(corr)
    with pd.option_context('display.width', 200, 'display.max_rows', 100):
        print(c.round(2).to_string(index=False))
        print('median rho:', c[list(METRICS)].median().round(2).to_dict())
        print(
            'fraction of fixtures with rho > 0.5:',
            (c[list(METRICS)] > 0.5).mean().round(2).to_dict(),
        )


if __name__ == '__main__':
    main()
