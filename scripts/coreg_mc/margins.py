#!/usr/bin/env python3
"""The yardsticks behind Task 3's margins (plan: "Decision gates", frozen in Task 3).

Development data only. Writes one table (``yardstick``, ``group``, ``n``, ``value``) with:

- ``equivalent_vox``: on the development fixtures, each D3 level's RMS error minus R1's, in
  acquisition voxels (median, 90th percentile, maximum). Arms of practically equal quality; the
  accuracy margin of gates 2 and 3 must sit above this.
- ``r1_error_vox``: R1's own error on the fixtures, in voxels (median, maximum).
- ``r1_inconsistency_mm``: R1's per-subject inconsistency (median over a subject's pairs, the
  conventional frame, pairs within the floor limit), by family and pair type.
- ``chance_mm``: the 95th percentile, over a bootstrap of subjects, of the family median of a
  zero-centred noise difference: R1's per-subject inconsistency measured in the conventional
  frame minus the perfusion frame. How far a family median moves by chance; the inconsistency
  margin of gate 4 should sit above it.
- ``worse_mm``: arms the review judged worse (P-mean; the opposite sign) minus R1, the family
  median of per-subject differences. What a margin can and cannot catch: consistency does not
  see an arm that is consistently wrong.

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/margins.py --freeze <freeze dir> --capture <comparators.tsv> \\
        --dev-real <dir> --out margins.tsv
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'hmc_probe'))

import dev_review  # noqa: E402
import freeze_comparators  # noqa: E402

#: Frozen 2026-11-02 (plan: "Frozen from the development review").
FLOOR_LIMIT_MM = 3.0
R1_DEFAULT, R1_OPPOSITE = 'FSL conventional default', 'FSL conventional opposite'
LEVELS = ('P-mean', 'P-robust', 'P-winsor', 'P-masked')
WORSE = ('FSL P-mean default', R1_OPPOSITE)
BOOTSTRAP = 2000
SEED = 20261102


def fixture_rows(freeze, capture):
    fix = pd.read_table(freeze / 'fixtures.tsv')
    voxel = pd.read_table(capture).drop_duplicates('fixture').set_index('fixture')['voxel_mm']
    fix['rms_vox'] = fix['rms_mm'] / fix['fixture'].map(voxel)
    w = fix.pivot(index='fixture', columns='arm', values='rms_vox')
    contrast = fix.drop_duplicates('fixture').set_index('fixture')['contrast']
    opposite = contrast.map(freeze_comparators.r1_sign) == 'opposite'
    r1 = w[R1_DEFAULT].where(~opposite, w[R1_OPPOSITE])
    rows = [
        ('r1_error_vox', 'median', len(r1), r1.median()),
        ('r1_error_vox', 'max', len(r1), r1.max()),
    ]
    diffs = pd.concat([w[f'FSL {level} default'] - r1 for level in LEVELS])
    rows += [
        ('equivalent_vox', 'median', len(diffs), diffs.median()),
        ('equivalent_vox', 'p90', len(diffs), diffs.quantile(0.9)),
        ('equivalent_vox', 'max', len(diffs), diffs.max()),
    ]
    return rows


def per_subject(dev_real):
    """Per-subject median inconsistency: (family, pair type, subject, frame) x arm."""
    cons = dev_review.run_real_tables(dev_real, 'consistency')
    cons = cons[cons['floor_mm'] <= FLOOR_LIMIT_MM]
    manifest = pd.read_table(HERE / 'manifest.tsv').set_index(['dataset', 'run'])['family']
    keys = zip(cons['dataset'], cons['run_a'], strict=True)
    cons = cons.assign(family=[manifest.get(k) for k in keys])
    index = ['family', 'pair_type', 'dataset', 'subject', 'frame', 'arm']
    return cons.groupby(index)['inconsistency_mm'].median().unstack('arm')


def chance(values, rng):
    """95th percentile of |median| of a bootstrap resample of the zero-centred values."""
    x = values - np.median(values)
    medians = [abs(np.median(rng.choice(x, len(x)))) for _ in range(BOOTSTRAP)]
    return float(np.quantile(medians, 0.95))


def real_rows(dev_real):
    subj = per_subject(dev_real)
    conv = subj.xs('conventional', level='frame')
    frames = subj[R1_DEFAULT].unstack('frame').dropna()
    noise = frames['conventional'] - frames['perfusion']
    rng = np.random.default_rng(SEED)
    rows = []
    for (family, pair_type), g in conv[R1_DEFAULT].dropna().groupby(['family', 'pair_type']):
        group = f'{family} / {pair_type}'
        rows.append(('r1_inconsistency_mm', group, len(g), g.median()))
        n = noise.xs((family, pair_type), level=('family', 'pair_type'), drop_level=False)
        if len(n) >= 3:
            rows.append(('chance_mm', group, len(n), chance(n.to_numpy(), rng)))

        for arm in WORSE:
            d = (conv[arm] - conv[R1_DEFAULT]).dropna()
            d = d.xs((family, pair_type), level=('family', 'pair_type'))
            rows.append((f'worse_mm: {arm}', group, len(d), d.median()))

    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--freeze', type=Path, required=True)
    p.add_argument('--capture', type=Path, required=True)
    p.add_argument('--dev-real', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args(argv)
    rows = fixture_rows(args.freeze, args.capture) + real_rows(args.dev_real)
    out = pd.DataFrame(rows, columns=['yardstick', 'group', 'n', 'value'])
    out.to_csv(args.out, sep='\t', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(out.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
