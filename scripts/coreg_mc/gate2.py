#!/usr/bin/env python3
"""Gate 2, development part, on Task 4's screen (plan: Task 4 rules, fixed 2026-11-02).

For each J1 variant and its matched ablations (S1 with the same masking; S2 and Q1 with the same
perfusion image, target and masking), on each family and form (``derivatives``, ``raw``):

- the variant is **better** on a stress family when its median error (voxels) is lower than
  each ablation's by at least ``MARGIN_VOX`` on every form the family has;
- it is **not worse** on any other family when its median is at most the best ablation's plus
  ``MARGIN_VOX`` on every form.

A variant passes with at least one stress family won and no family failed. At most ``KEEP``
passing variants go on: most families won, then lowest median error over all fixtures, then the
simpler variant (P-winsor, GM target, weight 0.4, masked: the frozen J1).

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/gate2.py --screen screen.tsv --out gate2.tsv
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

#: Frozen 2026-11-02 (plan: "Frozen from the development review", margins).
MARGIN_VOX = 0.1
KEEP = 2
#: fixture -> (family, form). Stress families first; ``other`` is "elsewhere".
FIXTURES = {
    'st_spiral': ('st', 'derivatives'),
    'st_vessels': ('st', 'derivatives'),
    'st_noperf': ('st', 'derivatives'),
    'st_spiral_raw': ('st', 'raw'),
    'st_vessels_raw': ('st', 'raw'),
    'st_noperf_raw': ('st', 'raw'),
    'rc_grase': ('grase', 'derivatives'),
    'rc_grase_off2': ('grase', 'derivatives'),
    'rc_grase_raw': ('grase', 'raw'),
    'mo_pcaslshort2_drift_off': ('short', 'derivatives'),
    'mo_pcaslshort4_drift_off': ('short', 'derivatives'),
    'rc_pcasl': ('other', 'derivatives'),
    'rc_pcasl_flip': ('other', 'derivatives'),
    'rc_pcasl_grid': ('other', 'derivatives'),
    'rc_pcasl_off2': ('other', 'derivatives'),
    'rc_pasl_multipld': ('other', 'derivatives'),
    'rc_pcasl_raw': ('other', 'raw'),
}
STRESS = ('st', 'grase', 'short')


def simplicity(row):
    """How many factors differ from the frozen J1 (fewer is simpler)."""
    return (
        (row['image'] != 'P-winsor')
        + (row['target'] != 'gm')
        + (row['weight'] != 0.4)
        + (not row['masked'])
    )


def judge(screen):
    df = screen.copy()
    missing = sorted(set(df['fixture']) ^ set(FIXTURES))
    if missing:
        raise ValueError(f'fixtures not in the gate 2 table or not screened: {missing}')

    if df['rms_voxels'].isna().any():
        raise ValueError('an arm has no transform; gate 2 needs every arm on every fixture')

    df[['family', 'form']] = pd.DataFrame(df['fixture'].map(FIXTURES).tolist(), index=df.index)
    med = df.groupby(['name', 'family', 'form'])['rms_voxels'].median()
    overall = df.groupby('name')['rms_voxels'].median()
    rows = []
    for v in df[df['arm'] == 'J1'].drop_duplicates('name').itertuples(index=False):
        tag = 'mask' if v.masked else 'nomask'
        ablations = [
            f'S1_{tag}',
            f'S2_{v.image}_{v.target}_{tag}',
            f'Q1_{v.image}_{v.target}_{tag}',
        ]
        record = {
            'name': v.name,
            'image': v.image,
            'target': v.target,
            'weight': v.weight,
            'masked': v.masked,
            'median_vox': overall[v.name],
        }
        won, failed = [], []
        for family in ('st', 'grase', 'short', 'other'):
            gaps = {}
            for form in med.loc[v.name, family].index:
                ablation = np.array([med.loc[a, family, form] for a in ablations])
                gaps[form] = med.loc[v.name, family, form] - ablation

            # A gap is the variant's median less an ablation's: the largest gap is to the best
            # ablation (lowest error), the smallest to the worst.
            # Better: below every ablation by the margin, on every form.
            better = all(g.max() <= -MARGIN_VOX for g in gaps.values())
            # Not worse: within the margin of the best ablation, on every form.
            ok = all(g.max() <= MARGIN_VOX for g in gaps.values())
            record[f'{family}_vs_best'] = max(float(g.max()) for g in gaps.values())
            record[f'{family}_vs_worst'] = max(float(g.min()) for g in gaps.values())
            if family in STRESS and better:
                won.append(family)
            elif not ok:
                failed.append(family)

        record |= {'won': ','.join(won), 'n_won': len(won), 'failed': ','.join(failed)}
        record['passes'] = bool(won) and not failed
        record['simplicity'] = simplicity(record)
        rows.append(record)

    out = pd.DataFrame(rows).sort_values(
        ['passes', 'n_won', 'median_vox', 'simplicity'],
        ascending=[False, False, True, True],
    )
    out['kept'] = False
    kept = out[out['passes']].index[:KEEP]
    out.loc[kept, 'kept'] = True
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--screen', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args(argv)
    out = judge(pd.read_table(args.screen))
    out.to_csv(args.out, sep='\t', index=False)
    with pd.option_context('display.width', 250, 'display.max_rows', 50):
        print(out.round(3).to_string(index=False))

    kept = list(out.loc[out['kept'], 'name'])
    print(f'{int(out["passes"].sum())} of {len(out)} variants pass; kept: {kept}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
