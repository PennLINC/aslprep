#!/usr/bin/env python3
"""Score ASLPrep's segmentation of each ``raw`` fixture against the truth (plan: "Design").

A ``raw`` fixture (``<name>_raw``, from ``make_recipes.py``) hands ASLPrep only the T1w; its
``derivatives`` counterpart (``<name>``) carries the truth brain mask and tissue maps of the same
T1w. Labels are the argmax of the three tissue probabilities inside the brain mask; ASLPrep's are
sampled at the truth's voxel centres (nearest neighbour, by world coordinates). A fixture is used
only if it passes every threshold below, which were fixed before any segmentation was seen.

Usage (in WSL)::

    python scripts/coreg_mc/check_raw.py --fixtures <aslscan fixtures> --runs <fixture runs> \\
        --names rc_pcasl_raw ... --out <tsv>
"""

import argparse
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
from nibabel.affines import apply_affine

#: Minimum Dice per structure (2026-10-29, Task 3); CSF is reported, not required.
THRESHOLDS = {'brain': 0.95, 'GM': 0.80, 'WM': 0.80}
TISSUES = ('GM', 'WM', 'CSF')


def labels(anat, glob_prefix=''):
    """(brain mask, tissue labels 0=outside 1=GM 2=WM 3=CSF, affine) from a derivatives folder."""
    mask_img = nb.load(next(anat.glob(f'{glob_prefix}*_desc-brain_mask.nii.gz')))
    mask = np.asanyarray(mask_img.dataobj) > 0
    probs = [
        np.asanyarray(nb.load(next(anat.glob(f'{glob_prefix}*_label-{t}_probseg.nii.gz'))).dataobj)
        for t in TISSUES
    ]
    tissue = np.argmax(np.stack(probs), axis=0) + 1
    return mask, np.where(mask, tissue, 0), mask_img.affine


def sample(volume, affine, points):
    """Nearest-neighbour values of ``volume`` at world ``points`` (0 outside its grid)."""
    ijk = np.rint(apply_affine(np.linalg.inv(affine), points)).astype(int)
    inside = np.all((ijk >= 0) & (ijk < volume.shape), axis=1)
    out = np.zeros(len(points), volume.dtype)
    out[inside] = volume[tuple(ijk[inside].T)]
    return out


def dice(a, b):
    return float(2 * (a & b).sum() / max(a.sum() + b.sum(), 1))


def score(truth_anat, run_anat):
    t_mask, t_labels, t_affine = labels(truth_anat)
    r_mask, r_labels, r_affine = labels(run_anat)
    ijk = np.argwhere(np.ones(t_mask.shape, bool))
    points = apply_affine(t_affine, ijk)
    got_mask = sample(r_mask, r_affine, points).reshape(t_mask.shape)
    got_labels = sample(r_labels, r_affine, points).reshape(t_mask.shape)
    out = {'brain': dice(t_mask, got_mask)}
    for i, t in enumerate(TISSUES, start=1):
        out[t] = dice(t_labels == i, got_labels == i)
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--fixtures', type=Path, required=True)
    p.add_argument('--runs', type=Path, required=True)
    p.add_argument('--names', nargs='+', required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args(argv)
    rows = []
    for name in args.names:
        base = name.removesuffix('_raw')
        truth = args.fixtures / base / 'derivatives' / 'anat' / 'sub-01' / 'anat'
        row = {'fixture': name} | score(truth, args.runs / name / 'sub-01' / 'anat')
        row['passes'] = all(row[k] >= v for k, v in THRESHOLDS.items())
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    print(df.round(3).to_string(index=False))
    return 0 if df['passes'].all() else 1


if __name__ == '__main__':
    sys.exit(main())
