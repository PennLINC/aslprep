#!/usr/bin/env python3
"""D3 perfusion-image builder and eligibility statistic (plan: "D3: the perfusion image").

The same functions are meant to run in the harness and, in Task 7, in the workflow. They use
no truth and no kinetic fit; inputs are the motion-corrected volumes on the raw grid, the
volume types, the post-labeling delays and the conventional brain mask.

Inputs by data type:

- control-label pairs: the pairwise differences (pairs matched in order of appearance); for
  multi-delay data only the pairs at delays at or above the run's median delay;
- repeated derived volumes (a delta-M or CBF series): those volumes;
- a single derived volume: that volume (the "D3-single" arm, never eligible by the rules below).

Levels: ``mean``, ``robust`` (voxelwise median; equal to the mean for one input), ``winsor``
(``robust`` clipped above the 99th percentile inside the mask) and ``masked`` (the ``winsor``
image; the mask is applied by the registration as an input weight, not here).

Eligibility: for pairs, at least ``min_pairs`` pairs after delay selection and an edge-content
statistic at or above ``limit``; for repeated derived volumes, the same statistic over halves of
the series. The statistic is the correlation, inside the mask, between the gradient magnitudes
of the two split-half images (alternate inputs, balanced within each delay); the split-half
correlation of the values is reported alongside. ``min_pairs`` and ``limit`` are placeholders
until Task 5b calibrates them.

Usage (on a harness output directory)::

    python perfusion.py --harness <dir> --metadata <asl.json> [--mask <mask.nii.gz>]
"""

import argparse
import json
from pathlib import Path

import nibabel as nb
import numpy as np

LEVELS = ('mean', 'robust', 'winsor', 'masked')
DERIVED = ('deltam', 'cbf')
#: Placeholders until Task 5b calibrates them.
MIN_PAIRS = 4
EDGE_LIMIT = 0.5


def select_inputs(data, context, plds):
    """The builder's inputs: (stack on the last axis, data type, delay of each input)."""
    controls = [i for i, t in enumerate(context) if t == 'control']
    labels = [i for i, t in enumerate(context) if t == 'label']
    if controls and labels:
        if len(controls) != len(labels):
            raise ValueError(f'{len(controls)} controls and {len(labels)} labels')
        pair_plds = np.array([plds[c] for c in controls], float)
        keep = pair_plds >= np.median(pair_plds)
        diffs = [
            data[..., c] - data[..., lab]
            for c, lab, k in zip(controls, labels, keep, strict=True)
            if k
        ]
        return np.stack(diffs, -1), 'pairs', pair_plds[keep]
    for kind in DERIVED:
        sel = [i for i, t in enumerate(context) if t == kind]
        if sel:
            dtype = 'derived_series' if len(sel) > 1 else 'derived_single'
            return data[..., sel], dtype, np.array([plds[i] for i in sel], float)
    raise ValueError('no control-label pairs and no derived volumes')


def combine(stack, level, mask=None):
    """One perfusion image from the inputs at ``level``."""
    if level not in LEVELS:
        raise ValueError(f'Unknown level: {level}')
    if level == 'mean':
        return stack.mean(-1)
    img = np.median(stack, axis=-1)
    if level in ('winsor', 'masked'):
        inside = img[mask] if mask is not None else img.ravel()
        img = np.minimum(img, np.percentile(inside, 99))
    return img


def split_halves(plds):
    """Alternate inputs within each delay: two index lists, balanced across delays."""
    halves = ([], [])
    for d in np.unique(plds):
        for j, i in enumerate(np.flatnonzero(plds == d)):
            halves[j % 2].append(int(i))
    return halves


def _grad(img, zooms):
    return np.linalg.norm(np.gradient(img, *zooms), axis=0)


def eligibility(
    stack, dtype, plds, mask, zooms, level='robust', min_pairs=MIN_PAIRS, limit=EDGE_LIMIT
):
    """The eligibility record of one run (see the module docstring)."""
    n = stack.shape[-1]
    out = {'data_type': dtype, 'n_inputs': int(n), 'min_pairs': min_pairs, 'limit': limit}
    if dtype == 'derived_single':
        return out | {'eligible': False, 'arm': 'D3-single', 'reason': 'no replicate'}
    # Both halves must hold every retained delay, or they differ in delay composition (as in
    # real_score.py's split-half rule).
    counts = {float(d): int((np.asarray(plds) == d).sum()) for d in np.unique(plds)}
    if min(counts.values()) < 2:
        return out | {'eligible': False, 'reason': 'a retained delay has fewer than two inputs'}
    a, b = split_halves(plds)
    ia, ib = combine(stack[..., a], level, mask), combine(stack[..., b], level, mask)
    edge = float(np.corrcoef(_grad(ia, zooms)[mask], _grad(ib, zooms)[mask])[0, 1])
    value = float(np.corrcoef(ia[mask], ib[mask])[0, 1])
    out |= {'edge_r': edge, 'value_r': value}
    if not np.isfinite(edge):
        return out | {'eligible': False, 'reason': 'non-finite statistic'}
    if dtype == 'pairs' and n < min_pairs:
        return out | {'eligible': False, 'reason': f'{n} pairs < {min_pairs}'}
    if edge < limit:
        return out | {'eligible': False, 'reason': f'edge_r {edge:.3f} < {limit}'}
    return out | {'eligible': True}


def epi_mask(ref_file):
    """A conventional brain mask of a reference (nilearn's EPI mask), for runs without one."""
    from nilearn.masking import compute_epi_mask

    ref = nb.load(ref_file)  # as float: nilearn fails on integer images
    ref = nb.Nifti1Image(np.asanyarray(ref.dataobj, np.float32), ref.affine)
    return np.asanyarray(compute_epi_mask(ref).dataobj).astype(bool)


def from_harness(harness_dir, metadata, mask=None):
    """Build every level and the eligibility record from a harness output directory."""
    import pandas as pd

    harness_dir = Path(harness_dir)
    corrected = nb.load(harness_dir / 'corrected.nii.gz')
    data = np.asanyarray(corrected.dataobj, np.float32)
    data = data[..., None] if data.ndim == 3 else data
    context = pd.read_table(harness_dir / 'hmc_aslcontext.tsv')['volume_type'].tolist()
    indices = json.loads((harness_dir / 'hmc_indices.json').read_text())
    pld = np.atleast_1d(np.asarray(metadata.get('PostLabelingDelay', 0.0), float))
    n_raw = max(indices) + 1
    pld = np.repeat(pld, n_raw) if pld.size == 1 else pld
    plds = [float(pld[i]) for i in indices]
    if mask is None:
        mask = epi_mask(harness_dir / 'aslref.nii.gz')
    stack, dtype, input_plds = select_inputs(data, context, plds)
    zooms = corrected.header.get_zooms()[:3]
    images = {lvl: combine(stack, lvl, mask) for lvl in LEVELS}
    return images, eligibility(stack, dtype, input_plds, mask, zooms), corrected.affine


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--harness', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, required=True, help='The ASL JSON sidecar.')
    parser.add_argument('--mask', type=Path)
    parser.add_argument('--out', type=Path, help='Directory for the images (default: --harness).')
    args = parser.parse_args(argv)
    mask = np.asanyarray(nb.load(args.mask).dataobj).astype(bool) if args.mask else None
    images, record, affine = from_harness(
        args.harness, json.loads(args.metadata.read_text()), mask
    )
    out = args.out or args.harness
    out.mkdir(parents=True, exist_ok=True)
    for lvl, img in images.items():
        if lvl != 'masked':  # the same image as winsor
            nb.Nifti1Image(img.astype(np.float32), affine).to_filename(
                out / f'perfusion_{lvl}.nii.gz'
            )
    (out / 'perfusion.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
