#!/usr/bin/env python3
"""Draw each run's masks over its ASL reference and mean perfusion image.

One PNG per job: the aslref (top row) and the mean perfusion image, or the true brain for
simulated runs (bottom row), in three axial, one coronal and one sagittal slice, with the
SynthStrip, legacy and anatomical (``anat_pv50``) mask outlines.

Usage::

    python montages.py --jobs jobs.json --masks masks/ --results results/ --out figures/
"""

import argparse
import json
from pathlib import Path

import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np

#: Categorical slots 1-3 of the report palette, in fixed order.
COLORS = {'synthstrip': '#2a78d6', 'legacy': '#eb6834', 'anat_pv50': '#1baf7a'}
LABELS = {
    'synthstrip': 'SynthStrip',
    'legacy': 'Legacy (BET + 3dAutomask)',
    'anat_pv50': 'Anatomical',
}


def _canonical(img):
    return nb.as_closest_canonical(img)


def _slices(mask):
    """Three axial slices through the brain, and the central coronal and sagittal slices."""
    idx = np.argwhere(mask)
    lo, hi = idx.min(0), idx.max(0)
    mid = (lo + hi) // 2
    z = np.linspace(lo[2] + 0.2 * (hi[2] - lo[2]), lo[2] + 0.8 * (hi[2] - lo[2]), 3).astype(int)
    return [('z', k) for k in z] + [('y', mid[1]), ('x', mid[0])]


def _cut(data, axis, k):
    return np.rot90(np.take(data, k, axis='xyz'.index(axis)))


def _aspect(zooms, axis):
    zx, zy, zz = zooms
    return {'z': zy / zx, 'y': zz / zx, 'x': zz / zy}[axis]


def draw(job, masks, bottom, bottom_label, out):
    ref = _canonical(nb.load(job['ref']))
    zooms = ref.header.get_zooms()[:3]
    ref_d = np.asanyarray(ref.dataobj, np.float32)
    mdata = {k: np.asanyarray(_canonical(nb.load(p)).dataobj) > 0 for k, p in masks.items()}
    bot = np.asanyarray(_canonical(bottom).dataobj, np.float32)
    cuts = _slices(mdata['synthstrip'] | mdata['anat_pv50'])

    fig, axes = plt.subplots(2, len(cuts), figsize=(2.2 * len(cuts), 4.9), facecolor='#1a1a19')
    for row, (img, label) in enumerate(((ref_d, 'aslref'), (bot, bottom_label))):
        finite = img[np.isfinite(img)]
        vmin, vmax = np.percentile(finite, [1, 99.5]) if label != 'truth' else (0, 1)
        if label == 'perfusion':
            vmin = 0
        for col, (axis, k) in enumerate(cuts):
            ax = axes[row, col]
            ax.imshow(
                _cut(img, axis, k),
                cmap='gray',
                vmin=vmin,
                vmax=vmax,
                aspect=_aspect(zooms, axis),
                interpolation='nearest',
            )
            for name, m in mdata.items():
                ax.contour(
                    _cut(m.astype(float), axis, k),
                    levels=[0.5],
                    colors=COLORS[name],
                    linewidths=1.2,
                )
            ax.set_axis_off()
            if col == 0:
                ax.text(
                    0.02,
                    0.98,
                    label,
                    transform=ax.transAxes,
                    va='top',
                    color='#f0efec',
                    fontsize=8,
                )
    handles = [plt.Line2D([], [], color=COLORS[k], lw=2, label=LABELS[k]) for k in COLORS]
    fig.legend(
        handles=handles,
        loc='lower center',
        ncol=3,
        frameon=False,
        fontsize=8,
        labelcolor='#f0efec',
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.07, wspace=0.03, hspace=0.03)
    fig.savefig(out, dpi=80, facecolor=fig.get_facecolor())
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--jobs', type=Path, required=True)
    parser.add_argument('--masks', type=Path, required=True)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    for job in json.loads(args.jobs.read_text()):
        anat = args.results / 'masks' / job['name'] / 'anat_pv50.nii.gz'
        if not anat.exists():
            continue
        masks = {
            'synthstrip': args.masks / job['name'] / 'synthstrip.nii.gz',
            'legacy': args.masks / job['name'] / 'legacy.nii.gz',
            'anat_pv50': anat,
        }
        if job['kind'] == 'real':
            bottom, label = nb.load(job['perfusion']), 'perfusion'
            sign = 1
            data = np.asanyarray(bottom.dataobj, np.float32)
            union = np.asanyarray(nb.load(anat).dataobj) > 0
            sign = 1 if np.median(data[union]) >= 0 else -1
            bottom = nb.Nifti1Image(sign * data, bottom.affine)
        else:
            bottom, label = nb.load(args.results / 'masks' / job['name'] / 'truth.nii.gz'), 'truth'
        draw(job, masks, bottom, label, args.out / f'{job["name"]}.png')
        print(job['name'], flush=True)


if __name__ == '__main__':
    main()
