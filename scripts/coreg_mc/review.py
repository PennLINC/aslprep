#!/usr/bin/env python3
"""Blinded review figures and batches (plan: "Blinded review figures").

One figure per pair of transforms on one run. The background is the run's own ASL images at
their native resolution (the conventional reference on top, the perfusion image P-robust below),
only reoriented to the closest canonical axes (a flip and permutation of voxel axes, never a
resampling). Over them, contours carried from T1w space through each side's pull (T1w point ->
aslref point): the white-matter boundary (the 0.5 isocontour of the run's WM mask) as a solid
line and the T1w brain mask's outline as a dashed line. Seven slices (three axial, two coronal,
two sagittal) are placed from the ASL brain mask's extent, which no arm changes, and are the
same for every figure of a run; each panel is cropped to that extent plus a margin.

Rows: side A, side B (one colour, the same for both), and an overlay of the two in two colours
assigned at random per figure. The intensity window is the 1-99th percentile inside the union of
the two carried brain masks. Each figure also has two flicker frames (``flicker``): the
reference and perfusion rows with only A's or only B's contours, in one colour, which the page
alternates so a small shift shows as motion.

A batch (``build_batch``) holds comparison pairs, about 10 % of them repeated later with sides
swapped, and known-answer pairs (one side an arm's solution, the other the same solution moved
by a rigid error of 1, 2 or 4 mm RMS over the brain). Sides and order are random from a recorded
seed. The page loads ``batch.json`` and the figures only; ``key.json`` (which run, which arm on
which side) is written beside them and never loaded by the page.

Usage: ``build_batch`` from Python (the pairs come from the scorer's tables); ``--dry-run``
builds a batch of every pair of the given arms on one work directory and checks it for leaks::

    python scripts/coreg_mc/review.py --dry-run --work <work dir> --runs <ASLPrep runs> \\
        --manifest scripts/coreg_mc/manifest.tsv --out <dir>
"""

import argparse
import json
import secrets
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
from scipy.ndimage import map_coordinates

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import scoring  # noqa: E402

#: Slice positions as fractions of the ASL brain mask's extent along each canonical axis.
SLICES = (('axial', 2, 0.3), ('axial', 2, 0.5), ('axial', 2, 0.7))
SLICES += (('coronal', 1, 0.4), ('coronal', 1, 0.6), ('sagittal', 0, 0.35), ('sagittal', 0, 0.65))
#: Contour upsampling within a slice (points per ASL voxel side).
UPSAMPLE = 4
SIDE_COLOUR = '#ffd400'
OVERLAY_COLOURS = ('#ff3bd4', '#00e676')
LINE_WIDTH = 1.2
#: Each panel is cropped to the ASL brain mask's extent plus this margin (mm), so contours
#: carried a few millimetres outside the mask stay in view.
CROP_MARGIN_MM = 8.0
#: Panel size (inches) and resolution: the page fits the figure to its width, and its
#: click-to-zoom shows it at full resolution.
PANEL_INCHES = 2.0
DPI = 200
KNOWN_MM = (1.0, 2.0, 4.0)
REPEAT_FRACTION = 0.1


class RunImages:
    """One run's ASL images (canonical axes) and T1w-space structures, from a work directory."""

    def __init__(self, work, perfusion='P-robust'):
        self.work = Path(work)
        ref = nb.as_closest_canonical(nb.load(self.work / 'conventional.nii.gz'))
        perf = nb.as_closest_canonical(nb.load(self.work / f'{perfusion}.nii.gz'))
        mask = nb.as_closest_canonical(nb.load(self.work / 'mask.nii.gz'))
        if not np.allclose(ref.affine, perf.affine, atol=1e-3):
            raise ValueError(f'{work}: the perfusion image is not on the reference grid')

        self.affine = ref.affine
        self.zooms = np.sqrt((ref.affine[:3, :3] ** 2).sum(0))
        self.images = [np.asanyarray(ref.dataobj, float), np.asanyarray(perf.dataobj, float)]
        self.asl_mask = np.asanyarray(mask.dataobj) > 0
        wm = nb.load(self.work / 'wm.nii.gz')
        t1_mask = nb.load(self.work / 't1w_mask.nii.gz')
        self.t1_affine = wm.affine
        self.wm = np.asanyarray(wm.dataobj, float)
        self.t1_mask = np.asanyarray(t1_mask.dataobj, float)
        idx = np.argwhere(self.asl_mask)
        self.extent = idx.min(0), idx.max(0)
        self.slices = self._place()
        self.brain_points = scoring.mask_points(self.work / 't1w_mask.nii.gz', every=20)

    def _place(self):
        """[(view, axis, index)] from the ASL brain mask's extent along each canonical axis."""
        lo, hi = self.extent
        span = hi - lo
        return [(view, axis, round(lo[axis] + f * span[axis])) for view, axis, f in SLICES]

    def crop(self, axis):
        """(xlim, ylim) in voxel coordinates for a slice across ``axis``: the ASL brain mask's
        extent in the slice's plane plus ``CROP_MARGIN_MM``. No arm changes it."""
        lo, hi = self.extent
        limits = []
        for a in (a for a in range(3) if a != axis):
            margin = 0.5 + CROP_MARGIN_MM / self.zooms[a]
            limits.append((lo[a] - margin, hi[a] + margin))
        return limits

    def plane(self, axis, index):
        """Fractional voxel coordinates (3, n, m) of an upsampled grid over one slice."""
        shape = list(self.images[0].shape[:3])
        inplane = [a for a in range(3) if a != axis]
        u = np.arange(shape[inplane[0]] * UPSAMPLE) / UPSAMPLE - 0.5 + 0.5 / UPSAMPLE
        v = np.arange(shape[inplane[1]] * UPSAMPLE) / UPSAMPLE - 0.5 + 0.5 / UPSAMPLE
        uu, vv = np.meshgrid(u, v, indexing='ij')
        coords = np.zeros((3, *uu.shape))
        coords[inplane[0]], coords[inplane[1]], coords[axis] = uu, vv, index
        return coords, uu, vv

    def carried(self, pull, axis, index):
        """WM and T1w-mask values on one slice's upsampled grid, carried through ``pull``."""
        coords, uu, vv = self.plane(axis, index)
        ijk = np.vstack([coords.reshape(3, -1), np.ones((1, coords[0].size))])
        world = self.affine @ ijk
        t1w = np.linalg.inv(pull) @ world
        t1_ijk = (np.linalg.inv(self.t1_affine) @ t1w)[:3]
        sample = [
            map_coordinates(v, t1_ijk, order=1, cval=0.0).reshape(uu.shape)
            for v in (self.wm, self.t1_mask)
        ]
        return uu, vv, sample

    def carried_mask_voxels(self, pull):
        """The T1w brain mask carried onto the ASL grid (nearest voxel centres), as a boolean."""
        ijk = np.argwhere(np.ones(self.images[0].shape[:3], bool))
        world = self.affine @ np.c_[ijk, np.ones(len(ijk))].T
        t1_ijk = (np.linalg.inv(self.t1_affine) @ np.linalg.inv(pull) @ world)[:3]
        values = map_coordinates(self.t1_mask, t1_ijk, order=0, cval=0.0)
        return values.reshape(self.images[0].shape[:3]) > 0.5


def _section(volume, axis, index):
    """A 2-D slice with the first in-plane axis horizontal."""
    return np.take(volume, index, axis=axis)


def _backgrounds(run, pull_a, pull_b):
    """The two background images (3-D) and their windows, the same for every rendering of a
    pair: the 1-99th percentile inside the union of the two carried brain masks."""
    union = run.carried_mask_voxels(pull_a) | run.carried_mask_voxels(pull_b)
    images = [im[..., 0] if im.ndim == 4 else im for im in run.images]
    windows = []
    for image in images:
        inside = image[union] if union.any() else image.ravel()
        windows.append(np.percentile(inside, [1, 99]))

    return images, windows


def _layout(run, rows, bottom_inches):
    """A figure of ``rows`` x slices panels, each column as wide as its cropped slice is (in
    mm) for a panel height of ``PANEL_INCHES``, so the panels fill it without gaps."""
    ratios = []
    for _view, axis, _index in run.slices:
        (x0, x1), (y0, y1) = run.crop(axis)
        u, v = (a for a in range(3) if a != axis)
        ratios.append((x1 - x0) * run.zooms[u] / ((y1 - y0) * run.zooms[v]))

    width = PANEL_INCHES * (sum(ratios) + 0.4)
    height = PANEL_INCHES * rows + bottom_inches
    fig, axes = plt.subplots(
        rows,
        len(ratios),
        figsize=(width, height),
        facecolor='black',
        gridspec_kw={'width_ratios': ratios},
        squeeze=False,
    )
    fig.subplots_adjust(
        left=PANEL_INCHES * 0.4 / width,
        right=0.995,
        top=1 - 0.01 * PANEL_INCHES / height,
        bottom=bottom_inches / height,
        wspace=0.02,
        hspace=0.02,
    )
    return fig, axes


def _panel(ax, run, images, windows, image_i, axis, index, contours, sides, title=None):
    """One slice of one background image with the given sides' contours, cropped."""
    lo, hi = windows[image_i]
    inplane = [a for a in range(3) if a != axis]
    ax.imshow(
        _section(images[image_i], axis, index).T,
        cmap='gray',
        vmin=lo,
        vmax=hi,
        origin='lower',
        interpolation='nearest',
        aspect=run.zooms[inplane[1]] / run.zooms[inplane[0]],
    )
    for side, colour in sides:
        uu, vv, (wm, mask) = contours[side]
        kwargs = {'colors': [colour], 'linewidths': LINE_WIDTH}
        if wm.max() > 0.5 > wm.min():
            ax.contour(uu, vv, wm, levels=[0.5], linestyles='solid', **kwargs)

        if mask.max() > 0.5 > mask.min():
            ax.contour(uu, vv, mask, levels=[0.5], linestyles='dashed', **kwargs)

    xlim, ylim = run.crop(axis)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_axis_off()
    if title:
        ax.text(
            -0.06,
            0.5,
            title,
            color='white',
            fontsize=14,
            ha='right',
            va='center',
            transform=ax.transAxes,
        )


def _save(fig, out, palette=None):
    """Save as a 256-colour palette PNG (a third of the size of RGBA); return the image, whose
    palette the second flicker frame reuses so the two frames' backgrounds match exactly."""
    import io

    from PIL import Image

    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=DPI, facecolor='black', metadata={'Software': None})
    plt.close(fig)
    image = Image.open(buf).convert('RGB')
    if palette is None:
        image = image.quantize(256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    else:
        image = image.quantize(palette=palette, dither=Image.Dither.NONE)

    image.save(out, 'PNG', optimize=True)
    return image


def _contours(run, pull_a, pull_b):
    """{slice: {side: carried contours}} for every slice of the run."""
    return {
        (axis, index): {
            side: run.carried(pull, axis, index) for side, pull in (('A', pull_a), ('B', pull_b))
        }
        for _view, axis, index in run.slices
    }


def figure(run, pull_a, pull_b, out, colours):
    """Write one blinded figure to ``out`` (PNG, no metadata)."""
    images, windows = _backgrounds(run, pull_a, pull_b)
    contours = _contours(run, pull_a, pull_b)
    fig, axes = _layout(run, 6, 0.2 * PANEL_INCHES)
    row_titles = ('A', 'A', 'B', 'B', 'A + B', 'A + B')
    for col, (_view, axis, index) in enumerate(run.slices):
        for row in range(6):
            if row < 4:
                sides = [('A' if row < 2 else 'B', SIDE_COLOUR)]
            else:
                sides = [('A', colours[0]), ('B', colours[1])]

            title = row_titles[row] if col == 0 else None
            c = contours[axis, index]
            _panel(axes[row, col], run, images, windows, row % 2, axis, index, c, sides, title)

    # The overlay's colour key, as coloured labels.
    for x, side, colour in ((0.45, 'A', colours[0]), (0.55, 'B', colours[1])):
        fig.text(x, 0.006, f'{side} (A + B rows)', color=colour, ha='center', fontsize=12)

    _save(fig, out)


def flicker(run, pull_a, pull_b, out_a, out_b):
    """Write the two flicker frames of a pair: the reference and perfusion rows with one side's
    contours each, in one colour, on identical backgrounds and windows, so the frames differ only
    where the contours do. The page alternates them and labels which side is showing."""
    images, windows = _backgrounds(run, pull_a, pull_b)
    contours = _contours(run, pull_a, pull_b)
    palette = None
    for side, out in (('A', out_a), ('B', out_b)):
        fig, axes = _layout(run, 2, 0.0)
        for col, (_view, axis, index) in enumerate(run.slices):
            for row in range(2):
                c = contours[axis, index]
                sides = [(side, SIDE_COLOUR)]
                _panel(axes[row, col], run, images, windows, row, axis, index, c, sides)

        palette = _save(fig, out, palette)


def known_answer(run, pull, size_mm, rng):
    """``pull`` with a random rigid error of ``size_mm`` RMS over the brain (T1w side)."""
    centre = run.brain_points.mean(0)
    error = scoring.scaled_rigid(rng.normal(size=6), size_mm, run.brain_points, centre)
    return pull @ error


def build_batch(pairs, known, out_dir, seed, repeat_fraction=REPEAT_FRACTION):
    """Write a batch: figures, ``batch.json`` (what the page loads) and ``key.json``.

    ``pairs``: dicts with ``run`` (a label for the key), ``work``, ``arm_a``, ``arm_b``, ``pull_a``
    and ``pull_b``. ``known``: dicts with ``run``, ``work``, ``arm``, ``pull`` and ``size_mm``.
    """
    rng = np.random.default_rng(seed)
    out_dir = Path(out_dir)
    (out_dir / 'figures').mkdir(parents=True, exist_ok=True)
    items = [{'kind': 'comparison', **p} for p in pairs]
    for k in known:
        run = RunImages(k['work'])
        moved = known_answer(run, k['pull'], k['size_mm'], rng)
        items.append(
            {
                'kind': 'known',
                'run': k['run'],
                'work': k['work'],
                'arm_a': k['arm'],
                'arm_b': f'{k["arm"]} moved {k["size_mm"]:g} mm',
                'pull_a': k['pull'],
                'pull_b': moved,
                'size_mm': k['size_mm'],
            },
        )

    order = list(rng.permutation(len(items)))
    comparisons = [i for i in order if items[i]['kind'] == 'comparison']
    n_repeat = round(repeat_fraction * len(comparisons))
    repeats = list(rng.choice(comparisons, size=n_repeat, replace=False)) if n_repeat else []
    for original in repeats:
        # Later in the batch than the original: somewhere in the second half after it.
        pos = order.index(original)
        lo = max(pos + 1, len(order) // 2)
        order.insert(int(rng.integers(lo, len(order) + 1)), ('repeat', original))

    batch, key, images, first = [], {}, {}, {}
    for entry in order:
        repeat_of = entry[1] if isinstance(entry, tuple) else None
        index = entry if repeat_of is None else repeat_of
        item = items[index]
        fid = secrets.token_hex(8)
        swap = bool(rng.integers(2))
        if repeat_of is not None:
            # Sides swapped relative to the original's presentation.
            swap = not key[first[repeat_of]]['swapped']

        left, right = ('b', 'a') if swap else ('a', 'b')
        colours = list(OVERLAY_COLOURS)
        if rng.integers(2):
            colours.reverse()

        work = str(item['work'])
        if work not in images:
            images[work] = RunImages(work)

        figure(
            images[work],
            item[f'pull_{left}'],
            item[f'pull_{right}'],
            out_dir / 'figures' / f'{fid}.png',
            colours,
        )
        frames = [f'figures/{secrets.token_hex(8)}.png' for _ in 'AB']
        flicker(
            images[work],
            item[f'pull_{left}'],
            item[f'pull_{right}'],
            *(out_dir / f for f in frames),
        )
        batch.append({'id': fid, 'image': f'figures/{fid}.png', 'flicker': frames})
        key[fid] = {
            'kind': 'repeat' if repeat_of is not None else item['kind'],
            'run': item['run'],
            'arm_A': item[f'arm_{left}'],
            'arm_B': item[f'arm_{right}'],
            'swapped': swap,
            'repeat_of': first[repeat_of] if repeat_of is not None else None,
            'size_mm': item.get('size_mm'),
            'overlay_colours': {'A': colours[0], 'B': colours[1]},
        }
        if repeat_of is None:
            first[index] = fid

    (out_dir / 'batch.json').write_text(json.dumps({'items': batch}, indent=1) + '\n')
    (out_dir / 'key.json').write_text(json.dumps({'seed': seed, 'figures': key}, indent=1) + '\n')
    return batch, key


def png_text(path):
    """The text a PNG carries outside its pixels: every chunk except IHDR, IDAT and IEND.

    Pixel data is compressed and can contain any short byte string by chance, so it is not
    searched; what a figure shows is checked by the figure's construction (no text but the side
    labels and the colour key).
    """
    data = Path(path).read_bytes()
    pos, out = 8, b''
    while pos < len(data):
        length = int.from_bytes(data[pos : pos + 4], 'big')
        kind = data[pos + 4 : pos + 8]
        if kind not in (b'IHDR', b'IDAT', b'IEND'):
            out += kind + data[pos + 8 : pos + 8 + length]

        pos += 12 + length
    return out


def leaks(out_dir, words):
    """Words from the key found in what the page loads: ``batch.json`` and the figures' text."""
    out_dir = Path(out_dir)
    texts = [(out_dir / 'batch.json').read_bytes()]
    texts += [png_text(f) for f in sorted((out_dir / 'figures').glob('*.png'))]
    found = set()
    for text in texts:
        found |= {w for w in words if w and w.lower().encode() in text.lower()}
    return sorted(found)


def main(argv=None):
    import itertools

    import run_real

    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--dry-run', action='store_true', required=True)
    p.add_argument('--work', type=Path, required=True, help='One run_real work directory.')
    p.add_argument('--out', type=Path, required=True)
    default_arms = ['FSL conventional default', 'FSL P-winsor default']
    p.add_argument('--arms', nargs='+', default=default_arms)
    p.add_argument('--seed', type=int, default=20261029)
    args = p.parse_args(argv)
    found, _ = run_real.collect(args.work, [])
    arms = {a: found[a] for a in args.arms}
    label = args.work.name
    pairs = [
        {
            'run': label,
            'work': args.work,
            'arm_a': a,
            'arm_b': b,
            'pull_a': arms[a],
            'pull_b': arms[b],
        }
        for a, b in itertools.combinations(arms, 2)
    ]
    first = next(iter(arms))
    known = [
        {'run': label, 'work': args.work, 'arm': first, 'pull': arms[first], 'size_mm': size}
        for size in KNOWN_MM
    ]
    batch, key = build_batch(pairs, known, args.out, args.seed, repeat_fraction=1.0)
    words = {label, *label.split('__'), 'FSL', 'ANTs', 'BBR', 'winsor', 'conventional', 'moved'}
    words |= {w for k in key.values() for w in (k['arm_A'], k['arm_B'], k['run'])}
    found_words = leaks(args.out, words)
    print(json.dumps({'figures': len(batch), 'kinds': [k['kind'] for k in key.values()]}))
    print('LEAKS: ' + ', '.join(found_words) if found_words else 'NO LEAKS')
    return 1 if found_words else 0


if __name__ == '__main__':
    sys.exit(main())
