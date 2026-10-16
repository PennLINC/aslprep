#!/usr/bin/env python3
"""Negative controls for the mean delta-M energy (``real_score.py``; Codex review, finding 6).

On the still fixtures, the mean delta-M of the raw pairs is compared with the same pairs after
a common shift applied to control and label alike (quarter-, half- and one-voxel, cubic
interpolation). A common shift creates no control-label misregistration, so an ideal metric
would not change; any change is interpolation smoothing of the noise and perfusion moving
relative to the fixed mask. Repeated with a fixed 1-voxel Gaussian pre-smoothing.

Result (2026-10-16): the energy falls by up to 4.5 % without pre-smoothing and up to 2 % with
it, while the short-listed candidates differ from production by about 0.2 % on real data. The
decision margin of ``real_score.py`` (``ENERGY_MARGIN``) is set above these effects.

Usage::

    python energy_controls.py
"""

import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, shift

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perfusion import epi_mask


def energy(img, mask, sigma):
    """Mean squared value inside ``mask``, after an optional Gaussian pre-smoothing."""
    img = gaussian_filter(img, sigma) if sigma else img
    return float((img[mask] ** 2).mean())


rows = []
for name in ('pcasl1pld', 'pcasl1pld_grase', 'pasl_multipld', 'pcasl_multipld'):
    fx = Path('/home/tsalo/aslscan-data-new/aslscan') / name / 'sub-01' / 'perf'
    asl = nb.load(next(fx.glob('*_asl.nii.gz')))
    d = np.asanyarray(asl.dataobj, np.float32)
    ctx = pd.read_table(next(fx.glob('*_aslcontext.tsv')))['volume_type'].to_numpy()
    c, lab = np.flatnonzero(ctx == 'control'), np.flatnonzero(ctx == 'label')
    mean_dm = (d[..., c] - d[..., lab]).mean(-1)
    # Energy is computed per pair-mean; the shift acts on the raw volumes, and by linearity on
    # the difference. Noise lives in each pair, so shift the pair differences, then average.
    pairs = d[..., c] - d[..., lab]
    ref = Path('/home/tsalo/hmc-probe/fixture-scores/c0') / name / 'aslref.nii.gz'
    mask = epi_mask(ref)
    for sigma in (0.0, 1.0):
        e0 = energy(mean_dm, mask, sigma)
        for label, s in (
            ('half-voxel shift', (0.5, 0.5, 0.0)),
            ('quarter-voxel', (0.25, 0.25, 0.25)),
            ('one-voxel shift', (1.0, 0.0, 0.0)),
        ):
            moved = np.stack(
                [shift(pairs[..., k], s, order=3, mode='nearest') for k in range(pairs.shape[-1])],
                -1,
            ).mean(-1)
            change = 100 * (energy(moved, mask, sigma) / e0 - 1)
            rows.append(
                {
                    'fixture': name,
                    'presmooth_vox': sigma,
                    'control': label,
                    'energy_change_pct': change,
                }
            )
df = pd.DataFrame(rows)
pd.set_option('display.width', 200)
print(
    df.pivot_table(
        index=['presmooth_vox', 'control'], columns='fixture', values='energy_change_pct'
    )
    .round(2)
    .to_string()
)
