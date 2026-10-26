#!/usr/bin/env python3
"""Why does registering to a blurred reference move an image? (Task 3)

Register an image to Gaussian-blurred copies of itself with MCFLIRT. There is no contrast
difference and no motion, so any transform is bias. Factors:

- blur: FWHM in mm (0 is a control);
- cost: MCFLIRT cost functions;
- field of view: as acquired; with the outermost slices dropped (edge slices that the blur
  fills from empty space); zero-padded by four slices; and shifted so the brain is centred;
- blur kernel: in-plane only, through-plane only, or isotropic; and ``both``, where the moving
  image is blurred like the reference (a candidate remedy for the GE path's smoothed
  reference).

Each condition reports the RMS displacement (80 mm sphere about the field-of-view centre) and
the translation along each axis, so a systematic drift towards one end of the field of view is
visible.

Usage (in the ASLPrep image, with FSL)::

    python blur_bias.py --image <3D image> --out <dir>
"""

import argparse
import itertools
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import nibabel as nb
import numpy as np
from scipy.ndimage import gaussian_filter

COSTS = ('mutualinfo', 'normcorr', 'corratio', 'leastsquares')
FWHMS = (0.0, 2.5, 5.0, 8.0)


def fsl_rms(mat, img):
    zooms = np.array(img.header.get_zooms()[:3])
    c = (np.array(img.shape[:3]) - 1) / 2 * zooms
    d = mat - np.eye(4)
    m, t = d[:3, :3], d[:3, 3] + d[:3, :3] @ c
    return float(np.sqrt(80.0**2 / 5 * np.trace(m.T @ m) + t @ t)), t.tolist()


def variants(img):
    data = np.asanyarray(img.dataobj, np.float32)
    data = data[..., 0] if data.ndim == 4 else data
    zooms = np.array(img.header.get_zooms()[:3])
    out = {'as_acquired': (data, img.affine)}
    out['drop_edge_slices'] = (data[..., 2:-2], img.affine @ _shift(0, 0, 2))
    pad = np.zeros(data.shape[:2] + (data.shape[2] + 8,), np.float32)
    pad[..., 4:-4] = data
    out['zero_padded'] = (pad, img.affine @ _shift(0, 0, -4))
    # Recentre: move the intensity centre of mass to the grid centre (integer voxel roll).
    com = np.array(np.nonzero(data > np.percentile(data, 80))).mean(axis=1)
    roll = np.round((np.array(data.shape) - 1) / 2 - com).astype(int)
    out['recentred'] = (np.roll(data, roll, axis=(0, 1, 2)), img.affine)
    return out, zooms


def _shift(i, j, k):
    m = np.eye(4)
    m[:3, 3] = [i, j, k]
    return m


def run(case, work):
    name, cost, fwhm, kernel, data, affine, zooms = case
    tag = f'{name}_{cost}_{fwhm:g}_{kernel}'
    d = work / tag
    d.mkdir(parents=True, exist_ok=True)
    sigma_mm = fwhm / 2.3548
    sig = np.array([sigma_mm] * 3) / zooms
    if kernel == 'inplane':
        sig[2] = 0
    elif kernel == 'throughplane':
        sig[:2] = 0
    ref = gaussian_filter(data, sig) if fwhm else data
    moving = ref if kernel == 'both' else data
    nb.Nifti1Image(ref, affine).to_filename(d / 'ref.nii.gz')
    nb.Nifti1Image(moving[..., None], affine).to_filename(d / 'in.nii.gz')
    cmd = ['mcflirt', '-in', str(d / 'in.nii.gz'), '-reffile', str(d / 'ref.nii.gz')]
    cmd += ['-out', str(d / 'out'), '-cost', cost, '-mats']
    subprocess.run(cmd, check=True, capture_output=True)
    mat = np.loadtxt(next((d / 'out.mat').glob('MAT_*')))
    rms, t = fsl_rms(mat, nb.Nifti1Image(data, affine))
    return {'fov': name, 'cost': cost, 'fwhm_mm': fwhm, 'kernel': kernel, 'rms_mm': rms, 't_mm': t}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--image', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    img = nb.load(args.image)
    fovs, zooms = variants(img)
    cases = []
    for (name, (data, affine)), cost, fwhm in itertools.product(fovs.items(), COSTS, FWHMS):
        kernels = ('iso', 'both')
        if name == 'as_acquired' and fwhm == 5.0:
            kernels += ('inplane', 'throughplane')
        cases += [(name, cost, fwhm, k, data, affine, zooms) for k in kernels]
    args.out.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(8) as pool:
        results = list(pool.map(lambda c: run(c, args.out / 'work'), cases))
    (args.out / 'blur_bias.json').write_text(json.dumps(results, indent=1) + '\n')
    for r in sorted(results, key=lambda r: (r['fov'], r['cost'], r['kernel'], r['fwhm_mm'])):
        t = ' '.join(f'{v:+.2f}' for v in r['t_mm'])
        print(
            f'{r["fov"]:17s} {r["cost"]:12s} {r["kernel"]:12s} fwhm {r["fwhm_mm"]:4.1f}  '
            f'rms {r["rms_mm"]:.3f}  t [{t}]'
        )


if __name__ == '__main__':
    main()
