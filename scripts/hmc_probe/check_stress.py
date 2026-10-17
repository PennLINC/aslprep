#!/usr/bin/env python3
"""Check that each Task 4 stress fixture shows its effect in the simulated images.

The plan requires each stress case to be checked before it is used (plan: "Stress cases"):

- ``st_vessels``: the mean delta-M along the arteries (phantom ``abv`` above 0.05) against
  grey matter's (the fixture's grey-matter partial volume at least 0.7, away from the arteries);
- ``st_noperf``: the mean delta-M in the territory against the mirror-image region;
- ``st_spiral``: the signal outside the head (the phantom's tissue labels dilated by 2 voxels,
  on the scan grid) as a fraction of the signal inside, in the M0 and control images.

Usage::

    python check_stress.py --fixtures <dir>/aslscan
"""

import argparse
import json
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
from nibabel.processing import resample_from_to
from scipy.ndimage import binary_dilation

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from fixtures import NO_PERFUSION  # noqa: E402


def _series(fx):
    perf = fx / 'sub-01' / 'perf'
    asl = nb.load(next(perf.glob('*_asl.nii.gz')))
    ctx = pd.read_table(next(perf.glob('*_aslcontext.tsv')))['volume_type'].to_numpy()
    return asl, np.asanyarray(asl.dataobj, np.float32), ctx


def _on_grid(phantom_file, ref, order=1):
    img = nb.load(phantom_file)
    return np.asanyarray(resample_from_to(img, (ref.shape[:3], ref.affine), order=order).dataobj)


def _mean_dm(data, ctx):
    return (data[..., ctx == 'control'] - data[..., ctx == 'label']).mean(-1)


def check_vessels(fx, phantom):
    asl, data, ctx = _series(fx)
    dm = _mean_dm(data, ctx)
    abv = _on_grid(phantom / 'abv.nii.gz', asl)
    gm_pv = np.asanyarray(
        nb.load(fx / 'sub-01' / 'perf' / 'ground-truth' / 'sub-01_desc-pvGM_gt.nii.gz').dataobj
    )
    vessel = abv > 0.05
    gm = (gm_pv >= 0.7) & ~binary_dilation(vessel, iterations=2)
    return {
        'vessel_voxels': int(vessel.sum()),
        'vessel_dm': float(dm[vessel].mean()),
        'gm_dm': float(dm[gm].mean()),
        'ratio': float(dm[vessel].mean() / dm[gm].mean()),
    }


def check_noperf(fx):
    asl, data, ctx = _series(fx)
    dm = _mean_dm(data, ctx)
    ijk = np.indices(asl.shape[:3]).reshape(3, -1).T
    xyz = nb.affines.apply_affine(asl.affine, ijk)
    centre, radius = np.asarray(NO_PERFUSION[0]), NO_PERFUSION[1]
    mirror = centre * np.array([-1.0, 1.0, 1.0])
    gm_pv = np.asanyarray(
        nb.load(fx / 'sub-01' / 'perf' / 'ground-truth' / 'sub-01_desc-pvGM_gt.nii.gz').dataobj
    )
    gm = (gm_pv >= 0.5).ravel()
    inner = radius - 8.0  # away from the territory's edge
    terr = (np.linalg.norm(xyz - centre, axis=1) <= inner) & gm
    contra = (np.linalg.norm(xyz - mirror, axis=1) <= inner) & gm
    flat = dm.ravel()
    return {
        'territory_voxels': int(terr.sum()),
        'territory_dm': float(flat[terr].mean()),
        'mirror_dm': float(flat[contra].mean()),
        'ratio': float(flat[terr].mean() / flat[contra].mean()),
    }


def check_spiral(fx, phantom):
    asl, data, ctx = _series(fx)
    dseg = _on_grid(phantom / 'dseg.nii.gz', asl, order=0)
    head = binary_dilation(dseg > 0, iterations=2)
    out = {}
    for kind in ('m0scan', 'control'):
        vol = data[..., ctx == kind].mean(-1)
        out[f'{kind}_outside_over_inside'] = float(
            np.abs(vol[~head]).mean() / np.abs(vol[head]).mean()
        )
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--fixtures', type=Path, required=True)
    args = parser.parse_args(argv)
    phantoms = args.fixtures / '_phantoms'
    report = {
        'st_vessels': check_vessels(args.fixtures / 'st_vessels', phantoms / 'tfmni-b0-vessels'),
        'st_noperf': check_noperf(args.fixtures / 'st_noperf'),
        'st_spiral': check_spiral(args.fixtures / 'st_spiral', phantoms / 'tfmni-b0'),
    }
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
