#!/usr/bin/env python3
"""Metric 5: a full run's CBF (and ATT) against the oracle run's, on a simulated fixture.

The target is the same quantification of the same data, with the true poses as the HMC
transforms (the ``oracle`` candidate, injected into a full ASLPrep run). Model, calibration and
masks are the same, so what differs is motion correction and what follows from it. Error
against the simulator's physical truth is a different question and is not scored here.

Frames as in ``score.py``: a run's ASL-space image ``y`` shows ``static(F^-1 y)``, with ``F`` the
polar mean of ``H_v^-1 P_v`` over its HMC transforms. The oracle's map is pulled into the
candidate's frame through ``F_o F^-1``, and both are compared inside the true brain seen in the
candidate's frame (one fixed support for every candidate, so an eroded mask cannot help). Grey
and white matter come from the truth's partial-volume maps (at least 0.5), the edge shell is the
outer 3 mm of the true brain.

Usage::

    python score_cbf.py --fixture <dir> --run <output dir> --hmc <xfm> \\
        --oracle-run <output dir> --oracle-hmc <xfm> [--out score_cbf.json]
"""

import argparse
import json
import sys
from pathlib import Path

import nibabel as nb
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

from score import mask_in_frame, polar_mean, shells  # noqa: E402
from score_full import _load_series  # noqa: E402

from aslprep.tests import truth_geometry as tg  # noqa: E402
from aslprep.tests import truth_scoring as ts  # noqa: E402

SCORE_CBF_VERSION = 1


def run_frame(fixture, hmc_file):
    """``F``: the polar mean of ``H_v^-1 P_v`` over the run's HMC transforms."""
    hmc = _load_series(hmc_file)
    poses = fixture.poses()
    keep = [
        i for i, t in enumerate(fixture.context) if t in ('control', 'label', 'm0scan', 'deltam')
    ]
    if len(keep) != len(hmc):
        keep = list(range(len(hmc)))
    return polar_mean([np.linalg.inv(h) @ poses[i] for h, i in zip(hmc, keep, strict=True)])


def pull_between_frames(src_img, src_frame, dst_img, dst_frame):
    """``src_img`` (frame ``src_frame``) resampled onto ``dst_img``'s grid in ``dst_frame``."""
    from scipy.ndimage import map_coordinates

    shape = dst_img.shape[:3]
    ijk = np.indices(shape).reshape(3, -1).T
    world = tg.apply_points(dst_img.affine, ijk)
    src = tg.apply_points(
        np.linalg.inv(src_img.affine) @ src_frame @ np.linalg.inv(dst_frame), world
    )
    data = np.asanyarray(src_img.dataobj, np.float64)
    data = data[..., 0] if data.ndim == 4 else data
    vals = map_coordinates(data, src.T, order=1, mode='constant', cval=np.nan)
    return vals.reshape(shape)


def _stats(diff, oracle, region, gm_mean):
    d = diff[region & np.isfinite(diff)]
    if not d.size:
        return {'n': 0}
    return {
        'n': int(d.size),
        'bias': float(d.mean()),
        'bias_rel_gm': float(d.mean() / gm_mean),
        'rmse': float(np.sqrt((d**2).mean())),
        'rmse_rel_gm': float(np.sqrt((d**2).mean()) / gm_mean),
        'p95_abs': float(np.percentile(np.abs(d), 95)),
        'oracle_mean': float(np.nanmean(oracle[region])),
    }


def _find(run, acq, suffix):
    hits = sorted((Path(run) / 'sub-01' / 'perf').glob(f'sub-01_acq-{acq}_{suffix}.nii.gz'))
    return hits[0] if len(hits) == 1 else None


def score_cbf(fixture_dir, run, hmc, oracle_run, oracle_hmc):
    fixture_dir = Path(fixture_dir)
    asl = next((fixture_dir / 'sub-01' / 'perf').glob('*_asl.nii.gz'))
    acq = asl.name.split('_acq-')[1].split('_')[0]
    fixture = ts.Fixture(fixture_dir, acq)
    frame = run_frame(fixture, hmc)
    frame_o = run_frame(fixture, oracle_hmc)
    gt = fixture_dir / 'sub-01' / 'perf' / 'ground-truth'
    out = {
        'score_cbf_version': SCORE_CBF_VERSION,
        'fixture': str(fixture_dir),
        'run': str(run),
        'oracle_run': str(oracle_run),
        'frame_offset_mm': tg.rms_displacement(frame, fixture.voxel_centers(fixture.brain)[1]),
        'oracle_frame_offset_mm': tg.rms_displacement(
            frame_o, fixture.voxel_centers(fixture.brain)[1]
        ),
    }
    for suffix in ('cbf', 'att'):
        cand_file, orac_file = _find(run, acq, suffix), _find(oracle_run, acq, suffix)
        if cand_file is None or orac_file is None:
            continue
        c_img = nb.load(cand_file)
        cand = np.asanyarray(c_img.dataobj, np.float64)
        cand = cand[..., 0] if cand.ndim == 4 else cand
        orac = pull_between_frames(nb.load(orac_file), frame_o, c_img, frame)
        brain = mask_in_frame(fixture.brain, fixture.affine, frame, c_img)
        regions = {'brain': brain, **shells(brain, c_img.affine)}
        for t in ('GM', 'WM'):
            pv_img = nb.load(gt / f'sub-01_desc-pv{t}_gt.nii.gz')
            pv = np.asanyarray(pv_img.dataobj) >= 0.5
            regions[t] = brain & mask_in_frame(pv, pv_img.affine, frame, c_img)
        gm_mean = float(np.nanmean(orac[regions['GM']]))
        diff = cand - orac
        res = {k: _stats(diff, orac, r, gm_mean) for k, r in regions.items()}
        if suffix == 'att':
            vals = cand[brain]
            lo, hi = np.nanmin(vals), np.nanmax(vals)
            res['fraction_at_bound'] = float(np.mean(np.isclose(vals, lo) | np.isclose(vals, hi)))
        out[suffix] = res
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--fixture', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--hmc', type=Path, required=True)
    parser.add_argument('--oracle-run', type=Path, required=True)
    parser.add_argument('--oracle-hmc', type=Path, required=True)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args(argv)
    text = json.dumps(
        score_cbf(args.fixture, args.run, args.hmc, args.oracle_run, args.oracle_hmc),
        indent=2,
        default=float,
    )
    if args.out:
        args.out.write_text(text + '\n')
    else:
        print(text)


if __name__ == '__main__':
    main()
