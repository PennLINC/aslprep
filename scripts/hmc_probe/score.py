#!/usr/bin/env python3
"""HMC probe scorer (Task 2): score one harness output on an aslscan fixture.

Frames follow ``aslprep/tests/truth_scoring.py`` (RAS point mappings):

- ``P_v``: the simulator's pose of volume ``v`` (static phantom -> moved);
- ``H_v``: the candidate's HMC transform, a pull mapping (aslref point -> volume-``v`` point),
  as written to ``hmc_xfms.txt``;
- ``A_v = H_v^-1 P_v``: the frame of corrected volume ``v``. The corrected image shows the
  static phantom as ``corrected_v(y) = static(A_v^-1 y)``; perfect relative correction makes
  every ``A_v`` the same;
- ``F``: the frame the reference image depicts. A separate M0 does not move in aslscan, and
  without motion every volume is static, so ``F`` is the identity then. With motion and a
  reference averaged from moved volumes, ``F`` has no single rigid value; the scorer then uses
  the polar mean of the selected volumes' poses and says so.

Metrics (numbers follow the plan's "Metrics" list):

1. ``relative_pose``: RMS displacement over the brain of ``A_v A_0^-1`` (relative to the
   first corrected volume) and of ``A_v Abar^-1`` (relative to the series' mean frame).
   A change of frame cancels out of both.
2. ``reference_consistency``: RMS displacement of ``A_v F^-1``, the offset between each
   corrected volume and what the reference depicts.
3. ``pair``: RMS displacement of ``A_c A_l^-1`` for each control-label pair.
4. ``deltam``: the candidate's mean delta-M against the oracle's, i.e. the same raw data
   resampled once with the true transforms. Two oracles: in the candidate's mean frame
   (``H*_v = P_v Abar^-1``, errors from relative misalignment only) and in the reference frame
   (``H*_v = P_v F^-1``, adding the common offset). Reported for the brain interior and for the
   outer shell (0-3 mm and 3-6 mm from the brain edge), as bias and RMSE, also relative to the
   oracle's grey-matter mean. The oracle's own error against the simulator truth is reported
   as the floor.
7. ``framewise``: estimated against true framewise displacement (RMS over the brain, between
   consecutive acquired volumes of the corrected series).

Usage::

    python score.py --fixture <fixture dir> --harness <harness output dir> [--out score.json]
"""

import argparse
import itertools
import json
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))  # the repository root, for aslprep.tests

from harness import load_itk_series, pair_deltam, resample_series  # noqa: E402

from aslprep.tests import truth_geometry as tg  # noqa: E402
from aslprep.tests.truth_scoring import Fixture  # noqa: E402

SCORE_VERSION = 2
#: Shell bounds (mm from the brain edge).
SHELLS = ((0.0, 3.0), (3.0, 6.0))


# ---------------------------------------------------------------------------------------------
# Transforms
# ---------------------------------------------------------------------------------------------
def polar_mean(mats):
    """Mean of rigid 4x4 matrices: averaged, with the linear part projected onto a rotation."""
    m = np.mean(np.asarray(mats), axis=0)
    u, _, vt = np.linalg.svd(m[:3, :3])
    r = u @ vt
    if np.linalg.det(r) < 0:
        u[:, -1] *= -1
        r = u @ vt
    out = np.eye(4)
    out[:3, :3] = r
    out[:3, 3] = m[:3, 3]
    return out


def frames(hmc, poses):
    """``A_v = H_v^-1 P_v`` for each corrected volume."""
    return np.stack([np.linalg.inv(h) @ p for h, p in zip(hmc, poses, strict=True)])


def _summary(values):
    values = np.asarray(values, dtype=float)
    if not values.size:
        return {'n': 0}
    return {
        'n': int(values.size),
        'median_mm': float(np.median(values)),
        'max_mm': float(np.max(values)),
        'by_volume_mm': [float(v) for v in values],
    }


def transform_metrics(a, context, points, ref_frame):
    """Metrics 1-3 from the frames ``a`` of the corrected volumes."""
    abar = polar_mean(a)
    out = {
        'relative_pose': {
            'to_first': _summary(
                [tg.rms_displacement(m @ np.linalg.inv(a[0]), points) for m in a]
            ),
            'to_mean': _summary([tg.rms_displacement(m @ np.linalg.inv(abar), points) for m in a]),
        },
        'reference_consistency': _summary(
            [tg.rms_displacement(m @ np.linalg.inv(ref_frame), points) for m in a]
        ),
        'mean_frame_offset_mm': tg.rms_displacement(abar @ np.linalg.inv(ref_frame), points),
    }
    controls = [i for i, t in enumerate(context) if t == 'control']
    labels = [i for i, t in enumerate(context) if t == 'label']
    if controls and len(controls) == len(labels):
        out['pair'] = _summary(
            [
                tg.rms_displacement(a[c] @ np.linalg.inv(a[lab]), points)
                for c, lab in zip(controls, labels, strict=True)
            ]
        )
    return out, abar


def framewise(hmc, poses, indices, points):
    """Estimated against true framewise displacement between consecutive corrected volumes."""
    order = np.argsort(indices)
    est, true = [], []
    for k0, k1 in itertools.pairwise(order):
        # The estimate of the head's motion between the two volumes is H_1 H_0^-1 (pull maps).
        est.append(tg.rms_displacement(hmc[k1] @ np.linalg.inv(hmc[k0]), points))
        true.append(tg.rms_displacement(poses[k1] @ np.linalg.inv(poses[k0]), points))
    est, true = np.asarray(est), np.asarray(true)
    if not est.size:
        return {'n': 0}
    corr = float(np.corrcoef(est, true)[0, 1]) if est.std() > 0 and true.std() > 0 else None
    return {
        'n': int(est.size),
        'mean_abs_error_mm': float(np.mean(np.abs(est - true))),
        'max_abs_error_mm': float(np.max(np.abs(est - true))),
        'estimated_mean_mm': float(est.mean()),
        'true_mean_mm': float(true.mean()),
        'correlation': corr,
    }


# ---------------------------------------------------------------------------------------------
# Delta-M
# ---------------------------------------------------------------------------------------------
def mask_in_frame(mask, mask_affine, frame, ref_img):
    """``mask`` (static frame) as seen on ``ref_img``'s grid in frame ``frame``.

    A corrected image in frame ``A`` shows ``static(A^-1 y)``, so the mask is pulled through
    ``A^-1``.
    """
    from scipy.ndimage import map_coordinates

    shape = ref_img.shape[:3]
    ijk = np.indices(shape).reshape(3, -1).T
    world = tg.apply_points(ref_img.affine, ijk)
    src = tg.apply_points(np.linalg.inv(mask_affine) @ np.linalg.inv(frame), world)
    vals = map_coordinates(mask.astype(np.float32), src.T, order=0, mode='constant', cval=0)
    return vals.reshape(shape) > 0.5


def shells(mask, affine):
    """Interior and edge shells of a mask, by distance (mm) from its edge."""
    from scipy.ndimage import distance_transform_edt

    sampling = np.sqrt((np.asarray(affine)[:3, :3] ** 2).sum(axis=0))
    # The EDT gives a voxel's distance to the nearest outside voxel centre, so edge voxels sit
    # at about one voxel; depth counts from there.
    depth = distance_transform_edt(mask, sampling=sampling) - sampling.min()
    out = {f'shell_{lo:g}_{hi:g}mm': mask & (depth >= lo) & (depth < hi) for lo, hi in SHELLS}
    out['interior'] = mask & (depth >= SHELLS[-1][1])
    return out


def _diff_stats(cand, oracle, regions, scale):
    out = {}
    for name, region in regions.items():
        d = (cand - oracle)[region]
        d = d[np.isfinite(d)]
        if not d.size:
            out[name] = {'n': 0}
            continue
        out[name] = {
            'n': int(d.size),
            'bias': float(np.mean(d)),
            'rmse': float(np.sqrt(np.mean(d**2))),
            'p95_abs': float(np.percentile(np.abs(d), 95)),
            'bias_rel_gm': float(np.mean(d) / scale) if scale else None,
            'rmse_rel_gm': float(np.sqrt(np.mean(d**2)) / scale) if scale else None,
        }
    return out


def deltam_metrics(fixture, harness_dir, hmc, poses, indices, context, abar, ref_frame):
    """Metric 4: candidate mean delta-M against the two oracles, by region."""
    cand_img = harness_dir / 'deltam.nii.gz'
    if not cand_img.exists():
        return {'available': False, 'reason': 'no control-label pairs'}
    ref_img = nb.load(harness_dir / 'aslref.nii.gz')
    # The harness may have run in a container; take the raw series from the fixture itself.
    asl_file = next((fixture.path / 'sub-01' / 'perf').glob('*_asl.nii.gz'))
    cand = np.asanyarray(nb.load(cand_img).dataobj).mean(-1)
    gm = fixture.pv['GM'] >= 0.7
    out = {'available': True}
    for name, target in (('oracle_mean_frame', abar), ('oracle_reference_frame', ref_frame)):
        pulls = [p @ np.linalg.inv(target) for p in poses]
        oracle = pair_deltam(
            resample_series(asl_file, indices, pulls, harness_dir / 'aslref.nii.gz'), context
        )
        oracle = np.asanyarray(oracle.dataobj)
        oracle = oracle[..., None] if oracle.ndim == 3 else oracle
        brain = mask_in_frame(fixture.brain, fixture.affine, target, ref_img)
        gm_t = mask_in_frame(gm, fixture.affine, target, ref_img)
        # The scale is the median over pairs of the oracle's GM mean: large spikes are not
        # rigid-correctable in the simulator (the B0 field stays with the scanner), and two
        # spiked pairs in ten can make the plain mean negative.
        scale = float(np.median(oracle[gm_t].mean(0))) if gm_t.any() else None
        oracle = oracle.mean(-1)
        regions = shells(brain, ref_img.affine)
        regions['brain'] = brain
        out[name] = _diff_stats(cand, oracle, regions, scale)
        out[name]['oracle_gm_mean'] = scale
        if name == 'oracle_reference_frame':
            out['oracle_floor'] = _oracle_floor(fixture, oracle, regions, scale, target, ref_img)
    return out


def _oracle_floor(fixture, oracle, regions, scale, target, ref_img):
    """The reference-frame oracle against the simulator's static delta-M truth."""
    gt_dir = fixture.path / 'sub-01' / 'perf' / 'ground-truth'
    gt_file = gt_dir / 'sub-01_desc-deltamStatic_gt.nii.gz'
    if not gt_file.exists():
        gt_file = gt_dir / 'sub-01_desc-deltam_gt.nii.gz'
    if not gt_file.exists():
        return {'available': False}
    from aslprep.tests import truth_models as tm

    gt = np.asanyarray(nb.load(gt_file).dataobj, dtype=np.float64)
    truth = tm.deltam_mean(gt, fixture.context) if gt.ndim == 4 else gt
    if not np.allclose(target, np.eye(4)) or ref_img.shape[:3] != fixture.shape:
        return {'available': False, 'reason': 'reference frame is not the static acquisition grid'}
    # The truth is in the simulator's units, before the readout and any suppression scaling:
    # fit one global scale over the brain, then report what that scale cannot explain.
    brain = regions['brain'] & np.isfinite(truth) & np.isfinite(oracle)
    t, o = truth[brain], oracle[brain]
    gain = float((o @ t) / (t @ t)) if t @ t > 0 else float('nan')
    return {
        'available': True,
        'gain': gain,
        **_diff_stats(oracle, gain * truth, regions, scale),
    }


# ---------------------------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------------------------
def reference_frame(fixture, harness, indices_all, a=None):
    """``F`` for the run's reference (see the module docstring).

    A reference averaged from corrected volumes (``reference_frame_kind`` ``corrected``, e.g.
    downstream policy D1) depicts those volumes' corrected frames, so ``F`` is the polar mean of
    their ``A_v``, not of their raw poses.
    """
    poses = fixture.poses()
    moving = fixture.motion is not None
    ref_type = harness.get('reference_type')
    ref_volumes = harness.get('reference_volumes')
    if harness.get('reference_frame_kind') == 'corrected' and a is not None and ref_volumes:
        pos = [indices_all.index(v) for v in ref_volumes]
        return polar_mean(a[pos]), f'corrected: polar mean of {len(pos)} corrected frame(s)'
    if not moving:
        return np.eye(4), 'static (no motion)'
    if ref_volumes == 'separate_m0':
        return np.eye(4), 'static (separate M0 does not move in aslscan)'
    if ref_volumes:
        note = 'exact: one volume' if len(ref_volumes) == 1 else 'approximate: polar mean of'
        return polar_mean(poses[ref_volumes]), f'{note} {len(ref_volumes)} volume pose(s)'
    if ref_type in ('separate_m0scan', 'separate m0scan'):
        return np.eye(4), 'static (separate M0 does not move in aslscan)'
    sel = [i for i, t in enumerate(fixture.context) if t == ref_type]
    if not sel:
        return np.eye(4), f'unknown ({ref_type}); identity assumed'
    return polar_mean(poses[sel]), f'approximate: polar mean of {len(sel)} {ref_type} poses'


def score(fixture_dir, harness_dir):
    fixture_dir, harness_dir = Path(fixture_dir), Path(harness_dir)
    manifest = json.loads((fixture_dir / 'manifest.json').read_text())
    # From the file name, so probe recipes and seed variants need no recipe lookup.
    asl = next((fixture_dir / 'sub-01' / 'perf').glob('*_asl.nii.gz'))
    acq = asl.name.split('_acq-')[1].split('_')[0]
    fixture = Fixture(fixture_dir, acq)
    harness = json.loads((harness_dir / 'harness.json').read_text())
    indices = json.loads((harness_dir / 'hmc_indices.json').read_text())
    context = pd.read_table(harness_dir / 'hmc_aslcontext.tsv')['volume_type'].tolist()
    hmc = load_itk_series(harness_dir / 'hmc_xfms.txt')
    poses = fixture.poses()[indices]
    _, points = fixture.voxel_centers(fixture.brain)
    a = frames(hmc, poses)
    ref_frame, ref_note = reference_frame(fixture, harness, indices, a)
    out = {
        'score_version': SCORE_VERSION,
        'fixture': str(fixture_dir),
        'recipe': manifest['recipe'],
        'harness': harness,
        'moving': fixture.motion is not None,
        'reference_frame': ref_note,
    }
    tm_out, abar = transform_metrics(a, context, points, ref_frame)
    out.update(tm_out)
    out['framewise'] = framewise(hmc, poses, indices, points)
    out['deltam'] = deltam_metrics(
        fixture, harness_dir, hmc, poses, indices, context, abar, ref_frame
    )
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--fixture', type=Path, required=True)
    parser.add_argument('--harness', type=Path, required=True)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args(argv)
    result = score(args.fixture, args.harness)
    text = json.dumps(result, indent=2)
    if args.out:
        args.out.write_text(text + '\n')
    else:
        print(text)


if __name__ == '__main__':
    main()
