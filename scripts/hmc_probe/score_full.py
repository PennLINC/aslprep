#!/usr/bin/env python3
"""Score the reference consumers of a full ASLPrep run on a simulated fixture (metrics 9-11).

Frames as in ``score.py`` and ``aslprep/tests/truth_scoring.py`` (RAS point mappings of the
static phantom). The aslref's frame is ``F = polar mean of A_v = H_v^-1 P_v`` over the
corrected volumes (the identity when HMC is perfect and the reference is static); every map
below is scored through it, in the scanner frame, so a compensating transform is not counted
as an error.

- 9 ``mask``: ASLPrep's ASL brain mask against the true brain seen in the aslref's frame:
  Dice, and missed and extra volume (ml), split into the outer 6 mm and the rest.
- 10 ``coregistration``: ``truth_scoring.score_frames`` (end to end and the registration step
  alone), plus the registration step measured against the aslref's actual frame ``F``:
  ``C F^-1 R^-1`` should be the identity (``C``: aslref -> T1w, pull from T1w).
- 11 ``m0`` and ``fieldmap``: the M0-to-aslref and aslref-to-fieldmap transforms against the
  known poses in ``ground-truth/probe_extras.json``.

Usage::

    python score_full.py --fixture <fixture dir> --aslprep <output dir> [--out score.json]
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

from aslprep.tests import truth_geometry as tg  # noqa: E402
from aslprep.tests import truth_scoring as ts  # noqa: E402

SCORE_FULL_VERSION = 1


def _load_series(path):
    from nitransforms.linear import load

    return np.asarray(load(path, fmt='itk').matrix).reshape(-1, 4, 4)


def _one(path):
    return _load_series(path)[0]


def _xfm(perf, pattern):
    hits = sorted(perf.glob(pattern))
    return hits[0] if len(hits) == 1 else None


def _transform_error(est_pull, truth_pull, points):
    """RMS displacement and rotation of ``est^-1 truth`` (identity when est == truth)."""
    d = np.linalg.inv(est_pull) @ truth_pull
    return {
        'rms_mm': tg.rms_displacement(d, points),
        'rot_deg': tg.rot_angle_deg(d),
        'rms_mm_if_inverted': tg.rms_displacement(est_pull @ truth_pull, points),
    }


def score_full(fixture_dir, aslprep_dir):
    fixture_dir, aslprep_dir = Path(fixture_dir), Path(aslprep_dir)
    asl = next((fixture_dir / 'sub-01' / 'perf').glob('*_asl.nii.gz'))
    acq = asl.name.split('_acq-')[1].split('_')[0]
    fixture = ts.Fixture(fixture_dir, acq)
    perf = aslprep_dir / 'sub-01' / 'perf'
    out = {'score_full_version': SCORE_FULL_VERSION, 'fixture': str(fixture_dir)}
    _, points = fixture.voxel_centers(fixture.brain)

    hmc_file = _xfm(perf, f'sub-01_acq-{acq}_from-orig_to-aslref*_xfm.txt')
    if hmc_file is None:
        return {**out, 'problems': ['no HMC transform']}
    hmc = _load_series(hmc_file)
    poses = fixture.poses()
    # HMC transforms cover the processed volumes only; ASLPrep keeps control/label/m0scan
    # (or deltam/m0scan) in acquisition order.
    keep = [
        i for i, t in enumerate(fixture.context) if t in ('control', 'label', 'm0scan', 'deltam')
    ]
    if len(keep) != len(hmc):
        keep = list(range(len(hmc)))
    frame = polar_mean([np.linalg.inv(h) @ poses[i] for h, i in zip(hmc, keep, strict=True)])
    out['aslref_frame_offset_mm'] = tg.rms_displacement(frame, points)

    # 9: brain mask
    mask_file = ts.find_output(aslprep_dir, acq, 'mask', desc='brain')
    if mask_file is not None:
        m_img = nb.load(mask_file)
        est = np.asanyarray(m_img.dataobj) > 0.5
        true = mask_in_frame(fixture.brain, fixture.affine, frame, m_img)
        vox_ml = abs(np.linalg.det(m_img.affine[:3, :3])) / 1000
        regions = shells(true, m_img.affine)
        edge = regions['shell_0_3mm'] | regions['shell_3_6mm']
        out['mask'] = {
            'dice': float(2 * (est & true).sum() / (est.sum() + true.sum())),
            'missed_ml': float((true & ~est).sum() * vox_ml),
            'missed_edge_ml': float((edge & ~est).sum() * vox_ml),
            'missed_interior_ml': float((regions['interior'] & ~est).sum() * vox_ml),
            'extra_ml': float((est & ~true).sum() * vox_ml),
            'true_ml': float(true.sum() * vox_ml),
        }

    # 10: coregistration
    frames = ts.score_frames(fixture, aslprep_dir)
    out['frames'] = frames
    coreg_file = ts.find_xfm(aslprep_dir, acq, 'from-aslref_to-T1w')
    if coreg_file is not None:
        c = _one(coreg_file)
        # aslref point y shows static(F^-1 y); T1w point x shows static(R^-1 x). The pull from
        # T1w to aslref is then C = F R^-1, so C^-1 F R^-1 is the registration step's error.
        out['coregistration_vs_aslref_frame'] = _transform_error(
            c, frame @ np.linalg.inv(fixture.R), points
        )

    # 11: M0 and fieldmap registration against known poses
    extras_file = fixture_dir / 'sub-01' / 'perf' / 'ground-truth' / 'probe_extras.json'
    extras = json.loads(extras_file.read_text()) if extras_file.exists() else {}
    m0_file = _xfm(perf, f'sub-01_acq-{acq}_from-m0scan_to-aslref*_xfm.txt')
    if m0_file is not None:
        m0_pose = np.asarray(extras.get('m0', {}).get('pose', np.eye(4)))
        # The M0 image shows static(Pm^-1 x); the pull from aslref to M0 is Pm F^-1.
        out['m0'] = _transform_error(_one(m0_file), m0_pose @ np.linalg.inv(frame), points)
    if 'fieldmap' in extras:
        fmap_files = [
            p
            for p in perf.glob(f'sub-01_acq-{acq}_from-aslref_to-*_xfm.txt')
            if 'T1w' not in p.name
        ]
        fmap_file = fmap_files[0] if len(fmap_files) == 1 else None
        if fmap_file is None:
            out['fieldmap'] = {'problems': ['no aslref-to-fieldmap transform found']}
        else:
            q = np.asarray(extras['fieldmap']['pose'])
            out['fieldmap'] = {
                'file': fmap_file.name,
                **_transform_error(_one(fmap_file), q @ np.linalg.inv(frame), points),
            }
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--fixture', type=Path, required=True)
    parser.add_argument('--aslprep', type=Path, required=True)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args(argv)
    text = json.dumps(score_full(args.fixture, args.aslprep), indent=2, default=float)
    if args.out:
        args.out.write_text(text + '\n')
    else:
        print(text)


if __name__ == '__main__':
    main()
