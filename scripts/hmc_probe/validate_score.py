#!/usr/bin/env python3
"""Validate the scorer (Task 2) with transforms whose correct scores are known.

For a fixture, write harness-style outputs for synthetic candidates and score them:

- ``identity``: no correction (``H_v = I``);
- ``oracle``: the true transforms (``H_v = P_v``);
- ``shifted``: the true transforms in a frame shifted by a known rigid ``S``
  (``H_v = P_v S^-1``), so every corrected volume sits at ``S``;
- ``one_wrong``: the true transforms with one volume off by a known rigid ``E``.

Expected: the oracle scores zero everywhere; ``shifted`` scores zero on relative pose, pairs
and the mean-frame delta-M, and ``|S|`` on reference consistency; ``one_wrong`` scores
``|E|`` on that volume's pose relative to the others; ``identity`` scores the true motion.
A failed expectation exits non-zero.

Usage::

    python validate_score.py --fixture <fixture dir> --out <dir>
"""

import argparse
import json
import shutil
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

from harness import pair_deltam, resample_series  # noqa: E402
from score import score  # noqa: E402

from aslprep.tests import truth_geometry as tg  # noqa: E402
from aslprep.tests.aslscan_fixtures import load_recipe  # noqa: E402
from aslprep.tests.truth_scoring import Fixture  # noqa: E402

S = tg.rigid_ras(0.7, -0.4, 0.3, 0.8, -0.5, 0.6)
E = tg.rigid_ras(0.5, 0.0, 0.0, 0.0, 0.0, 1.0)
#: ITK files carry single-precision round-off, about 1e-6 mm here.
TOL_MM = 1e-4


def write_itk_series(path, mats):
    """Write RAS pull matrices as an ITK multi-transform file (MCFLIRT2ITK's layout)."""
    lines = ['#Insight Transform File V1.0']
    for k, m in enumerate(mats):
        lps = tg.ras_to_lps(np.asarray(m, dtype=float))
        params = [*lps[:3, :3].ravel(), *lps[:3, 3]]
        lines += [
            f'#Transform {k}',
            'Transform: AffineTransform_double_3_3',
            f'Parameters: {" ".join(repr(float(v)) for v in params)}',
            'FixedParameters: 0 0 0',
        ]
    Path(path).write_text('\n'.join(lines) + '\n')


def make_candidate(fixture, asl_file, out, pulls):
    out.mkdir(parents=True, exist_ok=True)
    context = fixture.context
    indices = [i for i, t in enumerate(context) if t in ('control', 'label', 'm0scan')]
    ref = fixture.asl_img.slicer[..., 0]
    nb.Nifti1Image(np.asanyarray(ref.dataobj, dtype=np.float32), ref.affine).to_filename(
        out / 'aslref.nii.gz'
    )
    sel = [pulls[i] for i in indices]
    write_itk_series(out / 'hmc_xfms.txt', sel)
    ctx = [context[i] for i in indices]
    pd.DataFrame({'volume_type': ctx}).to_csv(out / 'hmc_aslcontext.tsv', sep='\t', index=False)
    (out / 'hmc_indices.json').write_text(json.dumps(indices))
    corrected = resample_series(asl_file, indices, sel, out / 'aslref.nii.gz')
    corrected.to_filename(out / 'corrected.nii.gz')
    dm = pair_deltam(corrected, ctx)
    if dm is not None:
        dm.to_filename(out / 'deltam.nii.gz')
    (out / 'harness.json').write_text(
        json.dumps(
            {
                'candidate': out.name,
                'inputs': {'asl_file': str(asl_file)},
                'reference_type': 'static-synthetic',
            }
        )
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--fixture', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    manifest = json.loads((args.fixture / 'manifest.json').read_text())
    fixture = Fixture(args.fixture, load_recipe(manifest['recipe']).acq)
    asl_file = next((args.fixture / 'sub-01' / 'perf').glob('*_asl.nii.gz'))
    poses = fixture.poses()
    _, points = fixture.voxel_centers(fixture.brain)
    n = len(poses)
    wrong = n // 2 if fixture.context[n // 2] != 'm0scan' else n // 2 + 1
    cases = {
        'identity': [np.eye(4)] * n,
        'oracle': list(poses),
        'shifted': [p @ np.linalg.inv(S) for p in poses],
        'one_wrong': [p @ E if i == wrong else p for i, p in enumerate(poses)],
    }
    shutil.rmtree(args.out, ignore_errors=True)
    results, failures = {}, []
    for name, pulls in cases.items():
        make_candidate(fixture, asl_file, args.out / name, pulls)
        res = score(args.fixture, args.out / name)
        # The synthetic reference is volume 0, so its frame is P_0, not the static frame.
        results[name] = res
    s_mm = tg.rms_displacement(S, points)
    e_mm = tg.rms_displacement(E, points)

    def check(cond, msg):
        if not cond:
            failures.append(msg)

    o = results['oracle']
    check(o['relative_pose']['to_mean']['max_mm'] < TOL_MM, 'oracle relative pose not zero')
    check(o['reference_consistency']['max_mm'] < TOL_MM, 'oracle reference consistency not zero')
    if 'pair' in o:
        check(o['pair']['max_mm'] < TOL_MM, 'oracle pair not zero')
    if o['deltam']['available']:
        for frame in ('oracle_mean_frame', 'oracle_reference_frame'):
            check(o['deltam'][frame]['brain']['rmse'] < 1e-4, f'oracle delta-M differs ({frame})')
    sh = results['shifted']
    check(sh['relative_pose']['to_mean']['max_mm'] < TOL_MM, 'shift changed relative pose')
    if 'pair' in sh:
        check(sh['pair']['max_mm'] < TOL_MM, 'shift changed pair error')
    check(
        abs(sh['reference_consistency']['median_mm'] - s_mm) < 1e-3,
        f'shift consistency {sh["reference_consistency"]["median_mm"]} != {s_mm}',
    )
    if sh['deltam']['available']:
        check(
            sh['deltam']['oracle_mean_frame']['brain']['rmse'] < 1e-3,
            'shift changed mean-frame delta-M',
        )
    ow = results['one_wrong']
    k = json.loads((args.out / 'one_wrong' / 'hmc_indices.json').read_text()).index(wrong)
    got = ow['relative_pose']['to_first']['by_volume_mm'][k]
    check(abs(got - e_mm) < 1e-3, f'one_wrong volume error {got} != {e_mm}')
    summary = {
        name: {
            'relative_to_mean_max': r['relative_pose']['to_mean']['max_mm'],
            'reference_consistency_median': r['reference_consistency']['median_mm'],
            'pair_max': r.get('pair', {}).get('max_mm'),
            'deltam_mean_frame_rmse': r['deltam']
            .get('oracle_mean_frame', {})
            .get('brain', {})
            .get('rmse'),
        }
        for name, r in results.items()
    }
    summary['expected'] = {
        'shift_mm': s_mm,
        'one_wrong_mm': e_mm,
        'moving': fixture.motion is not None,
    }
    print(json.dumps(summary, indent=2))
    (args.out / 'validation.json').write_text(
        json.dumps({'summary': summary, 'failures': failures}, indent=2)
    )
    if failures:
        print('FAILED:', *failures, sep='\n  ')
        sys.exit(1)
    print('all expectations met')


if __name__ == '__main__':
    main()
