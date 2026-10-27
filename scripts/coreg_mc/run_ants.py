#!/usr/bin/env python3
"""ANTs registration arms of the coregistration probe, scored against aslscan truth.

Plan: "Candidates" (S1, S2, Q1, J1) and Task 2's engine pilot. For each fixture, the HMC
probe's ``coreg_experiment.run_fixture`` (in ``scripts/hmc_probe``) prepares the work directory
from an ASLPrep run on that fixture (T1w brain and mask, white matter, the conventional reference
and its ASL brain mask, the perfusion images) and runs production's FSL path (``mri_coreg``, then
FLIRT-BBR with both signs). A grey-matter target is added here, and every ANTs arm starts from
the same ``mri_coreg`` estimate (``conventional_init.mat``), so arms differ only in their
objective.

Transform conventions. ASLPrep's ``from-aslref_to-T1w`` matrix, read as RAS by nitransforms, is
a pull map: T1w point -> aslref point. ANTs with the anatomy as the fixed image estimates the
same pull map (fixed point -> moving point), so its initial transform is that RAS pull written
in ITK form, and its output is read back the same way. Truth is ``F R^-1`` as in the HMC probe's
``coreg_experiment.run_fixture``.

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/run_ants.py --fixtures <aslscan fixtures> --runs <fixture runs> \\
        --license <license.txt> --work <dir> --out scores.tsv --names rc_pcasl rc_grase ...
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
from nibabel.affines import apply_affine
from scipy.ndimage import gaussian_filter

HMC_PROBE = Path(__file__).resolve().parents[1] / 'hmc_probe'
sys.path.insert(0, str(HMC_PROBE))

import coreg_experiment  # noqa: E402
from coreg_experiment import _pull, estimates  # noqa: E402
from run_panel import IMAGE  # noqa: E402 - the pinned image
from score import polar_mean  # noqa: E402
from score_full import _load_series, _one, _transform_error  # noqa: E402

from aslprep.tests import truth_scoring as ts  # noqa: E402

SEED = 20261022

#: Stage schedule, frozen by the engine pilot (Task 2) before Task 3. The fixed images are on
#: the ASL-resolution grid (``lowres``), so two levels suffice.
SCHEDULE = '--convergence [500x250,1e-6,10] --shrink-factors 2x1 --smoothing-sigmas 1x0vox'
#: A registration that takes longer than this has stalled (the first pilot's dense CC on the
#: 1 mm T1w grid ran for 19 hours without finishing one arm).
TIMEOUT_S = 1800


def mi(fixed, moving, weight):
    return f'MI[{fixed},{moving},{weight},32,Regular,0.25]'


def cc(fixed, moving, weight):
    return f'CC[{fixed},{moving},{weight},2]'


#: arm -> stages; each stage is a list of metric terms (fixed = anatomy, moving = ASL).
ARMS = {
    'S1': [[mi('/w/t1w_brain_lr.nii.gz', '/w/conventional.nii.gz', 1)]],
    'S2': [[cc('/w/gm_target_lr.nii.gz', '/w/P-winsor.nii.gz', 1)]],
    'Q1': [
        [mi('/w/t1w_brain_lr.nii.gz', '/w/conventional.nii.gz', 1)],
        [cc('/w/gm_target_lr.nii.gz', '/w/P-winsor.nii.gz', 1)],
    ],
    'J1': [
        [
            mi('/w/t1w_brain_lr.nii.gz', '/w/conventional.nii.gz', 0.6),
            cc('/w/gm_target_lr.nii.gz', '/w/P-winsor.nii.gz', 0.4),
        ],
    ],
}
#: Engine-pilot variants of the static term (Task 2), and a repeat of S2 for determinism.
PILOT_ARMS = {
    'S1-mattes': [['Mattes[/w/t1w_brain_lr.nii.gz,/w/conventional.nii.gz,1,32,Regular,0.25]']],
    'S1-dense': [['MI[/w/t1w_brain_lr.nii.gz,/w/conventional.nii.gz,1,32,None]']],
    'S1-mattes-dense': [['Mattes[/w/t1w_brain_lr.nii.gz,/w/conventional.nii.gz,1,32,None]']],
    'S1-mattes64-dense': [['Mattes[/w/t1w_brain_lr.nii.gz,/w/conventional.nii.gz,1,64,None]']],
    'S1-fullres': [[mi('/w/t1w_brain.nii.gz', '/w/conventional.nii.gz', 1)]],
    'S2-repeat': ARMS['S2'],
}
ALL_ARMS = ARMS | PILOT_ARMS


def gm_target(work, gm_probseg):
    """The GM probability map smoothed to the ASL resolution, on the T1w grid."""
    gm = nb.load(gm_probseg)
    t1 = nb.load(work / 't1w_brain.nii.gz')
    asl = nb.load(work / 'conventional.nii.gz')
    if gm.shape != t1.shape or not np.allclose(gm.affine, t1.affine, atol=1e-3):
        raise ValueError(f'{gm_probseg} is not on the T1w grid')
    t1_vox = np.sqrt((t1.affine[:3, :3] ** 2).sum(0))
    # The ASL voxel size along the T1w axes (rotation-aware, a conservative per-axis estimate).
    asl_vox = np.sqrt((asl.affine[:3, :3] ** 2).sum(0)).max()
    fwhm = np.sqrt(np.maximum(asl_vox**2 - t1_vox**2, 0))
    sigma_vox = fwhm / 2.3548 / t1_vox
    data = gaussian_filter(np.asanyarray(gm.dataobj, np.float32), sigma_vox)
    nb.Nifti1Image(data, t1.affine).to_filename(work / 'gm_target.nii.gz')


def lowres(work):
    """The anatomical fixed images on an isotropic grid at the ASL resolution.

    ANTs evaluates its metrics on the fixed image's grid; the 1 mm T1w grid gives no
    information the ASL image can use and makes dense CC far too slow. World coordinates are
    unchanged, so transforms are unaffected.
    """
    from nibabel.processing import resample_to_output

    asl = nb.load(work / 'conventional.nii.gz')
    vox = float(np.sqrt((asl.affine[:3, :3] ** 2).sum(0)).max())
    for name, order in (('t1w_brain', 1), ('t1w_mask', 0), ('gm_target', 1)):
        img = nb.load(work / f'{name}.nii.gz')
        img = nb.Nifti1Image(np.asanyarray(img.dataobj, np.float32), img.affine)
        out = resample_to_output(img, voxel_sizes=(vox, vox, vox), order=order)
        out.to_filename(work / f'{name}_lr.nii.gz')


def write_init(work):
    """Production's mri_coreg estimate as an ITK pull transform for ANTs."""
    import nitransforms as nt

    pull = _pull(
        work / 'conventional_init.mat',
        nb.load(work / 'conventional.nii.gz'),
        nb.load(work / 't1w_brain.nii.gz'),
    )
    nt.linear.Affine(pull).to_filename(work / 'init_itk.txt', fmt='itk')
    return pull


def register(work, arms):
    lines = []
    for arm in arms:
        cmd = [
            f'timeout {TIMEOUT_S} antsRegistration --verbose 1 --dimensionality 3 --float 0',
            f'--collapse-output-transforms 1 --output [/w/ants_{arm}_] --interpolation Linear',
            '--winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0',
            f'--random-seed {SEED} --initial-moving-transform /w/init_itk.txt',
        ]
        for terms in ALL_ARMS[arm]:
            cmd.append('--transform Rigid[0.1]')
            cmd += [f'--metric {t}' for t in terms]
            cmd.append(SCHEDULE)
            cmd.append('--masks [/w/t1w_mask_lr.nii.gz,/w/mask.nii.gz]')
        lines.append(' '.join(cmd) + f' > /w/ants_{arm}.log 2>&1')
        lines.append(f'echo "{arm} exit $?" >> /w/ants_status.txt')
    # No ``set -e``: an arm that fails or times out has no transform and is reported as such.
    script = ['export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1', ': > /w/ants_status.txt', *lines]
    (work / 'ants.sh').write_text('\n'.join(script) + '\n')
    cmd = [
        'docker',
        'run',
        '--rm',
        '-u',
        f'{os.getuid()}:{os.getgid()}',
        '-e',
        'HOME=/tmp',
        '-v',
        f'{work}:/w',
        '--entrypoint',
        'bash',
        IMAGE,
        '/w/ants.sh',
    ]
    subprocess.run(cmd, check=True, capture_output=True)


def read_ants(work, arm):
    import nitransforms as nt

    f = work / f'ants_{arm}_0GenericAffine.mat'
    return np.asarray(nt.linear.load(f, fmt='itk').matrix) if f.exists() else None


def fixture_truth(fixture_dir, run_dir):
    """Truth pull and scoring points, exactly as ``coreg_experiment.run_fixture`` defines them."""
    perf = run_dir / 'sub-01' / 'perf'
    fperf = fixture_dir / 'sub-01' / 'perf'
    asl = next(fperf.glob('*_asl.nii.gz'))
    acq = asl.name.split('_acq-')[1].split('_')[0]
    fixture = ts.Fixture(fixture_dir, acq)
    _, static_points = fixture.voxel_centers(fixture.brain)
    points = apply_affine(fixture.R, static_points)
    hmc = _load_series(next(perf.glob(f'*_acq-{acq}_from-orig_to-aslref*_xfm.txt')))
    poses = fixture.poses()
    keep = [i for i, t in enumerate(fixture.context) if t in ('control', 'label', 'm0scan')]
    keep = keep if len(keep) == len(hmc) else list(range(len(hmc)))
    frames = [np.linalg.inv(h) @ poses[i] for h, i in zip(hmc, keep, strict=True)]
    series_frame = polar_mean(frames)
    truth = series_frame @ np.linalg.inv(fixture.R)
    production = _one(ts.find_xfm(run_dir, acq, 'from-aslref_to-T1w'))
    voxel = float(np.sqrt((fixture.affine[:3, :3] ** 2).sum(axis=0)).max())
    return truth, points, production, voxel


def run_fixture(name, args):
    work = args.work / name
    if args.ants_only:
        # Reuse the prepared work directory (and its FSL results); register the ANTs arms only.
        register(work, args.arms)
    elif not args.rescore:
        shutil.rmtree(work, ignore_errors=True)
        # Prepares the work directory and runs production's FSL path (both BBR signs).
        coreg_experiment.run_fixture(args.fixtures / name, args.runs / name, work, args.license)
        anat = args.fixtures / name / 'derivatives' / 'anat' / 'sub-01' / 'anat'
        gm_target(work, next(anat.glob('*_label-GM_probseg.nii.gz')))
        lowres(work)
        write_init(work)
        register(work, args.arms)
    truth, points, production, voxel = fixture_truth(args.fixtures / name, args.runs / name)
    init = _pull(
        work / 'conventional_init.mat',
        nb.load(work / 'conventional.nii.gz'),
        nb.load(work / 't1w_brain.nii.gz'),
    )
    found = {'production': production, 'mri_coreg init': init}
    for (source, sign), (final, _bbr, _init, fallback) in estimates(work).items():
        if source in ('conventional', 'P-winsor'):
            found[f'FSL {source} BBR {sign}' + (' (fallback)' if fallback else '')] = final
    for arm in args.arms:
        found[f'ANTs {arm}'] = read_ants(work, arm)
    rows = []
    for arm, est in found.items():
        row = {'fixture': name, 'arm': arm, 'voxel_mm': voxel}
        if est is None:
            row['problem'] = 'no transform'
        else:
            err = _transform_error(est, truth, points)
            row |= {
                'rms_mm': err['rms_mm'],
                'rot_deg': err['rot_deg'],
                'rms_voxels': err['rms_mm'] / voxel,
                'inverted_rms_mm': err['rms_mm_if_inverted'],
            }
        rows.append(row)
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--fixtures', type=Path, required=True)
    p.add_argument('--runs', type=Path, required=True)
    p.add_argument('--license', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--names', nargs='+', required=True)
    p.add_argument('--arms', nargs='+', default=list(ARMS), choices=list(ALL_ARMS))
    p.add_argument('--rescore', action='store_true')
    p.add_argument('--ants-only', action='store_true')
    args = p.parse_args(argv)
    rows = []
    for name in args.names:
        rows += run_fixture(name, args)
        print(f'{name}: done', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    table = df.pivot_table(index='arm', columns='fixture', values='rms_mm', aggfunc='first')
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(table.round(2))
    print(json.dumps({'schedule': SCHEDULE, 'seed': SEED}))


if __name__ == '__main__':
    main()
