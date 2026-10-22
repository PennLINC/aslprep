#!/usr/bin/env python3
"""ANTs registration arms of the coregistration probe, scored against aslscan truth.

Plan: "Candidates" (S1, S2, Q1, J1) and Task 2's engine pilot. Every arm works on a work
directory prepared by the HMC probe's ``coreg_experiment.prepare`` (T1w brain and mask, white
matter, the conventional reference and its ASL brain mask, and the perfusion images), plus a
grey-matter target written here. Every arm starts from the same ``mri_coreg`` estimate as
production (``conventional_init.mat`` from the same work directory), so arms differ only in
their objective.

Transform conventions. ASLPrep's ``from-aslref_to-T1w`` matrix, read as RAS by nitransforms, is
a pull map: T1w point -> aslref point. ANTs with the anatomy as the fixed image estimates the
same pull map (fixed point -> moving point), so its initial transform is that RAS pull written
in ITK form, and its output is read back the same way. Truth is ``F R^-1`` as in the HMC probe's
``coreg_experiment.run_fixture``.

Usage (in WSL, ``hmc-probe`` environment)::

    python run_ants.py --fixtures <aslscan fixtures> --runs <fixture-runs/baseline> \\
        --prepared <coreg-5b> --work <dir> --out scores.tsv --names rc_pcasl rc_grase ...
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

HMC_PROBE = Path('/mnt/c/Users/tsalo/Documents/linc/aslprep-hmc-probe/scripts/hmc_probe')
sys.path.insert(0, str(HMC_PROBE))

from coreg_experiment import _pull, estimates  # noqa: E402
from score import polar_mean  # noqa: E402
from score_full import _load_series, _one, _transform_error  # noqa: E402

from aslprep.tests import truth_scoring as ts  # noqa: E402

IMAGE = 'pennlinc/aslprep@sha256:eacb03704072b411ca47ee89cbf17ea819b44640650ee3386beca8d29d9f78ff'
SEED = 20261022

#: Stage schedule, frozen by the engine pilot (Task 2) before Task 3.
SCHEDULE = '--convergence [1000x500x250,1e-7,10] --shrink-factors 4x2x1 --smoothing-sigmas 2x1x0vox'


def mi(fixed, moving, weight):
    return f'MI[{fixed},{moving},{weight},32,Regular,0.25]'


def cc(fixed, moving, weight):
    return f'CC[{fixed},{moving},{weight},4]'


#: arm -> stages; each stage is a list of metric terms (fixed = anatomy, moving = ASL).
ARMS = {
    'S1': [[mi('/w/t1w_brain.nii.gz', '/w/conventional.nii.gz', 1)]],
    'S2': [[cc('/w/gm_target.nii.gz', '/w/P-winsor.nii.gz', 1)]],
    'Q1': [
        [mi('/w/t1w_brain.nii.gz', '/w/conventional.nii.gz', 1)],
        [cc('/w/gm_target.nii.gz', '/w/P-winsor.nii.gz', 1)],
    ],
    'J1': [[
        mi('/w/t1w_brain.nii.gz', '/w/conventional.nii.gz', 0.6),
        cc('/w/gm_target.nii.gz', '/w/P-winsor.nii.gz', 0.4),
    ]],
}


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
            'antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1',
            f'--output [/w/ants_{arm}_] --interpolation Linear',
            '--winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0',
            f'--random-seed {SEED} --initial-moving-transform /w/init_itk.txt',
        ]
        for terms in ARMS[arm]:
            cmd.append('--transform Rigid[0.1]')
            cmd += [f'--metric {t}' for t in terms]
            cmd.append(SCHEDULE)
            cmd.append('--masks [/w/t1w_mask.nii.gz,/w/mask.nii.gz]')
        lines.append(' '.join(cmd) + f' > /w/ants_{arm}.log 2>&1')
    (work / 'ants.sh').write_text('export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1\n' + '\n'.join(lines) + '\n')
    cmd = ['docker', 'run', '--rm', '-u', f'{os.getuid()}:{os.getgid()}', '-e', 'HOME=/tmp',
           '-v', f'{work}:/w', '--entrypoint', 'bash', IMAGE, '/w/ants.sh']
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
    series_frame = polar_mean([np.linalg.inv(h) @ poses[i] for h, i in zip(hmc, keep, strict=True)])
    truth = series_frame @ np.linalg.inv(fixture.R)
    production = _one(ts.find_xfm(run_dir, acq, 'from-aslref_to-T1w'))
    voxel = float(np.sqrt((fixture.affine[:3, :3] ** 2).sum(axis=0)).max())
    return truth, points, production, voxel


def run_fixture(name, args):
    src = args.prepared / name
    work = args.work / name
    if not args.rescore:
        shutil.rmtree(work, ignore_errors=True)
        shutil.copytree(src, work)
        gm = next((args.fixtures / name / 'derivatives' / 'anat' / 'sub-01' / 'anat').glob(
            '*_label-GM_probseg.nii.gz'))
        gm_target(work, gm)
        write_init(work)
        register(work, args.arms)
    truth, points, production, voxel = fixture_truth(args.fixtures / name, args.runs / name)
    init = _pull(work / 'conventional_init.mat', nb.load(work / 'conventional.nii.gz'),
                 nb.load(work / 't1w_brain.nii.gz'))
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
            row |= {'rms_mm': err['rms_mm'], 'rot_deg': err['rot_deg'],
                    'rms_voxels': err['rms_mm'] / voxel,
                    'inverted_rms_mm': err['rms_mm_if_inverted']}
        rows.append(row)
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--fixtures', type=Path, required=True)
    p.add_argument('--runs', type=Path, required=True)
    p.add_argument('--prepared', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--names', nargs='+', required=True)
    p.add_argument('--arms', nargs='+', default=list(ARMS))
    p.add_argument('--rescore', action='store_true')
    args = p.parse_args(argv)
    rows = []
    for name in args.names:
        rows += run_fixture(name, args)
        print(f'{name}: done', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(df.pivot_table(index='arm', columns='fixture', values='rms_mm', aggfunc='first').round(2))
    print(json.dumps({'schedule': SCHEDULE, 'seed': SEED}))


if __name__ == '__main__':
    main()
