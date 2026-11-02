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
    # The fixed-image direction: ASL fixed, anatomy moving (metrics sampled on the ASL grid).
    'S1-inv': [[mi('/w/conventional.nii.gz', '/w/t1w_brain.nii.gz', 1)]],
    'S2-inv': [[cc('/w/P-winsor.nii.gz', '/w/gm_target.nii.gz', 1)]],
    'J1-inv': [
        [
            mi('/w/conventional.nii.gz', '/w/t1w_brain.nii.gz', 0.6),
            cc('/w/P-winsor.nii.gz', '/w/gm_target.nii.gz', 0.4),
        ],
    ],
}
ALL_ARMS = ARMS | PILOT_ARMS
#: Arms with the ASL image fixed: ANTs estimates the push map (aslref point -> T1w point), so
#: its initial transform is the inverse of production's and its result is inverted on reading.
INVERSE = {'S1-inv', 'S2-inv', 'J1-inv'}


#: ExploreASL's pseudo-tissue weights (``xASL_im_CreatePseudoCBF.m``: pGM * 80 + pWM * 26.7),
#: normalized to GM; CC does not depend on scale.
PSEUDO_CBF_WM = 26.7 / 80


def _smoothed_target(work, probsegs, weights, name):
    """A weighted sum of T1w-grid probability maps smoothed to the ASL resolution."""
    t1 = nb.load(work / 't1w_brain.nii.gz')
    asl = nb.load(work / 'conventional.nii.gz')
    data = np.zeros(t1.shape, np.float32)
    for probseg, weight in zip(probsegs, weights, strict=True):
        img = nb.load(probseg)
        if img.shape != t1.shape or not np.allclose(img.affine, t1.affine, atol=1e-3):
            raise ValueError(f'{probseg} is not on the T1w grid')

        data += weight * np.asanyarray(img.dataobj, np.float32)

    t1_vox = np.sqrt((t1.affine[:3, :3] ** 2).sum(0))
    # The ASL voxel size along the T1w axes (rotation-aware, a conservative per-axis estimate).
    asl_vox = np.sqrt((asl.affine[:3, :3] ** 2).sum(0)).max()
    fwhm = np.sqrt(np.maximum(asl_vox**2 - t1_vox**2, 0))
    sigma_vox = fwhm / 2.3548 / t1_vox
    nb.Nifti1Image(gaussian_filter(data, sigma_vox), t1.affine).to_filename(work / name)


def gm_target(work, gm_probseg):
    """The GM probability map smoothed to the ASL resolution, on the T1w grid."""
    _smoothed_target(work, [gm_probseg], [1.0], 'gm_target.nii.gz')


def pseudo_cbf_target(work, gm_probseg, wm_probseg):
    """ExploreASL's pseudo-tissue image (GM + 0.334 WM), smoothed as the GM target is."""
    _smoothed_target(
        work,
        [gm_probseg, wm_probseg],
        [1.0, PSEUDO_CBF_WM],
        'pseudocbf_target.nii.gz',
    )


def lowres(work, names=(('t1w_brain', 1), ('t1w_mask', 0), ('gm_target', 1))):
    """The anatomical fixed images on an isotropic grid at the ASL resolution.

    ANTs evaluates its metrics on the fixed image's grid; the 1 mm T1w grid gives no
    information the ASL image can use and makes dense CC far too slow. World coordinates are
    unchanged, so transforms are unaffected. ``names`` are (image, interpolation order) pairs.
    """
    from nibabel.processing import resample_to_output

    asl = nb.load(work / 'conventional.nii.gz')
    vox = float(np.sqrt((asl.affine[:3, :3] ** 2).sum(0)).max())
    for name, order in names:
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
    nt.linear.Affine(np.linalg.inv(pull)).to_filename(work / 'init_itk_inv.txt', fmt='itk')
    return pull


def ants_command(arm, init, prefix, schedule=SCHEDULE, stages=None, masked=True):
    """One arm's ``antsRegistration`` call (container paths), writing ``<prefix>0GenericAffine.mat``.

    ``init`` is the pull (or, for an ``INVERSE`` arm, the push) as an ITK file; every stage uses
    ``schedule``. ``stages`` replaces the arm's own (a screen variant); with ``masked`` false the
    stages have no masks.
    """
    masks = '[/w/t1w_mask_lr.nii.gz,/w/mask.nii.gz]'
    if arm in INVERSE:
        masks = '[/w/mask.nii.gz,/w/t1w_mask.nii.gz]'

    cmd = [
        f'timeout {TIMEOUT_S} antsRegistration --verbose 1 --dimensionality 3 --float 0',
        f'--collapse-output-transforms 1 --output [{prefix}] --interpolation Linear',
        '--winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0',
        f'--random-seed {SEED} --initial-moving-transform {init}',
    ]
    for terms in stages or ALL_ARMS[arm]:
        cmd.append('--transform Rigid[0.1]')
        cmd += [f'--metric {t}' for t in terms]
        cmd.append(schedule)
        if masked:
            cmd.append(f'--masks {masks}')

    return ' '.join(cmd)


def run_script(work, name, lines):
    """Run ``lines`` as ``<work>/<name>`` in the pinned image, single-threaded, ``work`` at /w."""
    script = ['export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1', *lines]
    (work / name).write_text('\n'.join(script) + '\n')
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
        f'/w/{name}',
    ]
    subprocess.run(cmd, check=True, capture_output=True)


def register(work, arms):
    if not arms:
        return

    lines = [': > /w/ants_status.txt']
    for arm in arms:
        init = '/w/init_itk_inv.txt' if arm in INVERSE else '/w/init_itk.txt'
        lines.append(ants_command(arm, init, f'/w/ants_{arm}_') + f' > /w/ants_{arm}.log 2>&1')
        lines.append(f'echo "{arm} exit $?" >> /w/ants_status.txt')
    # No ``set -e``: an arm that fails or times out has no transform and is reported as such.
    run_script(work, 'ants.sh', lines)


def read_ants(work, arm, f=None):
    """An arm's result as a pull; ``f`` is the output file if not the primary run's."""
    import nitransforms as nt

    f = f or work / f'ants_{arm}_0GenericAffine.mat'
    if not f.exists():
        return None

    matrix = np.asarray(nt.linear.load(f, fmt='itk').matrix)
    return np.linalg.inv(matrix) if arm in INVERSE else matrix


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


def anat_dir(fixture_dir, run_dir):
    """The anatomy the arms use: the fixture's derivatives, or for a ``raw`` fixture ASLPrep's."""
    given = fixture_dir / 'derivatives' / 'anat' / 'sub-01' / 'anat'
    return given if given.exists() else run_dir / 'sub-01' / 'anat'


def run_fixture(name, args):
    work = args.work / name
    if args.ants_only:
        # Reuse the prepared work directory (and its FSL results); register the ANTs arms only.
        register(work, args.arms)
    elif not args.rescore:
        shutil.rmtree(work, ignore_errors=True)
        anat = anat_dir(args.fixtures / name, args.runs / name)
        # Prepares the work directory and runs production's FSL path (both BBR signs).
        coreg_experiment.run_fixture(
            args.fixtures / name,
            args.runs / name,
            work,
            args.license,
            anat_dir=anat,
        )
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
    p.add_argument('--arms', nargs='*', default=list(ARMS), choices=list(ALL_ARMS))
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
