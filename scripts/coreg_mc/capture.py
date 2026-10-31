#!/usr/bin/env python3
"""Capture range of the coregistration arms on fixtures (plan: "Metrics", metric 3).

Every arm starts from ``C0 = C* P`` instead of its own initialization: ``C*`` is the truth pull
(T1w point -> aslref point, ``run_ants.fixture_truth``) and ``P`` a random rigid transform about
the brain's centre, applied on the T1w side, at 2, 5, 10 and 15 mm and degrees (a shift of that
many mm along a random direction, then a rotation of that many degrees about a random axis), 20
draws per level. The draws are drawn once from a fixed seed, so every arm and fixture gets the
same ones. FSL arms run FLIRT-BBR from ``C0`` with ``mri_coreg`` skipped (fMRIPrep's fallback
check then compares BBR with ``C0``); ANTs arms take ``C0`` as ``--initial-moving-transform``.
A draw converges if the result is within 0.5 acquisition voxel of ``C*``. The table also scores
the optimizer's result before the fallback check (``optimizer_*``): the check keeps ``C0`` when
BBR moved more than 15 mm from it, so with it no FSL arm can recover from a start that far off.
``--rescore`` scores results already on disk without running the arms.

``--identity-test`` checks the input conventions instead (gate 1, "the unperturbed start scores
zero"): each arm is given ``C0 = C*`` and run without optimizing (ANTs with zero iterations,
FLIRT with ``-applyxfm``), and its result, read back as the arm's own result is read, must equal
``C*``.

Arms are named as in the probe's tables: ``ANTs <arm>`` (``run_ants.ARMS``) and
``FSL <source> BBR <sign>`` (``coreg_experiment.SOURCES``; sign ``default`` or ``opposite``).
The fixture's work directory must already be prepared by ``run_ants.py``.

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/capture.py --fixtures <aslscan fixtures> --runs <fixture runs> \\
        --work <work dir> --out <tsv> --names rc_pcasl ... --arms 'ANTs S1' \\
        'FSL conventional BBR default' [--identity-test] [--parallel 4]
"""

import argparse
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'hmc_probe'))

import run_ants  # noqa: E402
import scoring  # noqa: E402
from candidates import _fsl_to_world  # noqa: E402
from coreg_experiment import SOURCES, _pull, fmriprep_fallback  # noqa: E402

LEVELS = (2, 5, 10, 15)
DRAWS = 20
SEED = 20261030
CONVERGED_VOXELS = 0.5
#: Zero iterations: ANTs composes and writes its initial transform unchanged.
NO_ITERATIONS = '--convergence [0,1e-6,10] --shrink-factors 1 --smoothing-sigmas 0vox'
#: The result of an identity test must match ``C*`` to this RMS (mm): file round-off only.
IDENTITY_TOLERANCE_MM = 1e-3


def unit(v):
    return v / np.linalg.norm(v)


def axis_rigid(translation, axis, angle_deg, centre):
    """A rotation of ``angle_deg`` about ``axis`` through ``centre``, then a shift."""
    k = unit(np.asarray(axis, float))
    kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    a = np.deg2rad(angle_deg)
    rot = np.eye(3) + np.sin(a) * kx + (1 - np.cos(a)) * kx @ kx
    out = np.eye(4)
    out[:3, :3] = rot
    centre = np.asarray(centre, float)
    out[:3, 3] = centre - rot @ centre + np.asarray(translation, float)
    return out


def draws(levels=LEVELS, n=DRAWS, seed=SEED):
    """[(level, draw, shift direction, rotation axis)], the same for every arm and fixture."""
    rng = np.random.default_rng(seed)
    return [
        (level, i, unit(rng.normal(size=3)), unit(rng.normal(size=3)))
        for level in levels
        for i in range(n)
    ]


def perturbation(level, direction, axis, centre):
    return axis_rigid(direction * level, axis, level, centre)


def fsl_matrix(pull, src_img, t1_img):
    """The FSL matrix (source -> T1w) of a RAS pull: the inverse of ``coreg_experiment._pull``."""
    push = np.linalg.inv(pull)
    return np.linalg.inv(_fsl_to_world(t1_img)) @ push @ _fsl_to_world(src_img)


def parse_arm(arm):
    """('ANTs', arm, None) or ('FSL', source, sign)."""
    parts = arm.split()
    if parts[0] == 'ANTs' and len(parts) == 2 and parts[1] in run_ants.ALL_ARMS:
        return 'ANTs', parts[1], None

    if parts[0] == 'FSL' and len(parts) == 4 and parts[1] in SOURCES and parts[2] == 'BBR':
        if parts[3] in ('default', 'opposite'):
            return 'FSL', parts[1], parts[3]

    raise ValueError(f'unknown arm {arm!r}')


def start_files(work, tag, c0):
    """Write ``C0`` for every arm's convention under ``<work>/capture``; container paths."""
    import nitransforms as nt

    out = work / 'capture'
    nt.linear.Affine(c0).to_filename(out / f'{tag}_itk.txt', fmt='itk')
    nt.linear.Affine(np.linalg.inv(c0)).to_filename(out / f'{tag}_itk_inv.txt', fmt='itk')
    t1_img = nb.load(work / 't1w_brain.nii.gz')
    for image in {v[0] for v in SOURCES.values()}:
        src_img = nb.load(work / f'{image}.nii.gz')
        np.savetxt(out / f'{tag}_{image}.mat', fsl_matrix(c0, src_img, t1_img))


def arm_lines(arm, tag, identity):
    """Shell lines that run one arm from the start ``tag``, writing ``capture/<tag>__<slug>``."""
    kind, name, sign = parse_arm(arm)
    slug = arm.replace(' ', '_')
    result = f'/w/capture/{tag}__{slug}'
    if kind == 'ANTs':
        init = f'/w/capture/{tag}_itk.txt'
        if name in run_ants.INVERSE:
            init = f'/w/capture/{tag}_itk_inv.txt'

        schedule = NO_ITERATIONS if identity else run_ants.SCHEDULE
        cmd = run_ants.ants_command(name, init, f'{result}_', schedule=schedule)
        return [f'{cmd} > {result}.log 2>&1']

    image, extra = SOURCES[name]
    init = f'/w/capture/{tag}_{image}.mat'
    src = f'/w/{image}.nii.gz'
    if identity:
        return [
            f'flirt -in {src} -ref /w/t1w_brain.nii.gz -applyxfm -init {init} '
            f'-omat {result}.mat > {result}.log 2>&1',
        ]

    slope = ['-bbrslope', '0.5'] if sign == 'opposite' else []
    args = ' '.join([*extra, *slope])
    return [
        f'flirt -in {src} -ref /w/t1w_brain.nii.gz -cost bbr -dof 6 -basescale 1 '
        f'-wmseg /w/wm.nii.gz -init {init} -omat {result}.mat -schedule $S {args} '
        f'> {result}.log 2>&1',
    ]


def read_result(work, arm, tag, c0):
    """(pull or None, fallback, optimizer's pull or None) for one arm and start.

    The first is the arm's own result, after fMRIPrep's fallback check for FSL arms; the last is
    the optimizer's result before that check.
    """
    kind, name, _ = parse_arm(arm)
    result = work / 'capture' / f'{tag}__{arm.replace(" ", "_")}'
    if kind == 'ANTs':
        pull = run_ants.read_ants(work, name, Path(f'{result}_0GenericAffine.mat'))
        return pull, False, pull

    mat = Path(f'{result}.mat')
    if not mat.exists():
        return None, False, None

    image = SOURCES[name][0]
    bbr = _pull(mat, nb.load(work / f'{image}.nii.gz'), nb.load(work / 't1w_brain.nii.gz'))
    fallback = fmriprep_fallback(c0, bbr)
    return (c0 if fallback else bbr), fallback, bbr


def run_fixture(name, args):
    work = args.work / name
    truth, points, _production, voxel = run_ants.fixture_truth(
        args.fixtures / name,
        args.runs / name,
    )
    centre = points.mean(0)
    starts = [('identity', 0, 0, np.eye(4))]
    if not args.identity_test:
        starts = [
            (f'L{level}_d{i:02d}', level, i, perturbation(level, direction, axis, centre))
            for level, i, direction, axis in draws()
        ]

    if not args.rescore:
        shutil.rmtree(work / 'capture', ignore_errors=True)
        (work / 'capture').mkdir()
        lines = ['S=$FSLDIR/etc/flirtsch/bbr.sch']
        for tag, _level, _draw, p in starts:
            start_files(work, tag, truth @ p)
            for arm in args.arms:
                lines += arm_lines(arm, tag, args.identity_test)

        run_ants.run_script(work, 'capture.sh', lines)

    rows = []
    for tag, level, draw, p in starts:
        c0 = truth @ p
        for arm in args.arms:
            est, fallback, optimizer = read_result(work, arm, tag, c0)
            row = {
                'fixture': name,
                'arm': arm,
                'level': level,
                'draw': draw,
                'voxel_mm': voxel,
                'start_error_mm': scoring.distance(c0, truth, points),
                'fallback': fallback,
            }
            if est is None:
                row |= {'error_mm': np.nan, 'converged': False, 'problem': 'no transform'}
                row |= {'optimizer_error_mm': np.nan, 'optimizer_converged': False}
            else:
                error = scoring.distance(est, truth, points)
                row |= {'error_mm': error, 'converged': error <= CONVERGED_VOXELS * voxel}
                error = scoring.distance(optimizer, truth, points)
                row |= {
                    'optimizer_error_mm': error,
                    'optimizer_converged': error <= CONVERGED_VOXELS * voxel,
                }
            rows.append(row)
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--fixtures', type=Path, required=True)
    p.add_argument('--runs', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--names', nargs='+', required=True)
    p.add_argument('--arms', nargs='+', required=True)
    p.add_argument('--identity-test', action='store_true')
    p.add_argument('--parallel', type=int, default=4)
    p.add_argument(
        '--rescore',
        action='store_true',
        help='score the results already under <work>/<fixture>/capture without running the arms',
    )
    args = p.parse_args(argv)
    for arm in args.arms:
        parse_arm(arm)

    with ThreadPoolExecutor(args.parallel) as pool:
        rows = [r for rows in pool.map(lambda n: run_fixture(n, args), args.names) for r in rows]
    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        if args.identity_test:
            worst = df.groupby('arm')['error_mm'].max()
            print(worst)
            ok = bool((df['error_mm'] <= IDENTITY_TOLERANCE_MM).all())
            print('IDENTITY PASS' if ok else 'IDENTITY FAIL')
            return 0 if ok else 1

        for column in ('converged', 'optimizer_converged'):
            print(column)
            print(df.pivot_table(index='arm', columns='level', values=column, aggfunc='mean'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
