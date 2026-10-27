#!/usr/bin/env python3
"""Coregistration arms on real panel runs, and their cross-run consistency (no truth).

Plan: Task 1 (harness, real data) and Task 2's engine pilot. For each selected run of
``manifest.tsv``, the HMC probe's ``coreg_experiment.run_real`` prepares the work directory from
the run's ASLPrep outputs and runs production's FSL path (both BBR signs, every D3 source); the
ANTs arms of ``run_ants`` then run from the same ``mri_coreg`` start, with the grey-matter target
taken from ASLPrep's own segmentation of the T1w. Every arm's pull transform (T1w point ->
aslref point, ASLPrep's convention) is saved to ``<work>/transforms.json``.

Consistency (plan: "Real data", design 1). For runs a and b of one session with the same
phase-encoding direction, an image of run a (moving) is registered directly to the same image of
run b (fixed) with ANTs (MI, rigid, from the identity), giving the pull ``D_ab`` (b point -> a
point). An arm with transforms ``C_a``, ``C_b`` is inconsistent by the RMS displacement of
``D_ab^-1 C_a C_b^-1`` over run b's brain mask. The direct registration's own error is bounded by
its noise floor: the RMS displacement of ``D_ab D_ba`` over run a's brain. It is done in two
frames (``FRAMES``): the conventional reference, which depicts the M0 (or first volume) at its own
pose, and the perfusion image built from the corrected series. An arm built from one image shares
that image's frame, so the score in the other frame carries the per-run offset between the two.

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/run_real.py --manifest scripts/coreg_mc/manifest.tsv \\
        --runs /home/tsalo/coreg-mc/runs/panel --adapted /home/tsalo/coreg-mc/adapted \\
        --license <license.txt> --work <dir> --out <prefix> --subjects ds000234/sub-01 ...
"""

import argparse
import csv
import itertools
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
from nibabel.affines import apply_affine

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'hmc_probe'))

import coreg_experiment  # noqa: E402
import run_ants  # noqa: E402
from score_full import _one, _transform_error  # noqa: E402

#: FSL sources kept as arms: production's reference and the D3 level chosen for R2 candidates.
FSL_SOURCES = ('conventional', 'P-winsor')

#: Consistency frames: {frame name: the work-directory image registered between runs}.
FRAMES = {'conventional': 'conventional', 'perfusion': 'P-winsor'}


def comparable_key(run):
    """Runs of one session that can be paired: the same phase-encoding direction."""
    pe = re.search(r'_dir-([A-Za-z]+)', run)
    return pe.group(1) if pe else ''


def record_for(row):
    """The record ``coreg_experiment.run_real`` expects, from a manifest row."""
    root = Path(row['root'])
    sub_dir = root / row['subject']
    candidates = sub_dir.rglob('perf/*_asl.nii*')
    asl = next(p for p in candidates if p.name.startswith(row['run'] + '_asl'))
    context = asl.with_name(re.sub(r'_asl\.nii(\.gz)?$', '_aslcontext.tsv', asl.name))
    return {
        'dataset': row['dataset'],
        'subject': row['subject'].removeprefix('sub-'),
        'run': row['run'],
        'family': row['family'],
        'root': str(root),
        'asl_path': str(asl.relative_to(root)),
        'aslcontext': str(context.relative_to(root)),
    }


def run_one(row, args):
    """Register every arm on one run; returns its work directory, {arm: pull}, FSL fallbacks."""
    record = record_for(row)
    work = args.work / f'{row["dataset"]}__{row["run"]}'
    if not args.rescore:
        shutil.rmtree(work, ignore_errors=True)
        coreg_experiment.run_real(record, args.runs, args.adapted, work, args.license)
        run_dir = args.runs / f'{record["dataset"]}_sub-{record["subject"]}'
        anat = sorted(run_dir.glob('sub-*/anat')) + sorted(run_dir.glob('sub-*/ses-*/anat'))
        gm = next(p for p in anat[0].glob('*_label-GM_probseg.nii.gz') if '_space-' not in p.name)
        run_ants.gm_target(work, gm)
        run_ants.lowres(work)
        run_ants.write_init(work)
        run_ants.register(work, args.arms)

    arms = {}
    fallbacks = set()
    run_dir = args.runs / f'{record["dataset"]}_sub-{record["subject"]}'
    perf = run_dir / Path(record['asl_path']).parent
    stem = re.sub(r'_echo-[^_]+', '', Path(record['asl_path']).name.split('_asl.')[0])
    xfm = sorted(perf.glob(f'{stem}_from-aslref_to-T1w*_xfm.txt'))
    if xfm:
        arms['R0 production'] = _one(xfm[0])

    for (source, sign), (final, _bbr, _init, fallback) in coreg_experiment.estimates(work).items():
        if source in FSL_SOURCES:
            arms[f'FSL {source} {sign}'] = final
            if fallback:
                fallbacks.add(f'FSL {source} {sign}')

    for arm in args.arms:
        est = run_ants.read_ants(work, arm)
        if est is not None:
            arms[f'ANTs {arm}'] = est

    record = {k: np.asarray(v).tolist() for k, v in arms.items()}
    record['fallbacks'] = sorted(fallbacks)
    (work / 'transforms.json').write_text(json.dumps(record, indent=1) + '\n')
    return work, arms, fallbacks


def brain_points(work, every=5):
    """World coordinates of every ``every``-th voxel of the run's ASL brain mask."""
    mask = nb.load(work / 'mask.nii.gz')
    ijk = np.argwhere(np.asanyarray(mask.dataobj) > 0)[::every]
    return apply_affine(mask.affine, ijk)


def direct(work_a, work_b, out, image, reuse=False):
    """ANTs rigid MI of run a's ``image`` (moving) to run b's (fixed): the pull b -> a."""
    import nitransforms as nt

    result = out / 'direct_0GenericAffine.mat'
    if reuse and result.exists():
        return np.asarray(nt.linear.load(result, fmt='itk').matrix)

    out.mkdir(parents=True, exist_ok=True)
    for name, src in (('a', work_a), ('b', work_b)):
        shutil.copyfile(src / f'{image}.nii.gz', out / f'{name}_image.nii.gz')
        shutil.copyfile(src / 'mask.nii.gz', out / f'{name}_mask.nii.gz')

    cmd_ants = (
        f'timeout {run_ants.TIMEOUT_S} antsRegistration --dimensionality 3 --float 0 '
        '--collapse-output-transforms 1 --output [/w/direct_] --interpolation Linear '
        '--winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 '
        f'--random-seed {run_ants.SEED} --transform Rigid[0.1] '
        '--metric MI[/w/b_image.nii.gz,/w/a_image.nii.gz,1,32,Regular,0.5] '
        '--convergence [500x250,1e-6,10] --shrink-factors 2x1 --smoothing-sigmas 1x0vox '
        '--masks [/w/b_mask.nii.gz,/w/a_mask.nii.gz] > /w/direct.log 2>&1'
    )
    (out / 'direct.sh').write_text(f'export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1\n{cmd_ants}\n')
    cmd = [
        'docker',
        'run',
        '--rm',
        '-u',
        f'{os.getuid()}:{os.getgid()}',
        '-e',
        'HOME=/tmp',
        '-v',
        f'{out}:/w',
        '--entrypoint',
        'bash',
        run_ants.IMAGE,
        '/w/direct.sh',
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    return np.asarray(nt.linear.load(result, fmt='itk').matrix)


def consistency(runs, args):
    """Inconsistency of every arm over every comparable pair of runs in one session."""
    rows = []
    for (row_a, work_a, arms_a, fb_a), (row_b, work_b, arms_b, fb_b) in itertools.permutations(
        runs,
        2,
    ):
        same_session = (row_a['dataset'], row_a['subject'], row_a['session']) == (
            row_b['dataset'],
            row_b['subject'],
            row_b['session'],
        )
        if not same_session or comparable_key(row_a['run']) != comparable_key(row_b['run']):
            continue

        points_a = brain_points(work_a)
        points_b = brain_points(work_b)
        for frame, image in FRAMES.items():
            frame_dir = args.work / 'direct' / frame
            pair_dir = frame_dir / f'{work_a.name}__to__{work_b.name}'
            d_ab = direct(work_a, work_b, pair_dir, image, reuse=args.rescore)
            reverse_dir = frame_dir / f'{work_b.name}__to__{work_a.name}'
            d_ba = direct(work_b, work_a, reverse_dir, image, reuse=args.rescore)
            floor = _transform_error(d_ab @ d_ba, np.eye(4), points_a)['rms_mm']
            for arm in sorted(set(arms_a) & set(arms_b)):
                implied = arms_a[arm] @ np.linalg.inv(arms_b[arm])
                err = _transform_error(d_ab, implied, points_b)['rms_mm']
                rows.append(
                    {
                        'dataset': row_a['dataset'],
                        'subject': row_a['subject'],
                        'session': row_a['session'],
                        'run_a': row_a['run'],
                        'run_b': row_b['run'],
                        'frame': frame,
                        'arm': arm,
                        'inconsistency_mm': err,
                        'floor_mm': floor,
                        'fallback': arm in fb_a or arm in fb_b,
                    },
                )
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--runs', type=Path, required=True)
    p.add_argument('--adapted', type=Path, required=True)
    p.add_argument('--license', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True, help='Output prefix for the TSVs.')
    p.add_argument('--subjects', nargs='+', required=True, help='dataset/sub-XX')
    p.add_argument('--arms', nargs='+', default=list(run_ants.ARMS))
    p.add_argument('--rescore', action='store_true')
    args = p.parse_args(argv)

    with args.manifest.open() as f:
        manifest = list(csv.DictReader(f, delimiter='\t'))
    wanted = {tuple(s.split('/')) for s in args.subjects}
    runs = []
    for row in manifest:
        if (row['dataset'], row['subject']) not in wanted:
            continue

        work, arms, fallbacks = run_one(row, args)
        runs.append((row, work, arms, fallbacks))
        print(f'{row["run"]}: {len(arms)} arms', flush=True)

    # Agreement of every arm with every other on each run (RMS over the run's brain, mm).
    agree = []
    for row, work, arms, _fallbacks in runs:
        points = brain_points(work)
        for a, b in itertools.combinations(sorted(arms), 2):
            dist = _transform_error(arms[a], arms[b], points)['rms_mm']
            agree.append({'run': row['run'], 'arm_a': a, 'arm_b': b, 'distance_mm': dist})
    pd.DataFrame(agree).to_csv(f'{args.out}_agreement.tsv', sep='\t', index=False)

    cons = pd.DataFrame(consistency(runs, args))
    cons.to_csv(f'{args.out}_consistency.tsv', sep='\t', index=False)
    if not cons.empty:
        table = cons.pivot_table(
            index='arm',
            columns=['dataset', 'frame'],
            values='inconsistency_mm',
            aggfunc='median',
        )
        floors = cons.groupby(['dataset', 'frame'])['floor_mm'].median()
        one_frame = cons[cons['frame'] == next(iter(FRAMES))]
        fell_back = one_frame[one_frame['fallback']].groupby('arm').size()
        with pd.option_context('display.width', 200, 'display.max_rows', 200):
            print(table.round(2))
            print('noise floor (median):', floors.round(2).to_dict())
            print('ordered pairs with an FSL fallback:', fell_back.to_dict())


if __name__ == '__main__':
    main()
