#!/usr/bin/env python3
"""Coregistration arms on real panel runs, and their truth-free scores.

Plan: Task 1 (harness, real data) and Task 2 (scorer). For each selected run of ``manifest.tsv``,
the HMC probe's ``coreg_experiment.run_real`` prepares the work directory from the run's ASLPrep
outputs and runs production's FSL path (both BBR signs, every D3 source); the ANTs arms of
``run_ants`` then run from the same ``mri_coreg`` start, with the grey-matter target taken from
ASLPrep's own segmentation of the T1w. Every arm's pull transform (T1w point -> aslref point,
ASLPrep's convention) is saved to ``<work>/transforms.json``. The scores are in ``scoring``.

Consistency (plan: "Real data", designs 1 and 2). For runs a and b of one subject with the same
phase-encoding direction, an image of run a (moving) is registered directly to the same image of
run b (fixed) with ANTs (MI, rigid), giving the pull ``D_ab``; pairs within a session start from
the identity, pairs across sessions (the head was repositioned) from the images' centres of mass.
Every session of a subject is coregistered to one anatomical reference (ASLPrep's subject-level
T1w template). The direct registration is done in two frames (``FRAMES``): the conventional
reference, which depicts the M0 (or first volume) at its own pose, and the perfusion image built
from the corrected series.

Gross failure (plan: "Gross failure"), per run and arm: the carried-mask rule of ``scoring``, and
inconsistency above ``scoring.GROSS_INCONSISTENCY_MM`` in the conventional frame, recorded two
ways: with any partner run (the plan's rule) and with every partner (which, with two or more
partners, points at the run that failed rather than at both).

Semi-synthetic recovery (design 3) is ``semi_synthetic``: the T1w-side images get a moved
header, production's chain and the ANTs arms run again, and each arm is scored against its own
solution on the original.

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

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'hmc_probe'))

import coreg_experiment  # noqa: E402
import run_ants  # noqa: E402
import scoring  # noqa: E402
from score_full import _one  # noqa: E402

#: FSL sources kept as arms: production's reference and the D3 level chosen for R2 candidates.
FSL_SOURCES = ('conventional', 'P-winsor')

#: Consistency frames: {frame name: the work-directory image registered between runs}.
FRAMES = {'conventional': 'conventional', 'perfusion': 'P-winsor'}

#: Images on the T1w side of the work directory, moved together by ``semi_synthetic``.
T1W_SIDE = (
    't1w_brain',
    't1w_mask',
    'wm',
    'gm_target',
    't1w_brain_lr',
    't1w_mask_lr',
    'gm_target_lr',
)
#: Images on the ASL side, copied unchanged.
ASL_SIDE = ('conventional', 'mask', 'P-mean', 'P-robust', 'P-winsor')


def comparable_key(run):
    """Runs of one subject that can be paired: the same phase-encoding direction."""
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


def collect(work, ants_arms):
    """Every FSL and ANTs arm's pull in a work directory, and the FSL arms that fell back."""
    arms = {}
    fallbacks = set()
    for (source, sign), (final, _bbr, _init, fallback) in coreg_experiment.estimates(work).items():
        if source in FSL_SOURCES:
            arms[f'FSL {source} {sign}'] = final
            if fallback:
                fallbacks.add(f'FSL {source} {sign}')

    for arm in ants_arms:
        est = run_ants.read_ants(work, arm)
        if est is not None:
            arms[f'ANTs {arm}'] = est
    return arms, fallbacks


def run_one(row, args):
    """Register every arm on one run; returns its work directory, {arm: pull}, FSL fallbacks."""
    record = record_for(row)
    work = args.work / f'{row["dataset"]}__{row["run"]}'
    run_dir = args.runs / f'{record["dataset"]}_sub-{record["subject"]}'
    if not args.rescore:
        shutil.rmtree(work, ignore_errors=True)
        coreg_experiment.run_real(record, args.runs, args.adapted, work, args.license)
        anat = sorted(run_dir.glob('sub-*/anat')) + sorted(run_dir.glob('sub-*/ses-*/anat'))
        gm = next(p for p in anat[0].glob('*_label-GM_probseg.nii.gz') if '_space-' not in p.name)
        run_ants.gm_target(work, gm)
        run_ants.lowres(work)
        run_ants.write_init(work)
        run_ants.register(work, args.arms)

    arms = {}
    perf = run_dir / Path(record['asl_path']).parent
    stem = re.sub(r'_echo-[^_]+', '', Path(record['asl_path']).name.split('_asl.')[0])
    xfm = sorted(perf.glob(f'{stem}_from-aslref_to-T1w*_xfm.txt'))
    if xfm:
        arms['R0 production'] = _one(xfm[0])

    found, fallbacks = collect(work, args.arms)
    arms |= found
    record = {k: np.asarray(v).tolist() for k, v in arms.items()}
    record['fallbacks'] = sorted(fallbacks)
    (work / 'transforms.json').write_text(json.dumps(record, indent=1) + '\n')
    return work, arms, fallbacks


def brain_points(work, every=5):
    """World coordinates of every ``every``-th voxel of the run's ASL brain mask."""
    return scoring.mask_points(work / 'mask.nii.gz', every)


def direct(work_a, work_b, out, image, reuse=False, centre_init=False):
    """ANTs rigid MI of run a's ``image`` (moving) to run b's (fixed): the pull b -> a.

    ``centre_init`` starts from the images' centres of mass instead of the identity.
    """
    import nitransforms as nt

    result = out / 'direct_0GenericAffine.mat'
    if reuse and result.exists():
        return np.asarray(nt.linear.load(result, fmt='itk').matrix)

    out.mkdir(parents=True, exist_ok=True)
    for name, src in (('a', work_a), ('b', work_b)):
        shutil.copyfile(src / f'{image}.nii.gz', out / f'{name}_image.nii.gz')
        shutil.copyfile(src / 'mask.nii.gz', out / f'{name}_mask.nii.gz')

    init = '--initial-moving-transform [/w/b_image.nii.gz,/w/a_image.nii.gz,1] '
    cmd_ants = (
        f'timeout {run_ants.TIMEOUT_S} antsRegistration --dimensionality 3 --float 0 '
        '--collapse-output-transforms 1 --output [/w/direct_] --interpolation Linear '
        '--winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 '
        f'--random-seed {run_ants.SEED} {init if centre_init else ""}--transform Rigid[0.1] '
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


def semi_synthetic(work, m, out, ants_arms, license_file):
    """Every arm on the run with its T1w side moved by ``m``; returns (arms, fallbacks).

    The ASL images are copied unchanged and the T1w-side images get the header ``m A``, so the
    expected pull of an arm whose solution on the original is ``C`` is ``C m^-1``. Production's
    chain runs again from its own ``mri_coreg`` initialization, as do the ANTs arms.
    """
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    for name in ASL_SIDE:
        shutil.copyfile(work / f'{name}.nii.gz', out / f'{name}.nii.gz')

    for name in T1W_SIDE:
        scoring.move_header(work / f'{name}.nii.gz', out / f'{name}.nii.gz', m)

    coreg_experiment.register(out, license_file)
    run_ants.write_init(out)
    run_ants.register(out, ants_arms)
    return collect(out, ants_arms)


def consistency(runs, args):
    """Inconsistency of every arm over every comparable pair of a subject's runs."""
    rows = []
    for (row_a, work_a, arms_a, fb_a), (row_b, work_b, arms_b, fb_b) in itertools.permutations(
        runs,
        2,
    ):
        same_subject = (row_a['dataset'], row_a['subject']) == (row_b['dataset'], row_b['subject'])
        if not same_subject or comparable_key(row_a['run']) != comparable_key(row_b['run']):
            continue

        pair_type = 'cross-run' if row_a['session'] == row_b['session'] else 'cross-session'
        points_a = brain_points(work_a)
        points_b = brain_points(work_b)
        options = {'reuse': args.rescore, 'centre_init': pair_type == 'cross-session'}
        for frame, image in FRAMES.items():
            forward = args.work / 'direct' / frame / f'{work_a.name}__to__{work_b.name}'
            reverse = args.work / 'direct' / frame / f'{work_b.name}__to__{work_a.name}'
            d_ab = direct(work_a, work_b, forward, image, **options)
            d_ba = direct(work_b, work_a, reverse, image, **options)
            floor = scoring.noise_floor(d_ab, d_ba, points_a)
            rows += [
                {
                    'dataset': row_a['dataset'],
                    'subject': row_a['subject'],
                    'session_a': row_a['session'],
                    'session_b': row_b['session'],
                    'run_a': row_a['run'],
                    'run_b': row_b['run'],
                    'pair_type': pair_type,
                    'frame': frame,
                    'arm': arm,
                    'inconsistency_mm': scoring.inconsistency(d_ab, c_a, arms_b[arm], points_b),
                    'floor_mm': floor,
                    'fallback': arm in fb_a or arm in fb_b,
                }
                for arm, c_a in sorted(arms_a.items())
                if arm in arms_b
            ]
    return rows


def partner_inconsistency(rows, dataset, run, arm):
    """{partner run: the larger inconsistency of the two ordered pairs}, conventional frame."""
    frame = next(iter(FRAMES))
    out = {}
    for r in rows:
        if (r['dataset'], r['arm'], r['frame']) != (dataset, arm, frame):
            continue

        if run in (r['run_a'], r['run_b']):
            partner = r['run_b'] if r['run_a'] == run else r['run_a']
            out[partner] = max(out.get(partner, 0.0), r['inconsistency_mm'])
    return out


def gross_failures(runs, cons_rows):
    """Per run and arm: the carried-mask rule and the two inconsistency rules."""
    limit = scoring.GROSS_INCONSISTENCY_MM
    rows = []
    for row, work, arms, _fallbacks in runs:
        carried = scoring.CarriedMask(work / 'mask.nii.gz', work / 't1w_mask.nii.gz')
        for arm, result in scoring.gross_by_mask(carried, arms).items():
            partners = partner_inconsistency(cons_rows, row['dataset'], row['run'], arm)
            values = list(partners.values())
            rows.append(
                {
                    'dataset': row['dataset'],
                    'run': row['run'],
                    'arm': arm,
                    **result,
                    'n_partners': len(values),
                    'max_inconsistency_mm': max(values, default=np.nan),
                    'gross_any_partner': any(v > limit for v in values),
                    'gross_every_partner': bool(values) and all(v > limit for v in values),
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
            dist = scoring.distance(arms[a], arms[b], points)
            agree.append({'run': row['run'], 'arm_a': a, 'arm_b': b, 'distance_mm': dist})
    pd.DataFrame(agree).to_csv(f'{args.out}_agreement.tsv', sep='\t', index=False)

    cons_rows = consistency(runs, args)
    cons = pd.DataFrame(cons_rows)
    cons.to_csv(f'{args.out}_consistency.tsv', sep='\t', index=False)
    gross = pd.DataFrame(gross_failures(runs, cons_rows))
    gross.to_csv(f'{args.out}_gross.tsv', sep='\t', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        if not cons.empty:
            table = cons.pivot_table(
                index='arm',
                columns=['pair_type', 'frame'],
                values='inconsistency_mm',
                aggfunc='median',
            )
            print(table.round(2))
            floors = cons.groupby(['pair_type', 'frame'])['floor_mm'].median()
            print('noise floor (median):', floors.round(2).to_dict())

        rules = ['gross_mask', 'gross_any_partner', 'gross_every_partner']
        print(gross.groupby('arm')[rules].sum())


if __name__ == '__main__':
    main()
