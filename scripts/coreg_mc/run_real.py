#!/usr/bin/env python3
"""Coregistration arms on real panel runs, and their truth-free scores.

Plan: Task 1 (harness, real data) and Task 2 (scorer). For each selected run of ``manifest.tsv``,
the HMC probe's ``coreg_experiment.run_real`` prepares the work directory from the run's ASLPrep
outputs and runs production's FSL path (both BBR signs, every D3 source); the ANTs arms of
``run_ants`` then run from the same ``mri_coreg`` start, with the grey-matter target taken from
ASLPrep's own segmentation of the T1w. Every arm's pull transform (T1w point -> aslref point,
ASLPrep's convention) is saved to ``<work>/transforms.json``. The scores are in ``scoring``.

Consistency (plan: "Real data", designs 1 and 2). For runs a and b of one subject with the same
acquisition (``comparable_key``) and anatomical reference (``pairable``), an image of run a
(moving) is registered directly to the same image of run b (fixed) with ANTs (MI, rigid), giving
the pull ``D_ab``, and with FLIRT as a second engine; pairs within a session start ANTs from the
identity, pairs across sessions (the head was repositioned) from the images' centres of mass.
Cross-session pairs need the subject-level T1w template, which every session's
``from-aslref_to-T1w`` refers to. The direct registration is done in two frames (``FRAMES``): the
conventional reference, which depicts the M0 (or first volume) at its own pose, and the perfusion
image built from the corrected series. Expected arms with no transform are failures to run
(``missing``) and count as gross failures.

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
import hashlib
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'hmc_probe'))

import coreg_experiment  # noqa: E402
import run_ants  # noqa: E402
import scoring  # noqa: E402
from score_full import _one  # noqa: E402

#: FSL sources kept as arms: production's reference and every D3 level (R2's level and
#: eligibility rule are frozen from them in Task 3).
FSL_SOURCES = tuple(coreg_experiment.SOURCES)

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


#: R0 (ASLPrep's stored transform) and the harness's reproduction of production's FSL path agree
#: to within this on every run checked; a larger distance means the harness is not looking at the
#: same anatomy or transform as ASLPrep, and the run's R0 is not comparable.
R0_REPRODUCTION_MM = 0.05


def comparable_key(row, work):
    """Runs of one subject that can be paired: the same acquisition, judged from the data.

    The phase-encoding direction and acquisition type from the ASL sidecar, the manifest's readout
    and background suppression, and the reference's grid (shape and voxel size). Runs that differ
    in any of these differ in distortion or contrast, which a direct registration would count as
    inconsistency.
    """
    record = record_for(row)
    sidecar = Path(record['root']) / record['asl_path']
    sidecar = sidecar.with_name(re.sub(r'\.nii(\.gz)?$', '.json', sidecar.name))
    meta = json.loads(sidecar.read_text()) if sidecar.exists() else {}
    ref = nb.load(work / 'conventional.nii.gz')
    zooms = tuple(round(float(z), 2) for z in ref.header.get_zooms()[:3])
    return (
        meta.get('PhaseEncodingDirection'),
        meta.get('MRAcquisitionType'),
        row['readout'],
        row['background_suppression'],
        ref.shape[:3],
        zooms,
    )


def anatomical_reference(run_dir):
    """The anatomical outputs the harness registers to: the subject-level T1w when there is one.

    ``coreg_experiment.run_real`` takes the first of these directories; with several sessions,
    ASLPrep writes one subject-level template (``sub-X/anat``) that every session's
    ``from-aslref_to-T1w`` refers to. Returns (directory, whether it is subject-level).
    """
    anat = sorted(run_dir.glob('sub-*/anat')) + sorted(run_dir.glob('sub-*/ses-*/anat'))
    if not anat:
        raise FileNotFoundError(f'no anat outputs in {run_dir}')
    return anat[0], not anat[0].parent.name.startswith('ses-')


def expected_arms(ants_arms):
    """Every arm a run should have: R0, the FSL arms and the ANTs arms."""
    fsl = [f'FSL {s} {sign}' for s in FSL_SOURCES for sign in ('default', 'opposite')]
    return ['R0 production', *fsl, *(f'ANTs {a}' for a in ants_arms)]


def fingerprint(paths, command):
    """A digest of input files and a command, to decide whether a stored result is current."""
    digest = hashlib.sha256(command.encode())
    for path in paths:
        digest.update(Path(path).read_bytes())
    return digest.hexdigest()


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
    """Register every arm on one run (``--rescore``: reuse the work directory as it is).

    Returns a dict: ``row``, ``work``, ``arms`` ({arm: pull}), ``fallbacks``, ``missing`` (expected
    arms with no transform: failures to run), ``excluded`` (R0, when its reproduction differs by
    more than ``R0_REPRODUCTION_MM``), ``key`` (``comparable_key``), ``anat`` and
    ``subject_level`` (``anatomical_reference``), and ``r0_vs_reproduction_mm``.
    """
    record = record_for(row)
    work = args.work / f'{row["dataset"]}__{row["run"]}'
    run_dir = args.runs / f'{record["dataset"]}_sub-{record["subject"]}'
    anat, subject_level = anatomical_reference(run_dir)
    if not args.rescore:
        shutil.rmtree(work, ignore_errors=True)
        coreg_experiment.run_real(record, args.runs, args.adapted, work, args.license)
        gm = next(p for p in anat.glob('*_label-GM_probseg.nii.gz') if '_space-' not in p.name)
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
    missing = sorted(set(expected_arms(args.arms)) - set(arms))
    r0_vs = np.nan
    reproduction = arms.get('FSL conventional default')
    if 'R0 production' in arms and reproduction is not None:
        r0_vs = scoring.distance(arms['R0 production'], reproduction, brain_points(work))

    # An R0 its reproduction does not match is not comparable with the other arms (the harness is
    # not looking at the same anatomy or transform): it is excluded, not scored and not a failure.
    excluded = []
    if r0_vs > R0_REPRODUCTION_MM:
        arms.pop('R0 production')
        excluded.append('R0 production')

    stored = {k: np.asarray(v).tolist() for k, v in arms.items()}
    stored |= {'fallbacks': sorted(fallbacks), 'missing': missing, 'excluded': excluded}
    stored['anat'] = str(anat)
    stored['r0_vs_reproduction_mm'] = r0_vs
    (work / 'transforms.json').write_text(json.dumps(stored, indent=1) + '\n')
    return {
        'row': row,
        'work': work,
        'arms': arms,
        'fallbacks': fallbacks,
        'missing': missing,
        'excluded': excluded,
        'key': comparable_key(row, work),
        'anat': anat,
        'subject_level': subject_level,
        'r0_vs_reproduction_mm': r0_vs,
    }


def brain_points(work, every=5):
    """World coordinates of every ``every``-th voxel of the run's ASL brain mask."""
    return scoring.mask_points(work / 'mask.nii.gz', every)


def direct(work_a, work_b, out, image, reuse=False, centre_init=False):
    """Run a's ``image`` (moving) registered directly to run b's (fixed) by two engines.

    Returns the pulls b -> a from ANTs (rigid MI, masked) and FLIRT (6 dof, normalized MI);
    their disagreement is part of the noise floor.
    ``centre_init`` starts ANTs from the images' centres of mass instead of the identity (FLIRT
    runs its own search). With ``reuse``, a stored result is used only if its fingerprint (the four
    input images and the commands) matches.
    """
    import nitransforms as nt

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
    # FLIRT's default +-90 degree search found mirrored solutions (117 mm from the identity for an
    # image registered to itself), and brain-mask weights made one pair worse (1.6 mm from ANTs,
    # against 0.6 without), so the search is limited to +-15 degrees and the images are unweighted.
    search = ' '.join(f'-searchr{axis} -15 15' for axis in 'xyz')
    cmd_flirt = (
        f'flirt -in /w/a_image.nii.gz -ref /w/b_image.nii.gz -dof 6 -cost normmi {search} '
        '-omat /w/flirt.mat > /w/flirt.log 2>&1'
    )
    script = f'export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1\n{cmd_ants}\n{cmd_flirt}\n'
    inputs = [src / f for src in (work_a, work_b) for f in (f'{image}.nii.gz', 'mask.nii.gz')]
    stamp = fingerprint(inputs, script + run_ants.IMAGE)
    ants_file, flirt_file = out / 'direct_0GenericAffine.mat', out / 'flirt.mat'

    def read():
        d_ants = np.asarray(nt.linear.load(ants_file, fmt='itk').matrix)
        a_img, b_img = nb.load(out / 'a_image.nii.gz'), nb.load(out / 'b_image.nii.gz')
        return d_ants, coreg_experiment._pull(flirt_file, a_img, b_img)

    stamp_file = out / 'fingerprint.txt'
    if reuse and stamp_file.exists() and stamp_file.read_text().strip() == stamp:
        return read()

    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    for name, src in (('a', work_a), ('b', work_b)):
        shutil.copyfile(src / f'{image}.nii.gz', out / f'{name}_image.nii.gz')
        shutil.copyfile(src / 'mask.nii.gz', out / f'{name}_mask.nii.gz')

    (out / 'direct.sh').write_text(script)
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
    stamp_file.write_text(stamp + '\n')
    return read()


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


def pairable(a, b):
    """Whether runs a and b form a pair, and of which type (``None`` when they do not).

    Same subject, the same acquisition (``comparable_key``) and the same anatomical reference;
    across sessions, that reference must be the subject-level T1w that every session refers to.
    """
    row_a, row_b = a['row'], b['row']
    if (row_a['dataset'], row_a['subject']) != (row_b['dataset'], row_b['subject']):
        return None

    if a['key'] != b['key'] or a['anat'] != b['anat']:
        return None

    if row_a['session'] == row_b['session']:
        return 'cross-run'
    return 'cross-session' if a['subject_level'] else None


def consistency(runs, args):
    """Inconsistency of every arm over every pair of a subject's runs (see ``pairable``).

    The noise floor is the larger of the forward-reverse closure of the ANTs direct registration
    and the disagreement between the ANTs and FLIRT direct registrations.
    """
    rows = []
    for a, b in itertools.permutations(runs, 2):
        pair_type = pairable(a, b)
        if pair_type is None:
            continue

        work_a, work_b = a['work'], b['work']
        points_a = brain_points(work_a)
        points_b = brain_points(work_b)
        options = {'reuse': args.rescore, 'centre_init': pair_type == 'cross-session'}
        for frame, image in FRAMES.items():
            forward = args.work / 'direct' / frame / f'{work_a.name}__to__{work_b.name}'
            reverse = args.work / 'direct' / frame / f'{work_b.name}__to__{work_a.name}'
            d_ab, d_ab_flirt = direct(work_a, work_b, forward, image, **options)
            d_ba, _ = direct(work_b, work_a, reverse, image, **options)
            closure = scoring.noise_floor(d_ab, d_ba, points_a)
            engines = scoring.distance(d_ab, d_ab_flirt, points_b)
            base = {
                'dataset': a['row']['dataset'],
                'subject': a['row']['subject'],
                'session_a': a['row']['session'],
                'session_b': b['row']['session'],
                'run_a': a['row']['run'],
                'run_b': b['row']['run'],
                'pair_type': pair_type,
                'frame': frame,
                'floor_closure_mm': closure,
                'floor_engines_mm': engines,
                'floor_mm': max(closure, engines),
            }
            rows += [
                base
                | {
                    'arm': arm,
                    'inconsistency_mm': scoring.inconsistency(d_ab, c_a, b['arms'][arm], points_b),
                    'fallback': arm in a['fallbacks'] or arm in b['fallbacks'],
                }
                for arm, c_a in sorted(a['arms'].items())
                if arm in b['arms']
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
    """Per run and arm: the carried-mask rule (both calibrations) and the inconsistency rules.

    ``verdict`` is the adopted rule (``scoring.gross_verdict``); the other rule columns are kept
    as diagnostics. An expected arm with no transform (it failed or timed out) is a gross failure
    under every rule (``no_transform``), so failures to run are counted, not dropped.
    """
    limit = scoring.GROSS_INCONSISTENCY_MM
    rows = []
    for run in runs:
        row, work, arms = run['row'], run['work'], run['arms']
        base = {
            'dataset': row['dataset'],
            'run': row['run'],
            'r0_vs_reproduction_mm': run['r0_vs_reproduction_mm'],
        }
        carried = scoring.CarriedMask(work / 'mask.nii.gz', work / 't1w_mask.nii.gz')
        kinds = ('axes', 'rigid')
        by_kind = {kind: scoring.gross_by_mask(carried, arms, kind=kind) for kind in kinds}
        for arm in arms:
            partners = partner_inconsistency(cons_rows, row['dataset'], row['run'], arm)
            values = list(partners.values())
            axes, rigid = by_kind['axes'][arm], by_kind['rigid'][arm]
            rows.append(
                base
                | {
                    'arm': arm,
                    'no_transform': False,
                    'excluded': False,
                    'dice': axes['dice'],
                    'best_dice': axes['best_dice'],
                    'best_arm': axes['best_arm'],
                    'limit_axes': axes['limit'],
                    'limit_rigid': rigid['limit'],
                    'gross_mask_axes': axes['gross_mask'],
                    'gross_mask_rigid': rigid['gross_mask'],
                    'n_partners': len(values),
                    'max_inconsistency_mm': max(values, default=np.nan),
                    'gross_any_partner': any(v > limit for v in values),
                    'gross_every_partner': bool(values) and all(v > limit for v in values),
                    'verdict': scoring.gross_verdict(values, rigid['gross_mask']),
                },
            )

        failed = {
            'no_transform': True,
            'excluded': False,
            'gross_mask_axes': True,
            'gross_mask_rigid': True,
            'gross_any_partner': True,
            'gross_every_partner': True,
            'verdict': 'fail',
        }
        rows += [base | {'arm': arm} | failed for arm in run['missing']]
        rows += [base | {'arm': arm, 'excluded': True} for arm in run['excluded']]
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
    p.add_argument(
        '--arms',
        nargs='*',
        default=list(run_ants.ARMS),
        help='ANTs arms; none (an empty list) scores the FSL arms and R0 only.',
    )
    p.add_argument(
        '--rescore',
        action='store_true',
        help=(
            'Score the work directories as they are, without registering again (they are not '
            'fingerprinted: rerun without this after any change to the inputs or arms); direct '
            'registrations are reused only when their fingerprint matches.'
        ),
    )
    args = p.parse_args(argv)

    with args.manifest.open() as f:
        manifest = list(csv.DictReader(f, delimiter='\t'))
    wanted = {tuple(s.split('/')) for s in args.subjects}
    runs = []
    failed = []
    for row in manifest:
        if (row['dataset'], row['subject']) not in wanted:
            continue

        try:
            run = run_one(row, args)
        except Exception as err:  # noqa: BLE001 - one bad run is recorded, not fatal to a batch
            failed.append({'dataset': row['dataset'], 'run': row['run'], 'error': repr(err)})
            print(f'{row["run"]}: FAILED {err!r}', flush=True)
            continue

        runs.append(run)
        r0_vs = run['r0_vs_reproduction_mm']
        note = ''
        if run['excluded']:
            note = f'; R0 excluded: {r0_vs:.2f} mm from its reproduction'

        print(f'{row["run"]}: {len(run["arms"])} arms, missing {run["missing"]}{note}', flush=True)

    pd.DataFrame(failed, columns=['dataset', 'run', 'error']).to_csv(
        f'{args.out}_failed.tsv',
        sep='\t',
        index=False,
    )
    # Agreement of every arm with every other on each run (RMS over the run's brain, mm).
    agree = []
    for run in runs:
        arms = run['arms']
        points = brain_points(run['work'])
        for a, b in itertools.combinations(sorted(arms), 2):
            dist = scoring.distance(arms[a], arms[b], points)
            agree.append({'run': run['row']['run'], 'arm_a': a, 'arm_b': b, 'distance_mm': dist})
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
            floor_columns = ['floor_closure_mm', 'floor_engines_mm', 'floor_mm']
            floors = cons.groupby(['pair_type', 'frame'])[floor_columns].median()
            print('noise floor (median):')
            print(floors.round(2))

        rules = [
            'no_transform',
            'gross_mask_axes',
            'gross_mask_rigid',
            'gross_any_partner',
            'gross_every_partner',
        ]
        print(gross.groupby('arm')[rules].sum())
        print('adopted rule (verdict per run):')
        print(pd.crosstab(gross['arm'], gross['verdict']))


if __name__ == '__main__':
    main()
