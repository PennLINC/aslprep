#!/usr/bin/env python3
"""Run the harness (``harness.py``) over a panel manifest, one container per run.

This is the cheap path for real data: reference building, HMC and subtraction only, with no
anatomical processing. It serves the motion screen and the HMC metrics that need no
coregistration. Each run's output goes to ``<out>/<label>/<dataset>/<run>/`` together with
``motion.json``, a summary of the estimated motion (no truth is needed).

The data roots come from the manifest (``root``); the adapted copies of
``panel_config.ADAPTATIONS`` replace the zoo's datasets unless ``--no-adapt`` is given.

Usage::

    python scripts/hmc_probe/batch_harness.py --manifest <manifest.json> --adapted <dir> \\
        --src <code dir containing aslprep/> --out <dir> --label c0 [--candidate c0] [--jobs 8]
"""

import argparse
import itertools
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from panel_config import ADAPTATIONS, EXCLUDED  # noqa: E402
from run_panel import IMAGE, SITE  # noqa: E402

RADIUS = 80.0


def rms(m, centre):
    """RMS displacement of a 4x4 RAS map over an 80 mm sphere about ``centre`` (rmsdiff)."""
    d = np.asarray(m) - np.eye(4)
    lin, t = d[:3, :3], d[:3, 3] + d[:3, :3] @ centre
    return float(np.sqrt(RADIUS**2 / 5 * np.trace(lin.T @ lin) + t @ t))


def motion_summary(run_dir):
    """Estimated motion from the harness transforms: framewise and from the series median."""
    import nibabel as nb
    from harness import load_itk_series

    h = load_itk_series(run_dir / 'hmc_xfms.txt')
    ref = nb.load(run_dir / 'aslref.nii.gz')
    centre = ref.affine @ np.r_[(np.array(ref.shape[:3]) - 1) / 2, 1]
    centre = centre[:3]
    indices = json.loads((run_dir / 'hmc_indices.json').read_text())
    context = pd.read_table(run_dir / 'hmc_aslcontext.tsv')['volume_type'].tolist()
    order = np.argsort(indices)
    fd = [rms(h[b] @ np.linalg.inv(h[a]), centre) for a, b in itertools.pairwise(order)]
    dev = [rms(m, centre) for m in h]
    pairs = []
    controls = [i for i, t in enumerate(context) if t == 'control']
    labels = [i for i, t in enumerate(context) if t == 'label']
    if controls and len(controls) == len(labels):
        pairs = [
            rms(h[c] @ np.linalg.inv(h[lab]), centre)
            for c, lab in zip(controls, labels, strict=True)
        ]
    return {
        'n_volumes': len(h),
        'fd_mean_mm': float(np.mean(fd)) if fd else 0.0,
        'fd_max_mm': float(np.max(fd)) if fd else 0.0,
        'offset_median_mm': float(np.median(dev)),
        'offset_max_mm': float(np.max(dev)),
        'pair_median_mm': float(np.median(pairs)) if pairs else None,
        'pair_max_mm': float(np.max(pairs)) if pairs else None,
    }


def run_one(row, args):
    name = row['dataset']
    root = Path(row['root'])
    if not args.no_adapt and name in ADAPTATIONS:
        root = args.adapted / name
    out = args.out / args.label / name / row['run']
    if (out / 'motion.json').exists() and not args.rerun:
        return row['run'], 'done'
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    cmd = ['docker', 'run', '--rm', '-u', f'{os.getuid()}:{os.getgid()}', '-e', 'HOME=/tmp']
    cmd += ['-e', 'FSLOUTPUTTYPE=NIFTI_GZ']
    for host, container, mode in (
        (HERE, '/probe', 'ro'),
        (args.src / 'aslprep', SITE, 'ro'),
        (root, '/data', 'ro'),
        (out, '/out', 'rw'),
    ):
        cmd += ['-v', f'{host}:{container}:{mode}']
    cmd += ['--entrypoint', '/app/.pixi/envs/aslprep/bin/python', IMAGE, '/probe/harness.py']
    cmd += ['--asl', f'/data/{row["asl_path"]}', '--candidate', args.candidate, '--out', '/out']
    m0 = row.get('m0_candidates') or []
    if m0 and row.get('m0type') == 'Separate':
        # Prefer the M0 with the run's phase-encoding direction (as the adapted copies pin it).
        direction = next((p for p in row['run'].split('_') if p.startswith('dir-')), None)
        same_dir = [m for m in m0 if direction and f'_{direction}_' in m]
        cmd += ['--m0scan', f'/data/{(same_dir or m0)[0]}']
    cmd += args.harness_args
    with (out / 'harness.log').open('w') as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)
    if proc.returncode:
        return row['run'], f'FAILED ({proc.returncode})'
    summary = {'dataset': name, 'run': row['run'], 'family': row['family'], **motion_summary(out)}
    (out / 'motion.json').write_text(json.dumps(summary, indent=2) + '\n')
    # Keep the transforms and the reference; the resampled series are regenerable.
    for f in ('corrected.nii.gz',):
        (out / f).unlink(missing_ok=True)
    return row['run'], 'ok'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--adapted', type=Path, required=True)
    parser.add_argument('--src', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--candidate', default='c0')
    parser.add_argument('--jobs', type=int, default=8)
    parser.add_argument('--no-adapt', action='store_true')
    parser.add_argument('--only', nargs='*')
    parser.add_argument('--rerun', action='store_true')
    parser.add_argument(
        '--harness-args',
        nargs=argparse.REMAINDER,
        default=[],
        help='Passed to harness.py (e.g. --opt engine=volreg --downstream d1); must come last.',
    )
    args = parser.parse_args(argv)
    args.out, args.src, args.adapted = (
        args.out.resolve(),
        args.src.resolve(),
        args.adapted.resolve(),
    )
    rows = [
        r
        for r in json.loads(args.manifest.read_text())
        if r['dataset'] not in EXCLUDED
        and not any(f.startswith('degenerate') for f in r['predicted_failures'])
        # e.g. a ds005752 series with 2 volumes under a 3-row dataset-level aslcontext
        and r.get('context_matches_volumes', True)
        and (not args.only or r['dataset'] in args.only)
    ]
    print(f'{len(rows)} runs, {args.jobs} at a time', flush=True)
    with ThreadPoolExecutor(args.jobs) as pool:
        for run, status in pool.map(lambda r: run_one(r, args), rows):
            print(f'{run}: {status}', flush=True)


if __name__ == '__main__':
    main()
