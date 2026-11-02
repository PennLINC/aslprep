#!/usr/bin/env python3
"""Run times of production's coregistration and the frozen ANTs arms (plan: metric 5, gate 5).

On a seeded sample of development real runs (``--per-family`` per family), each step runs
single-threaded (``ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1``, ``OMP_NUM_THREADS=1``) in the
pinned image, one run at a time, with the work directory mounted read-only and every output
written to a tmpfs (``/out``) inside the container:

- ``mri_coreg``: production's initialization, which every arm starts from;
- ``bbr``: FLIRT-BBR from it (R1's chain with the default sign);
- ``S1``, ``S2``, ``Q1``, ``J1``: the frozen ANTs arms (``run_ants.ARMS``) from the same start.

Writes one row per run and step (``seconds``, wall time inside the container).

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/time_arms.py --dev-real <dir> --license <file> --out <tsv>
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'hmc_probe'))

import run_ants  # noqa: E402
from run_panel import IMAGE  # noqa: E402

SEED = 20261102
STEPS = ('mri_coreg', 'bbr', *run_ants.ARMS)


def step_lines():
    """Shell lines that time each step and print ``<step> <seconds>``."""
    lines = [
        'export FS_LICENSE=/lic/license.txt ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1',
        'export OMP_NUM_THREADS=1',
        'S=$FSLDIR/etc/flirtsch/bbr.sch',
        't() { local s=$EPOCHREALTIME; "${@:2}" > /out/$1.log 2>&1; local rc=$?; '
        'awk -v n="$1" -v a="$s" -v b="$EPOCHREALTIME" -v r="$rc" '
        '\'BEGIN { printf "%s %.3f %d\\n", n, b - a, r }\'; }',
        't mri_coreg mri_coreg --mov /w/conventional.nii.gz --ref /w/t1w_brain.nii.gz '
        '--reg /out/init.lta --dof 6 --sep 4 --ftol 0.0001 --linmintol 0.01',
        'lta_convert --inlta /out/init.lta --outfsl /out/init.mat > /dev/null',
        't bbr flirt -in /w/conventional.nii.gz -ref /w/t1w_brain.nii.gz -cost bbr -dof 6 '
        '-basescale 1 -wmseg /w/wm.nii.gz -init /out/init.mat -omat /out/bbr.mat -schedule $S',
    ]
    for arm in run_ants.ARMS:
        cmd = run_ants.ants_command(arm, '/w/init_itk.txt', f'/out/ants_{arm}_')
        lines.append(f't {arm} {cmd}')

    return lines


def time_run(work, license_file):
    script = '\n'.join(step_lines()) + '\n'
    cmd = ['docker', 'run', '--rm', '-u', f'{os.getuid()}:{os.getgid()}', '-e', 'HOME=/tmp']
    cmd += ['-v', f'{work}:/w:ro', '-v', f'{license_file.parent}:/lic:ro', '--tmpfs', '/out']
    cmd += ['--entrypoint', 'bash', IMAGE, '-c', script]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        step, seconds, status = line.split()
        rows.append({'step': step, 'seconds': float(seconds), 'status': int(status)})

    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--dev-real', type=Path, required=True)
    p.add_argument('--license', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--per-family', type=int, default=2)
    args = p.parse_args(argv)
    manifest = pd.read_table(HERE / 'manifest.tsv')
    dev = manifest[manifest['side'] == 'development']
    rng = np.random.default_rng(SEED)
    sample = pd.concat(
        [g.iloc[rng.permutation(len(g))[: args.per_family]] for _, g in dev.groupby('family')],
    )
    rows = []
    for r in sample.itertuples(index=False):
        work = args.dev_real / 'work' / f'{r.dataset}__{r.run}'
        for row in time_run(work, args.license):
            rows.append({'dataset': r.dataset, 'run': r.run, 'family': r.family, **row})
            print(r.dataset, r.run, row, flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    print(df.pivot_table(index='step', values='seconds', aggfunc=['median', 'max']).round(1))


if __name__ == '__main__':
    main()
