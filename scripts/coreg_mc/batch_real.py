#!/usr/bin/env python3
"""Run ``run_real.py`` over every subject of one side of the split, several processes at a time.

Subjects are dealt round-robin into ``--groups`` groups (each a ``run_real`` process with its own
output prefix and log); pairs never cross subjects, so the groups share one work directory
without overlap. Panel subjects read ``<runs>/panel``, extension subjects ``<runs>/extension``.

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/batch_real.py --side development --runs /home/tsalo/coreg-mc/runs \\
        --adapted <dir> --license <file> --work <dir> --out <dir> [--groups 8] [--arms ...]
"""

import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--manifest', type=Path, default=HERE / 'manifest.tsv')
    p.add_argument('--side', required=True, choices=('development', 'held_out'))
    p.add_argument('--runs', type=Path, required=True)
    p.add_argument('--adapted', type=Path, required=True)
    p.add_argument('--license', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--groups', type=int, default=8)
    p.add_argument('--arms', nargs='*', default=[])
    p.add_argument('--rescore', action='store_true')
    args = p.parse_args(argv)
    manifest = pd.read_table(args.manifest)
    rows = manifest[manifest['side'] == args.side]
    subjects = rows[['source', 'dataset', 'subject']].drop_duplicates().to_dict('records')
    args.out.mkdir(parents=True, exist_ok=True)
    jobs = []
    for source in sorted({s['source'] for s in subjects}):
        mine = [s for s in subjects if s['source'] == source]
        for i in range(args.groups):
            group = mine[i :: args.groups]
            if group:
                jobs.append((f'{source}-{i}', args.runs / source, group))

    def run(job):
        name, runs, group = job
        cmd = [
            sys.executable,
            str(HERE / 'run_real.py'),
            '--manifest',
            str(args.manifest),
            '--runs',
            str(runs),
            '--adapted',
            str(args.adapted),
            '--license',
            str(args.license),
            '--work',
            str(args.work),
            '--out',
            str(args.out / name),
            '--arms',
            *args.arms,
            '--subjects',
            *(f'{s["dataset"]}/{s["subject"]}' for s in group),
        ]
        if args.rescore:
            cmd.append('--rescore')

        with (args.out / f'{name}.log').open('w') as log:
            proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)
        return name, len(group), proc.returncode

    print(f'{len(subjects)} subjects in {len(jobs)} groups', flush=True)
    with ThreadPoolExecutor(args.groups) as pool:
        for name, n, code in pool.map(run, jobs):
            print(f'{name}: {n} subjects, exit {code}', flush=True)


if __name__ == '__main__':
    main()
