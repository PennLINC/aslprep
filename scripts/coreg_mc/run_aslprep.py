#!/usr/bin/env python3
"""Run ASLPrep on panel subjects, one container per subject, several at a time.

Adapted from the HMC probe's ``scripts/hmc_probe/run_panel.py``. The pinned image runs the
revision it was built from; ``--src`` (a directory containing ``aslprep/``) optionally mounts
another revision over the package. Every run is fresh: its own output and working directories
and no precomputed derivatives. A successful run's working directory is deleted.

Jobs come from a TSV with columns ``tag``, ``root`` (BIDS directory), ``subject`` (without
``sub-``), ``extra`` (space-separated ASLPrep arguments, may be empty) and ``filters`` (a BIDS
filter JSON, may be empty).

Usage (in WSL)::

    python run_aslprep.py --jobs-file jobs.tsv --src <dir> --out <dir> --label <label> \\
        [--only tag ...] [--parallel 4]
"""

import argparse
import base64
import csv
import json
import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

#: pennlinc/aslprep:unstable as of 2026-10-26, built from 072ee0f (typewise_offset HMC, #773).
IMAGE = 'pennlinc/aslprep@sha256:7b34cc5066181c4f23d42602376569613ee59e9c68de0c9899076ffa9c8450f6'
SITE = '/app/.pixi/envs/aslprep/lib/python3.12/site-packages/aslprep'
#: The FreeSurfer license used by ASLPrep's own tests (aslprep/tests/conftest.py).
LICENSE = (
    'bWF0dGhldy5jaWVzbGFrQHBzeWNoLnVjc2IuZWR1CjIwNzA2CipDZmVWZEg1VVQ4clkKRlNCWVouVWtlVElDdwo='
)


def run_job(job, args):
    tag = job['tag']
    out = args.out / args.label / tag
    work = args.out / 'work' / args.label / tag
    log = args.out / args.label / f'{tag}.log'
    done = out / 'DONE.json'
    if done.exists() and not args.rerun:
        return tag, json.loads(done.read_text())['returncode']
    for d in (out, work):
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True)
    (work / 'home').mkdir()
    mounts = [
        (Path(job['root']), '/data', 'ro'),
        (out, '/out', 'rw'),
        (work, '/work', 'rw'),
        (args.out / 'license.txt', '/license.txt', 'ro'),
    ]
    if args.src:
        mounts.append((args.src / 'aslprep', SITE, 'ro'))

    extra = job['extra'].split()
    if job.get('filters'):
        (work / 'bids_filters.json').write_text(job['filters'] + '\n')
        extra += ['--bids-filter-file', '/work/bids_filters.json']

    cmd = ['docker', 'run', '--rm', '-u', f'{os.getuid()}:{os.getgid()}', '-e', 'HOME=/work/home']
    for host, container, mode in mounts:
        cmd += ['-v', f'{host}:{container}:{mode}']
    cmd += [
        IMAGE,
        '/data',
        '/out',
        'participant',
        '--participant-label',
        job['subject'],
        '-w',
        '/work',
        '--fs-license-file',
        '/license.txt',
        '--fs-no-reconall',
        '--output-spaces',
        'asl',
        'T1w',
        '--nprocs',
        str(args.nprocs),
        '--omp-nthreads',
        str(args.nprocs),
        '--mem',
        str(args.mem_mb),
        '--notrack',
        *extra,
    ]
    start = time.time()
    with log.open('w') as fobj:
        fobj.write(' '.join(cmd) + '\n\n')
        fobj.flush()
        proc = subprocess.run(cmd, stdout=fobj, stderr=subprocess.STDOUT, check=False)
    record = {
        **job,
        'src': str(args.src) if args.src else None,
        'image': IMAGE,
        'label': args.label,
        'returncode': proc.returncode,
        'seconds': round(time.time() - start),
        'cmd': cmd,
    }
    done.write_text(json.dumps(record, indent=2) + '\n')
    if proc.returncode == 0:
        shutil.rmtree(work, ignore_errors=True)
    return tag, proc.returncode


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--jobs-file', type=Path, required=True)
    parser.add_argument('--src', type=Path, help='Directory containing aslprep/ to mount.')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--only', nargs='*')
    parser.add_argument('--parallel', type=int, default=4)
    parser.add_argument('--nprocs', type=int, default=6)
    parser.add_argument('--mem-mb', type=int, default=6500)
    parser.add_argument('--rerun', action='store_true')
    args = parser.parse_args(argv)
    args.out = args.out.resolve()
    (args.out / args.label).mkdir(parents=True, exist_ok=True)
    (args.out / 'license.txt').write_text(base64.b64decode(LICENSE).decode())

    with args.jobs_file.open() as f:
        jobs = list(csv.DictReader(f, delimiter='\t'))
    if args.only:
        jobs = [j for j in jobs if j['tag'] in args.only]
    print(f'{len(jobs)} jobs, {args.parallel} at a time', flush=True)
    with ThreadPoolExecutor(args.parallel) as pool:
        futures = [pool.submit(run_job, j, args) for j in jobs]
        for fut in futures:
            tag, code = fut.result()
            print(f'{tag}: {"ok" if code == 0 else f"FAILED ({code})"}', flush=True)


if __name__ == '__main__':
    main()
