#!/usr/bin/env python3
"""Run ASLPrep on the real-data panel, one container per subject, several at a time.

The ASLPrep code is mounted from ``--src`` over the package in the image, so a pinned revision
(or a candidate) runs without rebuilding the image. Every run is fresh: its own output and
working directories, and no reuse of precomputed derivatives.

Usage::

    python scripts/hmc_probe/run_panel.py --manifest <manifest.json> --zoo <aslprep-zoo> \\
        --adapted <adapted dir> --src <code dir containing aslprep/> --out <dir> --label baseline
"""

import argparse
import base64
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from panel_config import (
    ADAPTATIONS,
    BIDS_FILTERS,
    EXCLUDED,
    EXTRA_ARGS,
    FULL_PIPELINE_EXCLUDED,
)

IMAGE = 'pennlinc/aslprep:unstable'
SITE = '/app/.pixi/envs/aslprep/lib/python3.12/site-packages/aslprep'
#: The FreeSurfer license used by ASLPrep's own tests (aslprep/tests/conftest.py).
LICENSE = (
    'bWF0dGhldy5jaWVzbGFrQHBzeWNoLnVjc2IuZWR1CjIwNzA2CipDZmVWZEg1VVQ4clkKRlNCWVouVWtlVElDdwo='
)


def _root(zoo, name):
    base = zoo / name
    return base / 'dset' if (base / 'dset' / 'dataset_description.json').exists() else base


def jobs_from(manifest, zoo, adapted, originals, no_adapt=False):
    """One job per (dataset, subject).

    ``originals`` runs the unadapted and excluded datasets. ``no_adapt`` takes every dataset
    from ``zoo`` as it is (for the extended panel, which ``extend_panel.py`` already adapted).
    """
    jobs = {}
    for row in manifest:
        name, sub = row['dataset'], row['subject']
        if any(f.startswith('degenerate affine') for f in row.get('predicted_failures', [])):
            print(f'skipping {name} sub-{sub}: {row["run"]} has a degenerate affine')
            continue
        if originals:
            if name not in ADAPTATIONS and name not in EXCLUDED:
                continue
            root = _root(zoo, name)
        elif name in FULL_PIPELINE_EXCLUDED:
            continue
        elif no_adapt:
            root = _root(zoo, name)
        else:
            if name in EXCLUDED:
                continue
            root = adapted / name if name in ADAPTATIONS else _root(zoo, name)
        jobs[(name, sub)] = root
    return jobs


def run_job(name, sub, root, args):
    tag = f'{name}_sub-{sub}'
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
    license_file = args.out / 'license.txt'
    filters = BIDS_FILTERS.get(name)
    extra = list(EXTRA_ARGS.get(name, []))
    # Copies made by adapt_panel.py or extend_panel.py carry probe notes and edits that the
    # image's legacy validator rejects; they are generated, so validation is skipped.
    if root != _root(args.zoo, name) or args.no_adapt:
        if '--skip-bids-validation' not in extra:
            extra.append('--skip-bids-validation')
    mounts = [
        (args.src / 'aslprep', SITE, 'ro'),
        (root, '/data', 'ro'),
        (out, '/out', 'rw'),
        (work, '/work', 'rw'),
        (license_file, '/license.txt', 'ro'),
    ]
    if filters:
        filter_file = work / 'bids_filters.json'
        filter_file.write_text(json.dumps(filters, indent=2))
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
        sub,
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
        *args.extra,
    ]
    start = time.time()
    with log.open('w') as fobj:
        fobj.write(' '.join(cmd) + '\n\n')
        fobj.flush()
        proc = subprocess.run(cmd, stdout=fobj, stderr=subprocess.STDOUT, check=False)
    record = {
        'dataset': name,
        'subject': sub,
        'data': str(root),
        'src': str(args.src),
        'label': args.label,
        'returncode': proc.returncode,
        'seconds': round(time.time() - start),
        'cmd': cmd,
    }
    done.write_text(json.dumps(record, indent=2) + '\n')
    # Disk space is limited: a successful run's working directory is not needed again.
    if proc.returncode == 0 and not args.keep_work:
        shutil.rmtree(work, ignore_errors=True)
    return tag, proc.returncode


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--zoo', type=Path, required=True)
    parser.add_argument('--adapted', type=Path, required=True)
    parser.add_argument('--src', type=Path, required=True, help='Directory containing aslprep/.')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--originals', action='store_true', help='Run unadapted originals.')
    parser.add_argument(
        '--no-adapt', action='store_true', help='Use every dataset in --zoo as it is.'
    )
    parser.add_argument('--only', nargs='*', help='Datasets to run.')
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--nprocs', type=int, default=6)
    parser.add_argument('--mem-mb', type=int, default=7000)
    parser.add_argument('--rerun', action='store_true')
    parser.add_argument(
        '--keep-work', action='store_true', help='Keep working directories of successful runs.'
    )
    parser.add_argument('extra', nargs='*', help='Extra ASLPrep arguments (after --).')
    args = parser.parse_args(argv)
    args.out, args.zoo = args.out.resolve(), args.zoo.resolve()
    (args.out / args.label).mkdir(parents=True, exist_ok=True)
    (args.out / 'license.txt').write_text(base64.b64decode(LICENSE).decode())

    manifest = json.loads(args.manifest.read_text())
    jobs = jobs_from(
        manifest, args.zoo.resolve(), args.adapted.resolve(), args.originals, args.no_adapt
    )
    if args.only:
        jobs = {k: v for k, v in jobs.items() if k[0] in args.only}
    print(f'{len(jobs)} jobs, {args.jobs} at a time', flush=True)
    with ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(run_job, n, s, r, args) for (n, s), r in sorted(jobs.items())]
        for fut in futures:
            tag, code = fut.result()
            print(f'{tag}: {"ok" if code == 0 else f"FAILED ({code})"}', flush=True)


if __name__ == '__main__':
    main()
