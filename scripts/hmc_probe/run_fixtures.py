#!/usr/bin/env python3
"""Run full ASLPrep on simulated fixtures, as ``run_panel.py`` does for the real panel.

Each fixture is a one-subject BIDS dataset made by ASLPrep's fixture machinery or by
``fixtures.py``. Its anatomical derivatives are passed with ``--derivatives`` and its M0
divisor with ``--m0_scale``, as ``aslprep/tests/aslscan_cli.py`` does.

Usage::

    python scripts/hmc_probe/run_fixtures.py --fixtures <dir>/aslscan --names rc_pcasl ... \\
        --src <code dir containing aslprep/> --out <dir> --label baseline [--sloppy]
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from run_panel import IMAGE, LICENSE, SITE


def run_one(name, args):
    import base64

    fixture = (args.fixtures / name).resolve()
    truth = json.loads((fixture / 'sub-01' / 'perf' / 'ground-truth' / 'truth.json').read_text())
    label = args.label + ('-sloppy' if args.sloppy else '')
    out = args.out / label / name
    work = args.out / 'work' / label / name
    done = out / 'DONE.json'
    if done.exists() and not args.rerun:
        return name, json.loads(done.read_text())['returncode']
    for d in (out, work):
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True)
    (work / 'home').mkdir()
    license_file = args.out / 'license.txt'
    license_file.write_text(base64.b64decode(LICENSE).decode())
    cmd = ['docker', 'run', '--rm', '-u', f'{os.getuid()}:{os.getgid()}', '-e', 'HOME=/work/home']
    for host, container, mode in (
        (args.src / 'aslprep', SITE, 'ro'),
        (fixture, '/data', 'ro'),
        (out, '/out', 'rw'),
        (work, '/work', 'rw'),
        (license_file, '/license.txt', 'ro'),
    ):
        cmd += ['-v', f'{host}:{container}:{mode}']
    hmc = (args.hmc_derivs / name).resolve() if args.hmc_derivs else None
    if hmc is not None:
        if not hmc.exists():
            return name, 'no injected HMC derivatives'
        cmd += ['-v', f'{hmc}:/hmcderiv:ro']
    cmd += [IMAGE, '/data', '/out', 'participant', '--participant-label', '01', '-w', '/work']
    cmd += ['--fs-license-file', '/license.txt', '--fs-no-reconall', '--notrack']
    # The image's legacy validator rejects the probe's BIDS fieldmaps (*_fieldmap with
    # *_magnitude); the fixtures are generated, so validation is skipped.
    cmd += ['--skip-bids-validation']
    cmd += ['--output-spaces', 'asl', 'T1w', '--smooth_kernel', '5']
    cmd += ['--nprocs', str(args.nprocs), '--omp-nthreads', str(args.nprocs)]
    cmd += ['--mem', str(args.mem_mb)]
    if hmc is not None:
        # One dataset holds the anatomical derivatives and the injected HMC products, and it must
        # be collected last: each derivatives dataset's functional cache replaces the
        # ``transforms`` entry, and ASLPrep appends its atlas datasets after the user's (cli/
        # parser.py), so by default they wipe any precomputed transform. Naming them first keeps
        # them out of the way (the parser only adds them when absent).
        cmd += [
            '--derivatives',
            'aslprepatlases=/home/aslprep/.cache/aslprep/XCPDAtlases',
            'aslprep4s=/home/aslprep/.cache/aslprep/AtlasPack',
            'anat=/hmcderiv',
        ]
    elif (fixture / 'derivatives' / 'anat').exists():
        cmd += ['--derivatives', 'anat=/data/derivatives/anat']
    divisor = float(truth.get('m0_divisor', 1.0))
    if divisor != 1:
        cmd += ['--m0_scale', str(divisor)]
    if args.sloppy:
        cmd.append('--sloppy')
    start = time.time()
    with (args.out / label / f'{name}.log').open('w') as fobj:
        fobj.write(' '.join(cmd) + '\n\n')
        fobj.flush()
        proc = subprocess.run(cmd, stdout=fobj, stderr=subprocess.STDOUT, check=False)
    done.write_text(
        json.dumps(
            {
                'fixture': str(fixture),
                'src': str(args.src),
                'label': label,
                'returncode': proc.returncode,
                'seconds': round(time.time() - start),
                'cmd': cmd,
            },
            indent=2,
        )
        + '\n'
    )
    if proc.returncode == 0 and not args.keep_work:
        shutil.rmtree(work, ignore_errors=True)
    return name, proc.returncode


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--fixtures', type=Path, required=True, help='<data dir>/aslscan')
    parser.add_argument('--names', nargs='+', required=True)
    parser.add_argument('--src', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--sloppy', action='store_true')
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--nprocs', type=int, default=6)
    parser.add_argument('--mem-mb', type=int, default=7000)
    parser.add_argument('--rerun', action='store_true')
    parser.add_argument('--keep-work', action='store_true')
    parser.add_argument(
        '--hmc-derivs',
        type=Path,
        help='Task 7: <dir>/<fixture> holds injected HMC derivatives (inject_hmc.py).',
    )
    args = parser.parse_args(argv)
    args.out = args.out.resolve()
    args.src = args.src.resolve()
    (args.out / (args.label + ('-sloppy' if args.sloppy else ''))).mkdir(
        parents=True, exist_ok=True
    )
    with ThreadPoolExecutor(args.jobs) as pool:
        for name, code in pool.map(lambda n: run_one(n, args), args.names):
            print(f'{name}: {"ok" if code == 0 else f"FAILED ({code})"}', flush=True)


if __name__ == '__main__':
    main()
