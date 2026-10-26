#!/usr/bin/env python3
"""Fit C6 on the host (PyTorch, GPU) for each fixture; the harness candidate ``c6`` reads the fits.

The ASLPrep image has no PyTorch, so the C6 fit runs here and writes one ``<asl stem>.npz`` per
run into ``--out``. ``batch_fixtures.py --mount <out>:/c6fit`` then runs the harness with
``--opt fit_dir=/c6fit``.

Usage::

    python scripts/hmc_probe/c6_fit.py --data-dir <dir> --out <fit dir> [--names ...]
        [--init zero|random --init-mm 1 --seed 0 --iters 1.0]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from c6_identifiability import tissue_mask  # noqa: E402
from candidates import c6_fit_file, c6_fit_inputs  # noqa: E402
from harness import resolve_inputs  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--names', nargs='*')
    parser.add_argument('--init', choices=('zero', 'random'), default='zero')
    parser.add_argument('--init-mm', type=float, default=1.0)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--iters', type=float, default=1.0, help='Scale on every stage.')
    parser.add_argument('--hires', type=int, default=2, help='1: the plain interpolation model.')
    parser.add_argument(
        '--support',
        choices=('image', 'truth'),
        default='image',
        help='The perfusion support: from the static image, or the true tissue (upper bound).',
    )
    parser.add_argument(
        '--p-bound', type=float, help='Bound |P| to this fraction of the static image.'
    )
    parser.add_argument('--rerun', action='store_true')
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    root = args.data_dir.resolve() / 'aslscan'
    fixtures = sorted(
        p
        for p in root.iterdir()
        if (p / 'manifest.json').exists() and (not args.names or p.name in args.names)
    )
    for fx in fixtures:
        asl = next((fx / 'sub-01' / 'perf').glob('*_asl.nii.gz'))
        out = c6_fit_file(args.out, asl)
        if out.exists() and not args.rerun:
            print(f'{fx.name}: done', flush=True)
            continue
        t0 = time.time()
        mask = tissue_mask(fx) if args.support == 'truth' else None
        res = c6_fit_inputs(
            resolve_inputs(str(asl)),
            args.init,
            args.init_mm,
            args.seed,
            args.iters,
            args.hires,
            mask,
            args.p_bound,
        )
        if res is None:
            np.savez(out, applied=False, meta=json.dumps({'reason': 'no pairs'}))
            print(f'{fx.name}: not applicable', flush=True)
            continue
        meta = {k: res[k] for k in ('init', 'init_mm', 'seed', 'hires', 'p_bound', 'w', 'history')}
        meta['seconds'] = round(time.time() - t0, 1)
        meta['start'] = None if res['start'] is None else np.asarray(res['start']).tolist()
        np.savez(
            out,
            applied=True,
            pulls=res['pulls'],
            theta=res['theta'],
            S=res['S'],
            P=res['P'],
            cl=np.asarray(res['cl']),
            meta=json.dumps(meta),
        )
        print(
            f'{fx.name}: {meta["seconds"]} s, final loss {res["history"][-1]["loss"]:.3g}',
            flush=True,
        )


if __name__ == '__main__':
    main()
