#!/usr/bin/env python3
"""Task 7: the D3 coregistration sources on every injected full run (metric 10 for D3).

For each ``<runs>/<label>/<fixture>`` (except the oracle), ``coreg_experiment.run_fixture`` builds
the perfusion-image levels from that run's corrected series and registers them, and the run's
conventional reference, by fMRIPrep's FSL path (both BBR signs), scored against the truth.
D3 keeps every other consumer as the run had it, so this is D3 on top of each arm's policy.

Usage::

    python task7_coreg.py --runs <dir> --derivs <dir> --fixtures <dir> [<dir> ...] \\
        --work <dir> --license <license.txt> --out <tsv> [--jobs 4]
"""

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from coreg_experiment import run_fixture  # noqa: E402
from task7_score import hmc_file  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--runs', type=Path, required=True)
    parser.add_argument('--derivs', type=Path, required=True)
    parser.add_argument('--fixtures', type=Path, nargs='+', required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--license', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args(argv)
    fixture_dirs = {p.name: p for root in args.fixtures for p in root.iterdir() if p.is_dir()}
    todo = []
    for done in sorted(args.runs.glob('*/*/DONE.json')):
        label, name = done.parent.parent.name, done.parent.name
        if label == 'oracle' or name not in fixture_dirs:
            continue
        if json.loads(done.read_text())['returncode'] == 0:
            todo.append((label, name))

    def one(item):
        label, name = item
        work = args.work / label / name
        cached = work / 'rows.json'
        if cached.exists():
            return json.loads(cached.read_text())
        try:
            rows = run_fixture(
                fixture_dirs[name],
                args.runs / label / name,
                work,
                args.license,
                hmc_file=hmc_file(args.derivs, label, name),
            )
        except Exception as exc:  # noqa: BLE001 - record and continue
            rows = [{'fixture': name, 'problems': repr(exc)}]
        rows = [{'label': label, **r} for r in rows]
        cached.write_text(json.dumps(rows, default=float) + '\n')
        print(f'{label}/{name}: {len(rows)} rows', flush=True)
        return rows

    with ThreadPoolExecutor(args.jobs) as pool:
        rows = [r for rs in pool.map(one, todo) for r in rs]
    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    print(df.pivot_table(index=['label'], columns=['source', 'sign'], values='rms_mm',
                         aggfunc='median').round(3).T.to_string())


if __name__ == '__main__':
    main()
