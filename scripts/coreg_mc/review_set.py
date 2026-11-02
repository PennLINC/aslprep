#!/usr/bin/env python3
"""Draw the held-out review set (plan: "The review set and applicability", frozen in Task 3).

One run per held-out subject, chosen at random within the subject from a recorded seed, before
any figure exists. Every held-out subject is kept, so every spiral and GRASE subject is in.
Subjects, not runs, are the review's sampling unit (the bootstrap is over subjects), so one run
per subject costs little evidence; it brings the review from up to 832 pairs (208 runs) to up to
212 (53 runs, two J-arms each against R1 and Q1), before the tie distance drops close pairs.

Writes ``review_set.tsv`` beside this script (the manifest's columns for the chosen runs).

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/review_set.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SEED = 20261102


def main():
    manifest = pd.read_table(HERE / 'manifest.tsv')
    held = manifest[manifest['side'] == 'held_out'].sort_values(['dataset', 'subject', 'run'])
    rng = np.random.default_rng(SEED)
    chosen = [g.iloc[rng.integers(len(g))] for _, g in held.groupby(['dataset', 'subject'])]
    out = pd.DataFrame(chosen).reset_index(drop=True)
    out.to_csv(HERE / 'review_set.tsv', sep='\t', index=False)
    summary = out.groupby('family').size().rename('runs').to_frame()
    summary['held_out_runs'] = held.groupby('family').size()
    print(summary.to_string())
    print(f'{len(out)} runs from {len(held)} held-out runs')
    return 0


if __name__ == '__main__':
    sys.exit(main())
