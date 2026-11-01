#!/usr/bin/env python3
"""Score a completed review batch against its key (plan: "Blinded review figures").

Reads the review page's choices (one JSON file per figure id, as ``ArtifactData list`` saves the
``choices`` collection with ``out_dir``), the batch's ``key.json`` and, for a development batch,
``dev_review.py``'s ``selection.tsv``. Every side the reviewer named is turned into the arm it
showed, so nothing here depends on which side an arm was drawn on.

- known-answer pairs: detected if the reviewer prefers the unmoved solution; a tie or the moved
  one is a miss. Detection rate per size.
- repeats: whether the repeat names the same arm (or a tie) as its original.
- comparisons: the arm preferred, or a tie, with the pair's distance and family; and which arms
  the reviewer flagged as misaligned.

Writes ``<out>/choices.tsv`` (one row per figure, unblinded) and prints the summaries.

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/score_review.py --choices <dir> --batch <review dir> --out <dir>
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

#: Known-answer detection at 2 mm below this means the figures are redesigned (plan, gate 4).
MIN_DETECTION_2MM = 0.8


def preferred(choice, key):
    """The arm the reviewer preferred, or ``'tie'``."""
    if choice == 'tie':
        return 'tie'

    return key[f'arm_{choice}']


def load(choices_dir, batch_dir):
    key = json.loads((batch_dir / 'key.json').read_text())['figures']
    rows = []
    for f in sorted(choices_dir.glob('*.json')):
        c = json.loads(f.read_text())
        k = key[f.stem]
        rows.append(
            {
                'id': f.stem,
                'position': c['position'],
                'kind': k['kind'],
                'run': k['run'],
                'arm_A': k['arm_A'],
                'arm_B': k['arm_B'],
                'repeat_of': k['repeat_of'],
                'size_mm': k['size_mm'],
                'choice': c['choice'],
                'preferred': preferred(c['choice'], k),
                'failed': ';'.join(k[f'arm_{s}'] for s in 'AB' if c[f'failed_{s}']),
                'note': c.get('note', ''),
            },
        )
    missing = sorted(set(key) - {r['id'] for r in rows})
    return pd.DataFrame(rows).sort_values('position'), missing


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--choices', type=Path, required=True)
    p.add_argument('--batch', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args(argv)
    df, missing = load(args.choices, args.batch)
    if missing:
        print(f'UNANSWERED: {", ".join(missing)}')
        return 1

    selection = args.batch / 'selection.tsv'
    if selection.exists():
        sel = pd.read_table(selection)
        sel = sel[sel['chosen']]
        sel = sel.assign(
            run=sel['dataset'] + '/' + sel['run'],
            pair=[frozenset(p) for p in zip(sel['arm_a'], sel['arm_b'], strict=True)],
        )
        df['pair'] = [frozenset(p) for p in zip(df['arm_A'], df['arm_B'], strict=True)]
        cols = ['run', 'pair', 'comparison', 'family', 'distance_mm', 'separation_mm']
        df = df.merge(sel[cols], on=['run', 'pair'], how='left').drop(columns='pair')

    known = df[df['kind'] == 'known'].copy()
    known['detected'] = known['preferred'] == known['arm_A'].where(
        ~known['arm_A'].str.contains(' moved '),
        known['arm_B'],
    )
    known['outcome'] = 'detected'
    known.loc[known['preferred'] == 'tie', 'outcome'] = 'tie'
    known.loc[~known['detected'] & (known['preferred'] != 'tie'), 'outcome'] = 'wrong'
    print('Known-answer pairs:')
    print(pd.crosstab(known['size_mm'], known['outcome'], margins=True).to_string())
    rate = known.groupby('size_mm')['detected'].mean()
    print(rate.to_string())
    ok = rate.get(2.0, 0.0) >= MIN_DETECTION_2MM
    print(f'2 mm detection {rate.get(2.0, 0.0):.0%}: {"PASS" if ok else "REDESIGN"}')

    by_id = df.set_index('id')
    repeats = df[df['kind'] == 'repeat']
    print('\nRepeats:')
    for r in repeats.itertuples(index=False):
        first = by_id.loc[r.repeat_of, 'preferred']
        same = 'agree' if first == r.preferred else 'DISAGREE'
        print(f'  {r.run}: first {first!r}, repeat {r.preferred!r}: {same}')

    comp = df[df['kind'] == 'comparison']
    print('\nComparisons:')
    show = ['comparison', 'family', 'run', 'distance_mm', 'preferred', 'failed', 'note']
    with pd.option_context('display.width', 250, 'display.max_colwidth', 70):
        print(comp.sort_values(['comparison', 'family'])[show].to_string(index=False))
        print(pd.crosstab([comp['comparison'], comp['family']], comp['preferred']).to_string())

    args.out.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out / 'choices.tsv', sep='\t', index=False)
    return 0


if __name__ == '__main__':
    sys.exit(main())
