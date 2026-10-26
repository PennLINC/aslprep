#!/usr/bin/env python3
"""Fix the development and held-out lists of the coregistration probe (plan Task 0).

Plan: "Development and held-out data". Three axes:

- dataset: whole datasets held out (``HELD_OUT_DATASETS``), including every dataset of the held-out
  spiral scanner group;
- subject: in each extension dataset, about half the subjects held out, drawn by ``SEED``;
- subjects that the HMC probe's extended panel also holds (``extended/split.json``) keep the HMC
  probe's assignment, so a subject is never development in one probe and held out in the other.

The panel's own exemplar subjects (the zoo) are development unless their dataset is held out.

Usage::

    python split.py --panel-manifest <hmc-probe panel/manifest.tsv> \\
        --hmc-split <hmc-probe extended/split.json> --subjects extension_subjects.tsv --out split.json
"""

import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path

SEED = 20261022
#: tarricone (multi-delay GRASE) and ds008719 (2D EPI) as in the HMC probe's plan; the spiral
#: scanner group with serial 66044 (ds006131, ds000240, ds008579).
HELD_OUT_DATASETS = {
    'tarricone': 'dataset (multi-delay GRASE)',
    'ds008719': 'dataset (2D EPI)',
    'ds006131': 'scanner group 66044 (stack-of-spirals)',
    'ds000240': 'scanner group 66044 (stack-of-spirals)',
    'ds008579': 'scanner group 66044 (stack-of-spirals)',
}


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--panel-manifest', type=Path, required=True)
    p.add_argument('--hmc-split', type=Path, required=True)
    p.add_argument('--subjects', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()

    rng = random.Random(SEED)
    hmc = json.loads(a.hmc_split.read_text())
    hmc_side = {
        (ds, sub): ('held_out' if sub in v['held_out'] else 'development')
        for ds, v in hmc.items()
        for sub in v['subjects']
    }

    rows = []
    with a.panel_manifest.open() as f:
        panel = {(r['dataset'], f'sub-{r["subject"]}') for r in csv.DictReader(f, delimiter='\t')}
    for ds, sub in sorted(panel):
        side, why = (
            ('held_out', HELD_OUT_DATASETS[ds])
            if ds in HELD_OUT_DATASETS
            else ('development', 'panel exemplar')
        )
        rows.append({'dataset': ds, 'subject': sub, 'source': 'panel', 'side': side, 'axis': why})

    with a.subjects.open() as f:
        ext = defaultdict(list)
        for r in csv.DictReader(f, delimiter='\t'):
            ext[r['dataset']].append(r['subject'])
    for ds, subs in sorted(ext.items()):
        fixed = {s: hmc_side[(ds, s)] for s in subs if (ds, s) in hmc_side}
        free = sorted(s for s in subs if s not in fixed)
        rng.shuffle(free)
        # Half of the dataset's extension subjects held out, counting the fixed ones.
        n_held = round(len(subs) / 2) - sum(v == 'held_out' for v in fixed.values())
        n_held = max(0, min(n_held, len(free)))
        for s in subs:
            if s in fixed:
                side, why = fixed[s], 'subject (HMC probe assignment)'
            else:
                side = 'held_out' if s in free[:n_held] else 'development'
                why = 'subject'
            row = {'dataset': ds, 'subject': s, 'source': 'extension', 'side': side, 'axis': why}
            rows.append(row)

    record = {'seed': SEED, 'held_out_datasets': HELD_OUT_DATASETS, 'subjects': rows}
    a.out.write_text(json.dumps(record, indent=1) + '\n')
    count = defaultdict(lambda: [0, 0])
    for r in rows:
        count[r['dataset']][r['side'] == 'held_out'] += 1
    for ds, (dev, held) in sorted(count.items()):
        print(f'{ds:28s} development {dev:3d}  held out {held:3d}')
    print('total', sum(c[0] for c in count.values()), sum(c[1] for c in count.values()))


if __name__ == '__main__':
    main()
