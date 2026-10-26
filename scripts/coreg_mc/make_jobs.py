#!/usr/bin/env python3
"""Write the ``run_aslprep.py`` job files for the panel and the panel extension.

Plan: Task 0. The panel is every ``source == 'panel'`` subject of ``split.json`` (the zoo's
exemplars), run as the HMC probe ran them: adapted copies, exclusions, extra arguments and BIDS
filters from ``scripts/hmc_probe/panel_config.py``. The extension is every subject of
``extension_subjects.tsv``, curated into ``--ext``.

Usage::

    python make_jobs.py --zoo <aslprep-zoo> --adapted <adapted copies> --ext <panel-ext> \\
        --split split.json --subjects extension_subjects.tsv --out <dir>
"""

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'hmc_probe'))
from panel_config import (
    ADAPTATIONS,
    BIDS_FILTERS,
    EXCLUDED,
    EXTRA_ARGS,
    FULL_PIPELINE_EXCLUDED,
)

#: Extension datasets' arguments (the curated copies passed the current BIDS validator, so the
#: image's legacy one is skipped).
EXTENSION_EXTRA = {'ds008020': ['--skull-strip-t1w', 'skip']}
FIELDS = ['tag', 'root', 'subject', 'extra', 'filters']


def zoo_root(zoo, name):
    base = zoo / name
    return base / 'dset' if (base / 'dset' / 'dataset_description.json').exists() else base


def write(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print(f'{path}: {len(rows)} jobs')


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--zoo', type=Path, required=True)
    p.add_argument('--adapted', type=Path, required=True)
    p.add_argument('--ext', type=Path, required=True)
    p.add_argument('--split', type=Path, required=True)
    p.add_argument('--subjects', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()

    panel = []
    split = json.loads(a.split.read_text())['subjects']
    for s in split:
        ds, sub = s['dataset'], s['subject'].removeprefix('sub-')
        if s['source'] != 'panel':
            continue

        if ds in EXCLUDED or ds in FULL_PIPELINE_EXCLUDED:
            print(f'skipping {ds} sub-{sub}: {(EXCLUDED | FULL_PIPELINE_EXCLUDED)[ds]}')
            continue

        root = a.adapted / ds if ds in ADAPTATIONS else zoo_root(a.zoo, ds)
        extra = list(EXTRA_ARGS.get(ds, []))
        if ds in ADAPTATIONS and '--skip-bids-validation' not in extra:
            extra.append('--skip-bids-validation')  # the adapted copies carry probe notes
        panel.append(
            {
                'tag': f'{ds}_sub-{sub}',
                'root': str(root),
                'subject': sub,
                'extra': ' '.join(extra),
                'filters': json.dumps(BIDS_FILTERS[ds]) if ds in BIDS_FILTERS else '',
            },
        )
    write(a.out / 'jobs_panel.tsv', panel)

    with a.subjects.open() as f:
        subjects = list(csv.DictReader(f, delimiter='\t'))
    ext = []
    for r in subjects:
        ds, sub = r['dataset'], r['subject'].removeprefix('sub-')
        ext.append(
            {
                'tag': f'{ds}_sub-{sub}',
                'root': str(a.ext / ds),
                'subject': sub,
                'extra': ' '.join(['--skip-bids-validation', *EXTENSION_EXTRA.get(ds, [])]),
                'filters': '',
            },
        )
    write(a.out / 'jobs_extension.tsv', ext)


if __name__ == '__main__':
    main()
