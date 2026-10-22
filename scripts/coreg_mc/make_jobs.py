#!/usr/bin/env python3
"""Write the ``run_aslprep.py`` job files for the panel sample and the panel extension.

Plan: Task 0. The sample is one run per family from the HMC probe's panel, rerun on the pinned
revision to decide whether the HMC probe's baseline runs can be reused. The extension is every
subject of ``extension_subjects.tsv``, curated into ``--ext``.

Usage::

    python make_jobs.py --zoo <aslprep-zoo> --adapted <hmc-probe panel/adapted> --ext <panel-ext> \\
        --subjects extension_subjects.tsv --out <dir>
"""

import argparse
import csv
from pathlib import Path

#: Per-dataset arguments, as the HMC probe's ``panel_config.EXTRA_ARGS``.
EXTRA = {
    'ds007522': '--skip-bids-validation',
    'ds008020': '--skull-strip-t1w skip',
}
#: One run per family: (dataset, subject, family). ds004146 runs from the HMC probe's adapted copy.
SAMPLE = [
    ('ds007116', '20721', 'stack-of-spirals'),
    ('ds003823', '5023', '3D pCASL, background-suppressed'),
    ('ds007522', 'HML0202', '3D GRASE, two pairs'),
    ('ds004146', '0337', '2D EPI, background suppression'),
    ('ds004856', '3131', '2D EPI, no background suppression'),
    ('ds008020', '19', '3D EPI at 7 T'),
    ('ds005752', 'ON13986', 'GE product output'),
]
ADAPTED = {'ds004146'}


def write(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['tag', 'root', 'subject', 'extra'], delimiter='\t',
                           lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print(f'{path}: {len(rows)} jobs')


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--zoo', type=Path, required=True)
    p.add_argument('--adapted', type=Path, required=True)
    p.add_argument('--ext', type=Path, required=True)
    p.add_argument('--subjects', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()

    sample = []
    for ds, sub, _family in SAMPLE:
        root = a.adapted / ds if ds in ADAPTED else a.zoo / ds
        sample.append({'tag': f'{ds}_sub-{sub}', 'root': str(root), 'subject': sub,
                       'extra': EXTRA.get(ds, '')})
    write(a.out / 'jobs_sample.tsv', sample)

    with a.subjects.open() as f:
        subjects = list(csv.DictReader(f, delimiter='\t'))
    ext = []
    for r in subjects:
        ds, sub = r['dataset'], r['subject'].removeprefix('sub-')
        # Curated and checked with the current BIDS validator; the image's legacy one is skipped.
        extra = ' '.join(x for x in ('--skip-bids-validation', EXTRA.get(ds, '')) if x)
        extra = ' '.join(dict.fromkeys(extra.split(' ')))
        ext.append({'tag': f'{ds}_sub-{sub}', 'root': str(a.ext / ds), 'subject': sub,
                    'extra': extra})
    write(a.out / 'jobs_extension.tsv', ext)


if __name__ == '__main__':
    main()
