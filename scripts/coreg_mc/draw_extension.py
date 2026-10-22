#!/usr/bin/env python3
"""Draw the repeat-scan subjects of the panel extension from the zoo's DataLad clones.

Plan: ``docs/plans/2026-10-22-joint-multicontrast-coregistration.md``, "Panel extension".
Nothing is downloaded: subjects, sessions and runs come from the clone's file tree, and sizes
from the git-annex keys. Each dataset has a file filter narrower than the zoo's
``02_download_subjects.sh``: only the files its curation keeps (``FILTERS``), from sessions with
ASL. The selected files are written to ``--files`` for the download.

A subject is eligible if it has ASL in at least two sessions, or at least two ASL runs in one
session that differ by more than the echo entity, has the anatomical image the zoo's curation
uses as its T1w, and is not already in the panel. Per dataset, subjects are drawn at random with
one seed, after the preference the plan states.

Usage::

    python draw_extension.py --cache ~/.cache/aslprep-zoo --out extension_subjects.tsv
"""

import argparse
import csv
import os
import random
import re
from pathlib import Path

SEED = 20261022

#: dataset -> (number to draw, glob of the image used as T1w, subjects already in the panel)
DATASETS = {
    'ds007116': (20, '*_rec-normrefaced*_T1w.nii.gz', {'sub-20721'}),
    'ds003823': (20, '*_acq-wholebrain_T1w.nii*', {'sub-5023'}),
    'ds004146': (20, '*_UNIT1_denoised.nii.gz', {'sub-0337'}),
    'ds008020': (None, '*_UNIT1.nii*', {'sub-01', 'sub-19'}),
    'ds007522': (10, '*_T1w.nii*', {'sub-HML0202', 'sub-HML0001'}),
    'ds005752': (None, '*_T1w.nii*', {'sub-ON13986'}),
}
#: dataset -> regular expression on ``<datatype>/<file name>`` of the files to download, from
#: sessions with ASL: what each curator in the zoo's ``03_curate_openneuro.py`` keeps.
FILTERS = {
    'ds007116': r'^(anat/.*_rec-normrefaced.*_T1w|anat/.*_rec-normdefaced.*_T2w|perf/.*)\.',
    'ds003823': r'^(anat/.*_acq-wholebrain_T1w|perf/.*)\.',
    # The separate M0 is not used (M0Type is Included; the HMC probe's adaptation removes it).
    'ds004146': r'^(anat/.*_UNIT1_denoised|anat/.*_T2w|perf/.*_(asl|aslcontext))\.',
    'ds008020': r'^(anat/.*_UNIT1|perf/.*)\.',
    'ds007522': r'^(anat/.*_(T1w|T2w|FLAIR)|perf/.*|fmap/.*)\.',
    'ds005752': r'^(anat/.*_(T1w|T2w)|perf/.*)\.',
}


def annex_size(path):
    if path.is_symlink():
        m = re.search(r'-s(\d+)--', os.readlink(path))
        return int(m.group(1)) if m else 0
    return path.stat().st_size if path.exists() else 0


def survey(root, t1w_glob, keep):
    keep = re.compile(keep)
    for sub in sorted(p for p in root.glob('sub-*') if p.is_dir()):
        sessions = sorted(p for p in sub.glob('ses-*') if p.is_dir()) or [sub]
        asl = {}
        for ses in sessions:
            names = sorted(p.name for p in (ses / 'perf').glob('*_asl.nii*'))
            if names:
                asl[ses] = names
        if not asl:
            continue
        runs = {s: {re.sub(r'_echo-[^_]+', '', n) for n in v} for s, v in asl.items()}
        multi_ses = len(asl) >= 2
        multi_run = any(len(v) >= 2 for v in runs.values())
        has_t1w = any(list((s / 'anat').glob(t1w_glob)) for s in asl)
        files = [
            p for ses in asl for dt in ('anat', 'perf', 'fmap') if (ses / dt).exists()
            for p in sorted((ses / dt).iterdir())
            if not p.is_dir() and keep.search(f'{dt}/{p.name}')
        ]
        nbytes = sum(annex_size(p) for p in files)
        variant = any('VARIANT' in n for v in asl.values() for n in v)
        yield {
            'subject': sub.name,
            'sessions': ','.join(s.name for s in asl if s != sub),
            'n_sessions': len(asl),
            'max_runs_per_session': max(len(v) for v in runs.values()),
            'multi_session': multi_ses,
            'multi_run': multi_run,
            'variant': variant,
            'has_t1w': has_t1w,
            'mb': round(nbytes / 1e6, 1),
            '_files': [str(p.relative_to(root)) for p in files],
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--files', type=Path, required=True)
    args = parser.parse_args()

    rng = random.Random(SEED)
    rows = []
    for ds, (n, t1w_glob, existing) in DATASETS.items():
        cands = [
            r for r in survey(args.cache / ds, t1w_glob, FILTERS[ds])
            if (r['multi_session'] or r['multi_run']) and r['has_t1w'] and r['subject'] not in existing
        ]
        if ds == 'ds007522':
            # The four-run pCASL protocol only (one subject mixes Philips and Siemens series).
            cands = [r for r in cands if r['max_runs_per_session'] >= 4]
        if ds == 'ds007116':
            # Subjects that keep one protocol across sessions first, then the variants.
            same = [r for r in cands if not r['variant']]
            rest = [r for r in cands if r['variant']]
            rng.shuffle(same)
            rng.shuffle(rest)
            ordered = same + rest
        else:
            ordered = sorted(cands, key=lambda r: r['subject'])
            rng.shuffle(ordered)
        chosen = ordered if n is None else ordered[:n]
        print(f'{ds}: {len(cands)} eligible, {len(chosen)} drawn, '
              f'{sum(r["mb"] for r in chosen) / 1e3:.2f} GB')
        rows += [{'dataset': ds, **r} for r in sorted(chosen, key=lambda r: r['subject'])]

    with args.files.open('w') as f:
        for r in rows:
            for rel in r['_files']:
                f.write(f'{r["dataset"]}\t{rel}\n')
    rows = [{k: v for k, v in r.items() if k != '_files'} for r in rows]
    with args.out.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    print(f'total: {len(rows)} subjects, {sum(r["mb"] for r in rows) / 1e3:.2f} GB')


if __name__ == '__main__':
    main()
