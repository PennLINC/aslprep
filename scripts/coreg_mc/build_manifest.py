#!/usr/bin/env python3
"""Task 0 manifest of the coregistration probe: one row per ASL run, with what applies to it.

Plan: Task 0 and "The review set and applicability". Each run's acquisition is read from its
sidecar and ``aslcontext.tsv``; the family comes from ``FAMILIES``; the split side from
``split.json``. Applicability is what the inputs allow before any ASLPrep run: a comparable
cross-run partner in the same session, a second session, a T2w, a separate or included M0. What
depends on ASLPrep outputs (perfusion-image eligibility, reference type) is added later.

Usage::

    python build_manifest.py --split split.json --zoo <aslprep-zoo> --adapted <panel/adapted> \\
        --ext <panel-ext> --out manifest.tsv
"""

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

#: Acquisition family per dataset, as in the plan's real-data panel table.
FAMILIES = {
    **dict.fromkeys(['ds000234', 'ds000235', 'ds000236', 'ds000240', 'ds007116', 'ds008579',
                     'ds006131', 'shanmugan_hup6'], 'stack-of-spirals'),
    **dict.fromkeys(['ds007908', 'ds006148', 'ds003823', 'kashyap', 'ds007522',
                     'jsr-aslprep-debug', 'lsa-aslprep-debug', 'tarricone',
                     'grotzinger_menopausal_asl'], '3D GRASE / 3D pCASL'),
    **dict.fromkeys(['ds008719', 'morys_bad_brainmask', 'ds004856', 'ds004146', 'ds004873',
                     '2d_mb_pcasl', 'ds000254'], '2D EPI'),
    'ds008020': '3D EPI 7T',
    'ds005752': 'GE product',
}
ADAPTED = {'ds004146', 'ds006148', 'ds000254'}


def root_for(row, zoo, adapted, ext):
    ds = row['dataset']
    if row['source'] == 'extension':
        return ext / ds
    if ds in ADAPTED and (adapted / ds).exists():
        return adapted / ds
    base = zoo / ds
    return base / 'dset' if (base / 'dset' / 'dataset_description.json').exists() else base


def comparable_key(name):
    """Runs of one session that can be paired: same phase-encoding direction, echo 1 only."""
    pe = re.search(r'_dir-([A-Za-z]+)_', name)
    return pe.group(1) if pe else ''


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--split', type=Path, required=True)
    p.add_argument('--zoo', type=Path, required=True)
    p.add_argument('--adapted', type=Path, required=True)
    p.add_argument('--ext', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()

    split = json.loads(a.split.read_text())['subjects']
    rows = []
    for s in split:
        root = root_for(s, a.zoo, a.adapted, a.ext)
        sub_dir = root / s['subject']
        if not sub_dir.exists():
            continue
        t2w = bool(list(sub_dir.rglob('*_T2w.nii*')))
        asl_files = sorted(p for p in sub_dir.rglob('perf/*_asl.nii*') if '_echo-2' not in p.name)
        sessions = defaultdict(list)
        for f in asl_files:
            ses = re.search(r'_(ses-[^_]+)_', f.name)
            sessions[ses.group(1) if ses else ''].append(f)
        for ses, files in sessions.items():
            for f in files:
                meta = json.loads(Path(re.sub(r'\.nii(\.gz)?$', '.json', str(f))).read_text())
                ctx_file = Path(re.sub(r'_asl\.nii(\.gz)?$', '_aslcontext.tsv', str(f)))
                ctx = ctx_file.read_text().split()[1:] if ctx_file.exists() else []
                partners = [g for g in files if g != f and comparable_key(g.name) == comparable_key(f.name)]
                fmap = any(
                    f.name in json.dumps(json.loads(j.read_text()).get('IntendedFor', ''))
                    for j in (f.parent.parent / 'fmap').glob('*.json')
                ) if (f.parent.parent / 'fmap').exists() else False
                pld = meta.get('PostLabelingDelay')
                rows.append({
                    'dataset': s['dataset'],
                    'subject': s['subject'],
                    'session': ses,
                    'run': f.name.split('_asl.')[0],
                    'family': FAMILIES.get(s['dataset'], 'unknown'),
                    'side': s['side'],
                    'axis': s['axis'],
                    'source': s['source'],
                    'manufacturer': meta.get('Manufacturer', ''),
                    'labeling': meta.get('ArterialSpinLabelingType', ''),
                    'readout': meta.get('MRAcquisitionType', ''),
                    'background_suppression': meta.get('BackgroundSuppression', ''),
                    'm0type': meta.get('M0Type', ''),
                    'n_delays': len(set(pld)) if isinstance(pld, list) else 1,
                    'volume_types': ','.join(sorted(set(ctx))),
                    'n_pairs': ctx.count('label'),
                    't2w': t2w,
                    'fieldmap': fmap,
                    'n_sessions': len(sessions),
                    'cross_run': bool(partners),
                    'cross_session': len(sessions) >= 2,
                    'root': str(root),
                })
    with a.out.open('w', newline='') as fobj:
        w = csv.DictWriter(fobj, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

    summary = defaultdict(lambda: defaultdict(set))
    for r in rows:
        key = (r['family'], r['side'])
        summary[key]['runs'].add((r['dataset'], r['subject'], r['session'], r['run']))
        summary[key]['subjects'].add((r['dataset'], r['subject']))
        if r['cross_run'] or r['cross_session']:
            summary[key]['repeat_subjects'].add((r['dataset'], r['subject']))
        if r['t2w']:
            summary[key]['t2w_subjects'].add((r['dataset'], r['subject']))
    print(f'{len(rows)} runs')
    for (fam, side), v in sorted(summary.items()):
        print(f'{fam:22s} {side:12s} runs {len(v["runs"]):4d}  subjects {len(v["subjects"]):3d}  '
              f'repeat subjects {len(v["repeat_subjects"]):3d}  T2w subjects {len(v["t2w_subjects"]):3d}')


if __name__ == '__main__':
    main()
