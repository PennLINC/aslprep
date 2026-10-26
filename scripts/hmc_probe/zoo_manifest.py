#!/usr/bin/env python3
"""Build the Task 0 manifest of the real-data panel (the aslprep-zoo collection).

For every ASL run, record its acquisition family, its volumes, the separate M0 scans and
fieldmaps associated with it, how each separate M0's grid relates to the ASL grid, and the
processing path and failures that ASLPrep's code predicts for it. The predictions mirror the
checks in ``aslprep/utils/bids.py`` (``collect_run_data``), ``aslprep/utils/asl.py``
(``infer_m0tr``), ``aslprep/interfaces/cbf.py`` (``ExtractDeltaM``) and
``aslprep/workflows/asl/base.py`` (the short-series GE path); the baseline run confirms them.

The panel is read in place and never modified.

Usage::

    python scripts/hmc_probe/zoo_manifest.py --zoo <aslprep-zoo> --out <dir>
"""

import argparse
import json
from collections import Counter
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd

#: Acquisition family of each dataset (keyed by dataset, or by dataset and subject).
FAMILIES = {
    'ds008719': '2D EPI, no BS',
    'morys_bad_brainmask': '2D EPI, no BS',
    'ds004856': 'Philips pCASL, no M0, no BS',
    'ds004146': '2D EPI, BS',
    'ds004873': '2D EPI, BS',
    '2d_mb_pcasl': '2D multiband, multi-delay',
    'ds000254': '2D multiband, GE',
    'ds007908': '3D GRASE, single delay',
    'ds006148': '3D GRASE, single delay',
    'ds003823': '3D GRASE, single delay',
    'kashyap': '3D GRASE, single delay',
    'ds007522': '3D GRASE, single delay',
    ('ds007522', 'HML0202'): '3D GRASE, single delay',
    'jsr-aslprep-debug': '3D GRASE, single delay',
    'lsa-aslprep-debug': '3D GRASE, multi-delay',
    'tarricone': '3D GRASE, multi-delay',
    'grotzinger_menopausal_asl': '3D GRASE, multi-delay',
    'ds000234': '3D stack-of-spirals',
    'ds000235': '3D stack-of-spirals',
    'ds000236': '3D stack-of-spirals',
    'ds000240': '3D stack-of-spirals',
    'ds007116': '3D stack-of-spirals',
    'ds008579': '3D stack-of-spirals',
    'ds006131': '3D stack-of-spirals',
    'shanmugan_hup6': '3D stack-of-spirals',
    'ds008020': '3D EPI, 7 T',
    'ds005752': 'GE product output',
}

#: ASLPrep drops these volume types before processing (``n/a``) or never models them.
DROPPED_TYPES = ('n/a', 'noRF')


def _dataset_roots(zoo):
    """Every BIDS root in the panel: ``<zoo>/<name>`` or ``<zoo>/<name>/dset``."""
    roots = []
    for child in sorted(zoo.iterdir()):
        if not child.is_dir() or child.name == 'code':
            continue
        root = child / 'dset' if (child / 'dset' / 'dataset_description.json').exists() else child
        if (root / 'dataset_description.json').exists():
            roots.append((child.name, root))
    return roots


def _entities(path):
    """BIDS entities of a file name, as a dict (``sub``, ``ses``, ``acq``...)."""
    stem = path.name.split('.')[0]
    parts = stem.split('_')
    ents = dict(p.split('-', 1) for p in parts[:-1] if '-' in p)
    ents['suffix'] = parts[-1]
    return ents


def _sidecar(path, root):
    """Merged JSON metadata, following BIDS inheritance from the root down."""
    ents = _entities(path)
    stem_suffix = ents['suffix']
    meta = {}
    candidates = []
    rel = path.relative_to(root)
    dirs = [root] + [root.joinpath(*rel.parts[: i + 1]) for i in range(len(rel.parts) - 1)]
    for d in dirs:
        for j in sorted(d.glob(f'*{stem_suffix}.json')):
            jents = _entities(j)
            if all(ents.get(k) == v for k, v in jents.items() if k != 'suffix'):
                candidates.append(j)
    for j in candidates:  # shallowest first, so deeper files override
        meta.update(json.loads(j.read_text()))
    return meta


def _aslcontext(path, root):
    """The nearest ``aslcontext.tsv`` whose entities are a subset of the run's."""
    ents = _entities(path)
    best = None
    for tsv in sorted(path.parent.glob('*_aslcontext.tsv')) + sorted(
        root.glob('*_aslcontext.tsv')
    ):
        tents = _entities(tsv)
        if all(ents.get(k) == v for k, v in tents.items() if k != 'suffix'):
            if best is None or len(tents) > len(_entities(best)):
                best = tsv
    return best


def _intended_for(meta, root, sub):
    """Absolute paths that a sidecar's ``IntendedFor`` names."""
    targets = meta.get('IntendedFor', [])
    targets = targets if isinstance(targets, list) else [targets]
    out = []
    for t in targets:
        if t.startswith('bids::'):
            out.append((root / t[len('bids::') :]).resolve())
        else:
            out.append((root / f'sub-{sub}' / t).resolve())
    return out


def _grid_relation(m0_img, asl_img):
    """How an M0 grid relates to an ASL grid.

    ``same``: same shape and affine (the #703 fix skips resampling).
    ``same-physical``: different storage (axis order or flips) of the same voxel centres.
    ``different``: a different physical sampling.
    """
    for img in (m0_img, asl_img):
        if abs(np.linalg.det(img.affine[:3, :3])) < 1e-6:
            return 'degenerate-affine'
    if m0_img.shape[:3] == asl_img.shape[:3] and np.allclose(
        m0_img.affine, asl_img.affine, rtol=0, atol=1e-4
    ):
        return 'same'
    m0_can = nb.as_closest_canonical(nb.Nifti1Image(np.zeros(m0_img.shape[:3]), m0_img.affine))
    asl_can = nb.as_closest_canonical(nb.Nifti1Image(np.zeros(asl_img.shape[:3]), asl_img.affine))
    if m0_can.shape == asl_can.shape and np.allclose(
        m0_can.affine, asl_can.affine, rtol=0, atol=1e-3
    ):
        return 'same-physical'
    return 'different'


def _pairs_by_delay(context, meta):
    """Number of control/label pairs at each post-labeling delay."""
    pld = np.atleast_1d(np.asarray(meta.get('PostLabelingDelay', np.nan), dtype=float))
    if pld.size == 1:
        pld = np.repeat(pld, len(context))
    controls = Counter(float(pld[i]) for i, t in enumerate(context) if t == 'control')
    labels = Counter(float(pld[i]) for i, t in enumerate(context) if t == 'label')
    return {d: min(controls[d], labels[d]) for d in sorted(set(controls) | set(labels))}


def _predict(row, meta, context):
    """Failures and path that ASLPrep's code predicts for a run."""
    failures = []
    if 'degenerate-affine' in row['m0_grid'] or 0 in row['zooms_mm']:
        failures.append('degenerate affine (a zero voxel size): the image cannot be resampled')
    m0type = meta.get('M0Type')
    if m0type == 'Separate' and not row['m0_candidates']:
        failures.append('M0Type Separate but no M0 scan found (collect_run_data)')
    if m0type != 'Separate' and row['m0_candidates']:
        failures.append(f'M0Type {m0type} but a separate M0 is associated (collect_run_data)')
    if (
        m0type == 'Absent'
        and 'control' in context
        and 'cbf' not in context
        and meta.get('BackgroundSuppression')
    ):
        failures.append('background-suppressed controls cannot calibrate (infer_m0tr)')
    if 'control' in context and context.count('control') != context.count('label'):
        failures.append(
            f'{context.count("control")} controls != {context.count("label")} labels (ExtractDeltaM)'
        )
    tr_prep = np.atleast_1d(np.asarray(meta.get('RepetitionTimePreparation', np.nan), dtype=float))
    if m0type == 'Included':
        m0_idx = [i for i, t in enumerate(context) if t == 'm0scan']
        m0tr = tr_prep[m0_idx] if tr_prep.size > 1 else tr_prep
        if np.any(m0tr == 0):
            failures.append('M0 RepetitionTimePreparation is 0: non-finite calibration')
    if meta.get('MRAcquisitionType') == '2D' and 'SliceTiming' not in meta:
        failures.append('2D without SliceTiming (fails BIDS validation)')
    if row['echo'] is not None and row.get('n_echoes', 1) > 1:
        failures.append('multi-echo: outputs drop the echo entity and overwrite each other')
    unknown = sorted(set(context) - {'control', 'label', 'm0scan', 'deltam', 'cbf', 'n/a', 'noRF'})
    if unknown:
        failures.append(f'unknown volume types {unknown}')

    path = 'GE (short series)' if row['n_volumes'] <= 5 else 'standard'
    if 'control' in context and 'label' in context:
        target = 'control'
    elif 'deltam' in context:
        target = 'deltam'
    elif 'cbf' in context:
        target = 'cbf'
    else:
        target = None
    # SelectHighestContrastVolumes, with prioritize_m0 on the GE path.
    ge = path.startswith('GE')
    if 'm0scan' in context and ge:
        ref = 'm0scan'
    elif row['m0_candidates'] and ge:
        ref = 'separate_m0scan'
    elif 'cbf' in context:
        ref = 'cbf'
    elif 'deltam' in context:
        ref = 'deltam'
    elif 'm0scan' in context:
        ref = 'm0scan'
    elif row['m0_candidates']:
        ref = 'separate_m0scan'
    else:
        ref = 'control'
    notes = []
    # ReduceASLFiles keeps the processing target's types and any in-series m0scan for HMC.
    corrected = {'control': ('control', 'label', 'm0scan'), 'deltam': ('deltam', 'm0scan')}
    if ref == 'separate_m0scan' or ref not in corrected.get(target, (target, 'm0scan')):
        notes.append(f'reference ({ref}) is not in the motion-corrected series')
    if set(context) & set(DROPPED_TYPES):
        notes.append(f'drops {sorted(set(context) & set(DROPPED_TYPES))} before processing')
    if m0type == 'Separate' and 'm0scan' in context:
        notes.append('in-series m0scan with M0Type Separate')
    if len(row['m0_candidates']) > 1:
        notes.append(f'{len(row["m0_candidates"])} M0 candidates; the first is used')
    return failures, path, target, ref, notes


def build(zoo):
    rows = []
    for name, root in _dataset_roots(zoo):
        fmap_targets = {}
        for fj in sorted(root.glob('sub-*/**/fmap/*.json')):
            sub = _entities(fj)['sub']
            for target in _intended_for(json.loads(fj.read_text()), root, sub):
                fmap_targets.setdefault(target, []).append(fj.name)
        m0_targets = {}
        for mj in sorted(root.glob('sub-*/**/perf/*_m0scan.json')):
            sub = _entities(mj)['sub']
            for target in _intended_for(json.loads(mj.read_text()), root, sub):
                m0_targets.setdefault(target, []).append(mj.with_name(mj.name[:-5] + '.nii.gz'))

        for asl in sorted(root.glob('sub-*/**/perf/*_asl.nii.gz')):
            ents = _entities(asl)
            meta = _sidecar(asl, root)
            ctx_file = _aslcontext(asl, root)
            context = (
                pd.read_table(ctx_file, keep_default_na=False)['volume_type'].tolist()
                if ctx_file
                else []
            )
            img = nb.load(asl)
            n_volumes = img.shape[3] if img.ndim == 4 else 1
            m0s = m0_targets.get(asl.resolve(), [])
            row = {
                'dataset': name,
                'run': asl.name.replace('_asl.nii.gz', ''),
                'asl_path': str(asl.relative_to(root)),
                'root': str(root),
                'subject': ents['sub'],
                'session': ents.get('ses'),
                'echo': ents.get('echo'),
                'n_echoes': len(
                    list(asl.parent.glob(asl.name.replace(f'echo-{ents.get("echo")}', 'echo-*')))
                )
                if ents.get('echo')
                else 1,
                'family': FAMILIES.get((name, ents['sub']), FAMILIES.get(name, 'unassigned')),
                'manufacturer': meta.get('Manufacturer'),
                'scanner_serial': meta.get('DeviceSerialNumber'),
                'labeling': meta.get('ArterialSpinLabelingType'),
                'acquisition': meta.get('MRAcquisitionType'),
                'background_suppression': meta.get('BackgroundSuppression'),
                'm0type': meta.get('M0Type'),
                'n_volumes': n_volumes,
                'shape': list(img.shape[:3]),
                'zooms_mm': [round(float(z), 3) for z in img.header.get_zooms()[:3]],
                'volume_types': dict(Counter(context)),
                'first_volume': context[0] if context else None,
                'pairs_by_delay': {str(k): v for k, v in _pairs_by_delay(context, meta).items()},
                'm0_candidates': [str(m.relative_to(root)) for m in m0s],
                'm0_grid': [_grid_relation(nb.load(m), img) for m in m0s],
                'fieldmaps': fmap_targets.get(asl.resolve(), []),
                'aslcontext': str(ctx_file.relative_to(root)) if ctx_file else None,
                'context_matches_volumes': len(context) == n_volumes,
            }
            failures, path, target, ref, notes = _predict(row, meta, context)
            row.update(
                {
                    'predicted_failures': failures,
                    'processing_path': path,
                    'processing_target': target,
                    'reference_type': ref,
                    'notes': notes,
                }
            )
            rows.append(row)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--zoo', type=Path, required=True, help='The aslprep-zoo directory.')
    parser.add_argument('--out', type=Path, required=True, help='Output directory.')
    args = parser.parse_args(argv)
    rows = build(args.zoo.resolve())
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'manifest.json').write_text(json.dumps(rows, indent=1) + '\n')
    table = pd.DataFrame(rows)
    for col in table.columns:
        table[col] = table[col].map(lambda v: json.dumps(v) if isinstance(v, list | dict) else v)
    table.to_csv(args.out / 'manifest.tsv', sep='\t', index=False)
    print(f'{len(rows)} runs in {table["dataset"].nunique()} datasets -> {args.out}')


if __name__ == '__main__':
    main()
