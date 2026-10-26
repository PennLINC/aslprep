#!/usr/bin/env python3
"""D3 eligibility statistics for every real-data harness run (input to Task 5b's calibration).

For each run of ``<runs>/<dataset>/<run>/`` with control-label pairs, the pairwise delta-M
volumes the harness kept (``deltam.nii.gz``) are the builder's inputs; delay selection,
split halves and the edge statistic follow ``perfusion.py``. The mask is nilearn's EPI mask of
the run's reference. Runs without pairs (derived volumes only) are listed with their data type
and no statistic, since the harness does not keep their corrected series.

Usage::

    python perfusion_screen.py --runs <harness-runs/c0> --manifest <manifest.json> [...] \\
        --adapted <panel/adapted> --out <perfusion_c0.tsv>
"""

import argparse
import json
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
from panel_config import ADAPTATIONS
from perfusion import EDGE_LIMIT, MIN_PAIRS, eligibility, epi_mask


def _sidecar(record, adapted):
    """The run's sidecar: the adapted copy's where one exists for this run, else the original."""
    asl = Path(record['root']) / record['asl_path']
    if record['dataset'] in ADAPTATIONS:
        copy = adapted / record['dataset'] / record['asl_path']
        asl = copy if copy.exists() else asl
    return json.loads(asl.with_name(asl.name.replace('.nii.gz', '.json')).read_text())


def run_record(run_dir, record, adapted):
    context = pd.read_table(run_dir / 'hmc_aslcontext.tsv')['volume_type'].tolist()
    indices = json.loads((run_dir / 'hmc_indices.json').read_text())
    out = {'dataset': record['dataset'], 'run': record['run'], 'family': record['family']}
    dm_file = run_dir / 'deltam.nii.gz'
    if not dm_file.exists():
        kinds = sorted(set(context) & {'deltam', 'cbf'})
        n = sum(t in kinds for t in context)
        dtype = 'derived_single' if n == 1 else 'derived_series'
        return out | {
            'data_type': dtype,
            'n_inputs': n,
            'eligible': False,
            'reason': 'single derived volume' if n == 1 else 'corrected series not kept',
        }
    meta = _sidecar(record, adapted)
    pld = np.atleast_1d(np.asarray(meta.get('PostLabelingDelay', 0.0), float))
    pld = np.repeat(pld, max(indices) + 1) if pld.size == 1 else pld
    controls = [indices[k] for k, t in enumerate(context) if t == 'control']
    pair_plds = np.array([pld[i] for i in controls], float)
    dm = nb.load(dm_file)
    data = np.asanyarray(dm.dataobj, np.float32)
    data = data[..., None] if data.ndim == 3 else data
    keep = pair_plds >= np.median(pair_plds)
    mask = epi_mask(run_dir / 'aslref.nii.gz')
    rec = eligibility(data[..., keep], 'pairs', pair_plds[keep], mask, dm.header.get_zooms()[:3])
    return out | rec | {'n_pairs_total': len(controls)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--runs', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, nargs='+', required=True)
    parser.add_argument('--adapted', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    records = {}
    for m in args.manifest:
        for r in json.loads(m.read_text()):
            records[(r['dataset'], r['run'])] = r
    rows = []
    for run_dir in sorted(p for p in args.runs.glob('*/*') if (p / 'harness.json').exists()):
        key = (run_dir.parent.name, run_dir.name)
        if key not in records:
            print(f'{key}: not in the manifests')
            continue
        try:
            rows.append(run_record(run_dir, records[key], args.adapted))
        except Exception as exc:  # noqa: BLE001 - report and continue over the panel
            print(f'{key}: {type(exc).__name__}: {exc}')
    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    print(
        f'{len(df)} runs -> {args.out} (placeholders: {MIN_PAIRS} pairs, edge_r >= {EDGE_LIMIT})'
    )
    print(df.groupby(['family', 'data_type'])[['edge_r', 'value_r']].median().round(3))
    print(df.groupby('family').eligible.mean().round(2))


if __name__ == '__main__':
    main()
