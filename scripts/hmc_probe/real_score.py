#!/usr/bin/env python3
"""Real-data scores for harness outputs (plan: "Real data"), no truth needed.

For each run in ``<runs>/<label>/<dataset>/<run>/`` (``batch_harness.py`` output):

Decision metrics (direction fixed in advance):

- ``split_half_r``: correlation, inside the brain, of the mean delta-M of two halves of the
  pairs, the halves balanced within each delay (alternate pairs of each delay). Runs with fewer
  than two pairs at any delay are not eligible (``None``).
- ``negative_fraction``: fraction of brain voxels whose mean delta-M is negative.

Diagnostics (no assumed good direction):

- ``edge_leakage_r``: correlation of absolute mean delta-M with the gradient magnitude of the
  reference, inside the brain;
- ``alternation_power``: for the six rigid parameters of the transforms, the share of their
  variance at the label/control frequency (alternate volumes), averaged.

The brain mask is the reference thresholded at 20 % of its 98th percentile, largest connected
component, so it is the same for every candidate only if the references agree; the fixed
evaluation support of the plan replaces it once full-pipeline masks exist.

Usage::

    python scripts/hmc_probe/real_score.py --runs <dir> --labels c0 c1 ... --manifest <json> ...
"""

import argparse
import json
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
from scipy import ndimage as ndi


def brain_mask(ref):
    data = np.asanyarray(ref.dataobj, np.float32)
    thr = 0.2 * np.percentile(data[data > 0], 98) if (data > 0).any() else 0
    mask = data > thr
    lab, n = ndi.label(mask)
    if n > 1:
        sizes = ndi.sum(mask, lab, range(1, n + 1))
        mask = lab == (np.argmax(sizes) + 1)
    return ndi.binary_fill_holes(mask)


def rigid_params(m):
    """Translation (mm) and ZYX Euler angles (degrees) of a 4x4 rigid matrix."""
    r = m[:3, :3]
    ry = np.arcsin(-np.clip(r[2, 0], -1, 1))
    rx = np.arctan2(r[2, 1], r[2, 2])
    rz = np.arctan2(r[1, 0], r[0, 0])
    return np.r_[m[:3, 3], np.degrees([rx, ry, rz])]


def score_run(run_dir, plds_by_volume):
    from harness import load_itk_series

    out = {}
    ctx = pd.read_table(run_dir / 'hmc_aslcontext.tsv')['volume_type'].tolist()
    indices = json.loads((run_dir / 'hmc_indices.json').read_text())
    ref = nb.load(run_dir / 'aslref.nii.gz')
    mask = brain_mask(ref)
    dm_file = run_dir / 'deltam.nii.gz'
    if dm_file.exists():
        dm = np.asanyarray(nb.load(dm_file).dataobj, np.float32)
        dm = dm[..., None] if dm.ndim == 3 else dm
        mean = dm.mean(-1)
        out['negative_fraction'] = float((mean[mask] < 0).mean())
        controls = [indices[k] for k, t in enumerate(ctx) if t == 'control']
        pair_pld = [plds_by_volume[i] for i in controls]
        halves = ([], [])
        eligible = True
        for p in sorted(set(pair_pld)):
            ks = [k for k, q in enumerate(pair_pld) if q == p]
            if len(ks) < 2:
                eligible = False
                break
            halves[0].extend(ks[0::2])
            halves[1].extend(ks[1::2])
        if eligible and halves[0] and halves[1]:
            a, b = dm[..., halves[0]].mean(-1)[mask], dm[..., halves[1]].mean(-1)[mask]
            out['split_half_r'] = float(np.corrcoef(a, b)[0, 1])
        else:
            out['split_half_r'] = None
        grad = np.linalg.norm(np.gradient(np.asanyarray(ref.dataobj, np.float32)), axis=0)
        out['edge_leakage_r'] = float(np.corrcoef(np.abs(mean[mask]), grad[mask])[0, 1])
    h = load_itk_series(run_dir / 'hmc_xfms.txt')
    order = np.argsort(indices)
    params = np.array([rigid_params(h[k]) for k in order])
    if len(params) >= 4:
        centred = params - params.mean(0)
        alt = np.array([(-1) ** k for k in range(len(params))])
        power = []
        for j in range(6):
            total = (centred[:, j] ** 2).sum()
            if total > 0:
                coef = centred[:, j] @ alt / len(alt)
                power.append(coef**2 * len(alt) / total)
        out['alternation_power'] = float(np.mean(power)) if power else None
    return out


def plds(manifest_row):
    root = Path(manifest_row['root'])
    side = root / manifest_row['asl_path'].replace('.nii.gz', '.json')
    meta = json.loads(side.read_text()) if side.exists() else {}
    pld = np.atleast_1d(np.asarray(meta.get('PostLabelingDelay', 0.0), float))
    n = manifest_row['n_volumes']
    return list(np.repeat(pld, n)) if pld.size == 1 else list(pld)


def main(argv=None):
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--runs', type=Path, required=True)
    parser.add_argument('--labels', nargs='+', required=True)
    parser.add_argument('--manifest', type=Path, nargs='+', required=True)
    parser.add_argument('--tsv', type=Path)
    args = parser.parse_args(argv)
    rows = {}
    for m in args.manifest:
        for r in json.loads(m.read_text()):
            rows[(r['dataset'], r['run'])] = r
    records = []
    for label in args.labels:
        for run_dir in sorted((args.runs / label).glob('*/*/hmc_xfms.txt')):
            run_dir = run_dir.parent
            key = (run_dir.parent.name, run_dir.name)
            if key not in rows:
                continue
            try:
                s = score_run(run_dir, plds(rows[key]))
            except Exception as exc:  # noqa: BLE001 - one bad run must not stop the table
                s = {'error': repr(exc)[:200]}
            records.append({'label': label, 'dataset': key[0], 'run': key[1],
                            'family': rows[key]['family'], **s})  # fmt: skip
    df = pd.DataFrame(records)
    cols = ['split_half_r', 'negative_fraction', 'edge_leakage_r', 'alternation_power']
    with pd.option_context('display.width', 200):
        print(df.groupby(['family', 'label'])[cols].median().round(4))
    if args.tsv:
        df.to_csv(args.tsv, sep='\t', index=False)


if __name__ == '__main__':
    main()
