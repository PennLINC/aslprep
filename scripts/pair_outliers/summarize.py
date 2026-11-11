"""Summarize the simulation and real-data results for the report.

Every number in the report comes from these functions, which read the JSON files written by
``evaluate.py`` and ``real_eval.py``.
"""

import json
from collections import defaultdict
from itertools import pairwise
from pathlib import Path

import numpy as np

#: Methods in display order, with their labels.
METHOD_LABELS = {
    'none': 'No exclusion',
    'zfrac': 'Per-delay voxel z (model-free)',
    'zfrac_loo': 'Per-delay voxel z, leave-one-out (model-free)',
    'empirical': 'Same-delay median (κ = 0)',
    'hybrid': 'Hybrid (κ = 2), proposed',
    'model': 'Model only (κ = ∞)',
    'hybrid_raw': 'Hybrid, no smoothed energy',
    'hybrid_noamp': 'Hybrid, no amplitude test',
    'oracle': 'Oracle (truly corrupted pairs)',
}

OUTCOMES = {
    'dm_rmse_pct': 'Mean ΔM error (% of GM signal)',
    'region_cbf_rmse': 'Regional CBF RMSE (mL/100 g/min)',
    'region_att_rmse': 'Regional ATT RMSE (s)',
    'voxel_cbf_mae': 'Voxel CBF median abs. error',
    'voxel_att_mae': 'Voxel ATT median abs. error (s)',
}


def load(directory):
    return [json.loads(p.read_text()) for p in sorted(Path(directory).glob('*.json'))]


def exclusion_rates(results, methods):
    """Per method: share of clean pairs excluded, and share of corrupted pairs excluded.

    Returns ``{method: {'clean_excl': ..., 'detect': ..., 'n_clean': ..., 'n_bad': ...}}``.
    """
    out = {}
    for method in methods:
        fp = tp = n_clean = n_bad = 0
        for r in results:
            keep = np.array(r['keep'][method])
            bad = np.array(r['bad'])
            fp += int((~keep & ~bad).sum())
            tp += int((~keep & bad).sum())
            n_clean += int((~bad).sum())
            n_bad += int(bad.sum())
        out[method] = {
            'clean_excl': fp / max(n_clean, 1),
            'detect': tp / max(n_bad, 1) if n_bad else np.nan,
            'n_clean': n_clean,
            'n_bad': n_bad,
        }
    return out


def clean_run_exclusions(results, methods):
    """Per method: fraction of runs without corruption that lose at least one pair."""
    clean = [r for r in results if not any(r['bad'])]
    return {
        m: float(np.mean([not all(r['keep'][m]) for r in clean])) if clean else np.nan
        for m in methods
    }


def detection_by_kind(results, methods, kinds):
    """Per method and artifact kind: fraction of corrupted pairs excluded."""
    out = {}
    for method in methods:
        out[method] = {}
        for kind in kinds:
            hits = [
                not r['keep'][method][i]
                for r in results
                for i in range(r['n_pairs'])
                if r['bad'][i] and r['kind'][i] == kind
            ]
            out[method][kind] = (float(np.mean(hits)) if hits else np.nan, len(hits))
    return out


def detection_by_severity(results, method, kind, edges):
    """Fraction of corrupted pairs of one kind excluded, in bins of severity."""
    sev, hit = [], []
    for r in results:
        for i in range(r['n_pairs']):
            if r['bad'][i] and r['kind'][i] == kind:
                sev.append(r['severity'][i])
                hit.append(not r['keep'][method][i])
    sev, hit = np.array(sev), np.array(hit)
    rates, counts = [], []
    for lo, hi in pairwise(edges):
        sel = (sev >= lo) & (sev < hi)
        rates.append(float(hit[sel].mean()) if sel.any() else np.nan)
        counts.append(int(sel.sum()))
    return rates, counts


def outcome_table(results, methods, outcome, key='design'):
    """Mean outcome per (key value, n_bad) and method, with the paired SE of the
    difference from no exclusion.

    Returns ``{(group, n_bad): {method: (mean, se_of_diff_from_none, n)}}``.
    """
    grouped = defaultdict(list)
    for r in results:
        grouped[(r[key], r['n_bad'])].append(r)
    table = {}
    for cell, runs in sorted(grouped.items()):
        base = np.array([r['metrics']['none'][outcome] for r in runs])
        table[cell] = {}
        for method in methods:
            values = np.array([r['metrics'][method][outcome] for r in runs], dtype=float)
            # The oracle's ΔM error is undefined where it excluded every pair at a delay.
            ok = np.isfinite(values) & np.isfinite(base)
            diff = values[ok] - base[ok]
            se = float(diff.std(ddof=1) / np.sqrt(diff.size)) if diff.size > 1 else np.nan
            mean = float(values[ok].mean()) if ok.any() else np.nan
            table[cell][method] = (mean, se, int(ok.sum()))
    return table


def threshold_tradeoff(results, method, thresholds):
    """Share of clean and corrupted pairs whose final severity exceeds each threshold.

    This is approximate: the severities are from the method's last scoring pass, which was
    scored against the fit without the pairs flagged at the default threshold, and the
    effect-size floors and the one-pair-per-delay rule are ignored.
    """
    clean, bad = [], []
    for r in results:
        if method not in r['scores']:
            continue
        severity = np.array(r['scores'][method]['severity'])
        flags = np.array(r['bad'])
        clean.extend(severity[~flags])
        bad.extend(severity[flags])
    clean, bad = np.array(clean), np.array(bad)
    return (
        [float(np.mean(clean > t)) for t in thresholds],
        [float(np.mean(bad > t)) if bad.size else np.nan for t in thresholds],
    )


def oracle_gain_recovered(table, method):
    """Per cell: (none - method) / (none - oracle), where the oracle improves on none."""
    out = {}
    for cell, row in table.items():
        gain = row['none'][0] - row['oracle'][0]
        out[cell] = (row['none'][0] - row[method][0]) / gain if gain > 0 else np.nan
    return out
