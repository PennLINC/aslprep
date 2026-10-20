"""Truth-free consistency tests of partial volume correction on real ASLPrep outputs.

Two tests, both on ASLPrep's parcellated CBF tables (``*_atlas-<atlas>_*cbf.tsv``):

``resolution``
    The same subject scanned at two resolutions in one session (kashyap: 2 mm ``hires`` and
    3.4 mm ``lores`` PCASL). Uncorrected parcel CBF should fall at the coarser resolution,
    because each parcel mixes in more white matter and CSF. A partial volume correction that
    works should make the two agree.
``retest``
    The same protocol in several sessions (ds007116: three sessions). Reports the
    within-parcel coefficient of variation across sessions, raw and with each session's
    global median removed, and the reliability of the parcel pattern (ICC(3,1), with
    parcels as targets and sessions as raters, on raw and log values).

``contrast``
    The GM/WM CBF contrast of every run, from the native-space maps. It is free of each
    run's calibration scale, which differs between protocols (kashyap's two resolutions also
    differ in labeling duration and readout), so it is the cleaner resolution test: after a
    correction that works, it should not depend on resolution, while the uncorrected contrast
    falls at the coarser resolution.

Every method is scored on the same parcels: those with a positive value, and at least 50%
coverage, for every method in every run compared.

Usage::

    python consistency.py resolution DERIV_DIR --subject 01 --low acq-lores --high acq-hires
    python consistency.py retest DERIV_DIR --subject 20721
    python consistency.py contrast DERIV_DIR --subject 01

Writes a JSON summary to ``--out`` (default: stdout).
"""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

#: Parcellated maps compared, by their ``desc`` entity (``None`` is the uncorrected CBF).
METHODS = {
    'uncorrected': None,
    'basil': 'basil',
    'lr_gm': 'lrGM',
    'basil_gm': 'basilGM',
}
#: Coverage table used for each method (BASIL's GM map has none of its own).
COVERAGE = {'uncorrected': None, 'basil': 'basil', 'lr_gm': 'lrGM', 'basil_gm': 'basil'}
MIN_COVERAGE = 0.5


def _run_key(path, atlas):
    """The run's entities before the atlas entity, e.g. ``sub-01_ses-20_acq-hires_dir-PA_run-1``."""
    return path.name.split(f'_atlas-{atlas}')[0]


def load_tables(deriv_dir, subject, atlas):
    """Parcel CBF (and coverage) for every run and method: ``{run: {method: Series}}``."""
    perf_files = sorted(Path(deriv_dir).glob(f'sub-{subject}/**/perf/*_atlas-{atlas}_*cbf.tsv'))
    runs = {}
    for path in perf_files:
        match = re.search(rf'_atlas-{atlas}_(?:desc-([A-Za-z0-9]+)_)?cbf\.tsv$', path.name)
        if not match:
            continue
        desc = match.group(1)
        method = next((m for m, d in METHODS.items() if d == desc), None)
        if method is None:
            continue
        values = pd.read_table(path).iloc[0].astype(float)
        runs.setdefault(_run_key(path, atlas), {})[method] = values
    coverage = {}
    for path in sorted(
        Path(deriv_dir).glob(f'sub-{subject}/**/perf/*_atlas-{atlas}_*coverage.tsv')
    ):
        match = re.search(rf'_atlas-{atlas}_(?:desc-([A-Za-z0-9]+)_)?coverage\.tsv$', path.name)
        if not match:
            continue
        table = pd.read_table(path, index_col=0).iloc[:, 0].astype(float)
        coverage.setdefault(_run_key(path, atlas), {})[match.group(1)] = table
    return runs, coverage


def _valid_one(runs, coverage, keys, method):
    """Parcels with finite, positive values and enough coverage in every run of ``keys``."""
    mask = None
    for key in keys:
        values = runs[key].get(method)
        if values is None:
            return None
        ok = np.isfinite(values) & (values > 0)
        cov = coverage.get(key, {}).get(COVERAGE[method])
        if cov is not None:
            ok &= cov.reindex(values.index).fillna(0) >= MIN_COVERAGE
        mask = ok if mask is None else mask & ok
    return mask


def _valid(runs, coverage, keys):
    """Parcels valid for every method present in every run, so methods share one parcel set."""
    masks = [_valid_one(runs, coverage, keys, m) for m in METHODS]
    masks = [m for m in masks if m is not None]
    if not masks:
        return None
    mask = masks[0]
    for other in masks[1:]:
        mask = mask & other.reindex(mask.index).fillna(False).astype(bool)
    return mask


def icc_3_1(table):
    """ICC(3,1) of an (n_targets, k_raters) array: consistency of single measurements."""
    data = np.asarray(table, float)
    n, k = data.shape
    grand = data.mean()
    ss_rows = k * ((data.mean(axis=1) - grand) ** 2).sum()
    ss_cols = n * ((data.mean(axis=0) - grand) ** 2).sum()
    ss_err = ((data - grand) ** 2).sum() - ss_rows - ss_cols
    ms_rows = ss_rows / (n - 1)
    ms_err = ss_err / ((n - 1) * (k - 1))
    return float((ms_rows - ms_err) / (ms_rows + (k - 1) * ms_err))


def resolution_test(runs, coverage, low, high):
    """Compare each session's low- and high-resolution runs, per method."""
    out = {}
    for low_key in sorted(k for k in runs if low in k):
        high_key = low_key.replace(low, high)
        if high_key not in runs:
            continue
        session = re.search(r'ses-[^_]+', low_key)
        name = session.group(0) if session else low_key
        out[name] = {}
        valid = _valid(runs, coverage, (low_key, high_key))
        if valid is None or valid.sum() < 10:
            continue
        for method in METHODS:
            if method not in runs[low_key] or method not in runs[high_key]:
                continue
            lo = runs[low_key][method][valid]
            hi = runs[high_key][method][valid]
            log_ratio = np.log(lo / hi)
            out[name][method] = {
                'n_parcels': int(valid.sum()),
                'median_low': float(lo.median()),
                'median_high': float(hi.median()),
                'median_ratio': float(np.exp(log_ratio.median())),
                'mad_log_ratio': float((log_ratio - log_ratio.median()).abs().median()),
                'r': float(np.corrcoef(lo, hi)[0, 1]),
            }
    return out


def retest_test(runs, coverage, group=None, exclude=()):
    """Across-session reliability of each method's parcel values.

    ``group`` restricts the runs to those whose key contains it (for example ``acq-lores``),
    and runs whose key contains any of ``exclude`` are left out.
    """
    keys = sorted(
        k for k in runs if (group is None or group in k) and not any(e in k for e in exclude)
    )
    out = {'runs': keys}
    valid = _valid(runs, coverage, keys) if len(keys) >= 2 else None
    if valid is None or valid.sum() < 10:
        return out
    for method in METHODS:
        if any(method not in runs[k] for k in keys):
            continue
        table = pd.concat([runs[k][method][valid] for k in keys], axis=1)
        cv = table.std(axis=1, ddof=1) / table.mean(axis=1)
        medians = table.median(axis=0)
        # Regional reproducibility, with each session's global scale removed
        scaled = table / medians
        cv_scaled = scaled.std(axis=1, ddof=1) / scaled.mean(axis=1)
        out[method] = {
            'n_parcels': int(valid.sum()),
            'median_cv': float(cv.median()),
            'median_cv_scaled': float(cv_scaled.median()),
            'icc_3_1': icc_3_1(table.to_numpy()),
            'icc_3_1_log': icc_3_1(np.log(table.to_numpy())),
            'session_medians': [float(v) for v in medians],
            'cv_of_session_medians': float(medians.std(ddof=1) / medians.mean()),
        }
    return out


#: Native maps for the contrast test: (GM map desc, WM map desc), ``None`` = uncorrected CBF.
CONTRAST_MAPS = {
    'uncorrected': (None, None),
    'basil': ('basil', 'basil'),
    'lr': ('lrGM', 'lrWM'),
    'basil_pvc': ('basilGM', 'basilWM'),
}


def _native(perf, stem, desc, suffix):
    name = f'{stem}_desc-{desc}_{suffix}.nii.gz' if desc else f'{stem}_{suffix}.nii.gz'
    path = perf / name
    return path if path.exists() else None


def contrast_test(deriv_dir, subject, pure_gm=0.7, pure_wm=0.9, min_pv=0.25):
    """GM/WM CBF contrast of every run, from the native-space maps.

    The contrast is free of the run's calibration scale, so it can be compared across
    protocols. Uncorrected maps (ASLPrep's and BASIL's) use voxels whose PV map is at least
    ``pure_gm`` GM or ``pure_wm`` WM; corrected maps use their estimates wherever the tissue's
    PV is at least ``min_pv``.
    """
    import nibabel as nb

    out = {}
    for mask_path in sorted(
        Path(deriv_dir).glob(f'sub-{subject}/**/perf/*_desc-brain_mask.nii.gz')
    ):
        if '_space-' in mask_path.name:
            continue
        perf = mask_path.parent
        stem = mask_path.name.replace('_desc-brain_mask.nii.gz', '')
        pv_paths = {t: perf / f'{stem}_label-{t}_desc-psf_probseg.nii.gz' for t in ('GM', 'WM')}
        if not all(p.exists() for p in pv_paths.values()):
            continue
        mask = np.asanyarray(nb.load(mask_path).dataobj) > 0
        pv = {t: nb.load(p).get_fdata() for t, p in pv_paths.items()}
        run = {}
        for method, (gm_desc, wm_desc) in CONTRAST_MAPS.items():
            gm_path = _native(perf, stem, gm_desc, 'cbf')
            wm_path = _native(perf, stem, wm_desc, 'cbf')
            if gm_path is None or wm_path is None:
                continue
            gm = nb.load(gm_path).get_fdata()
            wm = nb.load(wm_path).get_fdata()
            gm = gm[..., 0] if gm.ndim == 4 else gm
            wm = wm[..., 0] if wm.ndim == 4 else wm
            corrected = gm_desc != wm_desc
            gm_sel = mask & (pv['GM'] >= (min_pv if corrected else pure_gm))
            wm_sel = mask & (pv['WM'] >= (min_pv if corrected else pure_wm))
            gm_sel &= np.isfinite(gm) & (gm != 0)
            wm_sel &= np.isfinite(wm) & (wm != 0)
            if gm_sel.sum() < 50 or wm_sel.sum() < 50:
                continue
            gm_med, wm_med = float(np.median(gm[gm_sel])), float(np.median(wm[wm_sel]))
            run[method] = {
                'gm_median': gm_med,
                'wm_median': wm_med,
                'gm_wm_ratio': gm_med / wm_med,
                'n_gm': int(gm_sel.sum()),
                'n_wm': int(wm_sel.sum()),
            }
        out[stem] = run
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('test', choices=('resolution', 'retest', 'contrast'))
    parser.add_argument('deriv_dir', type=Path)
    parser.add_argument('--subject', required=True)
    parser.add_argument('--atlas', default='4S156Parcels')
    parser.add_argument('--low', default='acq-lores', help='resolution: low-resolution entity')
    parser.add_argument('--high', default='acq-hires', help='resolution: high-resolution entity')
    parser.add_argument('--group', default=None, help='retest: restrict runs to this entity')
    parser.add_argument(
        '--exclude', nargs='*', default=(), help='retest: leave out runs containing these'
    )
    parser.add_argument('--out', type=Path, default=None)
    args = parser.parse_args(argv)

    runs, coverage = load_tables(args.deriv_dir, args.subject, args.atlas)
    if args.test == 'contrast':
        result = contrast_test(args.deriv_dir, args.subject)
    elif not runs:
        parser.error(f'no parcellated CBF tables for sub-{args.subject} and {args.atlas}')
    elif args.test == 'resolution':
        result = resolution_test(runs, coverage, args.low, args.high)
    else:
        result = retest_test(runs, coverage, args.group, args.exclude)
    result = {'test': args.test, 'atlas': args.atlas, 'subject': args.subject, 'result': result}
    text = json.dumps(result, indent=2)
    if args.out:
        args.out.write_text(text + '\n')
    else:
        sys.stdout.write(text + '\n')


if __name__ == '__main__':
    main()
