#!/usr/bin/env python3
"""Task 7: score every injected full run (metrics 5, 6 and 9-11) and tabulate.

For each ``<runs>/<label>/<fixture>`` whose ``DONE.json`` reports success: ``score_full`` with the
injected HMC transforms (mask, coregistration, M0 and fieldmap registration, end to end) and
``score_cbf`` against ``<runs>/oracle/<fixture>`` (CBF and ATT against the same quantification
without motion-correction error). Writes ``<out>/<label>/<fixture>/{full,cbf}.json`` and a table,
``<out>/task7.tsv``, one row per run.

Usage::

    python task7_score.py --runs <dir> --derivs <dir> --fixtures <dir> [<dir> ...] --out <dir>
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from score_cbf import score_cbf  # noqa: E402
from score_full import score_full  # noqa: E402

ROW = {
    'mask_dice': ('full', 'mask', 'dice'),
    'mask_extra_ml': ('full', 'mask', 'extra_ml'),
    'coreg_end_to_end_mm': ('full', 'coregistration_vs_aslref_frame', 'rms_mm'),
    'coreg_step_mm': ('full', 'coregistration_vs_reference_image', 'rms_mm'),
    'coreg_only_rms_voxels': ('full', 'frames', 'coreg', 'coreg_only_rms_voxels'),
    'motion_rms_median_mm': ('full', 'frames', 'motion', 'rms_error_median_mm'),
    'm0_mm': ('full', 'm0', 'rms_mm'),
    'fieldmap_mm': ('full', 'fieldmap', 'rms_mm'),
    'cbf_gm_bias_rel': ('cbf', 'cbf', 'GM', 'bias_rel_gm'),
    'cbf_wm_bias_rel': ('cbf', 'cbf', 'WM', 'bias_rel_gm'),
    'cbf_brain_rmse_rel': ('cbf', 'cbf', 'brain', 'rmse_rel_gm'),
    'cbf_edge_rmse_rel': ('cbf', 'cbf', 'shell_0_3mm', 'rmse_rel_gm'),
    'cbf_brain_p95': ('cbf', 'cbf', 'brain', 'p95_abs'),
    'att_brain_rmse_rel': ('cbf', 'att', 'brain', 'rmse_rel_gm'),
    'att_at_bound': ('cbf', 'att', 'fraction_at_bound'),
    'oracle_gm_cbf': ('cbf', 'cbf', 'GM', 'oracle_mean'),
}
#: The oracle's grey-matter CBF must be positive for metric 5 to mean anything. On
#: ``mo_pcasl_spikes`` it is not: the simulated head moves through a B0 field fixed to the
#: scanner, so 2 mm spikes change the distortion and even the true poses leave the spike pairs
#: corrupt (the harness reports the same for its oracle). Such fixtures are dropped from metric 5.
CBF_COLUMNS = [k for k in ROW if k.startswith(('cbf_', 'att_'))]
#: The CBF target must be in the run's own frame, or the comparison resamples the target and not
#: the run. Production and C9 use production's reference (a separate M0 can sit millimetres from
#: the volumes), so they are compared with the oracle built in production's frame; the template
#: arms' frames are within a few hundredths of a millimetre of the static one.
ORACLE_FOR = {'c0': 'oracle-c0', 'c9': 'oracle-c0', 'c9-late': 'oracle-c0'}


def _get(d, path):
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return float('nan')
        d = d[k]
    return d


def hmc_file(derivs, label, fixture):
    hits = sorted((derivs / label / fixture / 'sub-01' / 'perf').glob('*_desc-hmc_xfm.txt'))
    return hits[0] if len(hits) == 1 else None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--runs', type=Path, required=True)
    parser.add_argument('--derivs', type=Path, required=True)
    parser.add_argument('--fixtures', type=Path, nargs='+', required=True, help='aslscan dirs')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--rescore', action='store_true')
    args = parser.parse_args(argv)
    fixture_dirs = {p.name: p for root in args.fixtures for p in root.iterdir() if p.is_dir()}
    rows = []
    for done in sorted(args.runs.glob('*/*/DONE.json')):
        run, label, name = done.parent, done.parent.parent.name, done.parent.name
        if json.loads(done.read_text())['returncode'] != 0 or name not in fixture_dirs:
            continue
        hmc = hmc_file(args.derivs, label, name)
        oracle = ORACLE_FOR.get(label, 'oracle')
        if not (args.runs / oracle / name / 'DONE.json').exists():
            oracle = 'oracle'
        oracle_hmc = hmc_file(args.derivs, oracle, name)
        oracle_run = args.runs / oracle / name
        out = args.out / label / name
        out.mkdir(parents=True, exist_ok=True)
        rec = {}
        for kind in ('full', 'cbf'):
            f = out / f'{kind}.json'
            if f.exists() and not args.rescore:
                rec[kind] = json.loads(f.read_text())
                continue
            try:
                if kind == 'full':
                    rec[kind] = score_full(fixture_dirs[name], run, hmc)
                elif oracle_hmc and (oracle_run / 'DONE.json').exists() and hmc:
                    rec[kind] = score_cbf(fixture_dirs[name], run, hmc, oracle_run, oracle_hmc)
                else:
                    continue
            except Exception as exc:  # noqa: BLE001 - record and continue
                rec[kind] = {'problems': [repr(exc)]}
            f.write_text(json.dumps(rec[kind], indent=2, default=float) + '\n')
        rows.append(
            {'label': label, 'fixture': name, 'oracle': oracle}
            | {k: _get(rec, p) for k, p in ROW.items()}
            | {'frame_vs_oracle_mm': _get(rec, ('cbf', 'frame_offset_mm'))}
        )
        print(f'{label}/{name}: scored', flush=True)
    df = pd.DataFrame(rows)
    invalid = ~(df.oracle_gm_cbf > 0)
    if invalid.any():
        print('metric 5 dropped (oracle GM CBF not positive):', sorted(set(df.fixture[invalid])))
    df.loc[invalid, CBF_COLUMNS] = float('nan')
    df.to_csv(args.out / 'task7.tsv', sep='\t', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 30):
        print(df.groupby('label').median(numeric_only=True).round(3).T.to_string())


if __name__ == '__main__':
    main()
