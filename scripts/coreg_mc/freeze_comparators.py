#!/usr/bin/env python3
"""Task 3 evidence for freezing R1's sign rule and R2's D3 level and eligibility rule.

Plan: "Candidates" and Task 3. Development data only: the development fixtures (truth) and the
development real runs (gross verdicts and inconsistency from ``run_real.py``, via
``batch_real.py``).

Per run (fixture or real) this collects:

- ``contrast``: the conventional reference's white-minus-grey contrast under the ``mri_coreg``
  initialization (the medians of the reference at eroded-WM and GM-core T1w voxels, pulled through
  the initialization, divided by the median inside the T1w brain). Available before BBR runs, so a
  rule on it can run in production. BBR's default sign expects grey brighter than white
  (``contrast < 0``).
- real runs: each FSL arm's adopted gross verdict, its median inconsistency over its partners
  (conventional frame), the pair noise floor, and the HMC probe's D3 eligibility record
  (``perfusion.eligibility``, recomputed by ``coreg_experiment.run_real`` with ``rescore``, which
  rewrites the work directory's inputs with identical content and registers nothing).
- fixtures: each FSL arm's RMS error against truth, and the same eligibility record.

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/freeze_comparators.py --fixtures <aslscan fixtures> \\
        --fixture-runs <fixture runs> --fixture-work <dir> --names rc_pcasl ... \\
        --dev-real <dir> --panel-runs <dir> --extension-runs <dir> --adapted <dir> \\
        --license <file> --out <dir>
"""

import argparse
import contextlib
import io
import json
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
from nibabel.affines import apply_affine
from scipy.ndimage import binary_erosion, map_coordinates

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'hmc_probe'))

import coreg_experiment  # noqa: E402
import run_ants  # noqa: E402
import run_real  # noqa: E402

#: Voxels of T1w-grid erosion applied to the WM mask (partial volume at the ASL resolution).
WM_EROSION = 2
#: GM core: the smoothed GM probability (``gm_target``) above this, outside the WM mask.
GM_CORE = 0.5


def contrast(work):
    """The conventional reference's WM-minus-GM contrast under the mri_coreg initialization."""
    est = coreg_experiment.estimates(work)
    if ('conventional', 'default') not in est:
        return np.nan
    init = est[('conventional', 'default')][2]
    ref = nb.load(work / 'conventional.nii.gz')
    data = np.asanyarray(ref.dataobj, np.float32)
    wm_img = nb.load(work / 'wm.nii.gz')
    wm = binary_erosion(np.asanyarray(wm_img.dataobj) > 0, iterations=WM_EROSION)
    gm = (np.asanyarray(nb.load(work / 'gm_target.nii.gz').dataobj) > GM_CORE) & ~(
        np.asanyarray(wm_img.dataobj) > 0
    )
    brain = np.asanyarray(nb.load(work / 't1w_mask.nii.gz').dataobj) > 0

    def sample(mask):
        world = apply_affine(wm_img.affine, np.argwhere(mask)[::3])
        aslref = apply_affine(init, world)
        ijk = apply_affine(np.linalg.inv(ref.affine), aslref).T
        values = map_coordinates(data, ijk, order=1, mode='constant', cval=np.nan)
        return np.nanmedian(values)

    return float((sample(wm) - sample(gm)) / sample(brain))


def fixtures(fixture_root, run_root, work_root, names, license_file):
    """One row per development fixture and FSL arm: RMS error against truth, and eligibility.

    ``coreg_experiment.run_fixture`` with ``rescore`` scores every source and sign already
    registered in the work directory (it rewrites the prepared inputs with identical content).
    """
    rows = []
    for name in names:
        fixture_dir, run_dir, work = fixture_root / name, run_root / name, work_root / name
        with contextlib.redirect_stdout(io.StringIO()):
            scored = coreg_experiment.run_fixture(
                fixture_dir,
                run_dir,
                work,
                license_file,
                rescore=True,
                anat_dir=run_ants.anat_dir(fixture_dir, run_dir),
            )
        c = contrast(work)
        for r in scored:
            arm = 'R0 production' if r['source'].startswith('production') else None
            rows.append(
                {
                    'fixture': name,
                    'arm': arm or f'FSL {r["source"]} {r["sign"]}',
                    'rms_mm': r['rms_mm'],
                    'fallback': r.get('fallback', False),
                    'contrast': c,
                }
                | {k: v for k, v in r.items() if k.startswith('elig_')},
            )
        print(f'{name}: contrast {c:+.3f}', flush=True)
    return pd.DataFrame(rows)


def eligibility(row, runs, adapted, work, license_file):
    """The D3 eligibility record of one real run (no registration; see the module docstring)."""
    record = run_real.record_for(row)
    with contextlib.redirect_stdout(io.StringIO()):
        rows = coreg_experiment.run_real(record, runs, adapted, work, license_file, rescore=True)
    if not rows:
        return {}
    return {k: v for k, v in rows[0].items() if k.startswith('elig_')}


def read_tables(folder, suffix):
    """Every ``*_<suffix>.tsv`` in ``folder``, concatenated (a group with no rows writes none)."""
    tables = []
    for f in sorted(folder.glob(f'*_{suffix}.tsv')):
        with contextlib.suppress(pd.errors.EmptyDataError):
            tables.append(pd.read_table(f))

    return pd.concat(tables, ignore_index=True)


def real_runs(dev_real, manifest, run_roots, adapted, license_file):
    """One row per development real run and FSL arm."""
    gross = read_tables(dev_real, 'gross')
    cons = read_tables(dev_real, 'consistency')
    cons = cons[cons['frame'] == 'conventional']
    # Every run's partners: each pair counts for both of its runs.
    keep = ['dataset', 'run', 'arm', 'inconsistency_mm', 'floor_mm']
    per_run = pd.concat(
        [cons.rename(columns={side: 'run'})[keep] for side in ('run_a', 'run_b')],
        ignore_index=True,
    )
    per_run = per_run.groupby(['dataset', 'run', 'arm'], as_index=False).agg(
        median_inconsistency_mm=('inconsistency_mm', 'median'),
        median_floor_mm=('floor_mm', 'median'),
    )
    out = gross[['dataset', 'run', 'arm', 'verdict', 'dice', 'n_partners', 'no_transform']].merge(
        per_run,
        on=['dataset', 'run', 'arm'],
        how='left',
    )
    out = out.drop_duplicates(['dataset', 'run', 'arm'], keep='last')
    out = out.merge(
        manifest[['dataset', 'run', 'family', 'source']],
        on=['dataset', 'run'],
        how='left',
    )
    runs = out[['dataset', 'run']].drop_duplicates()
    extra = []
    for ds, run in runs.itertuples(index=False):
        row = manifest[(manifest['dataset'] == ds) & (manifest['run'] == run)].iloc[0]
        work = dev_real / 'work' / f'{ds}__{run}'
        rec = {'dataset': ds, 'run': run, 'contrast': contrast(work)}
        rec |= eligibility(row, run_roots[row['source']], adapted, work, license_file)
        extra.append(rec)
        print(
            f'{ds} {run}: contrast {rec["contrast"]:+.3f} eligible {rec.get("elig_eligible")}',
            flush=True,
        )

    return out.merge(pd.DataFrame(extra), on=['dataset', 'run'], how='left')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--manifest', type=Path, default=HERE / 'manifest.tsv')
    p.add_argument('--dev-real', type=Path, required=True)
    p.add_argument('--fixtures', type=Path, required=True)
    p.add_argument('--fixture-runs', type=Path, required=True)
    p.add_argument('--fixture-work', type=Path, required=True)
    p.add_argument('--names', nargs='+', required=True)
    p.add_argument('--panel-runs', type=Path, required=True)
    p.add_argument('--extension-runs', type=Path, required=True)
    p.add_argument('--adapted', type=Path, required=True)
    p.add_argument('--license', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--reuse-fixtures', action='store_true', help='Read --out/fixtures.tsv.')
    args = p.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.reuse_fixtures and (args.out / 'fixtures.tsv').exists():
        fix = pd.read_table(args.out / 'fixtures.tsv')
    else:
        fix = fixtures(
            args.fixtures,
            args.fixture_runs,
            args.fixture_work,
            args.names,
            args.license,
        )
        fix.to_csv(args.out / 'fixtures.tsv', sep='\t', index=False)
    wide = fix.pivot(index='fixture', columns='arm', values='rms_mm')
    with pd.option_context('display.width', 300, 'display.max_columns', 30):
        print(wide.round(2).to_string())

    manifest = pd.read_table(args.manifest)
    roots = {'panel': args.panel_runs, 'extension': args.extension_runs}
    real = real_runs(args.dev_real, manifest, roots, args.adapted, args.license)
    real.to_csv(args.out / 'real.tsv', sep='\t', index=False)
    n_real = int(real[['dataset', 'run']].drop_duplicates().shape[0])
    print(json.dumps({'fixtures': len(fix), 'real_runs': n_real}))


if __name__ == '__main__':
    main()
