#!/usr/bin/env python3
"""List the ASL references to evaluate, from earlier ASLPrep runs, as ``jobs.json``.

Two sources, both left by the head-motion-correction probe (``scripts/hmc_probe`` on the
``worktree-hmc-probe`` branch):

- **real**: ``coreg-5b-real/<dataset>__<run>/``, one per aslprep-zoo run. Each holds the run's
  ``desc-coreg_aslref`` (``conventional.nii.gz``), its production SynthStrip mask
  (``mask.nii.gz``), the T1w brain mask, the mean perfusion image (``P-mean.nii.gz``, unmasked)
  and FSL-BBR registrations of several sources to the T1w (``<source>_bbr_<sign>.mat``), scored
  in ``coreg_5b_real.tsv``. Distortion correction is taken to have run when the panel manifest
  lists fieldmaps for the run.
- **sim**: aslscan fixtures (``fixtures/aslscan/<name>``) with their baseline ASLPrep run
  (``fixture-runs/baseline/<name>``), which carry the ground truth.

Usage::

    python jobs.py --probe /home/tsalo/hmc-probe --out jobs.json
"""

import argparse
import json
from pathlib import Path

import pandas as pd


def real_jobs(probe):
    root = probe / 'coreg-5b-real'
    manifest = pd.read_csv(probe / 'panel' / 'manifest.tsv', sep='\t', dtype=str).fillna('')
    manifest = manifest.set_index(manifest['dataset'] + '__' + manifest['run'])
    scores = pd.read_csv(root / 'coreg_5b_real.tsv', sep='\t')
    jobs = []
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        if not (d / 'conventional.nii.gz').exists() or d.name not in manifest.index:
            print(f'skipping {d.name}: no reference or no manifest row')
            continue
        row = manifest.loc[d.name]
        dataset, run = d.name.split('__', 1)
        sc = scores[(scores['dataset'] == dataset) & (scores['run'] == run)]
        fallback = {f'{r.source}/{r.sign}': bool(r.fallback) for r in sc.itertuples(index=False)}
        jobs.append(
            {
                'name': d.name,
                'kind': 'real',
                'dataset': dataset,
                'run': run,
                'family': row['family'],
                'acquisition': row['acquisition'],
                'zooms_mm': json.loads(row['zooms_mm']),
                'sdc': bool(json.loads(row['fieldmaps'] or '[]')),
                'ref': str(d / 'conventional.nii.gz'),
                'prod_mask': str(d / 'mask.nii.gz'),
                't1w_mask': str(d / 't1w_mask.nii.gz'),
                'perfusion': str(d / 'P-mean.nii.gz'),
                'dir': str(d),
                'fallback': fallback,
                'perfusion_eligible': bool(sc['elig_eligible'].iloc[0]) if len(sc) else None,
            }
        )
    return jobs


def sim_jobs(probe):
    jobs = []
    for fx in sorted((probe / 'fixtures' / 'aslscan').iterdir()):
        if fx.name.startswith('_') or not fx.is_dir():
            continue
        perf = probe / 'fixture-runs' / 'baseline' / fx.name / 'sub-01' / 'perf'
        refs = list(perf.glob('*_desc-coreg_aslref.nii.gz'))
        if len(refs) != 1:
            print(f'skipping {fx.name}: {len(refs)} coreg references')
            continue
        stem = refs[0].name.split('_desc-coreg_aslref')[0]
        acq = stem.split('_acq-')[1]
        anat = fx / 'derivatives' / 'anat' / 'sub-01' / 'anat'
        jobs.append(
            {
                'name': f'sim__{fx.name}',
                'kind': 'sim',
                'fixture': str(fx),
                'acq': acq,
                'sdc': (fx / 'sub-01' / 'fmap').exists(),
                'ref': str(refs[0]),
                'prod_mask': str(perf / f'{stem}_desc-brain_mask.nii.gz'),
                't1w_mask': str(anat / 'sub-01_desc-brain_mask.nii.gz'),
                'aslprep_dir': str(perf.parents[1]),
            }
        )
    return jobs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    jobs = real_jobs(args.probe) + sim_jobs(args.probe)
    args.out.write_text(json.dumps(jobs, indent=1))
    print(
        f'{sum(j["kind"] == "real" for j in jobs)} real, {sum(j["kind"] == "sim" for j in jobs)} sim'
    )


if __name__ == '__main__':
    main()
