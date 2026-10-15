#!/usr/bin/env python3
"""Task 5b: register each reference candidate to the T1w and score it against truth.

For each probe fixture with a finished baseline run (``run_fixtures.py``):

- sources: the conventional reference (``desc-coreg_aslref``: production's unwarped,
  bias-corrected reference, D0/D1 for C0) and the D3 perfusion levels, built by
  ``perfusion.py`` from the run's native-space preprocessed series (``desc-preproc_asl``,
  already motion- and distortion-corrected) with the run's brain mask;
- registration as fMRIPrep's FSL path does it (``init_fsl_bbr_wf``): ``mri_coreg``
  (6 dof, ``--sep 4 --ftol 0.0001 --linmintol 0.01``) to the brain-masked T1w, converted to an
  FSL initialization, then ``flirt -cost bbr -dof 6 -basescale 1`` with FSL's ``bbr.sch`` and
  the white-matter mask. ``P-masked`` adds ``-inweight`` (the brain mask); the opposite BBR
  sign (``-bbrslope 0.5``) is run for every source as a diagnostic;
- fMRIPrep's fallback check is reproduced approximately (Frobenius norm of the difference of
  the BBR and initial RAS matrices above 15 keeps the initial transform);
- scoring: the estimated T1w -> reference pull against the truth ``F R^-1`` (``score_full.py``),
  RMS displacement over the brain; and the source's grey/white contrast in the truth tissue map.

Usage::

    python coreg_experiment.py --fixtures <dir>/aslscan --runs <fixture-runs/baseline> \\
        --names rc_pcasl ... --work <dir> --license <license.txt> --out <tsv>
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

from candidates import _fsl_to_world  # noqa: E402
from perfusion import LEVELS, combine, eligibility, select_inputs  # noqa: E402
from run_panel import IMAGE  # noqa: E402
from score import polar_mean  # noqa: E402
from score_full import _load_series, _one, _transform_error  # noqa: E402

from aslprep.tests import truth_scoring as ts  # noqa: E402

FALLBACK_NORM = 15.0


def _wm_mask(anat, out):
    """White matter from the dseg: the label that best matches the WM probability map."""
    dseg = nb.load(next(anat.glob('*_dseg.nii.gz')))
    labels = np.asanyarray(dseg.dataobj).astype(int)
    wm_p = np.asanyarray(nb.load(next(anat.glob('*_label-WM_probseg.nii.gz'))).dataobj) > 0.5
    best = max((v for v in np.unique(labels) if v), key=lambda v: ((labels == v) & wm_p).sum())
    nb.Nifti1Image((labels == best).astype(np.uint8), dseg.affine).to_filename(out)
    return int(best)


def _native(perf, suffix):
    """The one native-space (no `space-` entity) output with `suffix`."""
    hits = [p for p in perf.glob(f'*_{suffix}') if '_space-' not in p.name]
    if len(hits) != 1:
        raise FileNotFoundError(f'{len(hits)} native {suffix} files in {perf}')
    return hits[0]


def prepare(fixture_dir, run_dir, work):
    """Write the T1w brain, WM mask, brain mask and every source image to ``work``."""
    work.mkdir(parents=True, exist_ok=True)
    anat = fixture_dir / 'derivatives' / 'anat' / 'sub-01' / 'anat'
    t1 = nb.load(next(anat.glob('*_desc-preproc_T1w.nii.gz')))
    t1_mask = np.asanyarray(nb.load(next(anat.glob('*_desc-brain_mask.nii.gz'))).dataobj) > 0
    brain = np.asanyarray(t1.dataobj, np.float32) * t1_mask
    nb.Nifti1Image(brain, t1.affine).to_filename(work / 't1w_brain.nii.gz')
    wm_label = _wm_mask(anat, work / 'wm.nii.gz')

    perf = run_dir / 'sub-01' / 'perf'
    shutil.copyfile(_native(perf, 'desc-coreg_aslref.nii.gz'), work / 'conventional.nii.gz')
    mask_img = nb.load(_native(perf, 'desc-brain_mask.nii.gz'))
    mask = np.asanyarray(mask_img.dataobj) > 0
    mask_img = nb.Nifti1Image(mask.astype(np.uint8), mask_img.affine)
    mask_img.to_filename(work / 'mask.nii.gz')
    pre = nb.load(_native(perf, 'desc-preproc_asl.nii.gz'))
    data = np.asanyarray(pre.dataobj, np.float32)
    fperf = fixture_dir / 'sub-01' / 'perf'
    context = pd.read_table(next(fperf.glob('*_aslcontext.tsv')))['volume_type'].tolist()
    meta = json.loads(next(fperf.glob('*_asl.json')).read_text())
    if len(context) != data.shape[-1]:
        raise ValueError(f'{len(context)} context rows for {data.shape[-1]} preprocessed volumes')
    pld = np.atleast_1d(np.asarray(meta.get('PostLabelingDelay', 0.0), float))
    pld = np.repeat(pld, len(context)) if pld.size == 1 else pld
    stack, dtype, input_plds = select_inputs(data, context, list(pld))
    record = eligibility(stack, dtype, input_plds, mask, pre.header.get_zooms()[:3])
    for level in LEVELS:
        if level != 'masked':  # the same image as winsor, registered with -inweight
            img = combine(stack, level, mask)
            nb.Nifti1Image(img.astype(np.float32), pre.affine).to_filename(
                work / f'P-{level}.nii.gz'
            )
    return {'wm_label': wm_label, 'eligibility': record}


SOURCES = {  # name: (image, extra FLIRT arguments)
    'conventional': ('conventional', []),
    'P-mean': ('P-mean', []),
    'P-robust': ('P-robust', []),
    'P-winsor': ('P-winsor', []),
    'P-masked': ('P-winsor', ['-inweight', '/w/mask.nii.gz']),
}


def register(work, license_file):
    """Run mri_coreg and FLIRT-BBR (both BBR signs) for every source, in the ASLPrep image."""
    lines = ['set -e', 'export FS_LICENSE=/lic/license.txt', 'S=$FSLDIR/etc/flirtsch/bbr.sch']
    for name, (image, extra) in SOURCES.items():
        src = f'/w/{image}.nii.gz'
        lines += [
            f'mri_coreg --mov {src} --ref /w/t1w_brain.nii.gz --reg /w/{name}_init.lta --dof 6 '
            f'--sep 4 --ftol 0.0001 --linmintol 0.01 > /w/{name}_mri_coreg.log 2>&1',
            f'lta_convert --inlta /w/{name}_init.lta --outfsl /w/{name}_init.mat > /dev/null',
        ]
        for sign, slope in (('default', []), ('opposite', ['-bbrslope', '0.5'])):
            args = ' '.join([*extra, *slope])
            lines.append(
                f'flirt -in {src} -ref /w/t1w_brain.nii.gz -cost bbr -dof 6 -basescale 1 '
                f'-wmseg /w/wm.nii.gz -init /w/{name}_init.mat -omat /w/{name}_bbr_{sign}.mat '
                f'-schedule $S {args}'
            )
    (work / 'register.sh').write_text('\n'.join(lines) + '\n')
    cmd = ['docker', 'run', '--rm', '-u', f'{os.getuid()}:{os.getgid()}', '-e', 'HOME=/tmp']
    cmd += ['-v', f'{work}:/w', '-v', f'{license_file.parent}:/lic:ro']
    cmd += ['--entrypoint', 'bash', IMAGE, '/w/register.sh']
    subprocess.run(cmd, check=True, capture_output=True)


def _pull(mat_file, src_img, t1_img):
    """FSL matrix (source -> T1w) as a RAS pull from T1w points to source points."""
    push = _fsl_to_world(t1_img) @ np.loadtxt(mat_file) @ np.linalg.inv(_fsl_to_world(src_img))
    return np.linalg.inv(push)


def score(fixture_dir, run_dir, work, prep):
    asl = next((fixture_dir / 'sub-01' / 'perf').glob('*_asl.nii.gz'))
    acq = asl.name.split('_acq-')[1].split('_')[0]
    fixture = ts.Fixture(fixture_dir, acq)
    _, points = fixture.voxel_centers(fixture.brain)
    perf = run_dir / 'sub-01' / 'perf'
    hmc = _load_series(next(perf.glob(f'*_acq-{acq}_from-orig_to-aslref*_xfm.txt')))
    poses = fixture.poses()
    keep = [i for i, t in enumerate(fixture.context) if t in ('control', 'label', 'm0scan')]
    keep = keep if len(keep) == len(hmc) else list(range(len(hmc)))
    frame = polar_mean([np.linalg.inv(h) @ poses[i] for h, i in zip(hmc, keep, strict=True)])
    truth = frame @ np.linalg.inv(fixture.R)
    production = _one(ts.find_xfm(run_dir, acq, 'from-aslref_to-T1w'))
    t1_img = nb.load(work / 't1w_brain.nii.gz')
    gm_wm = _contrast(fixture, work)
    rows = [
        {
            'fixture': fixture_dir.name,
            'source': 'production (ASLPrep run)',
            'sign': 'default',
            **_transform_error(production, truth, points),
        }
    ]
    for name, (image, _) in SOURCES.items():
        src_img = nb.load(work / f'{image}.nii.gz')
        init = _pull(work / f'{name}_init.mat', src_img, t1_img)
        for sign in ('default', 'opposite'):
            bbr = _pull(work / f'{name}_bbr_{sign}.mat', src_img, t1_img)
            fallback = bool(np.linalg.norm(bbr - init) > FALLBACK_NORM)
            rows.append(
                {
                    'fixture': fixture_dir.name,
                    'source': name,
                    'sign': sign,
                    **_transform_error(init if fallback else bbr, truth, points),
                    'bbr_rms_mm': _transform_error(bbr, truth, points)['rms_mm'],
                    'init_rms_mm': _transform_error(init, truth, points)['rms_mm'],
                    'fallback': fallback,
                    'gm_wm_ratio': gm_wm.get(image),
                }
            )
    for r in rows:
        r.update({f'elig_{k}': v for k, v in prep['eligibility'].items()})
    return rows


def _contrast(fixture, work):
    """Grey/white mean ratio of each source in the truth tissue maps (still fixtures: F ~ I)."""
    gm = fixture.pv['GM'] >= 0.7
    wm = fixture.pv['WM'] >= 0.7
    out = {}
    for image in {v[0] for v in SOURCES.values()}:
        img = nb.load(work / f'{image}.nii.gz')
        data = np.asanyarray(img.dataobj, np.float32)
        if data.shape != gm.shape:
            out[image] = None
            continue
        out[image] = float(data[gm].mean() / data[wm].mean())
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--fixtures', type=Path, required=True)
    parser.add_argument('--runs', type=Path, required=True)
    parser.add_argument('--names', nargs='+', required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--license', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    rows = []
    for name in args.names:
        work = args.work / name
        shutil.rmtree(work, ignore_errors=True)
        prep = prepare(args.fixtures / name, args.runs / name, work)
        register(work, args.license)
        rows += score(args.fixtures / name, args.runs / name, work, prep)
        print(f'{name}: done', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(
            df.pivot_table(
                index=['source', 'sign'], columns='fixture', values='rms_mm', aggfunc='first'
            ).round(2)
        )
        print(
            df.query("sign == 'default'")[
                ['fixture', 'source', 'init_rms_mm', 'bbr_rms_mm', 'fallback', 'gm_wm_ratio']
            ]
            .round(2)
            .to_string()
        )


if __name__ == '__main__':
    main()
