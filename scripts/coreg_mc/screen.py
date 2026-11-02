#!/usr/bin/env python3
"""Task 4's primary screen: J1 over its factors, with matched ablations, on development fixtures.

Plan: "Bounding the combinations" and the Task 4 rules (fixed 2026-11-02 before any J1 result).
J1 runs over perfusion image (P-winsor, P-robust) x perfusion target (smoothed GM, pseudo-CBF) x
perfusion weight (0.2, 0.4, 0.6) x masking (the stage masks, none), anatomy fixed: 24 variants.
Each is compared with the ablations built from the same factors: S1 (masking only), S2 and Q1
(perfusion image, target, masking). The frozen arms of ``run_ants.ARMS`` are the variants
``P-winsor``, ``gm``, 0.4 and masked.

Every arm starts from production's ``mri_coreg`` estimate in the fixture's work directory
(prepared by ``run_ants.py``); the pseudo-CBF target (``run_ants.pseudo_cbf_target``) is added
to it. Writes one row per fixture and arm with the factors and the error against truth.

Usage (in WSL, ``aslprep-pvc`` environment, from the repository root)::

    python scripts/coreg_mc/screen.py --fixtures <aslscan fixtures> --runs <fixture runs> \\
        --work <work dir> --out screen.tsv --names rc_pcasl ... [--parallel 6]
"""

import argparse
import itertools
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'hmc_probe'))

import freeze_comparators  # noqa: E402
import run_ants  # noqa: E402
from run_ants import cc, mi  # noqa: E402
from score_full import _transform_error  # noqa: E402

IMAGES = ('P-winsor', 'P-robust')
TARGETS = {'gm': 'gm_target_lr', 'pseudocbf': 'pseudocbf_target_lr'}
WEIGHTS = (0.2, 0.4, 0.6)
MASKING = (True, False)
STATIC = ('/w/t1w_brain_lr.nii.gz', '/w/conventional.nii.gz')


def variants():
    """Every arm of the screen: (name, arm, factors, stages, masked)."""
    out = []
    for masked in MASKING:
        tag = 'mask' if masked else 'nomask'
        out.append((f'S1_{tag}', 'S1', {'masked': masked}, [[mi(*STATIC, 1)]], masked))
        for image, target in itertools.product(IMAGES, TARGETS):
            perf = (f'/w/{TARGETS[target]}.nii.gz', f'/w/{image}.nii.gz')
            factors = {'image': image, 'target': target, 'masked': masked}
            base = f'{image}_{target}_{tag}'
            out.append((f'S2_{base}', 'S2', factors, [[cc(*perf, 1)]], masked))
            out.append((f'Q1_{base}', 'Q1', factors, [[mi(*STATIC, 1)], [cc(*perf, 1)]], masked))
            for w in WEIGHTS:
                stages = [[mi(*STATIC, round(1 - w, 1)), cc(*perf, w)]]
                name = f'J1_{base}_w{round(w * 10)}'
                out.append((name, 'J1', factors | {'weight': w}, stages, masked))

    return out


def prepare(name, args):
    """Add the pseudo-CBF target (full and ASL-resolution) to the fixture's work directory."""
    work = args.work / name
    anat = run_ants.anat_dir(args.fixtures / name, args.runs / name)
    gm = next(anat.glob('*_label-GM_probseg.nii.gz'))
    wm = next(anat.glob('*_label-WM_probseg.nii.gz'))
    run_ants.pseudo_cbf_target(work, gm, wm)
    run_ants.lowres(work, names=(('pseudocbf_target', 1),))


def run_fixture(name, args):
    work = args.work / name
    arms = variants()
    if not args.rescore:
        prepare(name, args)
        lines = []
        for arm_name, arm, _factors, stages, masked in arms:
            cmd = run_ants.ants_command(
                arm,
                '/w/init_itk.txt',
                f'/w/screen/{arm_name}_',
                stages=stages,
                masked=masked,
            )
            lines.append(f'{cmd} > /w/screen/{arm_name}.log 2>&1')

        (work / 'screen').mkdir(exist_ok=True)
        run_ants.run_script(work, 'screen.sh', lines)

    truth, points, _production, voxel = run_ants.fixture_truth(
        args.fixtures / name,
        args.runs / name,
    )
    rows = []
    for arm_name, arm, factors, _stages, _masked in arms:
        est = run_ants.read_ants(work, arm, work / 'screen' / f'{arm_name}_0GenericAffine.mat')
        row = {'fixture': name, 'name': arm_name, 'arm': arm, **factors, 'voxel_mm': voxel}
        if est is None:
            row['problem'] = 'no transform'
        else:
            err = _transform_error(est, truth, points)
            row |= {'rms_mm': err['rms_mm'], 'rms_voxels': err['rms_mm'] / voxel}

        rows.append(row)

    print(f'{name}: done', flush=True)
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--fixtures', type=Path, required=True)
    p.add_argument('--runs', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--names', nargs='+', required=True)
    p.add_argument('--parallel', type=int, default=4)
    p.add_argument('--rescore', action='store_true', help='score results already on disk')
    args = p.parse_args(argv)
    names = freeze_comparators.development_fixtures(args.names)
    with ThreadPoolExecutor(args.parallel) as pool:
        rows = [r for rows in pool.map(lambda n: run_fixture(n, args), names) for r in rows]

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    missing = df[df['rms_mm'].isna()] if 'rms_mm' in df else df
    print(f'{len(df)} rows, {len(missing)} without a transform')


if __name__ == '__main__':
    main()
