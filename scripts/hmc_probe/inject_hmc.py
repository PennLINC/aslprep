#!/usr/bin/env python3
"""Task 7: hand a candidate's HMC products to a fresh full ASLPrep run.

ASLPrep reuses a precomputed HMC reference and HMC transforms (``fit.py``, stages 1 and 2) and
builds everything after them anew: the separate M0's registration, the fieldmap registration,
the coregistration reference, the brain mask and the coregistration. Writing a harness output
(``aslref.nii.gz``, ``hmc_xfms.txt``) as a one-run derivatives dataset with ASLPrep's own names
therefore wires the candidate in, under the downstream policy the harness applied (D2 by
default: the candidate's own template for every consumer), with no change to the workflow.
Nothing else is precomputed, and every run gets fresh output and working directories.

``run_fixtures.py --hmc-derivs <out>/<label>`` then passes ``<out>/<label>/<fixture>`` as the
user dataset, after the atlas datasets. It holds a copy of the fixture's anatomical derivatives
as well. ASLPrep merges one functional cache per derivatives dataset, each replacing the
``transforms`` entry, and appends its atlas datasets after the user's, so by default the
injected transforms are wiped and HMC is recomputed (an ASLPrep bug: precomputed transforms are
never reused while atlases are enabled). ``run_fixtures.py`` names the atlas datasets first.

Usage::

    python inject_hmc.py --scores <fixture-scores> --fixtures <data dir>/aslscan \\
        --labels c0 c1-pld ... --names rc_pcasl ... --out <dir>
"""

import argparse
import json
import shutil
from pathlib import Path


def inject(harness_dir, fixture_dir, out_dir):
    harness_dir, fixture_dir, out_dir = Path(harness_dir), Path(fixture_dir), Path(out_dir)
    asl = next((fixture_dir / 'sub-01' / 'perf').glob('*_asl.nii.gz'))
    stem = asl.name.removesuffix('_asl.nii.gz')
    perf = out_dir / 'sub-01' / 'perf'
    if out_dir.exists():
        shutil.rmtree(out_dir)
    anat = fixture_dir / 'derivatives' / 'anat'
    if anat.exists():
        shutil.copytree(anat, out_dir, ignore=shutil.ignore_patterns('dataset_description.json'))
    perf.mkdir(parents=True, exist_ok=True)
    harness = json.loads((harness_dir / 'harness.json').read_text())
    (out_dir / 'dataset_description.json').write_text(
        json.dumps(
            {
                'Name': f'HMC probe: {harness["candidate"]} on {fixture_dir.name}',
                'BIDSVersion': '1.10.0',
                'DatasetType': 'derivative',
                'GeneratedBy': [{'Name': 'hmc_probe harness', 'Description': json.dumps(harness)}],
            },
            indent=2,
        )
        + '\n'
    )
    shutil.copyfile(harness_dir / 'aslref.nii.gz', perf / f'{stem}_desc-hmc_aslref.nii.gz')
    xfm = perf / f'{stem}_from-orig_to-aslref_mode-image_desc-hmc_xfm.txt'
    shutil.copyfile(harness_dir / 'hmc_xfms.txt', xfm)
    return xfm


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--scores', type=Path, required=True)
    parser.add_argument('--fixtures', type=Path, required=True, help='<data dir>/aslscan')
    parser.add_argument('--labels', nargs='+', required=True)
    parser.add_argument('--names', nargs='+', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    for label in args.labels:
        for name in args.names:
            src = args.scores / label / name
            if not (src / 'hmc_xfms.txt').exists():
                print(f'{label}/{name}: no harness output')
                continue
            print(f'{label}/{name}: {inject(src, args.fixtures / name, args.out / label / name)}')


if __name__ == '__main__':
    main()
