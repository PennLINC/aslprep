#!/usr/bin/env python3
"""Make the adapted copies of panel datasets listed in ``panel_config.ADAPTATIONS``.

Each copy gets an ``ADAPTED.md`` that lists every change and why. The zoo is never modified.

Usage::

    python scripts/hmc_probe/adapt_panel.py --zoo <aslprep-zoo> --out <dir>
"""

import argparse
import json
import shutil
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from panel_config import ADAPTATIONS


def _root(zoo, name):
    base = zoo / name
    return base / 'dset' if (base / 'dset' / 'dataset_description.json').exists() else base


def _remove_files(root, step):
    removed = sorted(root.glob(step['glob']))
    for f in removed:
        f.unlink()
    return [f'removed {f.relative_to(root)}' for f in removed]


def _unlink_m0(root, step):
    done = []
    for j in sorted(root.glob(step['glob'])):
        meta = json.loads(j.read_text())
        if meta.pop('IntendedFor', None) is not None:
            j.write_text(json.dumps(meta, indent=2) + '\n')
            done.append(f'removed IntendedFor from {j.relative_to(root)}')
    return done


def _drop_unpaired(root, step):
    """Drop the last unpaired control or label of every run whose counts differ."""
    done = []
    for ctx_file in sorted(root.glob('sub-*/**/perf/*_aslcontext.tsv')):
        ctx = pd.read_table(ctx_file, keep_default_na=False)
        types = ctx['volume_type'].tolist()
        n_c, n_l = types.count('control'), types.count('label')
        if n_c == n_l:
            continue
        extra = 'control' if n_c > n_l else 'label'
        drop = max(i for i, t in enumerate(types) if t == extra)
        asl_file = ctx_file.with_name(ctx_file.name.replace('_aslcontext.tsv', '_asl.nii.gz'))
        img = nb.load(asl_file)
        keep = [i for i in range(img.shape[3]) if i != drop]
        data = np.asanyarray(img.dataobj)[..., keep]
        nb.Nifti1Image(data, img.affine, img.header).to_filename(asl_file)
        ctx.drop(index=drop).to_csv(ctx_file, sep='\t', index=False)
        side = asl_file.with_name(asl_file.name.replace('.nii.gz', '.json'))
        meta = json.loads(side.read_text())
        for key, value in list(meta.items()):
            if isinstance(value, list) and len(value) == len(types) and key != 'SliceTiming':
                meta[key] = [v for i, v in enumerate(value) if i != drop]
        side.write_text(json.dumps(meta, indent=2) + '\n')
        done.append(f'dropped volume {drop} ({extra}) from {asl_file.relative_to(root)}')
    return done


OPS = {'remove_files': _remove_files, 'unlink_m0': _unlink_m0, 'drop_unpaired': _drop_unpaired}


def adapt(zoo, out, names=None):
    for name, steps in ADAPTATIONS.items():
        if names and name not in names:
            continue
        dest = out / name
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(
            _root(zoo, name), dest, ignore=shutil.ignore_patterns('sourcedata', '.git*')
        )
        lines = [
            f'# {name}: adapted copy for the HMC probe',
            '',
            f'Source: `{_root(zoo, name)}`',
            '',
        ]
        for step in steps:
            lines.append(f'- {step["why"]}')
            lines += [f'  - {d}' for d in OPS[step['op']](dest, step)]
        (dest / 'ADAPTED.md').write_text('\n'.join(lines) + '\n')
        print(f'{name} -> {dest}')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--zoo', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--only', nargs='*')
    args = parser.parse_args(argv)
    adapt(args.zoo.resolve(), args.out.resolve(), args.only)


if __name__ == '__main__':
    main()
