#!/usr/bin/env python3
"""Write PCASL multi-delay probe recipes into ``scripts/hmc_probe/recipes/``.

All derive from ASLPrep's ``pcasl_multipld`` fixture (2D PCASL, six delays from 0.25 to 1.5 s,
background suppression, arterial transit 0.7 s in grey matter and 1.0 s in white matter), the
protocol on which the label offset step's delay cutoff failed. Each varies one thing:

- ``mo_pcaslmulti_{drift,withinpairs,spikes}``: motion (``motion_recipes`` trajectories);
- ``st_pcaslmulti_attshort`` / ``st_pcaslmulti_attlong``: arterial transit 0.5 / 0.7 s and
  1.2 / 1.6 s (grey / white matter), no motion; ``mo_pcaslmulti_attlong_drift`` with drift;
- ``st_pcaslmulti_nobs``: no background suppression, no motion;
- ``st_pcaslmulti_sparse``: three delays (0.5, 1.0, 1.5 s, 14 pairs each), no motion.

Motion recipes are written for the development and held-out seeds (``_ho``).

Usage::

    PCASL_MULTI_BASE=<aslprep>/aslprep/tests/data/aslscan/pcasl_multipld \\
        python scripts/hmc_probe/pcasl_multi_recipes.py
"""

import json
import os
import re
import shutil
import zlib
from pathlib import Path

import numpy as np
from motion_recipes import OUT, REPO, SEEDS, check, trajectory

#: The recipe as on ASLPrep's main branch (after #763, for the current aslscan), where this
#: checkout's copy may be older; PCASL_MULTI_BASE names another ASLPrep checkout's copy.
BASE = Path(
    os.environ.get(
        'PCASL_MULTI_BASE', REPO / 'aslprep' / 'tests' / 'data' / 'aslscan' / 'pcasl_multipld'
    )
)
ATT_LINE = re.compile(r'(?m)^arterial_compartment_att = .*$')


def _att(overlay, gm, wm):
    new, n = ATT_LINE.subn(
        f'arterial_compartment_att = {{ grey_matter = {gm}, white_matter = {wm}, csf = 0.0 }}',
        overlay,
    )
    if n != 1:
        raise ValueError('expected one arterial_compartment_att line')
    return new


def _no_bs(sidecar, overlay):
    sidecar = dict(sidecar, BackgroundSuppression=False)
    for key in ('BackgroundSuppressionNumberPulses', 'BackgroundSuppressionPulseTime'):
        sidecar.pop(key, None)
    lines, out, skip = overlay.splitlines(), [], False
    for line in lines:
        if line.startswith('['):
            skip = line.strip() == '[background_suppression]'
        if not skip:
            out.append(line)
    return sidecar, '\n'.join(out) + '\n'


def _sparse(sidecar):
    plds = [p for p in (0.5, 1.0, 1.5) for _ in range(28)]
    sidecar = dict(sidecar, PostLabelingDelay=plds, TotalAcquiredPairs=len(plds) // 2)
    return sidecar, ['label', 'control'] * (len(plds) // 2)


#: name -> (motion kind or None, edit)
VARIANTS = {
    'mo_pcaslmulti_drift': ('drift', None),
    'mo_pcaslmulti_withinpairs': ('withinpairs', None),
    'mo_pcaslmulti_spikes': ('spikes', None),
    'st_pcaslmulti_attshort': (None, 'attshort'),
    'st_pcaslmulti_attlong': (None, 'attlong'),
    'mo_pcaslmulti_attlong_drift': ('drift', 'attlong'),
    'st_pcaslmulti_nobs': (None, 'nobs'),
    'st_pcaslmulti_sparse': (None, 'sparse'),
}


def write(name, kind, edit, split):
    sidecar = json.loads((BASE / 'asl.json').read_text())
    context = (BASE / 'aslcontext.tsv').read_text().split()[1:]
    overlay = (BASE / 'overlay.toml').read_text()
    recipe = (BASE / 'recipe.toml').read_text()
    if edit == 'attshort':
        overlay = _att(overlay, 0.5, 0.7)
    elif edit == 'attlong':
        overlay = _att(overlay, 1.2, 1.6)
    elif edit == 'nobs':
        sidecar, overlay = _no_bs(sidecar, overlay)
    elif edit == 'sparse':
        sidecar, context = _sparse(sidecar)
    plds = sidecar['PostLabelingDelay']
    full = name + {'heldout': '_ho'}.get(split, '')
    d = OUT / full
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    note = f'Derived from pcasl_multipld ({edit or "unchanged"})'
    if kind:
        seed = SEEDS[split] + zlib.crc32(f'pcaslmulti/{kind}'.encode()) % 1000
        traj = trajectory(kind, context, plds, np.random.default_rng(seed))
        rows = ['trans_x\ttrans_y\ttrans_z\trot_x\trot_y\trot_z']
        rows += ['\t'.join(f'{v:.6f}' for v in [*r[:3], *np.radians(r[3:])]) for r in traj]
        (d / 'motion.tsv').write_text('\n'.join(rows) + '\n')
        overlay = (
            overlay.rstrip('\n') + '\n[motion]\nmode = "trajectory"\ntrajectory = "motion.tsv"\n'
        )
        note += f', trajectory {kind} ({split}, seed {seed})'
    (d / 'asl.json').write_text(json.dumps(sidecar, indent=2) + '\n')
    (d / 'aslcontext.tsv').write_text('volume_type\n' + '\n'.join(context) + '\n')
    (d / 'overlay.toml').write_text(overlay)
    acq = full.replace('_', '')
    lines = [ln for ln in recipe.splitlines() if not ln.startswith('acq')]
    (d / 'recipe.toml').write_text('\n'.join([f'acq = "{acq}"', *lines]) + '\n')
    (d / 'SOURCE.md').write_text(
        f'# {full}\n\n{note}; written by scripts/hmc_probe/pcasl_multi_recipes.py.\n'
    )
    if kind:
        check(full, kind)
    return full


def main():
    names = [
        write(name, kind, edit, split)
        for name, (kind, edit) in VARIANTS.items()
        for split in (('dev', 'heldout') if kind else ('dev',))
    ]
    print('\n'.join(names))


if __name__ == '__main__':
    main()
