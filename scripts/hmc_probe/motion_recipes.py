#!/usr/bin/env python3
"""Write the probe's motion recipes (Task 4) into ``scripts/hmc_probe/recipes/``.

Each recipe is an existing protocol with a known rigid trajectory (aslscan's ``trajectory``
motion mode, translations in mm and rotations in radians, one row per volume). Trajectories
with random parts are drawn from ``SEEDS['dev']`` for development and ``SEEDS['heldout']`` for
the held-out set (recipe names ending ``_ho``); nothing else differs between the two.

Trajectory kinds:

- ``drift``: a slow drift to about 1.5 mm and 1.5 degrees;
- ``spikes``: still, with two 2 mm / 2 degree spikes;
- ``throughplane``: drift along z with rotation about x (through the slices);
- ``betweenpairs``: a random walk that changes only between control-label pairs;
- ``withinpairs``: as ``betweenpairs``, plus a 0.5 mm / 0.5 degree jump inside each pair;
- ``synchronous``: every label volume displaced by the same 0.4 mm / 0.4 degree (motion
  locked to label/control, which zig-zag regression would remove);
- ``pldconfounded``: the volumes of one delay displaced by 1 mm / 1 degree (multi-delay only);
- ``constant``: every volume at the same pose (the constant-pose control).

Protocols: ``pcasl`` (2D PCASL single delay, separate M0, from ``pcasl1pld``), ``pcaslinc``
(the same with two M0 volumes in the series, so the M0 moves with the head), ``pasl`` (2D PASL,
8 delays x 1 pair, separate M0, from ``pasl_multipld``), ``paslrep`` (2D PASL, 4 delays x 3
pairs), ``pcaslshort2`` and ``pcaslshort4`` (2 and 4 pairs), ``deltamrep`` (an in-series M0
and 8 delta-M volumes).

Usage::

    python scripts/hmc_probe/motion_recipes.py
"""

import json
import shutil
import zlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OUT = HERE / 'recipes'
SEEDS = {'dev': 20261015, 'heldout': 20261016}


# ---------------------------------------------------------------------------------------------
# Protocols
# ---------------------------------------------------------------------------------------------
def _base(name):
    d = REPO / 'aslprep' / 'tests' / 'data' / 'aslscan' / name
    sidecar = json.loads((d / 'asl.json').read_text())
    context = (d / 'aslcontext.tsv').read_text().split()[1:]
    overlay = (d / 'overlay.toml').read_text()
    recipe = (d / 'recipe.toml').read_text()
    return sidecar, context, overlay, recipe


def protocol(kind):
    """``(sidecar, context, overlay, recipe text, plds)`` for a protocol."""
    if kind in ('pcasl', 'pcaslinc', 'pcaslshort2', 'pcaslshort4', 'deltamrep'):
        sidecar, context, overlay, recipe = _base('pcasl1pld')
        if kind == 'pcaslinc':
            context = ['m0scan', 'm0scan', *context]
            sidecar['M0Type'] = 'Included'
            overlay = overlay.split('[m0]')[0]
        elif kind.startswith('pcaslshort'):
            n = int(kind[-1])
            context = context[: 2 * n]
            sidecar['TotalAcquiredPairs'] = n
        elif kind == 'deltamrep':
            context = ['m0scan', *(['deltam'] * 8)]
            sidecar['M0Type'] = 'Included'
            sidecar.pop('TotalAcquiredPairs', None)
            overlay = overlay.split('[m0]')[0]
        plds = [sidecar['PostLabelingDelay']] * len(context)
        return sidecar, context, overlay, recipe, plds
    if kind in ('pasl', 'paslrep'):
        sidecar, context, overlay, recipe = _base('pasl_multipld')
        if kind == 'paslrep':
            delays = [0.9, 1.5, 2.1, 2.7]
            context, plds = [], []
            for d in delays:
                context += ['label', 'control'] * 3
                plds += [d] * 6
            sidecar['PostLabelingDelay'] = plds
            sidecar['TotalAcquiredPairs'] = len(context) // 2
        plds = sidecar['PostLabelingDelay']
        return sidecar, context, overlay, recipe, plds
    raise ValueError(kind)


# ---------------------------------------------------------------------------------------------
# Trajectories: (n, 6) arrays of [tx ty tz mm, rx ry rz degrees]
# ---------------------------------------------------------------------------------------------
def _pairs(context):
    """Pair index of each volume (control and label in order of appearance; others -1)."""
    pair, out, seen = -1, [], {'control': 0, 'label': 0}
    for t in context:
        if t in seen:
            seen[t] += 1
            pair = max(seen['control'], seen['label']) - 1
            out.append(pair)
        else:
            out.append(-1)
    return out


def _jump(rng, trans_mm, rot_deg):
    """A random rigid jump: translation of norm ``trans_mm`` and rotation of norm ``rot_deg``.

    Translation and rotation directions are drawn and scaled independently, so each amplitude
    is the vector norm of its part (a 2 mm / 2 degree spike is 2 mm and 2 degrees, not a
    rotation scaled by the translation's norm).
    """
    t = rng.normal(size=3)
    r = rng.normal(size=3)
    return np.r_[trans_mm * t / np.linalg.norm(t), rot_deg * r / np.linalg.norm(r)]


def trajectory(kind, context, plds, rng):
    n = len(context)
    t = np.zeros((n, 6))
    ramp = np.linspace(0, 1, n)[:, None]
    if kind == 'drift':
        t = ramp * np.array([0.8, -0.6, 1.0, 0.8, -0.5, 0.7]) * rng.uniform(0.8, 1.2)
    elif kind == 'spikes':
        for v in sorted(rng.choice(np.arange(2, n - 1), size=min(2, n - 3), replace=False)):
            t[v] = _jump(rng, 2.0, 2.0)
    elif kind == 'throughplane':
        t[:, 2] = ramp[:, 0] * 2.0 * rng.uniform(0.8, 1.2)
        t[:, 3] = ramp[:, 0] * 2.0 * rng.uniform(0.8, 1.2)
    elif kind in ('betweenpairs', 'withinpairs'):
        pairs = _pairs(context)
        steps = rng.normal(scale=0.25, size=(max(pairs) + 1, 6))
        walk = np.cumsum(steps, axis=0)
        for v, p in enumerate(pairs):
            t[v] = walk[max(p, 0)]
        if kind == 'withinpairs':
            for v, ty in enumerate(context):
                if ty == 'label':
                    t[v] += _jump(rng, 0.5, 0.5)
    elif kind == 'synchronous':
        offset = np.array([0.4, 0.0, 0.0, 0.0, 0.0, 0.4]) * rng.choice([-1, 1])
        t[[v for v, ty in enumerate(context) if ty == 'label']] = offset
    elif kind == 'pldconfounded':
        unique = sorted(set(plds))
        target = unique[len(unique) // 2]
        t[[v for v, p in enumerate(plds) if p == target]] = np.array([1.0, 0, 0, 0, 0, 1.0])
    elif kind == 'constant':
        t[:] = np.array([1.0, -0.8, 0.6, 1.0, -1.0, 0.8])
    else:
        raise ValueError(kind)
    return t


# ---------------------------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------------------------
#: (protocol, trajectory kind) pairs to build.
PLAN = [
    *[('pcasl', k) for k in ('drift', 'spikes', 'throughplane', 'betweenpairs', 'withinpairs',
                             'synchronous', 'constant')],
    ('pcaslinc', 'constant'),
    ('pcaslinc', 'drift'),
    *[('pasl', k) for k in ('drift', 'betweenpairs', 'withinpairs', 'pldconfounded')],
    *[('paslrep', k) for k in ('drift', 'pldconfounded', 'withinpairs')],
    ('pcaslshort2', 'drift'),
    ('pcaslshort4', 'drift'),
    ('deltamrep', 'drift'),
    ('deltamrep', 'spikes'),
]  # fmt: skip


def write(protocol_name, kind, split):
    sidecar, context, overlay, recipe, plds = protocol(protocol_name)
    # crc32, not hash(): Python randomizes string hashes per process.
    seed = SEEDS[split] + zlib.crc32(f'{protocol_name}/{kind}'.encode()) % 1000
    rng = np.random.default_rng(seed)
    traj = trajectory(kind, context, plds, rng)
    name = f'mo_{protocol_name}_{kind}' + ('_ho' if split == 'heldout' else '')
    acq = name.replace('_', '')
    d = OUT / name
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    (d / 'asl.json').write_text(json.dumps(sidecar, indent=2) + '\n')
    (d / 'aslcontext.tsv').write_text('volume_type\n' + '\n'.join(context) + '\n')
    rows = ['trans_x\ttrans_y\ttrans_z\trot_x\trot_y\trot_z']
    rows += ['\t'.join(f'{v:.6f}' for v in [*r[:3], *np.radians(r[3:])]) for r in traj]
    (d / 'motion.tsv').write_text('\n'.join(rows) + '\n')
    overlay = overlay.rstrip('\n') + '\n[motion]\nmode = "trajectory"\ntrajectory = "motion.tsv"\n'
    (d / 'overlay.toml').write_text(overlay)
    lines = [ln for ln in recipe.splitlines() if not ln.startswith('acq')]
    (d / 'recipe.toml').write_text('\n'.join([f'acq = "{acq}"', *lines]) + '\n')
    (d / 'SOURCE.md').write_text(
        f'# {name}\n\nProtocol {protocol_name}, trajectory {kind} ({split}, seed {seed}); '
        'written by scripts/hmc_probe/motion_recipes.py.\n'
    )
    return name


#: Largest translation (mm) and rotation (degrees) per kind, as vector norms of each volume's
#: pose; checked for every generated trajectory.
BOUNDS = {
    'drift': (2.0, 1.6),
    'spikes': (2.0, 2.0),
    'throughplane': (2.5, 2.5),
    'betweenpairs': (None, None),  # a random walk; no fixed bound
    'withinpairs': (None, None),
    'synchronous': (0.4, 0.4),
    'pldconfounded': (1.0, 1.0),
    'constant': (1.5, 1.7),  # the fixed pose (1, -0.8, 0.6 mm; 1, -1, 0.8 deg): 1.41 mm, 1.62 deg
}


def check(name, kind):
    rows = np.loadtxt(OUT / name / 'motion.tsv', skiprows=1, ndmin=2)
    trans = np.linalg.norm(rows[:, :3], axis=1).max()
    rot = np.degrees(np.linalg.norm(rows[:, 3:], axis=1)).max()
    t_max, r_max = BOUNDS[kind]
    # motion.tsv keeps six decimals (radians for rotations), so allow 1e-3 of rounding.
    if t_max is not None and (trans > t_max + 1e-3 or rot > r_max + 1e-3):
        raise ValueError(f'{name}: {trans:.2f} mm / {rot:.2f} deg exceeds {t_max} / {r_max}')
    return trans, rot


def main():
    names = []
    for split in ('dev', 'heldout'):
        for protocol_name, kind in PLAN:
            name = write(protocol_name, kind, split)
            trans, rot = check(name, kind)
            names.append(f'{name}\t{trans:.2f} mm\t{rot:.2f} deg')
    print('\n'.join(names))


if __name__ == '__main__':
    main()
