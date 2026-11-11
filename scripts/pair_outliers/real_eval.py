"""Evaluate pair exclusion on real multi-delay runs, as acquired and with injected artifacts.

Reads the runs written by ``real_prep.py``. For each run:

``observed``
    Every method is applied to the run as acquired. There is no ground truth, so this shows
    which pairs each method excludes and how much the fits change.
``injected``
    Artifacts of the kinds in :mod:`simulate` are injected into randomly chosen pairs of the
    real data (one pair, or 15% of the pairs), and every method is applied. Detection is scored
    against the injected pairs, and the fits are compared with the fit of the uninjected run
    with every pair, so the real noise, physiology and model misfit are all present.

Usage::

    python real_eval.py PREP_DIR OUT_DIR [--reps 12] [--workers 20]
"""

import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import ndimage, sparse

sys.path.insert(0, str(Path(__file__).parent))

import evaluate
import simulate

from aslprep.utils.cbf import fit_deltam_multipld
from aslprep.utils.pair_outliers import _delay_groups, detect_outlier_pairs_multipld

#: Fraction of pairs corrupted in the heavier injection condition.
HEAVY_FRACTION = 0.15

#: Model-based methods run on real data (the ablations run in simulation only, to save time).
REAL_METHODS = ('model', 'hybrid')

#: Edge of the cubic parcels, in mm.
PARCEL_MM = 12.0

_RUNS = {}
_REFS = {}


class Run:
    """A prepared real run, restricted to its brain mask."""

    def __init__(self, path):
        data = np.load(path)
        self.name = path.stem
        self.mask = data['mask']
        self.affine = data['affine']
        self.voxel_mm = np.sqrt((self.affine[:3, :3] ** 2).sum(axis=0))
        self.control = data['control']
        self.label = data['label']
        self.m0 = data['m0'][self.mask].astype(float)
        self.pld = data['plds']
        self.tau = data['tau']
        self.readout = str(data['readout'])
        self.multiband = int(data['multiband'])
        slice_time = data['slice_time']
        shift = np.broadcast_to(slice_time[None, None, :], self.mask.shape)[self.mask]
        self.plds = self.pld[None, :] + shift[:, None]
        self.groups = _delay_groups(self.plds, self.tau, True)
        self.n_pairs = self.pld.size
        self.common = {
            'scaled_m0data': self.m0,
            'labeleff': simulate.LABELEFF,
            't1blood': simulate.T1BLOOD,
            'partition_coefficient': simulate.PARTITION,
            'is_casl': True,
        }
        self.smoother = evaluate.make_smoother(self.mask, evaluate.SMOOTH_SIGMA)
        # Parcels: 12 mm cubes of brain voxels, at least half full
        size = np.maximum(np.round(PARCEL_MM / self.voxel_mm).astype(int), 1)
        idx = np.indices(self.mask.shape) // size[:, None, None, None]
        labels = np.ravel_multi_index(tuple(idx), tuple(i.max() + 1 for i in idx))[self.mask]
        ids, counts = np.unique(labels, return_counts=True)
        ids = ids[counts >= 0.5 * np.prod(size)]
        members = (labels[None, :] == ids[:, None]).astype(float)
        self.parcels = sparse.csr_matrix(members / members.sum(axis=1, keepdims=True))

    def deltam(self, control=None, label=None):
        control = self.control if control is None else control
        label = self.label if label is None else label
        return (control - label)[self.mask].astype(float)

    def fit(self, deltam, keep):
        return fit_deltam_multipld(
            deltam_arr=deltam[:, keep], plds=self.plds[:, keep], tau=self.tau[keep], **self.common
        )

    def regional_fit(self, deltam, keep):
        common = dict(self.common, scaled_m0data=self.parcels @ self.m0)
        return fit_deltam_multipld(
            deltam_arr=self.parcels @ deltam[:, keep],
            plds=self.parcels @ self.plds[:, keep],
            tau=self.tau[keep],
            **common,
        )


def _get_run(prep_dir, name):
    if name not in _RUNS:
        # Keep at most two runs in memory per worker.
        while len(_RUNS) >= 2:
            _RUNS.pop(next(iter(_RUNS)))
        _RUNS[name] = Run(prep_dir / f'{name}.npz')
    return _RUNS[name]


def _reference(run):
    """Fits of the uninjected run with every pair, and its mean delta-M per delay."""
    if run.name not in _REFS:
        deltam = run.deltam()
        keep = np.ones(run.n_pairs, dtype=bool)
        mean_dm = np.column_stack(
            [deltam[:, run.groups == g].mean(axis=1) for g in np.unique(run.groups)]
        )
        _REFS[run.name] = (run.regional_fit(deltam, keep), mean_dm)
    return _REFS[run.name]


def _white_noise_sd(deltam):
    """Thermal noise SD of each pair, from differences between neighboring voxels."""
    diffs = np.diff(deltam, axis=0)
    return 1.4826 * np.median(np.abs(diffs - np.median(diffs, axis=0)), axis=0) / np.sqrt(2)


def inject(run, n_bad, rng):
    """Inject artifacts into ``n_bad`` random pairs. Returns delta-M and the injection table."""
    control = run.control.copy()
    label = run.label.copy()
    clean = run.deltam()
    noise_sd = _white_noise_sd(clean)
    bad = np.zeros(run.n_pairs, dtype=bool)
    bad[rng.choice(run.n_pairs, size=n_bad, replace=False)] = True
    kinds = np.full(run.n_pairs, '', dtype=object)
    severity = np.zeros(run.n_pairs)
    extra = np.zeros_like(clean)
    n_slices = run.mask.shape[2]
    for pair in np.flatnonzero(bad):
        kind = simulate.ARTIFACTS[rng.integers(len(simulate.ARTIFACTS))]
        kinds[pair] = kind
        if kind == 'motion':
            sev = rng.uniform(0.3, 1.5)
            label[..., pair] = _rigid_mm(label[..., pair], rng, sev, run.voxel_mm)
        elif kind == 'intensity':
            image = control if rng.random() < 0.5 else label
            if run.readout == '2D':
                sev = rng.uniform(0.02, 0.15)
                # One slice in every multiband group
                n_groups = n_slices // run.multiband
                start = rng.integers(2, n_groups - 2)
                slices = start + n_groups * np.arange(run.multiband)
                image[:, :, slices, pair] *= 1 - sev
            else:
                sev = rng.uniform(0.01, 0.05)
                image[..., pair] *= 1 + sev * rng.choice([-1, 1])
        elif kind == 'labeling':
            sev = rng.uniform(0.3, 1.0)
            members = run.groups == run.groups[pair]
            extra[:, pair] -= sev * clean[:, members].mean(axis=1)
        else:  # spike
            sev = rng.uniform(2.0, 5.0)
            extra[:, pair] += rng.normal(0, np.sqrt(sev**2 - 1) * noise_sd[pair], clean.shape[0])
        severity[pair] = sev
    return run.deltam(control, label) + extra, bad, kinds, severity


def _rigid_mm(volume, rng, magnitude, voxel_mm):
    """Move a volume by a random rigid transform of ``magnitude`` mm and degrees."""
    angles = rng.normal(0, 1, 3)
    angles *= np.deg2rad(magnitude) / np.linalg.norm(angles)
    shift = rng.normal(0, 1, 3)
    shift *= magnitude / np.linalg.norm(shift)
    cx, cy, cz = np.cos(angles)
    sx, sy, sz = np.sin(angles)
    rot = (
        np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
        @ np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
        @ np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    )
    # Rotate in mm about the volume center, then convert to voxels.
    scale = np.diag(voxel_mm)
    matrix = np.linalg.inv(scale) @ rot @ scale
    center = (np.array(volume.shape) - 1) / 2
    offset = center - matrix @ center + shift / voxel_mm
    return ndimage.affine_transform(volume, matrix, offset=offset, order=1, mode='nearest')


def run_methods(run, deltam):
    """Select pairs and fit with every method. Returns keeps, fits, scores and fit counts."""
    keeps = {
        'none': np.ones(run.n_pairs, dtype=bool),
        'zfrac': evaluate.zfrac_keep(deltam, run.groups, np.ones(deltam.shape[0], dtype=bool)),
        'zfrac_loo': evaluate.zfrac_loo_keep(
            deltam, run.groups, np.ones(deltam.shape[0], dtype=bool)
        ),
    }
    scores, n_fits = {}, {}
    for name in REAL_METHODS:
        kwargs = dict(evaluate.METHODS[name])
        kwargs['smooth'] = run.smoother if kwargs['smooth'] else None
        result = detect_outlier_pairs_multipld(
            deltam_arr=deltam, plds=run.plds, tau=run.tau, **run.common, **kwargs
        )
        keeps[name] = result['keep']
        n_fits[name] = result['n_fits']
        if result['scores'] is not None:
            scores[name] = {
                k: np.asarray(v).tolist() for k, v in result['scores'].items() if k != 'noise'
            }
    return keeps, scores, n_fits


def _metrics(run, deltam, keep):
    ref_region, ref_dm = _reference(run)
    region = run.regional_fit(deltam, keep)
    ok = np.isfinite(region[0]) & np.isfinite(ref_region[0])
    scale = float(np.median(ref_region[0][ok]))
    dm_err = []
    for i_group, group in enumerate(np.unique(run.groups)):
        members = (run.groups == group) & keep
        if not members.any():
            # Every pair at this delay was excluded (only the oracle can do this).
            continue
        err = deltam[:, members].mean(axis=1) - ref_dm[:, i_group]
        dm_err.append(np.sqrt(np.mean(err**2)))
    signal = float(np.mean(np.abs(ref_dm)))
    cbf_err = (region[0] - ref_region[0])[ok]
    att_err = (region[1] - ref_region[1])[ok]
    return {
        'n_kept': int(keep.sum()),
        'n_delays': len(dm_err),
        'dm_rmse_pct': 100 * float(np.mean(dm_err)) / signal,
        'region_cbf_rmse_pct': 100 * float(np.sqrt(np.mean(cbf_err**2))) / scale,
        'region_cbf_mae_pct': 100 * float(np.median(np.abs(cbf_err))) / scale,
        'region_att_rmse': float(np.sqrt(np.mean(att_err**2))),
        'region_att_mae': float(np.median(np.abs(att_err))),
        'region_cbf_median': float(np.median(region[0][ok])),
    }


def task_data(run, n_bad, rep):
    """The delta-M of one task: the run as acquired, or with reproducible injected artifacts."""
    if n_bad == 0:
        bad = np.zeros(run.n_pairs, dtype=bool)
        kinds = np.full(run.n_pairs, '', dtype=object)
        return run.deltam(), bad, kinds, np.zeros(run.n_pairs)
    rng = np.random.default_rng(10_000 * rep + 97 * n_bad + len(run.name))
    return inject(run, n_bad, rng)


def run_task(task):
    prep_dir, name, n_bad, rep = task
    start = time.time()
    run = _get_run(Path(prep_dir), name)
    deltam, bad, kinds, severity = task_data(run, n_bad, rep)

    keeps, scores, n_fits = run_methods(run, deltam)
    keeps['oracle'] = ~bad
    metrics = {method: _metrics(run, deltam, keep) for method, keep in keeps.items()}
    return {
        'run': name,
        'n_bad': n_bad,
        'rep': rep,
        'n_pairs': run.n_pairs,
        'n_voxels': int(run.mask.sum()),
        'pld': run.pld.tolist(),
        'bad': bad.tolist(),
        'kind': list(kinds),
        'severity': severity.tolist(),
        'keep': {k: v.tolist() for k, v in keeps.items()},
        'scores': scores,
        'n_fits': n_fits,
        'metrics': metrics,
        'seconds': time.time() - start,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('prep_dir', type=Path)
    parser.add_argument('out_dir', type=Path)
    parser.add_argument('--reps', type=int, default=6)
    parser.add_argument('--workers', type=int, default=20)
    parser.add_argument('--runs', nargs='+', default=None)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    names = args.runs or sorted(p.stem for p in args.prep_dir.glob('*.npz'))

    tasks = []
    for name in names:
        n_pairs = np.load(args.prep_dir / f'{name}.npz')['plds'].size
        heavy = max(1, round(HEAVY_FRACTION * n_pairs))
        conditions = [(0, 0)]
        for n_bad in sorted({1, heavy}):
            conditions += [(n_bad, rep) for rep in range(args.reps)]
        for n_bad, rep in conditions:
            if not (args.out_dir / f'{name}_bad-{n_bad}_rep-{rep}.json').exists():
                tasks.append((str(args.prep_dir), name, n_bad, rep))
    # Group tasks by run, so that each worker loads few runs.
    tasks.sort(key=lambda t: (t[1], t[2], t[3]))
    print(f'{len(tasks)} tasks', flush=True)

    with Pool(args.workers) as pool:
        for result in pool.imap_unordered(run_task, tasks):
            out = f'{result["run"]}_bad-{result["n_bad"]}_rep-{result["rep"]}.json'
            (args.out_dir / out).write_text(json.dumps(result))
            print(out, f'{result["seconds"]:.0f}s', flush=True)


if __name__ == '__main__':
    main()
