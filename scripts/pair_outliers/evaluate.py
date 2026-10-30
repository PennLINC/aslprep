"""Compare pair-exclusion methods on simulated multi-delay runs.

Each task simulates one run (see :mod:`simulate`), selects pairs with every method, and scores
the result against references computed from the noise-free signal:

``dm``
    The mean delta-M of the retained pairs at each delay, against the noise-free delta-M,
    over the brain: RMS error, as a percentage of the mean grey-matter signal, averaged over
    delays. This is what exclusion changes directly.
``voxel``
    The voxelwise general kinetic model fit, in grey matter, against the fit to noise-free
    data (so the model's misfit to the simulated physics is not counted as error):
    median absolute error of CBF and ATT, and the error of the grey-matter median CBF.
``region``
    The model fit to the grey-matter mean of each 12 mm cubic parcel, against the same fit to
    noise-free data: RMS error of CBF and ATT over parcels.

Results go to one JSON file per task, which ``report.py`` reads.

Usage::

    python evaluate.py OUT_DIR [--seeds 10] [--workers 20] [--designs hcpa rep4 rep2 rep1]

Set ``OMP_NUM_THREADS=1`` (and the other BLAS thread variables) so that the workers do not
oversubscribe the CPU.
"""

import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import ndimage, sparse, stats

sys.path.insert(0, str(Path(__file__).parent))

import simulate

from aslprep.utils.cbf import fit_deltam_multipld
from aslprep.utils.pair_outliers import _delay_groups, detect_outlier_pairs_multipld

#: Number of corrupted pairs simulated for each design.
CONTAMINATION = {
    'hcpa': [0, 1, 4, 9],
    'rep4': [0, 1, 2, 5],
    'rep2': [0, 1, 2],
    'rep1': [0, 1],
}

#: Gaussian smoothing of residuals before scoring, in voxels (4 mm), for the smoothed methods.
SMOOTH_SIGMA = 1.0

#: Model-based detection methods: keyword arguments to detect_outlier_pairs_multipld.
#: ``smooth: True`` is replaced by the smoothing function of the phantom grid.
METHODS = {
    'empirical': {'shrinkage': 0.0, 'smooth': True},
    'hybrid': {'shrinkage': 2.0, 'smooth': True},
    'model': {'shrinkage': np.inf, 'smooth': True},
    'hybrid_raw': {'shrinkage': 2.0, 'smooth': False},
    'hybrid_noamp': {'shrinkage': 2.0, 'smooth': True, 'use_amplitude': False},
}

#: Edge of the cubic parcels used for regional fits, in voxels (12 mm).
PARCEL_SIZE = 3

#: Fewest grey-matter voxels in a parcel.
PARCEL_MIN_VOXELS = 10

_PHANTOM = None
_PARCELS = None
_REFERENCE = {}


def _init_worker():
    global _PHANTOM, _PARCELS  # noqa: PLW0603
    _PHANTOM = simulate.Phantom(seed=0)
    idx = np.indices(_PHANTOM.shape) // PARCEL_SIZE
    labels = np.ravel_multi_index(tuple(idx), tuple(i.max() + 1 for i in idx))[_PHANTOM.head]
    gm = _PHANTOM.gm_mask[_PHANTOM.head]
    ids, counts = np.unique(labels[gm], return_counts=True)
    ids = ids[counts >= PARCEL_MIN_VOXELS]
    # (n_parcels, n_head_voxels) averaging matrix over each parcel's GM voxels
    members = (labels[None, :] == ids[:, None]) & gm[None, :]
    _PARCELS = sparse.csr_matrix(members / members.sum(axis=1, keepdims=True))


def make_smoother(mask, sigma):
    """Smooth columns of arrays of the voxels in ``mask`` with a Gaussian kernel."""

    def smooth(values):
        out = np.empty_like(values)
        volume = np.zeros(mask.shape)
        for i_col in range(values.shape[1]):
            volume[mask] = values[:, i_col]
            out[:, i_col] = ndimage.gaussian_filter(volume, sigma, mode='constant')[mask]
        return out

    return smooth


def zfrac_keep(deltam, groups, voxels, z=1.5, fraction=0.2):
    """Model-free comparator: reject pairs that differ from the other repeats at their delay.

    Within each delay with at least three pairs, a pair is rejected if more than ``fraction``
    of the brain voxels lie more than ``z`` standard deviations from that voxel's mean over
    the delay's pairs. At least one pair is kept at each delay.
    This is the rule the ``asl-outlier-eval`` study attributes to Harteveld et al. (2022).
    """
    keep = np.ones(deltam.shape[1], dtype=bool)
    data = deltam[voxels]
    for group in np.unique(groups):
        members = np.flatnonzero(groups == group)
        if members.size < 3:
            continue

        values = data[:, members]
        mean = values.mean(axis=1, keepdims=True)
        sd = values.std(axis=1, ddof=1, keepdims=True)
        far = np.abs(values - mean) > z * np.maximum(sd, np.finfo(float).tiny)
        outlying = far.mean(axis=0) > fraction
        if outlying.all():
            outlying[np.argmin(far.mean(axis=0))] = False
        keep[members[outlying]] = False
    return keep


def zfrac_loo_keep(deltam, groups, voxels, z=1.5, fraction=0.2):
    """Leave-one-out version of :func:`zfrac_keep`, calibrated to the same null exceedance.

    Each pair is compared with the mean and SD of the *other* pairs at its delay, which needs
    at least two others. Under normal noise, (x - mean) / (sd * sqrt(1 + 1/k)) follows a t
    distribution with k - 1 degrees of freedom for k others, so the threshold is the t quantile
    that voxels of a clean pair exceed as often as they exceed ``z`` in a large sample
    (13.4% for ``z = 1.5``). A pair is rejected if more than ``fraction`` of the voxels exceed it.
    """
    exceed = 2 * stats.norm.sf(z)
    keep = np.ones(deltam.shape[1], dtype=bool)
    data = deltam[voxels]
    for group in np.unique(groups):
        members = np.flatnonzero(groups == group)
        if members.size < 3:
            continue

        far_fraction = np.zeros(members.size)
        for i, pair in enumerate(members):
            others = data[:, members[members != pair]]
            k = others.shape[1]
            limit = stats.t.isf(exceed / 2, k - 1) * np.sqrt(1 + 1 / k)
            sd = np.maximum(others.std(axis=1, ddof=1), np.finfo(float).tiny)
            far = np.abs(data[:, pair] - others.mean(axis=1)) > limit * sd
            far_fraction[i] = far.mean()
        outlying = far_fraction > fraction
        if outlying.all():
            outlying[np.argmin(far_fraction)] = False
        keep[members[outlying]] = False
    return keep


def compute_metrics(run, data, keep, fit, groups, brain, gm, common, tau, refs):
    """Score one method's pair selection and fit against the noise-free references."""
    ref_voxel, ref_region = refs
    gm_signal = np.mean(run['truth'][gm])
    dm_err = []
    for group in np.unique(groups):
        members = (groups == group) & keep
        if not members.any():
            # Every pair at this delay was excluded (only the oracle can do this).
            continue
        err = data[brain][:, members].mean(axis=1) - run['truth'][brain][:, groups == group][:, 0]
        dm_err.append(np.sqrt(np.mean(err**2)))

    voxels = gm & np.isfinite(fit[0]) & np.isfinite(ref_voxel[0])
    cbf_err = fit[0][voxels] - ref_voxel[0][voxels]
    att_err = fit[1][voxels] - ref_voxel[1][voxels]
    region = _regional_fit(data[:, keep], run['plds'][:, keep], tau[keep], common)
    ok = np.isfinite(region[0]) & np.isfinite(ref_region[0])
    return {
        'n_kept': int(keep.sum()),
        'n_delays': len(dm_err),
        'dm_rmse_pct': 100 * float(np.mean(dm_err)) / gm_signal,
        'voxel_cbf_mae': float(np.median(np.abs(cbf_err))),
        'voxel_att_mae': float(np.median(np.abs(att_err))),
        'gm_cbf_bias': float(np.median(fit[0][voxels]) - np.median(ref_voxel[0][voxels])),
        'region_cbf_rmse': float(np.sqrt(np.mean((region[0] - ref_region[0])[ok] ** 2))),
        'region_att_rmse': float(np.sqrt(np.mean((region[1] - ref_region[1])[ok] ** 2))),
    }


def _reference(design, run, common, tau):
    """Fits to the noise-free signal: voxelwise and regional. Cached per design."""
    if design.name not in _REFERENCE:
        _, first = np.unique(run['pld_per_pair'], return_index=True)
        truth = run['truth'][:, first]
        plds = run['plds'][:, first]
        voxel = fit_deltam_multipld(deltam_arr=truth, plds=plds, tau=tau[first], **common)
        region = _regional_fit(truth, plds, tau[first], common)
        _REFERENCE[design.name] = (voxel, region)
    return _REFERENCE[design.name]


def _regional_fit(deltam, plds, tau, common):
    region_common = dict(common, scaled_m0data=_PARCELS @ common['scaled_m0data'])
    return fit_deltam_multipld(
        deltam_arr=_PARCELS @ deltam, plds=_PARCELS @ plds, tau=tau, **region_common
    )


def run_task(task):
    design_name, n_bad, seed = task
    design = simulate.DESIGNS[design_name]
    phantom = _PHANTOM
    start = time.time()
    run = simulate.simulate_run(phantom, design, n_bad, seed=1000 * seed + 17 * n_bad + 1)

    head = phantom.head
    brain = phantom.brain[head]
    gm = phantom.gm_mask[head]
    common = {
        'scaled_m0data': phantom.m0[head],
        'labeleff': simulate.LABELEFF,
        't1blood': simulate.T1BLOOD,
        'partition_coefficient': simulate.PARTITION,
        'is_casl': True,
    }
    tau = np.full(run['plds'].shape[1], design.tau)
    groups = _delay_groups(run['plds'], tau, True)
    ref_voxel, ref_region = _reference(design, run, common, tau)

    def plain_fit(data, keep):
        return fit_deltam_multipld(
            deltam_arr=data[:, keep], plds=run['plds'][:, keep], tau=tau[keep], **common
        )

    all_pairs = np.ones(tau.size, dtype=bool)
    keeps = {
        'none': all_pairs,
        'oracle': ~run['bad'],
        'zfrac': zfrac_keep(run['dirty'], groups, brain),
        'zfrac_loo': zfrac_loo_keep(run['dirty'], groups, brain),
    }
    fits = {name: plain_fit(run['dirty'], keep) for name, keep in keeps.items()}
    fits['clean'] = plain_fit(run['clean'], all_pairs)
    keeps['clean'] = all_pairs
    scores = {}
    n_fits = {}
    smoother = make_smoother(head, SMOOTH_SIGMA)
    for name, kwargs in METHODS.items():
        kwargs = dict(kwargs)
        kwargs['smooth'] = smoother if kwargs['smooth'] else None
        result = detect_outlier_pairs_multipld(
            deltam_arr=run['dirty'],
            plds=run['plds'],
            tau=tau,
            mask=brain,
            **common,
            **kwargs,
        )
        fits[name] = result['fit']
        keeps[name] = result['keep']
        n_fits[name] = result['n_fits']
        if result['scores'] is not None:
            scores[name] = {
                key: np.asarray(value).tolist()
                for key, value in result['scores'].items()
                if key != 'noise'
            }

    metrics = {}
    for name, fit in fits.items():
        data = run['clean'] if name == 'clean' else run['dirty']
        metrics[name] = compute_metrics(
            run, data, keeps[name], fit, groups, brain, gm, common, tau,
            (ref_voxel, ref_region),
        )  # fmt: skip

    return {
        'design': design_name,
        'n_bad': n_bad,
        'seed': seed,
        'n_pairs': int(tau.size),
        'pld': run['pld_per_pair'].tolist(),
        'bad': run['bad'].tolist(),
        'kind': list(run['kinds']),
        'severity': run['severity'].tolist(),
        'keep': {name: keep.tolist() for name, keep in keeps.items()},
        'scores': scores,
        'n_fits': n_fits,
        'metrics': metrics,
        'ref_gm_cbf': float(np.median(ref_voxel[0][gm & np.isfinite(ref_voxel[0])])),
        'ref_region_cbf': float(np.median(ref_region[0][np.isfinite(ref_region[0])])),
        'seconds': time.time() - start,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('out_dir', type=Path)
    parser.add_argument('--seeds', type=int, default=10)
    parser.add_argument('--workers', type=int, default=20)
    parser.add_argument('--designs', nargs='+', default=list(CONTAMINATION))
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    tasks = []
    for seed in range(args.seeds):
        for design in args.designs:
            for n_bad in CONTAMINATION[design]:
                out_file = args.out_dir / f'{design}_bad-{n_bad}_seed-{seed}.json'
                if not out_file.exists():
                    tasks.append((design, n_bad, seed))
    print(f'{len(tasks)} tasks', flush=True)

    with Pool(args.workers, initializer=_init_worker) as pool:
        for result in pool.imap_unordered(run_task, tasks):
            name = f'{result["design"]}_bad-{result["n_bad"]}_seed-{result["seed"]}.json'
            (args.out_dir / name).write_text(json.dumps(result))
            print(name, f'{result["seconds"]:.0f}s', flush=True)


if __name__ == '__main__':
    main()
