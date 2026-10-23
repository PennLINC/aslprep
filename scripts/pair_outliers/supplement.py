"""Update finished results without rerunning every method.

Simulated runs and injected artifacts are reproducible from their seeds, so each task is
regenerated and checked against its stored no-exclusion result before anything is added.

``sim`` and ``real``
    Add the leave-one-out model-free comparator (``zfrac_loo``, see
    :func:`evaluate.zfrac_loo_keep`), which was added after the runs started, in response to
    review.
``redo``
    Recompute the methods that the three-other-pairs rule for the same-delay correction
    changed (:data:`CHANGED`) for first-run results, writing complete results to ``--out-dir``.

Usage::

    python supplement.py redo FIRST_RUN_DIR --out-dir RESULTS_DIR [--workers 20]
    python supplement.py sim RESULTS_DIR [--workers 20]
    python supplement.py real RESULTS_DIR --prep-dir PREP_DIR [--workers 10]
"""

import argparse
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

import evaluate
import real_eval
import simulate

from aslprep.utils.cbf import fit_deltam_multipld
from aslprep.utils.pair_outliers import _delay_groups

METHOD = 'zfrac_loo'


def _check(stored, recomputed, path):
    if not np.isclose(stored, recomputed, rtol=1e-6):
        raise RuntimeError(f'{path}: regenerated data differ ({stored} vs {recomputed})')


def supplement_sim(path):
    path = Path(path)
    result = json.loads(path.read_text())
    if METHOD in result['keep']:
        return path.name, 'done'
    phantom = evaluate._PHANTOM
    design = simulate.DESIGNS[result['design']]
    n_bad, seed = result['n_bad'], result['seed']
    run = simulate.simulate_run(phantom, design, n_bad, seed=1000 * seed + 17 * n_bad + 1)
    if run['bad'].tolist() != result['bad']:
        raise RuntimeError(f'{path}: regenerated corruption differs')

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
    refs = evaluate._reference(design, run, common, tau)

    def metrics_for(keep):
        fit = fit_deltam_multipld(
            deltam_arr=run['dirty'][:, keep], plds=run['plds'][:, keep], tau=tau[keep], **common
        )
        return evaluate.compute_metrics(
            run, run['dirty'], keep, fit, groups, brain, gm, common, tau, refs
        )

    all_pairs = np.ones(tau.size, dtype=bool)
    # The no-exclusion fit is deterministic, so it must reproduce the stored result.
    _check(result['metrics']['none']['dm_rmse_pct'], metrics_for(all_pairs)['dm_rmse_pct'], path)
    keep = evaluate.zfrac_loo_keep(run['dirty'], groups, brain)
    result['keep'][METHOD] = keep.tolist()
    result['metrics'][METHOD] = (
        dict(result['metrics']['none'], n_kept=int(keep.sum()))
        if keep.all()
        else metrics_for(keep)
    )
    path.write_text(json.dumps(result))
    return path.name, int((~keep).sum())


#: Methods changed by the rule that the same-delay correction needs three other pairs.
#: The others (no exclusion, oracle, the z rules and the model alone) are unaffected by it.
CHANGED = ('empirical', 'hybrid', 'hybrid_raw', 'hybrid_noamp')


def redo_sim(args):
    """Recompute the methods in CHANGED for one first-run result, keeping the others."""
    src, out_dir = Path(args[0]), Path(args[1])
    out = out_dir / src.name
    if out.exists():
        return src.name, 'exists'
    result = json.loads(src.read_text())
    phantom = evaluate._PHANTOM
    design = simulate.DESIGNS[result['design']]
    n_bad, seed = result['n_bad'], result['seed']
    run = simulate.simulate_run(phantom, design, n_bad, seed=1000 * seed + 17 * n_bad + 1)
    if run['bad'].tolist() != result['bad']:
        raise RuntimeError(f'{src}: regenerated corruption differs')

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
    refs = evaluate._reference(design, run, common, tau)
    all_pairs = np.ones(tau.size, dtype=bool)
    fit = fit_deltam_multipld(deltam_arr=run['dirty'], plds=run['plds'], tau=tau, **common)
    none = evaluate.compute_metrics(
        run, run['dirty'], all_pairs, fit, groups, brain, gm, common, tau, refs
    )
    _check(result['metrics']['none']['dm_rmse_pct'], none['dm_rmse_pct'], src)
    _check(result['metrics']['none']['region_cbf_rmse'], none['region_cbf_rmse'], src)

    smoother = evaluate.make_smoother(head, evaluate.SMOOTH_SIGMA)
    for name in CHANGED:
        kwargs = dict(evaluate.METHODS[name])
        kwargs['smooth'] = smoother if kwargs['smooth'] else None
        detected = evaluate.detect_outlier_pairs_multipld(
            deltam_arr=run['dirty'], plds=run['plds'], tau=tau, mask=brain, **common, **kwargs
        )
        result['keep'][name] = detected['keep'].tolist()
        result['n_fits'][name] = detected['n_fits']
        result['metrics'][name] = evaluate.compute_metrics(
            run, run['dirty'], detected['keep'], detected['fit'], groups, brain, gm, common,
            tau, refs,
        )  # fmt: skip
        if detected['scores'] is None:
            result['scores'].pop(name, None)
        else:
            result['scores'][name] = {
                k: np.asarray(v).tolist() for k, v in detected['scores'].items() if k != 'noise'
            }
    result['redone'] = list(CHANGED)
    out.write_text(json.dumps(result))
    return src.name, 'redone'


def supplement_real(args):
    prep_dir, path = args
    path = Path(path)
    result = json.loads(path.read_text())
    if METHOD in result['keep']:
        return path.name, 'done'
    run = real_eval._get_run(Path(prep_dir), result['run'])
    deltam, bad, _, _ = real_eval.task_data(run, result['n_bad'], result['rep'])
    if bad.tolist() != result['bad']:
        raise RuntimeError(f'{path}: regenerated injection differs')
    all_pairs = np.ones(run.n_pairs, dtype=bool)
    _check(
        result['metrics']['none']['dm_rmse_pct'],
        real_eval._metrics(run, deltam, all_pairs)['dm_rmse_pct'],
        path,
    )
    keep = evaluate.zfrac_loo_keep(deltam, run.groups, np.ones(deltam.shape[0], dtype=bool))
    result['keep'][METHOD] = keep.tolist()
    result['metrics'][METHOD] = real_eval._metrics(run, deltam, keep)
    path.write_text(json.dumps(result))
    return path.name, int((~keep).sum())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('kind', choices=['sim', 'real', 'redo'])
    parser.add_argument('results_dir', type=Path)
    parser.add_argument('--prep-dir', type=Path)
    parser.add_argument('--out-dir', type=Path, help='for redo: where the redone results go')
    parser.add_argument('--workers', type=int, default=10)
    args = parser.parse_args()
    paths = sorted(args.results_dir.glob('*.json'))
    if args.kind == 'redo':
        tasks = [(str(p), str(args.out_dir)) for p in paths]
        with Pool(args.workers, initializer=evaluate._init_worker) as pool:
            for name, status in pool.imap_unordered(redo_sim, tasks):
                print(name, status, flush=True)
    elif args.kind == 'sim':
        with Pool(args.workers, initializer=evaluate._init_worker) as pool:
            for name, n in pool.imap_unordered(supplement_sim, paths):
                print(name, n, flush=True)
    else:
        # Group by run so that each worker loads few runs.
        tasks = [(str(args.prep_dir), str(p)) for p in paths]
        with Pool(args.workers) as pool:
            for name, n in pool.imap_unordered(supplement_real, tasks, chunksize=4):
                print(name, n, flush=True)


if __name__ == '__main__':
    main()
