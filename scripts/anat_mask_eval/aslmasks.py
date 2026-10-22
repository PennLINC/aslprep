#!/usr/bin/env python3
"""Compute ASLPrep's own SynthStrip and legacy brain masks for a set of ASL references.

Runs inside the ASLPrep image (needs ``--skull-strip-asl``, i.e. PennLINC/aslprep#753), on the
jobs written by ``jobs.py``. Each reference is a final coregistration reference
(``desc-coreg_aslref``: bias-corrected, and distortion-corrected when the run had a fieldmap), so
both methods run as ASLPrep runs them on such a reference:

- ``synthstrip``: :func:`~aslprep.workflows.asl.reference.init_synthstrip_aslref_wf` without N4
  (what ``init_skullstrip_aslref_wf`` uses);
- ``legacy``: niworkflows' ``init_skullstrip_bold_wf`` (BET and ``3dAutomask``) for references
  that were distortion-corrected, as ASLPrep does; otherwise
  :func:`~aslprep.workflows.asl.reference.init_enhance_and_skullstrip_asl_wf` without N4 (the
  reference is already bias-corrected).

Usage (in the container)::

    python aslmasks.py --jobs jobs.json --out /out --work /work --nprocs 8
"""

import argparse
import json
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def _workflows(sdc):
    from niworkflows.func.util import init_skullstrip_bold_wf

    from aslprep.workflows.asl.reference import (
        init_enhance_and_skullstrip_asl_wf,
        init_synthstrip_aslref_wf,
    )

    legacy = (
        init_skullstrip_bold_wf(name='legacy')
        if sdc
        else init_enhance_and_skullstrip_asl_wf(disable_n4=True, name='legacy')
    )
    return {
        'synthstrip': init_synthstrip_aslref_wf(disable_n4=True, name='synthstrip'),
        'legacy': legacy,
    }


def _copy(in_file, out_file):
    import shutil

    shutil.copyfile(in_file, out_file)
    return out_file


def run_job(job, out_root, work_root):
    from nipype.interfaces import utility as niu
    from nipype.pipeline import engine as pe

    from aslprep import config

    config.nipype.omp_nthreads = 1
    out = Path(out_root) / job['name']
    out.mkdir(parents=True, exist_ok=True)
    done = []
    for method, wf in _workflows(job['sdc']).items():
        target = out / f'{method}.nii.gz'
        if target.exists():
            continue
        # Identity output nodes are dropped from the executed graph: copy the mask out instead.
        outer = pe.Workflow(name=f'{method}_run', base_dir=str(Path(work_root) / job['name']))
        copy = pe.Node(
            niu.Function(function=_copy, input_names=['in_file', 'out_file']), name='copy'
        )
        copy.inputs.out_file = str(target)
        wf.inputs.inputnode.in_file = job['ref']
        outer.connect(wf, 'outputnode.mask_file', copy, 'in_file')
        outer.run()
        done.append(method)
    shutil.rmtree(Path(work_root) / job['name'], ignore_errors=True)
    return job['name'], done


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--jobs', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--nprocs', type=int, default=4)
    args = parser.parse_args(argv)
    jobs = json.loads(args.jobs.read_text())
    with ProcessPoolExecutor(args.nprocs) as pool:
        futures = {pool.submit(run_job, j, args.out, args.work): j['name'] for j in jobs}
        for fut in as_completed(futures):
            try:
                name, done = fut.result()
                print(f'{name}: {done or "already done"}', flush=True)
            except Exception as err:  # noqa: BLE001 - one failed run must not stop the rest
                print(f'{futures[fut]}: FAILED {type(err).__name__}: {err}', flush=True)


if __name__ == '__main__':
    main()
