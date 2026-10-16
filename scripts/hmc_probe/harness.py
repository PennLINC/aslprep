#!/usr/bin/env python3
"""HMC probe harness (Task 1): reference building, head-motion correction and subtraction only.

One run of one candidate on one ASL series writes, to ``--out``:

- ``aslref.nii.gz``: the reference the candidate built;
- ``hmc_xfms.txt``: one ITK affine per volume of the motion-corrected series, mapping
  the original volume to ``aslref`` (the convention of ASLPrep's ``from-orig_to-aslref`` file);
- ``hmc_aslcontext.tsv`` and ``hmc_indices.json``: which volumes of the raw series those are;
- ``corrected.nii.gz``: the raw volumes resampled once through their transforms onto ``aslref``;
- ``deltam.nii.gz``: control minus label for each pair, from ``corrected.nii.gz``;
- ``harness.json``: what was run (candidate, reference type, processing path, code location).

Candidates:

- ``c0``: production. It builds ASLPrep's own workflows (``init_raw_aslref_wf``, with its
  reference selection, ``RobustAverage``, M0 regridding and GE smoothing; ``ReduceASLFiles``;
  ``init_asl_hmc_wf``), so what it runs is whatever ``aslprep`` package is importable. The
  harness container mounts a pinned revision over the installed package.
- ``compat``: the simplified 2026-10-12 probe (plain-mean reference, raw MCFLIRT per volume
  type, mutual information), kept only to check the harness against ``results.json``.

The harness runs inside the ASLPrep image, which provides FSL and ANTs.

Usage::

    python harness.py --asl <..._asl.nii.gz> --candidate c0 --out <dir> [--m0scan <file>]
"""

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd

HARNESS_VERSION = 1


# ---------------------------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------------------------
def _stem(asl_file):
    name = Path(asl_file).name
    for suffix in ('_asl.nii.gz', '_asl.nii'):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f'{asl_file} is not an *_asl.nii[.gz] file')


def _strip_entity(stem, entity):
    return '_'.join(p for p in stem.split('_') if not p.startswith(f'{entity}-'))


def resolve_inputs(asl_file, aslcontext=None, metadata=None, m0scan=None):
    """Find the sidecar, ``aslcontext.tsv`` and separate M0 of an ASL run.

    The sidecar and context are looked for next to the run (also without its ``echo`` and
    ``part`` entities, as aslscan writes ``part-mag``). A separate M0 is taken from ``m0scan``
    or, if absent and ``M0Type`` is ``Separate``, from ``<stem>_m0scan.nii.gz`` next to the run.
    """
    asl_file = Path(asl_file).resolve()
    stem = _stem(asl_file)
    stems = [stem, _strip_entity(stem, 'echo'), _strip_entity(stem, 'part')]
    folder = asl_file.parent

    def find(suffix):
        for s in stems:
            for cand in (folder / f'{s}_{suffix}', folder.parent / f'{s}_{suffix}'):
                if cand.exists():
                    return cand
        return None

    meta = json.loads(Path(metadata).read_text()) if metadata else {}
    if not metadata:
        side = find('asl.json')
        if side is None:
            raise FileNotFoundError(f'no sidecar for {asl_file}')
        meta = json.loads(side.read_text())
    ctx = Path(aslcontext) if aslcontext else find('aslcontext.tsv')
    if ctx is None:
        raise FileNotFoundError(f'no aslcontext.tsv for {asl_file}')
    if m0scan is None and meta.get('M0Type') == 'Separate':
        m0scan = find('m0scan.nii.gz')
    return {
        'asl_file': str(asl_file),
        'aslcontext': str(ctx.resolve()),
        'metadata': meta,
        'm0scan': str(Path(m0scan).resolve()) if m0scan else None,
    }


def _context(aslcontext):
    return pd.read_table(aslcontext)['volume_type'].tolist()


def _n_volumes(asl_file):
    img = nb.load(asl_file)
    return img.shape[3] if img.ndim == 4 else 1


# ---------------------------------------------------------------------------------------------
# Candidates
# ---------------------------------------------------------------------------------------------
def _collect(out_dir, aslref, xforms, aslcontext, asl_file):
    """Copy a workflow's outputs to fixed names (runs as a Nipype Function node)."""
    import shutil
    from pathlib import Path

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(aslref, out / 'aslref.nii.gz')
    shutil.copyfile(xforms, out / 'hmc_xfms.txt')
    shutil.copyfile(aslcontext, out / 'hmc_aslcontext.tsv')
    return out_dir


def run_c0(inputs, work, smooth_kernel=5.0, force_ge=None):
    """Production reference building and HMC, built from ASLPrep's own workflows."""
    from nipype.interfaces import utility as niu
    from nipype.pipeline import engine as pe

    from aslprep import config
    from aslprep.interfaces.utility import ReduceASLFiles
    from aslprep.utils.asl import select_processing_target
    from aslprep.workflows.asl.fit import select_reference_volume_type
    from aslprep.workflows.asl.hmc import init_asl_hmc_wf
    from aslprep.workflows.asl.reference import init_raw_aslref_wf

    config.workflow.smooth_kernel = smooth_kernel
    metadata = inputs['metadata']
    n_vols = _n_volumes(inputs['asl_file'])
    # aslprep/workflows/asl/base.py: the short-series ("GE") path.
    use_ge = (n_vols <= 5) if force_ge is None else force_ge
    ref_type = select_reference_volume_type(inputs['aslcontext'], metadata, prioritize_m0=use_ge)
    target = select_processing_target(inputs['aslcontext'])

    wf = pe.Workflow(name='harness_c0', base_dir=str(work))
    ref_wf = init_raw_aslref_wf(
        name='hmc_aslref_wf',
        asl_file=inputs['asl_file'],
        reference_volume_type=ref_type,
        m0scan=metadata['M0Type'] == 'Separate',
        use_ge=use_ge,
    )
    ref_wf.inputs.inputnode.aslcontext = inputs['aslcontext']
    ref_wf.inputs.inputnode.m0scan = inputs['m0scan']
    reduce = pe.Node(
        ReduceASLFiles(
            metadata=metadata,
            processing_target=target,
            aslcontext=inputs['aslcontext'],
        ),
        name='reduce_asl_file',
    )
    hmc_wf = init_asl_hmc_wf(name='asl_hmc_wf', mem_gb=1, omp_nthreads=1)
    collected = work / 'collected'
    sink = pe.Node(
        niu.Function(
            function=_collect,
            input_names=['out_dir', 'aslref', 'xforms', 'aslcontext', 'asl_file'],
            output_names=['out_dir'],
        ),
        name='sink',
    )
    sink.inputs.out_dir = str(collected)
    wf.connect([
        (ref_wf, reduce, [('outputnode.asl_file', 'asl_file')]),
        (ref_wf, hmc_wf, [('outputnode.aslref', 'inputnode.raw_ref_image')]),
        (reduce, hmc_wf, [
            ('asl_file', 'inputnode.asl_file'),
            ('aslcontext', 'inputnode.aslcontext'),
        ]),
        (ref_wf, sink, [('outputnode.aslref', 'aslref')]),
        (hmc_wf, sink, [('outputnode.xforms', 'xforms')]),
        (reduce, sink, [('aslcontext', 'aslcontext'), ('asl_file', 'asl_file')]),
    ])  # fmt:skip
    wf.run()
    full = _context(inputs['aslcontext'])
    keep = {'control': ('control', 'label', 'm0scan'), 'deltam': ('deltam', 'm0scan')}.get(
        target, (target, 'm0scan')
    )
    indices = [i for i, t in enumerate(full) if t in keep]
    return {
        'aslref': str(collected / 'aslref.nii.gz'),
        'xforms': str(collected / 'hmc_xfms.txt'),
        'aslcontext': str(collected / 'hmc_aslcontext.tsv'),
        'indices': indices,
        'reference_type': ref_type,
        # RobustAverage registers the selected volumes (all of ``ref_type``) to the first of them
        # (3dvolreg's default base), so the reference depicts that volume's pose.
        'reference_volumes': 'separate_m0'
        if ref_type == 'separate_m0scan'
        else [i for i, t in enumerate(full) if t == ref_type][:1],
        'processing_target': target,
        'use_ge': use_ge,
    }


def _mcflirt_to_itk(mat_files, ref, out_file, work):
    from niworkflows.interfaces.itk import MCFLIRT2ITK

    node = MCFLIRT2ITK(in_files=mat_files, in_reference=ref, in_source=ref, num_threads=1)
    res = node.run(cwd=str(work))
    shutil.copyfile(res.outputs.out_file, out_file)
    return out_file


def run_compat(inputs, work):
    """The 2026-10-12 probe: plain-mean reference, raw MCFLIRT per type, mutual information."""
    asl = nb.load(inputs['asl_file'])
    data = np.asanyarray(asl.dataobj, dtype=np.float32)
    if data.ndim == 3:
        data = data[..., None]
    context = _context(inputs['aslcontext'])
    if inputs['m0scan']:
        m0 = np.asanyarray(nb.load(inputs['m0scan']).dataobj, dtype=np.float32)
        ref = m0.mean(-1) if m0.ndim == 4 else m0
        ref_type = 'separate m0scan'
    else:
        ref_type = next(t for t in ('cbf', 'm0scan', 'deltam', 'control') if t in context)
        ref = data[..., [i for i, t in enumerate(context) if t == ref_type]].mean(-1)
    ref_file = work / 'common_ref.nii.gz'
    nb.Nifti1Image(ref, asl.affine, asl.header).to_filename(ref_file)
    mats = [None] * len(context)
    for vt in dict.fromkeys(context):
        idx = [i for i, t in enumerate(context) if t == vt]
        f = work / f'{vt}.nii.gz'
        nb.Nifti1Image(data[..., idx], asl.affine, asl.header).to_filename(f)
        base = work / f'{vt}_mcf'
        cmd = ['mcflirt', '-in', str(f), '-reffile', str(ref_file), '-out', str(base)]
        cmd += ['-cost', 'mutualinfo', '-mats']
        subprocess.run(cmd, check=True, capture_output=True)
        got = sorted(Path(f'{base}.mat').glob('MAT_*'))
        for i, m in zip(idx, got, strict=True):
            mats[i] = str(m)
    xforms = _mcflirt_to_itk(mats, str(ref_file), work / 'hmc_xfms.txt', work)
    ctx_file = work / 'hmc_aslcontext.tsv'
    pd.DataFrame({'volume_type': context}).to_csv(ctx_file, sep='\t', index=False)
    return {
        'aslref': str(ref_file),
        'xforms': str(xforms),
        'aslcontext': str(ctx_file),
        'indices': list(range(len(context))),
        'reference_type': ref_type,
        'reference_volumes': 'separate_m0'
        if inputs['m0scan']
        else [i for i, t in enumerate(context) if t == ref_type],
        'processing_target': 'all',
        'use_ge': None,
        'fsl_mats': mats,
    }


CANDIDATES = {'c0': run_c0, 'compat': run_compat}


def _register_candidates():
    """Add the candidates of ``candidates.py`` (c1-c4)."""
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from candidates import CANDIDATES as extra

    CANDIDATES.update(extra)


_register_candidates()


# ---------------------------------------------------------------------------------------------
# Resampling and subtraction
# ---------------------------------------------------------------------------------------------
def load_itk_series(xforms):
    """The per-volume 4x4 RAS matrices (orig -> aslref, as points) of an ITK series file."""
    from nitransforms.linear import load

    return np.asarray(load(xforms, fmt='itk').matrix).reshape(-1, 4, 4)


def resample_series(asl_file, indices, matrices, reference, order=3):
    """Resample each listed raw volume once onto ``reference`` through its matrix."""
    from nitransforms.linear import Affine
    from nitransforms.resampling import apply

    asl = nb.load(asl_file)
    ref = nb.load(reference)
    out = np.zeros(ref.shape[:3] + (len(indices),), dtype=np.float32)
    for k, (i, mat) in enumerate(zip(indices, matrices, strict=True)):
        # A 3D file is a one-volume series (index 0), not a stack of slices.
        raw = asl.dataobj if asl.ndim == 3 else asl.dataobj[..., i]
        vol = nb.Nifti1Image(np.asanyarray(raw, dtype=np.float32), asl.affine)
        moved = apply(Affine(mat, reference=ref), vol, reference=ref, order=order)
        out[..., k] = np.asanyarray(moved.dataobj)
    return nb.Nifti1Image(out, ref.affine)


def pair_deltam(corrected, context):
    """Control minus label for each pair, matched in order of appearance."""
    data = np.asanyarray(corrected.dataobj)
    controls = [i for i, t in enumerate(context) if t == 'control']
    labels = [i for i, t in enumerate(context) if t == 'label']
    if len(controls) != len(labels):
        raise ValueError(f'{len(controls)} controls and {len(labels)} labels')
    if not controls:
        return None
    dm = np.stack(
        [data[..., c] - data[..., lab] for c, lab in zip(controls, labels, strict=True)], -1
    )
    return nb.Nifti1Image(dm.astype(np.float32), corrected.affine)


def apply_downstream(inputs, res, work, policy):
    """Replace a candidate's reference by downstream policy D0 or D1 (plan: "Downstream policy").

    - ``d1``: production's own volume selection, each volume resampled through the candidate's
      transform and averaged (plain mean; ``RobustAverage`` would re-register the volumes and
      move the frame). Defined only when every selected volume is in the corrected series;
      otherwise the run falls back to D0 and says so.
    - ``d0``: production's reference exactly, with the candidate's transforms re-expressed in
      its frame. The change of frame ``T`` (production aslref point -> candidate aslref point)
      is the polar mean of ``H_v^-1 H0_v`` over production's selected volumes that both runs
      corrected; with a separate M0 (no shared volume), a rigid FLIRT of the candidate's
      reference to production's.

    Returns the updated result and a record of what was done.
    """
    from candidates import _flirt, _fsl_to_world, _write_itk, polar_mean

    prod = run_c0(inputs, _mkdir(work / 'downstream_c0'))
    sel = prod['reference_volumes']
    h = load_itk_series(res['xforms'])
    idx = list(res['indices'])
    note = {'policy': policy, 'production_reference_volumes': sel}
    if policy == 'd1':
        if sel != 'separate_m0' and all(s in idx for s in sel):
            vols = resample_series(
                inputs['asl_file'], sel, [h[idx.index(s)] for s in sel], res['aslref']
            )
            ref = nb.Nifti1Image(np.asanyarray(vols.dataobj).mean(-1), vols.affine)
            ref.to_filename(work / 'aslref_d1.nii.gz')
            note['d1'] = f'mean of {len(sel)} corrected volume(s)'
            return {
                **res,
                'aslref': str(work / 'aslref_d1.nii.gz'),
                'reference_volumes': sel,
                'reference_frame_kind': 'corrected',
            }, note
        note['d1_fallback_to_d0'] = 'production selects volumes outside the corrected series'
    h0 = load_itk_series(prod['xforms'])
    idx0 = list(prod['indices'])
    shared = [s for s in (sel if sel != 'separate_m0' else []) if s in idx and s in idx0]
    if shared:
        t = polar_mean([np.linalg.inv(h[idx.index(s)]) @ h0[idx0.index(s)] for s in shared])
        note['frame_change'] = f'polar mean over {len(shared)} shared volume(s)'
    else:
        cand_img, prod_img = nb.load(res['aslref']), nb.load(prod['aslref'])
        m = _flirt(res['aslref'], prod['aslref'], work / 'd0_frame.mat', 'mutualinfo')
        world = _fsl_to_world(prod_img) @ m @ np.linalg.inv(_fsl_to_world(cand_img))
        t = np.linalg.inv(world)  # production point -> candidate point
        note['frame_change'] = 'rigid FLIRT (mutual information) of the candidate reference'
    xforms = _write_itk(work / 'hmc_xfms_d0.txt', [m @ t for m in h])
    return {
        **res,
        'aslref': prod['aslref'],
        'xforms': str(xforms),
        'reference_volumes': sel,
        'reference_type': f'{res["reference_type"]}; D0 ({prod["reference_type"]})',
        'reference_frame_kind': 'raw',
    }, note


def _mkdir(path):
    path.mkdir(parents=True, exist_ok=True)
    return path


def run(inputs, candidate, out, downstream='d2', **kwargs):
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out) as tmp:
        work = Path(tmp)
        res = CANDIDATES[candidate](inputs, _mkdir(work / 'candidate'), **kwargs)
        downstream_note = {'policy': downstream}
        if downstream in ('d0', 'd1'):
            res, downstream_note = apply_downstream(inputs, res, work, downstream)
        elif downstream != 'd2':
            raise ValueError(f'Unknown downstream policy: {downstream}')
        shutil.copyfile(res['aslref'], out / 'aslref.nii.gz')
        shutil.copyfile(res['xforms'], out / 'hmc_xfms.txt')
        shutil.copyfile(res['aslcontext'], out / 'hmc_aslcontext.tsv')
        artifacts = res.pop('artifacts', None)
        if artifacts:  # a candidate's audit trail (e.g. C8's regression), kept with the outputs
            shutil.copytree(artifacts, out / 'artifacts', dirs_exist_ok=True)
        fsl_mats = res.pop('fsl_mats', None)
        if fsl_mats:
            mat_dir = out / 'fsl_mats'
            mat_dir.mkdir(exist_ok=True)
            for i, m in enumerate(fsl_mats):
                shutil.copyfile(m, mat_dir / f'MAT_{i:04d}')
    matrices = load_itk_series(out / 'hmc_xfms.txt')
    if len(matrices) != len(res['indices']):
        raise RuntimeError(f'{len(matrices)} transforms for {len(res["indices"])} volumes')
    corrected = resample_series(
        inputs['asl_file'], res['indices'], matrices, out / 'aslref.nii.gz'
    )
    corrected.to_filename(out / 'corrected.nii.gz')
    context = _context(out / 'hmc_aslcontext.tsv')
    dm = pair_deltam(corrected, context)
    if dm is not None:
        dm.to_filename(out / 'deltam.nii.gz')
    (out / 'hmc_indices.json').write_text(json.dumps(res['indices']) + '\n')
    import aslprep

    record = {
        'harness_version': HARNESS_VERSION,
        'candidate': candidate,
        'inputs': {k: v for k, v in inputs.items() if k != 'metadata'},
        'reference_type': res['reference_type'],
        'processing_target': res['processing_target'],
        'use_ge': res['use_ge'],
        'groups': res.get('groups'),
        # Raw volumes the reference was built from, or 'separate_m0' (scoring uses their poses).
        'reference_volumes': res.get('reference_volumes'),
        # 'corrected': the reference is an average of volumes after this run's correction, so it
        # depicts the corrected frame; otherwise it depicts raw volumes (as production's does).
        'reference_frame_kind': res.get('reference_frame_kind', 'raw'),
        'downstream': downstream_note,
        'aslprep': aslprep.__file__,
        'aslprep_version': getattr(aslprep, '__version__', None),
        'options': kwargs,
    }
    (out / 'harness.json').write_text(json.dumps(record, indent=2) + '\n')
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--asl', required=True, help='The *_asl.nii.gz file.')
    parser.add_argument('--candidate', choices=sorted(CANDIDATES), required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--aslcontext')
    parser.add_argument('--metadata', help='JSON sidecar, if not next to the run.')
    parser.add_argument('--m0scan', help='Separate M0 scan, if not next to the run.')
    parser.add_argument('--smooth-kernel', type=float, default=5.0)
    parser.add_argument('--force-ge', choices=('yes', 'no'))
    parser.add_argument(
        '--downstream',
        choices=('d0', 'd1', 'd2'),
        default='d2',
        help="Reference policy: production's (d0), rebuilt from production's selection after "
        "correction (d1), or the candidate's own (d2, the default).",
    )
    parser.add_argument(
        '--opt', action='append', default=[], help='Candidate option, key=value (repeatable).'
    )
    args = parser.parse_args(argv)
    inputs = resolve_inputs(args.asl, args.aslcontext, args.metadata, args.m0scan)
    kwargs = dict(o.split('=', 1) for o in args.opt)
    if args.candidate == 'c0':
        kwargs['smooth_kernel'] = args.smooth_kernel
        if args.force_ge:
            kwargs['force_ge'] = args.force_ge == 'yes'
    record = run(inputs, args.candidate, args.out, downstream=args.downstream, **kwargs)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
