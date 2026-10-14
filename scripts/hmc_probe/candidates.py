"""HMC probe candidates beyond production (plan: "Candidates").

Each candidate is a function ``(inputs, work, **options) -> dict`` with the same result keys as
``harness.run_c0``, registered in ``harness.CANDIDATES``. Options come from the harness command
line as ``--opt key=value``.

- ``c1``: own-type templates. Volumes are grouped by type (and by delay where a group would
  have at least ``min_group`` volumes); each group is registered to its own template (the
  group mean, or with ``loo=1`` the mean of the other volumes); each template is then
  registered rigidly to one common frame (``frame=control`` or ``frame=m0``), unless
  ``cross=identity``. Derived volume types (delta-M, CBF) are never registered to another
  contrast (policy P1): a single derived volume gets the identity, a repeated series is
  registered within its own type only.
- ``c2``: as ``c1``, with each template rebuilt from the registered volumes ``iterations``
  times (3 by default), as QSIPrep builds its b=0 template.
- ``c3``: production's transforms, re-expressed relative to one acquired volume (``ref_volume``,
  by default the first control), whose raw image becomes the reference.
- ``c4``: one transform per control-label pair (strict pair constraint): each pair gets the
  polar mean of its two transforms from ``base`` (``c0`` or ``c1``).

All registrations use FSL (MCFLIRT for series, FLIRT for single images) with ``cost``
(``mutualinfo`` by default). Transforms are composed as FSL matrices on the raw grid and
converted to ITK once.
"""

import shutil
import subprocess
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd

#: Volume types never registered to another contrast (policy P1).
DERIVED = ('deltam', 'cbf')
#: What production keeps for HMC, by processing target (aslprep/interfaces/utility.py).
KEPT = {'control': ('control', 'label', 'm0scan'), 'deltam': ('deltam', 'm0scan')}


def _opt(options, key, default, cast=str):
    value = options.get(key, default)
    return cast(value) if value is not None else None


def _kept(inputs):
    from aslprep.utils.asl import select_processing_target

    target = select_processing_target(inputs['aslcontext'])
    context = pd.read_table(inputs['aslcontext'])['volume_type'].tolist()
    keep = KEPT.get(target, (target, 'm0scan'))
    idx = [i for i, t in enumerate(context) if t in keep]
    return target, idx, [context[i] for i in idx]


def _plds(inputs, n_total, idx):
    pld = np.atleast_1d(np.asarray(inputs['metadata'].get('PostLabelingDelay', 0.0), float))
    pld = np.repeat(pld, n_total) if pld.size == 1 else pld
    return [float(pld[i]) for i in idx]


def _save(data, img, path):
    nb.Nifti1Image(np.asarray(data, np.float32), img.affine, img.header).to_filename(path)
    return path


def _mcflirt(series, ref, out_base, cost):
    """MCFLIRT a 4D file to a reference; returns FSL matrices (input -> reference)."""
    cmd = ['mcflirt', '-in', str(series), '-reffile', str(ref), '-out', str(out_base)]
    cmd += ['-cost', cost, '-mats', '-spline_final']
    subprocess.run(cmd, check=True, capture_output=True)
    mats = sorted(Path(f'{out_base}.mat').glob('MAT_*'))
    return [np.loadtxt(m) for m in mats], Path(f'{out_base}.nii.gz')


def _flirt(moving, ref, out_mat, cost):
    """Rigid FLIRT of one image to another; returns the FSL matrix (moving -> reference)."""
    cmd = ['flirt', '-in', str(moving), '-ref', str(ref), '-omat', str(out_mat), '-dof', '6']
    cmd += ['-cost', cost, '-searchrx', '-10', '10', '-searchry', '-10', '10']
    cmd += ['-searchrz', '-10', '10']
    subprocess.run(cmd, check=True, capture_output=True)
    return np.loadtxt(out_mat)


def _to_itk(mats, ref_file, work):
    from harness import _mcflirt_to_itk

    files = []
    for k, m in enumerate(mats):
        f = work / f'final_{k:04d}.mat'
        np.savetxt(f, m)
        files.append(str(f))
    return _mcflirt_to_itk(files, str(ref_file), work / 'hmc_xfms.txt', work)


def _regrid(img_file, ref_img, out):
    """Put an image (e.g. a separate M0) on the raw ASL grid, linear interpolation."""
    from nibabel.processing import resample_from_to

    img = nb.load(img_file)
    data = np.asanyarray(img.dataobj, np.float32)
    if data.ndim == 4:
        data = data.mean(-1)
    img3 = nb.Nifti1Image(data, img.affine)
    target = nb.Nifti1Image(np.zeros(ref_img.shape[:3], np.float32), ref_img.affine)
    same = img3.shape == target.shape and np.allclose(img3.affine, target.affine, atol=1e-4)
    out_img = img3 if same else resample_from_to(img3, target, order=1)
    nb.Nifti1Image(np.asanyarray(out_img.dataobj, np.float32), ref_img.affine).to_filename(out)
    return out


def _groups(types, plds, by_pld, min_group):
    groups = {}
    for k, (t, p) in enumerate(zip(types, plds, strict=True)):
        groups.setdefault((t, p), []).append(k)
    if not by_pld:
        merged = {}
        for (t, _), ks in groups.items():
            merged.setdefault((t, None), []).extend(ks)
        return {g: sorted(ks) for g, ks in merged.items()}
    out = {}
    for (t, p), ks in groups.items():
        key = (t, p) if len(ks) >= min_group else (t, None)
        out.setdefault(key, []).extend(ks)
    return {g: sorted(ks) for g, ks in out.items()}


def run_c1(inputs, work, **options):
    """Own-type templates (c1) and iterated templates (c2, ``iterations`` > 1)."""
    cost = _opt(options, 'cost', 'mutualinfo')
    by_pld = _opt(options, 'grouping', 'type') == 'type_pld'
    min_group = _opt(options, 'min_group', 3, int)
    frame = _opt(options, 'frame', 'control')
    cross = _opt(options, 'cross', 'flirt')
    iterations = _opt(options, 'iterations', 1, int)
    loo = _opt(options, 'loo', 0, int) == 1

    target, idx, types = _kept(inputs)
    img = nb.load(inputs['asl_file'])
    data = np.asanyarray(img.dataobj, np.float32)
    data = data[..., None] if data.ndim == 3 else data
    plds = _plds(inputs, data.shape[-1], idx)
    groups = _groups(types, plds, by_pld, min_group)

    mats = [None] * len(idx)
    templates = {}
    for g, ks in groups.items():
        tag = f'{g[0]}_{g[1] if g[1] is not None else "all"}'.replace('.', 'p').replace('/', '')
        vols = data[..., [idx[k] for k in ks]]
        if len(ks) == 1:
            templates[g] = _save(vols[..., 0], img, work / f'tpl_{tag}.nii.gz')
            mats[ks[0]] = np.eye(4)
            continue
        series = _save(vols, img, work / f'series_{tag}.nii.gz')
        if loo:
            for j, k in enumerate(ks):
                others = np.delete(vols, j, axis=-1).mean(-1)
                tpl = _save(others, img, work / f'loo_{tag}_{j}.nii.gz')
                one = _save(vols[..., j], img, work / f'one_{tag}_{j}.nii.gz')
                mats[k] = _flirt(one, tpl, work / f'loo_{tag}_{j}.mat', cost)
            templates[g] = _save(vols.mean(-1), img, work / f'tpl_{tag}.nii.gz')
            continue
        tpl = _save(vols.mean(-1), img, work / f'tpl_{tag}_it0.nii.gz')
        for it in range(iterations):
            got, moved = _mcflirt(series, tpl, work / f'mcf_{tag}_it{it}', cost)
            if it < iterations - 1:
                moved_data = np.asanyarray(nb.load(moved).dataobj, np.float32)
                tpl = _save(moved_data.mean(-1), img, work / f'tpl_{tag}_it{it + 1}.nii.gz')
        for k, m in zip(ks, got, strict=True):
            mats[k] = m
        templates[g] = tpl

    # The common frame.
    controls = [g for g in templates if g[0] == 'control']
    m0_groups = [g for g in templates if g[0] == 'm0scan']
    if frame == 'm0' and inputs.get('m0scan'):
        frame_file = _regrid(inputs['m0scan'], img, work / 'frame_m0.nii.gz')
    elif frame == 'm0' and m0_groups:
        frame_file = templates[m0_groups[0]]
    elif controls:
        # With several control groups (one per delay), the longest-delay group's template.
        frame_file = templates[max(controls, key=lambda g: (g[1] is None, g[1] or 0))]
    else:
        frame_file = templates[next(iter(templates))]

    cross_mats = {}
    for g, tpl in templates.items():
        if Path(tpl) == Path(frame_file) or cross == 'identity' or g[0] in DERIVED:
            cross_mats[g] = np.eye(4)
        else:
            tag = f'{g[0]}_{g[1]}'.replace('.', 'p')
            cross_mats[g] = _flirt(tpl, frame_file, work / f'cross_{tag}.mat', cost)
    final = [None] * len(idx)
    for g, ks in groups.items():
        for k in ks:
            final[k] = cross_mats[g] @ mats[k]
    aslref = work / 'aslref.nii.gz'
    shutil.copyfile(frame_file, aslref)
    ctx_file = work / 'hmc_aslcontext.tsv'
    pd.DataFrame({'volume_type': types}).to_csv(ctx_file, sep='\t', index=False)
    return {
        'aslref': str(aslref),
        'xforms': str(_to_itk(final, aslref, work)),
        'aslcontext': str(ctx_file),
        'indices': idx,
        'reference_type': f'{frame} frame ({Path(frame_file).name})',
        'processing_target': target,
        'use_ge': None,
        'groups': {f'{g[0]}@{g[1]}': len(ks) for g, ks in groups.items()},
    }


def run_c2(inputs, work, **options):
    options = {'iterations': 3, **options}
    return run_c1(inputs, work, **options)


def _write_itk(path, mats):
    """Write RAS pull matrices as an ITK multi-transform file (MCFLIRT2ITK's layout)."""
    flip = np.diag([-1.0, -1.0, 1.0, 1.0])
    lines = ['#Insight Transform File V1.0']
    for k, m in enumerate(mats):
        lps = flip @ np.asarray(m, dtype=float) @ flip
        params = [*lps[:3, :3].ravel(), *lps[:3, 3]]
        lines += [
            f'#Transform {k}',
            'Transform: AffineTransform_double_3_3',
            f'Parameters: {" ".join(repr(float(v)) for v in params)}',
            'FixedParameters: 0 0 0',
        ]
    Path(path).write_text('\n'.join(lines) + '\n')
    return path


def _base(inputs, work, base, options):
    from harness import CANDIDATES, load_itk_series

    sub = work / f'base_{base}'
    sub.mkdir()
    res = CANDIDATES[base](inputs, sub, **{k: v for k, v in options.items() if k != 'base'})
    return res, load_itk_series(res['xforms'])


def run_c3(inputs, work, **options):
    """Production's transforms re-expressed relative to one acquired volume."""
    res, h = _base(inputs, work, 'c0', {})
    types = pd.read_table(res['aslcontext'])['volume_type'].tolist()
    ref_k = _opt(options, 'ref_volume', None, int)
    if ref_k is None:
        ref_k = types.index('control') if 'control' in types else 0
    # A point y of volume r's grid is the aslref point H_r^-1 y, so H'_v = H_v H_r^-1.
    new = [m @ np.linalg.inv(h[ref_k]) for m in h]
    img = nb.load(inputs['asl_file'])
    aslref = _save(
        np.asanyarray(img.dataobj[..., res['indices'][ref_k]], np.float32),
        img,
        work / 'aslref.nii.gz',
    )
    return {
        **res,
        'aslref': str(aslref),
        'xforms': str(_write_itk(work / 'hmc_xfms.txt', new)),
        'reference_type': f're-referenced to volume {res["indices"][ref_k]}',
    }


def run_c4(inputs, work, **options):
    """One transform per control-label pair: the polar mean of the pair's two transforms."""
    from score import polar_mean

    base = _opt(options, 'base', 'c0')
    res, h = _base(inputs, work, base, options)
    types = pd.read_table(res['aslcontext'])['volume_type'].tolist()
    controls = [i for i, t in enumerate(types) if t == 'control']
    labels = [i for i, t in enumerate(types) if t == 'label']
    new = list(h)
    for c, lab in zip(controls, labels, strict=True):
        new[c] = new[lab] = polar_mean([h[c], h[lab]])
    return {
        **res,
        'xforms': str(_write_itk(work / 'hmc_xfms.txt', new)),
        'reference_type': f'{res["reference_type"]}; pairs shared ({base})',
    }


CANDIDATES = {'c1': run_c1, 'c2': run_c2, 'c3': run_c3, 'c4': run_c4}
