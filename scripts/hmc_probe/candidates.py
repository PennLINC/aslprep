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
- ``c8``: zig-zag regression (Wang 2012). Control and label volumes are registered together to
  one reference (``reference=mean``, their mean, or ``reference=production``, the separate M0);
  the label/control alternation (one regressor per delay) is regressed out of the six rigid
  parameters about the field-of-view centre; in-series M0 volumes are registered separately,
  as their own group, and derived volumes follow P1. ``convention=spm`` reproduces ASLtbx and
  ExploreASL's ``cleanandsave_parameters`` instead: SPM parameters of the world-space
  transforms relative to the first volume, a positional +/-1 regressor projected out with no
  intercept, then re-referenced to the first volume (``fit`` and ``regressors`` are ignored).

Comparators from #513: ``fixvol`` registers every kept volume to one acquired raw volume;
``perfref`` re-registers each pair's corrected delta-M to a preliminary perfusion image built
after a first pass (``base``) and applies the residual to both volumes of the pair.

C1 and C2 take ``engine=mcflirt|volreg|ants`` for the series step (AFNI 3dvolreg, least
squares; ANTs rigid, mutual information); everything else is FSL.

All registrations use FSL (MCFLIRT for series, FLIRT for single images) with ``cost``
(``mutualinfo`` by default). Transforms are composed as FSL matrices on the raw grid and
converted to ITK once.
"""

import json
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


def _pull_ras_to_fsl(pulls, img):
    """RAS pull matrices (reference point -> volume point) -> FSL matrices (volume -> reference)."""
    f2w = _fsl_to_world(img)
    return [np.linalg.inv(f2w) @ np.linalg.inv(p) @ f2w for p in pulls]


def _volreg(series, ref, out_base):
    """AFNI 3dvolreg (least squares, heptic) of a 4D file to a reference; FSL matrices.

    ``-1Dmatrix_save`` writes, per volume, the 3x4 matrix from base to input coordinates in
    AFNI's DICOM order (x right-to-left, y anterior-to-posterior: LPS).
    """
    from nitransforms.io.afni import AFNILinearTransformArray

    aff = f'{out_base}.aff12.1D'
    cmd = ['3dvolreg', '-base', str(ref), '-prefix', f'{out_base}.nii.gz', '-1Dmatrix_save', aff]
    subprocess.run([*cmd, str(series)], check=True, capture_output=True)
    # AFNI deobliques implicitly: its matrices are in the cardinal (plumb) DICOM frames of the
    # two headers, so a plain LPS/RAS sign flip is wrong for oblique images (until revision 13
    # it was used, which was right only for the simulated, non-oblique fixtures).
    series_img, ref_img = nb.load(series), nb.load(ref)
    pulls = AFNILinearTransformArray.from_filename(aff).to_ras(
        moving=series_img, reference=ref_img
    )
    return _pull_ras_to_fsl(list(np.asarray(pulls).reshape(-1, 4, 4)), series_img)


def _ants(series, ref, out_base):
    """ANTs rigid (mutual information, three levels) of each volume to a reference; FSL matrices.

    The ITK transform maps fixed (reference) points to moving (volume) points, in LPS;
    nitransforms reads it as a RAS pull matrix.
    """
    import nitransforms as nt

    img = nb.load(series)
    data = np.asanyarray(img.dataobj, np.float32)
    pulls = []
    for k in range(data.shape[-1]):
        vol = _save(data[..., k], img, Path(f'{out_base}_v{k:03d}.nii.gz'))
        prefix = f'{out_base}_v{k:03d}_'
        cmd = ['antsRegistration', '--dimensionality', '3', '--float', '1']
        cmd += ['--output', prefix, '--interpolation', 'Linear']
        cmd += ['--transform', 'Rigid[0.1]', '--metric', f'MI[{ref},{vol},1,32,Regular,0.25]']
        cmd += ['--convergence', '[100x50x25,1e-6,10]', '--shrink-factors', '4x2x1']
        cmd += ['--smoothing-sigmas', '2x1x0vox']
        subprocess.run(cmd, check=True, capture_output=True)
        xfm = nt.linear.load(f'{prefix}0GenericAffine.mat', fmt='itk')
        pulls.append(np.asarray(xfm.matrix))
    return _pull_ras_to_fsl(pulls, img)


def _nearest_rigid(m, centre):
    """The rigid matrix nearest an affine ``m`` (polar factor), agreeing with it at ``centre``."""
    u, _, vt = np.linalg.svd(m[:3, :3])
    r = u @ vt
    if np.linalg.det(r) < 0:
        u[:, -1] *= -1
        r = u @ vt
    out = np.eye(4)
    out[:3, :3] = r
    out[:3, 3] = m[:3, :3] @ centre + m[:3, 3] - r @ centre
    return out


def _niimath(series, ref, out_base, rigid_ls=False):
    """niimath ``-allineate`` (a port of AFNI 3dAllineate) of each volume to a reference.

    The default (fast) engine fits a full affine only, so its result is projected onto the
    nearest rigid transform about the reference's centre; ``rigid_ls`` uses the classic engine
    with a rigid warp (``-warp shr``) and least squares instead. ``-savemat`` writes the
    world-space pull matrix (``fixed_to_moving``).
    """
    import os

    exe = os.environ.get('NIIMATH', 'niimath')
    img = nb.load(series)
    data = np.asanyarray(img.dataobj, np.float32)
    ref_img = nb.load(ref)
    centre = ref_img.affine @ np.r_[(np.array(ref_img.shape[:3]) - 1) / 2, 1]
    pulls = []
    for k in range(data.shape[-1]):
        vol = _save(data[..., k], img, Path(f'{out_base}_v{k:03d}.nii.gz'))
        mat = Path(f'{out_base}_v{k:03d}.json')
        classic = ['-cost', 'ls', '-warp', 'shr']
        opts = classic if rigid_ls else ['-cost', 'fast']
        cmd = [exe, str(vol), '-allineate', str(ref), '-savemat', str(mat)]
        out = f'{out_base}_v{k:03d}_out.nii.gz'
        proc = subprocess.run([*cmd, *opts, out], capture_output=True, text=True, check=False)
        used_classic = rigid_ls
        if proc.returncode and not rigid_ls and 'pyramid' in proc.stderr:
            # The fast engine cannot build its pyramid on coarse grids (e.g. 4 x 4 x 6 mm).
            print(
                f'niimath: fast engine failed on volume {k} ({proc.stderr.strip()}); '
                'using the classic rigid least-squares engine',
                flush=True,
            )
            proc = subprocess.run(
                [*cmd, *classic, out], capture_output=True, text=True, check=False
            )
            used_classic = True
        if proc.returncode:
            raise RuntimeError(f'niimath failed on volume {k}: {proc.stderr.strip()}')
        pull = np.asarray(json.loads(mat.read_text())['fixed_to_moving'], float)
        pulls.append(pull if used_classic else _nearest_rigid(pull, centre[:3]))
    return _pull_ras_to_fsl(pulls, img)


def _register_series(series, ref, out_base, cost, engine):
    """Register each volume of ``series`` to ``ref`` with ``engine``.

    Returns FSL matrices and, for MCFLIRT, its own resampled series (``None`` otherwise).
    ``cost`` applies to MCFLIRT only: 3dvolreg is least squares, ANTs uses mutual information.
    """
    if engine == 'mcflirt':
        return _mcflirt(series, ref, out_base, cost)
    if engine == 'volreg':
        return _volreg(series, ref, out_base), None
    if engine == 'ants':
        return _ants(series, ref, out_base), None
    if engine == 'niimath':
        return _niimath(series, ref, out_base), None
    if engine == 'niimath_ls':
        return _niimath(series, ref, out_base, rigid_ls=True), None
    raise ValueError(f'Unknown engine: {engine}')


def _apply_fsl(data, mats, img):
    """Resample each volume of ``data`` into the reference grid (linear), from FSL matrices."""
    from scipy.ndimage import map_coordinates

    fsl2vox = np.linalg.inv(img.affine) @ _fsl_to_world(img)
    ijk = np.indices(img.shape[:3]).reshape(3, -1)
    hom = np.vstack([ijk, np.ones(ijk.shape[1])])
    out = np.empty_like(data)
    for k, m in enumerate(mats):
        # Reference voxel -> reference FSL -> volume FSL (m^-1) -> volume voxel.
        src = fsl2vox @ np.linalg.inv(m) @ np.linalg.inv(fsl2vox) @ hom
        out[..., k] = map_coordinates(data[..., k], src[:3], order=1, mode='nearest').reshape(
            img.shape[:3]
        )
    return out


def _flirt(moving, ref, out_mat, cost):
    """Rigid registration of one image to another; the FSL matrix (moving -> reference).

    MCFLIRT on the single volume (local optimization from the identity), as the plan specifies
    for the cross-type step. FLIRT's global search was used until revision 13 and is not used:
    between near-identical images (label and control templates, correlation 0.9998 on real
    data) its cost is lower off the voxel grid, where interpolation smooths the noise, so it
    returned spurious 1.5-2 degree rotations on real data (``ds008020``, ``ds004856``,
    ``ds006148``) while MCFLIRT and 3dvolreg returned the identity.
    """
    # MCFLIRT writes its matrices into ``<base>.mat/``, so the base must not be ``out_mat``'s stem.
    base = Path(out_mat).with_name(f'{Path(out_mat).stem}_mcf')
    mats, _ = _mcflirt(moving, ref, base, cost)
    np.savetxt(out_mat, mats[0])
    return mats[0]


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
    engine = _opt(options, 'engine', 'mcflirt')

    target, idx, types = _kept(inputs)
    img = nb.load(inputs['asl_file'])
    data = np.asanyarray(img.dataobj, np.float32)
    data = data[..., None] if data.ndim == 3 else data
    plds = _plds(inputs, data.shape[-1], idx)
    groups = _groups(types, plds, by_pld, min_group)

    mats = [None] * len(idx)
    templates = {}
    corrected_groups = set()
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
            got, moved = _register_series(series, tpl, work / f'reg_{tag}_it{it}', cost, engine)
            if iterations > 1:  # the next template, and after the last pass the final one
                if moved is not None:  # MCFLIRT's own (spline) resampling, as before
                    moved_data = np.asanyarray(nb.load(moved).dataobj, np.float32)
                else:
                    moved_data = _apply_fsl(vols, got, img)
                tpl = _save(moved_data.mean(-1), img, work / f'tpl_{tag}_it{it + 1}.nii.gz')
        for k, m in zip(ks, got, strict=True):
            mats[k] = m
        templates[g] = tpl
        if iterations > 1:
            # The last template is rebuilt from the final transforms, so it is an average of
            # corrected volumes: it depicts their corrected frame (score.py's 'corrected').
            corrected_groups.add(g)

    # The common frame.
    controls = [g for g in templates if g[0] == 'control']
    m0_groups = [g for g in templates if g[0] == 'm0scan']
    frame_kind = 'raw'
    if frame == 'm0' and inputs.get('m0scan'):
        frame_file = _regrid(inputs['m0scan'], img, work / 'frame_m0.nii.gz')
        ref_volumes = 'separate_m0'
    else:
        if frame == 'm0' and m0_groups:
            g = m0_groups[0]
        elif controls:
            # With several control groups (one per delay), the longest-delay group's template.
            g = max(controls, key=lambda g: (g[1] is None, g[1] or 0))
        else:
            g = next(iter(templates))
        frame_file = templates[g]
        ref_volumes = [idx[k] for k in groups[g]]
        frame_kind = 'corrected' if g in corrected_groups else 'raw'

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
        'reference_volumes': ref_volumes,
        'reference_frame_kind': frame_kind,
        'processing_target': target,
        'use_ge': None,
        'groups': {f'{g[0]}@{g[1]}': len(ks) for g, ks in groups.items()},
    }


def run_c2(inputs, work, **options):
    options = {'iterations': 3, **options}
    return run_c1(inputs, work, **options)


def polar_mean(mats):
    """Mean of rigid 4x4 matrices, with the linear part projected onto a rotation.

    The same as ``score.polar_mean``, repeated here because the harness runs in a container
    without the repository.
    """
    m = np.mean(np.asarray(mats), axis=0)
    u, _, vt = np.linalg.svd(m[:3, :3])
    r = u @ vt
    if np.linalg.det(r) < 0:
        u[:, -1] *= -1
        r = u @ vt
    out = np.eye(4)
    out[:3, :3] = r
    out[:3, 3] = m[:3, 3]
    return out


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
        np.asanyarray(
            img.dataobj if img.ndim == 3 else img.dataobj[..., res['indices'][ref_k]], np.float32
        ),
        img,
        work / 'aslref.nii.gz',
    )
    return {
        **res,
        'aslref': str(aslref),
        'xforms': str(_write_itk(work / 'hmc_xfms.txt', new)),
        'reference_type': f're-referenced to volume {res["indices"][ref_k]}',
        'reference_volumes': [res['indices'][ref_k]],
    }


def run_c4(inputs, work, **options):
    """One transform per control-label pair: the polar mean of the pair's two transforms."""
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


def _about_centre(m, c):
    """FSL matrix -> (rotation vector in degrees, translation in mm about centre ``c``)."""
    from scipy.spatial.transform import Rotation

    r, t = m[:3, :3], m[:3, 3]
    return np.r_[np.degrees(Rotation.from_matrix(r).as_rotvec()), r @ c + t - c]


def _from_centre(p, c):
    from scipy.spatial.transform import Rotation

    r = Rotation.from_rotvec(np.radians(p[:3])).as_matrix()
    out = np.eye(4)
    out[:3, :3] = r
    out[:3, 3] = p[3:] - r @ c + c
    return out


def _fsl_to_world(img):
    """FSL scaled-voxel coordinates -> the image's world (RAS mm) coordinates."""
    zooms = np.array(img.header.get_zooms()[:3], float)
    d = np.diag(np.r_[zooms, 1.0])
    if np.linalg.det(img.affine[:3, :3]) > 0:  # FSL flips x: x_fsl = (N - 1 - i) * zx
        d[0, 0] = -zooms[0]
        d[0, 3] = (img.shape[0] - 1) * zooms[0]
    return img.affine @ np.linalg.inv(d)


def spm_matrix(p):
    """SPM's ``spm_matrix`` for the six rigid parameters (mm, radians): T @ R1 @ R2 @ R3."""
    c, s = np.cos(p[3:6]), np.sin(p[3:6])
    t = np.eye(4)
    t[:3, 3] = p[:3]
    r1, r2, r3 = np.eye(4), np.eye(4), np.eye(4)
    r1[1:3, 1:3] = [[c[0], s[0]], [-s[0], c[0]]]
    r2[np.ix_([0, 2], [0, 2])] = [[c[1], s[1]], [-s[1], c[1]]]
    r3[:2, :2] = [[c[2], s[2]], [-s[2], c[2]]]
    return t @ r1 @ r2 @ r3


def spm_imatrix(m):
    """SPM's ``spm_imatrix`` for a rigid matrix: the six parameters (mm, radians)."""
    r = m[:3, :3]
    p5 = np.arcsin(np.clip(r[0, 2], -1.0, 1.0))
    if abs(p5) - np.pi / 2 >= -1e-12:  # gimbal lock, as SPM handles it
        p4, p6 = 0.0, np.arctan2(-r[1, 0], -r[2, 0] / r[0, 2])
    else:
        c = np.cos(p5)
        p4, p6 = np.arctan2(r[1, 2] / c, r[2, 2] / c), np.arctan2(r[0, 1] / c, r[0, 0] / c)
    return np.r_[m[:3, 3], p4, p5, p6]


def spm_zigzag_clean(params, n_first_label=True):
    """ASLtbx/ExploreASL ``cleanandsave_parameters``, transcribed.

    ``params`` are SPM parameters of each volume relative to the first. The regressor is a
    hard-coded +/-1 alternation by position (-1 first, as when label comes first; the sign does
    not change the projection), projected out with no intercept, then every column is
    re-referenced to the first volume.
    """
    n = len(params)
    ref = -np.ones(n) if n_first_label else np.ones(n)
    ref[1::2] *= -1
    cleaned = params - np.outer(ref, ref @ params) / (ref @ ref)
    cleaned = cleaned - cleaned[0]
    return cleaned, params - cleaned, ref


def _huber(y, x, k=1.345, iterations=50):
    """Huber regression by iteratively reweighted least squares (MAD scale)."""
    b, *_ = np.linalg.lstsq(x, y, rcond=None)
    for _ in range(iterations):
        r = y - x @ b
        scale = np.median(np.abs(r - np.median(r))) / 0.6745 + 1e-9
        w = np.sqrt(np.minimum(1.0, k * scale / np.maximum(np.abs(r), 1e-12)))
        b, *_ = np.linalg.lstsq(x * w[:, None], y * w, rcond=None)
    return b


def zigzag_design(types, plds, regressors='pooled'):
    """The design: an intercept plus label/control alternation regressor(s).

    ``pooled`` (Wang 2012): one +1/-1 (control/label) regressor over the whole series.
    ``per_delay`` (an extension for multi-delay data): one regressor per delay, non-zero only
    for that delay's volumes. With one pair at a delay, its regressor fits that pair's whole
    difference, which makes the pair's two cleaned transforms identical: a strict pair
    constraint, not zig-zag regression. So ``per_delay`` is used only when every delay has at
    least two pairs; otherwise the pooled design is used and ``fell_back`` is True.
    """
    n = len(types)
    alt = np.array([1.0 if t == 'control' else -1.0 for t in types])
    fell_back = False
    if regressors == 'per_delay':
        delays = sorted(set(plds))
        pairs_per_delay = [
            sum(1 for t, p in zip(types, plds, strict=True) if p == d and t == 'control')
            for d in delays
        ]
        if min(pairs_per_delay) < 2:
            regressors, fell_back = 'pooled', True
    cols = [np.ones(n)]
    if regressors == 'pooled':
        cols.append(alt - alt.mean())
    else:
        for d in sorted(set(plds)):
            col = np.where(np.asarray(plds) == d, alt, 0.0)
            nz = col != 0
            col[nz] -= col[nz].mean()
            cols.append(col)
    return np.column_stack(cols), regressors, fell_back


def zigzag_clean(params, types, plds, robust=False, regressors='pooled'):
    """Remove the label/control alternation from motion parameters.

    ``params`` is (n, 6) for the control and label volumes in acquisition order. Only the fitted
    alternation is subtracted; the intercept, drift and spikes stay. Returns the cleaned
    parameters, the removed component, and an audit record (design, coefficients, fallback).
    """
    x, used, fell_back = zigzag_design(types, plds, regressors)
    if robust:
        beta = np.column_stack([_huber(params[:, j], x) for j in range(params.shape[1])])
    else:
        beta, *_ = np.linalg.lstsq(x, params, rcond=None)
    removed = x[:, 1:] @ beta[1:]
    audit = {'design': x, 'beta': beta, 'regressors': used, 'fell_back_to_pooled': fell_back}
    return params - removed, removed, audit


def run_c8(inputs, work, **options):
    """Zig-zag regression (C8): joint estimation, alternation regressed out, M0 separate."""
    cost = _opt(options, 'cost', 'mutualinfo')
    reference = _opt(options, 'reference', 'mean')
    target, idx, types = _kept(inputs)
    cl = [k for k, t in enumerate(types) if t in ('control', 'label')]
    if not cl:
        return run_c1(inputs, work, cost=cost)
    img = nb.load(inputs['asl_file'])
    data = np.asanyarray(img.dataobj, np.float32)
    data = data[..., None] if data.ndim == 3 else data
    plds = _plds(inputs, data.shape[-1], idx)
    cl_data = data[..., [idx[k] for k in cl]]
    series = _save(cl_data, img, work / 'cl_series.nii.gz')
    if reference == 'production':
        # Production's own reference: its selection, RobustAverage, M0 regridding and GE
        # smoothing, built by the c0 path (its HMC is discarded).
        from harness import run_c0

        sub = work / 'production_reference'
        sub.mkdir()
        c0 = run_c0(inputs, sub)
        ref_file, ref_volumes = c0['aslref'], c0['reference_volumes']
    else:
        ref_file = _save(cl_data.mean(-1), img, work / 'ref_cl_mean.nii.gz')
        ref_volumes = [idx[k] for k in cl]
    mats, _ = _mcflirt(series, ref_file, work / 'mcf_cl', cost)
    zooms = np.array(img.header.get_zooms()[:3])
    # FSL coordinates flip x for a positive-determinant affine (i -> N-1-i), which leaves the
    # field-of-view midpoint unchanged, so this centre is right in either case.
    centre = (np.array(img.shape[:3]) - 1) / 2 * zooms
    convention = _opt(options, 'convention', 'centre')
    cl_types = [types[k] for k in cl]
    if convention == 'spm':
        # ASLtbx/ExploreASL: SPM parameters of world-space transforms relative to the first
        # volume (Euler angles about the world origin), positional alternation, no intercept.
        f2w = _fsl_to_world(img)
        world = [f2w @ m @ np.linalg.inv(f2w) for m in mats]
        params = np.array([spm_imatrix(w @ np.linalg.inv(world[0])) for w in world])
        cleaned, removed, ref = spm_zigzag_clean(params)
        is_control = np.array([t == 'control' for t in cl_types])
        audit = {
            'design': ref[:, None],
            'beta': (ref @ params / (ref @ ref))[None],
            'regressors': 'spm positional',
            'fell_back_to_pooled': False,
            'alternation_matches_context': bool(
                np.all(is_control == (ref > 0)) or np.all(is_control == (ref < 0))
            ),
        }
    elif convention == 'centre':
        params = np.array([_about_centre(m, centre) for m in mats])
        cleaned, removed, audit = zigzag_clean(
            params,
            cl_types,
            [plds[k] for k in cl],
            robust=_opt(options, 'fit', 'ols') == 'huber',
            regressors=_opt(options, 'regressors', 'pooled'),
        )
    else:
        raise ValueError(f'Unknown convention: {convention}')
    if _opt(options, 'regress', 1, int) == 0:  # joint estimation only (a diagnostic variant)
        cleaned, removed = params, np.zeros_like(params)
        audit['regressors'] = 'none'
    # Audit trail: what was estimated, the parameter convention, and what was removed.
    art = work / 'artifacts'
    art.mkdir()
    np.save(art / 'initial_fsl_mats.npy', np.asarray(mats))
    np.save(art / 'params_initial.npy', params)
    np.save(art / 'params_cleaned.npy', cleaned)
    np.save(art / 'removed.npy', removed)
    np.save(art / 'design.npy', audit['design'])
    np.save(art / 'beta.npy', audit['beta'])
    (art / 'zigzag.json').write_text(
        json.dumps(
            {
                'convention': (
                    'SPM spm_imatrix of world-space transforms relative to the first volume '
                    '(translation mm, Euler radians about the world origin)'
                    if convention == 'spm'
                    else 'rotation vector (degrees) and translation (mm) about the '
                    'field-of-view centre, in FSL scaled-voxel coordinates'
                )
                + ', per volume of the control-label series (acquisition order)',
                'centre_fsl_mm': centre.tolist(),
                'volumes': [idx[k] for k in cl],
                'regressors': audit['regressors'],
                'fell_back_to_pooled': audit['fell_back_to_pooled'],
                'alternation_matches_context': audit.get('alternation_matches_context'),
                'fit': 'ols, no intercept' if convention == 'spm' else _opt(options, 'fit', 'ols'),
                'removed_max': float(np.abs(removed).max()) if removed.size else 0.0,
            },
            indent=2,
        )
    )
    final = [None] * len(idx)
    for k, p in zip(cl, cleaned, strict=True):
        if convention == 'spm':  # Vo(j).mat = spm_matrix(clQ(j,:)) * V(1).mat
            final[k] = np.linalg.inv(f2w) @ spm_matrix(p) @ world[0] @ f2w
        else:
            final[k] = _from_centre(p, centre)

    # M0 volumes: their own group, registered to its mean, then the mean to the reference.
    m0 = [k for k, t in enumerate(types) if t == 'm0scan']
    if m0:
        m0_data = data[..., [idx[k] for k in m0]]
        tpl = _save(m0_data.mean(-1), img, work / 'm0_mean.nii.gz')
        if len(m0) > 1:
            m0_mats, _ = _mcflirt(
                _save(m0_data, img, work / 'm0_series.nii.gz'), tpl, work / 'mcf_m0', cost
            )
        else:
            m0_mats = [np.eye(4)]
        cross = _flirt(tpl, ref_file, work / 'cross_m0.mat', cost)
        for k, m in zip(m0, m0_mats, strict=True):
            final[k] = cross @ m
    # Derived volumes (policy P1): identity.
    final = [np.eye(4) if f is None else f for f in final]
    aslref = work / 'aslref.nii.gz'
    shutil.copyfile(ref_file, aslref)
    ctx_file = work / 'hmc_aslcontext.tsv'
    pd.DataFrame({'volume_type': types}).to_csv(ctx_file, sep='\t', index=False)
    return {
        'aslref': str(aslref),
        'xforms': str(_to_itk(final, aslref, work)),
        'aslcontext': str(ctx_file),
        'indices': idx,
        'reference_type': f'zig-zag ({reference} reference, {audit["regressors"]})',
        'reference_volumes': ref_volumes,
        'processing_target': target,
        'use_ge': None,
        'artifacts': str(art),
    }


def _m0_into(final, data, idx, types, img, ref_file, work, cost):
    """Fill the M0 volumes' FSL matrices as in C8 (their own group, then their mean to the
    reference), and the identity for any volume still unset (derived volumes)."""
    final = list(final)
    m0 = [k for k, t in enumerate(types) if t == 'm0scan']
    if m0:
        m0_data = data[..., [idx[k] for k in m0]]
        tpl = _save(m0_data.mean(-1), img, work / 'm0_mean.nii.gz')
        m0_mats = (
            _mcflirt(_save(m0_data, img, work / 'm0_series.nii.gz'), tpl, work / 'mcf_m0', cost)[0]
            if len(m0) > 1
            else [np.eye(4)]
        )
        cross = _flirt(tpl, ref_file, work / 'cross_m0.mat', cost)
        for k, m in zip(m0, m0_mats, strict=True):
            final[k] = cross @ m
    return [np.eye(4) if f is None else f for f in final]


def c5_eligibility(metadata, plds_cl, presat_delay):
    """Reasons C5 cannot apply (#673): an empty list means eligible."""
    reasons = []
    if metadata.get('MRAcquisitionType') != '2D':
        reasons.append('not a 2D readout')
    if not metadata.get('SliceTiming'):
        reasons.append('no SliceTiming')
    if metadata.get('BackgroundSuppression') and metadata.get(
        'BackgroundSuppressionNumberPulses', 0
    ):
        reasons.append('background-suppression inversion pulses alter the recovery')
    if len(set(plds_cl)) < 2:
        reasons.append('fewer than two delays sample the recovery')
    if presat_delay is None:
        reasons.append('pre-saturation timing unknown (give presat_delay)')
    if metadata.get('ArterialSpinLabelingType') != 'PCASL':
        reasons.append('only (P)CASL timing is modelled')
    return reasons


def fit_recovery(data, times, t1_grid=None):
    """Per-voxel saturation recovery ``S = M0 (1 - exp(-t / T1eff))`` by grid search over T1eff.

    ``data`` is (x, y, z, n); ``times`` is (z, n), the recovery time of each volume at each slice.
    For each candidate T1eff the best M0 is closed-form (least squares); the T1eff with the
    smallest residual is kept. Returns M0 and T1eff maps.
    """
    t1_grid = np.linspace(0.3, 4.0, 75) if t1_grid is None else t1_grid
    m0 = np.zeros(data.shape[:3], np.float32)
    t1 = np.full(data.shape[:3], np.nan, np.float32)
    for z in range(data.shape[2]):
        s = data[:, :, z, :].reshape(-1, data.shape[-1]).astype(np.float64)
        best_res = np.full(s.shape[0], np.inf)
        best_m0 = np.zeros(s.shape[0])
        best_t1 = np.full(s.shape[0], np.nan)
        for t1_c in t1_grid:
            b = 1.0 - np.exp(-times[z] / t1_c)
            m0_c = s @ b / (b @ b)
            res = ((s - np.outer(m0_c, b)) ** 2).sum(axis=1)
            better = res < best_res
            best_res[better], best_m0[better], best_t1[better] = res[better], m0_c[better], t1_c
        m0[:, :, z] = best_m0.reshape(data.shape[:2])
        t1[:, :, z] = best_t1.reshape(data.shape[:2])
    return m0, t1


def run_c5(inputs, work, **options):
    """Saturation-recovery normalization for motion estimation (C5; #673 ``motion-only``).

    Eligible runs only (see ``c5_eligibility``); an ineligible run is processed as production
    (C0) and says so. Pass 1 fits ``S = M0 (1 - exp(-t/T1eff))`` to the unaligned controls,
    with ``t = presat_delay + LabelingDuration + PLD + SliceTiming[slice]`` (pre-saturation
    ``presat_delay`` seconds before labeling starts; default 0, at labeling start), divides
    controls and labels alike by ``1 - exp(-t/T1eff)``, and registers the normalized volumes
    jointly to their mean (MCFLIRT). Pass 2 refits on the aligned controls, maps the fitted
    T1eff back into each volume's frame, normalizes again and re-registers. The transforms are
    applied to the original intensities by the harness. M0 volumes and derived volumes are
    handled as in C8.

    ``t1_model``: ``slice`` (default) fits one T1eff per slice to the brain-mean control signal,
    normalizes once and registers twice (to the normalized mean, then to the aligned mean).
    ``voxel`` is the voxelwise two-pass fit above; on the drift fixtures it is biased toward
    its own pass-1 poses, because the T1eff map's edges, placed at the estimated pose, are
    imprinted on the normalized images.
    """
    cost = _opt(options, 'cost', 'mutualinfo')
    presat_delay = _opt(options, 'presat_delay', 0.0, float)
    t1_model = _opt(options, 't1_model', 'slice')
    if t1_model not in ('slice', 'voxel'):
        raise ValueError(f't1_model must be slice or voxel, not {t1_model}')
    target, idx, types = _kept(inputs)
    meta = inputs['metadata']
    img = nb.load(inputs['asl_file'])
    data = np.asanyarray(img.dataobj, np.float32)
    data = data[..., None] if data.ndim == 3 else data
    plds = _plds(inputs, data.shape[-1], idx)
    cl = [k for k, t in enumerate(types) if t in ('control', 'label')]
    reasons = c5_eligibility(meta, [plds[k] for k in cl], presat_delay)
    art = work / 'artifacts'
    art.mkdir()
    if reasons or not cl:
        from harness import run_c0

        sub = work / 'refused_c0'
        sub.mkdir()
        res = run_c0(inputs, sub)
        (art / 'c5.json').write_text(json.dumps({'eligible': False, 'reasons': reasons}, indent=2))
        return {
            **res,
            'reference_type': f'{res["reference_type"]}; C5 refused',
            'artifacts': str(art),
        }

    slice_times = np.asarray(meta['SliceTiming'], float)
    if slice_times.size != data.shape[2]:
        raise ValueError(f'{slice_times.size} slice times for {data.shape[2]} slices')
    tau = float(meta['LabelingDuration'])
    cl_data = data[..., [idx[k] for k in cl]]
    cl_types = [types[k] for k in cl]
    times = np.array(
        [[presat_delay + tau + plds[k] + slice_times[z] for k in cl] for z in range(data.shape[2])]
    )
    is_control = np.array([t == 'control' for t in cl_types])

    def normalize(t1):
        factor = 1.0 - np.exp(
            -times[None, None, :, :] / np.where(np.isfinite(t1), t1, 1.0)[..., None]
        )
        return cl_data / np.maximum(factor, 0.05)

    # Pass 1: fit on the unaligned controls, normalize, register jointly to the normalized mean.
    if t1_model == 'slice':
        # One T1eff per slice, fitted to the mean control signal of the brain voxels in it: the
        # factor depends only on (slice, volume), in the acquisition frame, so it cannot carry
        # anatomy at an assumed pose (the voxelwise map does, and biases registration toward it).
        ctrl = cl_data[..., is_control]
        brain = ctrl.mean(-1) > 0.2 * np.percentile(ctrl.mean(-1), 99)
        per_slice = np.array(
            [
                ctrl[:, :, z][brain[:, :, z]].mean(0)
                if brain[:, :, z].sum() >= 20
                else np.full(ctrl.shape[-1], np.nan)
                for z in range(ctrl.shape[2])
            ]
        )
        ok = np.isfinite(per_slice).all(1)
        _, t1_slice = fit_recovery(per_slice[ok][None, None], times[ok][:, is_control])
        t1_z = np.full(ctrl.shape[2], float(np.nanmedian(t1_slice)))
        t1_z[ok] = t1_slice[0, 0]
        t1 = np.broadcast_to(t1_z[None, None, :], cl_data.shape[:3]).astype(np.float32)
        norm2 = normalize(t1)
        mats = None
        t1_ref = t1
    else:
        _, t1 = fit_recovery(cl_data[..., is_control], times[:, is_control])
        norm = normalize(t1)
        series = _save(norm, img, work / 'norm1.nii.gz')
        ref1 = _save(norm.mean(-1), img, work / 'norm1_mean.nii.gz')
        mats, _ = _mcflirt(series, ref1, work / 'mcf1', cost)
        # Pass 2: refit on the aligned controls, map T1eff back to each volume, normalize.
        aligned = _apply_fsl(cl_data, mats, img)
        _, t1_ref = fit_recovery(aligned[..., is_control], times[:, is_control])
        t1_vols = _apply_fsl(
            np.repeat(np.nan_to_num(t1_ref, nan=1.0)[..., None], len(cl), -1),
            [np.linalg.inv(m) for m in mats],
            img,
        )
        factor = 1.0 - np.exp(-times[None, None, :, :] / np.maximum(t1_vols, 0.05))
        norm2 = cl_data / np.maximum(factor, 0.05)
    if mats is None:  # slice model: one normalization, then the same two-step registration
        mats, _ = _mcflirt(
            _save(norm2, img, work / 'norm1.nii.gz'),
            _save(norm2.mean(-1), img, work / 'norm1_mean.nii.gz'),
            work / 'mcf1',
            cost,
        )
    moved2 = _apply_fsl(norm2, mats, img)
    ref_file = _save(moved2.mean(-1), img, work / 'norm2_mean.nii.gz')
    mats2, _ = _mcflirt(_save(norm2, img, work / 'norm2.nii.gz'), ref_file, work / 'mcf2', cost)
    np.save(art / 't1eff_pass1.npy', t1)
    np.save(art / 't1eff_pass2.npy', t1_ref)
    (art / 'c5.json').write_text(
        json.dumps(
            {
                'eligible': True,
                't1_model': t1_model,
                'presat_delay_s': presat_delay,
                'labeling_duration_s': tau,
                'recovery_time_range_s': [float(times.min()), float(times.max())],
                't1eff_median_s': float(np.nanmedian(t1_ref[t1_ref > 0])),
            },
            indent=2,
        )
    )

    final = [None] * len(idx)
    for k, m in zip(cl, mats2, strict=True):
        final[k] = m
    final = _m0_into(final, data, idx, types, img, ref_file, work, cost)
    aslref = work / 'aslref.nii.gz'
    shutil.copyfile(ref_file, aslref)
    ctx_file = work / 'hmc_aslcontext.tsv'
    pd.DataFrame({'volume_type': types}).to_csv(ctx_file, sep='\t', index=False)
    return {
        'aslref': str(aslref),
        'xforms': str(_to_itk(final, aslref, work)),
        'aslcontext': str(ctx_file),
        'indices': idx,
        'reference_type': 'saturation-recovery normalized mean (two passes)',
        'reference_volumes': [idx[k] for k in cl],
        'reference_frame_kind': 'corrected',
        'processing_target': target,
        'use_ge': None,
        'artifacts': str(art),
    }


def c6_fit_inputs(
    inputs,
    init='zero',
    init_mm=1.0,
    seed=0,
    scale_iters=1.0,
    hires=2,
    support_mask=None,
    p_bound=None,
):
    """Fit the C6 model to a run's control and label volumes (needs PyTorch and, ideally, a GPU).

    The perfusion support is the static-image mask (mean control above 20% of its 99th
    percentile), the form that passed the identifiability test. ``init``: ``zero`` or ``random``
    (uniform within ``init_mm`` mm and degrees, for the initialization test). ``hires``: the
    acquisition model (maps on an in-plane grid ``hires`` times finer, k-space cropped, with a
    learned weight per phase-encoding line on the axis ``PhaseEncodingDirection`` names); 1 is
    the plain interpolation model. ``support_mask`` replaces the static-image support (a boolean
    array on the run's grid; the probe passes the true tissue mask as an upper bound). Returns
    None for a run without both controls and labels.
    """
    from c6_model import fit

    _, idx, types = _kept(inputs)
    img = nb.load(inputs['asl_file'])
    data = np.asanyarray(img.dataobj, np.float32)
    data = data[..., None] if data.ndim == 3 else data
    plds = _plds(inputs, data.shape[-1], idx)
    cl = [k for k, t in enumerate(types) if t in ('control', 'label')]
    if not {'control', 'label'} <= {types[k] for k in cl}:
        return None
    cl_data = data[..., [idx[k] for k in cl]].astype(np.float64)
    ctrl_mean = cl_data[..., [types[k] == 'control' for k in cl]].mean(-1)
    support = ctrl_mean > 0.2 * np.percentile(ctrl_mean, 99)
    if support_mask is not None:
        support = np.asarray(support_mask, bool)
    if data.shape[0] % 2 or data.shape[1] % 2:
        hires = 1  # the centred k-space crop needs even in-plane sizes
    rng = np.random.default_rng(seed)
    start = None if init == 'zero' else rng.uniform(-init_mm, init_mm, (len(cl), 6))
    stages = tuple(
        (fwhm, max(1, round(it * scale_iters)))
        for fwhm, it in ((None, 200), (4.0, 300), (2.0, 300), (0.0, 400))
    )
    res = fit(
        cl_data,
        [types[k] for k in cl],
        [plds[k] for k in cl],
        img.affine,
        support,
        init=start,
        stages=stages,
        seed=seed,
        hires=hires,
        pe_axis={'i': 0, 'j': 1}.get(
            str(inputs['metadata'].get('PhaseEncodingDirection', 'j'))[0]
        ),
        p_bound=p_bound,
    )
    return {
        **res,
        'cl': cl,
        'init': init,
        'init_mm': init_mm,
        'seed': seed,
        'start': start,
        'hires': hires,
        'p_bound': p_bound,
    }


def c6_fit_file(fit_dir, asl_file):
    return Path(fit_dir) / f'{Path(asl_file).name.removesuffix(".nii.gz")}.npz'


def run_c6(inputs, work, **options):
    """Differentiable contrast model fitted jointly with the motion (C6; ``c6_model``).

    The fit needs PyTorch, which the ASLPrep image does not have: ``c6_fit.py`` fits each run on
    the host and ``fit_dir`` (mounted read-only) holds the results; without ``fit_dir`` the fit
    runs here. Control and label volumes only; a run without both is processed as production
    and says so. The static image is the reference handed on; M0 volumes are handled as in C8.
    """
    cost = _opt(options, 'cost', 'mutualinfo')
    fit_dir = _opt(options, 'fit_dir', None)
    target, idx, types = _kept(inputs)
    img = nb.load(inputs['asl_file'])
    data = np.asanyarray(img.dataobj, np.float32)
    data = data[..., None] if data.ndim == 3 else data
    art = work / 'artifacts'
    art.mkdir()
    if fit_dir is not None:
        f = c6_fit_file(fit_dir, inputs['asl_file'])
        z = np.load(f, allow_pickle=False)
        res = (
            {k: z[k] for k in ('pulls', 'theta', 'S', 'P', 'cl')} | json.loads(str(z['meta']))
            if z['applied']
            else None
        )
    else:
        res = c6_fit_inputs(
            inputs,
            _opt(options, 'init', 'zero'),
            _opt(options, 'init_mm', 1.0, float),
            _opt(options, 'seed', 0, int),
            _opt(options, 'iters', 1.0, float),
            _opt(options, 'hires', 2, int),
        )
    if res is None:
        from harness import run_c0

        sub = work / 'refused_c0'
        sub.mkdir()
        res = run_c0(inputs, sub)
        (art / 'c6.json').write_text(json.dumps({'applied': False, 'reason': 'no pairs'}))
        return {
            **res,
            'reference_type': f'{res["reference_type"]}; C6 not applicable',
            'artifacts': str(art),
        }

    ref_file = _save(res['S'], img, work / 'c6_static.nii.gz')
    _save(res['P'], img, art / 'c6_perfusion.nii.gz')
    np.save(art / 'theta.npy', res['theta'])
    meta = {k: res.get(k) for k in ('init', 'init_mm', 'seed', 'hires', 'p_bound', 'w', 'history')}
    (art / 'c6.json').write_text(json.dumps({'applied': True, **meta}, indent=2))
    final = [None] * len(idx)
    for k, m in zip(res['cl'], _pull_ras_to_fsl(res['pulls'], img), strict=True):
        final[k] = m
    final = _m0_into(final, data, idx, types, img, ref_file, work, cost)
    aslref = work / 'aslref.nii.gz'
    shutil.copyfile(ref_file, aslref)
    ctx_file = work / 'hmc_aslcontext.tsv'
    pd.DataFrame({'volume_type': types}).to_csv(ctx_file, sep='\t', index=False)
    return {
        'aslref': str(aslref),
        'xforms': str(_to_itk(final, aslref, work)),
        'aslcontext': str(ctx_file),
        'indices': idx,
        'reference_type': 'C6 static image (joint contrast and motion model)',
        'reference_volumes': [idx[int(k)] for k in res['cl']],
        'reference_frame_kind': 'corrected',
        'processing_target': target,
        'use_ge': None,
        'artifacts': str(art),
    }


def run_fixvol(inputs, work, **options):
    """Comparator (#513): every kept volume registered by MCFLIRT to one acquired raw volume.

    ``ref_volume`` is a position in the kept series (default: the first control, else the
    first volume), as MCFLIRT's own ``-refvol`` would use a raw volume.
    """
    cost = _opt(options, 'cost', 'mutualinfo')
    target, idx, types = _kept(inputs)
    img = nb.load(inputs['asl_file'])
    data = np.asanyarray(img.dataobj, np.float32)
    data = data[..., None] if data.ndim == 3 else data
    ref_k = _opt(options, 'ref_volume', None, int)
    if ref_k is None:
        ref_k = types.index('control') if 'control' in types else 0
    series = _save(data[..., idx], img, work / 'series.nii.gz')
    ref_file = _save(data[..., idx[ref_k]], img, work / 'aslref.nii.gz')
    mats, _ = _mcflirt(series, ref_file, work / 'mcf', cost)
    ctx_file = work / 'hmc_aslcontext.tsv'
    pd.DataFrame({'volume_type': types}).to_csv(ctx_file, sep='\t', index=False)
    return {
        'aslref': str(ref_file),
        'xforms': str(_to_itk(mats, ref_file, work)),
        'aslcontext': str(ctx_file),
        'indices': idx,
        'reference_type': f'fixed raw volume {idx[ref_k]}',
        'reference_volumes': [idx[ref_k]],
        'processing_target': target,
        'use_ge': None,
    }


def run_perfref(inputs, work, **options):
    """Comparator (#513): a preliminary perfusion image as a second-pass HMC target.

    Pass 1 is ``base`` (``c0`` by default). Its corrected pairs give a robust perfusion image
    (``perfusion.combine``, level ``robust``, delays at or above the median); each pair's
    corrected delta-M is then registered rigidly to that image (same contrast, normalized
    correlation, MCFLIRT) and the residual is applied to both volumes of the pair.
    Needs ``min_pairs`` pairs; otherwise the run is pass 1 unchanged and says so.
    """
    from harness import resample_series
    from perfusion import combine, select_inputs

    base = _opt(options, 'base', 'c0')
    min_pairs = _opt(options, 'min_pairs', 4, int)
    res, h = _base(inputs, work, base, {})
    types = pd.read_table(res['aslcontext'])['volume_type'].tolist()
    idx = list(res['indices'])
    controls = [i for i, t in enumerate(types) if t == 'control']
    labels = [i for i, t in enumerate(types) if t == 'label']
    if len(controls) < min_pairs or len(controls) != len(labels):
        return {**res, 'reference_type': f'{res["reference_type"]}; perfref skipped (pairs)'}
    corrected = resample_series(inputs['asl_file'], idx, h, res['aslref'])
    cdata = np.asanyarray(corrected.dataobj, np.float32)
    raw_plds = _plds(inputs, nb.load(inputs['asl_file']).shape[-1], idx)
    stack, _, _ = select_inputs(cdata, types, raw_plds)
    ref_img = nb.load(res['aslref'])
    perf = _save(combine(stack, 'robust'), ref_img, work / 'perfusion_robust.nii.gz')
    f2w = _fsl_to_world(ref_img)
    new = list(h)
    for p, (c, lab) in enumerate(zip(controls, labels, strict=True)):
        dm = _save(cdata[..., c] - cdata[..., lab], ref_img, work / f'pair_{p:03d}.nii.gz')
        mat = _flirt(dm, perf, work / f'pair_{p:03d}.mat', 'normcorr')
        push = f2w @ mat @ np.linalg.inv(f2w)  # pair frame -> perfusion frame
        for v in (c, lab):
            new[v] = h[v] @ np.linalg.inv(push)
    return {
        **res,
        'xforms': str(_write_itk(work / 'hmc_xfms_perfref.txt', new)),
        'reference_type': f'{res["reference_type"]}; pairs re-registered to a perfusion image',
    }


CANDIDATES = {
    'c5': run_c5,
    'c6': run_c6,
    'c1': run_c1,
    'c2': run_c2,
    'c3': run_c3,
    'c4': run_c4,
    'c8': run_c8,
    'fixvol': run_fixvol,
    'perfref': run_perfref,
}
