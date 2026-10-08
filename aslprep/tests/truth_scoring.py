"""Score an ASLPrep run on an aslscan fixture against the ground truth.

See Sections 4.6 and 4.7 of ``docs/specs/2026-10-08-aslscan-ground-truth-tests-design.md``.
:func:`score_run` returns a nested dict of numbers (and lists of problems); it never raises
because an output is missing, and records the problem instead. Assertions live in the tests
(:mod:`aslprep.tests.truth_bounds`).

Frames (RAS point mappings, see :mod:`aslprep.tests.truth_geometry`):

- ``P_v``: the simulator's pose of volume ``v`` (static phantom -> moved);
- ``H_v``: ASLPrep's HMC transform, a pull mapping (aslref point -> volume-``v`` point);
- ``C``: ASLPrep's coregistration, a pull mapping (T1w point -> aslref point);
- ``R``: the anatomical offset (phantom point -> T1w point).

The estimated phantom-to-T1w mapping through volume ``v`` is ``E_v = C^-1 H_v^-1 P_v``, and
its error is ``D_v = E_v R^-1``. The aslref's pose is ``A_v = H_v^-1 P_v``.
"""

import json
from pathlib import Path

import nibabel as nb
import numpy as np

from aslprep.tests import truth_geometry as tg
from aslprep.tests import truth_models as tm

#: Denominator floor (ml/100 g/min) for ratios (spec Section 4.7).
FLOOR = 5.0
#: Partial-volume threshold of the tissue-dominant masks.
DOMINANT = 0.7
#: Displacement (mm) of the aslref pose below which native outputs compare voxel to voxel.
NATIVE_TOLERANCE_MM = 0.5


# ---------------------------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------------------------
class Fixture:
    """A generated fixture: raw inputs, ground truth, and generation records."""

    def __init__(self, path, acq):
        self.path = Path(path)
        self.acq = acq
        perf = self.path / 'sub-01' / 'perf'
        stem = f'sub-01_acq-{acq}'
        self.asl_img = nb.load(perf / f'{stem}_asl.nii.gz')
        self.metadata = json.loads((perf / f'{stem}_asl.json').read_text())
        self.context = (perf / f'{stem}_aslcontext.tsv').read_text().split()[1:]
        m0 = perf / f'{stem}_m0scan.nii.gz'
        self.m0scan = nb.load(m0).get_fdata() if m0.exists() else None
        self.m0_metadata = (
            json.loads((perf / f'{stem}_m0scan.json').read_text()) if m0.exists() else None
        )
        gt = perf / 'ground-truth'
        self.truth = json.loads((gt / 'truth.json').read_text())
        self.simulation = json.loads((gt / 'simulation.json').read_text())
        self.perfusion = nb.load(gt / 'sub-01_desc-perfusion_gt.nii.gz').get_fdata()
        self.att = nb.load(gt / 'sub-01_desc-att_gt.nii.gz').get_fdata()
        self.pv = {
            t: nb.load(gt / f'sub-01_desc-pv{t}_gt.nii.gz').get_fdata()
            for t in ('GM', 'WM', 'CSF')
        }
        self.motion = _read_tsv(gt / 'sub-01_desc-motion_gt.tsv')
        self.affine = self.asl_img.affine
        self.shape = self.asl_img.shape[:3]
        self.R = np.asarray(self.truth['R'])
        self.m0_scale = float(self.truth['m0_divisor'])

    @property
    def brain(self):
        pv = self.pv
        return pv['GM'] + pv['WM'] + pv['CSF'] >= 0.5

    def voxel_centers(self, mask):
        ijk = np.argwhere(mask)
        return ijk, tg.apply_points(self.affine, ijk)

    def poses(self):
        """``P_v`` for every volume (identity without motion)."""
        n = self.asl_img.shape[3]
        if self.motion is None:
            return np.repeat(np.eye(4)[None], n, axis=0)
        center = self.truth['motion_center']
        return np.stack(
            [
                tg.pose_matrix(
                    [row['trans_x'], row['trans_y'], row['trans_z']],
                    [row['rot_x'], row['rot_y'], row['rot_z']],
                    center,
                )
                for row in self.motion
            ]
        )


def _read_tsv(path):
    path = Path(path)
    if not path.exists():
        return None
    lines = path.read_text().strip().splitlines()
    header = lines[0].split('\t')
    rows = []
    for line in lines[1:]:
        values = line.split('\t')
        rows.append({k: _number(v) for k, v in zip(header, values, strict=True)})
    return rows


def _number(value):
    try:
        return float(value)
    except ValueError:
        return value


def find_output(aslprep_dir, acq, suffix, space=None, desc=None, ext='.nii.gz'):
    """An ASLPrep output by entities, or None. ``space=None`` means no space entity."""
    perf = Path(aslprep_dir) / 'sub-01' / 'perf'
    name = f'sub-01_acq-{acq}'
    if space:
        name += f'_space-{space}'
    pattern = name + (f'*_desc-{desc}_{suffix}{ext}' if desc else f'*_{suffix}{ext}')
    matches = [
        p
        for p in sorted(perf.glob(pattern))
        if (desc or '_desc-' not in p.name) and (space or '_space-' not in p.name)
    ]
    return matches[0] if len(matches) == 1 else None


def find_xfm(aslprep_dir, acq, name):
    perf = Path(aslprep_dir) / 'sub-01' / 'perf'
    matches = sorted(perf.glob(f'sub-01_acq-{acq}_{name}*_xfm.txt'))
    return matches[0] if len(matches) == 1 else None


def on_grid(img, affine, shape, order=1):
    """Data of ``img`` on the given grid (identity-copied when the grids already match)."""
    if img.ndim == 4 and img.shape[3] == 1:
        img = img.slicer[..., 0]
    if img.shape[:3] == tuple(shape) and np.allclose(img.affine, affine, atol=1e-4):
        return np.asanyarray(img.dataobj, dtype=np.float64)
    from nibabel.processing import resample_from_to

    target = nb.Nifti1Image(np.zeros(shape, np.int8), affine)
    return resample_from_to(img, target, order=order, cval=np.nan).get_fdata()


def sample_at(img, world_points, order=1):
    """Sample ``img`` at RAS world points (NaN outside the image)."""
    from scipy.ndimage import map_coordinates

    ijk = tg.apply_points(np.linalg.inv(img.affine), world_points)
    data = np.asanyarray(img.dataobj, dtype=np.float64)
    if data.ndim == 4 and data.shape[3] == 1:  # standard-space CBF carries a singleton time axis
        data = data[..., 0]
    return map_coordinates(data, ijk.T, order=order, mode='constant', cval=np.nan)


# ---------------------------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------------------------
def ratio_stats(values, reference, mask):
    """Median, tails and correlation of ``values / reference`` over ``mask``."""
    n = int(mask.sum())
    if n == 0:
        return {'n': 0}
    v, r = values[mask], reference[mask]
    finite = np.isfinite(v)
    ratio = v[finite] / r[finite]
    out = {
        'n': n,
        'finite': float(finite.mean()),
        'median': float(np.median(ratio)) if ratio.size else float('nan'),
        'median_dev': float(abs(np.median(ratio) - 1)) if ratio.size else float('nan'),
        'p05': float(np.percentile(ratio, 5)) if ratio.size else float('nan'),
        'p95': float(np.percentile(ratio, 95)) if ratio.size else float('nan'),
        'median_abs_dev': float(np.median(np.abs(ratio - 1))) if ratio.size else float('nan'),
        'p95_abs_dev': float(np.percentile(np.abs(ratio - 1), 95)) if ratio.size else float('nan'),
    }
    if ratio.size > 2 and np.std(v[finite]) > 0 and np.std(r[finite]) > 0:
        out['r'] = float(np.corrcoef(v[finite], r[finite])[0, 1])
    return out


def tissue_masks(fixture, extra=None, reference=None):
    """``valid`` and the tissue-dominant masks (spec Section 4.7)."""
    pv = fixture.pv
    reference = fixture.perfusion if reference is None else reference
    valid = fixture.brain & (pv['CSF'] < 0.5) & np.isfinite(reference) & (reference > FLOOR)
    if extra is not None:
        valid &= extra
    return {
        'valid': valid,
        'GM': valid & (pv['GM'] >= DOMINANT),
        'WM': valid & (pv['WM'] >= DOMINANT),
    }


# ---------------------------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------------------------
def expected_native(fixture, fwhm, mutation=None):
    """Tier A expectation for a single-delay fixture on the native grid."""
    return tm.expected_cbf(
        fixture.asl_img.get_fdata(),
        fixture.context,
        fixture.metadata,
        fixture.affine,
        m0scan=fixture.m0scan,
        m0_metadata=fixture.m0_metadata,
        m0_scale=fixture.m0_scale,
        fwhm=fwhm,
        mutation=mutation,
    )


def is_multi_delay(metadata):
    return len(set(np.atleast_1d(metadata['PostLabelingDelay']).tolist())) > 1


def score_native(fixture, aslprep_dir, fwhm, subset=300):
    """Tier A, Tier B, physical agreement and coverage for native-space ("asl") outputs."""
    out = {'problems': []}
    cbf_file = find_output(aslprep_dir, fixture.acq, 'cbf')
    mask_file = find_output(aslprep_dir, fixture.acq, 'mask', desc='brain')
    if cbf_file is None or mask_file is None:
        out['problems'].append('native cbf or brain mask not found')
        return out
    cbf_img = nb.load(cbf_file)
    out['grid_matches_input'] = bool(
        cbf_img.shape[:3] == fixture.shape
        and np.allclose(cbf_img.affine, fixture.affine, atol=1e-4)
    )
    cbf = on_grid(cbf_img, fixture.affine, fixture.shape)
    out_mask = on_grid(nb.load(mask_file), fixture.affine, fixture.shape, order=0) > 0.5
    out['coverage'] = {
        'mask_coverage': float((out_mask & fixture.brain).sum() / fixture.brain.sum()),
    }

    masks = tissue_masks(fixture, extra=out_mask)
    out['coverage']['n_GM'] = int(masks['GM'].sum())
    out['coverage']['n_WM'] = int(masks['WM'].sum())
    out['coverage']['finite'] = float(np.isfinite(cbf[masks['valid']]).mean())

    # Tier B against the physical truth
    out['tier_b'] = {t: ratio_stats(cbf, fixture.perfusion, masks[t]) for t in ('GM', 'WM')}
    gm, wm = masks['GM'], masks['WM']
    contrast_out = np.nanmedian(cbf[gm]) / np.nanmedian(cbf[wm])
    contrast_truth = np.median(fixture.perfusion[gm]) / np.median(fixture.perfusion[wm])
    out['tier_b']['gm_wm_contrast'] = float(contrast_out / contrast_truth)
    # how pure the tissue-dominant voxels are (report-only, to interpret Tier B)
    out['purity'] = {
        t: float(np.median(fixture.pv[t][masks[t]])) for t in ('GM', 'WM') if masks[t].any()
    }

    # ASLPrep's own preprocessed series (motion-corrected, in the aslref frame), if written
    preproc_file = find_output(aslprep_dir, fixture.acq, 'asl', desc='preproc')
    preproc = (
        on_grid(nb.load(preproc_file), fixture.affine, fixture.shape)
        if preproc_file is not None
        else None
    )
    if preproc is None:
        out['problems'].append('native desc-preproc_asl not found')
    interior = _interior(out_mask, fixture.affine, INTERIOR_MM)
    edge = out_mask & ~_interior(out_mask, fixture.affine, EDGE_MM)

    if is_multi_delay(fixture.metadata):
        out.update(
            _score_native_multi(
                fixture, preproc, cbf, att_masks_from=masks, subset=subset, aslprep_dir=aslprep_dir
            )
        )
        return out

    raw = fixture.asl_img.get_fdata()
    expected = expected_native(fixture, fwhm)
    a_masks = tissue_masks(fixture, extra=out_mask, reference=expected)
    # End to end: raw data through the documented model vs ASLPrep's CBF
    out['tier_a'] = ratio_stats(cbf, expected, a_masks['valid'])
    out['tier_a']['interior'] = ratio_stats(cbf, expected, a_masks['valid'] & interior)
    out['tier_a']['edge'] = ratio_stats(cbf, expected, a_masks['valid'] & edge)
    if preproc is not None:
        # Quantification given ASLPrep's preprocessed series: isolates the CBF model, the M0
        # handling and the calibration from motion-correction resampling.
        quant = _expected_from(fixture, preproc, fwhm)
        out['tier_a_quant'] = ratio_stats(cbf, quant, a_masks['valid'] & interior)
        out['hmc_deltam'] = _hmc_deltam(fixture, raw, preproc, a_masks['valid'] & interior)
    out['expected_ratio'] = {
        t: ratio_stats(expected, fixture.perfusion, masks[t])['median'] for t in ('GM', 'WM')
    }
    physical = _expected_from(
        fixture, raw, fwhm, alpha=tm.labeling_efficiency_physical(fixture.simulation)
    )
    out['physical'] = ratio_stats(cbf, physical, a_masks['valid'])
    return out


#: Distance (mm) from ASLPrep's brain-mask edge beyond which voxels count as interior, and
#: within which they count as edge, for the end-to-end Tier A breakdown (report-only).
INTERIOR_MM = 10.0
EDGE_MM = 5.0


def _interior(mask, affine, mm):
    from scipy.ndimage import distance_transform_edt

    sampling = np.sqrt((np.asarray(affine)[:3, :3] ** 2).sum(axis=0))
    return distance_transform_edt(mask, sampling=sampling) > mm


def _expected_from(fixture, asl, fwhm, alpha=None):
    return tm.expected_cbf(
        asl,
        fixture.context,
        fixture.metadata,
        fixture.affine,
        m0scan=fixture.m0scan,
        m0_metadata=fixture.m0_metadata,
        m0_scale=fixture.m0_scale,
        fwhm=fwhm,
        alpha=alpha,
    )


def _hmc_deltam(fixture, raw, preproc, mask):
    """How much motion correction changed delta-M (preprocessed vs raw), report-only."""
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = tm.deltam_mean(preproc, fixture.context) / tm.deltam_mean(raw, fixture.context)
    ratio = ratio[mask & np.isfinite(ratio)]
    if not ratio.size:
        return {'n': 0}
    return {
        'n': int(ratio.size),
        'median_dev': float(abs(np.median(ratio) - 1)),
        'p95_abs_dev': float(np.percentile(np.abs(ratio - 1), 95)),
    }


def _fit_subset(fixture, asl, voxels):
    return tm.expected_multi_delay(
        asl,
        fixture.context,
        fixture.metadata,
        fixture.affine,
        voxels,
        m0scan=fixture.m0scan,
        m0_metadata=fixture.m0_metadata,
        m0_scale=fixture.m0_scale,
    )


def _fit_stats(cbf, att, fit, voxels):
    i, j, k = voxels.T
    ok = np.isfinite(fit[:, 0]) & (fit[:, 0] > FLOOR) & np.isfinite(cbf[i, j, k])
    ratio = cbf[i, j, k][ok] / fit[ok, 0]
    return {
        'n': int(ok.sum()),
        'median_dev': float(abs(np.median(ratio) - 1)),
        'p95_abs_dev': float(np.percentile(np.abs(ratio - 1), 95)),
        'att_median_abs_diff': float(np.median(np.abs(att[i, j, k][ok] - fit[ok, 1]))),
    }, ok


def _score_native_multi(fixture, preproc, cbf, att_masks_from, subset, aslprep_dir):
    """Multi-delay Tier A on a voxel subset (the reference fit is per voxel and slow)."""
    masks = att_masks_from
    out = {}
    att_file = find_output(aslprep_dir, fixture.acq, 'att')
    if att_file is None:
        return {'problems': ['att not found']}
    att = on_grid(nb.load(att_file), fixture.affine, fixture.shape)
    voxels = np.argwhere(masks['valid'])
    voxels = voxels[:: max(1, len(voxels) // subset)][:subset]
    fit = _fit_subset(fixture, fixture.asl_img.get_fdata(), voxels)
    out['tier_a'], ok = _fit_stats(cbf, att, fit, voxels)
    if preproc is not None:
        out['tier_a_quant'], _ = _fit_stats(
            cbf, att, _fit_subset(fixture, preproc, voxels), voxels
        )
    i, j, k = voxels.T
    out['expected_ratio'] = {}
    for t in ('GM', 'WM'):
        sel = masks[t][i, j, k] & ok
        out['expected_ratio'][t] = (
            float(np.median(fit[sel, 0] / fixture.perfusion[i, j, k][sel])) if sel.any() else None
        )
    att_masks = {t: masks[t] & (fixture.att > 0) & (fixture.att < 100) for t in ('GM', 'WM')}
    out['tier_b_att'] = {
        t: {
            'median_abs_error': float(np.nanmedian(np.abs(att[m] - fixture.att[m]))),
            'r': float(np.corrcoef(np.nan_to_num(att[m]), fixture.att[m])[0, 1]),
        }
        for t, m in att_masks.items()
        if m.sum() > 2
    }
    return out


def score_frames(fixture, aslprep_dir):
    """Coregistration and head-motion errors from ASLPrep's transforms (Section 4.6)."""
    from nitransforms.linear import load

    out = {'problems': []}
    hmc = find_xfm(aslprep_dir, fixture.acq, 'from-orig_to-aslref')
    coreg = find_xfm(aslprep_dir, fixture.acq, 'from-aslref_to-T1w')
    if hmc is None or coreg is None:
        out['problems'].append('hmc or coreg transform not found')
        return out
    h = np.asarray(load(hmc, fmt='itk').matrix)
    h = h.reshape(-1, 4, 4)
    c = np.asarray(load(coreg, fmt='itk').matrix).reshape(4, 4)
    poses = fixture.poses()
    if len(h) != len(poses):
        out['problems'].append(f'{len(h)} HMC transforms for {len(poses)} volumes')
        return out

    _, points = fixture.voxel_centers(fixture.brain)
    a = np.stack([np.linalg.inv(h[v]) @ poses[v] for v in range(len(poses))])
    rel = [tg.rms_displacement(a[v] @ np.linalg.inv(a[0]), points) for v in range(len(a))]
    out['motion'] = {
        'rms_error_median_mm': float(np.median(rel)),
        'rms_error_max_mm': float(np.max(rel)),
    }
    native = [tg.rms_displacement(a[v], points) for v in range(len(a))]
    out['aslref_pose_rms_mm'] = float(np.max(native))
    out['native_comparable'] = bool(np.max(native) < NATIVE_TOLERANCE_MM)

    errors = [np.linalg.inv(c) @ a[v] @ np.linalg.inv(fixture.R) for v in range(len(a))]
    t1w_points = tg.apply_points(fixture.R, points)
    out['coreg'] = {
        'rot_deg': float(np.median([tg.rot_angle_deg(d) for d in errors])),
        'rms_mm': float(np.median([tg.rms_displacement(d, t1w_points) for d in errors])),
    }
    return out


def score_confounds(fixture, aslprep_dir):
    """Report-only motion parameter agreement: |r| per axis and mean framewise displacement."""
    path = find_output(aslprep_dir, fixture.acq, 'timeseries', desc='confounds', ext='.tsv')
    if path is None or fixture.motion is None:
        return {}
    conf = _read_tsv(path)
    out = {}
    for axis in ('trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'):
        est = np.array([row.get(axis, np.nan) for row in conf], dtype=float)
        true = np.array([row[axis] for row in fixture.motion], dtype=float)
        if len(est) == len(true) and np.std(est) > 0 and np.std(true) > 0:
            out[f'abs_r_{axis}'] = float(abs(np.corrcoef(est, true)[0, 1]))
    fd = np.array([row.get('framewise_displacement', np.nan) for row in conf], dtype=float)
    out['fd_mean'] = float(np.nanmean(fd))
    return out


def score_space(fixture, aslprep_dir, space, desc=None):
    """Tier B of a standard-space or T1w-space CBF map, sampled at truth voxel centres.

    T1w space: a phantom point ``p`` lies at ``R p``. MNI space (TemplateFlow phantoms): the
    phantom's world is the template's, so ``p`` is sampled as is.
    """
    path = find_output(aslprep_dir, fixture.acq, 'cbf', space=space, desc=desc)
    if path is None:
        return {'problems': [f'space-{space} cbf not found']}
    img = nb.load(path)
    masks = tissue_masks(fixture)

    def at(mask, shift_mm=0.0):
        ijk, points = fixture.voxel_centers(mask)
        points = points + [shift_mm, 0.0, 0.0]
        if space == 'T1w':
            points = tg.apply_points(fixture.R, points)
        return ijk, sample_at(img, points)

    out = {}
    for tissue in ('GM', 'WM'):
        ijk, values = at(masks[tissue])
        truth = fixture.perfusion[tuple(ijk.T)]
        finite = np.isfinite(values)
        out[tissue] = {
            'n': len(values),
            'finite': float(finite.mean()),
            'median': float(np.median(values[finite] / truth[finite])) if finite.any() else None,
        }

    # Alignment: correlation with ASLPrep's own native CBF at the same anatomical points
    # (independent of resampling blur), and with the truth. A 4 mm shift of the sampling points
    # gives the sensitivity reference.
    native_file = find_output(aslprep_dir, fixture.acq, 'cbf')
    if native_file is not None:
        native = on_grid(nb.load(native_file), fixture.affine, fixture.shape)
        ijk, values = at(masks['valid'])
        _, shifted = at(masks['valid'], shift_mm=4.0)
        ref, truth = native[tuple(ijk.T)], fixture.perfusion[tuple(ijk.T)]
        ok = np.isfinite(values) & np.isfinite(ref) & np.isfinite(shifted)
        out['r_native'] = float(np.corrcoef(values[ok], ref[ok])[0, 1])
        out['r_native_shifted_4mm'] = float(np.corrcoef(shifted[ok], ref[ok])[0, 1])
        out['r_truth'] = float(np.corrcoef(values[ok], truth[ok])[0, 1])
    return out


def score_clean(aslprep_dir):
    """Crash files and the report's error section."""
    aslprep_dir = Path(aslprep_dir)
    crashes = sorted(p.name for p in aslprep_dir.glob('sub-*/log/**/crash*'))
    report = aslprep_dir / 'sub-01.html'
    report_ok = report.exists() and 'No errors to report!' in report.read_text(errors='replace')
    return {'crash_files': crashes, 'report_ok': bool(report_ok)}


def score_run(fixture_path, aslprep_dir, acq, fwhm=5.0, spaces=(), extra_cbf=()):
    """Score one run (spec Section 4.7). Never raises for a missing output."""
    fixture = Fixture(fixture_path, acq)
    score = {'recipe': fixture.truth['recipe'], 'acq': acq}
    score['clean'] = score_clean(aslprep_dir)
    score['frames'] = score_frames(fixture, aslprep_dir)
    score['native'] = score_native(fixture, aslprep_dir, fwhm)
    score['confounds'] = score_confounds(fixture, aslprep_dir)
    for space in spaces:
        score[f'space-{space}'] = score_space(fixture, aslprep_dir, space)
    for desc in extra_cbf:
        path = find_output(aslprep_dir, acq, 'cbf', desc=desc)
        if path is None:
            score[f'desc-{desc}'] = {'problems': [f'desc-{desc} cbf not found']}
            continue
        cbf = on_grid(nb.load(path), fixture.affine, fixture.shape)
        masks = tissue_masks(fixture)
        score[f'desc-{desc}'] = {
            t: ratio_stats(cbf, fixture.perfusion, masks[t]) for t in ('GM', 'WM')
        }
    return score


def lookup(score, path):
    """``score['a']['b']`` for ``path`` ``('a', 'b')``; None when absent."""
    node = score
    for key in path:
        if not isinstance(node, dict) or key not in node:
            return None
        node = node[key]
    return node
