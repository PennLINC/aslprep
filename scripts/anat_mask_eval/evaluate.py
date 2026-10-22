#!/usr/bin/env python3
"""Compare ASL brain masks propagated from the T1w with SynthStrip and legacy masks.

For every job in ``jobs.json`` (see ``jobs.py``), with the SynthStrip and legacy masks from
``aslmasks.py``, this builds the anatomical masks in aslref space and scores all of them.

Anatomical masks (the T1w brain mask pulled onto the aslref grid through the coregistration):

- ``anat_nn``: nearest neighbour at the voxel centre (what PennLINC/aslprep#601 did, with ANTs'
  ``GenericLabel``);
- ``anat_pv50`` / ``anat_pv10``: the fraction of each aslref voxel inside the T1w mask
  (``SUPERSAMPLE``\\ :sup:`3` points per voxel), thresholded at 0.5 / 0.1;
- ``hybrid``: ``anat_pv50``, unless its Dice with SynthStrip is below ``GUARD_DICE`` (a sign
  that the coregistration failed), then SynthStrip's;
- ``hybrid_vol``: the same with a volume check (:func:`guarded_volume`), which keeps the
  anatomical mask when SynthStrip's is the smaller one.

Real data (no ground truth):

- agreement between masks (Dice, volume);
- perfusion edge metrics on the unmasked mean perfusion image ``P``: the share of the summed
  perfusion signal outside the mask (noise sums to about zero, so this is the brain signal,
  including the readout's blurring, that the mask leaves out), and the mean of ``P`` in a
  one-voxel shell just outside (``outer``) and just inside (``inner``) the mask, relative to the
  mean inside the mask. A mask that is too tight has a high ``outer``; one that is too loose has
  a low ``inner``;
- sensitivity of the anatomical mask to the registration: Dice with the masks from the other
  registration candidates of the coreg probe, and with the production transform perturbed by
  rigid translations of known size.

Simulated data (aslscan fixtures): Dice and errors against the true brain (tissue partial
volumes summing to 0.5 or more) in the aslref frame, the share of true grey matter each mask
keeps, and the anatomical mask through the true transform (the method's ceiling).

Usage::

    python evaluate.py --jobs jobs.json --masks masks/ --out results/
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
from scipy import ndimage as ndi

SUPERSAMPLE = 3
PERTURB_MM = (1.0, 2.0, 3.0, 4.0, 6.0)
PERTURB_DRAWS = 8
#: Registration sources of the coreg probe; ``conventional/default`` is ASLPrep's own.
SOURCES = ('conventional', 'P-mean', 'P-robust', 'P-winsor', 'P-masked')
SOURCE_IMAGE = {'P-masked': 'P-winsor'}
METHODS = ('synthstrip', 'legacy', 'anat_nn', 'anat_pv50', 'anat_pv10', 'hybrid', 'hybrid_vol')
#: Below this Dice with SynthStrip, the hybrid mask takes SynthStrip's (failed coregistration).
GUARD_DICE = 0.85
#: The volume-aware guard (see :func:`guarded_volume`); set after looking at these data.
GUARD_DICE_VOLUME = 0.90
GUARD_VOLUME = 0.95


# ---------------------------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------------------------
def fsl_to_world(img):
    """FSL scaled-voxel coordinates -> the image's world (RAS mm) coordinates."""
    zooms = np.array(img.header.get_zooms()[:3], float)
    d = np.diag(np.r_[zooms, 1.0])
    if np.linalg.det(img.affine[:3, :3]) > 0:  # FSL flips x
        d[0, 0] = -zooms[0]
        d[0, 3] = (img.shape[0] - 1) * zooms[0]
    return img.affine @ np.linalg.inv(d)


def fsl_pull(mat_file, src_img, t1_img):
    """An FSL matrix (source -> T1w) as a RAS pull from T1w points to source points."""
    push = fsl_to_world(t1_img) @ np.loadtxt(mat_file) @ np.linalg.inv(fsl_to_world(src_img))
    return np.linalg.inv(push)


def subvoxel_points(ref_img, n=SUPERSAMPLE):
    """World coordinates of ``n``\\ :sup:`3` points inside every voxel: (n**3, X, Y, Z, 3)."""
    offsets = (np.arange(n) + 0.5) / n - 0.5
    shape = ref_img.shape[:3]
    grid = np.stack(np.meshgrid(*[np.arange(s) for s in shape], indexing='ij'), -1).astype(float)
    pts = []
    for dx in offsets:
        for dy in offsets:
            for dz in offsets:
                ijk = grid + np.array([dx, dy, dz])
                pts.append(nb.affines.apply_affine(ref_img.affine, ijk).astype(np.float32))
    return np.stack(pts)


def sample(img_data, img_affine, world, order):
    """Sample an image at world points (any leading shape), zero outside."""
    vox = nb.affines.apply_affine(np.linalg.inv(img_affine), world.reshape(-1, 3))
    out = ndi.map_coordinates(img_data, vox.T, order=order, mode='constant', cval=0.0)
    return out.reshape(world.shape[:-1])


class AnatSampler:
    """Pulls a T1w mask onto an aslref grid, through any T1w -> aslref pull."""

    def __init__(self, t1_mask_img, ref_img, n=SUPERSAMPLE):
        self.data = (np.asanyarray(t1_mask_img.dataobj) > 0).astype(np.float32)
        self.t1_affine = t1_mask_img.affine
        self.ref_affine = ref_img.affine
        self.shape = ref_img.shape[:3]
        ijk = np.indices(self.shape).reshape(3, -1).T.astype(np.float32)
        offsets = (np.arange(n) + 0.5) / n - 0.5
        offsets = np.stack(np.meshgrid(offsets, offsets, offsets, indexing='ij'), -1)
        self.centres = ijk
        self.sub = ijk[None] + offsets.reshape(-1, 1, 3).astype(np.float32)  # (n**3, V, 3)

    def _vox(self, pull):
        """aslref voxel indices -> T1w voxel indices."""
        return np.linalg.inv(self.t1_affine) @ np.linalg.inv(pull) @ self.ref_affine

    def _at(self, m, pts):
        vox = pts @ m[:3, :3].T.astype(np.float32) + m[:3, 3].astype(np.float32)
        flat = vox.reshape(-1, 3).T
        out = ndi.map_coordinates(self.data, flat, order=0, mode='constant', cval=0.0)
        return out.reshape(vox.shape[:-1])

    def frac(self, pull, voxels=None):
        """Fraction of each aslref voxel inside the T1w mask (flat ``voxels`` only, if given)."""
        pts = self.sub if voxels is None else self.sub[:, voxels]
        return self._at(self._vox(pull), pts).mean(0)

    def masks(self, pull):
        m = self._vox(pull)
        nn = (self._at(m, self.centres) > 0.5).reshape(self.shape)
        frac = self.frac(pull).reshape(self.shape)
        return {'anat_nn': nn, 'anat_pv50': frac >= 0.5, 'anat_pv10': frac >= 0.1}, frac


def anat_masks(t1_mask_img, ref_img, pull):
    """The anatomical masks on the aslref grid; ``pull`` maps T1w points to aslref points."""
    return AnatSampler(t1_mask_img, ref_img).masks(pull)


def boundary_band(mask, voxel_mm, mm):
    """Flat indices of the voxels within ``mm`` (plus one voxel) of the mask's boundary."""
    edge = mask ^ ndi.binary_erosion(mask)
    steps = int(np.ceil(mm / voxel_mm)) + 1
    return np.flatnonzero(ndi.binary_dilation(edge, iterations=steps))


def translation(mm, rng):
    t = np.eye(4)
    v = rng.normal(size=3)
    t[:3, 3] = mm * v / np.linalg.norm(v)
    return t


# ---------------------------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------------------------
def dice(a, b):
    s = a.sum() + b.sum()
    return float(2 * (a & b).sum() / s) if s else np.nan


def shells(mask):
    inner = mask & ~ndi.binary_erosion(mask)
    outer = ndi.binary_dilation(mask) & ~mask
    return inner, outer


def perfusion_metrics(mask, perf):
    inner, outer = shells(mask)
    inside = perf[mask].mean()
    total = perf.sum()
    return {
        'perf_out_frac': float(perf[~mask].sum() / total),
        'outer_ratio': float(perf[outer].mean() / inside),
        'inner_ratio': float(perf[inner].mean() / inside),
    }


def significant_perfusion(perf, union, z=3.0):
    """Voxels whose smoothed perfusion is clearly above the noise.

    The noise (median absolute deviation) is measured far from the brain: outside ``union``
    dilated by 5 voxels.
    """
    smooth = ndi.gaussian_filter(perf, 1.0)
    far = ~ndi.binary_dilation(union, iterations=5)
    bg = smooth[far & np.isfinite(smooth)]
    if bg.size < 100:
        return np.zeros(perf.shape, bool)
    centre = np.median(bg)
    mad = 1.4826 * np.median(np.abs(bg - centre))
    return smooth > centre + z * mad


def guarded(masks):
    """The anatomical mask, unless it disagrees with SynthStrip (then SynthStrip's)."""
    if dice(masks['anat_pv50'], masks['synthstrip']) >= GUARD_DICE:
        return masks['anat_pv50']
    return masks['synthstrip']


def guarded_volume(masks):
    """Like :func:`guarded`, but only when SynthStrip's mask is not the smaller one.

    A SynthStrip failure removes brain, so its mask is smaller than the anatomical one; a failed
    coregistration moves the anatomical mask off the brain (and partly out of the field of
    view) without shrinking SynthStrip's. Below ``GUARD_DICE_VOLUME`` Dice, the anatomical mask
    is kept only when SynthStrip's has less than ``GUARD_VOLUME`` of its volume.
    """
    anat, ss = masks['anat_pv50'], masks['synthstrip']
    if dice(anat, ss) >= GUARD_DICE_VOLUME or ss.sum() < GUARD_VOLUME * anat.sum():
        return anat
    return ss


def load_mask(path, like=None):
    img = nb.load(path)
    data = np.asanyarray(img.dataobj) > 0
    if like is not None and data.shape != like.shape[:3]:
        raise ValueError(f'{path}: shape {data.shape} != {like.shape[:3]}')
    return data


# ---------------------------------------------------------------------------------------------
# Real data
# ---------------------------------------------------------------------------------------------
def evaluate_real(job, mask_dir, out_dir, rng):
    d = Path(job['dir'])
    ref = nb.load(job['ref'])
    t1 = nb.load(d / 't1w_brain.nii.gz')
    t1_mask = nb.load(job['t1w_mask'])
    vox_ml = float(np.prod(ref.header.get_zooms()[:3])) / 1000

    def pull_for(source, sign='default'):
        image = SOURCE_IMAGE.get(source, source)
        src = nb.load(d / f'{image}.nii.gz')
        kind = 'init' if job['fallback'].get(f'{source}/{sign}') else f'bbr_{sign}'
        mat = d / f'{source}_{kind}.mat'
        return fsl_pull(mat, src, t1) if mat.exists() else None

    pull = pull_for('conventional')
    masks = {
        'production': load_mask(job['prod_mask'], ref),
        'synthstrip': load_mask(mask_dir / job['name'] / 'synthstrip.nii.gz', ref),
        'legacy': load_mask(mask_dir / job['name'] / 'legacy.nii.gz', ref),
    }
    sampler = AnatSampler(t1_mask, ref)
    anat, frac = sampler.masks(pull)
    masks.update(anat)
    masks['hybrid'] = guarded(masks)
    masks['hybrid_vol'] = guarded_volume(masks)

    perf = np.asanyarray(nb.load(job['perfusion']).dataobj, np.float32)
    union = masks['synthstrip'] | masks['anat_pv50']
    sign = 1.0 if np.median(perf[union]) >= 0 else -1.0
    perf = sign * perf
    signif = significant_perfusion(perf, union)

    rows = []
    for name, m in masks.items():
        row = {
            'name': job['name'],
            'method': name,
            'vol_ml': float(m.sum() * vox_ml),
            'dice_anat_pv50': dice(m, masks['anat_pv50']),
            'dice_synthstrip': dice(m, masks['synthstrip']),
            'dice_legacy': dice(m, masks['legacy']),
        }
        row.update(perfusion_metrics(m, perf))
        row['signif_kept'] = float((signif & m).sum() / signif.sum()) if signif.any() else np.nan
        rows.append(row)

    # Registration sensitivity: other registration candidates, then known perturbations.
    alt = {}
    for source in SOURCES[1:]:
        p = pull_for(source)
        if p is not None:
            alt[source] = dice(sampler.masks(p)[0]['anat_pv50'], masks['anat_pv50'])
    perturb = []
    base = masks['anat_pv50'].ravel()
    # A translation of at most max(PERTURB_MM) changes voxels near the boundary only.
    band = boundary_band(masks['anat_pv50'], min(ref.header.get_zooms()[:3]), max(PERTURB_MM))
    for mm in PERTURB_MM:
        for _ in range(PERTURB_DRAWS):
            m = base.copy()
            m[band] = sampler.frac(pull @ translation(mm, rng), band) >= 0.5
            perturb.append({'name': job['name'], 'shift_mm': mm, 'dice': dice(m, base)})

    run_dir = out_dir / 'masks' / job['name']
    run_dir.mkdir(parents=True, exist_ok=True)
    for name in ('anat_nn', 'anat_pv50', 'anat_pv10'):
        nb.Nifti1Image(masks[name].astype(np.uint8), ref.affine).to_filename(
            run_dir / f'{name}.nii.gz'
        )
    nb.Nifti1Image(frac.astype(np.float32), ref.affine).to_filename(run_dir / 'anat_frac.nii.gz')

    info = {
        'name': job['name'],
        'dataset': job['dataset'],
        'family': job['family'],
        'acquisition': job['acquisition'],
        'sdc': job['sdc'],
        'voxel_mm': float(max(job['zooms_mm'])),
        # Identical references (several runs can share one) count once in the summaries.
        'ref_sha1': hashlib.sha1(Path(job['ref']).read_bytes(), usedforsecurity=False).hexdigest(),
        'perfusion_sign': sign,
        'perfusion_eligible': job['perfusion_eligible'],
        'coreg_fallback': bool(job['fallback'].get('conventional/default')),
        # Share of the T1w brain mask's volume that lands inside the ASL field of view.
        'anat_in_fov': float(
            frac.sum()
            * vox_ml
            / (
                (np.asanyarray(t1_mask.dataobj) > 0).sum()
                * np.prod(t1_mask.header.get_zooms()[:3])
                / 1000
            )
        ),
        'alt_dice_min': min(alt.values()) if alt else np.nan,
        **{f'alt_dice_{k}': v for k, v in alt.items()},
    }
    return rows, perturb, info


# ---------------------------------------------------------------------------------------------
# Simulated data
# ---------------------------------------------------------------------------------------------
def evaluate_sim(job, mask_dir, out_dir):
    from nitransforms.linear import load

    from aslprep.tests import truth_scoring as ts

    fixture = ts.Fixture(job['fixture'], job['acq'])
    ref = nb.load(job['ref'])
    t1_mask = nb.load(job['t1w_mask'])
    vox_ml = float(np.prod(ref.header.get_zooms()[:3])) / 1000
    perf_dir = Path(job['aslprep_dir']) / 'sub-01' / 'perf'

    h = np.asarray(
        load(ts.find_xfm(perf_dir.parents[1], job['acq'], 'from-orig_to-aslref'), fmt='itk').matrix
    )
    h = h.reshape(-1, 4, 4)
    c = np.asarray(
        load(ts.find_xfm(perf_dir.parents[1], job['acq'], 'from-aslref_to-T1w'), fmt='itk').matrix
    )
    c = c.reshape(4, 4)
    poses = fixture.poses()
    # The aslref's pose (phantom point -> aslref point), from the volume closest to the median.
    a = np.stack([np.linalg.inv(h[v]) @ poses[v] for v in range(len(poses))])
    _, points = fixture.voxel_centers(fixture.brain)
    cost = [
        np.median(
            [ts.tg.rms_displacement(a[v] @ np.linalg.inv(a[w]), points) for w in range(len(a))]
        )
        for v in range(len(a))
    ]
    pose = a[int(np.argmin(cost))]

    # Truth on the aslref grid: aslref point y shows phantom point pose^-1 y.
    sub = subvoxel_points(ref)
    to_phantom = nb.affines.apply_affine(np.linalg.inv(pose), sub)
    pv = {
        t: sample(fixture.pv[t].astype(np.float32), fixture.affine, to_phantom, 1).mean(0)
        for t in fixture.pv
    }
    brain_frac = pv['GM'] + pv['WM'] + pv['CSF']
    truth = brain_frac >= 0.5

    masks = {
        'production': load_mask(job['prod_mask'], ref),
        'synthstrip': load_mask(mask_dir / job['name'] / 'synthstrip.nii.gz', ref),
        'legacy': load_mask(mask_dir / job['name'] / 'legacy.nii.gz', ref),
    }
    anat, _ = anat_masks(t1_mask, ref, c)
    masks.update(anat)
    masks['hybrid'] = guarded(masks)
    masks['hybrid_vol'] = guarded_volume(masks)
    # Ceiling: the T1w mask through the true transform (T1w point R p -> aslref point pose p).
    true_pull = pose @ np.linalg.inv(fixture.R)
    masks['anat_truexfm'] = anat_masks(t1_mask, ref, true_pull)[0]['anat_pv50']

    t1w_points = ts.tg.apply_points(fixture.R, points)
    coreg_err = ts.tg.rms_displacement(np.linalg.inv(c) @ true_pull, t1w_points)

    gm_total = pv['GM'].sum()
    rows = []
    for name, m in masks.items():
        rows.append(
            {
                'name': job['name'],
                'method': name,
                'vol_ml': float(m.sum() * vox_ml),
                'truth_vol_ml': float(truth.sum() * vox_ml),
                'dice_truth': dice(m, truth),
                'false_pos_ml': float((m & ~truth).sum() * vox_ml),
                'false_neg_ml': float((~m & truth).sum() * vox_ml),
                'brain_pv_kept': float(brain_frac[m].sum() / brain_frac.sum()),
                'gm_kept': float(pv['GM'][m].sum() / gm_total),
                'nonbrain_pv_ml': float((1 - brain_frac[m]).clip(0).sum() * vox_ml),
                'coreg_err_mm': float(coreg_err),
                'sdc': job['sdc'],
            }
        )
    run_dir = out_dir / 'masks' / job['name']
    run_dir.mkdir(parents=True, exist_ok=True)
    for name in ('anat_nn', 'anat_pv50', 'anat_pv10', 'anat_truexfm'):
        nb.Nifti1Image(masks[name].astype(np.uint8), ref.affine).to_filename(
            run_dir / f'{name}.nii.gz'
        )
    nb.Nifti1Image(truth.astype(np.uint8), ref.affine).to_filename(run_dir / 'truth.nii.gz')
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--jobs', type=Path, required=True)
    parser.add_argument('--masks', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--only', nargs='*', help='Job names to run (default: all).')
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261022)
    jobs = json.loads(args.jobs.read_text())
    if args.only:
        jobs = [j for j in jobs if j['name'] in args.only]

    real_rows, perturb_rows, infos, sim_rows, problems = [], [], [], [], []
    for job in jobs:
        if (
            not (args.masks / job['name'] / 'synthstrip.nii.gz').exists()
            or not (args.masks / job['name'] / 'legacy.nii.gz').exists()
        ):
            problems.append({'name': job['name'], 'problem': 'masks missing'})
            continue
        try:
            if job['kind'] == 'real':
                rows, perturb, info = evaluate_real(job, args.masks, args.out, rng)
                real_rows += rows
                perturb_rows += perturb
                infos.append(info)
            else:
                sim_rows += evaluate_sim(job, args.masks, args.out)
        except Exception as err:  # noqa: BLE001 - record, keep going
            problems.append({'name': job['name'], 'problem': f'{type(err).__name__}: {err}'})
        print(job['name'], flush=True, file=sys.stderr)

    pd.DataFrame(real_rows).to_csv(args.out / 'real_masks.tsv', sep='\t', index=False)
    pd.DataFrame(perturb_rows).to_csv(args.out / 'real_perturb.tsv', sep='\t', index=False)
    pd.DataFrame(infos).to_csv(args.out / 'real_runs.tsv', sep='\t', index=False)
    pd.DataFrame(sim_rows).to_csv(args.out / 'sim_masks.tsv', sep='\t', index=False)
    pd.DataFrame(problems, columns=['name', 'problem']).to_csv(
        args.out / 'problems.tsv', sep='\t', index=False
    )
    print(
        f'{len(infos)} real runs, {len(sim_rows) and len({r["name"] for r in sim_rows})} sims, '
        f'{len(problems)} problems'
    )


if __name__ == '__main__':
    main()
