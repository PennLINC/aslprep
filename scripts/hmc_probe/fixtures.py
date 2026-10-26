#!/usr/bin/env python3
"""Probe-only simulated fixtures (Task 3a: the steps that consume the reference).

The probe's recipes live in ``scripts/hmc_probe/recipes/``, not with ASLPrep's test fixtures,
so they never enter the CI fixture digest. They are generated with ASLPrep's own fixture
machinery (``aslprep/tests/aslscan_fixtures.py``), pointed at that directory, and then given
the inputs aslscan does not produce, from each recipe's ``probe.toml``:

- ``[m0]``: the separate M0 moved by a known rigid ``pose`` (tx ty tz in mm, rx ry rz in
  degrees, about the field-of-view centre), and/or put on another grid (``zooms`` in mm, and
  ``flip_x`` to store the first axis reversed). A different grid is not a different pose: the
  two are set and recorded separately.
- ``[fieldmap]``: a BIDS fieldmap (``*_fieldmap.nii.gz`` in Hz and ``*_magnitude.nii.gz``)
  made from the phantom's B0 map, on its own grid (``zooms``) and at a known ``pose``, with
  ``IntendedFor`` the ASL run.

Moving a finished M0 image is an approximation of moving the head before the M0 is acquired
(as aslscan's own motion is for the ASL series). The truth goes to
``ground-truth/probe_extras.json``: every pose as a 4x4 RAS point mapping of the static phantom.

Phantoms: ``tfmni-b0`` is ASLPrep's ``tfmni`` phantom plus a synthetic B0 map (Gaussian
off-resonance blobs over the frontal sinus and the temporal lobes), so that aslscan simulates
EPI distortion.

Usage::

    python scripts/hmc_probe/fixtures.py --data-dir <dir> [--only rc_pcasl ...]
"""

import argparse
import hashlib
import json
import os
import shutil
import sys
import tomllib
from pathlib import Path

import nibabel as nb
import numpy as np
from scipy.ndimage import map_coordinates

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))

from aslprep.tests import aslscan_fixtures as af  # noqa: E402
from aslprep.tests import truth_geometry as tg  # noqa: E402

RECIPES_DIR = HERE / 'recipes'
#: The aslscan binary: the pinned CI build, unless ``ASLSCAN_BINARY`` names another (``st_spiral``
#: needs a newer build, for spiral gridding reconstruction). The fixture manifest records which.
BINARY = Path(
    os.environ.get(
        'ASLSCAN_BINARY', Path.home() / 'aslscan-build' / 'target' / 'release' / 'aslscan'
    )
)

#: Probe phantoms: base phantom plus a B0 map, as (centre in MNI mm, sigma in mm, peak Hz).
B0_BLOBS = (
    ((0.0, 50.0, -25.0), 14.0, 120.0),  # frontal sinus / orbitofrontal cortex
    ((-45.0, 0.0, -35.0), 12.0, -80.0),  # left temporal lobe (ear canal)
    ((45.0, 0.0, -35.0), 12.0, -80.0),  # right temporal lobe
)
#: Arterial segments for the macrovascular stress phantom, as (start, end) in MNI mm and a
#: radius in mm: the M1 segments of both middle cerebral arteries, the proximal anterior cerebral
#: arteries and their pericallosal course, and the P2 segments of the posterior cerebral arteries.
VESSELS = (
    ((-18.0, 4.0, -16.0), (-44.0, -2.0, -8.0), 2.5),
    ((18.0, 4.0, -16.0), (44.0, -2.0, -8.0), 2.5),
    ((-4.0, 14.0, -12.0), (0.0, 30.0, 18.0), 2.0),
    ((0.0, 30.0, 18.0), (0.0, -10.0, 34.0), 1.5),
    ((-10.0, -30.0, -14.0), (-26.0, -68.0, 0.0), 2.0),
    ((10.0, -30.0, -14.0), (26.0, -68.0, 0.0), 2.0),
)
#: Peak arterial blood volume fraction on a vessel's axis, and the arterial transit time (s)
#: there. Label is in the arteries only while aatt < t < aatt + tau; with rc_pcasl's 1.8 s label,
#: 1.8 s delay and 2D slices read 0-1.14 s later (t = 3.6-4.74 s), aatt must lie in 2.94-3.6 s:
#: a very slow-flow case (stenosis, older subjects), chosen so every slice shows the vessels.
VESSEL_ABV_PEAK = 0.3
VESSEL_AATT_S = 3.0
#: A territory with no perfusion: a sphere (MNI mm centre, radius mm) over the left middle
#: cerebral artery territory.
NO_PERFUSION = ((-40.0, -20.0, 20.0), 30.0)
PROBE_PHANTOMS = {
    'tfmni-b0': {'base': 'tfmni', 'b0_blobs': B0_BLOBS},
    # Task 4 stress phantoms for the D3 perfusion image (plan: "Stress cases").
    'tfmni-b0-vessels': {'base': 'tfmni-b0', 'vessels': VESSELS},
    'tfmni-b0-noperf': {'base': 'tfmni-b0', 'no_perfusion': NO_PERFUSION},
}


# ---------------------------------------------------------------------------------------------
# Pointing ASLPrep's fixture machinery at the probe's recipes and phantoms
# ---------------------------------------------------------------------------------------------
_orig_build_phantom = af.build_phantom
_orig_phantom_digest = af.phantom_digest


def _phantom_digest(name):
    if name in PROBE_PHANTOMS:
        spec = PROBE_PHANTOMS[name]
        if set(spec) == {'base', 'b0_blobs'}:  # the original form, so existing fixtures stay valid
            text = json.dumps({'base': _orig_phantom_digest(spec['base']), 'b0': spec['b0_blobs']})
            return hashlib.sha256(text.encode()).hexdigest()
        mods = {k: v for k, v in spec.items() if k != 'base'}
        if 'vessels' in spec:
            # Revision 2: arteries confined to the tissue labels (aslscan refuses background blood).
            mods['vessel_profile'] = [VESSEL_ABV_PEAK, VESSEL_AATT_S, 'rev3']
        text = json.dumps({'base': _phantom_digest(spec['base']), 'mods': mods})
        return hashlib.sha256(text.encode()).hexdigest()
    return _orig_phantom_digest(name)


def _segment_distance(xyz, p0, p1):
    """Distance (mm) from each point to the segment p0-p1."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0
    t = np.clip(((xyz - p0) @ d) / (d @ d), 0.0, 1.0)
    return np.linalg.norm(xyz - (p0 + t[:, None] * d), axis=1)


def _build_phantom(name, data_dir):
    if name not in PROBE_PHANTOMS:
        return _orig_build_phantom(name, data_dir)
    spec = PROBE_PHANTOMS[name]
    out = af._phantom_dir(name, data_dir)
    digest = _phantom_digest(name)
    prov = out / 'provenance.json'
    if prov.exists() and json.loads(prov.read_text()).get('digest') == digest:
        return out
    base = _build_phantom(spec['base'], data_dir)  # a probe phantom may build on another
    shutil.rmtree(out, ignore_errors=True)
    shutil.copytree(base, out)
    (out / 'provenance.json').unlink(missing_ok=True)
    img = nb.load(out / 'dseg.nii.gz')
    ijk = np.indices(img.shape).reshape(3, -1).T
    xyz = tg.apply_points(img.affine, ijk)
    if 'b0_blobs' in spec:
        field = np.zeros(len(xyz), np.float32)
        for centre, sigma, peak in spec['b0_blobs']:
            d2 = ((xyz - np.asarray(centre)) ** 2).sum(axis=1)
            field += peak * np.exp(-d2 / (2 * sigma**2))
        nb.Nifti1Image(field.reshape(img.shape), img.affine).to_filename(out / 'fieldmap.nii.gz')
        (out / 'fieldmap.json').write_text(json.dumps({'Units': 'Hz'}) + '\n')
    if 'vessels' in spec:
        # Arterial blood volume: a Gaussian profile across each vessel (sigma = radius), the
        # largest over vessels; aslscan switches its arterial compartment on from these maps.
        abv = np.zeros(len(xyz), np.float32)
        for p0, p1, radius in spec['vessels']:
            dist = _segment_distance(xyz, p0, p1)
            abv = np.maximum(abv, VESSEL_ABV_PEAK * np.exp(-(dist**2) / (2 * radius**2)))
        abv[abv < 1e-3] = 0.0
        # aslscan refuses blood in background voxels: keep the arteries inside the tissue labels.
        abv[np.asanyarray(img.dataobj).ravel() == 0] = 0.0
        aatt = np.full(len(xyz), VESSEL_AATT_S, np.float32)
        for fname, data, units in (('abv', abv, 'fraction'), ('aatt', aatt, 's')):
            nb.Nifti1Image(data.reshape(img.shape), img.affine).to_filename(
                out / f'{fname}.nii.gz'
            )
            (out / f'{fname}.json').write_text(json.dumps({'Units': units}) + '\n')
    if 'no_perfusion' in spec:
        centre, radius = spec['no_perfusion']
        perf_img = nb.load(out / 'perfusion.nii.gz')
        perf = np.asanyarray(perf_img.dataobj, np.float32).copy()
        inside = (((xyz - np.asarray(centre)) ** 2).sum(axis=1) <= radius**2).reshape(img.shape)
        perf[inside] = 0.0
        nb.Nifti1Image(perf, perf_img.affine, perf_img.header).to_filename(
            out / 'perfusion.nii.gz'
        )
    prov.write_text(json.dumps({'name': name, 'digest': digest, 'base': spec['base']}) + '\n')
    return out


def use_probe_recipes():
    af.RECIPES_DIR = RECIPES_DIR
    af.RECIPES = af._discover_recipes()
    af.build_phantom = _build_phantom
    af.phantom_digest = _phantom_digest
    if 'ASLSCAN_BINARY' in os.environ:
        # An explicitly named, non-pinned build: accept the revisions its own provenance record
        # states (check_aslscan still verifies the binary's hash), so the fixture's manifest
        # records the build that actually made it.
        stamp = json.loads((BINARY.parent / af.STAMP_NAME).read_text())
        keys = ('aslscan', 'mrsim-acq', 'toolchain', 'features', 'cargo_lock')
        af._expected_stamp = lambda: {k: stamp[k] for k in keys}


# ---------------------------------------------------------------------------------------------
# Extras
# ---------------------------------------------------------------------------------------------
def rigid(params, centre):
    """A rigid RAS point mapping (mm, degrees) about ``centre``."""
    t, r = params[:3], np.radians(params[3:])
    return tg.pose_matrix(t, r, centre)


def resample(data, affine, out_affine, out_shape, pose=None, order=3):
    """``data`` (on ``affine``) seen on a new grid after moving the head by ``pose``.

    A moved head shows the static image at ``pose^-1 x``.
    """
    pose = np.eye(4) if pose is None else pose
    ijk = np.indices(out_shape).reshape(3, -1).T
    world = tg.apply_points(out_affine, ijk)
    src = tg.apply_points(np.linalg.inv(affine) @ np.linalg.inv(pose), world)
    vals = map_coordinates(np.asarray(data, np.float64), src.T, order=order, mode='constant')
    return vals.reshape(out_shape).astype(np.float32)


def new_grid(affine, shape, zooms, flip_x=False):
    """A grid covering the same field of view, with other voxel sizes and/or a reversed x axis."""
    zooms = np.asarray(zooms, float)
    old_zooms = np.sqrt((affine[:3, :3] ** 2).sum(axis=0))
    extent = np.asarray(shape) * old_zooms
    new_shape = np.maximum(np.round(extent / zooms).astype(int), 1)
    centre = tg.fov_center(affine, shape)
    direction = affine[:3, :3] / old_zooms
    lin = direction * zooms
    if flip_x:
        lin[:, 0] *= -1
    out = np.eye(4)
    out[:3, :3] = lin
    out[:3, 3] = centre - lin @ ((new_shape - 1) / 2)
    return out, tuple(int(s) for s in new_shape)


def add_extras(fixture_dir, recipe_dir):
    """Apply ``probe.toml`` to a generated fixture."""
    cfg_file = recipe_dir / 'probe.toml'
    if not cfg_file.exists():
        return {}
    cfg = tomllib.loads(cfg_file.read_text())
    perf = fixture_dir / 'sub-01' / 'perf'
    asl = next(perf.glob('*_asl.nii.gz'))
    stem = asl.name[: -len('_asl.nii.gz')]
    asl_img = nb.load(asl)
    centre = tg.fov_center(asl_img.affine, asl_img.shape[:3])
    truth = {'centre': centre.tolist()}

    if 'm0' in cfg:
        m0_file = perf / f'{stem}_m0scan.nii.gz'
        m0_img = nb.load(m0_file)
        data = np.asanyarray(m0_img.dataobj, np.float32)
        pose = rigid(cfg['m0'].get('pose', [0.0] * 6), centre)
        affine, shape = m0_img.affine, m0_img.shape[:3]
        if 'zooms' in cfg['m0'] or cfg['m0'].get('flip_x'):
            zooms = cfg['m0'].get('zooms', np.sqrt((affine[:3, :3] ** 2).sum(axis=0)))
            affine, shape = new_grid(m0_img.affine, shape, zooms, cfg['m0'].get('flip_x', False))
        vols = data[..., None] if data.ndim == 3 else data
        out = np.stack(
            [
                resample(vols[..., k], m0_img.affine, affine, shape, pose)
                for k in range(vols.shape[-1])
            ],
            axis=-1,
        )
        out = out[..., 0] if data.ndim == 3 else out
        nb.Nifti1Image(out, affine).to_filename(m0_file)
        truth['m0'] = {'pose': pose.tolist(), 'affine': affine.tolist(), 'shape': list(shape)}

    if 'fieldmap' in cfg:
        phantom = af.build_phantom(af.load_recipe(recipe_dir.name).phantom, fixture_dir.parents[1])
        f_img = nb.load(phantom / 'fieldmap.nii.gz')
        mag_img = nb.load(phantom / 'M0map.nii.gz')
        pose = rigid(cfg['fieldmap'].get('pose', [0.0] * 6), centre)
        affine, shape = new_grid(asl_img.affine, asl_img.shape[:3], cfg['fieldmap']['zooms'])
        fmap_dir = fixture_dir / 'sub-01' / 'fmap'
        fmap_dir.mkdir(exist_ok=True)
        fstem = f'sub-01_acq-{cfg["fieldmap"].get("acq", "b0")}'
        field = resample(np.asanyarray(f_img.dataobj), f_img.affine, affine, shape, pose, order=1)
        mag = resample(
            np.asanyarray(mag_img.dataobj), mag_img.affine, affine, shape, pose, order=1
        )
        nb.Nifti1Image(field, affine).to_filename(fmap_dir / f'{fstem}_fieldmap.nii.gz')
        nb.Nifti1Image(mag, affine).to_filename(fmap_dir / f'{fstem}_magnitude.nii.gz')
        meta = {'Units': 'Hz', 'IntendedFor': [f'bids::sub-01/perf/{asl.name}']}
        (fmap_dir / f'{fstem}_fieldmap.json').write_text(json.dumps(meta, indent=2) + '\n')
        (fmap_dir / f'{fstem}_magnitude.json').write_text(json.dumps({}) + '\n')
        truth['fieldmap'] = {
            'pose': pose.tolist(),
            'affine': affine.tolist(),
            'shape': list(shape),
        }

    (perf / 'ground-truth' / 'probe_extras.json').write_text(json.dumps(truth, indent=2) + '\n')
    # The manifest's file hashes no longer describe the fixture; record that it was changed.
    manifest_file = fixture_dir / 'manifest.json'
    manifest = json.loads(manifest_file.read_text())
    manifest['probe_extras'] = cfg
    manifest['files'] = af._file_hashes(fixture_dir)
    manifest_file.write_text(json.dumps(manifest, indent=1) + '\n')
    return truth


def seed_variants(names, seeds, staging):
    """Copy recipes (ASLPrep's or the probe's) into ``staging`` as ``<name>__s<seed>``.

    Each copy has its overlay ``seed`` replaced, so it is another noise realization of the
    same protocol. The original recipe directory is not changed.
    """
    import re

    staging.mkdir(parents=True, exist_ok=True)
    out = []
    for name in names:
        src = RECIPES_DIR / name
        if not src.exists():
            src = af.TESTS_DIR / 'data' / 'aslscan' / name
        for seed in seeds:
            dest = staging / f'{name}__s{seed}'
            shutil.rmtree(dest, ignore_errors=True)
            shutil.copytree(src, dest)
            overlay = (dest / 'overlay.toml').read_text()
            overlay, n = re.subn(r'(?m)^seed\s*=.*$', f'seed = {seed}', overlay)
            if n != 1:
                raise ValueError(f'{name}: expected one seed line in overlay.toml')
            (dest / 'overlay.toml').write_text(overlay)
            out.append(dest.name)
    return out


def _robust_snr(raw_noisy, raw_clean, context, gm_mask):
    """As ``aslscan_fixtures._snr``, with the signal a per-voxel median over pair estimates.

    The clean series carries the recipe's motion, so its plain mean delta-M includes motion
    contamination; with 2 mm / 2 degree spikes in a 10-pair series it is negative and the
    calibration fails. The median ignores the spiked pairs and equals the mean when the head
    is still.
    """
    noisy = af.deltam_estimates(nb.load(af._raw_paths(raw_noisy)['asl']).get_fdata(), context)
    clean = af.deltam_estimates(nb.load(af._raw_paths(raw_clean)['asl']).get_fdata(), context)
    signal = float(np.mean(np.median(clean, axis=-1)[gm_mask]))
    noise = float(np.std((noisy - clean)[gm_mask]))
    return signal, noise


def generate(names, data_dir, threads=8, recipes_dir=None):
    use_probe_recipes()
    af._snr = _robust_snr  # probe-only: the CI fixtures keep the plain mean
    if recipes_dir is not None:
        af.RECIPES_DIR = recipes_dir
        af.RECIPES = af._discover_recipes()
    for name in names:
        out = af.generate(name, data_dir, BINARY, threads=threads)
        add_extras(out, af.RECIPES_DIR / name)
        print(f'{name} -> {out}', flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--only', nargs='*')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument(
        '--seeds',
        type=int,
        nargs='*',
        help='Generate these noise realizations of each recipe (as <name>__s<seed>).',
    )
    args = parser.parse_args(argv)
    use_probe_recipes()
    names = args.only or sorted(af.RECIPES)
    if args.seeds:
        staging = args.data_dir / 'recipes-staging'
        names = seed_variants(names, args.seeds, staging)
        generate(names, args.data_dir, args.threads, recipes_dir=staging)
    else:
        generate(names, args.data_dir, args.threads)


if __name__ == '__main__':
    main()
