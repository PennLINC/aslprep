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
BINARY = Path.home() / 'aslscan-build' / 'target' / 'release' / 'aslscan'

#: Probe phantoms: base phantom plus a B0 map, as (centre in MNI mm, sigma in mm, peak Hz).
B0_BLOBS = (
    ((0.0, 50.0, -25.0), 14.0, 120.0),  # frontal sinus / orbitofrontal cortex
    ((-45.0, 0.0, -35.0), 12.0, -80.0),  # left temporal lobe (ear canal)
    ((45.0, 0.0, -35.0), 12.0, -80.0),  # right temporal lobe
)
PROBE_PHANTOMS = {'tfmni-b0': {'base': 'tfmni', 'b0_blobs': B0_BLOBS}}


# ---------------------------------------------------------------------------------------------
# Pointing ASLPrep's fixture machinery at the probe's recipes and phantoms
# ---------------------------------------------------------------------------------------------
_orig_build_phantom = af.build_phantom
_orig_phantom_digest = af.phantom_digest


def _phantom_digest(name):
    if name in PROBE_PHANTOMS:
        spec = PROBE_PHANTOMS[name]
        text = json.dumps({'base': _orig_phantom_digest(spec['base']), 'b0': spec['b0_blobs']})
        return hashlib.sha256(text.encode()).hexdigest()
    return _orig_phantom_digest(name)


def _build_phantom(name, data_dir):
    if name not in PROBE_PHANTOMS:
        return _orig_build_phantom(name, data_dir)
    spec = PROBE_PHANTOMS[name]
    out = af._phantom_dir(name, data_dir)
    digest = _phantom_digest(name)
    prov = out / 'provenance.json'
    if prov.exists() and json.loads(prov.read_text()).get('digest') == digest:
        return out
    base = _orig_build_phantom(spec['base'], data_dir)
    shutil.rmtree(out, ignore_errors=True)
    shutil.copytree(base, out)
    img = nb.load(out / 'dseg.nii.gz')
    ijk = np.indices(img.shape).reshape(3, -1).T
    xyz = tg.apply_points(img.affine, ijk)
    field = np.zeros(len(xyz), np.float32)
    for centre, sigma, peak in spec['b0_blobs']:
        d2 = ((xyz - np.asarray(centre)) ** 2).sum(axis=1)
        field += peak * np.exp(-d2 / (2 * sigma**2))
    nb.Nifti1Image(field.reshape(img.shape), img.affine).to_filename(out / 'fieldmap.nii.gz')
    (out / 'fieldmap.json').write_text(json.dumps({'Units': 'Hz'}) + '\n')
    prov.write_text(json.dumps({'name': name, 'digest': digest, 'base': spec['base']}) + '\n')
    return out


def use_probe_recipes():
    af.RECIPES_DIR = RECIPES_DIR
    af.RECIPES = af._discover_recipes()
    af.build_phantom = _build_phantom
    af.phantom_digest = _phantom_digest


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


def generate(names, data_dir, threads=8, recipes_dir=None):
    use_probe_recipes()
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
