"""Tests for the aslscan fixture registry, geometry helpers and phantoms.

None of these need a generated fixture or the aslscan binary unless marked otherwise.
"""

import json

import nibabel as nb
import numpy as np
import pytest

from aslprep.tests import aslscan_fixtures as af
from aslprep.tests import truth_geometry as tg


def test_spec_file_matches():
    """The committed spec file (the CI cache key) must match the registry."""
    committed = af.SPEC_FILE.read_text()
    assert committed == af.spec_text(), (
        f'{af.SPEC_FILE.name} is out of date with the fixture registry. Regenerate it with:\n'
        f'    {af.SPEC_COMMAND}'
    )


def test_hashed_modules_exist():
    for mod in af.HASHED_MODULES:
        assert (af.TESTS_DIR / mod).is_file(), mod


def test_digest_dir_line_endings(tmp_path):
    """Line endings do not change a digest; names and lengths do."""
    lf, crlf = tmp_path / 'lf', tmp_path / 'crlf'
    for d, sep in ((lf, b'\n'), (crlf, b'\r\n')):
        (d / 'sub').mkdir(parents=True)
        (d / 'a.json').write_bytes(b'{' + sep + b'"x": 1' + sep + b'}' + sep)
        (d / 'sub' / 'b.tsv').write_bytes(b'volume_type' + sep + b'control' + sep)
    assert af._digest_dir(lf) == af._digest_dir(crlf)

    before = af._digest_dir(lf)
    (lf / 'sub' / 'b.tsv').rename(lf / 'sub' / 'c.tsv')
    assert af._digest_dir(lf) != before

    before = af._digest_dir(lf)
    (lf / 'a.json').write_bytes(b'{\n"x": 10\n}\n')
    assert af._digest_dir(lf) != before


def test_check_aslscan_requires_matching_stamp(tmp_path):
    binary = tmp_path / 'aslscan'
    binary.write_bytes(b'not really a binary')
    with pytest.raises(af.AslscanUnavailable, match=r'no aslscan\.build\.json'):
        af.check_aslscan(binary)

    stamp = {**af._expected_stamp(), 'sha256': af._sha256_file(binary)}
    (tmp_path / af.STAMP_NAME).write_text(json.dumps(stamp))
    assert af.check_aslscan(binary) == binary

    binary.write_bytes(b'rebuilt from something else')
    with pytest.raises(af.AslscanUnavailable, match='does not match the hash'):
        af.check_aslscan(binary)

    stamp = {**af._expected_stamp(), 'aslscan': '0' * 40, 'sha256': af._sha256_file(binary)}
    (tmp_path / af.STAMP_NAME).write_text(json.dumps(stamp))
    with pytest.raises(af.AslscanUnavailable, match='was built from'):
        af.check_aslscan(binary)


def test_find_aslscan_reports_every_candidate(tmp_path, monkeypatch):
    monkeypatch.setenv('ASLSCAN', str(tmp_path / 'missing'))
    monkeypatch.setenv('PATH', str(tmp_path))
    with pytest.raises(af.AslscanUnavailable, match='missing does not exist'):
        af.find_aslscan()


def test_spec_digest_tracks_inputs(monkeypatch):
    before = af.spec_digest()
    monkeypatch.setattr(af, 'CACHE_EPOCH', af.CACHE_EPOCH + 1)
    assert af.spec_digest() != before


# ---------------------------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------------------------
def test_rigid_ras_single_axes():
    def rot(rx, ry, rz, v):
        return tg.rigid_ras(0, 0, 0, rx, ry, rz)[:3, :3] @ v

    np.testing.assert_allclose(rot(0, 0, 90, [1, 0, 0]), [0, 1, 0], atol=1e-12)
    np.testing.assert_allclose(rot(90, 0, 0, [0, 1, 0]), [0, 0, 1], atol=1e-12)
    np.testing.assert_allclose(rot(0, 90, 0, [0, 0, 1]), [1, 0, 0], atol=1e-12)
    np.testing.assert_allclose(tg.rigid_ras(1, 2, 3, 0, 0, 0)[:3, 3], [1, 2, 3])
    # Rz Ry Rx: x is applied first, so x stays put under Rx and Rz takes it to y
    # (Rx Rz would give z instead)
    np.testing.assert_allclose(rot(90, 0, 90, [1, 0, 0]), [0, 1, 0], atol=1e-12)


def test_pose_matrix_rotates_about_center():
    center = [10.0, -20.0, 5.0]
    m = tg.pose_matrix([1, 2, 3], [0.1, -0.2, 0.3], center)
    np.testing.assert_allclose(tg.apply_points(m, [center])[0], [11, -18, 8])
    np.testing.assert_allclose(m[:3, :3], tg.rotation_zyx(0.1, -0.2, 0.3))
    assert tg.rot_angle_deg(m) > 0


def test_fov_center():
    affine = np.diag([2.0, 3.0, 4.0, 1.0])
    affine[:3, 3] = [-10, -20, -30]
    np.testing.assert_allclose(tg.fov_center(affine, (11, 21, 31)), [0, 10, 30])


def test_itk_affine_round_trip_through_nitransforms(tmp_path):
    """Our writer and nitransforms' reader agree on the RAS matrix (independent LPS handling)."""
    from nitransforms.linear import load

    m = tg.rigid_ras(3, -4, 2, 5, -3, 4)
    tg.write_itk_affine(tmp_path / 'xfm.txt', m)
    np.testing.assert_allclose(load(tmp_path / 'xfm.txt', fmt='itk').matrix, m, atol=1e-9)


def test_overlap_matrix_hand_cases():
    # aligned 1 mm cells into 2 mm cells
    w = tg.overlap_matrix(np.arange(5.0), np.array([0.0, 2.0, 4.0]))
    np.testing.assert_allclose(w, [[0.5, 0.5, 0, 0], [0, 0, 0.5, 0.5]])
    # half-cell offset
    w = tg.overlap_matrix(np.arange(4.0), np.array([0.5, 2.5]))
    np.testing.assert_allclose(w, [[0.25, 0.5, 0.25]])
    # target extends past the source: weights sum to the covered fraction
    w = tg.overlap_matrix(np.array([0.0, 1.0, 2.0]), np.array([1.0, 4.0]))
    np.testing.assert_allclose(w.sum(), 1 / 3)
    # fractional 3.5 mm cells over 1 mm cells
    w = tg.overlap_matrix(np.arange(8.0), np.array([0.0, 3.5, 7.0]))
    np.testing.assert_allclose(w.sum(axis=1), [1, 1])
    np.testing.assert_allclose(w[0, 3], 0.5 / 3.5)


def test_overlap_fractions_sum_to_coverage():
    rng = np.random.default_rng(0)
    labels = rng.integers(0, 4, size=(20, 22, 18))
    src = np.eye(4)
    dst = np.diag([3.5, 3.5, 5.0, 1.0])
    dst[:3, 3] = [1.25, 1.25, 2.0]  # cell edges at -0.5, the source's first edge
    frac = tg.overlap_fractions(labels, src, dst, (5, 6, 3), label_values=(0, 1, 2, 3))
    total = sum(frac.values())
    np.testing.assert_allclose(total[:5, :6, :3], 1.0, atol=1e-12)


def test_overlap_rejects_oblique_grids():
    oblique = tg.rigid_ras(0, 0, 0, 0, 0, 10)
    with pytest.raises(ValueError, match='axis-aligned'):
        tg.overlap_mean(np.zeros((3, 3, 3)), oblique, np.eye(4), (3, 3, 3))


def _blob_image(points, affine, shape, sigma=1.5):
    ijk = np.stack(np.meshgrid(*[np.arange(n) for n in shape], indexing='ij'), axis=-1)
    world = ijk @ affine[:3, :3].T + affine[:3, 3]
    data = np.zeros(shape)
    for p in points:
        data += np.exp(-np.sum((world - p) ** 2, axis=-1) / (2 * sigma**2))
    return nb.Nifti1Image(data.astype(np.float32), affine), world


def test_itk_pull_convention_by_landmarks(tmp_path):
    """Spec 9.1: resampling through a written transform moves landmarks to M^-1 p.

    An ITK image transform maps reference points to moving points (a pull mapping), so a blob
    at world point ``p`` in the moving image appears at ``M^-1 p`` on the reference grid.
    """
    from nitransforms.linear import load
    from nitransforms.resampling import apply

    affine = np.diag([2.0, 2.0, 2.0, 1.0])
    affine[:3, 3] = -39
    shape = (40, 40, 40)
    points = np.array([[x, y, z] for x in (-12, 12) for y in (-14, 10) for z in (-8, 13)], float)
    moving, world = _blob_image(points, affine, shape)

    m = tg.rigid_ras(3, -4, 2, 5, -3, 4)
    tg.write_itk_affine(tmp_path / 'xfm.txt', m)
    moved = apply(load(tmp_path / 'xfm.txt', fmt='itk'), moving, reference=moving, order=1)
    data = np.asarray(moved.dataobj)

    expected = tg.apply_points(np.linalg.inv(m), points)
    for p in expected:
        near = np.sum((world - p) ** 2, axis=-1) < 5**2
        w = data * near
        centroid = (world * w[..., None]).sum(axis=(0, 1, 2)) / w.sum()
        np.testing.assert_allclose(centroid, p, atol=0.1)


# ---------------------------------------------------------------------------------------------
# Phantoms
# ---------------------------------------------------------------------------------------------
def _write_synthetic_phantom(path, shape=(12, 10, 8)):
    """A small phantom satisfying aslscan's contract: three label slabs and background."""
    path.mkdir(parents=True)
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    affine[:3, 3] = [-6, -5, -4]
    labels = np.zeros(shape, np.int16)
    labels[2:10, 2:8, 1:3] = 1
    labels[2:10, 2:8, 3:5] = 2
    labels[2:10, 2:8, 5:7] = 3
    for quantity, units in af.PHANTOM_MAPS.items():
        data = np.zeros(shape, np.float32)
        for value, tissue in af.TISSUES.items():
            data[labels == value] = tissue[quantity]
        nb.Nifti1Image(data, affine).to_filename(path / f'{quantity}.nii.gz')
        (path / f'{quantity}.json').write_text(json.dumps({'Units': units}))
    nb.Nifti1Image(labels, affine).to_filename(path / 'dseg.nii.gz')
    (path / 'dseg.json').write_text(
        json.dumps({'Units': 'label indices', 'LabelMap': {'1': 'grey_matter'}})
    )
    (path / 'anat').mkdir()
    nb.Nifti1Image(labels.astype(np.float32), affine).to_filename(path / 'anat' / 'T1w.nii.gz')
    return path


def test_phantom_contract_accepts_valid_and_rejects_violations(tmp_path):
    good = _write_synthetic_phantom(tmp_path / 'good')
    af.check_phantom_contract(good)

    def corrupt(name, edit):
        img = nb.load(good / f'{name}.nii.gz')
        data = np.asanyarray(img.dataobj).copy()
        edit(data)
        nb.Nifti1Image(data, img.affine).to_filename(good / f'{name}.nii.gz')

    corrupt('M0map', lambda d: d.__setitem__((0, 0, 0), 5.0))
    with pytest.raises(af.PhantomContractError, match='zero in the background'):
        af.check_phantom_contract(good)
    corrupt('M0map', lambda d: d.__setitem__((0, 0, 0), 0.0))

    corrupt('T2starmap', lambda d: d.__setitem__((5, 5, 2), 1.0))
    with pytest.raises(af.PhantomContractError, match='below T2map'):
        af.check_phantom_contract(good)
    corrupt('T2starmap', lambda d: d.__setitem__((5, 5, 2), 0.066))

    (good / 'att.json').write_text(json.dumps({'Units': 'ms'}))
    with pytest.raises(af.PhantomContractError, match='Units'):
        af.check_phantom_contract(good)


def test_sanitize_maps_records_every_change():
    dseg = np.array([0, 1, 1, 2, 2])
    maps = {
        'perfusion': np.array([0, 50, -1, 20, np.nan]),
        'att': np.array([0, 0.8, 0.9, 1.2, 1.3]),
        'T1map': np.array([0, 1.3, 0.0, 0.8, 0.8]),
        'T2map': np.array([0, 0.08, 0.08, 0.1, 0.1]),
        'T2starmap': np.array([0, 0.09, 0.05, 0.05, 0.05]),
        'M0map': np.array([7, 70, 70, 60, 60]),
    }
    out, changes = af.sanitize_maps(maps, dseg)
    assert out['M0map'][0] == 0
    assert out['T1map'][2] == pytest.approx(0.05)
    assert out['T2starmap'][1] == pytest.approx(0.95 * 0.08)
    np.testing.assert_array_equal(out['perfusion'][[2, 4]], 0)
    np.testing.assert_array_equal(out['att'][[2, 4]], 1000)
    assert set(changes) == {
        'M0map zeroed outside the segmentation',
        'T1map raised to 0.05 in foreground',
        'T2starmap clipped to 0.95 x T2map',
        'perfusion set to 0 where not positive',
        'att set to 1000 where unperfused',
    }


def test_crop_keeps_world_coordinates(tmp_path):
    base = _write_synthetic_phantom(tmp_path / 'base')
    out = tmp_path / 'crop'
    out.mkdir()
    af._crop_phantom(base, '2:10,1:9,1:7', out)
    af.check_phantom_contract(out)
    full, crop = nb.load(base / 'perfusion.nii.gz'), nb.load(out / 'perfusion.nii.gz')
    assert crop.shape == (8, 8, 6)
    # voxel (0, 0, 0) of the crop is voxel (2, 1, 1) of the base, at the same world point
    np.testing.assert_allclose(crop.affine @ [0, 0, 0, 1], full.affine @ [2, 1, 1, 1])
    np.testing.assert_array_equal(crop.get_fdata(), full.get_fdata()[2:10, 1:9, 1:7])
    assert (out / 'anat' / 'T1w.nii.gz').exists()


def test_modulation_field_range_and_variation():
    affine = np.diag([2.0, 2.0, 2.0, 1.0])
    affine[:3, 3] = -60
    params = af.PHANTOM_PARAMS['tfmni']['modulation']['perfusion']
    field = af.modulation_field(params, affine, (60, 60, 60))
    assert field.min() >= 1 - params['amplitude'] - 1e-9
    assert field.max() <= 1 + params['amplitude'] + 1e-9
    assert field.std() > 0.02
