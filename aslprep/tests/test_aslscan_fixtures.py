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


def test_rot_angle_is_accurate_near_identity():
    for deg in (1e-4, 0.01, 2.0, 90.0, 179.0):
        assert tg.rot_angle_deg(tg.rigid_ras(0, 0, 0, deg, 0, 0)) == pytest.approx(deg, rel=1e-6)
    assert tg.rot_angle_deg(np.eye(4)) == 0.0


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


# ---------------------------------------------------------------------------------------------
# Recipes and generation
# ---------------------------------------------------------------------------------------------
def _aslscan_or_skip():
    try:
        return af.find_aslscan()
    except af.AslscanUnavailable as exc:
        pytest.skip(str(exc))


def test_every_recipe_loads():
    for name in af.RECIPES:
        recipe = af.load_recipe(name)
        assert recipe.acq
        assert recipe.context


@pytest.mark.parametrize(
    ('extent', 'voxel', 'expected'),
    [(193.0, 5.0, 39), (190.0, 5.0, 38), (190.0000000001, 5.0, 38), (24.0, 4.0, 6), (1.0, 4.0, 1)],
)
def test_expected_slices_matches_aslscan_rule(extent, voxel, expected):
    assert af.expected_slices(extent, voxel) == expected


def _recipe_in(tmp_path, monkeypatch, sidecar, context=('control', 'label'), **settings):
    monkeypatch.setattr(af, 'RECIPES_DIR', tmp_path / 'recipes')
    d = tmp_path / 'recipes' / 'r'
    d.mkdir(parents=True)
    (d / 'asl.json').write_text(json.dumps(sidecar))
    (d / 'aslcontext.tsv').write_text('volume_type\n' + '\n'.join(context) + '\n')
    (d / 'overlay.toml').write_text(
        'seed = 1\n[acquisition]\nmatrix = [12, 10]\nnoise_variance = 0.0\n'
    )
    lines = [f'{k} = {json.dumps(v)}' for k, v in {'acq': 'test', **settings}.items()]
    (d / 'recipe.toml').write_text('\n'.join(lines) + '\n')
    return af.Recipe('r')


def test_check_recipe_slice_rules(tmp_path, monkeypatch):
    phantom = _write_synthetic_phantom(tmp_path / 'phantom')  # 8 mm in z
    side = {'MRAcquisitionType': '2D', 'AcquisitionVoxelSize': [2, 2, 3]}
    af.check_recipe(
        _recipe_in(tmp_path / 'a', monkeypatch, {**side, 'SliceTiming': [0, 0.1, 0.2]}), phantom
    )

    with pytest.raises(af.RecipeError, match='gives 3 slices'):
        af.check_recipe(
            _recipe_in(tmp_path / 'b', monkeypatch, {**side, 'SliceTiming': [0, 1]}), phantom
        )

    side4 = {**side, 'AcquisitionVoxelSize': [2, 2, 2], 'MultibandAccelerationFactor': 2}
    ok = _recipe_in(tmp_path / 'c', monkeypatch, {**side4, 'SliceTiming': [0, 0.1, 0, 0.1]})
    af.check_recipe(ok, phantom)
    bad = _recipe_in(tmp_path / 'd', monkeypatch, {**side4, 'SliceTiming': [0, 0, 0.1, 0.1]})
    with pytest.raises(af.RecipeError, match='must share a slice time'):
        af.check_recipe(bad, phantom)

    three_d = {'MRAcquisitionType': '3D', 'AcquisitionVoxelSize': [2, 2, 2], 'SliceTiming': [0]}
    with pytest.raises(af.RecipeError, match='3D acquisitions'):
        af.check_recipe(_recipe_in(tmp_path / 'e', monkeypatch, three_d), phantom)

    plds = {**side, 'SliceTiming': [0, 0.1, 0.2], 'PostLabelingDelay': [1, 2, 3]}
    with pytest.raises(af.RecipeError, match='PostLabelingDelay has 3 values'):
        af.check_recipe(_recipe_in(tmp_path / 'f', monkeypatch, plds), phantom)


def test_recipe_rejects_unknown_keys(tmp_path, monkeypatch):
    with pytest.raises(af.RecipeError, match=r'unknown recipe\.toml keys'):
        _recipe_in(tmp_path, monkeypatch, {}, colour='blue')


def test_set_noise_replaces_exactly_one_line():
    text = 'seed = 1\n[acquisition]\nnoise_variance = 0.0  # set by the generator\n'
    assert 'noise_variance = 2.5  #' in af._set_noise(text, 2.5)
    with pytest.raises(af.RecipeError):
        af._set_noise('seed = 1\n', 1.0)


def test_deltam_estimates_pairs_by_order():
    asl = np.array([10.0, 9.0, 11.0, 9.5, 3.0])  # label-first pairs, then a deltam row
    context = ['label', 'control', 'label', 'control', 'deltam']
    # control 9 - label 10, control 9.5 - label 11, then the deltam row
    np.testing.assert_allclose(af.deltam_estimates(asl, context), [-1.0, -1.5, 3.0])


def test_fixture_dir_refuses_generation_when_required(tmp_path, monkeypatch):
    monkeypatch.setenv('ASLPREP_REQUIRE_FIXTURES', '1')
    with pytest.raises(af.FixtureUnavailable, match='forbids generating'):
        af.fixture_dir('fast_pcasl_seq', tmp_path)


def test_generated_fixture_matches_aslscan_truth(tmp_path):
    """Our exact-overlap partial volumes reproduce aslscan's perfusion truth (spec 4.4)."""
    binary = _aslscan_or_skip()
    from aslprep.tests import truth_geometry as tg

    out = af.generate('fast_pcasl_seq', tmp_path, binary)
    assert af.verify(tmp_path, ['fast_pcasl_seq']) == {}
    phantom = af.build_phantom(af.load_recipe('fast_pcasl_seq').phantom, tmp_path)
    perf = nb.load(phantom / 'perfusion.nii.gz')
    truth = nb.load(out / 'sub-01/perf/ground-truth/sub-01_desc-perfusion_gt.nii.gz')
    mine = tg.overlap_mean(perf.get_fdata(), perf.affine, truth.affine, truth.shape)
    np.testing.assert_allclose(mine, truth.get_fdata(), atol=1e-4)

    # tampering is detected
    (out / 'README').write_text('edited')
    assert af.verify(tmp_path, ['fast_pcasl_seq']) == {
        'fast_pcasl_seq': 'files differ from their manifest'
    }


def test_thread_count_does_not_change_output():
    binary = _aslscan_or_skip()
    assert 'identical' in af.check_threads('fast_pcasl_seq', binary)


# ---------------------------------------------------------------------------------------------
# Scoring (no ASLPrep run: outputs are fabricated from a fast fixture)
# ---------------------------------------------------------------------------------------------
def _write_itk_array(path, matrices):
    """An ITK multi-transform text file, as niworkflows' MCFLIRT2ITK writes for HMC."""
    blocks = ['#Insight Transform File V1.0']
    for i, m in enumerate(matrices):
        lps = tg.ras_to_lps(m)
        params = ' '.join(repr(float(v)) for v in [*lps[:3, :3].ravel(), *lps[:3, 3]])
        blocks += [
            f'#Transform {i}',
            'Transform: AffineTransform_double_3_3',
            f'Parameters: {params}',
            'FixedParameters: 0 0 0',
        ]
    path.write_text('\n'.join(blocks) + '\n')


def _fake_run(fixture, out, scale=1.0, shrink=False, nans=False):
    """Fabricate ASLPrep's native outputs: CBF = scale x the Tier A expectation."""
    from aslprep.tests import truth_scoring as ts

    fx = ts.Fixture(fixture, af.load_recipe('fast_pcasl_seq').acq)
    perf = out / 'sub-01' / 'perf'
    perf.mkdir(parents=True)
    cbf = np.nan_to_num(ts.expected_native(fx, fwhm=0.0)) * scale
    mask = fx.brain.copy()
    if shrink:
        mask[: mask.shape[0] // 2] = False
    if nans:
        cbf[fx.brain] = np.nan
    stem = f'sub-01_acq-{fx.acq}'
    nb.Nifti1Image(cbf.astype(np.float32), fx.affine).to_filename(perf / f'{stem}_cbf.nii.gz')
    nb.Nifti1Image(mask.astype(np.uint8), fx.affine).to_filename(
        perf / f'{stem}_desc-brain_mask.nii.gz'
    )
    return fx


@pytest.fixture
def fast_fixture(data_dir):
    try:
        return af.fixture_dir('fast_pcasl_seq', data_dir)
    except af.FixtureUnavailable as exc:
        pytest.skip(str(exc))


def test_scoring_exact_output_passes_and_errors_fail(fast_fixture, tmp_path):
    from aslprep.tests import truth_bounds as tb
    from aslprep.tests import truth_scoring as ts

    _fake_run(fast_fixture, tmp_path / 'exact')
    score = {
        'native': ts.score_native(
            ts.Fixture(fast_fixture, 'fastpcaslseq'), tmp_path / 'exact', 0.0
        )
    }
    assert score['native']['grid_matches_input']
    assert score['native']['tier_a']['median_abs_dev'] < 1e-6
    tb.check(score, 'tier_a_median', 'fast_pcasl_seq')
    tb.check(score, 'mask_coverage', 'fast_pcasl_seq')
    tb.check(
        score, 'tier_b', 'fast_pcasl_seq', reference=('native', 'expected_ratio', '{t}'), t='GM'
    )

    for kwargs, name, message in (
        ({'scale': 1.1}, 'tier_a_median', 'expected <= 0.05'),
        ({'shrink': True}, 'mask_coverage', 'expected >= 0.95'),
        ({'nans': True}, 'finite', 'expected >= 0.999'),
    ):
        out = tmp_path / name
        fx = _fake_run(fast_fixture, out, **kwargs)
        bad = {'native': ts.score_native(fx, out, 0.0)}
        with pytest.raises(AssertionError, match=message):
            tb.check(bad, name, 'fast_pcasl_seq')


def test_bounds_fail_on_missing_metric_or_ceiling():
    from aslprep.tests import truth_bounds as tb

    with pytest.raises(AssertionError, match='missing or not finite'):
        tb.check({'native': {}}, 'tier_a_median', 'any')
    with pytest.raises(AssertionError, match='no ceiling named'):
        tb.check({}, 'not_a_metric', 'any')
    with pytest.raises(AssertionError, match='missing or not finite'):
        tb.check({'native': {'tier_a': {'median_abs_dev': float('nan')}}}, 'tier_a_median', 'any')


def test_regression_bands_are_narrow():
    """A band wider than 0.1 in ratio units could hide a 10 % scaling error."""
    from aslprep.tests import truth_bounds as tb

    for (recipe, name), (lo, hi, _) in tb.BANDS.items():
        if tb.CEILINGS[name][3] == 'ratio':
            assert hi - lo <= 0.1, (recipe, name)


def test_frames_algebra(tmp_path):
    """Spec 9.3: consistent HMC and coregistration give no error; a 2 degree error is measured."""
    from types import SimpleNamespace

    from aslprep.tests import truth_scoring as ts

    r = tg.rigid_ras(3, -4, 2, 5, -3, 4)
    poses = np.stack(
        [tg.pose_matrix([0.5 * v, 0, 0], [0, 0, 0.01 * v], [0, 0, 0]) for v in range(4)]
    )
    affine = np.diag([3.0, 3.0, 3.0, 1.0])
    affine[:3, 3] = -30
    brain = np.zeros((20, 20, 20), bool)
    brain[4:16, 4:16, 4:16] = True
    fx = SimpleNamespace(
        acq='test',
        R=r,
        brain=brain,
        affine=affine,
        motion=[{}] * len(poses),  # a moving recipe: the aslref pose comes from the HMC transforms
        poses=lambda: poses,
        voxel_centers=lambda m: (np.argwhere(m), tg.apply_points(affine, np.argwhere(m))),
    )
    perf = tmp_path / 'sub-01' / 'perf'
    perf.mkdir(parents=True)
    # aslref in the static frame: H_v = P_v. Coregistration pulls T1w points: C = R^-1.
    _write_itk_array(
        perf / 'sub-01_acq-test_from-orig_to-aslref_mode-image_desc-hmc_xfm.txt', poses
    )
    coreg = perf / 'sub-01_acq-test_from-aslref_to-T1w_mode-image_desc-coreg_xfm.txt'
    tg.write_itk_affine(coreg, np.linalg.inv(r))
    good = ts.score_frames(fx, tmp_path)
    assert good['coreg']['rot_deg'] < 1e-3
    assert good['coreg']['rms_mm'] < 1e-3
    assert good['motion']['rms_error_max_mm'] < 1e-3

    tg.write_itk_affine(coreg, np.linalg.inv(r) @ tg.rigid_ras(0, 0, 0, 0, 0, 2))
    bad = ts.score_frames(fx, tmp_path)
    assert bad['coreg']['rot_deg'] == pytest.approx(2.0, abs=1e-6)
    assert bad['coreg']['rms_mm'] > 0.1


def test_derivative_transforms_undo_the_offset(tmp_path):
    """Spec 9.2: resampling the offset T1w through from-T1w_to-MNI recovers the phantom T1w."""
    from nitransforms.linear import load
    from nitransforms.resampling import apply

    phantom = _write_synthetic_phantom(tmp_path / 'phantom', shape=(30, 32, 28))
    for name in ('brainmask', 'probseg-GM', 'probseg-WM', 'probseg-CSF'):
        nb.load(phantom / 'anat' / 'T1w.nii.gz').to_filename(phantom / 'anat' / f'{name}.nii.gz')
    t1w = nb.load(phantom / 'anat' / 'T1w.nii.gz')
    smooth = nb.Nifti1Image(
        np.asarray(
            __import__('scipy.ndimage', fromlist=['gaussian_filter']).gaussian_filter(
                t1w.get_fdata(), 1.0
            ),
            np.float32,
        ),
        t1w.affine,
    )
    smooth.to_filename(phantom / 'anat' / 'T1w.nii.gz')

    r = tg.rigid_ras(3, -4, 2, 5, -3, 4)
    out = tmp_path / 'fixture'
    af._write_anatomy(phantom, out, r, derivatives=True)
    danat = out / 'derivatives' / 'anat' / 'sub-01' / 'anat'
    moved = nb.load(danat / 'sub-01_desc-preproc_T1w.nii.gz')
    np.testing.assert_allclose(moved.affine, r @ smooth.affine, atol=1e-6)

    xfm = load(danat / f'sub-01_from-T1w_to-{af.TEMPLATE}_mode-image_xfm.txt', fmt='itk')
    back = np.asarray(apply(xfm, moved, reference=smooth, order=1).dataobj)
    interior = (slice(5, -5),) * 3
    a, b = back[interior].ravel(), smooth.get_fdata()[interior].ravel()
    assert np.corrcoef(a, b)[0, 1] > 0.99


def test_simulator_motion_convention(data_dir):
    """Spec 9.4: each moved volume matches the static one pushed through ``pose_matrix``.

    The prediction is compared with alternatives (the inverse pose; rotation about the world
    origin), which must fit worse, so a convention error cannot pass.
    """
    from scipy.ndimage import map_coordinates

    from aslprep.tests import truth_scoring as ts

    try:
        fixture = af.fixture_dir('geom_motion', data_dir)
    except af.FixtureUnavailable as exc:
        pytest.skip(str(exc))
    fx = ts.Fixture(fixture, 'geommotion')
    data = fx.asl_img.get_fdata()
    static = data[..., 0]  # volume 0 (control) has the identity pose
    ijk = np.argwhere(np.ones(fx.shape, bool))
    world = tg.apply_points(fx.affine, ijk)
    interior = np.zeros(fx.shape, bool)
    interior[3:-3, 3:-3, 1:-1] = True
    center = fx.truth['motion_center']

    def predicted(pose):
        # the moved head shows, at x, what the static head had at pose^-1 x
        src = tg.apply_points(np.linalg.inv(pose), world)
        vox = tg.apply_points(np.linalg.inv(fx.affine), src)
        return map_coordinates(static, vox.T, order=1, cval=np.nan).reshape(fx.shape)

    def fit(volume, pose):
        pred = predicted(pose)
        ok = interior & np.isfinite(pred)
        return np.corrcoef(pred[ok], volume[ok])[0, 1]

    for v, row in enumerate(fx.motion):
        if v == 0 or fx.context[v] != 'control':
            continue
        trans = [row['trans_x'], row['trans_y'], row['trans_z']]
        rot = [row['rot_x'], row['rot_y'], row['rot_z']]
        pose = tg.pose_matrix(trans, rot, center)
        r_true = fit(data[..., v], pose)
        # interpolating an ASL-resolution image limits r (0.89 measured for a pure shift);
        # the discriminating checks are the comparisons with the wrong conventions below
        assert r_true > 0.8, f'volume {v}: r = {r_true:.3f}'
        assert r_true > fit(data[..., v], np.linalg.inv(pose))
        if any(rot):  # the centre only matters for rotations
            assert r_true > fit(data[..., v], tg.pose_matrix(trans, rot, [0, 0, 0]))


def test_fixture_digest_is_per_recipe(tmp_path, monkeypatch):
    """Editing one recipe invalidates only its own fixture (spec Section 4.5)."""
    side = {'MRAcquisitionType': '2D', 'AcquisitionVoxelSize': [2, 2, 3]}
    _recipe_in(tmp_path, monkeypatch, side)  # creates recipes/r and points RECIPES_DIR there
    other = tmp_path / 'recipes' / 'other'
    other.mkdir()
    for f in (tmp_path / 'recipes' / 'r').iterdir():
        (other / f.name).write_bytes(f.read_bytes())
    monkeypatch.setattr(af, 'RECIPES', {'r': '', 'other': ''})
    monkeypatch.setattr(af, 'phantom_digest', lambda name: 'fixed')
    before = {name: af.fixture_digest(name) for name in ('r', 'other')}
    (other / 'aslcontext.tsv').write_text('volume_type\nlabel\ncontrol\n')
    assert af.fixture_digest('r') == before['r']
    assert af.fixture_digest('other') != before['other']


def test_frames_without_motion_separate_hmc_from_coregistration(tmp_path):
    """Without motion, a spurious HMC transform is a motion error, not a coregistration error."""
    from types import SimpleNamespace

    from aslprep.tests import truth_scoring as ts

    r = tg.rigid_ras(3, -4, 2, 5, -3, 4)
    affine = np.diag([3.0, 3.0, 3.0, 1.0])
    affine[:3, 3] = -30
    brain = np.zeros((20, 20, 20), bool)
    brain[4:16, 4:16, 4:16] = True
    fx = SimpleNamespace(
        acq='test',
        R=r,
        brain=brain,
        affine=affine,
        motion=None,
        poses=lambda: np.repeat(np.eye(4)[None], 3, axis=0),
        voxel_centers=lambda m: (np.argwhere(m), tg.apply_points(affine, np.argwhere(m))),
    )
    perf = tmp_path / 'sub-01' / 'perf'
    perf.mkdir(parents=True)
    spurious = [np.eye(4), np.eye(4), tg.rigid_ras(1.0, 0, 0, 0, 0, 2)]
    _write_itk_array(
        perf / 'sub-01_acq-test_from-orig_to-aslref_mode-image_desc-hmc_xfm.txt', spurious
    )
    tg.write_itk_affine(
        perf / 'sub-01_acq-test_from-aslref_to-T1w_mode-image_desc-coreg_xfm.txt', np.linalg.inv(r)
    )
    frames = ts.score_frames(fx, tmp_path)
    assert frames['coreg']['rms_mm'] < 1e-3  # coregistration is exact
    assert frames['motion']['rms_error_median_mm'] < 1e-3
    assert frames['motion']['rms_error_max_mm'] > 1.0  # the third volume moved spuriously
    assert frames['motion']['reference_offset_mm'] < 1e-3  # most volumes stayed in place
