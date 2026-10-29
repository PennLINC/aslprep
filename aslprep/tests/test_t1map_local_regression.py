"""Tests for per-voxel tissue T1s from a T1 map, by local regression (gh-283)."""

import json
from pathlib import Path

import nibabel as nb
import numpy as np
import pytest

from aslprep import config
from aslprep.tests.tests import mock_config
from aslprep.utils import t1map as t1m

SHAPE = (24, 24, 24)
AFFINE = np.diag([1.0, 1.0, 1.0, 1.0])
R1 = {'gm': 0.75, 'wm': 1.25, 'csf': 0.25}
INTERIOR = (slice(3, -3),) * 3


def _anatomy(gm_gradient=0.0):
    """Gray matter, white matter and CSF that vary over a few voxels, as at 1 mm, with a voxel
    R1 mixed by their fractions.

    The gray matter R1 may vary along the first axis, by ``gm_gradient`` per voxel.
    """
    x, y, z = np.meshgrid(*(np.arange(n, dtype=float) for n in SHAPE), indexing='ij')
    gm = np.clip(0.5 + 0.7 * np.sin(x / 1.3 + y / 2.1), 0, 1)
    csf = 0.6 * np.clip(np.cos(z / 1.1 + x / 3), 0, 1) * (1 - gm)
    wm = 1 - gm - csf
    gm_r1 = R1['gm'] + gm_gradient * (x - SHAPE[0] / 2)
    r1 = gm * gm_r1 + wm * R1['wm'] + csf * R1['csf']
    return 1 / r1, gm, wm, csf, gm_r1


def _medians(gm_t1=1 / R1['gm'], wm_t1=1 / R1['wm']):
    return {
        'gray': {'t1': gm_t1, 'iqr': 0.0, 'n': 1000},
        'white': {'t1': wm_t1, 'iqr': 0.0, 'n': 1000},
    }


@pytest.mark.parametrize(
    ('kernel_mm', 'voxel_size', 'expected'),
    [(7.0, 1.0, 7), (7.0, 0.8, 9), (7.0, 1.5, 5), (7.0, 2.0, 3), (7.0, 4.0, 3), (3.0, 4.0, 3)],
)
def test_kernel_voxels(kernel_mm, voxel_size, expected):
    """The cube is the odd number of voxels nearest its edge in mm, and at least 3."""
    affine = np.diag([voxel_size, voxel_size, voxel_size, 1.0])
    assert t1m.kernel_voxels(kernel_mm, affine) == (expected,) * 3
    with pytest.raises(ValueError, match='finite and positive'):
        t1m.kernel_voxels(np.inf, affine)


def test_local_regression_r1_recovers_tissues():
    """With R1 mixed linearly, the regression separates each tissue's R1, CSF included."""
    t1, gm, wm, csf, _ = _anatomy()
    mask = np.ones(SHAPE, bool)
    result = t1m.local_regression_r1(t1, gm, wm, csf, mask, (7, 7, 7), _medians())

    status = result['status']
    assert np.isin(status[mask], list(t1m.R1_FIT_STATUS)).all()
    # Away from the image's edges, where the box is cut off
    assert (status[INTERIOR] == 1).mean() > 0.95
    for tissue in ('gm', 'wm'):
        estimated = result[f'{tissue}_estimated']
        assert estimated[INTERIOR].mean() > 0.95
        # Exact where all three tissues are fit, and within the contamination limit otherwise
        full = estimated & (status == 1)
        np.testing.assert_allclose(result[f'{tissue}_r1'][full], R1[tissue], rtol=1e-6)
        bias = t1m.MAX_CONTAMINATION * (max(R1.values()) - min(R1.values())) * 2
        np.testing.assert_allclose(result[f'{tissue}_r1'][estimated], R1[tissue], atol=bias)

    assert result['estimated_fraction']['gray'] > 0.7


def test_local_regression_r1_follows_regional_variation():
    """A gray matter R1 that varies across the image is followed, unlike its median."""
    gradient = 0.005
    t1, gm, wm, csf, gm_r1 = _anatomy(gm_gradient=gradient)
    mask = np.ones(SHAPE, bool)
    result = t1m.local_regression_r1(t1, gm, wm, csf, mask, (7, 7, 7), _medians())
    estimated = result['gm_estimated']
    assert estimated.mean() > 0.8
    # The kernel averages the gradient, so the error is a small fraction of the range
    error = np.abs(result['gm_r1'][estimated] - gm_r1[estimated])
    assert np.median(error) < 0.1 * gradient * SHAPE[0]
    assert np.corrcoef(result['gm_r1'][estimated], gm_r1[estimated])[0, 1] > 0.95


def test_local_regression_r1_fills_with_medians():
    """Outside the mask, where a tissue is absent, and for implausible estimates, the median."""
    t1, gm, wm, csf, _ = _anatomy()
    # No white matter in the first slab, so its R1 cannot be estimated there
    wm[:8] = 0.0
    gm[:8] = 1 - csf[:8]
    t1[:8] = 1 / (gm[:8] * R1['gm'] + csf[:8] * R1['csf'])
    mask = np.ones(SHAPE, bool)
    mask[-3:] = False
    medians = _medians(gm_t1=1.5, wm_t1=0.9)
    result = t1m.local_regression_r1(t1, gm, wm, csf, mask, (7, 7, 7), medians)

    assert not result['wm_estimated'][:5].any()
    np.testing.assert_allclose(result['wm_r1'][:5], 1 / 0.9)
    assert (result['status'][-3:] == 0).all()
    assert not result['gm_estimated'][-3:].any()
    np.testing.assert_allclose(result['gm_r1'][-3:], 1 / 1.5)

    # A median far from the local values makes every local estimate implausible
    result = t1m.local_regression_r1(t1, gm, wm, csf, mask, (7, 7, 7), _medians(gm_t1=10.0))
    assert not result['gm_estimated'].any()
    assert result['estimated_fraction']['gray'] == 0.0
    np.testing.assert_allclose(result['gm_r1'], 0.1)


def test_local_regression_r1_slabs_match_whole_volume(monkeypatch):
    """Unmixing in slabs, with halos, gives the same maps as unmixing the volume at once."""
    t1, gm, wm, csf, _ = _anatomy(gm_gradient=0.003)
    mask = np.zeros(SHAPE, bool)
    mask[2:-2, 3:-3, 1:-1] = True
    # Slabs of 5 slices put seams inside the mask's bounding box (which starts at slice 2)
    monkeypatch.setattr(t1m, '_SLAB', 5)
    sliced = t1m.local_regression_r1(t1, gm, wm, csf, mask, (5, 5, 5), _medians())
    monkeypatch.setattr(t1m, '_SLAB', SHAPE[0])
    whole = t1m.local_regression_r1(t1, gm, wm, csf, mask, (5, 5, 5), _medians())
    assert (whole['status'][mask] == 1).mean() > 0.5
    # The box filter's running sums round differently over different extents
    for key in ('gm_r1', 'wm_r1'):
        np.testing.assert_allclose(sliced[key], whole[key], rtol=1e-10)

    for key in ('status', 'gm_estimated', 'wm_estimated'):
        np.testing.assert_array_equal(sliced[key], whole[key])


def test_build_pv_maps_clip():
    """Maps of quantities other than probabilities are averaged without clipping."""
    from aslprep.utils.pvc import build_pv_maps

    shape = (8, 8, 8)
    field = np.full(shape, 1.5)
    maps, _ = build_pv_maps(
        [field, field],
        AFFINE,
        shape,
        AFFINE,
        np.eye(4),
        fwhm=(1.0, 1.0, 1.0),
        clip=[True, False],
    )
    np.testing.assert_allclose(maps[0][2:-2, 2:-2, 2:-2], 1.0)
    np.testing.assert_allclose(maps[1][2:-2, 2:-2, 2:-2], 1.5)
    with pytest.raises(ValueError, match='clip values'):
        build_pv_maps(
            [field],
            AFFINE,
            shape,
            AFFINE,
            np.eye(4),
            fwhm=(1, 1, 1),
            clip=[True, False],
        )


def _identity_xfm(path):
    path.write_text(
        '#Insight Transform File V1.0\n'
        '#Transform 0\n'
        'Transform: MatrixOffsetTransformBase_double_3_3\n'
        'Parameters: 1 0 0 0 1 0 0 0 1 0 0 0\n'
        'FixedParameters: 0 0 0\n'
    )
    return str(path)


def test_partial_volume_maps_tissue_t1(tmp_path):
    """Tissue T1s on the ASL grid are 1 over each tissue's R1, weighted as its PV map is."""
    from aslprep.interfaces.pvc import PartialVolumeMaps

    _, gm, wm, csf, gm_r1 = _anatomy(gm_gradient=0.004)
    files = {}
    for name, data in {'gm': gm, 'wm': wm, 'csf': csf, 'gm_r1': gm_r1}.items():
        files[name] = str(tmp_path / f'{name}.nii.gz')
        nb.Nifti1Image(data.astype(np.float32), AFFINE).to_filename(files[name])

    files['wm_r1'] = str(tmp_path / 'wm_r1.nii.gz')
    nb.Nifti1Image(np.full(SHAPE, R1['wm'], np.float32), AFFINE).to_filename(files['wm_r1'])
    # A coarser ASL grid
    ref_affine = np.diag([3.0, 3.0, 3.0, 1.0])
    reference = str(tmp_path / 'ref.nii.gz')
    nb.Nifti1Image(np.ones((8, 8, 8), np.uint8), ref_affine).to_filename(reference)
    out = (
        PartialVolumeMaps(
            gm_tpm=files['gm'],
            wm_tpm=files['wm'],
            csf_tpm=files['csf'],
            gm_r1=files['gm_r1'],
            wm_r1=files['wm_r1'],
            reference=reference,
            aslref2anat_xfm=_identity_xfm(tmp_path / 'xfm.txt'),
            metadata={'MRAcquisitionType': '3D'},
        )
        .run(cwd=tmp_path)
        .outputs
    )
    gm_pv = nb.load(out.gm_pv).get_fdata()
    gm_t1 = nb.load(out.gm_t1).get_fdata()
    wm_t1 = nb.load(out.wm_t1).get_fdata()
    inside = gm_pv > 0.05
    assert inside.sum() > 100
    # A constant R1 gives exactly its T1; a varying one gives a T1 within its range
    wm_inside = nb.load(out.wm_pv).get_fdata() > 0.05
    np.testing.assert_allclose(wm_t1[wm_inside], 1 / R1['wm'], rtol=1e-5)
    assert np.all(gm_t1[inside] >= 1 / gm_r1.max() - 1e-6)
    assert np.all(gm_t1[inside] <= 1 / gm_r1.min() + 1e-6)
    # Exactly the ratio of the responses to the PV-weighted R1 and to the PV
    from aslprep.utils.pvc import build_pv_maps

    (pv_ref, weighted_ref), _ = build_pv_maps(
        [gm.astype(np.float32), gm.astype(np.float32) * gm_r1.astype(np.float32)],
        AFFINE,
        (8, 8, 8),
        ref_affine,
        np.eye(4),
        fwhm=(3.0, 3.0, 3.0),
        clip=[True, False],
    )
    np.testing.assert_allclose(gm_t1[inside], (pv_ref / weighted_ref)[inside], rtol=1e-5)
    # And it follows the gradient along the first axis
    assert np.mean(gm_t1[:3][inside[:3]]) > np.mean(gm_t1[-3:][inside[-3:]])

    # Without R1 maps, there are no T1 outputs
    plain = PartialVolumeMaps(
        gm_tpm=files['gm'],
        wm_tpm=files['wm'],
        csf_tpm=files['csf'],
        reference=reference,
        aslref2anat_xfm=_identity_xfm(tmp_path / 'xfm.txt'),
        metadata={'MRAcquisitionType': '3D'},
    )
    plain_dir = tmp_path / 'plain'
    plain_dir.mkdir()
    plain_out = plain.run(cwd=plain_dir).outputs
    np.testing.assert_allclose(nb.load(plain_out.gm_pv).get_fdata(), gm_pv)
    assert not plain_out.gm_t1


def _calibration():
    from aslprep.utils.calibration import resolve_calibration

    metadata = {
        'ArterialSpinLabelingType': 'PCASL',
        'MagneticFieldStrength': 3,
        'MRAcquisitionType': '3D',
        'LabelingDuration': 1.8,
        'PostLabelingDelay': [0.5, 1.0, 1.5, 2.0, 2.5],
        'M0Type': 'Estimate',
        'M0Estimate': 1000.0,
        'BackgroundSuppression': False,
    }
    return metadata, resolve_calibration(metadata, kinetic_model='tissue-t1')


def test_voxelwise_tissue_t1prime():
    """Per-voxel T1s give per-voxel T1's, and missing ones fall back to the specification."""
    from aslprep.utils.calibration import (
        calibration_value,
        tissue_t1_kinetics,
        voxelwise_tissue_t1prime,
    )

    _, calibration = _calibration()
    scalar = tissue_t1_kinetics(calibration)['t1prime']
    defaults = [
        calibration_value(calibration, 'GrayMatterT1'),
        calibration_value(calibration, 'WhiteMatterT1'),
    ]
    t1s = np.array([defaults, [0.0, np.nan], [1.0, 0.7]])
    t1prime = voxelwise_tissue_t1prime(calibration, t1s)
    np.testing.assert_allclose(t1prime[0], scalar)
    np.testing.assert_allclose(t1prime[1], scalar)
    assert np.all(t1prime[2] < scalar)


def test_compute_cbf_tissue_t1_maps(tmp_path):
    """Maps of the specification's T1s change nothing; other T1s change CBF where they differ."""
    from aslprep.interfaces.cbf import ComputeCBF
    from aslprep.utils.calibration import tissue_t1_kinetics
    from aslprep.utils.pvc import gkm_predict

    metadata, calibration = _calibration()
    plds = np.repeat(metadata['PostLabelingDelay'], 2)
    metadata = {**metadata, 'PostLabelingDelay': plds.tolist()}
    shape = (4, 4, 2)
    affine = np.diag([3.0, 3.0, 5.0, 1.0])
    gm_pv = np.full(shape, 0.6)
    wm_pv = np.full(shape, 0.4)
    t1prime = tissue_t1_kinetics(calibration)['t1prime']
    deltam = gkm_predict(
        np.full(32, 50.0),
        np.full(32, 1.2),
        np.full(32, 0.5),
        np.zeros(32),
        np.full(32, 1000.0),
        np.tile(plds, (32, 1)),
        labeleff=0.85 * 1.0,
        t1blood=1.65,
        is_casl=True,
        tau=1.8,
        t1prime=t1prime[None, :1],
    ).reshape(*shape, -1)

    def write(name, data):
        path = str(tmp_path / f'{name}.nii.gz')
        nb.Nifti1Image(np.asarray(data, np.float32), affine).to_filename(path)
        return path

    files = {
        'deltam': write('deltam', deltam),
        'm0_file': write('m0', np.full(shape, 1000.0)),
        'gm_pv': write('gm_pv', gm_pv),
        'wm_pv': write('wm_pv', wm_pv),
    }
    gm_t1 = calibration['Parameters']['GrayMatterT1']['Value']
    wm_t1 = calibration['Parameters']['WhiteMatterT1']['Value']
    half = np.full(shape, gm_t1)
    half[:2] = 1.0

    def run(name, **t1_files):
        run_dir = tmp_path / name
        run_dir.mkdir()
        out = (
            ComputeCBF(
                cbf_only=False,
                metadata=metadata,
                calibration=calibration,
                **files,
                **t1_files,
            )
            .run(cwd=run_dir)
            .outputs
        )
        return nb.load(out.mean_cbf).get_fdata()

    scalar = run('scalar')
    same = run(
        'same',
        gm_t1=write('gm_same', np.full(shape, gm_t1)),
        wm_t1=write('wm_same', np.zeros(shape)),
    )
    np.testing.assert_allclose(same, scalar, rtol=1e-6)
    changed = run(
        'changed',
        gm_t1=write('gm_half', half),
        wm_t1=write('wm_t1', np.full(shape, wm_t1)),
    )
    np.testing.assert_allclose(changed[2:], scalar[2:], rtol=1e-6)
    assert not np.allclose(changed[:2], scalar[:2], rtol=1e-3)


def test_multidelay_lr_pvc_tissue_t1_maps(tmp_path):
    """The multi-delay PVC fits each tissue curve with that tissue's T1' in the voxel."""
    from aslprep.interfaces.pvc import MultiDelayLocalRegressionPVC
    from aslprep.utils.calibration import calibration_value

    metadata, calibration = _calibration()
    plds = np.repeat(metadata['PostLabelingDelay'], 2)
    metadata = {**metadata, 'PostLabelingDelay': plds.tolist()}
    rng = np.random.default_rng(3)
    shape = (12, 12, 2)
    affine = np.diag([3.0, 3.0, 5.0, 1.0])
    pv_gm = rng.uniform(0.2, 0.7, shape)
    pv_wm = 1 - pv_gm
    gm_t1 = calibration_value(calibration, 'GrayMatterT1')
    # A simple rising delta-M, shared by both tissues: only the T1s differ between the runs
    deltam = (pv_gm + 0.5 * pv_wm)[..., None] * np.linspace(5, 20, plds.size)

    def write(name, data, run_dir):
        path = str(run_dir / f'{name}.nii.gz')
        nb.Nifti1Image(np.asarray(data, np.float32), affine).to_filename(path)
        return path

    def run(name, gm_map=None):
        run_dir = tmp_path / name
        run_dir.mkdir()
        files = {
            key: write(key, data, run_dir)
            for key, data in {
                'deltam': deltam,
                'm0_calib': np.full(shape, 1000.0),
                'quant_valid': np.ones(shape),
                'mask': np.ones(shape),
                'gm_pv': pv_gm,
                'wm_pv': pv_wm,
                'csf_pv': np.zeros(shape),
                'mean_cbf': np.full(shape, 40.0),
                'att': np.full(shape, 1.2),
            }.items()
        }
        if gm_map is not None:
            files['gm_t1'] = write('gm_t1', gm_map, run_dir)
            files['wm_t1'] = write('wm_t1', np.zeros(shape), run_dir)

        out = (
            MultiDelayLocalRegressionPVC(
                metadata=metadata,
                calibration=calibration,
                kernel_size=[3, 3, 1],
                **files,
            )
            .run(cwd=run_dir)
            .outputs
        )
        return nb.load(out.gm_cbf).get_fdata(), nb.load(out.wm_cbf).get_fdata()

    gm_scalar, wm_scalar = run('scalar')
    assert np.count_nonzero(gm_scalar) > 50
    gm_same, wm_same = run('same', np.full(shape, gm_t1))
    np.testing.assert_allclose(gm_same, gm_scalar, rtol=1e-6)
    np.testing.assert_allclose(wm_same, wm_scalar, rtol=1e-6)

    # A shorter gray matter T1 in half the image changes only its gray matter CBF
    shorter = np.full(shape, gm_t1)
    shorter[:6] = 0.8 * gm_t1
    gm_changed, wm_changed = run('changed', shorter)
    np.testing.assert_allclose(wm_changed, wm_scalar, rtol=1e-6)
    np.testing.assert_allclose(gm_changed[6:], gm_scalar[6:], rtol=1e-6)
    fitted = gm_scalar[:6] > 0
    assert fitted.sum() > 20
    assert not np.allclose(gm_changed[:6][fitted], gm_scalar[:6][fitted], rtol=1e-3)


def test_fit_tissue_curves_per_curve_t1prime():
    """One T1' per curve is subset with the curves that fit, when some do not."""
    from aslprep.utils.pvc import fit_tissue_curves, gkm_predict

    plds = np.tile(np.repeat([0.5, 1.0, 1.5, 2.0, 2.5], 2), (4, 1))
    t1prime = np.array([[1.2], [1.3], [1.4], [1.5]])
    kinetics = {'labeleff': 0.85, 't1blood': 1.65, 'is_casl': True, 'tau': 1.8}
    curves = gkm_predict(
        np.full(4, 50.0),
        np.full(4, 1.0),
        np.full(4, 0.5),
        np.zeros(4),
        np.full(4, 1000.0),
        plds,
        t1prime=t1prime,
        **kinetics,
    )
    curves[1] = np.nan
    result = fit_tissue_curves(curves, np.full(4, 1000.0), plds, t1prime=t1prime, **kinetics)
    assert result['status'][1] == 4
    fitted = result['status'] != 4
    np.testing.assert_allclose(result['cbf'][fitted], 50.0, rtol=1e-3)
    np.testing.assert_allclose(result['kinres'][fitted], 0.0, atol=1e-6)


def test_tissue_t1map_sink_names(tmp_path):
    """The per-voxel tissue T1 maps are written as perf derivatives with a label entity."""
    from aslprep.interfaces.bids import DerivativesDataSink
    from aslprep.workflows.asl.base import _tissue_t1_metadata

    source = tmp_path / 'sub-01' / 'perf' / 'sub-01_acq-x_asl.nii.gz'
    source.parent.mkdir(parents=True)
    nb.Nifti1Image(np.zeros((2, 2, 2, 2), np.float32), AFFINE).to_filename(source)
    in_file = tmp_path / 'gm_t1.nii.gz'
    nb.Nifti1Image(np.ones((2, 2, 2), np.float32), AFFINE).to_filename(in_file)
    summary = {'Description': 'Estimated in each voxel.', 'Method': 'local-regression'}
    sink = DerivativesDataSink(
        base_directory=str(tmp_path / 'out'),
        source_file=str(source),
        in_file=str(in_file),
        label='GM',
        suffix='T1map',
        compress=True,
        Units='s',
        dismiss_entities=('echo',),
        meta_dict=_tissue_t1_metadata(summary, 'gm'),
    )
    out = sink.run(cwd=tmp_path).outputs.out_file
    assert out.endswith('sub-01/perf/sub-01_acq-x_label-GM_T1map.nii.gz')
    sidecar = json.loads(Path(out.replace('.nii.gz', '.json')).read_text())
    assert sidecar['Units'] == 's'
    assert sidecar['TissueT1Method']['Method'] == 'local-regression'
    assert 'gray matter' in sidecar['Description']


def test_unmix_tissue_r1_interface(tmp_path):
    """The interface writes the maps and completes the provenance of the medians."""
    from aslprep.interfaces.t1map import UnmixTissueR1
    from aslprep.utils.calibration import measured_tissue_t1_updates

    t1, gm, wm, csf, _ = _anatomy()
    files = {}
    for name, data in {'t1': t1, 'gm': gm, 'wm': wm, 'csf': csf, 'mask': np.ones(SHAPE)}.items():
        files[name] = str(tmp_path / f'{name}.nii.gz')
        nb.Nifti1Image(data.astype(np.float32), AFFINE).to_filename(files[name])

    medians = _medians()
    updates = measured_tissue_t1_updates(1 / R1['gm'], 1 / R1['wm'], {'gray': 'GM', 'white': 'WM'})
    out = (
        UnmixTissueR1(
            t1=files['t1'],
            tpms=[files['gm'], files['wm'], files['csf']],
            mask=files['mask'],
            kernel_mm=7.0,
            medians=medians,
            updates=updates,
        )
        .run(cwd=tmp_path)
        .outputs
    )
    status = nb.load(out.status).get_fdata()
    full = status == 1
    assert full.mean() > 0.8
    np.testing.assert_allclose(nb.load(out.wm_r1).get_fdata()[full], R1['wm'], rtol=1e-4)
    assert out.updates['GrayMatterT1']['Value'] == updates['GrayMatterT1']['Value']
    assert out.updates['GrayMatterT1']['Origin'].startswith('GM. ')
    assert 'local linear regression' in out.updates['WhiteMatterT1']['Origin']
    # The inputs' updates are not changed
    assert updates['GrayMatterT1']['Origin'] == 'GM'
    assert out.summary['KernelVoxels'] == [7, 7, 7]
    assert out.summary['EstimatedFraction']['GrayMatter'] > 0.7


def test_init_t1map_wf_local_regression():
    """The local-regression method adds the unmixing, whose updates the workflow outputs."""
    from aslprep.workflows.t1map import init_t1map_wf

    t1map = {
        't1map': 'sub-01_T1map.nii.gz',
        'metadata': {'Units': 's'},
        'uri': 'bids:raw:sub-01/anat/sub-01_T1map.nii.gz',
        'companion': 'sub-01_UNIT1.nii.gz',
    }
    with mock_config():
        median_wf = init_t1map_wf(t1map=t1map, identity=True, output_dir='.')
        config.workflow.t1map_method = 'local-regression'
        config.workflow.t1map_kernel_size = 5.0
        wf = init_t1map_wf(t1map=t1map, identity=True, output_dir='.')

    assert median_wf.get_node('unmix_tissue_r1') is None
    unmix = wf.get_node('unmix_tissue_r1')
    assert unmix.inputs.kernel_mm == 5.0
    edge = wf._graph.get_edge_data(unmix, wf.get_node('outputnode'))['connect']
    assert ('updates', 'updates') in edge
    assert ('gm_r1', 'gm_r1') in edge
    measure = wf.get_node('measure_tissue_t1')
    measure_edge = wf._graph.get_edge_data(measure, wf.get_node('outputnode'))
    assert ('updates', 'updates') not in measure_edge['connect']
    assert '5 mm cube' in wf.__desc__
    assert 'sub-01' not in wf.__desc__


def test_init_cbf_wf_tissue_t1_maps(tmp_path):
    """With the local-regression method, per-voxel T1s reach the fit and the multi-delay PVC."""
    import pandas as pd

    from aslprep.tests.tests import reset_config
    from aslprep.workflows.asl.cbf import init_cbf_wf

    asl_file = tmp_path / 'sub-01_asl.nii.gz'
    asl_file.touch()
    metadata, _ = _calibration()
    plds = metadata['PostLabelingDelay']
    metadata = {**metadata, 'PostLabelingDelay': [delay for delay in plds for _ in range(2)]}
    aslcontext = tmp_path / 'sub-01_aslcontext.tsv'
    pd.DataFrame({'volume_type': ['control', 'label'] * len(plds)}).to_csv(
        aslcontext,
        sep='\t',
        index=False,
    )
    t1map = {'t1map': 'x.nii.gz', 'metadata': {}, 'uri': 'x', 'companion': None}
    reset_config()
    with mock_config():
        config.workflow.kinetic_model = 'tissue-t1'
        config.workflow.t1map_method = 'local-regression'
        wf = init_cbf_wf(
            name_source=str(asl_file),
            processing_target='control',
            aslcontext=str(aslcontext),
            metadata=metadata,
            t1map=t1map,
            pvc_method='local-regression',
            pvc_kernel_size=(3, 3, 1),
        )
        desc = wf.__desc__

    reset_config()
    pv_wf = wf.get_node('asl_pv_wf')
    compute_cbf = wf.get_node('compute_cbf')
    lr_pvc_wf = wf.get_node('lr_pvc_wf')
    to_pv = wf._graph.get_edge_data(wf.get_node('inputnode'), pv_wf)['connect']
    assert ('tissue_gm_r1', 'inputnode.gm_r1') in to_pv
    to_cbf = wf._graph.get_edge_data(pv_wf, compute_cbf)['connect']
    assert ('outputnode.gm_t1', 'gm_t1') in to_cbf
    to_lr = wf._graph.get_edge_data(pv_wf, lr_pvc_wf)['connect']
    assert ('outputnode.wm_t1', 'inputnode.wm_t1') in to_lr
    assert 'by\nlocal regression' in desc or 'by local regression' in desc
    assert pv_wf.get_node('pv_maps') is not None
