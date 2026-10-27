"""Tests for aslprep.interfaces.pvc."""

import nibabel as nb
import numpy as np
import pytest
from nitransforms.linear import load as load_xfm
from nitransforms.resampling import apply as apply_xfm

from aslprep.interfaces import pvc
from aslprep.tests import truth_geometry as tg


def _sphere(shape, centre_vox, radius_vox, width=2.0):
    grid = np.indices(shape).astype(float)
    dist = np.sqrt(sum((g - c) ** 2 for g, c in zip(grid, centre_vox, strict=True)))
    return np.clip((radius_vox + width / 2 - dist) / width, 0, 1).astype(np.float32)


def _centre_of_mass_ras(data, affine):
    idx = np.indices(data.shape).reshape(3, -1)
    weights = data.ravel()
    com = (idx * weights).sum(axis=1) / weights.sum()
    return (affine @ np.append(com, 1))[:3]


@pytest.fixture
def pv_inputs(tmp_path):
    """An anatomical sphere, a 4 mm ASL grid, and an aslref-to-anat transform with a shift."""
    anat_affine = np.eye(4)
    anat_affine[:3, 3] = [-30, -30, -30]
    gm = _sphere((60, 60, 60), (30, 30, 30), 10)
    wm = _sphere((60, 60, 60), (30, 30, 30), 5)
    gm = np.clip(gm - wm, 0, 1)
    csf = np.clip(_sphere((60, 60, 60), (30, 30, 30), 14) - gm - wm, 0, 1)
    gm_file = tmp_path / 'gm.nii.gz'
    wm_file = tmp_path / 'wm.nii.gz'
    csf_file = tmp_path / 'csf.nii.gz'
    nb.Nifti1Image(gm, anat_affine).to_filename(gm_file)
    nb.Nifti1Image(wm, anat_affine).to_filename(wm_file)
    nb.Nifti1Image(csf, anat_affine).to_filename(csf_file)

    ref_affine = np.diag([4.0, 4.0, 4.0, 1.0])
    ref_affine[:3, 3] = [-28, -28, -28]
    ref_file = tmp_path / 'ref.nii.gz'
    nb.Nifti1Image(np.ones((15, 15, 15), np.uint8), ref_affine).to_filename(ref_file)

    # The file maps anatomical points to aslref points (it pulls the T1w onto the ASL grid).
    shift = np.eye(4)
    shift[:3, 3] = [5.0, -3.0, 2.0]
    xfm_file = tmp_path / 'from-aslref_to-T1w_xfm.txt'
    tg.write_itk_affine(xfm_file, shift)
    return {
        'gm': gm_file,
        'wm': wm_file,
        'csf': csf_file,
        'ref': ref_file,
        'xfm': xfm_file,
        'shift': shift,
        'anat_affine': anat_affine,
        'ref_affine': ref_affine,
    }


def test_partial_volume_maps_transform_direction(pv_inputs, tmp_path):
    """PV maps land where resampling the TPM through the transform puts it."""
    results = pvc.PartialVolumeMaps(
        gm_tpm=str(pv_inputs['gm']),
        wm_tpm=str(pv_inputs['wm']),
        csf_tpm=str(pv_inputs['csf']),
        reference=str(pv_inputs['ref']),
        aslref2anat_xfm=str(pv_inputs['xfm']),
        metadata={'MRAcquisitionType': '2D'},
        effective_resolution='nominal',
    ).run(cwd=tmp_path)

    pv_gm = nb.load(results.outputs.gm_pv).get_fdata()
    assert pv_gm.max() > 0.5
    pv_csf = nb.load(results.outputs.csf_pv).get_fdata()
    csf_volume = nb.load(pv_inputs['csf']).get_fdata().sum()
    assert pv_csf.sum() * 64 == pytest.approx(csf_volume, rel=0.02)

    # Oracle: nitransforms pulls the anatomical image onto the ASL grid through the inverse
    # of the aslref-to-anat point mapping, as antsApplyTransforms does with the inverse flag.
    xfm = load_xfm(str(pv_inputs['xfm']))
    oracle = apply_xfm(
        ~xfm, nb.load(pv_inputs['gm']), reference=nb.load(pv_inputs['ref']), order=1
    ).get_fdata()

    com = _centre_of_mass_ras(pv_gm, pv_inputs['ref_affine'])
    com_oracle = _centre_of_mass_ras(oracle, pv_inputs['ref_affine'])
    np.testing.assert_allclose(com, com_oracle, atol=0.5)
    # And the shift moved the sphere away from its anatomical centre (the origin).
    assert np.linalg.norm(com) > 3

    # PV maps are tissue fractions: their sums match the tissue volumes.
    gm_volume = nb.load(pv_inputs['gm']).get_fdata().sum() * 1.0
    assert pv_gm.sum() * 64 == pytest.approx(gm_volume, rel=0.02)

    meta = results.outputs.pv_metadata
    assert meta['EffectiveResolutionSource'] == 'nominal'
    assert meta['EffectiveResolution'] == [4.0, 4.0, 4.0]
    assert meta['SliceEncodingAxis'] == 'k'


def test_partial_volume_maps_nominal_ignores_header_zooms(pv_inputs, tmp_path):
    """'nominal' uses the affine's voxel size, so header zooms that disagree add no blur."""
    ref_img = nb.load(pv_inputs['ref'])
    bad_ref = tmp_path / 'ref_bad_zooms.nii.gz'
    img = nb.Nifti1Image(np.asanyarray(ref_img.dataobj), ref_img.affine)
    for zooms, name in (((3.0, 3.0, 3.0), 'smaller'), ((5.0, 5.0, 5.0), 'larger')):
        img.header.set_zooms(zooms)
        img.to_filename(bad_ref)
        assert nb.load(bad_ref).header.get_zooms()[:3] == zooms
        (tmp_path / name).mkdir()
        results = pvc.PartialVolumeMaps(
            gm_tpm=str(pv_inputs['gm']),
            wm_tpm=str(pv_inputs['wm']),
            csf_tpm=str(pv_inputs['csf']),
            reference=str(bad_ref),
            aslref2anat_xfm=str(pv_inputs['xfm']),
            metadata={'MRAcquisitionType': '2D'},
            effective_resolution='nominal',
        ).run(cwd=tmp_path / name)
        meta = results.outputs.pv_metadata
        assert meta['EffectiveResolution'] == [4.0, 4.0, 4.0]
        assert meta['GaussianSigma'] == [0.0, 0.0, 0.0]


def test_partial_volume_maps_blur(pv_inputs, tmp_path):
    """A wider effective resolution spreads the tissue but keeps its volume."""
    kwargs = {
        'gm_tpm': str(pv_inputs['gm']),
        'wm_tpm': str(pv_inputs['wm']),
        'csf_tpm': str(pv_inputs['csf']),
        'reference': str(pv_inputs['ref']),
        'aslref2anat_xfm': str(pv_inputs['xfm']),
        'metadata': {'MRAcquisitionType': '3D', 'PulseSequenceType': 'GRASE'},
    }
    (tmp_path / 'nominal').mkdir()
    (tmp_path / 'auto').mkdir()
    nominal = pvc.PartialVolumeMaps(effective_resolution='nominal', **kwargs).run(
        cwd=tmp_path / 'nominal'
    )
    auto = pvc.PartialVolumeMaps(effective_resolution='auto', **kwargs).run(cwd=tmp_path / 'auto')
    a = nb.load(nominal.outputs.wm_pv).get_fdata()
    b = nb.load(auto.outputs.wm_pv).get_fdata()
    assert b.sum() == pytest.approx(a.sum(), rel=1e-3)
    assert b.max() < a.max()
    assert auto.outputs.pv_metadata['EffectiveResolutionSource'] == 'auto:grase-single-shot'
    np.testing.assert_allclose(
        auto.outputs.pv_metadata['EffectiveResolution'], [4.4, 4.4, 4 * 1.38]
    )


def test_local_regression_pvc_interface(tmp_path):
    """The interface unmixes CBF * M0 and calibrates with the centre voxel's M0."""
    rng = np.random.default_rng(1)
    shape = (16, 16, 4)
    affine = np.diag([3.0, 3.0, 5.0, 1.0])
    pv_gm = rng.uniform(0.1, 0.8, shape)
    pv_wm = (1 - pv_gm) * rng.uniform(0.5, 1.0, shape)
    pv_csf = 1 - pv_gm - pv_wm
    m0 = rng.uniform(800, 1200, shape)
    m0[0] = 0  # not quantifiable
    # Tissue amplitudes d_t = M0_t * f_t; uncorrected CBF = (sum_t pv_t d_t) / M0.
    # CSF carries a small spurious signal, which the CSF QC fit must find.
    d_gm, d_wm, d_csf = 60.0 * 1000, 20.0 * 800, 5.0 * 1500
    with np.errstate(divide='ignore', invalid='ignore'):
        cbf = np.nan_to_num((pv_gm * d_gm + pv_wm * d_wm + pv_csf * d_csf) / m0)
    quant_valid = m0 > 0
    mask = np.ones(shape, bool)

    files = {}
    for name, data in {
        'cbf': cbf,
        'm0_calib': m0,
        'quant_valid': quant_valid.astype(np.uint8),
        'mask': mask.astype(np.uint8),
        'gm_pv': pv_gm,
        'wm_pv': pv_wm,
        'csf_pv': pv_csf,
    }.items():
        files[name] = str(tmp_path / f'{name}.nii.gz')
        nb.Nifti1Image(data, affine).to_filename(files[name])

    results = pvc.LocalRegressionPVC(
        metadata={'MRAcquisitionType': '2D', 'SliceTiming': [0, 0.1, 0.2, 0.3]},
        **files,
    ).run(cwd=tmp_path)
    out = results.outputs
    status = nb.load(out.status).get_fdata()
    csf_status = nb.load(out.csf_status).get_fdata()
    csf_cbf = nb.load(out.csf_cbf).get_fdata()

    assert np.all(status[0] == 4)
    assert np.all(csf_status[0] == 4)
    # The three-tissue QC fit recovers the CSF signal, calibrated by the voxel's M0.
    full = (csf_status == 1) & (pv_csf >= 0.1)
    assert full.sum() > 100
    np.testing.assert_allclose(csf_cbf[full], d_csf / m0[full], rtol=1e-6)
    # The two-tissue fit assumes CSF is zero, so its GM estimate absorbs the CSF signal:
    # the QC fit is what exposes it, and the GM map must not come from the QC fit.
    gm = nb.load(out.gm_cbf).get_fdata()
    gm_valid = nb.load(out.gm_valid).get_fdata() > 0
    two = (status == 1) & gm_valid
    assert two.sum() > 100
    assert np.all(gm[~gm_valid] == 0)
    assert not np.allclose(gm[two], d_gm / m0[two], rtol=1e-3)
    assert sum(out.summary.values()) == pytest.approx(1.0)
    assert out.pvc_metadata['KernelSize'] == [5, 5, 1]
    assert out.qc['lr_csf_n_voxels'] == 0  # no voxel is half CSF in this phantom
    assert set(out.qc) == {
        'lr_csf_cbf_median',
        'lr_csf_cbf_iqr',
        'lr_csf_n_voxels',
        'lr_fraction_two_tissue',
    }


def test_local_regression_pvc_interface_exact_two_tissue(tmp_path):
    """Without CSF signal, the GM map is exact and the CSF QC estimate is zero."""
    rng = np.random.default_rng(3)
    shape = (16, 16, 2)
    affine = np.diag([3.0, 3.0, 5.0, 1.0])
    pv_gm = rng.uniform(0.1, 0.6, shape)
    pv_csf = rng.uniform(0.3, 0.6, shape) * (1 - pv_gm)
    pv_wm = 1 - pv_gm - pv_csf
    m0 = rng.uniform(800, 1200, shape)
    d_gm, d_wm = 60.0 * 1000, 20.0 * 800
    cbf = (pv_gm * d_gm + pv_wm * d_wm) / m0
    files = {}
    for name, data in {
        'cbf': cbf,
        'm0_calib': m0,
        'quant_valid': np.ones(shape, np.uint8),
        'mask': np.ones(shape, np.uint8),
        'gm_pv': pv_gm,
        'wm_pv': pv_wm,
        'csf_pv': pv_csf,
    }.items():
        files[name] = str(tmp_path / f'{name}.nii.gz')
        nb.Nifti1Image(data, affine).to_filename(files[name])
    out = (
        pvc.LocalRegressionPVC(metadata={'MRAcquisitionType': '3D'}, **files)
        .run(cwd=tmp_path)
        .outputs
    )
    status = nb.load(out.status).get_fdata()
    gm = nb.load(out.gm_cbf).get_fdata()
    two = (status == 1) & (nb.load(out.gm_valid).get_fdata() > 0)
    np.testing.assert_allclose(gm[two], d_gm / m0[two], rtol=1e-6)
    csf = nb.load(out.csf_cbf).get_fdata()
    csf_valid = nb.load(out.csf_status).get_fdata() == 1
    assert csf_valid.sum() > 100
    np.testing.assert_allclose(csf[csf_valid & (pv_csf >= 0.1)], 0, atol=1e-6)


def test_pvc_diagnostics_plot(tmp_path):
    """The diagnostics figure renders from the LR outputs.

    Its text keeps the default size when another plot in the same process left a larger
    global font size behind (as nireports' carpet plot does).
    """
    import matplotlib.pyplot as plt

    from aslprep.interfaces.plotting import PVCDiagnosticsPlot

    rng = np.random.default_rng(0)
    shape = (12, 12, 4)
    files = {}
    data = {
        'mean_cbf': rng.normal(50, 10, shape),
        'gm_cbf': rng.normal(70, 10, shape),
        'gm_valid': (rng.random(shape) > 0.3).astype(np.uint8),
        'gm_pv': rng.random(shape),
        'status': rng.integers(0, 5, shape).astype(np.uint8),
        'gm_noise': rng.random(shape),
    }
    for name, arr in data.items():
        files[name] = str(tmp_path / f'{name}.nii.gz')
        nb.Nifti1Image(arr, np.eye(4)).to_filename(files[name])

    with plt.rc_context({'font.size': 24}):
        results = PVCDiagnosticsPlot(
            summary={'two-tissue fit': 0.7, 'not identifiable': 0.3},
            pv_metadata={
                'EffectiveResolution': [3.0, 3.0, 4.0],
                'EffectiveResolutionSource': 'auto:2d',
            },
            **files,
        ).run(cwd=tmp_path)
        assert plt.rcParams['font.size'] == 24

    assert results.outputs.out_file.endswith('.svg')
    svg = open(results.outputs.out_file).read()
    assert '<svg' in svg[:2000]
    # Matplotlib draws text as glyph paths scaled by font size / 100
    title = svg.split('<!-- Fit status -->', 1)[1][:200]
    assert 'scale(0.12 -0.12)' in title


def test_divide_by_weight(tmp_path):
    """Values are divided by the resampled validity, and dropped where it is low."""
    affine = np.eye(4)
    data = np.array([[[45.0, 30.0, 0.0]]])
    weight = np.array([[[0.75, 0.5, 0.25]]])
    in_file = str(tmp_path / 'data.nii.gz')
    weight_file = str(tmp_path / 'weight.nii.gz')
    nb.Nifti1Image(data, affine).to_filename(in_file)
    nb.Nifti1Image(weight, affine).to_filename(weight_file)
    results = pvc.DivideByWeight(in_file=in_file, weight_file=weight_file).run(cwd=tmp_path)
    out = nb.load(results.outputs.out_file).get_fdata()
    np.testing.assert_allclose(out, [[[60.0, 60.0, 0.0]]])


def test_multidelay_local_regression_pvc_uses_calibration(tmp_path):
    """The tissue fits and the CSF diagnostic use the given calibration, not the metadata.

    The data are simulated with a labeling efficiency, blood T1 and partition coefficient that
    differ from what the metadata would give, so the true CBF values are only recovered if
    the interface quantifies with the specification.
    """
    from aslprep.utils.calibration import resolve_calibration
    from aslprep.utils.pvc import gkm_predict

    rng = np.random.default_rng(7)
    shape = (16, 16, 2)
    affine = np.diag([3.0, 3.0, 5.0, 1.0])
    plds = np.repeat(np.arange(0.25, 3.01, 0.25), 2)
    metadata = {
        'ArterialSpinLabelingType': 'PCASL',
        'MagneticFieldStrength': 3,
        'MRAcquisitionType': '3D',
        'LabelingDuration': 1.8,
        'PostLabelingDelay': plds.tolist(),
        'LabelingEfficiency': 0.85,
        'BackgroundSuppression': False,
    }
    calibration = resolve_calibration(metadata)
    params = calibration['Parameters']
    params['EffectiveLabelingEfficiency']['Value'] = 0.6
    params['BloodT1']['Value'] = 1.9
    params['PartitionCoefficient']['Value'] = 0.8
    kinetics = {
        'labeleff': 0.6,
        't1blood': 1.9,
        'partition_coefficient': 0.8,
        'is_casl': True,
        'tau': 1.8,
    }

    def curve(cbf, att, m0):
        return gkm_predict(
            np.array([cbf]),
            np.array([att]),
            np.array([0.5]),
            np.array([0.0]),
            np.array([m0]),
            plds[None, :],
            **kinetics,
        )[0]

    pv_gm = rng.uniform(0.2, 0.7, shape)
    pv_wm = (1 - pv_gm) * rng.uniform(0.6, 0.95, shape)
    pv_csf = 1 - pv_gm - pv_wm
    tissue = pv_gm[..., None] * curve(60.0, 1.0, 1000.0) + pv_wm[..., None] * curve(
        20.0, 1.6, 800.0
    )
    csf_signal = pv_csf[..., None] * curve(-10.0, 1.2, 1500.0)
    m0 = rng.uniform(700.0, 1100.0, shape)

    def run(deltam, name):
        run_dir = tmp_path / name
        run_dir.mkdir()
        files = {}
        for key, data in {
            'deltam': deltam,
            'm0_calib': m0,
            'quant_valid': np.ones(shape, np.uint8),
            'mask': np.ones(shape, np.uint8),
            'gm_pv': pv_gm,
            'wm_pv': pv_wm,
            'csf_pv': pv_csf,
            'mean_cbf': np.full(shape, 45.0),
            'att': np.full(shape, 1.2),
        }.items():
            files[key] = str(run_dir / f'{key}.nii.gz')
            nb.Nifti1Image(data, affine).to_filename(files[key])
        interface = pvc.MultiDelayLocalRegressionPVC(
            metadata=metadata, calibration=calibration, **files
        )
        return interface.run(cwd=run_dir).outputs

    # Tissue fits, without CSF signal (the two-tissue regression does not model it)
    out = run(tissue, 'tissue')
    status = nb.load(out.status).get_fdata()
    good = (status == 1) & (nb.load(out.gm_kinstatus).get_fdata() == 1)
    assert good.sum() > 100
    gm = nb.load(out.gm_cbf).get_fdata()
    np.testing.assert_allclose(gm[good], 60.0 * 1000.0 / m0[good], rtol=0.01)

    # The CSF diagnostic, whose basis is predicted with the same parameters
    out = run(tissue + csf_signal, 'csf')
    csf = nb.load(out.csf_cbf).get_fdata()
    csf_ok = (nb.load(out.csf_status).get_fdata() == 1) & (pv_csf >= 0.1)
    assert csf_ok.sum() > 50
    np.testing.assert_allclose(csf[csf_ok], -10.0 * 1500.0 / m0[csf_ok], rtol=0.01)


def test_multidelay_local_regression_pvc_interface(tmp_path, monkeypatch):
    """Multi-delay LR PVC unmixes each delay and fits each tissue curve (spec 5.7).

    M0 varies over the image, so each estimate must be calibrated by its own voxel's M0. A
    second run adds a negative CSF curve, which the signed CSF diagnostic must recover, also
    where the uncorrected CBF is zero (its ATT then falls back to the median). A third run
    forces each kinetic status, to check which estimates are written.
    """
    from aslprep.utils.pvc import gkm_predict

    rng = np.random.default_rng(5)
    shape = (16, 16, 2)
    affine = np.diag([3.0, 3.0, 5.0, 1.0])
    plds = np.repeat(np.arange(0.25, 3.01, 0.25), 2)  # two repeats of each delay
    kinetics = {'labeleff': 0.85, 't1blood': 1.65, 'is_casl': True, 'tau': 1.8}
    metadata = {
        'ArterialSpinLabelingType': 'PCASL',
        'MagneticFieldStrength': 3,
        'MRAcquisitionType': '3D',
        'LabelingDuration': 1.8,
        'PostLabelingDelay': plds.tolist(),
        'LabelingEfficiency': 0.85,
    }
    pv_gm = rng.uniform(0.2, 0.7, shape)
    pv_wm = (1 - pv_gm) * rng.uniform(0.6, 0.95, shape)
    pv_csf = 1 - pv_gm - pv_wm
    m0_gm, m0_wm, m0_csf = 1000.0, 800.0, 1500.0

    def curve(cbf, att, m0):
        return gkm_predict(
            np.array([cbf]),
            np.array([att]),
            np.array([0.5]),
            np.array([0.0]),
            np.array([m0]),
            plds[None, :],
            **kinetics,
        )[0]

    tissue = pv_gm[..., None] * curve(60.0, 1.0, m0_gm) + pv_wm[..., None] * curve(
        20.0, 1.6, m0_wm
    )
    csf_signal = pv_csf[..., None] * curve(-10.0, 1.2, m0_csf)
    m0 = rng.uniform(700.0, 1100.0, shape)
    m0[0, 0] = 0  # not quantifiable
    mean_cbf = np.full(shape, 45.0)
    att = np.full(shape, 1.2)
    mean_cbf[8:] = 0  # no uncorrected CBF, so no uncorrected ATT either
    att[8:] = 0

    def run(deltam, name):
        run_dir = tmp_path / name
        run_dir.mkdir()
        files = {}
        for key, data in {
            'deltam': deltam,
            'm0_calib': m0,
            'quant_valid': (m0 > 0).astype(np.uint8),
            'mask': np.ones(shape, np.uint8),
            'gm_pv': pv_gm,
            'wm_pv': pv_wm,
            'csf_pv': pv_csf,
            'mean_cbf': mean_cbf,
            'att': att,
        }.items():
            files[key] = str(run_dir / f'{key}.nii.gz')
            nb.Nifti1Image(data, affine).to_filename(files[key])
        interface = pvc.MultiDelayLocalRegressionPVC(metadata=metadata, **files)
        return interface.run(cwd=run_dir).outputs

    out = run(tissue, 'tissue')
    status = nb.load(out.status).get_fdata()
    gm_kinstatus = nb.load(out.gm_kinstatus).get_fdata()
    gm = nb.load(out.gm_cbf).get_fdata()
    gm_att = nb.load(out.gm_att).get_fdata()
    wm_att = nb.load(out.wm_att).get_fdata()
    good = (status == 1) & (gm_kinstatus == 1)
    assert good.sum() > 100
    # The estimand: tissue CBF calibrated by the centre voxel's M0
    np.testing.assert_allclose(gm[good], 60.0 * m0_gm / m0[good], rtol=0.01)
    np.testing.assert_allclose(gm_att[good], 1.0, atol=0.01)
    wm_good = (status == 1) & (nb.load(out.wm_kinstatus).get_fdata() == 1)
    assert wm_good.sum() > 100
    np.testing.assert_allclose(wm_att[wm_good], 1.6, atol=0.01)
    assert np.all(status[0, 0] == 4)
    assert np.all(gm[0, 0] == 0)
    # No CSF signal: the signed CSF diagnostic is zero where it is estimated.
    csf = nb.load(out.csf_cbf).get_fdata()
    csf_status = nb.load(out.csf_status).get_fdata()
    csf_ok = (csf_status == 1) & (pv_csf >= 0.1)
    assert csf_ok.sum() > 50
    np.testing.assert_allclose(csf[csf_ok], 0, atol=1e-6)
    assert out.pvc_metadata['ResidualUnits'] == 'delta-M'
    assert out.pvc_metadata['CSFQualityControlDefaultATT'] == pytest.approx(1.2)
    assert set(out.qc) == {
        'lr_csf_cbf_median',
        'lr_csf_cbf_iqr',
        'lr_csf_n_voxels',
        'lr_fraction_two_tissue',
    }

    # A negative CSF curve is recovered, with its sign, in both halves of the image.
    out = run(tissue + csf_signal, 'csf')
    csf = nb.load(out.csf_cbf).get_fdata()
    csf_status = nb.load(out.csf_status).get_fdata()
    csf_ok = (csf_status == 1) & (pv_csf >= 0.1)
    assert csf_ok[:8].sum() > 25
    assert csf_ok[8:].sum() > 25
    np.testing.assert_allclose(csf[csf_ok], -10.0 * m0_csf / m0[csf_ok], rtol=0.01)

    # Kinetic statuses: CBF is written for 1 and 2 only, ATT for 1 only.
    real_fit = pvc.fit_tissue_curves

    def cycled_fit(*args, **kwargs):
        result = real_fit(*args, **kwargs)
        result['status'] = np.resize(np.array([1, 2, 3, 5], np.uint8), result['status'].shape)
        return result

    monkeypatch.setattr(pvc, 'fit_tissue_curves', cycled_fit)
    out = run(tissue, 'statuses')
    for tissue_name in ('gm', 'wm'):
        kinstatus = nb.load(getattr(out, f'{tissue_name}_kinstatus')).get_fdata()
        cbf = nb.load(getattr(out, f'{tissue_name}_cbf')).get_fdata()
        cbf_valid = nb.load(getattr(out, f'{tissue_name}_valid')).get_fdata()
        att_out = nb.load(getattr(out, f'{tissue_name}_att')).get_fdata()
        att_valid = nb.load(getattr(out, f'{tissue_name}_att_valid')).get_fdata()
        assert set(np.unique(kinstatus)) >= {1, 2, 3, 5}
        rejected = np.isin(kinstatus, (3, 5))
        assert np.all(cbf[rejected] == 0)
        assert np.all(cbf_valid[rejected] == 0)
        np.testing.assert_array_equal(cbf_valid, np.isin(kinstatus, (1, 2)))
        np.testing.assert_array_equal(att_valid, kinstatus == 1)
        assert np.all(att_out[kinstatus != 1] == 0)
        assert np.all(att_out[kinstatus == 1] > 0)
