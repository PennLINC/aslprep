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
    """The diagnostics figure renders from the LR outputs."""
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
    results = PVCDiagnosticsPlot(
        summary={'two-tissue fit': 0.7, 'not identifiable': 0.3},
        pv_metadata={
            'EffectiveResolution': [3.0, 3.0, 4.0],
            'EffectiveResolutionSource': 'auto:2d',
        },
        **files,
    ).run(cwd=tmp_path)
    assert results.outputs.out_file.endswith('.svg')
    assert '<svg' in open(results.outputs.out_file).read(2000)


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
