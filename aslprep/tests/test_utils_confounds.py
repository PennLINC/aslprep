"""Tests for aslprep.utils.confounds.compute_cbf_threshold_stats."""

import json
import os

import nibabel as nb
import numpy as np
import pandas as pd
import pytest

from aslprep.utils.confounds import (
    compute_cbf_threshold_stats,
    compute_clipped_voxel_percentage,
)


def _save_img(arr, tmpdir, fname):
    """Persist a numpy array as a NIfTI-1 image and return the file path."""
    img = nb.Nifti1Image(arr, affine=np.eye(4))
    fpath = os.path.join(str(tmpdir), fname)
    img.to_filename(fpath)
    return fpath


def test_threshold_stats_basic():
    """Test compute_cbf_threshold_stats with a simple 3D array."""
    # 10 voxels: 10, 20, 30, ..., 100
    cbf = np.arange(10, 110, 10, dtype=float).reshape(2, 5, 1)
    mask = np.ones_like(cbf, dtype=np.uint8)

    result = compute_cbf_threshold_stats(cbf, mask, thresholds=(50,))
    # 60, 70, 80, 90, 100 are > 50, so 5 out of 10
    assert result['perc_voxels_cbf_gt_50'] == pytest.approx(50.0)


def test_threshold_stats_4d_input():
    """A 4D CBF image should be averaged across time before thresholding."""
    vol_lo = np.full((2, 2, 2), 80.0)
    vol_hi = np.full((2, 2, 2), 120.0)
    cbf_4d = np.stack([vol_lo, vol_hi], axis=3)  # mean = 100.0 everywhere
    mask = np.ones((2, 2, 2), dtype=np.uint8)

    result = compute_cbf_threshold_stats(cbf_4d, mask, thresholds=(100,))
    assert result['perc_voxels_cbf_gt_100'] == pytest.approx(0.0)  # 100 is NOT > 100

    result = compute_cbf_threshold_stats(cbf_4d, mask, thresholds=(99,))
    assert result['perc_voxels_cbf_gt_99'] == pytest.approx(100.0)


def test_threshold_stats_nan_voxels():
    """NaN voxels should be excluded from both numerator and denominator."""
    cbf = np.array([200.0, np.nan, 50.0, np.nan, 150.0]).reshape(5, 1, 1)
    mask = np.ones((5, 1, 1), dtype=np.uint8)

    result = compute_cbf_threshold_stats(cbf, mask, thresholds=(100,))
    # 3 valid voxels (200, 50, 150), 2 of them > 100
    assert result['perc_voxels_cbf_gt_100'] == pytest.approx(200.0 / 3.0)


def test_threshold_stats_empty_mask():
    """When no voxels are inside the mask, all results should be NaN."""
    cbf = np.array([100.0, 200.0, 300.0]).reshape(3, 1, 1)
    mask = np.zeros((3, 1, 1), dtype=np.uint8)

    result = compute_cbf_threshold_stats(cbf, mask)
    for val in result.values():
        assert np.isnan(val)


def test_threshold_stats_all_above():
    """If every voxel is above every threshold, all percentages should be 100."""
    cbf = np.full((3, 3, 3), 250.0)
    mask = np.ones_like(cbf, dtype=np.uint8)

    result = compute_cbf_threshold_stats(cbf, mask)
    assert result['perc_voxels_cbf_gt_100'] == pytest.approx(100.0)
    assert result['perc_voxels_cbf_gt_150'] == pytest.approx(100.0)
    assert result['perc_voxels_cbf_gt_200'] == pytest.approx(100.0)


def test_threshold_stats_all_below():
    """If every voxel is below every threshold, all percentages should be 0."""
    cbf = np.full((3, 3, 3), 10.0)
    mask = np.ones_like(cbf, dtype=np.uint8)

    result = compute_cbf_threshold_stats(cbf, mask)
    assert result['perc_voxels_cbf_gt_100'] == pytest.approx(0.0)
    assert result['perc_voxels_cbf_gt_150'] == pytest.approx(0.0)
    assert result['perc_voxels_cbf_gt_200'] == pytest.approx(0.0)


def test_threshold_stats_custom_thresholds():
    """Non-default threshold values should work the same way."""
    cbf = np.array([10.0, 30.0, 50.0, 70.0]).reshape(2, 2, 1)
    mask = np.ones((2, 2, 1), dtype=np.uint8)

    result = compute_cbf_threshold_stats(cbf, mask, thresholds=(20, 60))
    assert result['perc_voxels_cbf_gt_20'] == pytest.approx(75.0)  # 30, 50, 70
    assert result['perc_voxels_cbf_gt_60'] == pytest.approx(25.0)  # 70


def test_threshold_stats_strict_gt():
    """Values exactly equal to the threshold must NOT be counted."""
    cbf = np.array([100.0, 100.0, 100.0, 101.0]).reshape(2, 2, 1)
    mask = np.ones((2, 2, 1), dtype=np.uint8)

    result = compute_cbf_threshold_stats(cbf, mask, thresholds=(100,))
    # only 101 is strictly > 100
    assert result['perc_voxels_cbf_gt_100'] == pytest.approx(25.0)


def test_threshold_stats_from_nifti_paths(tmp_path):
    """The function should also accept NIfTI file paths, not just arrays."""
    cbf_data = np.full((5, 5, 5), 160.0)
    mask_data = np.ones((5, 5, 5), dtype=np.uint8)

    cbf_file = _save_img(cbf_data, tmp_path, 'cbf.nii.gz')
    mask_file = _save_img(mask_data, tmp_path, 'mask.nii.gz')

    result = compute_cbf_threshold_stats(cbf_file, mask_file)
    assert result['perc_voxels_cbf_gt_100'] == pytest.approx(100.0)
    assert result['perc_voxels_cbf_gt_150'] == pytest.approx(100.0)
    assert result['perc_voxels_cbf_gt_200'] == pytest.approx(0.0)


def test_cbf_qc_interface_threshold_columns(tmp_path):
    """ComputeCBFQC should write threshold columns to the output TSV and metadata JSON."""
    from aslprep.interfaces.confounds import ComputeCBFQC

    shape = (10, 10, 10)

    # half the brain at 50, half at 160
    cbf_data = np.full(shape, 50.0, dtype=np.float32)
    cbf_data[5:, :, :] = 160.0
    mean_cbf = _save_img(cbf_data, tmp_path, 'mean_cbf.nii.gz')

    # tissue probability maps
    gm = np.zeros(shape, dtype=np.float32)
    gm[2:8, 2:8, 2:8] = 0.9
    wm = np.zeros(shape, dtype=np.float32)
    wm[0:2, :, :] = 0.9
    csf = np.zeros(shape, dtype=np.float32)
    csf[8:, :, :] = 0.9

    mask = np.ones(shape, dtype=np.uint8)

    # minimal confounds file
    confounds_df = pd.DataFrame(
        {
            'framewise_displacement': [0.1, 0.2, 0.15],
            'rmsd': [0.05, 0.06, 0.04],
        }
    )
    confounds_file = os.path.join(str(tmp_path), 'confounds.tsv')
    confounds_df.to_csv(confounds_file, sep='\t', index=False)

    # name_source must look like a BIDS filename (used for entity extraction)
    name_source = _save_img(np.zeros(shape, dtype=np.float32), tmp_path, 'sub-01_asl.nii.gz')

    interface = ComputeCBFQC(
        name_source=name_source,
        mean_cbf=mean_cbf,
        gm_tpm=_save_img(gm, tmp_path, 'gm.nii.gz'),
        wm_tpm=_save_img(wm, tmp_path, 'wm.nii.gz'),
        csf_tpm=_save_img(csf, tmp_path, 'csf.nii.gz'),
        asl_mask=_save_img(mask, tmp_path, 'asl_mask.nii.gz'),
        t1w_mask=_save_img(mask, tmp_path, 't1w_mask.nii.gz'),
        confounds_file=confounds_file,
    )
    results = interface.run(cwd=str(tmp_path))

    # --- check the TSV ---
    qc_df = pd.read_csv(results.outputs.qc_file, sep='\t')

    # mean CBF columns should be present and non-NaN
    for thresh in (100, 150, 200):
        col = f'perc_voxels_cbf_gt_{thresh}'
        assert col in qc_df.columns, f'{col} missing from QC TSV'
        assert not pd.isna(qc_df[col].iloc[0]), f'{col} should not be NaN'

    # SCORE/SCRUB/BASIL columns should exist (NaN since those maps were not provided)
    for variant in ('score', 'scrub', 'basil', 'basil_gm'):
        for thresh in (100, 150, 200):
            col = f'perc_voxels_cbf_{variant}_gt_{thresh}'
            assert col in qc_df.columns, f'{col} missing from QC TSV'

    # spot-check the actual values (half at 50, half at 160)
    assert qc_df['perc_voxels_cbf_gt_100'].iloc[0] == pytest.approx(50.0)
    assert qc_df['perc_voxels_cbf_gt_150'].iloc[0] == pytest.approx(50.0)
    assert qc_df['perc_voxels_cbf_gt_200'].iloc[0] == pytest.approx(0.0)

    # --- check the metadata JSON ---
    with open(results.outputs.qc_metadata) as f:
        metadata = json.load(f)

    for thresh in (100, 150, 200):
        key = f'perc_voxels_cbf_gt_{thresh}'
        assert key in metadata, f'{key} missing from metadata JSON'
        assert metadata[key]['Units'] == 'percent'


@pytest.mark.parametrize(
    ('name', 'expected'),
    [
        ('FramewiseDisplacement', 'framewise_displacement'),
        ('StdDVARS', 'std_dvars'),
        ('NonStdDVARS', 'non_std_dvars'),
        ('vx-wisestdDVARS', 'vx-wisestd_dvars'),
        ('DVARS', 'dvars'),
        ('RMSD', 'rmsd'),
        ('aCompCor00', 'a_comp_cor00'),
        ('WhiteMatter', 'white_matter'),
        ('GlobalSignal', 'global_signal'),
        ('RotX', 'rot_x'),
        ('X', 'x'),
        ('global_signal', 'global_signal'),
        ('trans_x', 'trans_x'),
    ],
)
def test_camel_to_snake(name, expected):
    """Confound column names are converted to snake_case."""
    from aslprep.utils.confounds import _camel_to_snake

    assert _camel_to_snake(name) == expected


def _clipping_series(tmp_path, clipped_fraction, shape=(10, 10, 10, 4), ceiling=4095):
    """Write a series with noise below a ceiling, and a block of voxels clipped at it."""
    rng = np.random.default_rng(0)
    data = rng.integers(100, 3000, size=shape).astype(np.int16)
    n_clipped = round(clipped_fraction * shape[0])
    # Clip the block in one volume only, as a voxel is unusable once any volume is clipped.
    data[:n_clipped, :, :, 1] = ceiling
    in_file = os.path.join(str(tmp_path), 'series.nii.gz')
    nb.Nifti1Image(data, np.eye(4)).to_filename(in_file)
    return in_file


@pytest.mark.parametrize(
    ('clipped_fraction', 'expected'),
    [(0, 0.1), (0.2, 20.0), (1, 100.0)],
)
def test_clipped_voxel_percentage(tmp_path, clipped_fraction, expected):
    """The percentage reflects the share of brain voxels at the maximum intensity."""
    in_file = _clipping_series(tmp_path, clipped_fraction)
    mask = _save_img(np.ones((10, 10, 10), dtype=np.uint8), tmp_path, 'mask.nii.gz')

    perc = compute_clipped_voxel_percentage(in_file, mask)
    if clipped_fraction == 0:
        # Without clipping, only the voxel(s) holding the maximum are counted.
        assert 0 < perc <= expected
    else:
        assert perc == pytest.approx(expected)


def test_clipped_voxel_percentage_outside_mask(tmp_path):
    """Voxels at the maximum intensity outside of the mask are not counted."""
    in_file = _clipping_series(tmp_path, 0.2)
    mask_data = np.zeros((10, 10, 10), dtype=np.uint8)
    mask_data[5:, :, :] = 1
    mask = _save_img(mask_data, tmp_path, 'mask.nii.gz')

    assert compute_clipped_voxel_percentage(in_file, mask) == 0


def test_clipped_voxel_percentage_3d_and_scaled(tmp_path):
    """A 3D image with a scale factor (for example, an M0 scan) is handled."""
    data = np.full((10, 10, 10), 1000, dtype=np.int16)
    data[:3] = 4095
    img = nb.Nifti1Image(data, np.eye(4))
    img.header.set_slope_inter(2.5, 0)
    in_file = os.path.join(str(tmp_path), 'm0scan.nii.gz')
    img.to_filename(in_file)
    mask = _save_img(np.ones((10, 10, 10), dtype=np.uint8), tmp_path, 'mask.nii.gz')

    assert compute_clipped_voxel_percentage(in_file, mask) == pytest.approx(30.0)


def test_clipped_voxel_percentage_resamples_mask(tmp_path):
    """A mask on a different grid is resampled to the image."""
    # The image has 2 mm voxels; the mask covers the same field of view with 4 mm voxels.
    data = np.full((10, 10, 10), 1000, dtype=np.int16)
    data[:5] = 4095
    in_file = os.path.join(str(tmp_path), 'm0scan.nii.gz')
    nb.Nifti1Image(data, np.diag([2, 2, 2, 1])).to_filename(in_file)
    mask_file = os.path.join(str(tmp_path), 'mask.nii.gz')
    mask_affine = np.diag([4.0, 4.0, 4.0, 1.0])
    mask_affine[:3, 3] = 1
    nb.Nifti1Image(np.ones((5, 5, 5), dtype=np.uint8), mask_affine).to_filename(mask_file)

    assert compute_clipped_voxel_percentage(in_file, mask_file) == pytest.approx(50.0)


def test_clipped_voxel_percentage_empty_mask(tmp_path):
    """An empty mask gives NaN, not a division by zero."""
    in_file = _clipping_series(tmp_path, 0.2)
    mask = _save_img(np.zeros((10, 10, 10), dtype=np.uint8), tmp_path, 'mask.nii.gz')

    assert np.isnan(compute_clipped_voxel_percentage(in_file, mask))


@pytest.mark.parametrize('with_inputs', [True, False])
def test_cbf_qc_interface_clipping_columns(tmp_path, with_inputs):
    """ComputeCBFQC reports clipping for the ASL file and the M0 scan when they are given."""
    from aslprep.interfaces.confounds import ComputeCBFQC

    shape = (10, 10, 10)
    cbf_data = np.full(shape, 50.0, dtype=np.float32)
    cbf_data[5:, :, :] = 60.0
    gm = np.zeros(shape, dtype=np.float32)
    gm[2:8, 2:8, 2:8] = 0.9
    wm = np.zeros(shape, dtype=np.float32)
    wm[0:2, :, :] = 0.9
    csf = np.zeros(shape, dtype=np.float32)
    csf[8:, :, :] = 0.9
    mask = np.ones(shape, dtype=np.uint8)
    confounds_file = os.path.join(str(tmp_path), 'confounds.tsv')
    pd.DataFrame({'framewise_displacement': [0.1, 0.2], 'rmsd': [0.05, 0.06]}).to_csv(
        confounds_file,
        sep='\t',
        index=False,
    )
    m0_data = np.full(shape, 1000, dtype=np.int16)
    m0_data[:4] = 4095

    kwargs = {}
    if with_inputs:
        kwargs = {
            'asl_file': _clipping_series(tmp_path, 0.2),
            'm0scan': _save_img(m0_data, tmp_path, 'sub-01_m0scan.nii.gz'),
        }

    results = ComputeCBFQC(
        name_source=_save_img(np.zeros(shape, dtype=np.float32), tmp_path, 'sub-01_asl.nii.gz'),
        mean_cbf=_save_img(cbf_data, tmp_path, 'mean_cbf.nii.gz'),
        gm_tpm=_save_img(gm, tmp_path, 'gm.nii.gz'),
        wm_tpm=_save_img(wm, tmp_path, 'wm.nii.gz'),
        csf_tpm=_save_img(csf, tmp_path, 'csf.nii.gz'),
        asl_mask=_save_img(mask, tmp_path, 'asl_mask.nii.gz'),
        t1w_mask=_save_img(mask, tmp_path, 't1w_mask.nii.gz'),
        confounds_file=confounds_file,
        **kwargs,
    ).run(cwd=str(tmp_path))

    qc_df = pd.read_csv(results.outputs.qc_file, sep='\t')
    if with_inputs:
        assert qc_df['perc_voxels_clipped'].iloc[0] == pytest.approx(20.0)
        assert qc_df['perc_voxels_clipped_m0scan'].iloc[0] == pytest.approx(40.0)
    else:
        assert pd.isna(qc_df['perc_voxels_clipped'].iloc[0])
        assert pd.isna(qc_df['perc_voxels_clipped_m0scan'].iloc[0])

    with open(results.outputs.qc_metadata) as f:
        metadata = json.load(f)

    for key in ('perc_voxels_clipped', 'perc_voxels_clipped_m0scan'):
        assert metadata[key]['Units'] == 'percent'


@pytest.mark.parametrize('with_pvc', [True, False])
def test_cbf_qc_interface_pvc_columns(tmp_path, with_pvc):
    """Local-regression PVC QC measures are added to the QC file only when it ran."""
    from aslprep.interfaces.confounds import PVC_QC_METADATA, ComputeCBFQC

    shape = (10, 10, 10)
    cbf_data = np.full(shape, 50.0, dtype=np.float32)
    gm = np.zeros(shape, dtype=np.float32)
    gm[2:8, 2:8, 2:8] = 0.9
    wm = np.zeros(shape, dtype=np.float32)
    wm[0:2, :, :] = 0.9
    csf = np.zeros(shape, dtype=np.float32)
    csf[8:, :, :] = 0.9
    mask = np.ones(shape, dtype=np.uint8)
    confounds_file = os.path.join(str(tmp_path), 'confounds.tsv')
    pd.DataFrame({'framewise_displacement': [0.1, 0.2], 'rmsd': [0.05, 0.06]}).to_csv(
        confounds_file, sep='\t', index=False
    )
    pvc_qc = {
        'lr_csf_cbf_median': 2.5,
        'lr_csf_cbf_iqr': 4.0,
        'lr_csf_n_voxels': 120,
        'lr_fraction_two_tissue': 0.9,
    }
    kwargs = {'pvc_qc': pvc_qc} if with_pvc else {}
    results = ComputeCBFQC(
        name_source=_save_img(np.zeros(shape, dtype=np.float32), tmp_path, 'sub-01_asl.nii.gz'),
        mean_cbf=_save_img(cbf_data, tmp_path, 'mean_cbf.nii.gz'),
        gm_tpm=_save_img(gm, tmp_path, 'gm.nii.gz'),
        wm_tpm=_save_img(wm, tmp_path, 'wm.nii.gz'),
        csf_tpm=_save_img(csf, tmp_path, 'csf.nii.gz'),
        asl_mask=_save_img(mask, tmp_path, 'asl_mask.nii.gz'),
        t1w_mask=_save_img(mask, tmp_path, 't1w_mask.nii.gz'),
        confounds_file=confounds_file,
        **kwargs,
    ).run(cwd=str(tmp_path))

    qc_df = pd.read_csv(results.outputs.qc_file, sep='\t')
    with open(results.outputs.qc_metadata) as f:
        metadata = json.load(f)
    for key, value in pvc_qc.items():
        if with_pvc:
            assert qc_df[key].iloc[0] == pytest.approx(value)
            assert metadata[key] == PVC_QC_METADATA[key]
        else:
            assert key not in qc_df.columns
            assert key not in metadata
