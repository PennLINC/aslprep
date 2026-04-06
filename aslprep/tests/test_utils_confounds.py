"""Tests for aslprep.utils.confounds.compute_cbf_threshold_stats."""

import json
import os

import nibabel as nb
import numpy as np
import pandas as pd
import pytest

from aslprep.utils.confounds import compute_cbf_threshold_stats


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
