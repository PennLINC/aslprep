"""Tests for helper functions in aslprep.workflows.asl.confounds."""

import nibabel as nb
import numpy as np

from aslprep.workflows.asl.confounds import _clip_mask_to_coverage


def _write(arr, path):
    nb.Nifti1Image(arr, np.eye(4)).to_filename(path)
    return str(path)


def test_clip_mask_to_coverage(tmp_path, monkeypatch):
    """Crown mask voxels outside the ASL coverage are removed."""
    monkeypatch.chdir(tmp_path)
    mask = np.ones((4, 4, 4), dtype=np.uint8)
    mask_file = _write(mask, tmp_path / 'mask.nii.gz')

    # 4D data: only the first two slices have signal that varies over time
    asl = np.zeros((4, 4, 4, 4), dtype=np.float32)
    asl[:, :, :2, :] = np.array([1, 2, 1, 2], dtype=np.float32)
    asl_file = _write(asl, tmp_path / 'asl.nii.gz')
    clipped = nb.load(_clip_mask_to_coverage(mask_file, asl_file)).get_fdata()
    assert clipped[:, :, :2].all()
    assert not clipped[:, :, 2:].any()

    # Single-volume data: fall back to nonzero signal
    asl_1vol = asl[..., :1]
    asl_1vol_file = _write(asl_1vol, tmp_path / 'asl_1vol.nii.gz')
    clipped = nb.load(_clip_mask_to_coverage(mask_file, asl_1vol_file)).get_fdata()
    assert clipped[:, :, :2].all()
    assert not clipped[:, :, 2:].any()
