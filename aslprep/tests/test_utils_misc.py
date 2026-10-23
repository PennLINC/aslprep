"""Tests for aslprep.utils.misc."""

import nibabel as nb
import numpy as np
import pytest

from aslprep.utils.misc import estimate_series_mem_gb


def test_estimate_series_mem_gb(tmp_path):
    """The estimate holds four 8-byte copies of the series, with a floor of 0.2 GB."""
    shape = (64, 64, 36, 50)
    asl_file = tmp_path / 'sub-01_asl.nii.gz'
    nb.Nifti1Image(np.zeros(shape, dtype=np.uint8), np.eye(4)).to_filename(asl_file)

    # Whatever the data type on disk, tools may work in 8-byte floats
    assert estimate_series_mem_gb(str(asl_file)) == pytest.approx(4 * 8 * np.prod(shape) / 1024**3)

    small_file = tmp_path / 'sub-02_asl.nii.gz'
    nb.Nifti1Image(np.zeros((4, 4, 4, 2), dtype=np.float32), np.eye(4)).to_filename(small_file)
    assert estimate_series_mem_gb(str(small_file)) == 0.2


def test_estimate_series_mem_gb_placeholder(tmp_path):
    """An unreadable series, such as an empty placeholder, gets the default."""
    asl_file = tmp_path / 'sub-01_asl.nii.gz'
    asl_file.touch()
    assert estimate_series_mem_gb(str(asl_file)) == 1.0
    assert estimate_series_mem_gb(str(asl_file), default=0.5) == 0.5
