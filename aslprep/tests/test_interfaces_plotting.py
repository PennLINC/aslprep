"""Tests for aslprep.interfaces.plotting."""

import nibabel as nb
import numpy as np
import pytest

from aslprep.interfaces.plotting import ATTStatusPlot
from aslprep.utils.cbf import ATT_STATUS


@pytest.mark.parametrize('slice_axis', [0, 2])
@pytest.mark.parametrize('with_seg', [True, False])
def test_att_status_plot(tmp_path, slice_axis, with_seg):
    """The ATT status figure renders, with or without a segmentation."""
    rng = np.random.default_rng(0)
    shape = (12, 10, 8)
    mask = np.zeros(shape, np.uint8)
    mask[2:-2, 2:-2, 1:-1] = 1
    data = {
        'att': rng.uniform(0.5, 2.0, shape),
        'earliest_delay': rng.uniform(0.5, 2.0, shape),
        'att_status': rng.choice(list(ATT_STATUS.values()), shape).astype(np.uint8),
        'ref_vol': rng.uniform(0, 1000, shape),
        'mask': mask,
    }
    if with_seg:
        data['seg_file'] = rng.integers(0, 4, shape).astype(np.uint8)

    files = {}
    for name, arr in data.items():
        files[name] = str(tmp_path / f'{name}.nii.gz')
        nb.Nifti1Image(arr, np.eye(4)).to_filename(files[name])

    results = ATTStatusPlot(slice_axis=slice_axis, **files).run(cwd=tmp_path)
    assert results.outputs.out_file.endswith('.svg')
    assert '<svg' in open(results.outputs.out_file).read(2000)
