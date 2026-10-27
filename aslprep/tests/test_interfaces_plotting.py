"""Tests for aslprep.interfaces.plotting."""

import nibabel as nb
import numpy as np
import pytest

from aslprep.interfaces.plotting import (
    ATTStatusPlot,
    CBFByTissueTypePlot,
    MultiPLDByTissueTypePlot,
)
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


def _write_maps(tmp_path, shape=(12, 12, 6)):
    """Write a segmentation and three CBF maps, one with zeros outside gray matter."""
    rng = np.random.default_rng(0)
    seg = rng.integers(0, 4, shape).astype(np.uint8)
    affine = np.eye(4)
    seg_file = str(tmp_path / 'dseg.nii.gz')
    nb.Nifti1Image(seg, affine).to_filename(seg_file)

    means = np.array([0, 60, 25, 5], dtype=float)[seg]
    maps = {
        'basic': means + rng.normal(0, 8, shape),
        'score': means + rng.normal(0, 6, shape),
        'gmpvc': np.where(seg == 1, means + rng.normal(0, 8, shape), 0),
    }
    files = []
    for name, data in maps.items():
        files.append(str(tmp_path / f'{name}.nii.gz'))
        nb.Nifti1Image(data, affine).to_filename(files[-1])

    return seg_file, files


def test_cbf_by_tissue_type_plot(tmp_path):
    """The ridge plot renders one row per map, including a map that is empty outside GM."""
    seg_file, files = _write_maps(tmp_path)
    results = CBFByTissueTypePlot(
        in_files=files,
        labels=['Basic', 'SCORE', 'BASIL GM-PVC'],
        seg_file=seg_file,
    ).run(cwd=tmp_path)
    out_file = results.outputs.out_file
    assert out_file.endswith('.svg')
    text = open(out_file).read()
    assert '<svg' in text[:2000]
    for label in ['Basic', 'SCORE', 'BASIL GM-PVC']:
        assert label in text


def test_cbf_by_tissue_type_plot_label_mismatch(tmp_path):
    """Each map needs a label."""
    seg_file, files = _write_maps(tmp_path)
    with pytest.raises(ValueError, match='labels'):
        CBFByTissueTypePlot(in_files=files, labels=['Basic'], seg_file=seg_file).run(cwd=tmp_path)


@pytest.mark.parametrize(
    ('fields', 'panels'),
    [
        # Basic model only: one row per panel
        (['att', 'abat', 'abv'], ['Transit', 'Arrival', 'Volume']),
        # With BASIL, which may not write its aBV map
        (['att', 'abat', 'abv', 'att_basil'], ['Transit', 'Arrival', 'Volume']),
        (['att', 'abat', 'abv', 'att_basil', 'abv_basil'], ['Transit', 'Arrival', 'Volume']),
        # A panel without maps is left out
        (['att', 'att_basil'], ['Transit']),
    ],
)
def test_multipld_by_tissue_type_plot(tmp_path, fields, panels):
    """The multi-delay figure has a panel per parameter with maps, and a row per map."""
    seg_file, files = _write_maps(tmp_path)
    results = MultiPLDByTissueTypePlot(
        seg_file=seg_file,
        **dict.fromkeys(fields, files[0]),
    ).run(cwd=tmp_path)
    text = open(results.outputs.out_file).read()
    for panel in ['Transit', 'Arrival', 'Volume']:
        assert (panel in text) == (panel in panels)

    assert ('BASIL' in text) == any(field.endswith('_basil') for field in fields)


def test_multipld_by_tissue_type_plot_no_maps(tmp_path):
    """The multi-delay figure needs at least one map."""
    seg_file, _ = _write_maps(tmp_path)
    with pytest.raises(ValueError, match='None of the multi-delay'):
        MultiPLDByTissueTypePlot(seg_file=seg_file).run(cwd=tmp_path)
