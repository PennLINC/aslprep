"""Tests for aslprep.interfaces.parcellation module."""

import os

import nibabel as nb
import numpy as np
import pandas as pd
import pytest

from aslprep.interfaces.parcellation import NiftiParcellate


def _read_tables(results):
    scalarstats = pd.read_table(results.outputs.scalarstats, index_col='node')
    coverage = pd.read_table(results.outputs.coverage, index_col='node')
    return scalarstats, coverage


def test_nifti_parcellate(tmp_path_factory):
    """Parcels missing from the atlas or with too little coverage are n/a."""
    tmpdir = tmp_path_factory.mktemp('test_nifti_parcellate')

    arr = np.zeros((4, 4, 4))
    arr[0, 0, :] = 3
    arr[0, 1, :] = 4
    arr[0, 2, :] = 5
    atlas_img = nb.Nifti1Image(arr.astype(np.int32), np.eye(4))
    lut = pd.DataFrame(
        columns=['index', 'name'],
        data=[[1, 'Region A'], [2, 'Region B'], [3, 'Region C'], [4, 'Region D'], [5, 'Region E']],
    )

    mask = np.ones((4, 4, 4))
    mask_img = nb.Nifti1Image(mask.astype(np.int32), np.eye(4))
    atlas_file = os.path.join(tmpdir, 'atlas_01.nii.gz')
    atlas_img.to_filename(atlas_file)
    mask_file = os.path.join(tmpdir, 'mask_01.nii.gz')
    mask_img.to_filename(mask_file)
    lut_file = os.path.join(tmpdir, 'lut_01.tsv')
    lut.to_csv(lut_file, sep='\t', index=False)

    # Some parcels are not present, but none are masked out
    results = NiftiParcellate(
        in_files=[atlas_file],
        masks=[mask_file],
        names=['cbf'],
        atlas=atlas_file,
        atlas_labels=lut_file,
        min_coverage=0.5,
    ).run(cwd=tmpdir)
    scalarstats, coverage = _read_tables(results)
    assert scalarstats.index.tolist() == lut['name'].tolist()
    assert scalarstats.columns.tolist() == coverage.columns.tolist() == ['cbf']
    np.testing.assert_array_equal(coverage['cbf'], [0, 0, 1, 1, 1])
    np.testing.assert_array_equal(scalarstats['cbf'], [np.nan, np.nan, 3, 4, 5])

    # Now let's mask out some voxels
    mask[0, 0, 0] = 0  # 1/4 of the third parcel
    mask[0, 1, :2] = 0  # 1/2 of the fourth parcel
    mask[0, 2, :3] = 0  # 3/4 of the fifth parcel
    mask_img = nb.Nifti1Image(mask.astype(np.int32), np.eye(4))
    mask_file = os.path.join(tmpdir, 'mask_02.nii.gz')
    mask_img.to_filename(mask_file)
    results = NiftiParcellate(
        in_files=[atlas_file],
        masks=[mask_file],
        names=['cbf'],
        atlas=atlas_file,
        atlas_labels=lut_file,
        min_coverage=0.5,
    ).run(cwd=tmpdir)
    scalarstats, coverage = _read_tables(results)
    np.testing.assert_array_equal(coverage['cbf'], [0, 0, 0.75, 0.5, 0.25])
    np.testing.assert_array_equal(scalarstats['cbf'], [np.nan, np.nan, 3, 4, np.nan])


def _write_img(path, data):
    nb.Nifti1Image(data, np.eye(4)).to_filename(path)
    return str(path)


def _atlas_inputs(tmp_path):
    atlas = np.zeros((4, 4, 4), np.int32)
    atlas[:2] = 1
    atlas[2:] = 2
    lut_file = str(tmp_path / 'lut.tsv')
    pd.DataFrame({'index': [1, 2], 'name': ['A', 'B']}).to_csv(lut_file, sep='\t', index=False)
    return {'atlas': _write_img(tmp_path / 'atlas.nii.gz', atlas), 'atlas_labels': lut_file}


def test_nifti_parcellate_validity_mask(tmp_path):
    """With a validity mask, failed estimates (stored as 0) are left out; genuine zeros are not."""
    values = np.zeros((4, 4, 4))
    values[:2, :2] = 60.0  # parcel A: half valid 60, half failed (0)
    values[2:] = 0.0  # parcel B: genuine zero estimates
    mask = np.zeros((4, 4, 4), np.int32)
    mask[:2, :2] = 1
    mask[2:] = 1
    results = NiftiParcellate(
        in_files=[_write_img(tmp_path / 'map.nii.gz', values)],
        masks=[_write_img(tmp_path / 'mask.nii.gz', mask)],
        names=['cbf_lrGM'],
        min_coverage=0.25,
        **_atlas_inputs(tmp_path),
    ).run(cwd=tmp_path)
    scalarstats, coverage = _read_tables(results)
    assert scalarstats.loc['A', 'cbf_lrGM'] == 60.0
    assert scalarstats.loc['B', 'cbf_lrGM'] == 0.0
    np.testing.assert_allclose(coverage['cbf_lrGM'], [0.5, 1.0])


def test_nifti_parcellate_empty_mask(tmp_path):
    """An empty mask (nothing estimated) gives zero coverage and no values, not an error."""
    results = NiftiParcellate(
        in_files=[_write_img(tmp_path / 'map.nii.gz', np.zeros((4, 4, 4)))],
        masks=[_write_img(tmp_path / 'mask.nii.gz', np.zeros((4, 4, 4), np.int32))],
        names=['cbf_lrGM'],
        min_coverage=0.5,
        **_atlas_inputs(tmp_path),
    ).run(cwd=tmp_path)
    scalarstats, coverage = _read_tables(results)
    assert scalarstats.shape == (2, 1)
    assert scalarstats['cbf_lrGM'].isna().all()
    np.testing.assert_array_equal(coverage['cbf_lrGM'], [0.0, 0.0])


def test_nifti_parcellate_several_maps(tmp_path):
    """Maps sharing a mask are parcellated together; each keeps its own mask and column."""
    full_mask = np.ones((4, 4, 4), np.int32)
    half_mask = np.zeros((4, 4, 4), np.int32)
    half_mask[:2, :2] = 1  # half of parcel A, none of parcel B
    cbf = np.full((4, 4, 4), 50.0)
    cbf[2:] = 40.0
    att = np.full((4, 4, 4), 1.5)
    sparse = np.zeros((4, 4, 4))
    sparse[:2, :2] = 70.0

    full_mask_file = _write_img(tmp_path / 'full_mask.nii.gz', full_mask)
    results = NiftiParcellate(
        in_files=[
            _write_img(tmp_path / 'cbf.nii.gz', cbf),
            _write_img(tmp_path / 'sparse.nii.gz', sparse),
            _write_img(tmp_path / 'att.nii.gz', att),
        ],
        masks=[
            full_mask_file,
            _write_img(tmp_path / 'half_mask.nii.gz', half_mask),
            full_mask_file,
        ],
        names=['cbf', 'cbf_lrGM', 'att'],
        min_coverage=0.5,
        **_atlas_inputs(tmp_path),
    ).run(cwd=tmp_path)
    scalarstats, coverage = _read_tables(results)

    # Columns keep the order of the inputs, not the order of the mask groups.
    assert scalarstats.columns.tolist() == coverage.columns.tolist() == ['cbf', 'cbf_lrGM', 'att']
    np.testing.assert_allclose(scalarstats['cbf'], [50.0, 40.0])
    np.testing.assert_allclose(scalarstats['att'], [1.5, 1.5])
    np.testing.assert_allclose(scalarstats.loc['A', 'cbf_lrGM'], 70.0)
    assert np.isnan(scalarstats.loc['B', 'cbf_lrGM'])
    np.testing.assert_allclose(coverage['cbf'], [1.0, 1.0])
    np.testing.assert_allclose(coverage['att'], [1.0, 1.0])
    np.testing.assert_allclose(coverage['cbf_lrGM'], [0.5, 0.0])


def test_nifti_parcellate_rejects_bad_inputs(tmp_path):
    """Mismatched inputs or repeated names are errors, not silently misaligned columns."""
    map_file = _write_img(tmp_path / 'map.nii.gz', np.ones((4, 4, 4)))
    mask_file = _write_img(tmp_path / 'mask.nii.gz', np.ones((4, 4, 4), np.int32))
    atlas_inputs = _atlas_inputs(tmp_path)

    with pytest.raises(ValueError, match='2 maps, 1 masks and 2 names'):
        NiftiParcellate(
            in_files=[map_file, map_file], masks=[mask_file], names=['a', 'b'], **atlas_inputs
        ).run(cwd=tmp_path)
    with pytest.raises(ValueError, match='unique'):
        NiftiParcellate(
            in_files=[map_file, map_file],
            masks=[mask_file, mask_file],
            names=['a', 'a'],
            **atlas_inputs,
        ).run(cwd=tmp_path)
