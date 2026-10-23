"""Tests for how ASLPrep's derivatives are written.

Every derivative niworkflows rewrites is compressed with ISA-L instead of zlib.
"""

import gzip

import nibabel as nb
import numpy as np
import pytest
from niworkflows.interfaces import bids as niworkflows_bids
from niworkflows.utils.images import unsafe_write_nifti_header_and_data as zlib_writer

from aslprep.interfaces import bids

needs_isal = pytest.mark.skipif(bids.IGzipFile is None, reason='isal is not installed')


def _image(shape=(6, 5, 4, 3), dtype='float32'):
    rng = np.random.default_rng(0)
    data = (rng.standard_normal(shape) * 100 + 1000).astype(dtype)
    return nb.Nifti1Image(data, np.diag([1.8, 1.8, 1.8, 1.0]))


@needs_isal
def test_niworkflows_sinks_write_with_isal():
    """Test that importing ASLPrep's sinks routes niworkflows' writer to ISA-L."""
    assert niworkflows_bids.unsafe_write_nifti_header_and_data is bids.write_nifti_isal


@needs_isal
@pytest.mark.parametrize('dtype', ['float32', 'float64', 'int16'])
def test_isal_writer_matches_the_zlib_writer(tmp_path, dtype):
    """Test that the ISA-L file decompresses to exactly the bytes niworkflows writes."""
    img = _image(dtype=dtype)
    data = np.asanyarray(img.dataobj)
    zlib_writer(tmp_path / 'zlib.nii.gz', img.header, data)
    bids.write_nifti_isal(tmp_path / 'isal.nii.gz', img.header, data)

    with gzip.open(tmp_path / 'zlib.nii.gz') as a, gzip.open(tmp_path / 'isal.nii.gz') as b:
        assert a.read() == b.read()
    out = nb.load(tmp_path / 'isal.nii.gz')
    assert out.get_data_dtype() == np.dtype(dtype)
    assert np.array_equal(np.asanyarray(out.dataobj), data)


@needs_isal
def test_isal_writer_is_deterministic(tmp_path):
    """Test that two writes of one image are byte-identical (no name or time in the header)."""
    img = _image()
    data = np.asanyarray(img.dataobj)
    bids.write_nifti_isal(tmp_path / 'a.nii.gz', img.header, data)
    bids.write_nifti_isal(tmp_path / 'b.nii.gz', img.header, data)
    assert (tmp_path / 'a.nii.gz').read_bytes() == (tmp_path / 'b.nii.gz').read_bytes()


@needs_isal
def test_isal_writer_leaves_uncompressed_files_uncompressed(tmp_path):
    """Test that a .nii file is written uncompressed."""
    img = _image()
    data = np.asanyarray(img.dataobj)
    bids.write_nifti_isal(tmp_path / 'x.nii', img.header, data)
    assert (tmp_path / 'x.nii').read_bytes()[:2] != b'\x1f\x8b'
    assert np.array_equal(np.asanyarray(nb.load(tmp_path / 'x.nii').dataobj), data)


def test_datasink_compresses_a_4d_series(tmp_path):
    """Test a datasink run end to end: an uncompressed input written as .nii.gz."""
    img = _image()
    img.to_filename(tmp_path / 'in.nii')
    result = bids.DerivativesDataSink(
        base_directory=str(tmp_path / 'out'),
        source_file='/data/sub-01/perf/sub-01_asl.nii.gz',
        desc='preproc',
        suffix='asl',
        extension='.nii.gz',
        compress=True,
        in_file=str(tmp_path / 'in.nii'),
    ).run(cwd=tmp_path)
    out = result.outputs.out_file
    out = out[0] if isinstance(out, list) else out
    assert out.endswith('sub-01_desc-preproc_asl.nii.gz')
    with open(out, 'rb') as fobj:
        assert fobj.read(2) == b'\x1f\x8b'
    assert np.array_equal(np.asanyarray(nb.load(out).dataobj), np.asanyarray(img.dataobj))
