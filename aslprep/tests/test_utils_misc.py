"""Tests for the aslprep.utils.misc module."""
import numpy as np

from aslprep.utils import misc


def test_estimate_att_pcasl_3d():
    """Smoke test aslprep.utils.misc.estimate_att_pcasl with 3D data (no slice timing)."""
    n_voxels, n_unique_plds = 1000, 6
    base_plds = np.linspace(0.25, 1.5, n_unique_plds)
    plds = np.repeat(base_plds, n_voxels, axis=0)
    deltam_arr = np.random.random((n_voxels, n_unique_plds))
    tau = np.full(n_unique_plds, 1.4)
    att = misc.estimate_att_pcasl(
        deltam_arr=deltam_arr,
        plds=plds,
        lds=tau,
        t1blood=1.6,
        t1tissue=1.3,
    )
    assert att.shape == (n_voxels,)


def test_estimate_att_pcasl_2d():
    """Smoke test aslprep.utils.misc.estimate_att_pcasl with slice-shifted PLDs."""
    n_voxels, n_unique_plds, n_slice_times = 1000, 6, 50
    base_plds = np.linspace(0.25, 1.5, n_unique_plds)
    slice_times = np.linspace(0, 3, n_slice_times)
    slice_shifted_plds = np.dot(slice_times[:, None], base_plds[None, :])
    plds = np.repeat(slice_shifted_plds, n_voxels // n_slice_times, axis=0)
    deltam_arr = np.random.random((n_voxels, n_unique_plds))
    tau = np.full(n_unique_plds, 1.4)
    att = misc.estimate_att_pcasl(
        deltam_arr=deltam_arr,
        plds=plds,
        lds=tau,
        t1blood=1.6,
        t1tissue=1.3,
    )
    assert att.shape == (n_voxels,)


def test_estimate_cbf_pcasl_multipld():
    """Smoke test aslprep.utils.misc.estimate_cbf_pcasl_multipld with slice-shifted PLDs."""
    n_voxels, n_unique_plds, n_slice_times = 1000, 6, 50
    base_plds = np.linspace(0.25, 1.5, n_unique_plds)
    slice_times = np.linspace(0, 3, n_slice_times)
    slice_shifted_plds = np.dot(slice_times[:, None], base_plds[None, :])
    plds = np.repeat(slice_shifted_plds, n_voxels // n_slice_times, axis=0)
    deltam_arr = np.random.random((n_voxels, n_unique_plds))
    scaled_m0data = np.random.random(n_voxels)

    misc.estimate_cbf_pcasl_multipld(
        deltam_arr=deltam_arr,
        scaled_m0data=scaled_m0data,
        plds=plds,
        tau=np.array(1.4),
        labeleff=0.7,
        t1blood=1.6,
        t1tissue=1.3,
        unit_conversion=6000,
        partition_coefficient=0.9,
    )
