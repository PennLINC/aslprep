"""Tests for the aslprep.interfaces.cbf submodule."""

import os

import nibabel as nb
import numpy as np
import pandas as pd
import pytest

from aslprep.interfaces import cbf


def test_computecbf_casl(datasets, tmp_path_factory):
    """Test aslprep.interfaces.cbf.ComputeCBF with (P)CASL."""
    tmpdir = tmp_path_factory.mktemp("test_computecbf_casl")
    aslcontext_file = os.path.join(datasets["test_001"], "sub-01/perf/sub-01_aslcontext.tsv")

    aslcontext = pd.read_table(aslcontext_file)
    n_deltam = aslcontext.loc[aslcontext["volume_type"] == "label"].shape[0]
    n_volumes = aslcontext.shape[0]

    # Simulate ASL data and a brain mask.
    asl_data = np.random.random((30, 30, 30, n_deltam)).astype(np.float32)
    asl_file = _save_img(asl_data, tmpdir, "asl.nii.gz")
    asl_mask = np.zeros((30, 30, 30), dtype=np.uint8)
    asl_mask[10:20, 10:20, 10:20] = 1
    mask_file = _save_img(asl_mask, tmpdir, "mask.nii.gz")
    m0_file = _save_img(asl_mask, tmpdir, "m0.nii.gz")

    single_pld = 1.5
    plds = np.zeros(n_volumes)
    temp_plds = np.linspace(0.5, 3.5, n_deltam)
    plds[aslcontext["volume_type"] == "m0scan"] = 0
    plds[aslcontext["volume_type"] == "label"] = temp_plds
    plds[aslcontext["volume_type"] == "control"] = temp_plds
    bad_multiple_plds = plds.tolist()
    good_multiple_plds = plds[aslcontext["volume_type"] == "control"]

    BASE_METADATA = {
        "MagneticFieldStrength": 3,
        "LabelingDuration": 1.6,
    }
    ACQ_DICTS = [
        {"MRAcquisitionType": "3D"},
        {
            "MRAcquisitionType": "2D",
            "SliceTiming": list(np.linspace(0.1, 0.5, 30)),
        },
    ]

    for asltype in ["PCASL", "CASL"]:
        for acq_dict in ACQ_DICTS:
            # Scenario 1: PCASL with a single PostLabelingDelay
            # This should produce CBF time series and mean CBF, but no ATT
            metadata = {
                "ArterialSpinLabelingType": asltype,
                "PostLabelingDelay": single_pld,
                **BASE_METADATA,
                **acq_dict,
            }
            interface = cbf.ComputeCBF(
                cbf_only=False,
                deltam=asl_file,
                metadata=metadata,
                m0_scale=1,
                m0_file=m0_file,
                mask=mask_file,
            )
            results = interface.run(cwd=tmpdir)
            assert os.path.isfile(results.outputs.cbf_ts)
            cbf_img = nb.load(results.outputs.cbf_ts)
            assert cbf_img.ndim == 4
            assert cbf_img.shape[3] == n_deltam
            assert os.path.isfile(results.outputs.mean_cbf)
            mean_cbf_img = nb.load(results.outputs.mean_cbf)
            assert mean_cbf_img.ndim == 3
            assert results.outputs.att is None

            # Scenario 2: PCASL with one PostLabelingDelay for each volume (bad)
            metadata = {
                "ArterialSpinLabelingType": asltype,
                "PostLabelingDelay": bad_multiple_plds,
                **BASE_METADATA,
                **acq_dict,
            }

            interface = cbf.ComputeCBF(
                cbf_only=False,
                deltam=asl_file,
                metadata=metadata,
                m0_scale=1,
                m0_file=m0_file,
                mask=mask_file,
            )
            with pytest.raises(ValueError, match="Number of PostLabelingDelays"):
                results = interface.run(cwd=tmpdir)

            # Scenario 3: PCASL with one PostLabelingDelay for each deltam volume (good)
            # This should produce ATT and mean CBF volumes, but no CBF time series
            metadata = {
                "ArterialSpinLabelingType": asltype,
                "PostLabelingDelay": good_multiple_plds,
                **BASE_METADATA,
                **acq_dict,
            }

            interface = cbf.ComputeCBF(
                cbf_only=False,
                deltam=asl_file,
                metadata=metadata,
                m0_scale=1,
                m0_file=m0_file,
                mask=mask_file,
            )
            results = interface.run(cwd=tmpdir)
            assert results.outputs.cbf_ts is None
            assert os.path.isfile(results.outputs.mean_cbf)
            mean_cbf_img = nb.load(results.outputs.mean_cbf)
            assert mean_cbf_img.ndim == 3
            assert os.path.isfile(results.outputs.att)
            att_img = nb.load(results.outputs.att)
            assert att_img.ndim == 3


def test_computecbf_pasl(datasets, tmp_path_factory):
    """Test aslprep.interfaces.cbf.ComputeCBF with PASL."""
    tmpdir = tmp_path_factory.mktemp("test_computecbf_pasl")
    aslcontext_file = os.path.join(datasets["test_001"], "sub-01/perf/sub-01_aslcontext.tsv")

    aslcontext = pd.read_table(aslcontext_file)
    n_deltam = aslcontext.loc[aslcontext["volume_type"] == "label"].shape[0]
    n_volumes = aslcontext.shape[0]

    # Simulate ASL data and a brain mask.
    asl_data = np.random.random((30, 30, 30, n_deltam)).astype(np.float32)
    asl_file = _save_img(asl_data, tmpdir, "asl.nii.gz")
    asl_mask = np.zeros((30, 30, 30), dtype=np.uint8)
    asl_mask[10:20, 10:20, 10:20] = 1
    mask_file = _save_img(asl_mask, tmpdir, "mask.nii.gz")
    m0_file = _save_img(asl_mask, tmpdir, "m0.nii.gz")

    single_pld = 1.5
    plds = np.zeros(n_volumes)
    temp_plds = np.linspace(0.5, 3.5, n_deltam)
    plds[aslcontext["volume_type"] == "m0scan"] = 0
    plds[aslcontext["volume_type"] == "label"] = temp_plds
    plds[aslcontext["volume_type"] == "control"] = temp_plds
    bad_multiple_plds = plds.tolist()
    good_multiple_plds = plds[aslcontext["volume_type"] == "control"]

    BASE_METADATA = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
    }
    ACQ_DICTS = [
        {"MRAcquisitionType": "3D"},
        {
            "MRAcquisitionType": "2D",
            "SliceTiming": list(np.linspace(0.1, 0.5, 30)),
        },
    ]

    for acq_dict in ACQ_DICTS:
        # Scenario 1: PASL without BolusCutOff (raises ValueError).
        metadata = {
            "BolusCutOffFlag": False,
            "PostLabelingDelay": single_pld,
            **BASE_METADATA,
            **acq_dict,
        }
        with pytest.raises(ValueError, match="not supported in ASLPrep."):
            interface = cbf.ComputeCBF(
                cbf_only=False,
                deltam=asl_file,
                metadata=metadata,
                m0_scale=1,
                m0_file=m0_file,
                mask=mask_file,
            )
            results = interface.run(cwd=tmpdir)

        # Scenario 2: QUIPSS PASL with a single PostLabelingDelay
        # This should produce CBF time series and mean CBF, but no ATT
        metadata = {
            "BolusCutOffFlag": True,
            "BolusCutOffTechnique": "QUIPSS",
            "BolusCutOffDelayTime": 0.5,
            "PostLabelingDelay": single_pld,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
            mask=mask_file,
        )
        results = interface.run(cwd=tmpdir)
        assert os.path.isfile(results.outputs.cbf_ts)
        cbf_img = nb.load(results.outputs.cbf_ts)
        assert cbf_img.ndim == 4
        assert cbf_img.shape[3] == n_deltam
        assert os.path.isfile(results.outputs.mean_cbf)
        mean_cbf_img = nb.load(results.outputs.mean_cbf)
        assert mean_cbf_img.ndim == 3

        # Scenario 3: QUIPSS PASL with one PostLabelingDelay for each volume (bad)
        metadata = {
            "BolusCutOffFlag": True,
            "BolusCutOffTechnique": "QUIPSS",
            "BolusCutOffDelayTime": 0.5,
            "PostLabelingDelay": bad_multiple_plds,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
            mask=mask_file,
        )
        with pytest.raises(ValueError, match="Multi-delay data are not supported"):
            results = interface.run(cwd=tmpdir)

        # Scenario 4: QUIPSS PASL with one PostLabelingDelay for each deltam volume (good)
        metadata = {
            "BolusCutOffFlag": True,
            "BolusCutOffTechnique": "QUIPSS",
            "BolusCutOffDelayTime": 0.5,
            "PostLabelingDelay": good_multiple_plds,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
            mask=mask_file,
        )
        with pytest.raises(ValueError, match="Multi-delay data are not supported"):
            results = interface.run(cwd=tmpdir)

        # Scenario 5: QUIPSSII PASL with one PostLabelingDelay
        # This should produce CBF time series and mean CBF, but no ATT
        metadata = {
            "BolusCutOffFlag": True,
            "BolusCutOffTechnique": "QUIPSSII",
            "BolusCutOffDelayTime": 0.5,
            "PostLabelingDelay": single_pld,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
            mask=mask_file,
        )
        results = interface.run(cwd=tmpdir)
        assert os.path.isfile(results.outputs.cbf_ts)
        cbf_img = nb.load(results.outputs.cbf_ts)
        assert cbf_img.ndim == 4
        assert cbf_img.shape[3] == n_deltam
        assert os.path.isfile(results.outputs.mean_cbf)
        mean_cbf_img = nb.load(results.outputs.mean_cbf)
        assert mean_cbf_img.ndim == 3

        # Scenario 6: QUIPSSII PASL with multiple PostLabelingDelays
        metadata = {
            "BolusCutOffFlag": True,
            "BolusCutOffTechnique": "QUIPSSII",
            "BolusCutOffDelayTime": 0.5,
            "PostLabelingDelay": good_multiple_plds,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
            mask=mask_file,
        )
        with pytest.raises(ValueError, match="Multi-delay data are not supported"):
            results = interface.run(cwd=tmpdir)

        # Scenario 7: Q2TIPS PASL with one PostLabelingDelay
        # This should produce CBF time series and mean CBF, but no ATT
        metadata = {
            "BolusCutOffFlag": True,
            "BolusCutOffTechnique": "Q2TIPS",
            "BolusCutOffDelayTime": [0.7, 1.6],
            "PostLabelingDelay": single_pld,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
            mask=mask_file,
        )
        results = interface.run(cwd=tmpdir)
        assert os.path.isfile(results.outputs.cbf_ts)
        cbf_img = nb.load(results.outputs.cbf_ts)
        assert cbf_img.ndim == 4
        assert cbf_img.shape[3] == n_deltam
        assert os.path.isfile(results.outputs.mean_cbf)
        mean_cbf_img = nb.load(results.outputs.mean_cbf)
        assert mean_cbf_img.ndim == 3

        # Scenario 8: Q2TIPS PASL with multiple PostLabelingDelays
        metadata = {
            "BolusCutOffFlag": True,
            "BolusCutOffTechnique": "Q2TIPS",
            "BolusCutOffDelayTime": [0.7, 1.6],
            "PostLabelingDelay": good_multiple_plds,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
            mask=mask_file,
        )
        with pytest.raises(ValueError, match="Multi-delay data are not supported"):
            results = interface.run(cwd=tmpdir)


def test_compare_slicetiming(datasets, tmp_path_factory):
    """CBF estimates should be the same with slice-timing modification as not.

    As long as the slice times are all zero, of course.
    """
    tmpdir = tmp_path_factory.mktemp("test_computecbf_casl")
    aslcontext_file = os.path.join(datasets["test_001"], "sub-01/perf/sub-01_aslcontext.tsv")

    aslcontext = pd.read_table(aslcontext_file)
    n_deltam = aslcontext.loc[aslcontext["volume_type"] == "label"].shape[0]

    # Simulate ASL data and a brain mask.
    asl_data = np.random.random((30, 30, 30, n_deltam)).astype(np.float32)
    asl_file = _save_img(asl_data, tmpdir, "asl.nii.gz")
    asl_mask = np.zeros((30, 30, 30), dtype=np.uint8)
    asl_mask[10:20, 10:20, 10:20] = 1
    mask_file = _save_img(asl_mask, tmpdir, "mask.nii.gz")
    m0_file = _save_img(asl_mask, tmpdir, "m0.nii.gz")

    ACQ_DICTS = [
        {"MRAcquisitionType": "3D"},
        {
            "MRAcquisitionType": "2D",
            "SliceTiming": list(np.zeros(30)),
        },
    ]

    cbf_data = []
    for acq_dict in ACQ_DICTS:
        metadata = {
            "ArterialSpinLabelingType": "PCASL",
            "MagneticFieldStrength": 3,
            "LabelingDuration": 1.6,
            "PostLabelingDelay": 1.5,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
            mask=mask_file,
        )
        results = interface.run(cwd=tmpdir)
        cbf_data.append(nb.load(results.outputs.cbf_ts).get_fdata())

    assert np.array_equal(cbf_data[0], cbf_data[1])


def _save_img(arr, tmpdir, fname):
    img = nb.Nifti1Image(arr, affine=np.eye(4))
    fname = os.path.join(tmpdir, fname)
    img.to_filename(fname)
    return fname
