"""Tests for the aslprep.interfaces.cbf_computation submodule."""
import os

import nibabel as nb
import numpy as np
import pandas as pd
import pytest

from aslprep.interfaces import cbf_computation


def test_computecbf_casl(datasets, tmp_path_factory):
    """Test aslprep.interfaces.cbf_computation.ComputeCBF with (P)CASL."""
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
    temp_plds = np.linspace(0, 1, n_deltam)
    plds[aslcontext["volume_type"] == "m0scan"] = 0
    plds[aslcontext["volume_type"] == "label"] = temp_plds
    plds[aslcontext["volume_type"] == "control"] = temp_plds
    bad_multiple_plds = plds.tolist()
    good_multiple_plds = plds[aslcontext["volume_type"] == "control"]

    # Scenario 1: PCASL with a single PostLabelingDelay
    metadata = {
        "ArterialSpinLabelingType": "PCASL",
        "MagneticFieldStrength": 3,
        "LabelingDuration": 1.6,
        "PostLabelingDelay": single_pld,
    }
    interface = cbf_computation.ComputeCBF(
        cbf_only=False,
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
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

    # Scenario 2: PCASL with one PostLabelingDelay for each volume (bad)
    metadata = {
        "ArterialSpinLabelingType": "PCASL",
        "MagneticFieldStrength": 3,
        "LabelingDuration": 1.6,
        "PostLabelingDelay": bad_multiple_plds,
    }

    interface = cbf_computation.ComputeCBF(
        cbf_only=False,
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    with pytest.raises(ValueError, match="Number of PostLabelingDelay values"):
        interface.run(cwd=tmpdir)

    # Scenario 3: PCASL with one PostLabelingDelay for each deltam volume (good)
    metadata = {
        "ArterialSpinLabelingType": "PCASL",
        "MagneticFieldStrength": 3,
        "LabelingDuration": 1.6,
        "PostLabelingDelay": good_multiple_plds,
    }

    interface = cbf_computation.ComputeCBF(
        cbf_only=False,
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
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

    # Scenario 4: CASL with a single PostLabelingDelay
    metadata = {
        "ArterialSpinLabelingType": "CASL",
        "MagneticFieldStrength": 3,
        "LabelingDuration": 1.6,
        "PostLabelingDelay": single_pld,
    }
    interface = cbf_computation.ComputeCBF(
        cbf_only=False,
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
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

    # Scenario 5: CASL with multiple PostLabelingDelays
    metadata = {
        "ArterialSpinLabelingType": "CASL",
        "MagneticFieldStrength": 3,
        "LabelingDuration": 1.6,
        "PostLabelingDelay": good_multiple_plds,
    }
    interface = cbf_computation.ComputeCBF(
        cbf_only=False,
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
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
    """Test aslprep.interfaces.cbf_computation.ComputeCBF with PASL."""
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
    temp_plds = np.linspace(0, 1, n_deltam)
    plds[aslcontext["volume_type"] == "m0scan"] = 0
    plds[aslcontext["volume_type"] == "label"] = temp_plds
    plds[aslcontext["volume_type"] == "control"] = temp_plds
    bad_multiple_plds = plds.tolist()
    good_multiple_plds = plds[aslcontext["volume_type"] == "control"]

    # Scenario 1: PASL without BolusCutOff (raises ValueError).
    metadata = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": False,
        "PostLabelingDelay": single_pld,
    }
    with pytest.raises(ValueError, match="not supported in ASLPrep."):
        interface = cbf_computation.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0scale=1,
            m0_file=m0_file,
            mask=mask_file,
        )
        results = interface.run(cwd=tmpdir)

    # Scenario 2: QUIPSS PASL with a single PostLabelingDelay
    metadata = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": True,
        "BolusCutOffTechnique": "QUIPSS",
        "BolusCutOffDelayTime": 0.5,
        "PostLabelingDelay": single_pld,
    }
    interface = cbf_computation.ComputeCBF(
        cbf_only=False,
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
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
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": True,
        "BolusCutOffTechnique": "QUIPSS",
        "BolusCutOffDelayTime": 0.5,
        "PostLabelingDelay": bad_multiple_plds,
    }
    interface = cbf_computation.ComputeCBF(
        cbf_only=False,
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    with pytest.raises(ValueError, match="Multi-PLD data are not supported for PASL"):
        results = interface.run(cwd=tmpdir)

    # Scenario 4: QUIPSS PASL with one PostLabelingDelay for each deltam volume (good)
    metadata = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": True,
        "BolusCutOffTechnique": "QUIPSS",
        "BolusCutOffDelayTime": 0.5,
        "PostLabelingDelay": good_multiple_plds,
    }
    interface = cbf_computation.ComputeCBF(
        cbf_only=False,
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    with pytest.raises(ValueError, match="Multi-PLD data are not supported for PASL"):
        results = interface.run(cwd=tmpdir)

    # Scenario 5: QUIPSSII PASL with one PostLabelingDelay
    metadata = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": True,
        "BolusCutOffTechnique": "QUIPSSII",
        "BolusCutOffDelayTime": 0.5,
        "PostLabelingDelay": single_pld,
    }
    interface = cbf_computation.ComputeCBF(
        cbf_only=False,
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
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
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": True,
        "BolusCutOffTechnique": "QUIPSSII",
        "BolusCutOffDelayTime": 0.5,
        "PostLabelingDelay": good_multiple_plds,
    }
    interface = cbf_computation.ComputeCBF(
        cbf_only=False,
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    with pytest.raises(ValueError, match="Multi-PLD data are not supported for PASL"):
        results = interface.run(cwd=tmpdir)

    # Scenario 7: Q2TIPS PASL with one PostLabelingDelay
    metadata = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": True,
        "BolusCutOffTechnique": "Q2TIPS",
        "BolusCutOffDelayTime": [0.7, 1.6],
        "PostLabelingDelay": single_pld,
    }
    interface = cbf_computation.ComputeCBF(
        cbf_only=False,
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    with pytest.raises(ValueError, match="Q2TIPS is not supported in ASLPrep."):
        results = interface.run(cwd=tmpdir)

    # Scenario 8: Q2TIPS PASL with multiple PostLabelingDelays
    metadata = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": True,
        "BolusCutOffTechnique": "Q2TIPS",
        "BolusCutOffDelayTime": [0.7, 1.6],
        "PostLabelingDelay": good_multiple_plds,
    }
    interface = cbf_computation.ComputeCBF(
        cbf_only=False,
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    with pytest.raises(ValueError, match="Multi-PLD data are not supported for PASL"):
        results = interface.run(cwd=tmpdir)


def _save_img(arr, tmpdir, fname):
    img = nb.Nifti1Image(arr, affine=np.eye(4))
    fname = os.path.join(tmpdir, fname)
    img.to_filename(fname)
    return fname
