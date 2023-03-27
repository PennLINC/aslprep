"""Tests for the aslprep.interfaces.cbf_computation submodule."""
import json
import os

import nibabel as nb
import numpy as np
import pandas as pd
import pytest

from aslprep.interfaces import cbf_computation


def test_computecbf_casl(datasets, tmp_path_factory):
    """Test aslprep.interfaces.cbf_computation.ComputeCBF with (P)CASL."""
    tmpdir = tmp_path_factory.mktemp("test_computecbf_casl")
    json_file = os.path.join(datasets["dset"], "sub-01/perf/sub-01_asl.json")
    aslcontext_file = os.path.join(datasets["dset"], "sub-01/perf/sub-01_aslcontext.tsv")

    aslcontext = pd.read_table(aslcontext_file)
    n_deltam = aslcontext.loc[aslcontext["volume_type"] == "label"].shape[0]
    n_volumes = aslcontext.shape[0]
    with open(json_file, "r") as fo:
        metadata = json.load(fo)

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
    multiple_plds = plds.tolist()

    # Scenario 1: PCASL with a single PostLabelingDelay
    metadata = {
        "ArterialSpinLabelingType": "PCASL",
        "MagneticFieldStrength": 3,
        "LabelingDuration": 1.6,
        "PostLabelingDelay": single_pld,
    }
    pcasl_singlepld = cbf_computation.ComputeCBF(
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    results = pcasl_singlepld.run(cwd=tmpdir)
    assert os.path.isfile(results.outputs.cbf)
    cbf_img = nb.load(results.outputs.cbf)
    assert cbf_img.ndim == 4
    assert cbf_img.shape[3] == n_deltam
    assert os.path.isfile(results.outputs.mean_cbf)
    mean_cbf_img = nb.load(results.outputs.mean_cbf)
    assert mean_cbf_img.ndim == 3

    # Scenario 2: PCASL with multiple PostLabelingDelays
    metadata = {
        "ArterialSpinLabelingType": "PCASL",
        "MagneticFieldStrength": 3,
        "LabelingDuration": 1.6,
        "PostLabelingDelay": multiple_plds,
    }

    pcasl_multipld = cbf_computation.ComputeCBF(
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    results = pcasl_multipld.run(cwd=tmpdir)
    assert os.path.isfile(results.outputs.cbf)
    cbf_img = nb.load(results.outputs.cbf)
    assert cbf_img.ndim == 4
    assert cbf_img.shape[3] == n_deltam
    assert os.path.isfile(results.outputs.mean_cbf)
    mean_cbf_img = nb.load(results.outputs.mean_cbf)
    assert mean_cbf_img.ndim == 3

    # Scenario 3: CASL with a single PostLabelingDelay
    metadata = {
        "ArterialSpinLabelingType": "CASL",
        "MagneticFieldStrength": 3,
        "LabelingDuration": 1.6,
        "PostLabelingDelay": single_pld,
    }
    pcasl_singlepld = cbf_computation.ComputeCBF(
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    results = pcasl_singlepld.run(cwd=tmpdir)
    assert os.path.isfile(results.outputs.cbf)
    cbf_img = nb.load(results.outputs.cbf)
    assert cbf_img.ndim == 4
    assert cbf_img.shape[3] == n_deltam
    assert os.path.isfile(results.outputs.mean_cbf)
    mean_cbf_img = nb.load(results.outputs.mean_cbf)
    assert mean_cbf_img.ndim == 3

    # Scenario 4: CASL with multiple PostLabelingDelays
    metadata = {
        "ArterialSpinLabelingType": "CASL",
        "MagneticFieldStrength": 3,
        "LabelingDuration": 1.6,
        "PostLabelingDelay": multiple_plds,
    }
    pcasl_multipld = cbf_computation.ComputeCBF(
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    results = pcasl_multipld.run(cwd=tmpdir)
    assert os.path.isfile(results.outputs.cbf)
    cbf_img = nb.load(results.outputs.cbf)
    assert cbf_img.ndim == 4
    assert cbf_img.shape[3] == n_deltam
    assert os.path.isfile(results.outputs.mean_cbf)
    mean_cbf_img = nb.load(results.outputs.mean_cbf)
    assert mean_cbf_img.ndim == 3


def test_computecbf_pasl(datasets, tmp_path_factory):
    """Test aslprep.interfaces.cbf_computation.ComputeCBF with PASL."""
    tmpdir = tmp_path_factory.mktemp("test_computecbf_pasl")
    aslcontext_file = os.path.join(datasets["dset"], "sub-01/perf/sub-01_aslcontext.tsv")

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
    multiple_plds = np.zeros(n_volumes)
    temp_plds = np.linspace(0, 1, n_deltam)
    multiple_plds[aslcontext["volume_type"] == "m0scan"] = 0
    multiple_plds[aslcontext["volume_type"] == "label"] = temp_plds
    multiple_plds[aslcontext["volume_type"] == "control"] = temp_plds
    multiple_plds = multiple_plds.tolist()

    # Scenario 1: PASL without BolusCutOff (raises ValueError).
    metadata = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": False,
        "PostLabelingDelay": single_pld,
    }
    with pytest.raises(ValueError, match="not supported in ASLPrep."):
        pasl_no_bcof = cbf_computation.ComputeCBF(
            deltam=asl_file,
            metadata=metadata,
            m0scale=1,
            m0_file=m0_file,
            mask=mask_file,
        )
        results = pasl_no_bcof.run(cwd=tmpdir)

    # Scenario 2: QUIPSS PASL with a single PostLabelingDelay
    metadata = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": True,
        "BolusCutOffTechnique": "QUIPSS",
        "BolusCutOffDelayTime": 0.5,
        "PostLabelingDelay": single_pld,
    }
    quipss_singlepld = cbf_computation.ComputeCBF(
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    results = quipss_singlepld.run(cwd=tmpdir)
    assert os.path.isfile(results.outputs.cbf)
    cbf_img = nb.load(results.outputs.cbf)
    assert cbf_img.ndim == 4
    assert cbf_img.shape[3] == n_deltam
    assert os.path.isfile(results.outputs.mean_cbf)
    mean_cbf_img = nb.load(results.outputs.mean_cbf)
    assert mean_cbf_img.ndim == 3

    # Scenario 3: QUIPSS PASL with multiple PostLabelingDelays
    metadata = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": True,
        "BolusCutOffTechnique": "QUIPSS",
        "BolusCutOffDelayTime": 0.5,
        "PostLabelingDelay": multiple_plds,
    }
    quipss_multipld = cbf_computation.ComputeCBF(
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    results = quipss_multipld.run(cwd=tmpdir)
    assert os.path.isfile(results.outputs.cbf)
    cbf_img = nb.load(results.outputs.cbf)
    assert cbf_img.ndim == 4
    assert cbf_img.shape[3] == n_deltam
    assert os.path.isfile(results.outputs.mean_cbf)
    mean_cbf_img = nb.load(results.outputs.mean_cbf)
    assert mean_cbf_img.ndim == 3

    # Scenario 4: QUIPSSII PASL with one PostLabelingDelay
    metadata = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": True,
        "BolusCutOffTechnique": "QUIPSSII",
        "BolusCutOffDelayTime": 0.5,
        "PostLabelingDelay": single_pld,
    }
    quipssii_singlepld = cbf_computation.ComputeCBF(
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    results = quipssii_singlepld.run(cwd=tmpdir)
    assert os.path.isfile(results.outputs.cbf)
    cbf_img = nb.load(results.outputs.cbf)
    assert cbf_img.ndim == 4
    assert cbf_img.shape[3] == n_deltam
    assert os.path.isfile(results.outputs.mean_cbf)
    mean_cbf_img = nb.load(results.outputs.mean_cbf)
    assert mean_cbf_img.ndim == 3

    # Scenario 5: QUIPSSII PASL with multiple PostLabelingDelays
    metadata = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": True,
        "BolusCutOffTechnique": "QUIPSSII",
        "BolusCutOffDelayTime": 0.5,
        "PostLabelingDelay": multiple_plds,
    }
    quipssii_multipld = cbf_computation.ComputeCBF(
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    results = quipssii_multipld.run(cwd=tmpdir)
    assert os.path.isfile(results.outputs.cbf)
    cbf_img = nb.load(results.outputs.cbf)
    assert cbf_img.ndim == 4
    assert cbf_img.shape[3] == n_deltam
    assert os.path.isfile(results.outputs.mean_cbf)
    mean_cbf_img = nb.load(results.outputs.mean_cbf)
    assert mean_cbf_img.ndim == 3

    # Scenario 6: Q2TIPS PASL with one PostLabelingDelay
    metadata = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": True,
        "BolusCutOffTechnique": "Q2TIPS",
        "BolusCutOffDelayTime": [0.7, 1.6],
        "PostLabelingDelay": single_pld,
    }
    q2tips_singlepld = cbf_computation.ComputeCBF(
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    results = q2tips_singlepld.run(cwd=tmpdir)
    assert os.path.isfile(results.outputs.cbf)
    cbf_img = nb.load(results.outputs.cbf)
    assert cbf_img.ndim == 4
    assert cbf_img.shape[3] == n_deltam
    assert os.path.isfile(results.outputs.mean_cbf)
    mean_cbf_img = nb.load(results.outputs.mean_cbf)
    assert mean_cbf_img.ndim == 3

    # Scenario 7: Q2TIPS PASL with multiple PostLabelingDelays
    metadata = {
        "ArterialSpinLabelingType": "PASL",
        "MagneticFieldStrength": 3,
        "BolusCutOffFlag": True,
        "BolusCutOffTechnique": "Q2TIPS",
        "BolusCutOffDelayTime": [0.7, 1.6],
        "PostLabelingDelay": multiple_plds,
    }
    q2tips_multipld = cbf_computation.ComputeCBF(
        deltam=asl_file,
        metadata=metadata,
        m0scale=1,
        m0_file=m0_file,
        mask=mask_file,
    )
    results = q2tips_multipld.run(cwd=tmpdir)
    assert os.path.isfile(results.outputs.cbf)
    cbf_img = nb.load(results.outputs.cbf)
    assert cbf_img.ndim == 4
    assert cbf_img.shape[3] == n_deltam
    assert os.path.isfile(results.outputs.mean_cbf)
    mean_cbf_img = nb.load(results.outputs.mean_cbf)
    assert mean_cbf_img.ndim == 3


def _save_img(arr, tmpdir, fname):
    img = nb.Nifti1Image(arr, affine=np.eye(4))
    fname = os.path.join(tmpdir, fname)
    img.to_filename(fname)
    return fname
