"""Test CBF."""
import json
import os

import nibabel as nb
import numpy as np
import pandas as pd
import pytest

from aslprep.utils import misc


def test_computecbf_casl(datasets, tmp_path_factory):
    """Test aslprep.interfaces.cbf_computation.ComputeCBF with (P)CASL."""
    tmpdir = tmp_path_factory.mktemp("test_computecbf_casl")
    json_file = os.path.join(datasets["dset"], "sub-01/perf/sub-01_asl.json")
    aslcontext_file = os.path.join(datasets["dset"], "sub-01/perf/sub-01_aslcontext.tsv")

    aslcontext = pd.read_table(aslcontext_file)
    aslcontext = aslcontext.loc[aslcontext["volume_type"] != "m0scan"]
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

    plds = np.zeros(n_volumes)
    temp_plds = np.linspace(0, 1, n_deltam)
    plds[aslcontext["volume_type"] == "label"] = temp_plds
    plds[aslcontext["volume_type"] == "control"] = temp_plds
    bad_multiple_plds = plds.tolist()

    # Scenario 1: PCASL with a single PostLabelingDelay
    metadata = {
        "ArterialSpinLabelingType": "PCASL",
        "MagneticFieldStrength": 3,
        "LabelingDuration": 1.6,
        "PostLabelingDelay": bad_multiple_plds,
    }
    tcbf, meancbf, att = misc.compute_cbf(
        metadata=metadata,
        mask=mask_file,
        m0file=m0_file,
        cbffile=asl_file,
        m0scale=1,
    )
    raise Exception(tcbf.shape)


def _save_img(arr, tmpdir, fname):
    img = nb.Nifti1Image(arr, affine=np.eye(4))
    fname = os.path.join(tmpdir, fname)
    img.to_filename(fname)
    return fname
