"""Tests for the aslprep.interfaces.utility submodule."""

import nibabel as nb
import numpy as np
import pandas as pd

from aslprep.interfaces.utility import ReduceASLFiles
from aslprep.utils.asl import ignore_volume_types, select_processing_target


def test_reduceaslfiles_drops_ignored_volumes(tmp_path):
    """Volumes ignored with --ignore are removed, along with their volume-wise metadata."""
    volume_types = ['m0scan', 'deltam', 'cbf', 'deltam']
    data = np.zeros((2, 2, 2, len(volume_types)), np.float32)
    data += np.arange(len(volume_types))
    asl_file = tmp_path / 'sub-01_asl.nii.gz'
    nb.Nifti1Image(data, np.eye(4)).to_filename(asl_file)
    aslcontext = tmp_path / 'sub-01_aslcontext.tsv'
    aslcontext.write_text('\n'.join(['volume_type', *volume_types]) + '\n')

    ignored_aslcontext = ignore_volume_types(
        str(aslcontext),
        ignore=['deltam'],
        out_dir=tmp_path / 'work',
    )
    processing_target = select_processing_target(ignored_aslcontext)
    assert processing_target == 'cbf'

    results = ReduceASLFiles(
        asl_file=str(asl_file),
        aslcontext=ignored_aslcontext,
        processing_target=processing_target,
        metadata={'PostLabelingDelay': [0, 1.5, 1.8, 2.0], 'M0Type': 'Included'},
    ).run(cwd=tmp_path)

    assert pd.read_table(results.outputs.aslcontext)['volume_type'].tolist() == ['m0scan', 'cbf']
    reduced_data = nb.load(results.outputs.asl_file).get_fdata()
    assert reduced_data.shape[3] == 2
    assert np.array_equal(reduced_data[0, 0, 0, :], [0, 2])
    assert results.outputs.metadata['PostLabelingDelay'] == [0, 1.8]
