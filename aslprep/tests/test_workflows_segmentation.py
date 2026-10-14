"""Tests for the aslprep.workflows.segmentation module."""

import shutil

import nibabel as nb
import numpy as np
import pytest
from nipype.interfaces import fsl

from aslprep import config
from aslprep.interfaces.freesurfer import FixHeaderSynthStrip, SynthStripPadding, Zeropad
from aslprep.tests.tests import mock_config, reset_config


def _edges(wf, dst_name):
    return {
        dst_field: (src.name, src_field)
        for src, dst, data in wf._graph.edges(data=True)
        if dst.name == dst_name
        for src_field, dst_field in data['connect']
    }


@pytest.mark.parametrize('shape', [(64, 64, 64), (150, 59, 100), (160, 191, 1)])
def test_synthstrip_padding(tmp_path, shape):
    """The padding brings each in-plane axis up to a multiple of 64, split across both sides."""
    in_file = tmp_path / 'image.nii.gz'
    nb.Nifti1Image(np.zeros(shape, dtype=np.float32), np.eye(4)).to_filename(in_file)

    out = SynthStripPadding(input_image=str(in_file)).run(cwd=str(tmp_path)).outputs

    for nvoxels, low, high in [(shape[0], out.R, out.L), (shape[1], out.A, out.P)]:
        assert (nvoxels + low + high) % 64 == 0
        assert 0 <= high - low <= 1
    assert (shape[2] + out.S) % 64 == 0
    if shape == (64, 64, 64):
        assert (out.S, out.R, out.L, out.A, out.P) == (0, 0, 0, 0, 0)


def test_synthstrip_padding_rejects_4d(tmp_path):
    """SynthStrip needs a 3D image."""
    in_file = tmp_path / 'image.nii.gz'
    nb.Nifti1Image(np.zeros((4, 4, 4, 2), dtype=np.float32), np.eye(4)).to_filename(in_file)

    with pytest.raises(ValueError, match='3D inputs'):
        SynthStripPadding(input_image=str(in_file)).run(cwd=str(tmp_path))


def test_init_dl_prep_wf():
    """The padding amounts are computed by one node and applied by an AFNI Zeropad node."""
    from aslprep.workflows.segmentation import init_dl_prep_wf

    wf = init_dl_prep_wf()

    assert isinstance(wf.get_node('synthstrip_padding').interface, SynthStripPadding)
    zeropad = wf.get_node('prepare_synthstrip_reference')
    assert isinstance(zeropad.interface, Zeropad)
    zeropad_inputs = _edges(wf, 'prepare_synthstrip_reference')
    assert zeropad_inputs['in_files'] == ('skulled_autobox', 'out_file')
    for side in 'SRLAP':
        assert zeropad_inputs[side] == ('synthstrip_padding', side)
    assert _edges(wf, 'resample_skulled_to_reference')['reference_image'] == (
        'prepare_synthstrip_reference',
        'out_file',
    )


@pytest.mark.parametrize('sloppy', [False, True])
def test_init_synthstrip_wf(sloppy):
    """Sloppy mode uses BET, connected through its own field names, in place of SynthStrip."""
    from aslprep.workflows.segmentation import init_synthstrip_wf

    reset_config()
    with mock_config():
        config.execution.sloppy = sloppy
        wf = init_synthstrip_wf(do_padding=False)

    name, interface, in_field, mask_field = (
        ('bet', fsl.BET, 'in_file', 'mask_file')
        if sloppy
        else ('synthstrip', FixHeaderSynthStrip, 'input_image', 'out_brain_mask')
    )
    node = wf.get_node(name)
    assert isinstance(node.interface, interface)
    assert _edges(wf, name)[in_field] == ('inputnode', 'padded_image')
    assert _edges(wf, 'mask_to_original_grid')['input_image'] == (name, mask_field)
    if sloppy:
        assert node.inputs.mask is True


@pytest.mark.skipif(
    shutil.which('mri_synthstrip') is None or shutil.which('3dZeropad') is None,
    reason='SynthStrip and AFNI are only installed in the ASLPrep image.',
)
def test_synthstrip_wf_smoke(tmp_path):
    """Run the real SynthStrip workflow on a heavily downsampled head image (gh-614).

    Every integration test runs in sloppy mode, where BET stands in for SynthStrip,
    so this is the only test that runs SynthStrip itself.
    """
    from nibabel.processing import resample_to_output
    from templateflow import api as tflow

    from aslprep.workflows.segmentation import init_synthstrip_wf

    template = 'MNI152NLin2009cAsym'
    head = nb.load(tflow.get(template, resolution=2, desc=None, suffix='T1w'))
    brain = nb.load(tflow.get(template, resolution=2, desc='brain', suffix='mask'))
    head_file = tmp_path / 'head.nii.gz'
    resample_to_output(head, voxel_sizes=4, order=1).to_filename(head_file)
    truth = np.asanyarray(resample_to_output(brain, voxel_sizes=4, order=0).dataobj) > 0

    reset_config()
    with mock_config():
        config.execution.sloppy = False
        wf = init_synthstrip_wf()
        wf.base_dir = str(tmp_path / 'work')
        wf.inputs.inputnode.original_image = str(head_file)
        graph = wf.run(plugin='Linear')

    node = next(n for n in graph.nodes() if n.name == 'mask_to_original_grid')
    mask_img = nb.load(node.result.outputs.output_image)
    mask = np.asanyarray(mask_img.dataobj) > 0

    assert mask_img.shape == truth.shape
    assert np.allclose(mask_img.affine, nb.load(head_file).affine)
    dice = 2 * np.sum(mask & truth) / (mask.sum() + truth.sum())
    assert dice > 0.85, f'Dice with the template brain mask is {dice:.3f}'
