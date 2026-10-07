"""Tests for aslprep.workflows.asl.base."""

import pytest

from aslprep import config
from aslprep.tests.tests import mock_config, reset_config


def _build_asl_wf(output_spaces, cifti_output=False, project_goodvoxels=False):
    from aslprep.workflows.asl.base import init_asl_wf

    config.workflow.level = 'full'
    config.workflow.cifti_output = cifti_output
    config.workflow.project_goodvoxels = project_goodvoxels
    config.execution.output_spaces = output_spaces
    config.init_spaces()
    asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
    return init_asl_wf(asl_file=str(asl_file))


@pytest.mark.parametrize('cifti_output', [False, '91k'])
def test_init_asl_wf_surface_templates(cifti_output):
    """Standard surface spaces other than fsaverage are resampled with the Workbench.

    Reproduces nipreps/fmriprep#3461.
    """
    reset_config()
    with mock_config():
        config.workflow.run_reconall = True
        wf = _build_asl_wf(
            'asl MNI152NLin2009cAsym fsLR:den-32k',
            cifti_output=cifti_output,
            project_goodvoxels=True,
        )
        graph = wf._create_flat_graph()
        names = {node.fullname for node in graph.nodes()}

        # One CBF datasink per derivative, template, and density
        assert any('asl_wb_surf_wf.ds_mean_cbf_fsLR_32k' in name for name in names)

        # A single goodvoxels mask feeds every surface resampling
        goodvoxels = [
            name for name in names if name.endswith('goodvoxels_bold_mask_wf.outputnode')
        ]
        assert len(goodvoxels) <= 1
        assert any(name.endswith('ds_goodvoxels_mask') for name in names)
        roi_targets = {
            dst.fullname
            for src, dst in graph.edges()
            if src.fullname.endswith('goodvoxels_bold_mask_wf.outputnode')
        }
        assert any('asl_wb_surf_wf' in name for name in roi_targets)
        if cifti_output:
            assert any('asl_cifti_resample_wf' in name for name in roi_targets)


def test_init_asl_wf_surface_templates_require_freesurfer():
    """Without FreeSurfer surfaces, surface templates are skipped instead of crashing."""
    reset_config()
    with mock_config():
        config.workflow.run_reconall = False
        wf = _build_asl_wf('asl MNI152NLin2009cAsym fsLR:den-32k')
        names = {node.fullname for node in wf._create_flat_graph().nodes()}
        assert not any('asl_wb_surf_wf' in name for name in names)
