"""Tests for aslprep.workflows.asl.base."""

import pytest

from aslprep import config
from aslprep.tests.tests import mock_config, reset_config


def _build_asl_wf(
    output_spaces,
    cifti_output=False,
    project_goodvoxels=False,
    output_timeseries=False,
):
    from aslprep.workflows.asl.base import init_asl_wf

    config.workflow.level = 'full'
    config.workflow.cifti_output = cifti_output
    config.workflow.project_goodvoxels = project_goodvoxels
    config.workflow.output_timeseries = output_timeseries
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
            output_timeseries=True,
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


@pytest.mark.parametrize('output_timeseries', [False, True])
@pytest.mark.parametrize('project_goodvoxels', [False, True])
def test_init_asl_wf_output_timeseries(output_timeseries, project_goodvoxels):
    """ASL and CBF time series are only resampled and written with --output-timeseries."""
    reset_config()
    with mock_config():
        config.workflow.run_reconall = True
        config.workflow.scorescrub = True
        wf = _build_asl_wf(
            'asl T1w MNI152NLin2009cAsym fsLR:den-32k',
            cifti_output='91k',
            project_goodvoxels=project_goodvoxels,
            output_timeseries=output_timeseries,
        )
        names = {node.fullname.split('.', 1)[1] for node in wf._create_flat_graph().nodes()}

        timeseries_nodes = [
            'ds_asl_native_wf.ds_asl',
            'ds_asl_native_wf.ds_cbf_ts',
            'ds_asl_native_wf.ds_cbf_ts_score',
            'ds_asl_t1_wf.ds_asl',
            'ds_asl_t1_wf.ds_cbf_ts',
            'ds_asl_std_wf.ds_asl',
            'ds_asl_std_wf.ds_cbf_ts_score',
            'asl_std_wf.resample',
            'ds_asl_cifti_wf.ds_asl_cifti',
            'asl_cifti_resample_wf.asl_MNI6_wf.resample',
        ]
        for name in timeseries_nodes:
            assert (name in names) is output_timeseries, name

        # The goodvoxels mask needs the T1w-space ASL series
        assert ('asl_anat_wf.resample' in names) is (output_timeseries or project_goodvoxels)

        # Mean CBF and the resampling references are always produced
        for name in [
            'ds_asl_native_wf.ds_mean_cbf',
            'ds_asl_t1_wf.ds_mean_cbf',
            'ds_asl_std_wf.ds_mean_cbf',
            'ds_asl_std_wf.ds_mask',
            'asl_anat_wf.gen_ref',
            'asl_std_wf.gen_ref',
        ]:
            assert name in names, name


def test_init_asl_wf_surface_templates_require_freesurfer():
    """Without FreeSurfer surfaces, surface templates are skipped instead of crashing."""
    reset_config()
    with mock_config():
        config.workflow.run_reconall = False
        wf = _build_asl_wf('asl MNI152NLin2009cAsym fsLR:den-32k')
        names = {node.fullname for node in wf._create_flat_graph().nodes()}
        assert not any('asl_wb_surf_wf' in name for name in names)
