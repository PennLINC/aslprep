"""Tests for ASL workflow wiring."""

import pytest

from aslprep import config
from aslprep.tests.tests import mock_config, reset_config


def test_remove_keys():
    """Patched metadata keys are removed from preprocessed ASL sidecars."""
    from aslprep.workflows.asl.outputs import _remove_keys

    metadata = {'RepetitionTime': 4.0, 'PostLabelingDelay': [1.8, 1.8]}
    assert _remove_keys(metadata, ['RepetitionTime']) == {'PostLabelingDelay': [1.8, 1.8]}
    assert _remove_keys(metadata, []) == metadata


@pytest.mark.parametrize('level', ['minimal', 'full'])
def test_init_asl_wf_output_metadata(level):
    """Preprocessed ASL sidecars use the reduced series' metadata."""
    reset_config()
    with mock_config():
        config.workflow.level = level
        config.workflow.cifti_output = False
        config.execution.output_spaces = 'asl T1w MNI152NLin2009cAsym'
        config.init_spaces()
        from aslprep.workflows.asl.base import init_asl_wf

        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        wf = init_asl_wf(asl_file=str(asl_file))
        graph = wf._create_flat_graph()
        targets = {
            dst.fullname
            for src, dst in graph.edges()
            if src.fullname.endswith('asl_output_metadata')
        }
        assert any('ds_asl_native_wf' in name for name in targets)
        if level == 'full':
            assert any('ds_asl_t1_wf' in name for name in targets)
            assert any('ds_asl_std_wf' in name for name in targets)


def test_init_asl_fit_wf_fallback_total_readout_time():
    """The TotalReadoutTime fallback reaches the ASL distortion parameters.

    Reproduces the ASL side of nipreps/fmriprep#3423.
    """
    from aslprep.workflows.asl.fit import init_asl_fit_wf, init_asl_native_wf

    reset_config()
    with mock_config():
        config.workflow.fallback_total_readout_time = 0.05
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        fit_wf = init_asl_fit_wf(
            asl_file=str(asl_file),
            aslcontext=str(asl_file).replace('.nii.gz', 'context.tsv'),
            m0scan=None,
            use_ge=False,
            fieldmap_id='auto_00000',
        )
        assert fit_wf.get_node('distortion_params').inputs.fallback == 0.05
        native_wf = init_asl_native_wf(asl_file=str(asl_file), fieldmap_id='auto_00000')
        assert native_wf.get_node('distortion_params').inputs.fallback == 0.05
