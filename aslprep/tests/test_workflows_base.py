"""Tests for aslprep.workflows.base."""

from aslprep import config
from aslprep.tests.tests import mock_config, reset_config


def test_init_aslprep_wf_preserves_config():
    """Building the workflow must not modify the global config.

    Reproduces nipreps/fmriprep#3689: ``asl2anat_init`` was overwritten per subject,
    so a T2w-based choice for one subject leaked into subsequent subjects.
    """
    from aslprep.workflows.base import init_aslprep_wf

    reset_config()
    with mock_config():
        config.workflow.asl2anat_init = 'auto'
        before = config.get(flat=True)
        init_aslprep_wf()
        assert config.get(flat=True) == before


def test_init_asl_wf_auto_asl2anat_init():
    """Building an ASL workflow directly resolves an 'auto' asl2anat_init.

    Without subject-level resolution, 'auto' must fall back to 't1w',
    since the FunctionalSummary interface does not accept 'auto'.
    """
    from aslprep.workflows.asl.base import init_asl_wf

    reset_config()
    with mock_config():
        config.workflow.asl2anat_init = 'auto'
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        wf = init_asl_wf(asl_file=str(asl_file))
        summary = wf.get_node('asl_fit_wf.summary')
        assert summary.inputs.registration_init == 't1w'


def test_init_single_subject_wf_datatype_in_path(tmp_path):
    """The error for missing ASL runs explains a dataset path that PyBIDS misreads (gh-704)."""
    import shutil

    import pytest

    from aslprep.workflows.base import init_single_subject_wf

    reset_config()
    with mock_config():
        bids_dir = tmp_path / 'motion' / 'ds000240'
        shutil.copytree(config.execution.bids_dir, bids_dir)
        config.execution.bids_dir = bids_dir
        config.execution.bids_database_dir = None
        config.execution.layout = None
        config.execution._layout = None
        config.execution.init()

        asl_files = config.execution.layout.get(subject='01', suffix='asl', extension='.nii.gz')
        if asl_files[0].entities.get('datatype') == 'perf':
            pytest.skip('This version of PyBIDS only parses entities within the dataset root.')

        with pytest.raises(RuntimeError, match="contains a directory named 'motion'"):
            init_single_subject_wf('01')
