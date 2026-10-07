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
