"""Tests for aslprep.config."""

from bids.layout import Query

from aslprep import config
from aslprep.tests.tests import mock_config, reset_config


def test_bids_filters_non_string_values():
    """Non-string BIDS filter values must not crash config initialization.

    Reproduces nipreps/fmriprep#3641.
    """
    reset_config()
    with mock_config():
        config.execution.bids_filters = {
            'asl': {'run': 1, 'echo': [1, 2], 'session': '<Query.NONE: 1>', 'task': None},
        }
        config.execution.init()
        filters = config.execution.bids_filters['asl']
        assert filters['run'] == 1
        assert filters['echo'] == [1, 2]
        assert filters['session'] is Query.NONE
        assert filters['task'] is None
        config.execution.bids_filters = None
