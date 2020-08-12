from aslprep.workflows.tests import mock_config
from aslprep.workflows.base import init_single_subject_wf
with mock_config():
    wf = init_single_subject_wf('01')