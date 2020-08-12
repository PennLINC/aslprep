from aslprep.workflows.tests import mock_config
from aslprep.workflows.base import init_aslprep_wf
with mock_config():
    wf = init_aslprep_wf()