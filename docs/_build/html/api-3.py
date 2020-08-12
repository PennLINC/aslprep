from aslprep.workflows.tests import mock_config
from aslprep import config
from aslprep.workflows.bold.base import init_func_preproc_wf
with mock_config():
    bold_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_task-restEyesOpen_asl.nii.gz'
    wf = init_func_preproc_wf(str(bold_file))