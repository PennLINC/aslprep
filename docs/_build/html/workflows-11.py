from pathlib import Path
from pkg_resources import resource_filename as pkgrf
bids_dir=Path(pkgrf('aslprep', 'data/tests/ds000240')).absolute()
from aslprep.workflows.bold.cbf import init_cbfqc_compt_wf
bold_file = bids_dir / 'sub-01' / 'perf'/ 'sub-01_task-restEyesOpen_asl.nii.gz'
metadata = bids_dir / 'sub-01' / 'perf'/ 'sub-01_task-restEyesOpen_asl.json'
wf = init_cbfqc_compt_wf(mem_gb=0.1,bold_file=str(bold_file),metadata=str(metadata),omp_nthreads=1)