from pathlib import Path
from pkg_resources import resource_filename as pkgrf
bids_dir=Path(pkgrf('aslprep', 'data/tests/ds000240')).absolute()
metadatafile = bids_dir / 'sub-01' / 'perf'/ 'sub-01_task-restEyesOpen_asl.json'
import json
with open(metadatafile) as f:
    metadata = json.load(f)
from aslprep.workflows.bold.cbf import init_cbf_compt_wf
wf = init_cbf_compt_wf(mem_gb=0.1,metadata=metadata, omp_nthreads=4,smooth_kernel=5,dummy_vols=0)