from niworkflows.utils.spaces import SpatialReferences
from aslprep.workflows.bold import init_bold_std_trans_wf
wf = init_bold_std_trans_wf(
    freesurfer=True,
    mem_gb=3,
    omp_nthreads=1,
    spaces=SpatialReferences(
        spaces=['MNI152Lin',
                ('MNIPediatricAsym', {'cohort': '6'})],
        checkpoint=True),
)