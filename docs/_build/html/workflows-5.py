from collections import namedtuple
from aslprep.niworkflows.utils.spaces import SpatialReferences
from aslprep.workflows.bold.base import init_func_preproc_wf
BIDSLayout = namedtuple('BIDSLayout', ['root'])
wf = init_func_preproc_wf(
    bold2t1w_dof=9,
    bold_file='/completely/made/up/path/sub-01_task-rest_asl.nii.gz',
    cifti_output=False,
    debug=False,
    dummy_scans=None,
    err_on_aroma_warn=False,
    fmap_bspline=True,
    fmap_demean=True,
    force_syn=True,
    freesurfer=True,
    ignore=[],
    low_mem=False,
    medial_surface_nan=False,
    omp_nthreads=1,
    output_dir='.',
    reportlets_dir='.',
    t2s_coreg=False,
    spaces=SpatialReferences(
        spaces=[('MNI152Lin', {}),
                ('fsaverage', {'density': '10k'}),
                ('T1w', {}),
                ('fsnative', {})],
        checkpoint=True),
    use_bbr=True,
    use_syn=True,
    layout=BIDSLayout('.'),
    num_bold=1,
)