import os
from collections import namedtuple, OrderedDict
BIDSLayout = namedtuple('BIDSLayout', ['root'])
from aslprep.workflows.base import init_aslprep_wf
from ..niworkflows.utils.spaces import Reference, SpatialReferences

os.environ['FREESURFER_HOME'] = os.getcwd()
wf = init_aslprep_wf(
    anat_only=False,
    aroma_melodic_dim=-200,
    bold2t1w_dof=9,
    cifti_output=False,
    debug=False,
    dummy_scans=None,
    echo_idx=None,
    err_on_aroma_warn=False,
    fmap_bspline=False,
    fmap_demean=True,
    force_syn=True,
    freesurfer=True,
    fs_subjects_dir=None,
    hires=True,
    ignore=[],
    layout=BIDSLayout('.'),
    longitudinal=False,
    low_mem=False,
    medial_surface_nan=False,
    omp_nthreads=1,
    output_dir='.',
    regressors_all_comps=False,
    regressors_dvars_th=1.5,
    regressors_fd_th=0.5,
    run_uuid='X',
    skull_strip_fixed_seed=False,
    skull_strip_template=Reference('OASIS30ANTs'),
    spaces=SpatialReferences(
        spaces=['MNI152Lin',
                ('fsaverage', {'density': '10k'}),
                'T1w',
                'fsnative'],
        checkpoint=True),
    subject_list=['aslpreptest'],
    t2s_coreg=False,
    task_id='',
    use_aroma=False,
    use_bbr=True,
    use_syn=True,
    work_dir='.',
    bids_filters=None,
)