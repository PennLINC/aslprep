from aslprep.niworkflows.utils.spaces import Reference, SpatialReferences
from aslprep.smriprep.workflows.anatomical import init_anat_preproc_wf
wf = init_anat_preproc_wf(
    bids_root='.',
    freesurfer=True,
    hires=True,
    longitudinal=False,
    num_t1w=1,
    omp_nthreads=1,
    output_dir='.',
    skull_strip_template=Reference('MNI152NLin2009cAsym'),
    spaces=SpatialReferences([
        ('MNI152Lin', {}),
        ('fsaverage', {'density': '10k'}),
        ('T1w', {}),
        ('fsnative', {})
    ]),
    reportlets_dir='.',
    skull_strip_fixed_seed=False,
)