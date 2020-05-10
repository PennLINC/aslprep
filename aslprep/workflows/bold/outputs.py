# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Writing out derivative files."""
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
from ...niworkflows.interfaces.cifti import CiftiNameSource
from ...niworkflows.interfaces.surf import GiftiNameSource
from ...niworkflows.interfaces.utility import KeySelect
from ...niworkflows.utils.spaces import format_reference as _fmt_space

from ...config import DEFAULT_MEMORY_MIN_GB
from ...interfaces import DerivativesDataSink


def init_func_derivatives_wf(
    bids_root,
    cifti_output,
    freesurfer,
    metadata,
    output_dir,
    spaces,
    #use_aroma,
    name='func_derivatives_wf',
):
    """
    Set up a battery of datasinks to store derivatives in the right location.

    Parameters
    ----------
    bids_root : :obj:`str`
        Original BIDS dataset path.
    cifti_output : :obj:`bool`
        Whether the ``--cifti-output`` flag was set.
    freesurfer : :obj:`bool`
        Whether FreeSurfer anatomical processing was run.
    metadata : :obj:`dict`
        Metadata dictionary associated to the BOLD run.
    output_dir : :obj:`str`
        Where derivatives should be written out to.
    spaces : :py:class:`~...niworkflows.utils.spaces.SpatialReferences`
        A container for storing, organizing, and parsing spatial normalizations. Composed of
        :py:class:`~...niworkflows.utils.spaces.Reference` objects representing spatial references.
        Each ``Reference`` contains a space, which is a string of either TemplateFlow template IDs
        (e.g., ``MNI152Lin``, ``MNI152NLin6Asym``, ``MNIPediatricAsym``), nonstandard references
        (e.g., ``T1w`` or ``anat``, ``sbref``, ``run``, etc.), or a custom template located in
        the TemplateFlow root directory. Each ``Reference`` may also contain a spec, which is a
        dictionary with template specifications (e.g., a specification of ``{'resolution': 2}``
        would lead to resampling on a 2mm resolution of the space).
    use_aroma : :obj:`bool`
        Whether ``--use-aroma`` flag was set.
    name : :obj:`str`
        This workflow's identifier (default: ``func_derivatives_wf``).

    """
    from  ...smriprep.workflows.outputs import _bids_relative
    nonstd_spaces = set(spaces.get_nonstandard())
    workflow = Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=[
        'bold_aparc_std', 'bold_aparc_t1', 'bold_aseg_std',
        'bold_aseg_t1', 'bold_cifti', 'bold_mask_std', 'bold_mask_t1', 'bold_std',
        'bold_std_ref', 'bold_t1', 'bold_t1_ref', 'bold_native', 'bold_native_ref',
        'bold_mask_native', 'cifti_variant', 'cifti_metadata', 'cifti_density',
        'confounds', 'confounds_metadata', 'source_file', 'surf_files', 'surf_refs', 
        'template', 'spatial_reference','cbf','meancbf','score','avgscore',
        'scrub','basil','pv','cbf_t1','meancbf_t1','att_t1','score_t1','avgscore_t1',
        'scrub_t1','basil_t1','pv_t1','cbf_std','meancbf_std','score_std',
        'avgscore_std','scrub_std','basil_std','pv_std','qc_file',
        'cbf_hvoxf','score_hvoxf','scrub_hvoxf','basil_hvoxf','pvc_hvoxf',
        'cbf_sc207','score_sc207','scrub_sc207','basil_sc207','pvc_sc207',
        'cbf_sc217','score_sc217','scrub_sc217','basil_sc217','pvc_sc217',
        'cbf_sc407','score_sc407','scrub_sc407','basil_sc407','pvc_sc407',
        'cbf_sc417','score_sc417','scrub_sc417','basil_sc417','pvc_sc417']),
        name='inputnode')

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name='raw_sources')
    raw_sources.inputs.bids_root = bids_root

    ds_confounds = pe.Node(DerivativesDataSink(
        base_directory=output_dir, desc='confounds', suffix='regressors'),
        name="ds_confounds", run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)
    
    
    workflow.connect([
        (inputnode, raw_sources, [('source_file', 'in_files')]),
        (inputnode, ds_confounds, [('source_file', 'source_file'),
                                   ('confounds', 'in_file'),
                                   ('confounds_metadata', 'meta_dict')]),
    ])

    qcfile = pe.Node(
            DerivativesDataSink(base_directory=output_dir,desc='quality_control',suffix='cbf', compress=False),
            name='qcfile', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
    workflow.connect([
        #(inputnode, raw_sources, [('source_file', 'in_files')]),
        (inputnode,qcfile,[('source_file', 'source_file'),
                                 ('qc_file', 'in_file')]),
    ])

    cbf_hvoxf = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='HavardOxford',suffix='mean_cbf',
            compress=False),name='cbf_hvoxf',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    score_hvoxf = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='HavardOxford',suffix='mean_score',
            compress=False),name='score_hvoxf',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    scrub_hvoxf = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='HavardOxford',suffix='mean_scrub',
            compress=False),name='scrub_hvoxf',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    basil_hvoxf = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='HavardOxford',suffix='mean_basil',
            compress=False),name='basil_hvoxf',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    pvc_hvoxf = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='HavardOxford',suffix='mean_pvc',
            compress=False),name='pvc_hvoxf',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    
    workflow.connect([
        (inputnode,cbf_hvoxf,[('source_file', 'source_file'),
                                ('cbf_hvoxf', 'in_file')]),
        (inputnode,score_hvoxf,[('source_file', 'source_file'),
                                ('score_hvoxf', 'in_file')]),
        (inputnode,scrub_hvoxf,[('source_file', 'source_file'),
                                ('scrub_hvoxf', 'in_file')]),
        (inputnode,basil_hvoxf,[('source_file', 'source_file'),
                                ('basil_hvoxf', 'in_file')]),
        (inputnode,pvc_hvoxf,[('source_file', 'source_file'),
                                ('pvc_hvoxf', 'in_file')]),
      ])

    
    cbf_sc207 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer200x7',suffix='mean_cbf',
            compress=False),name='cbf_sc207',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    score_sc207 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer200x7',suffix='mean_score',
            compress=False),name='score_sc207',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    scrub_sc207 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer200x7',suffix='mean_scrub',
            compress=False),name='scrub_sc207',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    basil_sc207 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer200x7',suffix='mean_basil',
            compress=False),name='basil_sc207',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    pvc_sc207 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer200x7',suffix='mean_pvc',
            compress=False),name='pvc_sc207',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    
    workflow.connect([
        (inputnode,cbf_sc207,[('source_file', 'source_file'),
                                ('cbf_sc207', 'in_file')]),
        (inputnode,score_sc207,[('source_file', 'source_file'),
                                ('score_sc207', 'in_file')]),
        (inputnode,scrub_sc207,[('source_file', 'source_file'),
                                ('scrub_sc207', 'in_file')]),
        (inputnode,basil_sc207,[('source_file', 'source_file'),
                                ('basil_sc207', 'in_file')]),
        (inputnode,pvc_sc207,[('source_file', 'source_file'),
                                ('pvc_sc207', 'in_file')]),
      ])
    
    cbf_sc217 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer200x17',suffix='mean_cbf',
            compress=False),name='cbf_sc217',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    score_sc217 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer200x17',suffix='mean_score',
            compress=False),name='score_sc217',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    scrub_sc217 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer200x17',suffix='mean_scrub',
            compress=False),name='scrub_sc217',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    basil_sc217 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer200x17',suffix='mean_basil',
            compress=False),name='basil_sc217',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    pvc_sc217 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer200x17',suffix='mean_pvc',
            compress=False),name='pvc_sc217',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    
    workflow.connect([
        (inputnode,cbf_sc217,[('source_file', 'source_file'),
                                ('cbf_sc217', 'in_file')]),
        (inputnode,score_sc217,[('source_file', 'source_file'),
                                ('score_sc217', 'in_file')]),
        (inputnode,scrub_sc217,[('source_file', 'source_file'),
                                ('scrub_sc217', 'in_file')]),
        (inputnode,basil_sc217,[('source_file', 'source_file'),
                                ('basil_sc217', 'in_file')]),
        (inputnode,pvc_sc217,[('source_file', 'source_file'),
                                ('pvc_sc217', 'in_file')]),
      ])
    

    cbf_sc407 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer400x7',suffix='mean_cbf',
            compress=False),name='cbf_sc407',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    score_sc407 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer400x7',suffix='mean_score',
            compress=False),name='score_sc407',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    scrub_sc407 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer400x7',suffix='mean_scrub',
            compress=False),name='scrub_sc407',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    basil_sc407 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer400x7',suffix='mean_basil',
            compress=False),name='basil_sc407',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    pvc_sc407 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer400x7',suffix='mean_pvc',
            compress=False),name='pvc_sc407',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    
    workflow.connect([
        (inputnode,cbf_sc407,[('source_file', 'source_file'),
                                ('cbf_sc407', 'in_file')]),
        (inputnode,score_sc407,[('source_file', 'source_file'),
                                ('score_sc407', 'in_file')]),
        (inputnode,scrub_sc407,[('source_file', 'source_file'),
                                ('scrub_sc407', 'in_file')]),
        (inputnode,basil_sc407,[('source_file', 'source_file'),
                                ('basil_sc407', 'in_file')]),
        (inputnode,pvc_sc407,[('source_file', 'source_file'),
                                ('pvc_sc407', 'in_file')]),
      ])
    

    cbf_sc417 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer400x17',suffix='mean_cbf',
            compress=False),name='cbf_sc417',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    score_sc417 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer400x17',suffix='mean_score',
            compress=False),name='score_sc417',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    scrub_sc417 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer400x17',suffix='mean_scrub',
            compress=False),name='scrub_sc417',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    basil_sc417 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer400x17',suffix='mean_basil',
            compress=False),name='basil_sc417',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    pvc_sc417 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,label='schaefer400x17',suffix='mean_pvc',
            compress=False),name='pvc_sc417',run_without_submitting=True,mem_gb=DEFAULT_MEMORY_MIN_GB)
    
    workflow.connect([
        (inputnode,cbf_sc417,[('source_file', 'source_file'),
                                ('cbf_sc417', 'in_file')]),
        (inputnode,score_sc417,[('source_file', 'source_file'),
                                ('score_sc417', 'in_file')]),
        (inputnode,scrub_sc417,[('source_file', 'source_file'),
                                ('scrub_sc417', 'in_file')]),
        (inputnode,basil_sc417,[('source_file', 'source_file'),
                                ('basil_sc417', 'in_file')]),
        (inputnode,pvc_sc417,[('source_file', 'source_file'),
                                ('pvc_sc417', 'in_file')]),
      ])


    
    if nonstd_spaces.intersection(('func', 'run', 'bold', 'boldref', 'sbref')):
        ds_bold_native = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='preproc',
                                keep_dtype=True, compress=True, SkullStripped=False,
                                RepetitionTime=metadata.get('RepetitionTime'),
                                TaskName=metadata.get('TaskName')),
            name='ds_bold_native', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_bold_native_ref = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='aslref', compress=True),
            name='ds_bold_native_ref', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_bold_mask_native = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='brain',
                                suffix='mask', compress=True),
            name='ds_bold_mask_native', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        cbfnative = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='cbf', compress=True),
            name='cbfnative', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        meancbfnative = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='mean_cbf', compress=True),
            name='meancbfnative', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        scorenative = pe.Node(
            DerivativesDataSink(base_directory=output_dir,desc='score',suffix='cbf', compress=True),
            name='scorenative', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        meanscorenative = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='score',suffix='mean_cbf', compress=True),
            name='meanscorenative', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        scrubnative = pe.Node(
            DerivativesDataSink(base_directory=output_dir,desc='scrub',suffix='cbf', compress=True),
            name='scrubnative', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        basilnative = pe.Node(
            DerivativesDataSink(base_directory=output_dir,desc='basil',suffix='cbf', compress=True),
            name='basilnative', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        pvnative = pe.Node(
            DerivativesDataSink(base_directory=output_dir,desc='pvc',suffix='cbf', compress=True),
            name='pvcnative', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
       
        

        workflow.connect([
            (inputnode, ds_bold_native, [('source_file', 'source_file'),
                                         ('bold_native', 'in_file')]),
            (inputnode, ds_bold_native_ref, [('source_file', 'source_file'),
                                             ('bold_native_ref', 'in_file')]),
            (inputnode, ds_bold_mask_native, [('source_file', 'source_file'),
                                              ('bold_mask_native', 'in_file')]),
            (inputnode,cbfnative,[('source_file', 'source_file'),
                                              ('cbf', 'in_file')]),
            (inputnode,meancbfnative,[('source_file', 'source_file'),
                                              ('meancbf', 'in_file')]),
            (inputnode,scorenative,[('source_file', 'source_file'),
                                              ('score', 'in_file')]),
            (inputnode,meanscorenative,[('source_file', 'source_file'),
                                              ('avgscore', 'in_file')]),
            (inputnode,scrubnative,[('source_file', 'source_file'),
                                              ('scrub', 'in_file')]),
            (inputnode,basilnative,[('source_file', 'source_file'),
                                              ('basil', 'in_file')]),
            (inputnode,pvnative,[('source_file', 'source_file'),
                                              ('pv', 'in_file')]),
            (raw_sources, ds_bold_mask_native, [('out', 'RawSources')]),
        ])
  

    # Resample to T1w space
    if nonstd_spaces.intersection(('T1w', 'anat')):
        ds_bold_t1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, space='T1w', desc='preproc',
                                keep_dtype=True, compress=True, SkullStripped=False,
                                RepetitionTime=metadata.get('RepetitionTime'),
                                TaskName=metadata.get('TaskName')),
            name='ds_bold_t1', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_bold_t1_ref = pe.Node(
            DerivativesDataSink(base_directory=output_dir, space='T1w',
                                suffix='aslref', compress=True),
            name='ds_bold_t1_ref', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)

        ds_bold_mask_t1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, space='T1w', desc='brain',
                                suffix='mask', compress=True),
            name='ds_bold_mask_t1', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        
        cbfnativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='cbf', space='T1w',compress=True),
            name='cbfnativet1', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        meancbfnativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,suffix='mean_cbf', space='T1w',compress=True),
            name='meancbfnativet1', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        scorenativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,desc='score',suffix='cbf',space='T1w', compress=True),
            name='scorenativet1', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        meanscorenativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='mean_cbf',desc='score',space='T1w', compress=True),
            name='meanscorenativet1', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        scrubnativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir,desc='scrub',suffix='cbf',space='T1w', compress=True),
            name='scrubnativet1', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        basilnativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='basil',suffix='cbf', space='T1w',compress=True),
            name='basilnativet1', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        pvnativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='pvc',suffix='cbf',space='T1w', compress=True),
            name='pvcnativet1', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)

        
        workflow.connect([
            (inputnode, ds_bold_t1, [('source_file', 'source_file'),
                                     ('bold_t1', 'in_file')]),
            (inputnode, ds_bold_t1_ref, [('source_file', 'source_file'),
                                         ('bold_t1_ref', 'in_file')]),
            (inputnode, ds_bold_mask_t1, [('source_file', 'source_file'),
                                          ('bold_mask_t1', 'in_file')]),
            (inputnode,cbfnativet1,[('source_file', 'source_file'),
                                              ('cbf_t1', 'in_file')]),
            (inputnode,meancbfnativet1,[('source_file', 'source_file'),
                                              ('meancbf_t1', 'in_file')]),
            (inputnode,scorenativet1,[('source_file', 'source_file'),
                                              ('score_t1', 'in_file')]),
            (inputnode,meanscorenativet1,[('source_file', 'source_file'),
                                              ('avgscore_t1', 'in_file')]),
            (inputnode,scrubnativet1,[('source_file', 'source_file'),
                                              ('scrub_t1', 'in_file')]),
            (inputnode,basilnativet1,[('source_file', 'source_file'),
                                              ('basil_t1', 'in_file')]),
            (inputnode,pvnativet1,[('source_file', 'source_file'),
                                              ('pv_t1', 'in_file')]),
            (raw_sources, ds_bold_mask_t1, [('out', 'RawSources')]),
        ])
 

        if freesurfer:
            ds_bold_aseg_t1 = pe.Node(DerivativesDataSink(
                base_directory=output_dir, space='T1w', desc='aseg', suffix='dseg'),
                name='ds_bold_aseg_t1', run_without_submitting=True,
                mem_gb=DEFAULT_MEMORY_MIN_GB)
            ds_bold_aparc_t1 = pe.Node(DerivativesDataSink(
                base_directory=output_dir, space='T1w', desc='aparcaseg', suffix='dseg'),
                name='ds_bold_aparc_t1', run_without_submitting=True,
                mem_gb=DEFAULT_MEMORY_MIN_GB)
            workflow.connect([
                (inputnode, ds_bold_aseg_t1, [('source_file', 'source_file'),
                                              ('bold_aseg_t1', 'in_file')]),
                (inputnode, ds_bold_aparc_t1, [('source_file', 'source_file'),
                                               ('bold_aparc_t1', 'in_file')]),
            ])
#"""
    #if use_aroma:
        #ds_aroma_noise_ics = pe.Node(DerivativesDataSink(
         #   base_directory=output_dir, suffix='AROMAnoiseICs'),
          #  name="ds_aroma_noise_ics", run_without_submitting=True,
           # mem_gb=DEFAULT_MEMORY_MIN_GB)
        #ds_melodic_mix = pe.Node(DerivativesDataSink(
         #   base_directory=output_dir, desc='MELODIC', suffix='mixing'),
          #  name="ds_melodic_mix", run_without_submitting=True,
          #  mem_gb=DEFAULT_MEMORY_MIN_GB)
        #ds_aroma_std = pe.Node(
        #    DerivativesDataSink(base_directory=output_dir, space='MNI152NLin6Asym',
        #                        desc='smoothAROMAnonaggr', keep_dtype=True),
        #    name='ds_aroma_std', run_without_submitting=True,
        #   mem_gb=DEFAULT_MEMORY_MIN_GB)

        #workflow.connect([
        #    (inputnode, ds_aroma_noise_ics, [('source_file', 'source_file'),
        #                                     ('aroma_noise_ics', 'in_file')]),
        #    (inputnode, ds_melodic_mix, [('source_file', 'source_file'),
        #                                 ('melodic_mix', 'in_file')]),
        #    (inputnode, ds_aroma_std, [('source_file', 'source_file'),
        #                               ('nonaggr_denoised_file', 'in_file')]),
        #])

    #if getattr(spaces, '_cached') is None:
    #    return workflow

    # Store resamplings in standard spaces when listed in --output-spaces
    if spaces.cached.references:
        itersource = pe.Node(niu.IdentityInterface(fields=['space_definition']),
                             name='itersource')

        itersource.iterables = (
            'space_definition',
            [(s.fullname, s.spec) for s in spaces.cached.get_standard(dim=(3,))]
        )

        select_std = pe.Node(KeySelect(
            fields=['template', 'bold_std','bold_std_ref','bold_mask_std',
            'cbf_std','meancbf_std','score_std','avgscore_std','scrub_std',
            'basil_std','pv_std']),
            name='select_std', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)

        ds_bold_std = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='preproc',
                                keep_dtype=True, compress=True, SkullStripped=False,
                                RepetitionTime=metadata.get('RepetitionTime'),
                                TaskName=metadata.get('TaskName')),
            name='ds_bold_std', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_bold_std_ref = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='aslref'),
            name='ds_bold_std_ref', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_bold_mask_std = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='brain',
                                suffix='mask'),
            name='ds_bold_mask_std', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)

        cbfstd = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='cbf', compress=True),
            name='cbfstd', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        meancbfstd = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='mean_cbf', compress=True),
            name='meancbfstd', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        scorestd = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='score',suffix='cbf', compress=True),
            name='scorestd', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        meanscorestd = pe.Node(
            DerivativesDataSink(base_directory=output_dir,desc='score', suffix='mean_cbf', compress=True),
            name='meanscorestd', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        scrubstd = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='scrub',suffix='cbf', compress=True),
            name='scrubstd', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        basilstd = pe.Node(
            DerivativesDataSink(base_directory=output_dir,desc='basil',suffix='cbf', compress=True),
            name='basilstd', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        pvstd = pe.Node(
            DerivativesDataSink(base_directory=output_dir,desc='pvc',suffix='cbf', compress=True),
            name='pvcstd', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
    

        workflow.connect([
            (inputnode, ds_bold_std, [('source_file', 'source_file')]),
            (inputnode, ds_bold_std_ref, [('source_file', 'source_file')]),
            (inputnode, ds_bold_mask_std,[('source_file', 'source_file')]),
            (inputnode, cbfstd, [('source_file', 'source_file')]),
            (inputnode, meancbfstd, [('source_file', 'source_file')]),
            (inputnode, scorestd, [('source_file', 'source_file')]),
            (inputnode, meanscorestd, [('source_file', 'source_file')]),
            (inputnode, scrubstd, [('source_file', 'source_file')]),
            (inputnode, basilstd, [('source_file', 'source_file')]),
            (inputnode, pvstd, [('source_file', 'source_file')]),
            (inputnode, select_std, [('bold_std', 'bold_std'),
                                     ('bold_std_ref', 'bold_std_ref'),
                                     ('bold_mask_std', 'bold_mask_std'),
                                     ('cbf_std', 'cbf_std'),
                                     ('meancbf_std', 'meancbf_std'),
                                     ('score_std', 'score_std'),
                                     ('avgscore_std', 'avgscore_std'),
                                     ('scrub_std', 'scrub_std'),
                                     ('basil_std', 'basil_std'),
                                     ('pv_std', 'pv_std'),
                                     ('template', 'template'),
                                     ('spatial_reference', 'keys')]),
            (itersource, select_std, [(('space_definition', _fmt_space), 'key')]),
            (select_std, ds_bold_std, [('bold_std', 'in_file'),
                                       ('key', 'space')]),
            (select_std, ds_bold_std_ref, [('bold_std_ref', 'in_file'),
                                           ('key', 'space')]),
            (select_std, ds_bold_mask_std, [('bold_mask_std', 'in_file'),
                                            ('key', 'space')]),
            (select_std, cbfstd, [('cbf_std', 'in_file'),
                                            ('key', 'space')]),
            (select_std, meancbfstd, [('meancbf_std', 'in_file'),
                                            ('key', 'space')]),
            (select_std, scorestd, [('score_std', 'in_file'),
                                            ('key', 'space')]),
            (select_std, meanscorestd, [('avgscore_std', 'in_file'),
                                            ('key', 'space')]),
            (select_std, scrubstd, [('scrub_std', 'in_file'),
                                            ('key', 'space')]),
            (select_std, basilstd, [('basil_std', 'in_file'),
                                            ('key', 'space')]),
            (select_std, pvstd, [('pv_std', 'in_file'),
                                            ('key', 'space')]),
            (raw_sources, ds_bold_mask_std, [('out', 'RawSources')]),
        ])

    

        if freesurfer:
            select_fs_std = pe.Node(KeySelect(
                fields=['bold_aseg_std', 'bold_aparc_std', 'template']),
                name='select_fs_std', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
            ds_bold_aseg_std = pe.Node(DerivativesDataSink(
                base_directory=output_dir, desc='aseg', suffix='dseg'),
                name='ds_bold_aseg_std', run_without_submitting=True,
                mem_gb=DEFAULT_MEMORY_MIN_GB)
            ds_bold_aparc_std = pe.Node(DerivativesDataSink(
                base_directory=output_dir, desc='aparcaseg', suffix='dseg'),
                name='ds_bold_aparc_std', run_without_submitting=True,
                mem_gb=DEFAULT_MEMORY_MIN_GB)
            workflow.connect([
                (itersource, select_fs_std, [
                    (('space_definition', _fmt_space), 'key')]),
                (inputnode, select_fs_std, [('bold_aseg_std', 'bold_aseg_std'),
                                            ('bold_aparc_std', 'bold_aparc_std'),
                                            ('template', 'template'),
                                            ('spatial_reference', 'keys')]),
                (select_fs_std, ds_bold_aseg_std, [('bold_aseg_std', 'in_file'),
                                                   ('key', 'space')]),
                (select_fs_std, ds_bold_aparc_std, [('bold_aparc_std', 'in_file'),
                                                    ('key', 'space')]),
                (inputnode, ds_bold_aseg_std, [('source_file', 'source_file')]),
                (inputnode, ds_bold_aparc_std, [('source_file', 'source_file')])
            ])

    fs_outputs = spaces.cached.get_fs_spaces()
    if freesurfer and fs_outputs:
        select_fs_surf = pe.Node(KeySelect(
            fields=['surfaces', 'surf_kwargs']), name='select_fs_surf',
            run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        select_fs_surf.iterables = [('key', fs_outputs)]
        select_fs_surf.inputs.surf_kwargs = [{'space': s} for s in fs_outputs]

        name_surfs = pe.MapNode(GiftiNameSource(
            pattern=r'(?P<LR>[lr])h.\w+',
            template='space-{space}_hemi-{LR}.func'),
            iterfield=['in_file'], name='name_surfs',
            mem_gb=DEFAULT_MEMORY_MIN_GB, run_without_submitting=True)

        ds_bold_surfs = pe.MapNode(DerivativesDataSink(base_directory=output_dir),
                                   iterfield=['in_file', 'suffix'], name='ds_bold_surfs',
                                   run_without_submitting=True,
                                   mem_gb=DEFAULT_MEMORY_MIN_GB)

        workflow.connect([
            (inputnode, select_fs_surf, [
                ('surf_files', 'surfaces'),
                ('surf_refs', 'keys')]),
            (select_fs_surf, name_surfs, [('surfaces', 'in_file'),
                                          ('surf_kwargs', 'template_kwargs')]),
            (inputnode, ds_bold_surfs, [('source_file', 'source_file')]),
            (select_fs_surf, ds_bold_surfs, [('surfaces', 'in_file')]),
            (name_surfs, ds_bold_surfs, [('out_name', 'suffix')]),
        ])

    # CIFTI output
    if cifti_output:
        name_cifti = pe.MapNode(
            CiftiNameSource(), iterfield=['variant', 'density'], name='name_cifti',
            mem_gb=DEFAULT_MEMORY_MIN_GB, run_without_submitting=True)
        cifti_bolds = pe.MapNode(
            DerivativesDataSink(base_directory=output_dir, compress=False),
            iterfield=['in_file', 'suffix'], name='cifti_bolds',
            run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        
        cifti_key = pe.MapNode(DerivativesDataSink(
            base_directory=output_dir), iterfield=['in_file', 'suffix'],
            name='cifti_key', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        workflow.connect([
            (inputnode, name_cifti, [('cifti_variant', 'variant'),
                                     ('cifti_density', 'density')]),
            (inputnode, cifti_bolds, [('bold_cifti', 'in_file'),
                                      ('source_file', 'source_file')]),
            (name_cifti, cifti_bolds, [('out_name', 'suffix')]),
            (name_cifti, cifti_key, [('out_name', 'suffix')]),
            (inputnode, cifti_key, [('source_file', 'source_file'),
                                    ('cifti_metadata', 'in_file')]),
        ])

    return workflow
     