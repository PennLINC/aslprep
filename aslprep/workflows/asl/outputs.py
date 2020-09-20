# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Writing out derivative files."""
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from ...config import DEFAULT_MEMORY_MIN_GB
from ...interfaces import DerivativesDataSink


def init_asl_derivatives_wf(
    bids_root,
    cifti_output,
    freesurfer,
    metadata,
    output_dir,
    spaces,
    name='asl_derivatives_wf',
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
        Metadata dictionary associated to the ASL run.
    output_dir : :obj:`str`
        Where derivatives should be written out to.
    spaces : :py:class:`~niworkflows.utils.spaces.SpatialReferences`
        A container for storing, organizing, and parsing spatial normalizations. Composed of
        :py:class:`~niworkflows.utils.spaces.Reference` objects representing spatial references.
        Each ``Reference`` contains a space, which is a string of either TemplateFlow template IDs
        (e.g., ``MNI152Lin``, ``MNI152NLin6Asym``, ``MNIPediatricAsym``), nonstandard references
        (e.g., ``T1w`` or ``anat``, ``sbref``, ``run``, etc.), or a custom template located in
        the TemplateFlow root directory. Each ``Reference`` may also contain a spec, which is a
        dictionary with template specifications (e.g., a specification of ``{'resolution': 2}``
        would lead to resampling on a 2mm resolution of the space).
    name : :obj:`str`
        This workflow's identifier (default: ``func_derivatives_wf``).

    """
    from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from ...niworkflows.interfaces.utility import KeySelect
    from ...smriprep.workflows.outputs import _bids_relative

    nonstd_spaces = set(spaces.get_nonstandard())
    workflow = Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=[
        'asl_aparc_std', 'asl_aparc_t1', 'asl_aseg_std',
        'asl_aseg_t1', 'asl_cifti', 'asl_mask_std', 'asl_mask_t1', 'asl_std',
        'asl_std_ref', 'asl_t1', 'asl_t1_ref', 'asl_native', 'asl_native_ref',
        'asl_mask_native', 'cifti_variant', 'cifti_metadata', 'cifti_density',
        'confounds', 'confounds_metadata', 'source_file', 'surf_files', 'surf_refs',
        'template', 'spatial_reference', 'cbf', 'meancbf', 'score', 'avgscore',
        'scrub', 'basil', 'pv', 'cbf_t1', 'meancbf_t1', 'att_t1', 'score_t1', 'avgscore_t1',
        'scrub_t1', 'basil_t1', 'pv_t1', 'cbf_std', 'meancbf_std', 'score_std',
        'avgscore_std', 'scrub_std', 'basil_std', 'pv_std','att','att_std','qc_file',
        'cbf_hvoxf', 'score_hvoxf', 'scrub_hvoxf', 'basil_hvoxf', 'pvc_hvoxf',
        'cbf_sc207', 'score_sc207', 'scrub_sc207', 'basil_sc207', 'pvc_sc207',
        'cbf_sc217', 'score_sc217', 'scrub_sc217', 'basil_sc217', 'pvc_sc217',
        'cbf_sc407', 'score_sc407', 'scrub_sc407', 'basil_sc407', 'pvc_sc407',
        'cbf_sc417', 'score_sc417', 'scrub_sc417', 'basil_sc417', 'pvc_sc417'
        ]),
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
            DerivativesDataSink(base_directory=output_dir,
                                desc='quality_control',
                                suffix='cbf', compress=False),
            name='qcfile', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
    workflow.connect([
        (inputnode, qcfile, [('source_file', 'source_file'),
                             ('qc_file', 'in_file')]),
    ])

    cbf_hvoxf = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='HavardOxford', suffix='mean_cbf',
                                compress=False),
            name='cbf_hvoxf', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    score_hvoxf = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='HavardOxford',
                                suffix='mean_score', compress=False),
            name='score_hvoxf', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    scrub_hvoxf = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='HavardOxford',
                                suffix='mean_scrub', compress=False),
            name='scrub_hvoxf', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    basil_hvoxf = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='HavardOxford',
                                suffix='mean_basil', compress=False),
            name='basil_hvoxf', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    pvc_hvoxf = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='HavardOxford', suffix='mean_pvc',
                                compress=False),
            name='pvc_hvoxf', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow.connect([
        (inputnode, cbf_hvoxf, [('source_file', 'source_file'),
                                ('cbf_hvoxf', 'in_file')]),
        (inputnode, score_hvoxf, [('source_file', 'source_file'),
                                  ('score_hvoxf', 'in_file')]),
        (inputnode, scrub_hvoxf, [('source_file', 'source_file'),
                                  ('scrub_hvoxf', 'in_file')]),
        (inputnode, basil_hvoxf, [('source_file', 'source_file'),
                                  ('basil_hvoxf', 'in_file')]),
        (inputnode, pvc_hvoxf, [('source_file', 'source_file'),
                                ('pvc_hvoxf', 'in_file')]),
      ])

    cbf_sc207 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer200x7',
                                suffix='mean_cbf', compress=False),
            name='cbf_sc207', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    score_sc207 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer200x7',
                                suffix='mean_score', compress=False),
            name='score_sc207', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    scrub_sc207 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer200x7',
                                suffix='mean_scrub', compress=False),
            name='scrub_sc207', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    basil_sc207 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer200x7',
                                suffix='mean_basil', compress=False),
            name='basil_sc207', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    pvc_sc207 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer200x7', suffix='mean_pvc',
                                compress=False),
            name='pvc_sc207', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow.connect([
        (inputnode, cbf_sc207, [('source_file', 'source_file'),
                                ('cbf_sc207', 'in_file')]),
        (inputnode, score_sc207, [('source_file', 'source_file'),
                                  ('score_sc207', 'in_file')]),
        (inputnode, scrub_sc207, [('source_file', 'source_file'),
                                  ('scrub_sc207', 'in_file')]),
        (inputnode, basil_sc207, [('source_file', 'source_file'),
                                  ('basil_sc207', 'in_file')]),
        (inputnode, pvc_sc207, [('source_file', 'source_file'),
                                ('pvc_sc207', 'in_file')]),
      ])

    cbf_sc217 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer200x17',
                                suffix='mean_cbf', compress=False),
            name='cbf_sc217', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    score_sc217 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer200x17',
                                suffix='mean_score', compress=False),
            name='score_sc217', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    scrub_sc217 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer200x17',
                                suffix='mean_scrub', compress=False),
            name='scrub_sc217', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    basil_sc217 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer200x17',
                                suffix='mean_basil', compress=False),
            name='basil_sc217', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    pvc_sc217 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer200x17',
                                suffix='mean_pvc', compress=False),
            name='pvc_sc217', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow.connect([
        (inputnode, cbf_sc217, [('source_file', 'source_file'),
                                ('cbf_sc217', 'in_file')]),
        (inputnode, score_sc217, [('source_file', 'source_file'),
                                  ('score_sc217', 'in_file')]),
        (inputnode, scrub_sc217, [('source_file', 'source_file'),
                                  ('scrub_sc217', 'in_file')]),
        (inputnode, basil_sc217, [('source_file', 'source_file'),
                                  ('basil_sc217', 'in_file')]),
        (inputnode, pvc_sc217, [('source_file', 'source_file'),
                                ('pvc_sc217', 'in_file')]),
      ])

    cbf_sc407 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer400x7', suffix='mean_cbf',
                                compress=False),
            name='cbf_sc407', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    score_sc407 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer400x7',
                                suffix='mean_score', compress=False),
            name='score_sc407', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    scrub_sc407 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer400x7',
                                suffix='mean_scrub', compress=False),
            name='scrub_sc407', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    basil_sc407 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer400x7',
                                suffix='mean_basil', compress=False),
            name='basil_sc407', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    pvc_sc407 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer400x7', suffix='mean_pvc',
                                compress=False),
            name='pvc_sc407', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow.connect([
        (inputnode, cbf_sc407, [('source_file', 'source_file'),
                                ('cbf_sc407', 'in_file')]),
        (inputnode, score_sc407, [('source_file', 'source_file'),
                                  ('score_sc407', 'in_file')]),
        (inputnode, scrub_sc407, [('source_file', 'source_file'),
                                  ('scrub_sc407', 'in_file')]),
        (inputnode, basil_sc407, [('source_file', 'source_file'),
                                  ('basil_sc407', 'in_file')]),
        (inputnode, pvc_sc407, [('source_file', 'source_file'),
                                ('pvc_sc407', 'in_file')]),
      ])

    cbf_sc417 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer400x17',
                                suffix='mean_cbf', compress=False),
            name='cbf_sc417', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    score_sc417 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer400x17',
                                suffix='mean_score', compress=False),
            name='score_sc417', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    scrub_sc417 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer400x17',
                                suffix='mean_scrub', compress=False),
            name='scrub_sc417', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    basil_sc417 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer400x17',
                                suffix='mean_basil', compress=False),
            name='basil_sc417', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
    pvc_sc417 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='schaefer400x17',
                                suffix='mean_pvc', compress=False),
            name='pvc_sc417', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow.connect([
        (inputnode, cbf_sc417, [('source_file', 'source_file'),
                                ('cbf_sc417', 'in_file')]),
        (inputnode, score_sc417, [('source_file', 'source_file'),
                                  ('score_sc417', 'in_file')]),
        (inputnode, scrub_sc417, [('source_file', 'source_file'),
                                  ('scrub_sc417', 'in_file')]),
        (inputnode, basil_sc417, [('source_file', 'source_file'),
                                  ('basil_sc417', 'in_file')]),
        (inputnode, pvc_sc417, [('source_file', 'source_file'),
                                ('pvc_sc417', 'in_file')]),
      ])

    if nonstd_spaces.intersection(('func', 'run', 'sbref')):
        ds_asl_native = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir, desc='preproc', compress=True, SkullStripped=False,
                RepetitionTime=metadata.get('RepetitionTime'), TaskName=metadata.get('TaskName')),
            name='ds_asl_native', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_asl_native_ref = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='aslref', compress=True,
                                dismiss_entities=("echo",)),
            name='ds_asl_native_ref', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_asl_mask_native = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='brain', suffix='mask',
                                compress=True, dismiss_entities=("echo",)),
            name='ds_asl_mask_native', run_without_submitting=True,
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
            DerivativesDataSink(base_directory=output_dir, desc='score', suffix='cbf',
                                compress=True),
            name='scorenative', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        meanscorenative = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='score', suffix='mean_cbf',
                                compress=True),
            name='meanscorenative', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        scrubnative = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='scrub', suffix='cbf',
                                compress=True),
            name='scrubnative', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        basilnative = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='basil', suffix='cbf',
                                compress=True),
            name='basilnative', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        pvnative = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='pvc', suffix='cbf',
                                compress=True),
            name='pvcnative', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        attnative = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='bat', suffix='cbf',
                                compress=True),
            name='attcnative', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)

        workflow.connect([
            (inputnode, ds_asl_native, [('source_file', 'source_file'),
                                         ('asl_native', 'in_file')]),
            (inputnode, ds_asl_native_ref, [('source_file', 'source_file'),
                                             ('asl_native_ref', 'in_file')]),
            (inputnode, ds_asl_mask_native, [('source_file', 'source_file'),
                                              ('asl_mask_native', 'in_file')]),
            (inputnode, cbfnative, [('source_file', 'source_file'),
                                    ('cbf', 'in_file')]),
            (inputnode, meancbfnative, [('source_file', 'source_file'),
                                        ('meancbf', 'in_file')]),
            (inputnode, scorenative, [('source_file', 'source_file'),
                                      ('score', 'in_file')]),
            (inputnode, meanscorenative, [('source_file', 'source_file'),
                                          ('avgscore', 'in_file')]),
            (inputnode, scrubnative, [('source_file', 'source_file'),
                                      ('scrub', 'in_file')]),
            (inputnode, basilnative, [('source_file', 'source_file'),
                                      ('basil', 'in_file')]),
            (inputnode, pvnative, [('source_file', 'source_file'),
                                   ('pv', 'in_file')]),
            (inputnode, attnative, [('source_file', 'source_file'),
                                   ('att', 'in_file')]),
            (raw_sources, ds_asl_mask_native, [('out', 'RawSources')]),
        ])

    # Resample to T1w space
    if nonstd_spaces.intersection(('T1w', 'anat')):
        ds_asl_t1 = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir, space='T1w', desc='preproc', compress=True,
                SkullStripped=False, RepetitionTime=metadata.get('RepetitionTime'),
                TaskName=metadata.get('TaskName'), dismiss_entities=("echo",)),
            name='ds_asl_t1', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_asl_t1_ref = pe.Node(
            DerivativesDataSink(base_directory=output_dir, space='T1w', suffix='aslref',
                                compress=True, dismiss_entities=("echo",)),
            name='ds_asl_t1_ref', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)

        ds_asl_mask_t1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, space='T1w', desc='brain',
                                suffix='mask', compress=True, dismiss_entities=("echo",)),
            name='ds_asl_mask_t1', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)

        cbfnativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='cbf', space='T1w',
                                compress=True),
            name='cbfnativet1', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        meancbfnativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='mean_cbf', space='T1w',
                                compress=True),
            name='meancbfnativet1', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        scorenativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='score', suffix='cbf',
                                space='T1w', compress=True),
            name='scorenativet1', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        meanscorenativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='mean_cbf', desc='score',
                                space='T1w', compress=True),
            name='meanscorenativet1', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        scrubnativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='scrub', suffix='cbf',
                                space='T1w', compress=True),
            name='scrubnativet1', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        basilnativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='basil', suffix='cbf',
                                space='T1w', compress=True),
            name='basilnativet1', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        pvnativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='pvc', suffix='cbf',
                                space='T1w', compress=True),
            name='pvcnativet1', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        attnativet1 = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='bat', suffix='cbf',
                                space='T1w', compress=True),
            name='attnativet1', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)

        workflow.connect([
            (inputnode, ds_asl_t1, [('source_file', 'source_file'),
                                     ('asl_t1', 'in_file')]),
            (inputnode, ds_asl_t1_ref, [('source_file', 'source_file'),
                                         ('asl_t1_ref', 'in_file')]),
            (inputnode, ds_asl_mask_t1, [('source_file', 'source_file'),
                                          ('asl_mask_t1', 'in_file')]),
            (inputnode, cbfnativet1, [('source_file', 'source_file'),
                                      ('cbf_t1', 'in_file')]),
            (inputnode, meancbfnativet1, [('source_file', 'source_file'),
                                          ('meancbf_t1', 'in_file')]),
            (inputnode, scorenativet1, [('source_file', 'source_file'),
                                        ('score_t1', 'in_file')]),
            (inputnode, meanscorenativet1, [('source_file', 'source_file'),
                                            ('avgscore_t1', 'in_file')]),
            (inputnode, scrubnativet1, [('source_file', 'source_file'),
                                        ('scrub_t1', 'in_file')]),
            (inputnode, basilnativet1, [('source_file', 'source_file'),
                                        ('basil_t1', 'in_file')]),
            (inputnode, pvnativet1, [('source_file', 'source_file'),
                                     ('pv_t1', 'in_file')]),
            (inputnode, attnativet1, [('source_file', 'source_file'),
                                     ('att_t1', 'in_file')]),
            (raw_sources, ds_asl_mask_t1, [('out', 'RawSources')]),
        ])
        if freesurfer:
            ds_asl_aseg_t1 = pe.Node(DerivativesDataSink(
                base_directory=output_dir, space='T1w', desc='aseg', suffix='dseg',
                compress=True, dismiss_entities=("echo",)),
                name='ds_asl_aseg_t1', run_without_submitting=True,
                mem_gb=DEFAULT_MEMORY_MIN_GB)
            ds_asl_aparc_t1 = pe.Node(DerivativesDataSink(
                base_directory=output_dir, space='T1w', desc='aparcaseg', suffix='dseg',
                compress=True, dismiss_entities=("echo",)),
                name='ds_asl_aparc_t1', run_without_submitting=True,
                mem_gb=DEFAULT_MEMORY_MIN_GB)
            workflow.connect([
                (inputnode, ds_asl_aseg_t1, [('source_file', 'source_file'),
                                              ('asl_aseg_t1', 'in_file')]),
                (inputnode, ds_asl_aparc_t1, [('source_file', 'source_file'),
                                               ('asl_aparc_t1', 'in_file')]),
            ])

    if getattr(spaces, '_cached') is None:
        return workflow

    # Store resamplings in standard spaces when listed in --output-spaces
    if spaces.cached.references:
        from ...niworkflows.interfaces.space import SpaceDataSource

        spacesource = pe.Node(SpaceDataSource(),
                              name='spacesource', run_without_submitting=True)
        spacesource.iterables = ('in_tuple', [
            (s.fullname, s.spec) for s in spaces.cached.get_standard(dim=(3,))
        ])

        select_std = pe.Node(KeySelect(
            fields=['template', 'asl_std', 'asl_std_ref', 'asl_mask_std',
                    'cbf_std', 'meancbf_std', 'score_std', 'avgscore_std', 'scrub_std',
                    'basil_std', 'pv_std','att_std']),
            name='select_std', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_asl_std = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir, desc='preproc', compress=True, SkullStripped=False,
                RepetitionTime=metadata.get('RepetitionTime'), TaskName=metadata.get('TaskName'),
                dismiss_entities=("echo",)),
            name='ds_asl_std', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_asl_std_ref = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='aslref', compress=True,
                                dismiss_entities=("echo",)),
            name='ds_asl_std_ref', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_asl_mask_std = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='brain', suffix='mask',
                                compress=True, dismiss_entities=("echo",)),
            name='ds_asl_mask_std', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        cbfstd = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='cbf', compress=True),
            name='cbfstd', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        meancbfstd = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='mean_cbf', compress=True),
            name='meancbfstd', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        scorestd = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='score', suffix='cbf',
                                compress=True),
            name='scorestd', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        meanscorestd = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='score', suffix='mean_cbf',
                                compress=True),
            name='meanscorestd', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        scrubstd = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='scrub', suffix='cbf',
                                compress=True),
            name='scrubstd', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        basilstd = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='basil', suffix='cbf',
                                compress=True),
            name='basilstd', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        pvstd = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='pvc', suffix='cbf',
                                compress=True),
            name='pvcstd', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        
        attstd = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='bat', suffix='cbf',
                                compress=True),
            name='attstd', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)


        workflow.connect([
            (inputnode, ds_asl_std, [('source_file', 'source_file')]),
            (inputnode, ds_asl_std_ref, [('source_file', 'source_file')]),
            (inputnode, ds_asl_mask_std, [('source_file', 'source_file')]),
            (inputnode, cbfstd, [('source_file', 'source_file')]),
            (inputnode, meancbfstd, [('source_file', 'source_file')]),
            (inputnode, scorestd, [('source_file', 'source_file')]),
            (inputnode, meanscorestd, [('source_file', 'source_file')]),
            (inputnode, scrubstd, [('source_file', 'source_file')]),
            (inputnode, basilstd, [('source_file', 'source_file')]),
            (inputnode, pvstd, [('source_file', 'source_file')]),
            (inputnode, attstd, [('source_file', 'source_file')]),
            (inputnode, select_std, [('asl_std', 'asl_std'),
                                     ('asl_std_ref', 'asl_std_ref'),
                                     ('asl_mask_std', 'asl_mask_std'),
                                     ('cbf_std', 'cbf_std'),
                                     ('meancbf_std', 'meancbf_std'),
                                     ('score_std', 'score_std'),
                                     ('avgscore_std', 'avgscore_std'),
                                     ('scrub_std', 'scrub_std'),
                                     ('basil_std', 'basil_std'),
                                     ('pv_std', 'pv_std'),
                                     ('att_std', 'att_std'),
                                     ('template', 'template'),
                                     ('spatial_reference', 'keys')]),
            (spacesource, select_std, [('uid', 'key')]),
            (select_std, ds_asl_std, [('asl_std', 'in_file')]),
            (spacesource, ds_asl_std, [('space', 'space'),
                                        ('cohort', 'cohort'),
                                        ('resolution', 'resolution'),
                                        ('density', 'density')]),
            (select_std, ds_asl_std_ref, [('asl_std_ref', 'in_file')]),
            (spacesource, ds_asl_std_ref, [('space', 'space'),
                                            ('cohort', 'cohort'),
                                            ('resolution', 'resolution'),
                                            ('density', 'density')]),
            (select_std, ds_asl_mask_std, [('asl_mask_std', 'in_file')]),
            (spacesource, ds_asl_mask_std, [('space', 'space'),
                                             ('cohort', 'cohort'),
                                             ('resolution', 'resolution'),
                                             ('density', 'density')]),
            (select_std, cbfstd, [('cbf_std', 'in_file')]),
            (spacesource, cbfstd, [('space', 'space'),
                                   ('cohort', 'cohort'),
                                   ('resolution', 'resolution'),
                                   ('density', 'density')]),
            (select_std, meancbfstd, [('meancbf_std', 'in_file')]),
            (spacesource, meancbfstd, [('space', 'space'),
                                       ('cohort', 'cohort'),
                                       ('resolution', 'resolution'),
                                       ('density', 'density')]),
            (select_std, scorestd, [('score_std', 'in_file')]),
            (spacesource, scorestd, [('space', 'space'),
                                     ('cohort', 'cohort'),
                                     ('resolution', 'resolution'),
                                     ('density', 'density')]),
            (select_std, meanscorestd, [('avgscore_std', 'in_file')]),
            (spacesource, meanscorestd, [('space', 'space'),
                                         ('cohort', 'cohort'),
                                         ('resolution', 'resolution'),
                                         ('density', 'density')]),
            (select_std, scrubstd, [('scrub_std', 'in_file')]),
            (spacesource, scrubstd, [('space', 'space'),
                                     ('cohort', 'cohort'),
                                     ('resolution', 'resolution'),
                                     ('density', 'density')]),
            (select_std, basilstd, [('basil_std', 'in_file')]),
            (spacesource, basilstd, [('space', 'space'),
                                     ('cohort', 'cohort'),
                                     ('resolution', 'resolution'),
                                     ('density', 'density')]),
            (select_std, pvstd, [('pv_std', 'in_file')]),
            (spacesource, pvstd, [('space', 'space'),
                                  ('cohort', 'cohort'),
                                  ('resolution', 'resolution'),
                                  ('density', 'density')]),
            (select_std, attstd, [('att_std', 'in_file')]),
            (spacesource, attstd, [('space', 'space'),
                                  ('cohort', 'cohort'),
                                  ('resolution', 'resolution'),
                                  ('density', 'density')]),
            (raw_sources, ds_asl_mask_std, [('out', 'RawSources')]),
        ])

        if freesurfer:
            select_fs_std = pe.Node(KeySelect(
                fields=['asl_aseg_std', 'asl_aparc_std', 'template']),
                name='select_fs_std', run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
            ds_asl_aseg_std = pe.Node(DerivativesDataSink(
                base_directory=output_dir, desc='aseg', suffix='dseg', compress=True,
                dismiss_entities=("echo",)),
                name='ds_asl_aseg_std', run_without_submitting=True,
                mem_gb=DEFAULT_MEMORY_MIN_GB)
            ds_asl_aparc_std = pe.Node(DerivativesDataSink(
                base_directory=output_dir, desc='aparcaseg', suffix='dseg', compress=True,
                dismiss_entities=("echo",)),
                name='ds_asl_aparc_std', run_without_submitting=True,
                mem_gb=DEFAULT_MEMORY_MIN_GB)
            workflow.connect([
                (spacesource, select_fs_std, [('uid', 'key')]),
                (inputnode, select_fs_std, [('asl_aseg_std', 'asl_aseg_std'),
                                            ('asl_aparc_std', 'asl_aparc_std'),
                                            ('template', 'template'),
                                            ('spatial_reference', 'keys')]),
                (select_fs_std, ds_asl_aseg_std, [('asl_aseg_std', 'in_file')]),
                (spacesource, ds_asl_aseg_std, [('space', 'space'),
                                                 ('cohort', 'cohort'),
                                                 ('resolution', 'resolution'),
                                                 ('density', 'density')]),
                (select_fs_std, ds_asl_aparc_std, [('asl_aparc_std', 'in_file')]),
                (spacesource, ds_asl_aparc_std, [('space', 'space'),
                                                  ('cohort', 'cohort'),
                                                  ('resolution', 'resolution'),
                                                  ('density', 'density')]),
                (inputnode, ds_asl_aseg_std, [('source_file', 'source_file')]),
                (inputnode, ds_asl_aparc_std, [('source_file', 'source_file')])
            ])

    fs_outputs = spaces.cached.get_fs_spaces()
    if freesurfer and fs_outputs:
        from ...niworkflows.interfaces.surf import Path2BIDS

        select_fs_surf = pe.Node(KeySelect(
            fields=['surfaces', 'surf_kwargs']), name='select_fs_surf',
            run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        select_fs_surf.iterables = [('key', fs_outputs)]
        select_fs_surf.inputs.surf_kwargs = [{'space': s} for s in fs_outputs]

        name_surfs = pe.MapNode(Path2BIDS(pattern=r'(?P<hemi>[lr])h.\w+'),
                                iterfield='in_file', name='name_surfs',
                                run_without_submitting=True)

        ds_asl_surfs = pe.MapNode(DerivativesDataSink(
            base_directory=output_dir, extension="func.gii", dismiss_entities=("echo",)),
            iterfield=['in_file', 'hemi'], name='ds_asl_surfs',
            run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)

        workflow.connect([
            (inputnode, select_fs_surf, [
                ('surf_files', 'surfaces'),
                ('surf_refs', 'keys')]),
            (select_fs_surf, name_surfs, [('surfaces', 'in_file')]),
            (inputnode, ds_asl_surfs, [('source_file', 'source_file')]),
            (select_fs_surf, ds_asl_surfs, [('surfaces', 'in_file'),
                                             ('key', 'space')]),
            (name_surfs, ds_asl_surfs, [('hemi', 'hemi')]),
        ])

    # CIFTI output
    if cifti_output:
        ds_asl_cifti = pe.Node(DerivativesDataSink(
            base_directory=output_dir, suffix='asl', compress=False, dismiss_entities=("echo",)),
            name='ds_asl_cifti', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        workflow.connect([
            (inputnode, ds_asl_cifti, [(('asl_cifti', _unlist), 'in_file'),
                                        ('source_file', 'source_file'),
                                        (('cifti_metadata', _get_surface), 'space'),
                                        ('cifti_density', 'density'),
                                        (('cifti_metadata', _read_json), 'meta_dict')])
        ])
    return workflow


def init_asl_preproc_report_wf(mem_gb, reportlets_dir, name='asl_preproc_report_wf'):
    """
    Generate a visual report.

    This workflow generates and saves a reportlet showing the effect of resampling
    the ASL signal using the standard deviation maps.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from fmriprep.workflows.asl.resampling import init_asl_preproc_report_wf
            wf = init_asl_preproc_report_wf(mem_gb=1, reportlets_dir='.')

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of ASL file in GB
    reportlets_dir : :obj:`str`
        Directory in which to save reportlets
    name : :obj:`str`, optional
        Workflow name (default: asl_preproc_report_wf)

    Inputs
    ------
    in_pre
        ASL time-series, before resampling
    in_post
        ASL time-series, after resampling
    name_source
        ASL series NIfTI file
        Used to recover original information lost during processing

    """
    from nipype.algorithms.confounds import TSNR
    from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from ...niworkflows.interfaces import SimpleBeforeAfter
    from ...interfaces import DerivativesDataSink

    workflow = Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_pre', 'in_post', 'name_source']), name='inputnode')

    pre_tsnr = pe.Node(TSNR(), name='pre_tsnr', mem_gb=mem_gb * 4.5)
    pos_tsnr = pe.Node(TSNR(), name='pos_tsnr', mem_gb=mem_gb * 4.5)

    asl_rpt = pe.Node(SimpleBeforeAfter(), name='asl_rpt',
                       mem_gb=0.1)
    ds_report_asl = pe.Node(
        DerivativesDataSink(base_directory=reportlets_dir, desc='preproc',
                            datatype="figures", dismiss_entities=("echo",)),
        name='ds_report_asl', mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True
    )

    workflow.connect([
        (inputnode, ds_report_asl, [('name_source', 'source_file')]),
        (inputnode, pre_tsnr, [('in_pre', 'in_file')]),
        (inputnode, pos_tsnr, [('in_post', 'in_file')]),
        (pre_tsnr, asl_rpt, [('stddev_file', 'before')]),
        (pos_tsnr, asl_rpt, [('stddev_file', 'after')]),
        (asl_rpt, ds_report_asl, [('out_report', 'in_file')]),
    ])

    return workflow


def _unlist(in_file):
    while isinstance(in_file, (list, tuple)) and len(in_file) == 1:
        in_file = in_file[0]
    return in_file


def _get_surface(in_file):
    from pathlib import Path
    from json import loads
    return loads(Path(in_file).read_text())["surface"]


def _read_json(in_file):
    from pathlib import Path
    from json import loads
    return loads(Path(in_file).read_text())
