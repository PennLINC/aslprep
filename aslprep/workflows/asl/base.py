# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Preprocessing workflows for ASL data."""

import numpy as np
from fmriprep.workflows.bold.apply import init_bold_volumetric_resample_wf
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from aslprep import config
from aslprep.interfaces.bids import DerivativesDataSink
from aslprep.utils.asl import determine_multi_pld, select_processing_target
from aslprep.utils.bids import collect_run_data
from aslprep.utils.misc import _create_mem_gb, _get_wf_name, get_n_volumes
from aslprep.workflows.asl.apply import init_asl_cifti_resample_wf
from aslprep.workflows.asl.cbf import init_cbf_wf, init_parcellate_cbf_wf
from aslprep.workflows.asl.confounds import (
    init_asl_confounds_wf,
    init_carpetplot_wf,
    init_cbf_confounds_wf,
)
from aslprep.workflows.asl.fit import init_asl_fit_wf, init_asl_native_wf
from aslprep.workflows.asl.outputs import (
    init_ds_asl_native_wf,
    init_ds_ciftis_wf,
    init_ds_volumes_wf,
)
from aslprep.workflows.asl.plotting import init_cbf_reporting_wf


def init_asl_wf(
    *,
    asl_file: str,
    precomputed: dict = None,
    fieldmap_id: str | None = None,
):
    """Perform the functional preprocessing stages of ASLPrep.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep import config
            from aslprep.workflows.asl.base import init_asl_wf

            with mock_config():
                asl_file = (
                    config.execution.bids_dir / "sub-01" / "perf" /
                    "sub-01_asl.nii.gz"
                )
                wf = init_asl_wf(asl_file=str(asl_file))

    Parameters
    ----------
    asl_file
        ASL NIfTI file
    precomputed
        Dictionary containing precomputed derivatives to reuse, if possible.
    fieldmap_id
        ID of the fieldmap to use to correct this BOLD series. If :obj:`None`,
        no correction will be applied.

    Inputs
    ------
    asl_file
        asl series NIfTI file
    t1w_preproc
        Bias-corrected structural template image
    t1w_mask
        Mask of the skull-stripped template image
    t1w_dseg
        Segmentation of preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF)
    t1w_tpms
        List of tissue probability maps in T1w space
    fmap_id
        Unique identifiers to select fieldmap files
    fmap
        List of estimated fieldmaps (collated with fmap_id)
    fmap_ref
        List of fieldmap reference files (collated with fmap_id)
    fmap_coeff
        List of lists of spline coefficient files (collated with fmap_id)
    fmap_mask
        List of fieldmap masks (collated with fmap_id)
    sdc_method
        List of fieldmap correction method names (collated with fmap_id)
    template
        List of templates to target
    anat2std_xfm
        List of transform files, collated with templates

    Outputs
    -------
    asl_t1
        asl series, resampled to T1w space
    asl_mask_t1
        asl series mask in T1w space
    asl_std
        asl series, resampled to template space
    asl_mask_std
        asl series mask in template space
    confounds
        TSV of confounds
    cbf_ts_t1
        cbf times series in T1w space
    mean_cbf_t1
        mean cbf   in T1w space
    cbf_ts_score_t1
        scorecbf times series in T1w space
    mean_cbf_score_t1
        mean score cbf  in T1w space
    mean_cbf_scrub_t1, mean_cbf_gm_basil_t1, mean_cbf_basil_t1
        scrub, partial volume corrected and basil cbf in T1w space
    cbf_ts_std
        cbf times series in template space
    mean_cbf_std
        mean cbf   in template space
    cbf_ts_score_std
        scorecbf times series in template space
    mean_cbf_score_std
        mean score cbf  in template space
    mean_cbf_scrub_std, mean_cbf_gm_basil_std, mean_cbf_basil_std
        scrub, partial volume corrected and basil cbf in template space
    qc_file
        quality control measures

    Notes
    -----
    Note that ``anat2std_xfm``, ``std_space``, ``std_resolution``,
    ``std_cohort``, ``std_t1w`` and ``std_mask`` are treated as single
    inputs. In order to resample to multiple target spaces, connect
    these fields to an iterable.

    1.  Brain-mask T1w.
    2.  Generate initial ASL reference image.
        I think this is just the BOLD reference workflow from niworkflows,
        without any modifications for ASL data. I could be missing something.
        -   Not in GE workflow.
        -   In the GE workflow, the reference image comes from the M0 scan.
    3.  Motion correction.
        -   Outputs the HMC transform and the motion parameters.
        -   Not in GE workflow.
    4.  Susceptibility distortion correction.
        -   Outputs the SDC transform.
        -   Not in GE workflow.
    5.  Register ASL to T1w.
    6.  Apply the ASL-to-T1w transforms to get T1w-space outputs
        (passed along to derivatives workflow).
    7.  Apply the ASL-to-ASLref transforms to get native-space outputs.
        -   Not in GE workflow.
    8.  Calculate confounds.
        -   Not in GE workflow.
    9.  Calculate CBF.
    10. Refine the brain mask.
    11. Generate distortion correction report.
        -   Not in GE workflow.
    12. Apply ASL-to-template transforms to get template-space outputs.
    13. CBF QC workflow.
    14. CBF plotting workflow.
    15. Generate carpet plots.
        -   Not in GE workflow.
    16. Parcellate CBF results.
    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    # Have some options handy
    omp_nthreads = config.nipype.omp_nthreads
    spaces = config.workflow.spaces
    dummy_scans = config.workflow.dummy_scans
    smooth_kernel = config.workflow.smooth_kernel
    m0_scale = config.workflow.m0_scale
    scorescrub = config.workflow.scorescrub
    basil = config.workflow.basil
    nonstd_spaces = set(spaces.get_nonstandard())
    freesurfer_spaces = spaces.get_fs_spaces()
    layout = config.execution.layout

    # If number of ASL volumes is less than 5, motion correction, etc. will be skipped.
    n_vols = get_n_volumes(asl_file)
    use_ge = config.workflow.use_ge if isinstance(config.workflow.use_ge, bool) else n_vols <= 5
    if use_ge:
        config.loggers.workflow.warning('Using GE-specific processing. HMC will be disabled.')
        if scorescrub:
            config.loggers.workflow.warning(
                'SCORE/SCRUB processing is disabled for GE-specific processing'
            )
            scorescrub = False

    asl_tlen, mem_gb = _create_mem_gb(asl_file)

    wf_name = _get_wf_name(asl_file)
    config.loggers.workflow.debug(
        (
            'Creating asl processing workflow for "%s" (%.2f GB / %d TRs). '
            'Memory resampled/largemem=%.2f/%.2f GB.'
        ),
        asl_file,
        mem_gb['filesize'],
        asl_tlen,
        mem_gb['resampled'],
        mem_gb['largemem'],
    )

    # Collect associated files
    run_data = collect_run_data(layout, asl_file)
    metadata = run_data['asl_metadata'].copy()
    # Patch RepetitionTimePreparation into RepetitionTime,
    # for the sake of BOLD-based interfaces and workflows.
    # This value shouldn't be used for anything except figures and reportlets.
    metadata['RepetitionTime'] = metadata.get(
        'RepetitionTime',
        np.mean(metadata['RepetitionTimePreparation']),
    )

    is_multi_pld = determine_multi_pld(metadata=metadata)
    if scorescrub and is_multi_pld:
        config.loggers.workflow.warning(
            f'SCORE/SCRUB processing will be disabled for multi-delay {asl_file}'
        )
        scorescrub = False

    # Determine which volumes to use in the pipeline
    processing_target = select_processing_target(aslcontext=run_data['aslcontext'])

    # Determine which CBF outputs to expect
    att_derivs = []
    cbf_3d_derivs = ['mean_cbf']
    cbf_4d_derivs = []

    if is_multi_pld:
        att_derivs += ['att']
    else:
        cbf_4d_derivs += ['cbf_ts']

    if scorescrub:
        cbf_4d_derivs += ['cbf_ts_score']
        cbf_3d_derivs += [
            'mean_cbf_score',
            'mean_cbf_scrub',
        ]

    if basil:
        cbf_3d_derivs += [
            'mean_cbf_basil',
            'mean_cbf_gm_basil',
            'mean_cbf_wm_basil',
        ]
        att_derivs += ['att_basil']

    cbf_derivs = att_derivs + cbf_3d_derivs + cbf_4d_derivs

    # Build workflow
    workflow = Workflow(name=wf_name)
    workflow.__postdesc__ = """\
All resampling in *ASLPrep* uses a single interpolation step that concatenates all transformations.
Gridded (volumetric) resampling was performed using `antsApplyTransforms`,
configured with *Lanczos* interpolation to minimize the smoothing effects of other kernels
[@lanczos].
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                # ASL-specific elements
                'asl_file',
                'asl_metadata',
                'aslcontext',
                'm0scan',
                'm0scan_metadata',
                # Anatomical coregistration
                't1w_preproc',
                't1w_mask',
                't1w_dseg',
                't1w_tpms',
                # FreeSurfer outputs
                'subjects_dir',
                'subject_id',
                'fsnative2t1w_xfm',
                'white',
                'midthickness',
                'pial',
                'sphere_reg_fsLR',
                'midthickness_fsLR',
                'cortex_mask',
                'anat_ribbon',
                # Fieldmap registration
                'fmap',
                'fmap_ref',
                'fmap_coeff',
                'fmap_mask',
                'fmap_id',
                'sdc_method',
                # Volumetric templates
                'anat2std_xfm',
                'std_t1w',
                'std_mask',
                'std_space',
                'std_resolution',
                'std_cohort',
                # MNI152NLin6Asym warp, for CIFTI use
                'anat2mni6_xfm',
                'mni6_mask',
                # MNI152NLin2009cAsym inverse warp, for carpetplotting and CBF QC
                'mni2009c2anat_xfm',
            ],
        ),
        name='inputnode',
    )
    inputnode.inputs.asl_file = asl_file
    inputnode.inputs.asl_metadata = metadata
    inputnode.inputs.aslcontext = run_data['aslcontext']
    inputnode.inputs.m0scan_metadata = run_data['m0scan_metadata']

    # Perform minimal preprocessing of the ASL data, including HMC and SDC
    asl_fit_wf = init_asl_fit_wf(
        asl_file=asl_file,
        m0scan=run_data['m0scan'],
        use_ge=use_ge,
        precomputed=precomputed,
        fieldmap_id=fieldmap_id,
        omp_nthreads=omp_nthreads,
    )

    workflow.connect([
        (inputnode, asl_fit_wf, [
            ('aslcontext', 'inputnode.aslcontext'),
            # Original inputs from fMRIPrep
            ('t1w_preproc', 'inputnode.t1w_preproc'),
            ('t1w_mask', 'inputnode.t1w_mask'),
            ('t1w_dseg', 'inputnode.t1w_dseg'),
            ('subjects_dir', 'inputnode.subjects_dir'),
            ('subject_id', 'inputnode.subject_id'),
            ('fsnative2t1w_xfm', 'inputnode.fsnative2t1w_xfm'),
            ('fmap', 'inputnode.fmap'),
            ('fmap_ref', 'inputnode.fmap_ref'),
            ('fmap_coeff', 'inputnode.fmap_coeff'),
            ('fmap_mask', 'inputnode.fmap_mask'),
            ('fmap_id', 'inputnode.fmap_id'),
            ('sdc_method', 'inputnode.sdc_method'),
        ]),
    ])  # fmt:skip

    # Resample to aslref space.
    # NOTE: This differs from fMRIPrep, which puts this step at the 'resampling' level,
    # because ASLPrep's main output is CBF and we need aslref-space data to calculate CBF.
    asl_native_wf = init_asl_native_wf(
        asl_file=asl_file,
        m0scan=run_data['m0scan'],
        fieldmap_id=fieldmap_id,
        omp_nthreads=omp_nthreads,
        name='asl_native_wf',
    )

    workflow.connect([
        (inputnode, asl_native_wf, [
            ('aslcontext', 'inputnode.aslcontext'),
            ('fmap_ref', 'inputnode.fmap_ref'),
            ('fmap_coeff', 'inputnode.fmap_coeff'),
            ('fmap_id', 'inputnode.fmap_id'),
        ]),
        (asl_fit_wf, asl_native_wf, [
            ('outputnode.coreg_aslref', 'inputnode.aslref'),
            ('outputnode.asl_mask', 'inputnode.asl_mask'),
            ('outputnode.motion_xfm', 'inputnode.motion_xfm'),
            ('outputnode.aslref2fmap_xfm', 'inputnode.aslref2fmap_xfm'),
            ('outputnode.dummy_scans', 'inputnode.dummy_scans'),
        ]),
    ])  # fmt:skip

    # Calculate CBF
    # Compute CBF from ASLRef-space ASL data.
    cbf_wf = init_cbf_wf(
        name_source=asl_file,
        processing_target=processing_target,
        dummy_scans=dummy_scans,
        m0_scale=m0_scale,
        scorescrub=scorescrub,
        basil=basil,
        smooth_kernel=smooth_kernel,
        metadata=metadata,
        name='cbf_wf',
    )
    workflow.connect([
        (inputnode, cbf_wf, [
            ('t1w_mask', 'inputnode.t1w_mask'),
            ('t1w_tpms', 'inputnode.t1w_tpms'),
            ('m0scan_metadata', 'inputnode.m0scan_metadata'),
        ]),
        (asl_fit_wf, cbf_wf, [
            ('outputnode.asl_mask', 'inputnode.asl_mask'),
            ('outputnode.aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
        ]),
        # Use post-HMC+SDC ASL as basis for CBF calculation
        (asl_native_wf, cbf_wf, [
            ('outputnode.asl_native', 'inputnode.asl_file'),
            ('outputnode.aslcontext', 'inputnode.aslcontext'),
            ('outputnode.m0scan_native', 'inputnode.m0scan'),
        ]),
    ])  # fmt:skip

    if config.workflow.level == 'minimal':
        return workflow

    #
    # Resampling outputs workflow:
    #   - Calculate ASL confounds and CBF QC metrics.
    #   - Generate plots for CBF.
    #   - Resample to anatomical space.
    #   - Save aslref-space outputs only if requested.
    #

    asl_confounds_wf = init_asl_confounds_wf(
        n_volumes=n_vols,
        mem_gb=mem_gb['largemem'],
        freesurfer=config.workflow.run_reconall,
        name='asl_confounds_wf',
    )

    ds_confounds = pe.Node(
        DerivativesDataSink(
            base_directory=config.execution.aslprep_dir,
            desc='confounds',
            suffix='timeseries',
            dismiss_entities=('echo',),
        ),
        name='ds_confounds',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    ds_confounds.inputs.source_file = asl_file

    workflow.connect([
        (inputnode, asl_confounds_wf, [
            ('t1w_tpms', 'inputnode.t1w_tpms'),
            ('t1w_mask', 'inputnode.t1w_mask'),
        ]),
        (asl_fit_wf, asl_confounds_wf, [
            ('outputnode.asl_mask', 'inputnode.asl_mask'),
            ('outputnode.movpar_file', 'inputnode.movpar_file'),
            ('outputnode.rmsd_file', 'inputnode.rmsd_file'),
            ('outputnode.aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
            ('outputnode.dummy_scans', 'inputnode.skip_vols'),
        ]),
        (asl_native_wf, asl_confounds_wf, [('outputnode.asl_native', 'inputnode.asl')]),
        (asl_confounds_wf, ds_confounds, [
            ('outputnode.confounds_file', 'in_file'),
            ('outputnode.confounds_metadata', 'meta_dict'),
        ]),
    ])  # fmt:skip

    # Generate QC metrics
    cbf_confounds_wf = init_cbf_confounds_wf(
        scorescrub=scorescrub,
        basil=basil,
        name='cbf_confounds_wf',
    )
    workflow.connect([
        (inputnode, cbf_confounds_wf, [
            ('asl_file', 'inputnode.name_source'),
            ('t1w_tpms', 'inputnode.t1w_tpms'),
            ('t1w_mask', 'inputnode.t1w_mask'),
            ('mni2009c2anat_xfm', 'inputnode.mni2009c2anat_xfm'),
        ]),
        (asl_fit_wf, cbf_confounds_wf, [
            ('outputnode.asl_mask', 'inputnode.asl_mask'),
            ('outputnode.aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
            ('outputnode.rmsd_file', 'inputnode.rmsd_file'),
        ]),
        (asl_confounds_wf, cbf_confounds_wf, [
            ('outputnode.confounds_file', 'inputnode.confounds_file'),
        ]),
    ])  # fmt:skip
    for cbf_deriv in cbf_3d_derivs:
        workflow.connect([
            (cbf_wf, cbf_confounds_wf, [(f'outputnode.{cbf_deriv}', f'inputnode.{cbf_deriv}')]),
        ])  # fmt:skip

    # Plot CBF outputs.
    # NOTE: CIFTI input won't be provided unless level is set to 'full'.
    cbf_reporting_wf = init_cbf_reporting_wf(
        metadata=metadata,
        plot_timeseries=not (is_multi_pld or use_ge or (config.workflow.level == 'resampling')),
        scorescrub=scorescrub,
        basil=basil,
        name='cbf_reporting_wf',
    )
    workflow.connect([
        (inputnode, cbf_reporting_wf, [
            ('t1w_dseg', 'inputnode.t1w_dseg'),
            ('mni2009c2anat_xfm', 'inputnode.std2anat_xfm'),
        ]),
        (asl_confounds_wf, cbf_reporting_wf, [
            ('outputnode.confounds_file', 'inputnode.confounds_file'),
        ]),
        (cbf_wf, cbf_reporting_wf, [
            ('outputnode.score_outlier_index', 'inputnode.score_outlier_index'),
        ]),
        (cbf_confounds_wf, cbf_reporting_wf, [('outputnode.qc_file', 'inputnode.qc_file')]),
        (asl_fit_wf, cbf_reporting_wf, [
            ('outputnode.coreg_aslref', 'inputnode.aslref'),
            ('outputnode.aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
            # XXX: Used to use the one from refine_mask/reduce_mask
            ('outputnode.asl_mask', 'inputnode.asl_mask'),
        ]),
        (asl_confounds_wf, cbf_reporting_wf, [
            ('outputnode.crown_mask', 'inputnode.crown_mask'),
            ('outputnode.acompcor_masks', 'inputnode.acompcor_masks'),
        ]),
    ])  # fmt:skip

    for cbf_deriv in cbf_derivs:
        workflow.connect([
            (cbf_wf, cbf_reporting_wf, [(f'outputnode.{cbf_deriv}', f'inputnode.{cbf_deriv}')]),
        ])  # fmt:skip

    # If we want aslref-space outputs, then call the appropriate workflow
    aslref_out = bool(nonstd_spaces.intersection(('func', 'run', 'asl', 'aslref', 'sbref')))
    aslref_out &= config.workflow.level == 'full'

    if aslref_out:
        ds_asl_native_wf = init_ds_asl_native_wf(
            bids_root=str(config.execution.bids_dir),
            output_dir=config.execution.aslprep_dir,
            asl_output=aslref_out,
            metadata=metadata,
            cbf_3d=cbf_3d_derivs,
            cbf_4d=cbf_4d_derivs,
            att=att_derivs,
        )
        ds_asl_native_wf.inputs.inputnode.source_files = [asl_file]

        workflow.connect([
            (asl_fit_wf, ds_asl_native_wf, [('outputnode.asl_mask', 'inputnode.asl_mask')]),
            (asl_native_wf, ds_asl_native_wf, [('outputnode.asl_native', 'inputnode.asl')]),
        ])  # fmt:skip

        for cbf_deriv in cbf_derivs:
            workflow.connect([
                (cbf_wf, ds_asl_native_wf, [
                    (f'outputnode.{cbf_deriv}', f'inputnode.{cbf_deriv}'),
                ]),
            ])  # fmt:skip

    if config.workflow.level == 'resampling':
        # Fill in datasinks of reportlets seen so far
        for node in workflow.list_node_names():
            if node.split('.')[-1].startswith('ds_report'):
                workflow.get_node(node).inputs.base_directory = config.execution.aslprep_dir
                workflow.get_node(node).inputs.source_file = asl_file

        return workflow

    # Resample ASL file to anatomical space.
    # This doesn't write out the resampled file to the derivatives.
    asl_anat_wf = init_bold_volumetric_resample_wf(
        metadata=metadata,
        fieldmap_id=fieldmap_id,
        omp_nthreads=omp_nthreads,
        mem_gb=mem_gb,
        jacobian='fmap-jacobian' not in config.workflow.ignore,
        name='asl_anat_wf',
    )
    asl_anat_wf.inputs.inputnode.resolution = 'native'

    workflow.connect([
        (inputnode, asl_anat_wf, [
            ('t1w_preproc', 'inputnode.target_ref_file'),
            ('t1w_mask', 'inputnode.target_mask'),
            ('fmap_ref', 'inputnode.fmap_ref'),
            ('fmap_coeff', 'inputnode.fmap_coeff'),
            ('fmap_id', 'inputnode.fmap_id'),
        ]),
        (asl_fit_wf, asl_anat_wf, [
            ('outputnode.coreg_aslref', 'inputnode.bold_ref_file'),
            ('outputnode.aslref2fmap_xfm', 'inputnode.boldref2fmap_xfm'),
            ('outputnode.aslref2anat_xfm', 'inputnode.boldref2anat_xfm'),
        ]),
        (asl_native_wf, asl_anat_wf, [
            ('outputnode.asl_minimal', 'inputnode.bold_file'),
            ('outputnode.motion_xfm', 'inputnode.motion_xfm'),
        ]),
    ])  # fmt:skip

    # Write out anatomical-space derivatives.
    if nonstd_spaces.intersection(('anat', 'T1w')):
        ds_asl_t1_wf = init_ds_volumes_wf(
            bids_root=str(config.execution.bids_dir),
            output_dir=config.execution.aslprep_dir,
            metadata=metadata,
            cbf_3d=cbf_3d_derivs,
            cbf_4d=cbf_4d_derivs,
            att=att_derivs,
            name='ds_asl_t1_wf',
        )
        ds_asl_t1_wf.inputs.inputnode.source_files = [asl_file]
        ds_asl_t1_wf.inputs.inputnode.space = 'T1w'

        workflow.connect([
            (inputnode, ds_asl_t1_wf, [('t1w_preproc', 'inputnode.ref_file')]),
            (asl_fit_wf, ds_asl_t1_wf, [
                ('outputnode.asl_mask', 'inputnode.asl_mask'),
                ('outputnode.coreg_aslref', 'inputnode.aslref'),
                ('outputnode.aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
            ]),
            (asl_anat_wf, ds_asl_t1_wf, [('outputnode.bold_file', 'inputnode.asl')]),
        ])  # fmt:skip

        for cbf_deriv in cbf_derivs:
            workflow.connect([
                (cbf_wf, ds_asl_t1_wf, [(f'outputnode.{cbf_deriv}', f'inputnode.{cbf_deriv}')]),
            ])  # fmt:skip

    # Resample derivatives to standard space and write them out.
    # This is different from having internally-produced standard-space data.
    if spaces.cached.get_spaces(nonstandard=False, dim=(3,)):
        # Missing:
        #  * Clipping BOLD after resampling
        #  * Resampling parcellations
        asl_std_wf = init_bold_volumetric_resample_wf(
            metadata=metadata,
            fieldmap_id=fieldmap_id,
            omp_nthreads=omp_nthreads,
            mem_gb=mem_gb,
            jacobian='fmap-jacobian' not in config.workflow.ignore,
            name='asl_std_wf',
        )
        ds_asl_std_wf = init_ds_volumes_wf(
            bids_root=str(config.execution.bids_dir),
            output_dir=config.execution.aslprep_dir,
            metadata=metadata,
            cbf_3d=cbf_3d_derivs,
            cbf_4d=cbf_4d_derivs,
            att=att_derivs,
            name='ds_asl_std_wf',
        )
        ds_asl_std_wf.inputs.inputnode.source_files = [asl_file]

        workflow.connect([
            (inputnode, asl_std_wf, [
                ('std_t1w', 'inputnode.target_ref_file'),
                ('std_mask', 'inputnode.target_mask'),
                ('anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('std_resolution', 'inputnode.resolution'),
                ('fmap_ref', 'inputnode.fmap_ref'),
                ('fmap_coeff', 'inputnode.fmap_coeff'),
                ('fmap_id', 'inputnode.fmap_id'),
            ]),
            (asl_fit_wf, asl_std_wf, [
                ('outputnode.coreg_aslref', 'inputnode.bold_ref_file'),
                ('outputnode.aslref2fmap_xfm', 'inputnode.boldref2fmap_xfm'),
                ('outputnode.aslref2anat_xfm', 'inputnode.boldref2anat_xfm'),
            ]),
            (asl_native_wf, asl_std_wf, [
                ('outputnode.asl_minimal', 'inputnode.bold_file'),
                ('outputnode.motion_xfm', 'inputnode.motion_xfm'),
            ]),
            (inputnode, ds_asl_std_wf, [
                ('anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('std_space', 'inputnode.space'),
                ('std_resolution', 'inputnode.resolution'),
                ('std_cohort', 'inputnode.cohort'),
            ]),
            (asl_fit_wf, ds_asl_std_wf, [
                ('outputnode.asl_mask', 'inputnode.asl_mask'),
                ('outputnode.coreg_aslref', 'inputnode.aslref'),
                ('outputnode.aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
            ]),
            (asl_std_wf, ds_asl_std_wf, [
                ('outputnode.bold_file', 'inputnode.asl'),
                ('outputnode.resampling_reference', 'inputnode.ref_file'),
            ]),
        ])  # fmt:skip

        for cbf_deriv in cbf_derivs:
            workflow.connect([
                (cbf_wf, ds_asl_std_wf, [(f'outputnode.{cbf_deriv}', f'inputnode.{cbf_deriv}')]),
            ])  # fmt:skip

    # GIFTI outputs
    if config.workflow.run_reconall and freesurfer_spaces:
        from aslprep.workflows.asl.resampling import init_asl_surf_wf

        workflow.__postdesc__ += """\
Non-gridded (surface) resamplings were performed using `mri_vol2surf` (FreeSurfer).
"""
        config.loggers.workflow.debug('Creating ASL surface-sampling workflow.')

        # init_bold_surf_wf uses prepare_timing_parameters, which uses the config object.
        # The uninitialized fMRIPrep config will have config.workflow.ignore set to None
        # instead of a list, which will raise an error.
        asl_surf_wf = init_asl_surf_wf(
            mem_gb=mem_gb['resampled'],
            metadata=metadata,
            surface_spaces=freesurfer_spaces,
            medial_surface_nan=config.workflow.medial_surface_nan,
            output_dir=config.execution.aslprep_dir,
            cbf_3d=cbf_3d_derivs,
            cbf_4d=cbf_4d_derivs,
            att=att_derivs,
            name='asl_surf_wf',
        )
        asl_surf_wf.inputs.inputnode.source_file = asl_file
        workflow.connect([
            (inputnode, asl_surf_wf, [
                ('t1w_preproc', 'inputnode.anat'),
                ('subjects_dir', 'inputnode.subjects_dir'),
                ('subject_id', 'inputnode.subject_id'),
                ('fsnative2t1w_xfm', 'inputnode.fsnative2t1w_xfm'),
            ]),
            (asl_fit_wf, asl_surf_wf, [
                ('outputnode.aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
            ]),
        ])  # fmt:skip
        for cbf_deriv in cbf_derivs:
            workflow.connect([
                (cbf_wf, asl_surf_wf, [(f'outputnode.{cbf_deriv}', f'inputnode.{cbf_deriv}')]),
            ])  # fmt:skip

    if config.workflow.cifti_output:
        asl_cifti_resample_wf = init_asl_cifti_resample_wf(
            metadata=metadata,
            mem_gb=mem_gb,
            fieldmap_id=fieldmap_id,
            omp_nthreads=omp_nthreads,
        )

        workflow.connect([
            (inputnode, asl_cifti_resample_wf, [
                ('mni6_mask', 'inputnode.mni6_mask'),
                ('anat2mni6_xfm', 'inputnode.anat2mni6_xfm'),
                ('fmap_ref', 'inputnode.fmap_ref'),
                ('fmap_coeff', 'inputnode.fmap_coeff'),
                ('fmap_id', 'inputnode.fmap_id'),
                ('white', 'inputnode.white'),
                ('pial', 'inputnode.pial'),
                ('midthickness', 'inputnode.midthickness'),
                ('midthickness_fsLR', 'inputnode.midthickness_fsLR'),
                ('sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR'),
                ('cortex_mask', 'inputnode.cortex_mask'),
                ('anat_ribbon', 'inputnode.anat_ribbon'),
            ]),
            (asl_fit_wf, asl_cifti_resample_wf, [
                ('outputnode.coreg_aslref', 'inputnode.coreg_aslref'),
                ('outputnode.aslref2fmap_xfm', 'inputnode.aslref2fmap_xfm'),
                ('outputnode.aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
            ]),
            (asl_native_wf, asl_cifti_resample_wf, [
                ('outputnode.asl_minimal', 'inputnode.asl_file'),
                ('outputnode.motion_xfm', 'inputnode.motion_xfm'),
            ]),
            (asl_anat_wf, asl_cifti_resample_wf, [('outputnode.bold_file', 'inputnode.asl_anat')]),
        ])  # fmt:skip

        ds_asl_cifti_wf = init_ds_ciftis_wf(
            bids_root=str(config.execution.bids_dir),
            output_dir=config.execution.aslprep_dir,
            metadata=metadata,
            cbf_3d=cbf_3d_derivs,
            cbf_4d=cbf_4d_derivs,
            att=att_derivs,
            omp_nthreads=omp_nthreads,
            name='ds_asl_cifti_wf',
        )
        ds_asl_cifti_wf.inputs.inputnode.source_files = [asl_file]
        workflow.connect([
            (inputnode, ds_asl_cifti_wf, [
                ('t1w_preproc', 'inputnode.anat'),
                ('mni6_mask', 'inputnode.mni6_mask'),
                ('anat2mni6_xfm', 'inputnode.anat2mni6_xfm'),
                ('white', 'inputnode.white'),
                ('pial', 'inputnode.pial'),
                ('midthickness', 'inputnode.midthickness'),
                ('midthickness_fsLR', 'inputnode.midthickness_fsLR'),
                ('sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR'),
                ('cortex_mask', 'inputnode.cortex_mask'),
            ]),
            (asl_fit_wf, ds_asl_cifti_wf, [
                ('outputnode.aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
            ]),
            (asl_cifti_resample_wf, ds_asl_cifti_wf, [
                ('outputnode.asl_cifti', 'inputnode.asl_cifti'),
                ('outputnode.goodvoxels_mask', 'inputnode.goodvoxels_mask'),
            ]),
        ])  # fmt:skip
        for cbf_deriv in cbf_derivs:
            workflow.connect([
                (cbf_wf, ds_asl_cifti_wf, [(f'outputnode.{cbf_deriv}', f'inputnode.{cbf_deriv}')]),
            ])  # fmt:skip

        # Feed CIFTI into CBF-reporting workflow
        if 'cbf_ts' in cbf_4d_derivs:
            workflow.connect([
                (ds_asl_cifti_wf, cbf_reporting_wf, [
                    ('outputnode.cbf_ts', 'inputnode.cifti_cbf_ts'),
                ]),
            ])  # fmt:skip

    # Should always be reached, since ASLPrep includes MNI152NLin2009cAsym automatically
    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # Don't make carpet plots for short or GE data
        if not use_ge:
            carpetplot_wf = init_carpetplot_wf(
                mem_gb=mem_gb['resampled'],
                confounds_list=[
                    ('global_signal', None, 'GS'),
                    ('csf', None, 'GSCSF'),
                    ('white_matter', None, 'GSWM'),
                    ('std_dvars', None, 'DVARS'),
                    ('framewise_displacement', 'mm', 'FD'),
                ],
                metadata=metadata,
                cifti_output=config.workflow.cifti_output,
                suffix='asl',
                name='carpetplot_wf',
            )

            if config.workflow.cifti_output:
                workflow.connect([
                    (asl_cifti_resample_wf, carpetplot_wf, [
                        ('outputnode.asl_cifti', 'inputnode.cifti_asl'),
                    ]),
                ])  # fmt:skip

            def _last(inlist):
                if not isinstance(inlist, list):
                    raise ValueError(f'_last: input is not list ({inlist})')

                return inlist[-1]

            workflow.connect([
                (inputnode, carpetplot_wf, [('mni2009c2anat_xfm', 'inputnode.std2anat_xfm')]),
                (asl_fit_wf, carpetplot_wf, [
                    ('outputnode.dummy_scans', 'inputnode.dummy_scans'),
                    ('outputnode.asl_mask', 'inputnode.asl_mask'),
                    ('outputnode.aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
                ]),
                (asl_native_wf, carpetplot_wf, [('outputnode.asl_native', 'inputnode.asl')]),
                (asl_confounds_wf, carpetplot_wf, [
                    ('outputnode.confounds_file', 'inputnode.confounds_file'),
                    ('outputnode.crown_mask', 'inputnode.crown_mask'),
                    (('outputnode.acompcor_masks', _last), 'inputnode.acompcor_mask'),
                ]),
            ])  # fmt:skip

        # Parcellate CBF maps and write out parcellated TSV files and atlases
        parcellate_cbf_wf = init_parcellate_cbf_wf(
            cbf_3d=cbf_3d_derivs,
            name='parcellate_cbf_wf',
        )
        workflow.connect([
            (inputnode, parcellate_cbf_wf, [
                ('asl_file', 'inputnode.source_file'),
                ('mni2009c2anat_xfm', 'inputnode.MNI152NLin2009cAsym_to_anat_xfm'),
            ]),
            (asl_fit_wf, parcellate_cbf_wf, [
                ('outputnode.asl_mask', 'inputnode.asl_mask'),
                ('outputnode.aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
            ]),
        ])  # fmt:skip

        for cbf_deriv in cbf_3d_derivs:
            workflow.connect([
                (cbf_wf, parcellate_cbf_wf, [
                    (f'outputnode.{cbf_deriv}', f'inputnode.{cbf_deriv}'),
                ]),
            ])  # fmt:skip

    # Fill in datasinks of reportlets seen so far
    for node in workflow.list_node_names():
        if node.split('.')[-1].startswith('ds_report'):
            workflow.get_node(node).inputs.base_directory = config.execution.aslprep_dir
            workflow.get_node(node).inputs.source_file = asl_file

    return workflow


def _read_json(in_file):
    from json import loads
    from pathlib import Path

    if not isinstance(in_file, str):
        raise ValueError(f'_read_json: input is not str ({in_file})')

    return loads(Path(in_file).read_text())
