# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Orchestrating the ASL-preprocessing workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_asl_preproc_wf
.. autofunction:: init_asl_derivatives_wf

"""
from ... import config

import os

import nibabel as nb
from nipype.interfaces.fsl import Split as FSLSplit
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from ...utils.meepi import combine_meepi_source

from ...interfaces import DerivativesDataSink
from ...interfaces.reports import FunctionalSummary

# asl workflows
from .confounds import init_asl_confs_wf, init_carpetplot_wf
from .hmc import init_asl_hmc_wf
from .stc import init_asl_stc_wf
from .t2s import init_asl_t2s_wf
from .registration import init_asl_t1_trans_wf, init_asl_reg_wf
from .resampling import (
    init_asl_surf_wf,
    init_asl_std_trans_wf,
    init_asl_preproc_trans_wf,
)

# cbf workflows
from .cbf import (
    init_cbf_compt_wf,
    init_cbfqc_compt_wf,
    init_cbfplot_wf,
    init_cbfroiquant_wf)
from .outputs import init_asl_derivatives_wf


from ...interfaces.cbf_computation import refinemask


def init_asl_preproc_wf(asl_file):
    """
    This workflow controls the functional preprocessing stages of *aslprep*.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.tests import mock_config
            from aslprep import config
            from aslprep.workflows.asl.base import init_asl_preproc_wf
            with mock_config():
                asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_task-restEyesOpen_asl.nii.gz'
                wf = init_asl_preproc_wf(str(asl_file))

    Parameters
    ----------
    asl_file
        asl series NIfTI file

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
    template
        List of templates to target
    anat2std_xfm
        List of transform files, collated with templates
    std2anat_xfm
        List of inverse transform files, collated with templates
    

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
    cbf_t1
        cbf times series in T1w space
    meancbf_t1
        mean cbf   in T1w space
    scorecbf_t1
        scorecbf times series in T1w space
    avgscorecbf_t1
        mean score cbf  in T1w space
    scrub_t1, pv_t1, basil_t1
        scrub, parital volume corrected and basil cbf   in T1w space
    cbf_std
        cbf times series in template space
    meancbf_std
        mean cbf   in template space
    scorecbf_std
        scorecbf times series in template space
    avgscorecbf_std
        mean score cbf  in template space
    scrub_std, pv_std, basil_std
        scrub, parital volume corrected and basil cbf   in template space
    qc_file
        quality control meausres 

    See Also
    --------

    * :py:func:`~aslprep.niworkflows.func.util.init_asl_reference_wf`
    * :py:func:`~aslprep.workflows.asl.stc.init_asl_stc_wf`
    * :py:func:`~aslprep.workflows.asl.hmc.init_asl_hmc_wf`
    * :py:func:`~aslprep.workflows.asl.t2s.init_asl_t2s_wf`
    * :py:func:`~aslprep.workflows.asl.registration.init_asl_t1_trans_wf`
    * :py:func:`~aslprep.workflows.asl.registration.init_asl_reg_wf`
    * :py:func:`~aslprep.workflows.asl.confounds.init_asl_confounds_wf`
    * :py:func:`~aslprep.workflows.asl.confounds.init_ica_aroma_wf`
    * :py:func:`~aslprep.workflows.asl.resampling.init_asl_std_trans_wf`
    * :py:func:`~aslprep.workflows.asl.resampling.init_asl_preproc_trans_wf`
    * :py:func:`~aslprep.workflows.asl.resampling.init_asl_surf_wf`
    * :py:func:`~sdcflows.workflows.fmap.init_fmap_wf`
    * :py:func:`~sdcflows.workflows.pepolar.init_pepolar_unwarp_wf`
    * :py:func:`~sdcflows.workflows.phdiff.init_phdiff_wf`
    * :py:func:`~sdcflows.workflows.syn.init_syn_sdc_wf`
    * :py:func:`~sdcflows.workflows.unwarp.init_sdc_unwarp_wf`

    """
    from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from ...niworkflows.func.util import init_asl_reference_wf
    from ...niworkflows.interfaces.nibabel import ApplyMask
    from ...niworkflows.interfaces.utility import KeySelect
    from sdcflows.workflows.base import init_sdc_estimate_wf, fieldmap_wrangler

    ref_file = asl_file
    mem_gb = {'filesize': 1, 'resampled': 1, 'largemem': 1}
    asl_tlen = 10
    multiecho = isinstance(asl_file, list)

    # Have some options handy
    layout = config.execution.layout
    omp_nthreads = config.nipype.omp_nthreads
    spaces = config.workflow.spaces
    output_dir = str(config.execution.output_dir)
    dummyvols = config.workflow.dummy_vols
    smoothkernel = config.workflow.smooth_kernel
    mscale = config.workflow.m0_scale
    scorescrub = config.workflow.scorescrub
    basil = config.workflow.basil
    

    if multiecho:
        tes = [layout.get_metadata(echo)['EchoTime'] for echo in asl_file]
        ref_file = dict(zip(tes, asl_file))[min(tes)]

    if os.path.isfile(ref_file):
        asl_tlen, mem_gb = _create_mem_gb(ref_file)

    wf_name = _get_wf_name(ref_file)
    config.loggers.workflow.debug(
        'Creating asl processing workflow for "%s" (%.2f GB / %d TRs). '
        'Memory resampled/largemem=%.2f/%.2f GB.',
        ref_file, mem_gb['filesize'], asl_tlen, mem_gb['resampled'], mem_gb['largemem'])

    sbref_file = None
    # Find associated sbref, if possible
    entities = layout.parse_file_entities(ref_file)
    entities['suffix'] = 'sbref'
    entities['extension'] = ['nii', 'nii.gz']  # Overwrite extensions
    files = layout.get(return_type='file', **entities)
    refbase = os.path.basename(ref_file)
    if 'sbref' in config.workflow.ignore:
        config.loggers.workflow.info("Single-band reference files ignored.")
    elif files and multiecho:
        config.loggers.workflow.warning(
            "Single-band reference found, but not supported in "
            "multi-echo workflows at this time. Ignoring.")
    elif files:
        sbref_file = files[0]
        sbbase = os.path.basename(sbref_file)
        if len(files) > 1:
            config.loggers.workflow.warning(
                "Multiple single-band reference files found for {}; using "
                "{}".format(refbase, sbbase))
        else:
            config.loggers.workflow.info("Using single-band reference file %s.",
                                         sbbase)
    else:
        config.loggers.workflow.info("No single-band-reference found for %s.",
                                     refbase)

    metadata = layout.get_metadata(ref_file)

    # Find fieldmaps. Options: (phase1|phase2|phasediff|epi|fieldmap|syn)
    fmaps = None
    if 'fieldmaps' not in config.workflow.ignore:
        fmaps = fieldmap_wrangler(layout, ref_file,
                                  use_syn=config.workflow.use_syn_sdc,
                                  force_syn=config.workflow.force_syn)
    elif config.workflow.use_syn_sdc or config.workflow.force_syn:
        # If fieldmaps are not enabled, activate SyN-SDC in unforced (False) mode
        fmaps = {'syn': False}

    # Short circuits: (True and True and (False or 'TooShort')) == 'TooShort'
    run_stc = (bool(metadata.get("SliceTiming")) and
               'slicetiming' not in config.workflow.ignore and
               (_get_series_len(ref_file) > 4 or "TooShort"))

    # Build workflow
    workflow = Workflow(name=wf_name)
    workflow.__postdesc__ = """\
All resamplings can be performed with *a single interpolation
step* by composing all the pertinent transformations (i.e. head-motion
transform matrices, susceptibility distortion correction when available,
and co-registrations to anatomical and output spaces).
Gridded (volumetric) resamplings were performed using `antsApplyTransforms` (ANTs),
configured with Lanczos interpolation to minimize the smoothing
effects of other kernels [@lanczos].
"""

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['asl_file',
                't1w_preproc', 't1w_mask', 't1w_dseg', 't1w_tpms',
                'anat2std_xfm', 'std2anat_xfm', 'template'
                ]),
        name='inputnode')
    inputnode.inputs.asl_file = asl_file
    subj_dir=str(config.execution.bids_dir) + '/sub-' + str(config.execution.participant_label[0])
    if sbref_file is not None:
        from ...niworkflows.interfaces.images import ValidateImage
        val_sbref = pe.Node(ValidateImage(in_file=sbref_file), name='val_sbref')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['asl_t1', 'asl_t1_ref', 'asl_mask_t1',
                'asl_std', 'asl_std_ref', 'asl_mask_std', 
                'asl_native','cbf_t1', 'cbf_std', 'meancbf_t1', 
                'meancbf_std', 'score_t1', 'score_std',
                'avgscore_t1', 'avgscore_std', ' scrub_t1', 'scrub_std',
                'basil_t1', 'basil_std', 'pv_t1', 'pv_std', 
                'pv_native','att','att_t1','att_std',
                'confounds', 'confounds_metadata', 'qc_file']),
        name='outputnode')
      
    # Generate a brain-masked conversion of the t1w
    t1w_brain = pe.Node(ApplyMask(), name='t1w_brain')

    # asl buffer: an identity used as a pointer to either the original asl
    # or the STC'ed one for further use.
    aslbuffer = pe.Node(niu.IdentityInterface(fields=['asl_file']), name='aslbuffer')

    summary = pe.Node(
        FunctionalSummary(
            slice_timing=run_stc,
            registration=('FSL'),
            registration_dof=config.workflow.asl2t1w_dof,
            registration_init=config.workflow.asl2t1w_init,
            pe_direction=metadata.get("PhaseEncodingDirection"),
            tr=metadata.get("RepetitionTime")),
        name='summary', mem_gb=config.DEFAULT_MEMORY_MIN_GB, run_without_submitting=True)
    #summary.inputs.dummy_scans = config.workflow.dummy_scans

    asl_derivatives_wf = init_asl_derivatives_wf(
        bids_root=layout.root,
        metadata=metadata,
        output_dir=output_dir,
        spaces=spaces,
        scorescrub=scorescrub,
        basil=basil
    )

    workflow.connect([
        (outputnode, asl_derivatives_wf, [
            ('asl_t1', 'inputnode.asl_t1'),
            ('asl_t1_ref', 'inputnode.asl_t1_ref'),
            ('asl_mask_t1', 'inputnode.asl_mask_t1'),
            ('asl_native', 'inputnode.asl_native'),
            ('confounds', 'inputnode.confounds'),
            ('confounds_metadata', 'inputnode.confounds_metadata'),
        ]),
    ])

    # Generate a tentative aslref
    asl_reference_wf = init_asl_reference_wf(omp_nthreads=omp_nthreads)
    asl_reference_wf.inputs.inputnode.dummy_scans = 0
    if sbref_file is not None:
        workflow.connect([
            (val_sbref, asl_reference_wf, [('out_file', 'inputnode.sbref_file')]),
        ])

    # Top-level asl splitter
    asl_split = pe.Node(FSLSplit(dimension='t'), name='asl_split',
                         mem_gb=mem_gb['filesize'] * 3)

    # HMC on the asl
    asl_hmc_wf = init_asl_hmc_wf(name='asl_hmc_wf',
                                   mem_gb=mem_gb['filesize'],
                                   omp_nthreads=omp_nthreads)

    # calculate asl registration to T1w
    asl_reg_wf = init_asl_reg_wf(
        asl2t1w_dof=config.workflow.asl2t1w_dof,
        asl2t1w_init=config.workflow.asl2t1w_init,
        mem_gb=mem_gb['resampled'],
        name='asl_reg_wf',
        omp_nthreads=omp_nthreads,
        sloppy=config.execution.debug,
        use_bbr=config.workflow.use_bbr,
        use_compression=False,
    )

    # apply asl registration to T1w
    nonstd_spaces = set(spaces.get_nonstandard())
    t1cbfspace=False
    if nonstd_spaces.intersection(('T1w', 'anat')):
        t1cbfspace=True
    
    asl_t1_trans_wf = init_asl_t1_trans_wf(name='asl_t1_trans_wf',
                                             use_fieldwarp=bool(fmaps),
                                             multiecho=multiecho,
                                             cbft1space=t1cbfspace,
                                             scorescrub=scorescrub,
                                             basil=basil,
                                             mem_gb=mem_gb['resampled'],
                                             omp_nthreads=omp_nthreads,
                                             use_compression=False)

    # get confounds
    asl_confounds_wf = init_asl_confs_wf(
        mem_gb=mem_gb['largemem'],
        metadata=metadata,
        name='asl_confounds_wf')
    asl_confounds_wf.get_node('inputnode').inputs.t1_transform_flags = [False]

    # Apply transforms in 1 shot
    # Only use uncompressed output if AROMA is to be run
    asl_asl_trans_wf = init_asl_preproc_trans_wf(
        mem_gb=mem_gb['resampled'],
        omp_nthreads=omp_nthreads,
        use_compression=not config.execution.low_mem,
        use_fieldwarp=bool(fmaps),
        name='asl_asl_trans_wf'
    )
    asl_asl_trans_wf.inputs.inputnode.name_source = ref_file

    #refinemaskj = pe.Node(refinemask(),mem_gb=0.2, 
                                       #run_without_submitting=True, 
                                       #name="refinemask")

    # SLICE-TIME CORRECTION (or bypass) #############################################
    if run_stc is True:  # bool('TooShort') == True, so check True explicitly
        asl_stc_wf = init_asl_stc_wf(name='asl_stc_wf', metadata=metadata)
        workflow.connect([
            (asl_reference_wf, asl_stc_wf, [
                ('outputnode.skip_vols', 'inputnode.skip_vols')]),
            (asl_stc_wf, aslbuffer, [('outputnode.stc_file', 'asl_file')]),
        ])
        if not multiecho:
            workflow.connect([
                (asl_reference_wf, asl_stc_wf, [
                    ('outputnode.asl_file', 'inputnode.asl_file')])])
        else:  # for meepi, iterate through stc_wf for all workflows
            meepi_echos = aslbuffer.clone(name='meepi_echos')
            meepi_echos.iterables = ('asl_file', asl_file)
            workflow.connect([
                (meepi_echos, asl_stc_wf, [('asl_file', 'inputnode.asl_file')])])
    elif not multiecho:  # STC is too short or False
        # bypass STC from original asl to the splitter through aslbuffer
        workflow.connect([
            (asl_reference_wf, aslbuffer, [('outputnode.asl_file', 'asl_file')])])
    else:
        # for meepi, iterate over all meepi echos to aslbuffer
        aslbuffer.iterables = ('asl_file', asl_file)

    # SDC (SUSCEPTIBILITY DISTORTION CORRECTION) or bypass ##########################
    asl_sdc_wf = init_sdc_estimate_wf(fmaps, metadata,
                                       omp_nthreads=omp_nthreads,
                                       debug=config.execution.debug)

    # MULTI-ECHO EPI DATA #############################################
    if multiecho:
        from ...niworkflows.func.util import init_skullstrip_asl_wf
        skullstrip_asl_wf = init_skullstrip_asl_wf(name='skullstrip_asl_wf')

        inputnode.inputs.asl_file = ref_file  # Replace reference w first echo

        join_echos = pe.JoinNode(niu.IdentityInterface(fields=['asl_files']),
                                 joinsource=('meepi_echos' if run_stc is True else 'aslbuffer'),
                                 joinfield=['asl_files'],
                                 name='join_echos')

        # create optimal combination, adaptive T2* map
        asl_t2s_wf = init_asl_t2s_wf(echo_times=tes,
                                       mem_gb=mem_gb['resampled'],
                                       omp_nthreads=omp_nthreads,
                                       name='asl_t2smap_wf')

        workflow.connect([
            (skullstrip_asl_wf, join_echos, [
                ('outputnode.skull_stripped_file', 'asl_files')]),
            (join_echos, asl_t2s_wf, [
                ('asl_files', 'inputnode.asl_file')]),
        ])

    # MAIN WORKFLOW STRUCTURE #######################################################
    workflow.connect([
        (inputnode, t1w_brain, [('t1w_preproc', 'in_file'),
                                ('t1w_mask', 'in_mask')]),
        # Generate early reference
        (inputnode, asl_reference_wf, [('asl_file', 'inputnode.asl_file')]),
        # asl buffer has slice-time corrected if it was run, original otherwise
        (aslbuffer, asl_split, [('asl_file', 'in_file')]),
        # HMC
        (asl_reference_wf, asl_hmc_wf, [
            ('outputnode.raw_ref_image', 'inputnode.raw_ref_image'),
            ('outputnode.asl_file', 'inputnode.asl_file')]),
        #(asl_reference_wf, summary, [
            #('outputnode.algo_dummy_scans', 'algo_dummy_scans')]),
        # EPI-T1 registration workflow
        (inputnode, asl_reg_wf, [
            ('t1w_dseg', 'inputnode.t1w_dseg'),]),
        (t1w_brain, asl_reg_wf, [
            ('out_file', 'inputnode.t1w_brain')]),
        (inputnode, asl_t1_trans_wf, [
            ('asl_file', 'inputnode.name_source'),
            ('t1w_mask', 'inputnode.t1w_mask'),
           ]),
        (t1w_brain, asl_t1_trans_wf, [
            ('out_file', 'inputnode.t1w_brain')]),
        # unused if multiecho, but this is safe
        (asl_hmc_wf, asl_t1_trans_wf, [('outputnode.xforms', 'inputnode.hmc_xforms')]),
        (asl_reg_wf, asl_t1_trans_wf, [
            ('outputnode.itk_asl_to_t1', 'inputnode.itk_asl_to_t1')]),
        (asl_t1_trans_wf, outputnode, [('outputnode.asl_t1', 'asl_t1'),
                                        ('outputnode.asl_t1_ref', 'asl_t1_ref'),]),
        (asl_reg_wf, summary, [('outputnode.fallback', 'fallback')]),
        # SDC (or pass-through workflow)
        (t1w_brain, asl_sdc_wf, [
            ('out_file', 'inputnode.t1w_brain')]),
        (asl_reference_wf, asl_sdc_wf, [
            ('outputnode.ref_image', 'inputnode.epi_file'),
            ('outputnode.ref_image_brain', 'inputnode.epi_brain'),
            ('outputnode.asl_mask', 'inputnode.epi_mask')]),
        (asl_sdc_wf, asl_t1_trans_wf, [
            ('outputnode.out_warp', 'inputnode.fieldwarp'),
            ('outputnode.epi_mask', 'inputnode.ref_asl_mask'),
            ('outputnode.epi_brain', 'inputnode.ref_asl_brain')]),
        (asl_sdc_wf, asl_asl_trans_wf, [
            ('outputnode.out_warp', 'inputnode.fieldwarp'),
            ('outputnode.epi_mask', 'inputnode.asl_mask')]),
        (asl_sdc_wf, asl_reg_wf, [
            ('outputnode.epi_brain', 'inputnode.ref_asl_brain')]),
        (asl_sdc_wf, summary, [('outputnode.method', 'distortion_correction')]),
        # Connect asl_confounds_wf
        (inputnode, asl_confounds_wf, [('t1w_tpms', 'inputnode.t1w_tpms'),
                                        ('t1w_mask', 'inputnode.t1w_mask')]),
        (asl_hmc_wf, asl_confounds_wf, [
            ('outputnode.movpar_file', 'inputnode.movpar_file')]),
        (asl_reg_wf, asl_confounds_wf, [
            ('outputnode.itk_t1_to_asl', 'inputnode.t1_asl_xform')]),
        (asl_reference_wf, asl_confounds_wf, [
            ('outputnode.skip_vols', 'inputnode.skip_vols')]),
        (asl_asl_trans_wf, asl_confounds_wf, [
            ('outputnode.asl_mask', 'inputnode.asl_mask'),
        ]),
        (asl_confounds_wf, outputnode, [
            ('outputnode.confounds_file', 'confounds'),
        ]),
        (asl_confounds_wf, outputnode, [
            ('outputnode.confounds_metadata', 'confounds_metadata'),
        ]),
        # Connect asl_asl_trans_wf
        (asl_split, asl_asl_trans_wf, [
            ('out_files', 'inputnode.asl_file')]),
        (asl_hmc_wf, asl_asl_trans_wf, [
            ('outputnode.xforms', 'inputnode.hmc_xforms')]),
        # Summary
        (outputnode, summary, [('confounds', 'confounds_file')]),
    ])

    # for standard EPI data, pass along correct file
    if not multiecho:
        workflow.connect([
            (inputnode, asl_derivatives_wf, [
                ('asl_file', 'inputnode.source_file')]),
            (asl_asl_trans_wf, asl_confounds_wf, [
                ('outputnode.asl', 'inputnode.asl')]),
            (asl_split, asl_t1_trans_wf, [
                ('out_files', 'inputnode.asl_split')]),
        ])
    else:  # for meepi, create and use optimal combination
        workflow.connect([
            # update name source for optimal combination
            (inputnode, asl_derivatives_wf, [
                (('asl_file', combine_meepi_source), 'inputnode.source_file')]),
            (asl_asl_trans_wf, skullstrip_asl_wf, [
                ('outputnode.asl', 'inputnode.in_file')]),
            (asl_t2s_wf, asl_confounds_wf, [
                ('outputnode.asl', 'inputnode.asl')]),
            (asl_t2s_wf, asl_t1_trans_wf, [
                ('outputnode.asl', 'inputnode.asl_split')]),
        ])

    # compute  the CBF here
    compt_cbf_wf = init_cbf_compt_wf(name='compt_cbf_wf',
                                     mem_gb=mem_gb['filesize'],
                                     omp_nthreads=omp_nthreads,
                                     dummy_vols=dummyvols,
                                     M0Scale=mscale,
                                     bids_dir=subj_dir,
                                     scorescrub=scorescrub,
                                     basil=basil,
                                     smooth_kernel=smoothkernel,
                                     metadata=metadata)

    # cbf computation workflow
    workflow.connect([
         (asl_asl_trans_wf, compt_cbf_wf, [('outputnode.asl', 'inputnode.asl_file'),
                                             ('outputnode.asl_mask', 'inputnode.asl_mask')]),
         (inputnode, compt_cbf_wf, [('t1w_tpms', 'inputnode.t1w_tpms'),
                                    ('asl_file', 'inputnode.in_file')]),
         (asl_reg_wf, compt_cbf_wf, [('outputnode.itk_t1_to_asl', 'inputnode.t1_asl_xform')]),
         (asl_reg_wf, compt_cbf_wf, [('outputnode.itk_asl_to_t1', 'inputnode.itk_asl_to_t1')]),
         (inputnode, compt_cbf_wf, [('t1w_mask', 'inputnode.t1w_mask')]),
     ])


    refine_mask = pe.Node(refinemask(), mem_gb=1.0, 
                                        run_without_submitting=True, 
                                        name="refinemask")
    workflow.connect([
        (asl_asl_trans_wf, refine_mask, 
                           [('outputnode.asl_mask', 'in_aslmask')]),
        (asl_reg_wf, refine_mask, 
                        [('outputnode.itk_t1_to_asl', 'transforms')]),
        (inputnode, refine_mask, 
                        [('t1w_mask', 'in_t1mask')]),
    ])

    if fmaps:
        from sdcflows.workflows.outputs import init_sdc_unwarp_report_wf
        # Report on asl correction
        fmap_unwarp_report_wf = init_sdc_unwarp_report_wf()
        workflow.connect([
            (inputnode, fmap_unwarp_report_wf, [
                ('t1w_dseg', 'inputnode.in_seg')]),
            (asl_reference_wf, fmap_unwarp_report_wf, [
                ('outputnode.ref_image', 'inputnode.in_pre')]),
            (asl_reg_wf, fmap_unwarp_report_wf, [
                ('outputnode.itk_t1_to_asl', 'inputnode.in_xfm')]),
            (asl_sdc_wf, fmap_unwarp_report_wf, [
                ('outputnode.epi_corrected', 'inputnode.in_post')]),
        ])

        # Overwrite ``out_path_base`` of unwarping DataSinks
        # And ensure echo is dropped from report
        for node in fmap_unwarp_report_wf.list_node_names():
            if node.split('.')[-1].startswith('ds_'):
                fmap_unwarp_report_wf.get_node(node).interface.out_path_base = 'aslprep'
                fmap_unwarp_report_wf.get_node(node).inputs.dismiss_entities = ("echo",)

        for node in asl_sdc_wf.list_node_names():
            if node.split('.')[-1].startswith('ds_'):
                asl_sdc_wf.get_node(node).interface.out_path_base = 'aslprep'
                asl_sdc_wf.get_node(node).inputs.dismiss_entities = ("echo",)

        if 'syn' in fmaps:
            sdc_select_std = pe.Node(
                KeySelect(fields=['std2anat_xfm']),
                name='sdc_select_std', run_without_submitting=True)
            sdc_select_std.inputs.key = 'MNI152NLin2009cAsym'
            workflow.connect([
                (inputnode, sdc_select_std, [('std2anat_xfm', 'std2anat_xfm'),
                                             ('template', 'keys')]),
                (sdc_select_std, asl_sdc_wf, [('std2anat_xfm', 'inputnode.std2anat_xfm')]),
            ])

        if fmaps.get('syn') is True:  # SyN forced
            syn_unwarp_report_wf = init_sdc_unwarp_report_wf(
                name='syn_unwarp_report_wf', forcedsyn=True)
            workflow.connect([
                (inputnode, syn_unwarp_report_wf, [
                    ('t1w_dseg', 'inputnode.in_seg')]),
                (asl_reference_wf, syn_unwarp_report_wf, [
                    ('outputnode.ref_image', 'inputnode.in_pre')]),
                (asl_reg_wf, syn_unwarp_report_wf, [
                    ('outputnode.itk_t1_to_asl', 'inputnode.in_xfm')]),
                (asl_sdc_wf, syn_unwarp_report_wf, [
                    ('outputnode.syn_ref', 'inputnode.in_post')]),
            ])

            # Overwrite ``out_path_base`` of unwarping DataSinks
            # And ensure echo is dropped from report
            for node in syn_unwarp_report_wf.list_node_names():
                if node.split('.')[-1].startswith('ds_'):
                    syn_unwarp_report_wf.get_node(node).interface.out_path_base = 'aslprep'
                    syn_unwarp_report_wf.get_node(node).inputs.dismiss_entities = ("echo",)

    # Map final asl mask into T1w space (if required)
    nonstd_spaces = set(spaces.get_nonstandard())
    if nonstd_spaces.intersection(('T1w', 'anat')):
        from ...niworkflows.interfaces.fixes import (
            FixHeaderApplyTransforms as ApplyTransforms
        )

        aslmask_to_t1w = pe.Node(ApplyTransforms(interpolation='MultiLabel'),
                                  name='aslmask_to_t1w', mem_gb=0.1)
        workflow.connect([
            (asl_reg_wf, aslmask_to_t1w, [
                ('outputnode.itk_asl_to_t1', 'transforms')]),
            (asl_t1_trans_wf, aslmask_to_t1w, [
                ('outputnode.asl_mask_t1', 'reference_image')]),
            (refine_mask, aslmask_to_t1w, [
                ('out_mask', 'input_image')]),
            (aslmask_to_t1w, outputnode, [
                ('output_image', 'asl_mask_t1')]),
        ])
        workflow.connect([
            (compt_cbf_wf, asl_t1_trans_wf, [('outputnode.out_cbf', 'inputnode.cbf'),
                                              ('outputnode.out_mean', 'inputnode.meancbf'),
                                              ]),
            (asl_t1_trans_wf, asl_derivatives_wf, [('outputnode.cbf_t1', 'inputnode.cbf_t1'),
                                            ('outputnode.meancbf_t1', 'inputnode.meancbf_t1'),
                                
                                            
                                            ]),
                          ])

        if scorescrub: 
            workflow.connect([
                (compt_cbf_wf, asl_t1_trans_wf, [('outputnode.out_score', 'inputnode.score'),
                                              ('outputnode.out_avgscore', 'inputnode.avgscore'),
                                              ('outputnode.out_scrub', 'inputnode.scrub'),]),
               (asl_t1_trans_wf, asl_derivatives_wf,[
                                            ('outputnode.scrub_t1', 'inputnode.scrub_t1'),
                                            ('outputnode.score_t1', 'inputnode.score_t1'),
                                            ('outputnode.avgscore_t1', 'inputnode.avgscore_t1'),])
            ])
        if basil:
            workflow.connect([
                (compt_cbf_wf, asl_t1_trans_wf, [('outputnode.out_cbfb', 'inputnode.basil'),
                                            ('outputnode.out_cbfpv', 'inputnode.pv'),
                                            ('outputnode.out_att', 'inputnode.att'),]),
               (asl_t1_trans_wf, asl_derivatives_wf,[
                                            ('outputnode.basil_t1', 'inputnode.basil_t1'),
                                            ('outputnode.pv_t1', 'inputnode.pv_t1'),
                                            ('outputnode.att_t1', 'inputnode.att_t1'),])
            ])


    if nonstd_spaces.intersection(('func', 'run', 'asl', 'aslref', 'sbref')):
        workflow.connect([
            (asl_asl_trans_wf, outputnode, [
                ('outputnode.asl', 'asl_native')]),
            (asl_asl_trans_wf, asl_derivatives_wf, [
                ('outputnode.asl_ref', 'inputnode.asl_native_ref')]),
            (refine_mask, asl_derivatives_wf, [
                ('out_mask', 'inputnode.asl_mask_native' )]),
            (compt_cbf_wf, asl_derivatives_wf, [
                ('outputnode.out_cbf', 'inputnode.cbf'),
                ('outputnode.out_mean', 'inputnode.meancbf'),
                ]),
        ])

        if scorescrub:
            workflow.connect([
            (compt_cbf_wf, asl_derivatives_wf, [
                ('outputnode.out_score', 'inputnode.score'),
                ('outputnode.out_avgscore', 'inputnode.avgscore'),
                ('outputnode.out_scrub', 'inputnode.scrub'),])
            ])
        if basil:
            workflow.connect([
            (compt_cbf_wf, asl_derivatives_wf, [
                ('outputnode.out_cbfb', 'inputnode.basil'),
                ('outputnode.out_cbfpv', 'inputnode.pv'),
                ('outputnode.out_att', 'inputnode.att')]),
            ])

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # Apply transforms in 1 shot
        # Only use uncompressed output if AROMA is to be run
        asl_std_trans_wf = init_asl_std_trans_wf(
            mem_gb=mem_gb['resampled'],
            omp_nthreads=omp_nthreads,
            spaces=spaces,
            scorescrub=scorescrub,
            basil=basil,
            name='asl_std_trans_wf',
            use_compression=not config.execution.low_mem,
            use_fieldwarp=bool(fmaps),
        )
        workflow.connect([
            (inputnode, asl_std_trans_wf, [
                ('template', 'inputnode.templates'),
                ('anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('asl_file', 'inputnode.name_source'),]),
            (asl_hmc_wf, asl_std_trans_wf, [
                ('outputnode.xforms', 'inputnode.hmc_xforms')]),
            (asl_reg_wf, asl_std_trans_wf, [
                ('outputnode.itk_asl_to_t1', 'inputnode.itk_asl_to_t1')]),
            (refine_mask, asl_std_trans_wf, [
                ('out_mask', 'inputnode.asl_mask')]),
            (asl_sdc_wf, asl_std_trans_wf, [
                ('outputnode.out_warp', 'inputnode.fieldwarp')]),
            (compt_cbf_wf, asl_std_trans_wf, [('outputnode.out_cbf', 'inputnode.cbf'),
                                              ('outputnode.out_mean', 'inputnode.meancbf'),
                                               ]),
            ])

        if scorescrub:
            workflow.connect([
                  (compt_cbf_wf, asl_std_trans_wf, [
                                               ('outputnode.out_score', 'inputnode.score'),
                                               ('outputnode.out_avgscore', 'inputnode.avgscore'),
                                               ('outputnode.out_scrub', 'inputnode.scrub'),]),  
                ])
        if basil:
            workflow.connect([
                  (compt_cbf_wf, asl_std_trans_wf, [
                                               ('outputnode.out_cbfb', 'inputnode.basil'),
                                               ('outputnode.out_cbfpv', 'inputnode.pv'),
                                               ('outputnode.out_att', 'inputnode.att'),]),  
                ])

        if not multiecho:
            workflow.connect([
                (asl_split, asl_std_trans_wf, [
                    ('out_files', 'inputnode.asl_split')])
            ])
        else:
            split_opt_comb = asl_split.clone(name='split_opt_comb')
            workflow.connect([
                (asl_t2s_wf, split_opt_comb, [
                    ('outputnode.asl', 'in_file')]),
                (split_opt_comb, asl_std_trans_wf, [
                    ('out_files', 'inputnode.asl_split')
                ])
            ])

        # asl_derivatives_wf internally parametrizes over snapshotted spaces.
        # asl_derivatives_wf internally parametrizes over snapshotted spaces.
        workflow.connect([
            (asl_std_trans_wf, asl_derivatives_wf, [
                ('outputnode.template', 'inputnode.template'),
                ('outputnode.spatial_reference', 'inputnode.spatial_reference'),
                ('outputnode.asl_std_ref', 'inputnode.asl_std_ref'),
                ('outputnode.asl_std', 'inputnode.asl_std'),
                ('outputnode.asl_mask_std', 'inputnode.asl_mask_std'),
                ('outputnode.cbf_std', 'inputnode.cbf_std'),
                ('outputnode.meancbf_std','inputnode.meancbf_std')]),
        ])
        if scorescrub:
            workflow.connect([
            (asl_std_trans_wf, asl_derivatives_wf, [
                ('outputnode.score_std', 'inputnode.score_std'),
                ('outputnode.avgscore_std', 'inputnode.avgscore_std'),
                ('outputnode.scrub_std', 'inputnode.scrub_std'),]),
            ])

        if basil:
            workflow.connect([
            (asl_std_trans_wf, asl_derivatives_wf, [
                ('outputnode.basil_std', 'inputnode.basil_std'),
                ('outputnode.pv_std', 'inputnode.pv_std'),
                ('outputnode.att_std', 'inputnode.att_std')]),
            ])

    

    compt_qccbf_wf = init_cbfqc_compt_wf(name='compt_qccbf_wf',
                                         mem_gb=mem_gb['filesize'],
                                         omp_nthreads=omp_nthreads,
                                         asl_file=asl_file,
                                         scorescrub=scorescrub,
                                         basil=basil,
                                         metadata=metadata)
    workflow.connect([
         (refine_mask, compt_qccbf_wf, [('out_mask', 'inputnode.asl_mask')]),
         (inputnode, compt_qccbf_wf, [('t1w_tpms', 'inputnode.t1w_tpms')]),
         (asl_reg_wf, compt_qccbf_wf, [('outputnode.itk_t1_to_asl', 'inputnode.t1_asl_xform')]),
         (inputnode, compt_qccbf_wf, [('t1w_mask', 'inputnode.t1w_mask')]),
         (compt_cbf_wf, compt_qccbf_wf, [('outputnode.out_mean', 'inputnode.meancbf')]),
         (asl_confounds_wf, compt_qccbf_wf, [
            ('outputnode.confounds_file', 'inputnode.confmat')]),
         (compt_qccbf_wf, outputnode, [('outputnode.qc_file', 'qc_file')]),
         (compt_qccbf_wf, asl_derivatives_wf, [('outputnode.qc_file', 'inputnode.qc_file')]),
         (compt_qccbf_wf, summary, [('outputnode.qc_file', 'qc_file')]),
    ])

    if scorescrub:
        workflow.connect([
        (compt_cbf_wf, compt_qccbf_wf,[('outputnode.out_avgscore', 'inputnode.avgscore'),
                                         ('outputnode.out_scrub', 'inputnode.scrub'),]),
        ])
    if basil:
        workflow.connect([
        (compt_cbf_wf, compt_qccbf_wf,[('outputnode.out_cbfb', 'inputnode.basil'),
                                         ('outputnode.out_cbfpv', 'inputnode.pv'),]),
        ])

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        workflow.connect([
                (asl_std_trans_wf, compt_qccbf_wf, [('outputnode.asl_mask_std',
                                                      'inputnode.asl_mask_std')]),
            ])

    cbf_plot = init_cbfplot_wf(mem_gb=mem_gb['filesize'], metadata=metadata,scorescrub=scorescrub,
                               basil=basil,omp_nthreads=omp_nthreads, name='cbf_plot')
    workflow.connect([
       (compt_cbf_wf, cbf_plot, [('outputnode.out_mean', 'inputnode.cbf'),
                                 ('outputnode.out_avgscore', 'inputnode.score'),
                                 ('outputnode.out_scrub', 'inputnode.scrub'),
                                 ('outputnode.out_cbfb', 'inputnode.basil'),
                                 ('outputnode.out_cbfpv', 'inputnode.pvc'),
                                 ('outputnode.out_score', 'inputnode.score_ts'),
                                 ('outputnode.out_cbf', 'inputnode.cbf_ts')]),
       (inputnode, cbf_plot, [('std2anat_xfm', 'inputnode.std2anat_xfm')]),
       (asl_reg_wf, cbf_plot, [('outputnode.itk_t1_to_asl', 'inputnode.t1_asl_xform')]),
       (refine_mask, cbf_plot, [('out_mask', 'inputnode.asl_mask')]),
       (asl_reference_wf, cbf_plot, [('outputnode.ref_image_brain', 'inputnode.asl_ref')]),
       (asl_confounds_wf, cbf_plot, [('outputnode.confounds_file', 'inputnode.confounds_file')]),
       (compt_cbf_wf, cbf_plot, [('outputnode.out_scoreindex', 'inputnode.scoreindex')]),
       ])

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        carpetplot_wf = init_carpetplot_wf(
            mem_gb=mem_gb['resampled'],
            metadata=metadata,
            # cifti_output=config.workflow.cifti_output,
            name='carpetplot_wf')

            # Xform to 'MNI152NLin2009cAsym' is always computed.
        carpetplot_select_std = pe.Node(
                KeySelect(fields=['std2anat_xfm'], key='MNI152NLin2009cAsym'),
                name='carpetplot_select_std', run_without_submitting=True)

        workflow.connect([
                (inputnode, carpetplot_select_std, [
                    ('std2anat_xfm', 'std2anat_xfm'),
                    ('template', 'keys')]),
                (carpetplot_select_std, carpetplot_wf, [
                    ('std2anat_xfm', 'inputnode.std2anat_xfm')]),
                (asl_asl_trans_wf if not multiecho else asl_t2s_wf, carpetplot_wf, [
                    ('outputnode.asl', 'inputnode.asl')]),
                (refine_mask, carpetplot_wf, [
                    ('out_mask', 'inputnode.asl_mask')]),
                (asl_reg_wf, carpetplot_wf, [
                    ('outputnode.itk_t1_to_asl', 'inputnode.t1_asl_xform')]),
            ])

        workflow.connect([
            (asl_confounds_wf, carpetplot_wf, [
                        ('outputnode.confounds_file', 'inputnode.confounds_file')])
        ])
    cbfroiqu = init_cbfroiquant_wf(mem_gb=mem_gb['filesize'],basil=basil,scorescrub=scorescrub,
                                   omp_nthreads=omp_nthreads,
                                   name='cbf_roiquant')
    workflow.connect([
            (asl_asl_trans_wf, cbfroiqu, [('outputnode.asl_mask', 'inputnode.aslmask')]),
            (inputnode, cbfroiqu, [('std2anat_xfm', 'inputnode.std2anat_xfm')]),
            (asl_reg_wf, cbfroiqu, [('outputnode.itk_t1_to_asl', 'inputnode.t1_asl_xform')]),
            (compt_cbf_wf, cbfroiqu, [('outputnode.out_mean', 'inputnode.cbf'),
                                      ]),
            (cbfroiqu, asl_derivatives_wf, [
                ('outputnode.cbf_hvoxf', 'inputnode.cbf_hvoxf'),
                ('outputnode.cbf_sc207', 'inputnode.cbf_sc207'),
                ('outputnode.cbf_sc217', 'inputnode.cbf_sc217'),
                ('outputnode.cbf_sc407', 'inputnode.cbf_sc407'),
                ('outputnode.cbf_sc417', 'inputnode.cbf_sc417'),
                
                ]),
    ])
     
    if scorescrub:
        workflow.connect([
            (compt_cbf_wf, cbfroiqu, [('outputnode.out_avgscore', 'inputnode.score'),
                                      ('outputnode.out_scrub', 'inputnode.scrub')]),
            (cbfroiqu, asl_derivatives_wf, [
                ('outputnode.score_hvoxf', 'inputnode.score_hvoxf'),
                ('outputnode.score_sc207', 'inputnode.score_sc207'),
                ('outputnode.score_sc217', 'inputnode.score_sc217'),
                ('outputnode.score_sc407', 'inputnode.score_sc407'),
                ('outputnode.score_sc417', 'inputnode.score_sc417'),
                ('outputnode.scrub_hvoxf', 'inputnode.scrub_hvoxf'),
                ('outputnode.scrub_sc207', 'inputnode.scrub_sc207'),
                ('outputnode.scrub_sc217', 'inputnode.scrub_sc217'),
                ('outputnode.scrub_sc407', 'inputnode.scrub_sc407'),
                ('outputnode.scrub_sc417', 'inputnode.scrub_sc417'),]),
        ])
    if basil:
        workflow.connect([
            (compt_cbf_wf, cbfroiqu, [('outputnode.out_cbfb', 'inputnode.basil'),
                                      ('outputnode.out_cbfpv', 'inputnode.pvc')]),
            (cbfroiqu, asl_derivatives_wf, [
                ('outputnode.basil_hvoxf', 'inputnode.basil_hvoxf'),
                ('outputnode.basil_sc207', 'inputnode.basil_sc207'),
                ('outputnode.basil_sc217', 'inputnode.basil_sc217'),
                ('outputnode.basil_sc407', 'inputnode.basil_sc407'),
                ('outputnode.basil_sc417', 'inputnode.basil_sc417'),
                ('outputnode.pvc_hvoxf', 'inputnode.pvc_hvoxf'),
                ('outputnode.pvc_sc207', 'inputnode.pvc_sc207'),
                ('outputnode.pvc_sc217', 'inputnode.pvc_sc217'),
                ('outputnode.pvc_sc407', 'inputnode.pvc_sc407'),
                ('outputnode.pvc_sc417', 'inputnode.pvc_sc417'),]),
        ])

        





    # REPORTING ############################################################
    ds_report_summary = pe.Node(
        DerivativesDataSink(desc='summary', datatype="figures", dismiss_entities=("echo",)),
        name='ds_report_summary', run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB)

    ds_report_validation = pe.Node(
        DerivativesDataSink(base_directory=output_dir, desc='validation', datatype="figures",
                            dismiss_entities=("echo",)),
        name='ds_report_validation', run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB)

    workflow.connect([
        (summary, ds_report_summary, [('out_report', 'in_file')]),
        (asl_reference_wf, ds_report_validation, [
            ('outputnode.validation_report', 'in_file')]),
    ])

    # Fill-in datasinks of reportlets seen so far
    for node in workflow.list_node_names():
        if node.split('.')[-1].startswith('ds_report'):
            workflow.get_node(node).inputs.base_directory = output_dir
            workflow.get_node(node).inputs.source_file = ref_file

    return workflow


def _get_series_len(asl_fname):
    from ...niworkflows.interfaces.registration import _get_vols_to_discard
    img = nb.load(asl_fname)
    if len(img.shape) < 4:
        return 1

    skip_vols = _get_vols_to_discard(img)

    return img.shape[3] - skip_vols


def _create_mem_gb(asl_fname):
    asl_size_gb = os.path.getsize(asl_fname) / (1024**3)
    asl_tlen = nb.load(asl_fname).shape[-1]
    mem_gb = {
        'filesize': asl_size_gb,
        'resampled': asl_size_gb * 4,
        'largemem': asl_size_gb * (max(asl_tlen / 100, 1.0) + 4),
    }

    return asl_tlen, mem_gb


def _get_wf_name(asl_fname):
    """
    Derive the workflow name for supplied asl file.

    >>> _get_wf_name('/completely/made/up/path/sub-01_task-nback_asl.nii.gz')
    'func_preproc_task_nback_wf'
    >>> _get_wf_name('/completely/made/up/path/sub-01_task-nback_run-01_echo-1_asl.nii.gz')
    'func_preproc_task_nback_run_01_echo_1_wf'

    """
    from nipype.utils.filemanip import split_filename
    fname = split_filename(asl_fname)[1]
    fname_nosub = '_'.join(fname.split("_")[1:])
    # if 'echo' in fname_nosub:
    #     fname_nosub = '_'.join(fname_nosub.split("_echo-")[:1]) + "_asl"
    name = "asl_preproc_" + fname_nosub.replace(
        ".", "_").replace(" ", "").replace("-", "_").replace("_asl", "_wf")

    return name


def _to_join(in_file, join_file):
    """Join two tsv files if the join_file is not ``None``."""
    from ...niworkflows.interfaces.utils import JoinTSVColumns
    if join_file is None:
        return in_file
    res = JoinTSVColumns(in_file=in_file, join_file=join_file).run()
    return res.outputs.out_file
