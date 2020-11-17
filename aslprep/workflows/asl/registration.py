# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Registration workflows
++++++++++++++++++++++

.. autofunction:: init_asl_reg_wf
.. autofunction:: init_asl_t1_trans_wf
.. autofunction:: init_bbreg_wf
.. autofunction:: init_fsl_bbr_wf

"""
from ... import config

import os
import os.path as op

import pkg_resources as pkgr

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl, c3

from ...interfaces import DerivativesDataSink

DEFAULT_MEMORY_MIN_GB = config.DEFAULT_MEMORY_MIN_GB
LOGGER = config.loggers.workflow


def init_asl_reg_wf(
        use_bbr,
        asl2t1w_dof,
        asl2t1w_init,
        mem_gb,
        omp_nthreads,
        name='asl_reg_wf',
        sloppy=False,
        use_compression=True,
        write_report=True,
):
    """
    Build a workflow to run same-subject, ASL-to-T1w image-registration.

    Calculates the registration between a reference ASL image and T1w-space
    using a boundary-based registration (BBR) cost function.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.registration import init_asl_reg_wf
            wf = init_asl_reg_wf( mem_gb=3,
                                  omp_nthreads=1,
                                  use_bbr=True,
                                  asl2t1w_dof=9,
                                  asl2t1w_init='register')

    Parameters
    ----------
    use_bbr : :obj:`bool` or None
        Enable/disable boundary-based registration refinement.
        If ``None``, test BBR result for distortion before accepting.
    asl2t1w_dof : 6, 9 or 12
        Degrees-of-freedom for ASL-T1w registration
    asl2t1w_init : str, 'header' or 'register'
        If ``'header'``, use header information for initialization of ASL and T1 images.
        If ``'register'``, align volumes by their centers.
    mem_gb : :obj:`float`
        Size of ASL file in GB
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    name : :obj:`str`
        Name of workflow (default: ``asl_reg_wf``)
    use_compression : :obj:`bool`
        Save registered ASL series as ``.nii.gz``
    use_fieldwarp : :obj:`bool`
        Include SDC warp in single-shot transform from ASL to T1
    write_report : :obj:`bool`
        Whether a reportlet should be stored

    Inputs
    ------
    ref_asl_brain
        Reference image to which ASL series is aligned
        If ``fieldwarp == True``, ``ref_asl_brain`` should be unwarped
    t1w_brain
        Skull-stripped ``t1w_preproc``
    t1w_dseg
        Segmentation of preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF)

    Outputs
    -------
    itk_asl_to_t1
        Affine transform from ``ref_asl_brain`` to T1 space (ITK format)
    itk_t1_to_asl
        Affine transform from T1 space to ASL space (ITK format)
    fallback
        Boolean indicating whether BBR was rejected (mri_coreg registration returned)

    See Also
    --------
      * :py:func:`~aslprep.workflows.asl.registration.init_bbreg_wf`
      * :py:func:`~aslprep.workflows.asl.registration.init_fsl_bbr_wf`

    """
    from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow

    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['ref_asl_brain', 't1w_brain', 't1w_dseg',
                    ]),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=[
            'itk_asl_to_t1', 'itk_t1_to_asl', 'fallback']),
        name='outputnode'
    )

   
   
    bbr_wf = init_fsl_bbr_wf(use_bbr=use_bbr, asl2t1w_dof=asl2t1w_dof,
                                 asl2t1w_init=asl2t1w_init, sloppy=sloppy)

    workflow.connect([
        (inputnode, bbr_wf, [
            ('ref_asl_brain', 'inputnode.in_file'),
            ('t1w_dseg', 'inputnode.t1w_dseg'),
            ('t1w_brain', 'inputnode.t1w_brain')]),
        (bbr_wf, outputnode, [('outputnode.itk_asl_to_t1', 'itk_asl_to_t1'),
                              ('outputnode.itk_t1_to_asl', 'itk_t1_to_asl'),
                              ('outputnode.fallback', 'fallback')]),
    ])

    if write_report:
        ds_report_reg = pe.Node(
            DerivativesDataSink(datatype="figures", dismiss_entities=("echo",)),
            name='ds_report_reg', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)

        def _asl_reg_suffix(fallback):
            return 'flirtbbr'

        workflow.connect([
            (bbr_wf, ds_report_reg, [
                ('outputnode.out_report', 'in_file'),
                (('outputnode.fallback', _asl_reg_suffix), 'desc')]),
        ])

    return workflow


def init_asl_t1_trans_wf( mem_gb, omp_nthreads,scorescrub=False,basil=False, cbft1space=False, 
                        multiecho=False, use_fieldwarp=False,
                          use_compression=True, name='asl_t1_trans_wf'):
    """
    Co-register the reference ASL image to T1w-space.

    The workflow uses :abbr:`BBR (boundary-based registration)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.registration import init_asl_t1_trans_wf
            wf = init_asl_t1_trans_wf(
                                       mem_gb=3,
                                       omp_nthreads=1)

    Parameters
    ----------
    use_fieldwarp : :obj:`bool`
        Include SDC warp in single-shot transform from ASLto T1
    multiecho : :obj:`bool`
        If multiecho data was supplied, HMC already performed
    mem_gb : :obj:`float`
        Size of ASL file in GB
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    use_compression : :obj:`bool`
        Save registered ASL series as ``.nii.gz``
    name : :obj:`str`
        Name of workflow (default: ``asl_reg_wf``)

    Inputs
    ------
    name_source
        ASL series NIfTI file
        Used to recover original information lost during processing
    ref_asl_brain
        Reference image to which ASL series is aligned
        If ``fieldwarp == True``, ``ref_asl_brain`` should be unwarped
    ref_asl_mask
        Skull-stripping mask of reference image
    t1w_brain
        Skull-stripped bias-corrected structural template image
    t1w_mask
        Mask of the skull-stripped template image
    asl_split
        Individual 3D ASL volumes, not motion corrected
    hmc_xforms
        List of affine transforms aligning each volume to ``ref_image`` in ITK format
    itk_asl_to_t1
        Affine transform from ``ref_asl_brain`` to T1 space (ITK format)
    fieldwarp
        a :abbr:`DFM (displacements field map)` in ITK format

    Outputs
    -------
    ASL_t1
        Motion-corrected ASL series in T1 space
    asl_t1_ref
        Reference, contrast-enhanced summary of the motion-corrected ASL series in T1w space
    asl_mask_t1
        ASL mask in T1 space

    See also
    --------
      * :py:func:`~aslprep.workflows.asl.registration.init_bbreg_wf`
      * :py:func:`~aslprep.workflows.asl.registration.init_fsl_bbr_wf`

    """
    from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from ...niworkflows.func.util import init_asl_reference_wf
    from ...niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
    from ...niworkflows.interfaces.itk import MultiApplyTransforms
    from ...niworkflows.interfaces.nilearn import Merge
    from ...niworkflows.interfaces.utils import GenerateSamplingReference

    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['name_source', 'ref_asl_brain', 'ref_asl_mask','t1w_brain', 't1w_mask', 
                    'asl_split', 'fieldwarp', 'hmc_xforms', 'cbf', 'meancbf','att',
                    'score', 'avgscore', 'scrub', 'basil', 'pv', 'itk_asl_to_t1']),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=[
            'asl_t1', 'asl_t1_ref', 'asl_mask_t1','att_t1',
            'cbf_t1', 'meancbf_t1', 'score_t1', 'avgscore_t1', 
            'scrub_t1', 'basil_t1', 'pv_t1']),
        name='outputnode'
    )

    gen_ref = pe.Node(GenerateSamplingReference(), name='gen_ref',
                      mem_gb=0.3)  # 256x256x256 * 64 / 8 ~ 150MB

    mask_t1w_tfm = pe.Node(ApplyTransforms(interpolation='MultiLabel'),
                           name='mask_t1w_tfm', mem_gb=0.1)

    workflow.connect([
        (inputnode, gen_ref, [('ref_asl_brain', 'moving_image'),
                              ('t1w_brain', 'fixed_image'),
                              ('t1w_mask', 'fov_mask')]),
        (inputnode, mask_t1w_tfm, [('ref_asl_mask', 'input_image')]),
        (gen_ref, mask_t1w_tfm, [('out_file', 'reference_image')]),
        (inputnode, mask_t1w_tfm, [('itk_asl_to_t1', 'transforms')]),
        (mask_t1w_tfm, outputnode, [('output_image', 'asl_mask_t1')]),
    ])

    asl_to_t1w_transform = pe.Node(
                  MultiApplyTransforms(interpolation="LanczosWindowedSinc", float=True, copy_dtype=True),
                  name='asl_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
    # merge 3D volumes into 4D timeseries
    merge = pe.Node(Merge(compress=use_compression), name='merge', mem_gb=mem_gb)

    # Generate a reference on the target T1w space
    gen_final_ref = init_asl_reference_wf(omp_nthreads, pre_mask=True)
    
    if not multiecho:
        # Merge transforms placing the head motion correction last
        nforms = 2 + int(use_fieldwarp)
        merge_xforms = pe.Node(niu.Merge(nforms), name='merge_xforms',
                               run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)
        if use_fieldwarp:
            workflow.connect([
                (inputnode, merge_xforms, [('fieldwarp', 'in2')])
            ])

        workflow.connect([
            # merge transforms
            (inputnode, merge_xforms, [
                ('hmc_xforms', 'in%d' % nforms),
                ('itk_asl_to_t1', 'in1')]),
            (merge_xforms, asl_to_t1w_transform, [('out', 'transforms')]),
            (inputnode, asl_to_t1w_transform, [('asl_split', 'input_image')]),
            (inputnode, merge, [('name_source', 'header_source')]),
            (gen_ref, asl_to_t1w_transform, [('out_file', 'reference_image')]),
            (asl_to_t1w_transform, merge, [('out_files', 'in_files')]),
            (merge, gen_final_ref, [('out_file', 'inputnode.asl_file')]),
            (mask_t1w_tfm, gen_final_ref, [('output_image', 'inputnode.asl_mask')]),
            (merge, outputnode, [('out_file', 'asl_t1')]),
        ])

    else:
        from nipype.interfaces.fsl import Split as FSLSplit
        asl_split = pe.Node(FSLSplit(dimension='t'), name='asl_split',
                             mem_gb=DEFAULT_MEMORY_MIN_GB)

        workflow.connect([
            (inputnode, asl_split, [('asl_split', 'in_file')]),
            (asl_split, asl_to_t1w_transform, [('out_files', 'input_image')]),
            (inputnode, asl_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
            (inputnode, merge, [('name_source', 'header_source')]),
            (gen_ref, asl_to_t1w_transform, [('out_file', 'reference_image')]),
            (asl_to_t1w_transform, merge, [('out_files', 'in_files')]),
            (merge, gen_final_ref, [('out_file', 'inputnode.asl_file')]),
            (mask_t1w_tfm, gen_final_ref, [('output_image', 'inputnode.asl_mask')]),
            (merge, outputnode, [('out_file', 'asl_t1')]),
        ])

    if cbft1space:
        cbf_to_t1w_transform = pe.Node(
               ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3,
                        dimension=3),
              name='cbf_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        meancbf_to_t1w_transform = pe.Node(
                         ApplyTransforms(interpolation="LanczosWindowedSinc", float=True,input_image_type=3),
                         name='meancbf_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)

        workflow.connect([
         
         (gen_final_ref, outputnode, [('outputnode.ref_image', 'asl_t1_ref')]),
         (inputnode, cbf_to_t1w_transform, [('cbf', 'input_image')]),
         (cbf_to_t1w_transform, outputnode, [('output_image', 'cbf_t1')]),
         (inputnode, cbf_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (gen_ref, cbf_to_t1w_transform, [('out_file', 'reference_image')]),
         (inputnode, meancbf_to_t1w_transform, [('meancbf', 'input_image')]),
         (meancbf_to_t1w_transform, outputnode, [('output_image', 'meancbf_t1')]),
         (inputnode, meancbf_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (gen_ref, meancbf_to_t1w_transform, [('out_file', 'reference_image')]),
            ])
    if cbft1space and scorescrub:
        score_to_t1w_transform = pe.Node(
             ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3,
                        dimension=3),
             name='score_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        avgscore_to_t1w_transform = pe.Node(
            ApplyTransforms(interpolation="LanczosWindowedSinc", float=True,input_image_type=3),
           name='avgscore_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        scrub_to_t1w_transform = pe.Node(
              ApplyTransforms(interpolation="LanczosWindowedSinc", float=True,input_image_type=3),
              name='scrub_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)

        workflow.connect([

         (inputnode, score_to_t1w_transform, [('score', 'input_image')]),
         (score_to_t1w_transform, outputnode, [('output_image', 'score_t1')]),
         (inputnode, score_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (gen_ref, score_to_t1w_transform, [('out_file', 'reference_image')]),

         (inputnode, avgscore_to_t1w_transform, [('avgscore', 'input_image')]),
         (avgscore_to_t1w_transform, outputnode, [('output_image', 'avgscore_t1')]),
         (inputnode, avgscore_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (gen_ref, avgscore_to_t1w_transform, [('out_file', 'reference_image')]),

         (inputnode, scrub_to_t1w_transform, [('scrub', 'input_image')]),
         (scrub_to_t1w_transform, outputnode, [('output_image', 'scrub_t1')]),
         (inputnode, scrub_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (gen_ref, scrub_to_t1w_transform, [('out_file', 'reference_image')]),

         ])
    
    if cbft1space and basil:
        basil_to_t1w_transform = pe.Node(
               ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
            name='basil_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        pv_to_t1w_transform = pe.Node(
               ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
               name='pv_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        att_to_t1w_transform = pe.Node(
               ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
               name='att_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        
        
        workflow.connect([
         (inputnode, basil_to_t1w_transform, [('basil', 'input_image')]),
         (basil_to_t1w_transform, outputnode, [('output_image', 'basil_t1')]),
         (inputnode, basil_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (gen_ref, basil_to_t1w_transform, [('out_file', 'reference_image')]),

         (inputnode, pv_to_t1w_transform, [('pv', 'input_image')]),
         (pv_to_t1w_transform, outputnode, [('output_image', 'pv_t1')]),
         (inputnode, pv_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (gen_ref, pv_to_t1w_transform, [('out_file', 'reference_image')]),

         (inputnode, att_to_t1w_transform, [('att', 'input_image')]),
         (att_to_t1w_transform, outputnode, [('output_image', 'att_t1')]),
         (inputnode, att_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (gen_ref, att_to_t1w_transform, [('out_file', 'reference_image')]),
          ])

    return workflow


def init_bbreg_wf(use_bbr, asl2t1w_dof, asl2t1w_init, omp_nthreads, name='bbreg_wf'):
    """
    Build a workflow to run FreeSurfer's ``bbregister``.

    This workflow uses FreeSurfer's ``bbregister`` to register a ASL image to
    a T1-weighted structural image.

    It is a counterpart to :py:func:`~aslprep.workflows.asl.registration.init_fsl_bbr_wf`,
    which performs the same task using FSL's FLIRT with a BBR cost function.
    The ``use_bbr`` option permits a high degree of control over registration.
    If ``False``, standard, affine coregistration will be performed using
    FreeSurfer's ``mri_coreg`` tool.
    If ``True``, ``bbregister`` will be seeded with the initial transform found
    by ``mri_coreg`` (equivalent to running ``bbregister --init-coreg``).
    If ``None``, after ``bbregister`` is run, the resulting affine transform
    will be compared to the initial transform found by ``mri_coreg``.
    Excessive deviation will result in rejecting the BBR refinement and
    accepting the original, affine registration.

    Workflow Graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.registration import init_bbreg_wf
            wf = init_bbreg_wf(use_bbr=True, asl2t1w_dof=9,
                               asl2t1w_init='register', omp_nthreads=1)


    Parameters
    ----------
    use_bbr : :obj:`bool` or None
        Enable/disable boundary-based registration refinement.
        If ``None``, test BBR result for distortion before accepting.
    asl2t1w_dof : 6, 9 or 12
        Degrees-of-freedom for ASL-T1w registration
    asl2t1w_init : str, 'header' or 'register'
        If ``'header'``, use header information for initialization of ASL and T1 images.
        If ``'register'``, align volumes by their centers.
    name : :obj:`str`, optional
        Workflow name (default: bbreg_wf)

    Inputs
    ------
    in_file
        Reference ASL image to be registered
    fsnative2t1w_xfm
        FSL-style affine matrix translating from FreeSurfer T1.mgz to T1w
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID (must have folder in SUBJECTS_DIR)
    t1w_brain
        Unused (see :py:func:`~aslprep.workflows.asl.registration.init_fsl_bbr_wf`)
    t1w_dseg
        Unused (see :py:func:`~aslprep.workflows.asl.registration.init_fsl_bbr_wf`)

    Outputs
    -------
    itk_asl_to_t1
        Affine transform from ``ref_asl_brain`` to T1 space (ITK format)
    itk_t1_to_asl
        Affine transform from T1 space to ASL space (ITK format)
    out_report
        Reportlet for assessing registration quality
    fallback
        Boolean indicating whether BBR was rejected (mri_coreg registration returned)

    """
    from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
    
    from ...niworkflows.interfaces.freesurfer import (
        PatchedBBRegisterRPT as BBRegisterRPT,
        PatchedMRICoregRPT as MRICoregRPT,
        PatchedLTAConvert as LTAConvert
    )
    from ...niworkflows.interfaces.nitransforms import ConcatenateXFMs

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The ASL reference was then co-registered to the T1w reference using
`bbregister` (FreeSurfer) which implements boundary-based registration [@bbr].
Co-registration was configured with {dof} degrees of freedom{reason}.
""".format(dof={6: 'six', 9: 'nine', 12: 'twelve'}[asl2t1w_dof],
           reason='' if asl2t1w_dof == 6 else
                  'to account for distortions remaining in the ASL reference')

    inputnode = pe.Node(
        niu.IdentityInterface([
            'in_file',
            'fsnative2t1w_xfm', 'subjects_dir', 'subject_id',  # BBRegister
            't1w_dseg', 't1w_brain']),  # FLIRT BBR
        name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(['itk_asl_to_t1', 'itk_t1_to_asl', 'out_report', 'fallback']),
        name='outputnode')

    if asl2t1w_init not in ("register", "header"):
        raise ValueError(f"Unknown ASL-T1w initialization option: {asl2t1w_init}")

    # For now make BBR unconditional - in the future, we can fall back to identity,
    # but adding the flexibility without testing seems a bit dangerous
    if asl2t1w_init == "header":
        if use_bbr is False:
            raise ValueError("Cannot disable BBR and use header registration")
        if use_bbr is None:
            LOGGER.warning("Initializing BBR with header; affine fallback disabled")
            use_bbr = True

    merge_ltas = pe.Node(niu.Merge(2), name='merge_ltas', run_without_submitting=True)
    concat_xfm = pe.Node(ConcatenateXFMs(inverse=True), name='concat_xfm')

    workflow.connect([
        # Output ITK transforms
        (inputnode, merge_ltas, [('fsnative2t1w_xfm', 'in2')]),
        (merge_ltas, concat_xfm, [('out', 'in_xfms')]),
        (concat_xfm, outputnode, [('out_xfm', 'itk_asl_to_t1')]),
        (concat_xfm, outputnode, [('out_inv', 'itk_t1_to_asl')]),
    ])

    # Define both nodes, but only connect conditionally
    mri_coreg = pe.Node(
        MRICoregRPT(dof=asl2t1w_dof, sep=[4], ftol=0.0001, linmintol=0.01,
                    generate_report=not use_bbr),
        name='mri_coreg', n_procs=omp_nthreads, mem_gb=5)

    bbregister = pe.Node(
        BBRegisterRPT(dof=asl2t1w_dof, contrast_type='t2', registered_file=True,
                      out_lta_file=True, generate_report=True),
        name='bbregister', mem_gb=12)

    # Use mri_coreg
    if asl2t1w_init == "register":
        workflow.connect([
            (inputnode, mri_coreg, [('subjects_dir', 'subjects_dir'),
                                    ('subject_id', 'subject_id'),
                                    ('in_file', 'source_file')]),
        ])

        # Short-circuit workflow building, use initial registration
        if use_bbr is False:
            workflow.connect([
                (mri_coreg, outputnode, [('out_report', 'out_report')]),
                (mri_coreg, merge_ltas, [('out_lta_file', 'in1')])])
            outputnode.inputs.fallback = True

            return workflow

    # Use bbregister
    workflow.connect([
        (inputnode, bbregister, [('subjects_dir', 'subjects_dir'),
                                 ('subject_id', 'subject_id'),
                                 ('in_file', 'source_file')]),
    ])

    if asl2t1w_init == "header":
        bbregister.inputs.init = "header"
    else:
        workflow.connect([(mri_coreg, bbregister, [('out_lta_file', 'init_reg_file')])])

    # Short-circuit workflow building, use boundary-based registration
    if use_bbr is True:
        workflow.connect([
            (bbregister, outputnode, [('out_report', 'out_report')]),
            (bbregister, merge_ltas, [('out_lta_file', 'in1')])])
        outputnode.inputs.fallback = False

        return workflow

    # Only reach this point if asl2t1w_init is "register" and use_bbr is None

    #transforms = pe.Node(niu.Merge(2), run_without_submitting=True, name='transforms')
    #reports = pe.Node(niu.Merge(2), run_without_submitting=True, name='reports')

    #lta_ras2ras = pe.MapNode(LTAConvert(out_lta=True), iterfield=['in_lta'],
                             #name='lta_ras2ras', mem_gb=2)
    #compare_transforms = pe.Node(niu.Function(function=compare_xforms), name='compare_transforms')

    #select_transform = pe.Node(niu.Select(), run_without_submitting=True, name='select_transform')
    #select_report = pe.Node(niu.Select(), run_without_submitting=True, name='select_report')

    #workflow.connect([
        #(bbregister, transforms, [('out_lta_file', 'in1')]),
        #(mri_coreg, transforms, [('out_lta_file', 'in2')]),
        # Normalize LTA transforms to RAS2RAS (inputs are VOX2VOX) and compare
        #(transforms, lta_ras2ras, [('out', 'in_lta')]),
        #(lta_ras2ras, compare_transforms, [('out_lta', 'lta_list')]),
        #(compare_transforms, outputnode, [('out', 'fallback')]),
        # Select output transform
        #(transforms, select_transform, [('out', 'inlist')]),
        #(compare_transforms, select_transform, [('out', 'index')]),
        #(select_transform, merge_ltas, [('out', 'in1')]),
        # Select output report
        #(bbregister, reports, [('out_report', 'in1')]),
        #(mri_coreg, reports, [('out_report', 'in2')]),
        #(reports, select_report, [('out', 'inlist')]),
        #(compare_transforms, select_report, [('out', 'index')]),
        #(select_report, outputnode, [('out', 'out_report')]),
    #])

    #return workflow


def init_fsl_bbr_wf(use_bbr, asl2t1w_dof, asl2t1w_init, sloppy=False, name='fsl_bbr_wf'):
    """
    Build a workflow to run FSL's ``flirt``.

    This workflow uses FSL FLIRT to register a ASL image to a T1-weighted
    structural image, using a boundary-based registration (BBR) cost function.
    It is a counterpart to :py:func:`~aslprep.workflows.asl.registration.init_bbreg_wf`,
    which performs the same task using FreeSurfer's ``bbregister``.

    The ``use_bbr`` option permits a high degree of control over registration.
    If ``False``, standard, rigid coregistration will be performed by FLIRT.
    If ``True``, FLIRT-BBR will be seeded with the initial transform found by
    the rigid coregistration.
    If ``None``, after FLIRT-BBR is run, the resulting affine transform
    will be compared to the initial transform found by FLIRT.
    Excessive deviation will result in rejecting the BBR refinement and
    accepting the original, affine registration.

    Workflow Graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.registration import init_fsl_bbr_wf
            wf = init_fsl_bbr_wf(use_bbr=True, asl2t1w_dof=9, asl2t1w_init='register')


    Parameters
    ----------
    use_bbr : :obj:`bool` or None
        Enable/disable boundary-based registration refinement.
        If ``None``, test BBR result for distortion before accepting.
    asl2t1w_dof : 6, 9 or 12
        Degrees-of-freedom for ASL-T1w registration
    asl2t1w_init : str, 'header' or 'register'
        If ``'header'``, use header information for initialization of ASL and T1 images.
        If ``'register'``, align volumes by their centers.
    name : :obj:`str`, optional
        Workflow name (default: fsl_bbr_wf)

    Inputs
    ------
    in_file
        Reference ASL image to be registered
    t1w_brain
        Skull-stripped T1-weighted structural image
    t1w_dseg
        FAST segmentation of ``t1w_brain``
    

    Outputs
    -------
    itk_asl_to_t1
        Affine transform from ``ref_asl_brain`` to T1w space (ITK format)
    itk_t1_to_asl
        Affine transform from T1 space to ASL space (ITK format)
    out_report
        Reportlet for assessing registration quality
    fallback
        Boolean indicating whether BBR was rejected (rigid FLIRT registration returned)

    """
    from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from ...niworkflows.utils.images import dseg_label as _dseg_label
    from ...niworkflows.interfaces.freesurfer import PatchedLTAConvert as LTAConvert
    from ...niworkflows.interfaces.registration import FLIRTRPT
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The ASL reference was then co-registered to the T1w reference using `flirt` [FSL {fsl_ver}, @flirt] with the boundary-based registration [@bbr]
cost-function. Co-registration was configured with nine degrees of freedom to account for distortions remaining in the ASL reference.
""".format(fsl_ver=FLIRTRPT().version or '<ver>')

    inputnode = pe.Node(
        niu.IdentityInterface([
            'in_file', 't1w_dseg', 't1w_brain']),  # FLIRT BBR
        name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(['itk_asl_to_t1', 'itk_t1_to_asl', 'out_report', 'fallback']),
        name='outputnode')

    wm_mask = pe.Node(niu.Function(function=_dseg_label), name='wm_mask')
    wm_mask.inputs.label = 2  # BIDS default is WM=2
    flt_bbr_init = pe.Node(FLIRTRPT(dof=6, generate_report=not use_bbr,
                                    uses_qform=True), name='flt_bbr_init')

    if asl2t1w_init not in ("register", "header"):
        raise ValueError(f"Unknown ASL-T1w initialization option: {asl2t1w_init}")

    if asl2t1w_init == "header":
        raise NotImplementedError("Header-based registration initialization not supported for FSL")

    invt_bbr = pe.Node(fsl.ConvertXFM(invert_xfm=True), name='invt_bbr',
                       mem_gb=DEFAULT_MEMORY_MIN_GB)

    # ASL to T1 transform matrix is from fsl, using c3 tools to convert to
    # something ANTs will like.
    fsl2itk_fwd = pe.Node(c3.C3dAffineTool(fsl2ras=True, itk_transform=True),
                          name='fsl2itk_fwd', mem_gb=DEFAULT_MEMORY_MIN_GB)
    fsl2itk_inv = pe.Node(c3.C3dAffineTool(fsl2ras=True, itk_transform=True),
                          name='fsl2itk_inv', mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow.connect([
        (inputnode, flt_bbr_init, [('in_file', 'in_file'),
                                   ('t1w_brain', 'reference')]),
        (inputnode, fsl2itk_fwd, [('t1w_brain', 'reference_file'),
                                  ('in_file', 'source_file')]),
        (inputnode, fsl2itk_inv, [('in_file', 'reference_file'),
                                  ('t1w_brain', 'source_file')]),
        (invt_bbr, fsl2itk_inv, [('out_file', 'transform_file')]),
        (fsl2itk_fwd, outputnode, [('itk_transform', 'itk_asl_to_t1')]),
        (fsl2itk_inv, outputnode, [('itk_transform', 'itk_t1_to_asl')]),
    ])

    # Short-circuit workflow building, use rigid registration
    if use_bbr is False:
        workflow.connect([
            (flt_bbr_init, invt_bbr, [('out_matrix_file', 'in_file')]),
            (flt_bbr_init, fsl2itk_fwd, [('out_matrix_file', 'transform_file')]),
            (flt_bbr_init, outputnode, [('out_report', 'out_report')]),
        ])
        outputnode.inputs.fallback = True

        return workflow

    flt_bbr = pe.Node(
        FLIRTRPT(cost_func='bbr', dof=asl2t1w_dof, generate_report=True),
        name='flt_bbr')

    FSLDIR = os.getenv('FSLDIR')
    if FSLDIR:
        flt_bbr.inputs.schedule = op.join(FSLDIR, 'etc/flirtsch/bbr.sch')
    else:
        # Should mostly be hit while building docs
        LOGGER.warning("FSLDIR unset - using packaged BBR schedule")
        flt_bbr.inputs.schedule = pkgr.resource_filename('aslprep', 'data/flirtsch/bbr.sch')

    workflow.connect([
        (inputnode, wm_mask, [('t1w_dseg', 'in_seg')]),
        (inputnode, flt_bbr, [('in_file', 'in_file')]),
        (flt_bbr_init, flt_bbr, [('out_matrix_file', 'in_matrix_file')]),
    ])

    if sloppy is True:
        downsample = pe.Node(niu.Function(
            function=_conditional_downsampling, output_names=["out_file", "out_mask"]),
            name='downsample')
        workflow.connect([
            (inputnode, downsample, [("t1w_brain", "in_file")]),
            (wm_mask, downsample, [("out", "in_mask")]),
            (downsample, flt_bbr, [('out_file', 'reference'),
                                   ('out_mask', 'wm_seg')]),
        ])
    else:
        workflow.connect([
            (inputnode, flt_bbr, [('t1w_brain', 'reference')]),
            (wm_mask, flt_bbr, [('out', 'wm_seg')]),
        ])

    # Short-circuit workflow building, use boundary-based registration
    if use_bbr is True:
        workflow.connect([
            (flt_bbr, invt_bbr, [('out_matrix_file', 'in_file')]),
            (flt_bbr, fsl2itk_fwd, [('out_matrix_file', 'transform_file')]),
            (flt_bbr, outputnode, [('out_report', 'out_report')]),
        ])
        outputnode.inputs.fallback = False

        return workflow

    #transforms = pe.Node(niu.Merge(2), run_without_submitting=True, name='transforms')
    #reports = pe.Node(niu.Merge(2), run_without_submitting=True, name='reports')

    #compare_transforms = pe.Node(niu.Function(function=compare_xforms,output_names="out"), name='compare_transforms')

    #select_transform = pe.Node(niu.Select(), run_without_submitting=True, name='select_transform')
    #select_report = pe.Node(niu.Select(), run_without_submitting=True, name='select_report')

    #fsl_to_lta = pe.MapNode(LTAConvert(out_lta=True), iterfield=['in_fsl'],
                            #name='fsl_to_lta')
    #select_transform.inputs.index= 0

    workflow.connect([
        #(flt_bbr, transforms, [('out_matrix_file', 'in1')]),
        #(flt_bbr_init, transforms, [('out_matrix_file', 'in2')]),
        #Convert FSL transforms to LTA (RAS2RAS) transforms and compare
        #(inputnode, fsl_to_lta, [('in_file', 'source_file'),
                                 #('t1w_brain', 'target_file')]),
        #(transforms, fsl_to_lta, [('out', 'in_fsl')]),
        #(fsl_to_lta, compare_transforms, [('out_lta', 'lta_list')]),
        #(compare_transforms, outputnode, [('out', 'fallback')]),
        # Select output transform
        #(transforms, select_transform, [('out', 'inlist')]),
        #(compare_transforms, select_transform, [('out', 'index')]),
        (flt_bbr, invt_bbr, [('out_matrix_file', 'in_file')]),
        (flt_bbr, fsl2itk_fwd, [('out_matrix_file', 'transform_file')]),
        #(flt_bbr, reports, [('out_report', 'in1')]),
        #(flt_bbr_init, reports, [('out_report', 'in2')]),
        #(reports, select_report, [('out', 'inlist')]),
        #(compare_transforms, select_report, [('out', 'index')]),
        (flt_bbr, outputnode, [('out_report', 'out_report')]),
    ])

    return workflow


def compare_xforms(lta_list, norm_threshold=10):
    """
    Computes a normalized displacement between two affine transforms as the
    maximum overall displacement of the midpoints of the faces of a cube, when
    each transform is applied to the cube.
    This combines displacement resulting from scaling, translation and rotation.

    Although the norm is in mm, in a scaling context, it is not necessarily
    equivalent to that distance in translation.

    We choose a default threshold of 15mm as a rough heuristic.
    Normalized displacement above 20mm showed clear signs of distortion, while
    "good" BBR refinements were frequently below 10mm displaced from the rigid
    transform.
    The 10-20mm range was more ambiguous, and 15mm chosen as a compromise.
    This is open to revisiting in either direction.

    

    Parameters
    ----------

      lta_list : :obj:`list` or :obj:`tuple` of :obj:`str`
          the two given affines in LTA format
      norm_threshold : :obj:`float`
          the upper bound limit to the normalized displacement caused by the
          second transform relative to the first (default: `15`)

    """
    from ...niworkflows.interfaces.surf import load_transform
    from nipype.algorithms.rapidart import _calc_norm_affine

    bbr_affine = load_transform(lta_list[0])
    fallback_affine = load_transform(lta_list[1])

    norm, _ = _calc_norm_affine([fallback_affine, bbr_affine], use_differences=True)

    return norm[1] > norm_threshold


def _conditional_downsampling(in_file, in_mask, zoom_th=4.0):
    """Downsamples the input dataset for sloppy mode."""
    from pathlib import Path
    import numpy as np
    import nibabel as nb
    import nitransforms as nt
    from scipy.ndimage.filters import gaussian_filter

    img = nb.load(in_file)

    zooms = np.array(img.header.get_zooms()[:3])
    if not np.any(zooms < zoom_th):
        return in_file, in_mask

    out_file = Path('desc-resampled_input.nii.gz').absolute()
    out_mask = Path('desc-resampled_mask.nii.gz').absolute()

    shape = np.array(img.shape[:3])
    scaling = zoom_th / zooms
    newrot = np.diag(scaling).dot(img.affine[:3, :3])
    newshape = np.ceil(shape / scaling).astype(int)
    old_center = img.affine.dot(np.hstack((0.5 * (shape - 1), 1.0)))[:3]
    offset = old_center - newrot.dot((newshape - 1) * 0.5)
    newaffine = nb.affines.from_matvec(newrot, offset)

    newref = nb.Nifti1Image(np.zeros(newshape, dtype=np.uint8), newaffine)
    nt.Affine(reference=newref).apply(img).to_filename(out_file)

    mask = nb.load(in_mask)
    mask.set_data_dtype(float)
    mdata = gaussian_filter(mask.get_fdata(dtype=float), scaling)
    floatmask = nb.Nifti1Image(mdata, mask.affine, mask.header)
    newmask = nt.Affine(reference=newref).apply(floatmask)
    hdr = newmask.header.copy()
    hdr.set_data_dtype(np.uint8)
    newmaskdata = (newmask.get_fdata(dtype=float) > 0.5).astype(np.uint8)
    nb.Nifti1Image(newmaskdata, newmask.affine, hdr).to_filename(out_mask)

    return str(out_file), str(out_mask)
