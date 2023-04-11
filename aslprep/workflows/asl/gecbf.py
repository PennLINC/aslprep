# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""CBF-processing workflows for GE data."""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.interfaces.cbf_computation import RefineMask
from aslprep.interfaces.reports import FunctionalSummary
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.interfaces.nibabel import ApplyMask
from aslprep.niworkflows.interfaces.utility import KeySelect
from aslprep.utils.bids import collect_run_data
from aslprep.utils.misc import _create_mem_gb, _get_wf_name
from aslprep.workflows.asl.cbf import init_gecbf_compt_wf, init_parcellate_cbf_wf
from aslprep.workflows.asl.ge_utils import (
    init_asl_geref_wf,
    init_asl_gereg_wf,
    init_asl_gestd_trans_wf,
    init_asl_t1_getrans_wf,
)
from aslprep.workflows.asl.outputs import init_asl_derivatives_wf
from aslprep.workflows.asl.plotting import init_gecbfplot_wf
from aslprep.workflows.asl.qc import init_compute_cbf_qc_wf


def init_asl_gepreproc_wf(asl_file):
    """Manage the functional preprocessing stages of ASLPrep, for GE data.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep import config
            from aslprep.workflows.asl.gecbf import init_asl_gepreproc_wf
            with mock_config():
                asl_file = (
                    config.execution.bids_dir / 'sub-01' / 'perf' /
                    'sub-01_task-restEyesOpen_asl.nii.gz'
                )
                wf = init_asl_gepreproc_wf(str(asl_file))

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
    t1w_asec
        Segmentation of structural image, done with FreeSurfer.
    t1w_aparc
        Parcellation of structural image, done with FreeSurfer.
    t1w_tpms
        List of tissue probability maps in T1w space
    template
        List of templates to target
    anat2std_xfm
        List of transform files, collated with templates
    std2anat_xfm
        List of inverse transform files, collated with templates
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    t1w2fsnative_xfm
        LTA-style affine matrix translating from T1w to FreeSurfer-conformed subject space
    fsnative2t1w_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed subject space to T1w

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
    surfaces
        asl series, resampled to FreeSurfer surfaces
    asl_cifti
        asl CIFTI image
    cifti_variant
        combination of target spaces for `asl_cifti`
    cbf_ts_t1
        cbf times series in T1w space
    mean_cbf_t1
        mean cbf   in T1w space
    cbf_ts_score_t1
        scorecbf times series in T1w space
    mean_cbf_score_t1
        mean score cbf  in T1w space
    mean_cbf_scrub_t1, mean_cbf_gm_basil_t1, mean_cbf_basil_t1
        scrub, parital volume corrected and basil cbf   in T1w space
    cbf_ts_std
        cbf times series in template space
    mean_cbf_std
        mean cbf   in template space
    cbf_ts_score_std
        scorecbf times series in template space
    mean_cbf_score_std
        mean score cbf  in template space
    mean_cbf_scrub_std, mean_cbf_gm_basil_std, mean_cbf_basil_std
        scrub, parital volume corrected and basil cbf   in template space
    qc_file
        quality control meausres

    Notes
    -----
    1.  Brain-mask T1w.
    2.  Extract averaged, smoothed M0 image and reference image (which is generally the M0 image).
    3.  Estimate transforms from ASL to T1w.
    4.  Compute CBF.
    5.  Co-register the reference ASL image to T1w-space. (what about Step 3?)
    6.  Refine the ASL brain mask.
    7.  Warp the ASL brain mask to T1w-space.
    8.  CBF plotting workflow.
    9.  CBF QC workflow.
    10. Parcellate CBF results.

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
    mem_gb = {"filesize": 1, "resampled": 1, "largemem": 1}
    asl_tlen = 10

    # Have some options handy
    layout = config.execution.layout
    omp_nthreads = config.nipype.omp_nthreads
    spaces = config.workflow.spaces
    output_dir = str(config.execution.output_dir)
    mscale = config.workflow.m0_scale
    scorescrub = config.workflow.scorescrub
    basil = config.workflow.basil
    smoothkernel = config.workflow.smooth_kernel

    if scorescrub:
        config.loggers.workflow.warning(f"SCORE/SCRUB processing will be disabled for {asl_file}")
        scorescrub = False

    ref_file = asl_file
    asl_tlen, mem_gb = _create_mem_gb(ref_file)

    wf_name = _get_wf_name(ref_file)
    config.loggers.workflow.debug(
        'Creating asl processing workflow for "%s" (%.2f GB / %d TRs). '
        "Memory resampled/largemem=%.2f/%.2f GB.",
        ref_file,
        mem_gb["filesize"],
        asl_tlen,
        mem_gb["resampled"],
        mem_gb["largemem"],
    )

    # Collect associated files
    run_data = collect_run_data(layout, ref_file, multiecho=False)
    metadata = run_data["asl_metadata"].copy()

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

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "m0scan",
                "m0scan_metadata",
                "t1w_preproc",
                "t1w_mask",
                "t1w_dseg",
                "t1w_tpms",
                "anat2std_xfm",
                "std2anat_xfm",
                "template",
            ]
        ),
        name="inputnode",
    )
    inputnode.inputs.asl_file = asl_file
    inputnode.inputs.m0scan = run_data["m0scan"]
    inputnode.inputs.m0scan_metadata = run_data["m0scan_metadata"]

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_t1",
                "aslref_t1",
                "asl_mask_t1",
                "asl_std",
                "aslref_std",
                "asl_mask_std",
                "asl_native",
                "cbf_ts_t1",
                "cbf_ts_std",
                "mean_cbf_t1",
                "mean_cbf_std",
                "cbf_ts_score_t1",
                "cbf_ts_score_std",
                "mean_cbf_score_t1",
                "mean_cbf_score_std",
                "mean_cbf_scrub_t1",
                "mean_cbf_scrub_std",
                "aslref_to_t1w_xfm",
                "itk_t1_to_asl",
                "mean_cbf_basil_t1",
                "mean_cbf_basil_std",
                "mean_cbf_gm_basil_t1",
                "mean_cbf_wm_basil_t1",
                "mean_cbf_gm_basil_std",
                "mean_cbf_wm_basil_std",
                "pvwm_native",
                "mean_cbf_gm_basil",
                "att",
                "att_t1",
                "att_std",
                "qc_file",
            ],
        ),
        name="outputnode",
    )

    # Generate a brain-masked conversion of the t1w
    t1w_brain = pe.Node(ApplyMask(), name="t1w_brain")

    # fmt:off
    workflow.connect([
        (inputnode, t1w_brain, [
            ("t1w_preproc", "in_file"),
            ("t1w_mask", "in_mask"),
        ]),
    ])
    # fmt:on

    summary = pe.Node(
        FunctionalSummary(
            registration=("FSL"),
            registration_dof=config.workflow.asl2t1w_dof,
            registration_init=config.workflow.asl2t1w_init,
            distortion_correction="No distortion correction",
            pe_direction=metadata.get("PhaseEncodingDirection"),
            tr=metadata.get("RepetitionTime", metadata["RepetitionTimePreparation"]),
        ),
        name="summary",
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )
    # summary.inputs.dummy_scans = 0

    asl_derivatives_wf = init_asl_derivatives_wf(
        bids_root=layout.root,
        metadata=metadata,
        output_dir=output_dir,
        scorescrub=scorescrub,
        basil=basil,
        spaces=spaces,
        output_confounds=False,  # GE workflow doesn't generate volume-wise confounds
    )

    # fmt:off
    workflow.connect([
        (inputnode, asl_derivatives_wf, [("asl_file", "inputnode.source_file")]),
        (outputnode, asl_derivatives_wf, [
            ("asl_native", "inputnode.asl_native"),
            ("asl_t1", "inputnode.asl_t1"),
            ("aslref_t1", "inputnode.aslref_t1"),
            ("asl_mask_t1", "inputnode.asl_mask_t1"),
        ]),
    ])
    # fmt:on

    # begin workflow
    # Extract averaged, smoothed M0 image and reference image (which is generally the M0 image).
    gen_ref_wf = init_asl_geref_wf(
        metadata=metadata,
        aslcontext=run_data["aslcontext"],
        smooth_kernel=smoothkernel,
        name="asl_gereference_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, gen_ref_wf, [
            ("asl_file", "inputnode.asl_file"),
            ("m0scan", "inputnode.m0scan"),
            ("m0scan_metadata", "inputnode.m0scan_metadata"),
        ]),
    ])
    # fmt:on

    # ASL-to-T1w registration
    asl_reg_wf = init_asl_gereg_wf(
        use_bbr=config.workflow.use_bbr,
        asl2t1w_dof=config.workflow.asl2t1w_dof,
        asl2t1w_init=config.workflow.asl2t1w_init,
        name="asl_reg_wf",
        sloppy=False,
        write_report=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, asl_reg_wf, [("t1w_dseg", "inputnode.t1w_dseg")]),
        (gen_ref_wf, asl_reg_wf, [("outputnode.ref_image_brain", "inputnode.ref_asl_brain")]),
        (t1w_brain, asl_reg_wf, [("out_file", "inputnode.t1w_brain")]),
        (asl_reg_wf, outputnode, [
            ("outputnode.aslref_to_t1w_xfm", "aslref_to_t1w_xfm"),
            ("outputnode.itk_t1_to_asl", "itk_t1_to_asl"),
        ]),
        (asl_reg_wf, asl_derivatives_wf, [
            ("outputnode.itk_t1_to_asl", "inputnode.itk_t1_to_asl"),
            ("outputnode.aslref_to_t1w_xfm", "inputnode.aslref_to_t1w_xfm"),
        ]),
        (asl_reg_wf, summary, [("outputnode.fallback", "fallback")]),
    ])
    # fmt:on

    cbf_compt_wf = init_gecbf_compt_wf(
        name_source=asl_file,
        aslcontext=run_data["aslcontext"],
        metadata=metadata,
        scorescrub=scorescrub,
        basil=basil,
        M0Scale=mscale,
        mem_gb=mem_gb["filesize"],
        name="cbf_compt_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, cbf_compt_wf, [
            ("asl_file", "inputnode.asl_file"),
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("t1w_mask", "inputnode.t1w_mask"),
        ]),
        (gen_ref_wf, cbf_compt_wf, [
            ("outputnode.asl_mask", "inputnode.asl_mask"),
            ("outputnode.m0_file", "inputnode.m0_file"),
            ("outputnode.m0tr", "inputnode.m0tr"),
        ]),
        (asl_reg_wf, cbf_compt_wf, [
            ("outputnode.aslref_to_t1w_xfm", "inputnode.aslref_to_t1w_xfm"),
            ("outputnode.itk_t1_to_asl", "inputnode.t1w_to_aslref_xfm"),
        ]),
    ])
    # fmt:on

    nonstd_spaces = set(spaces.get_nonstandard())
    t1cbfspace = False
    if nonstd_spaces.intersection(("T1w", "anat")):
        t1cbfspace = True

    t1w_gereg_wf = init_asl_t1_getrans_wf(
        mem_gb=3,
        omp_nthreads=omp_nthreads,
        cbft1space=t1cbfspace,
        scorescrub=scorescrub,
        basil=basil,
        name="asl_t1_trans_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, t1w_gereg_wf, [
            ("asl_file", "inputnode.name_source"),
            ("t1w_mask", "inputnode.t1w_mask"),
        ]),
        (gen_ref_wf, t1w_gereg_wf, [
            ("outputnode.ref_image_brain", "inputnode.ref_asl_brain"),
            ("outputnode.asl_mask", "inputnode.ref_asl_mask"),
        ]),
        (t1w_brain, t1w_gereg_wf, [("out_file", "inputnode.t1w_brain")]),
        # unused if multiecho, but this is safe
        (asl_reg_wf, t1w_gereg_wf, [
            ("outputnode.aslref_to_t1w_xfm", "inputnode.aslref_to_t1w_xfm"),
        ]),
        (t1w_gereg_wf, outputnode, [
            ("outputnode.asl_t1", "asl_t1"),
            ("outputnode.aslref_t1", "aslref_t1"),
        ]),
        (cbf_compt_wf, t1w_gereg_wf, [
            ("outputnode.cbf_ts", "inputnode.cbf"),
            ("outputnode.mean_cbf", "inputnode.mean_cbf"),
        ]),
        (t1w_gereg_wf, asl_derivatives_wf, [
            ("outputnode.cbf_ts_t1", "inputnode.cbf_ts_t1"),
            ("outputnode.mean_cbf_t1", "inputnode.mean_cbf_t1"),
        ]),
    ])
    # fmt:on

    refine_mask = pe.Node(
        RefineMask(),
        mem_gb=1.0,
        run_without_submitting=True,
        name="refine_mask",
    )

    # fmt:off
    workflow.connect([
        (inputnode, refine_mask, [("t1w_mask", "t1w_mask")]),
        (gen_ref_wf, refine_mask, [("outputnode.asl_mask", "asl_mask")]),
        (asl_reg_wf, refine_mask, [("outputnode.itk_t1_to_asl", "transforms")]),
    ])
    # fmt:on

    cbf_plot = init_gecbfplot_wf(
        scorescrub=scorescrub,
        basil=basil,
        name="cbf_plot",
    )

    # fmt:off
    workflow.connect([
        (cbf_compt_wf, cbf_plot, [("outputnode.mean_cbf", "inputnode.cbf")]),
        (gen_ref_wf, cbf_plot, [("outputnode.ref_image_brain", "inputnode.aslref")]),
    ])
    # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (cbf_compt_wf, cbf_plot, [
                ("outputnode.mean_cbf_score", "inputnode.mean_cbf_score"),
                ("outputnode.mean_cbf_scrub", "inputnode.mean_cbf_scrub"),
            ]),
        ])
        # fmt:on

    if basil:
        # fmt:off
        workflow.connect([
            (cbf_compt_wf, cbf_plot, [
                ("outputnode.mean_cbf_basil", "inputnode.basil"),
                ("outputnode.mean_cbf_gm_basil", "inputnode.pvc"),
            ]),
        ])
        # fmt:on

    if nonstd_spaces.intersection(("T1w", "anat")):
        t1cbfspace = True
        from aslprep.interfaces.ants import ApplyTransforms

        aslmask_to_t1w = pe.Node(
            ApplyTransforms(interpolation="MultiLabel"),
            name="aslmask_to_t1w",
            mem_gb=0.1,
        )
        # fmt:off
        workflow.connect([
            (asl_reg_wf, aslmask_to_t1w, [("outputnode.aslref_to_t1w_xfm", "transforms")]),
            (t1w_gereg_wf, aslmask_to_t1w, [("outputnode.asl_mask_t1", "reference_image")]),
            (refine_mask, aslmask_to_t1w, [("out_mask", "input_image")]),
            (aslmask_to_t1w, outputnode, [("output_image", "asl_mask_t1")]),
        ])
        # fmt:on

        if scorescrub:
            # fmt:off
            workflow.connect([
                (cbf_compt_wf, t1w_gereg_wf, [
                    ("outputnode.cbf_ts_score", "inputnode.cbf_ts_score"),
                    ("outputnode.mean_cbf_score", "inputnode.mean_cbf_score"),
                    ("outputnode.mean_cbf_scrub", "inputnode.mean_cbf_scrub"),
                ]),
                (t1w_gereg_wf, asl_derivatives_wf, [
                    ("outputnode.cbf_ts_score_t1", "inputnode.cbf_ts_score_t1"),
                    ("outputnode.mean_cbf_score_t1", "inputnode.mean_cbf_score_t1"),
                    ("outputnode.mean_cbf_scrub_t1", "inputnode.mean_cbf_scrub"),
                ]),
            ])
            # fmt:on

        if basil:
            # fmt:off
            workflow.connect([
                (cbf_compt_wf, t1w_gereg_wf, [
                    ("outputnode.mean_cbf_basil", "inputnode.mean_cbf_basil"),
                    ("outputnode.mean_cbf_gm_basil", "inputnode.mean_cbf_gm_basil"),
                    ("outputnode.mean_cbf_wm_basil", "inputnode.mean_cbf_wm_basil"),
                    ("outputnode.att", "inputnode.att"),
                ]),
                (t1w_gereg_wf, asl_derivatives_wf, [
                    ("outputnode.mean_cbf_basil_t1", "inputnode.mean_cbf_basil_t1"),
                    ("outputnode.mean_cbf_gm_basil_t1", "inputnode.mean_cbf_gm_basil"),
                    ("outputnode.mean_cbf_wm_basil_t1", "inputnode.mean_cbf_wm_basil"),
                    ("outputnode.att_t1", "inputnode.att_t1"),
                ]),
            ])
            # fmt:on

    compute_cbf_qc_wf = init_compute_cbf_qc_wf(
        is_ge=True,
        scorescrub=scorescrub,
        basil=basil,
        name="compute_cbf_qc_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, compute_cbf_qc_wf, [
            ("asl_file", "inputnode.name_source"),
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("t1w_mask", "inputnode.t1w_mask"),
        ]),
        (refine_mask, compute_cbf_qc_wf, [("out_mask", "inputnode.asl_mask")]),
        (asl_reg_wf, compute_cbf_qc_wf, [
            ("outputnode.itk_t1_to_asl", "inputnode.t1w_to_aslref_xfm"),
        ]),
        (cbf_compt_wf, compute_cbf_qc_wf, [("outputnode.mean_cbf", "inputnode.mean_cbf")]),
        (compute_cbf_qc_wf, outputnode, [("outputnode.qc_file", "qc_file")]),
        (compute_cbf_qc_wf, asl_derivatives_wf, [("outputnode.qc_file", "inputnode.qc_file")]),
        (compute_cbf_qc_wf, summary, [("outputnode.qc_file", "qc_file")]),
    ])
    # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (cbf_compt_wf, compute_cbf_qc_wf, [
                ("outputnode.mean_cbf_score", "inputnode.mean_cbf_score"),
                ("outputnode.mean_cbf_scrub", "inputnode.mean_cbf_scrub"),
            ]),
        ])
        # fmt:on

    if basil:
        # fmt:off
        workflow.connect([
            (cbf_compt_wf, compute_cbf_qc_wf, [
                ("outputnode.mean_cbf_basil", "inputnode.mean_cbf_basil"),
                ("outputnode.mean_cbf_gm_basil", "inputnode.mean_cbf_gm_basil"),
            ]),
        ])
        # fmt:on

    if nonstd_spaces.intersection(("func", "run", "asl")):
        # fmt:off
        workflow.connect([
            (inputnode, outputnode, [("asl_file", "asl_native")]),
            (gen_ref_wf, asl_derivatives_wf, [
                ("outputnode.raw_ref_image", "inputnode.aslref_native"),
            ]),
            (refine_mask, asl_derivatives_wf, [("out_mask", "inputnode.asl_mask_native")]),
            (cbf_compt_wf, asl_derivatives_wf, [
                ("outputnode.cbf_ts", "inputnode.cbf_ts_native"),
                ("outputnode.mean_cbf", "inputnode.mean_cbf_native"),
            ]),
        ])
        # fmt:on

        if scorescrub:
            # fmt:off
            workflow.connect([
                (cbf_compt_wf, asl_derivatives_wf, [
                    ("outputnode.cbf_ts_score", "inputnode.cbf_ts_score_native"),
                    ("outputnode.mean_cbf_score", "inputnode.mean_cbf_score_native"),
                    ("outputnode.mean_cbf_scrub", "inputnode.mean_cbf_scrub_native"),
                ]),
            ])
            # fmt:on

        if basil:
            # fmt:off
            workflow.connect([
                (cbf_compt_wf, asl_derivatives_wf, [
                    ("outputnode.mean_cbf_basil", "inputnode.mean_cbf_basil_native"),
                    ("outputnode.mean_cbf_gm_basil", "inputnode.mean_cbf_gm_basil_native"),
                    ("outputnode.mean_cbf_wm_basil", "inputnode.mean_cbf_wm_basil_native"),
                    ("outputnode.att", "inputnode.att_native"),
                ]),
            ])
            # fmt:on

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # Apply transforms in 1 shot
        std_gereg_wf = init_asl_gestd_trans_wf(
            mem_gb=4,
            omp_nthreads=omp_nthreads,
            spaces=spaces,
            scorescrub=scorescrub,
            basil=basil,
            name="asl_gestd_trans_wf",
        )

        # fmt:off
        workflow.connect([
            (inputnode, std_gereg_wf, [
                ("template", "inputnode.templates"),
                ("anat2std_xfm", "inputnode.anat2std_xfm"),
                ("asl_file", "inputnode.name_source"),
                ("asl_file", "inputnode.asl_file"),
            ]),
            (asl_reg_wf, std_gereg_wf, [
                ("outputnode.aslref_to_t1w_xfm", "inputnode.aslref_to_t1w_xfm"),
            ]),
            (refine_mask, std_gereg_wf, [("out_mask", "inputnode.asl_mask")]),
            (std_gereg_wf, outputnode, [
                ("outputnode.asl_std", "asl_std"),
                ("outputnode.aslref_std", "aslref_std"),
                ("outputnode.asl_mask_std", "asl_mask_std"),
            ]),
            (cbf_compt_wf, std_gereg_wf, [
                ("outputnode.cbf_ts", "inputnode.cbf_ts"),
                ("outputnode.mean_cbf", "inputnode.mean_cbf"),
            ]),
        ])
        # fmt:on

        if scorescrub:
            # fmt:off
            workflow.connect([
                (cbf_compt_wf, std_gereg_wf, [
                    ("outputnode.cbf_ts_score", "inputnode.cbf_ts_score"),
                    ("outputnode.mean_cbf_score", "inputnode.mean_cbf_score"),
                    ("outputnode.mean_cbf_scrub", "inputnode.mean_cbf_scrub"),
                ]),
            ])
            # fmt:on

        if basil:
            # fmt:off
            workflow.connect([
                (cbf_compt_wf, std_gereg_wf, [
                    ("outputnode.mean_cbf_basil", "inputnode.mean_cbf_basil"),
                    ("outputnode.mean_cbf_gm_basil", "inputnode.mean_cbf_gm_basil"),
                    ("outputnode.mean_cbf_wm_basil", "inputnode.mean_cbf_wm_basil"),
                    ("outputnode.att", "inputnode.att"),
                ]),
            ])
            # fmt:on

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # fmt:off
        workflow.connect([
            (std_gereg_wf, compute_cbf_qc_wf, [
                ("outputnode.asl_mask_std", "inputnode.asl_mask_std"),
            ]),
        ])
        # fmt:on

        # asl_derivatives_wf internally parametrizes over snapshotted spaces.
        # fmt:off
        workflow.connect([
            (std_gereg_wf, asl_derivatives_wf, [
                ("outputnode.template", "inputnode.template"),
                ("outputnode.spatial_reference", "inputnode.spatial_reference"),
                ("outputnode.aslref_std", "inputnode.aslref_std"),
                ("outputnode.asl_std", "inputnode.asl_std"),
                ("outputnode.asl_mask_std", "inputnode.asl_mask_std"),
                ("outputnode.cbf_ts_std", "inputnode.cbf_ts_std"),
                ("outputnode.mean_cbf_std", "inputnode.mean_cbf_std"),
            ]),
        ])
        # fmt:on

        if scorescrub:
            # fmt:off
            workflow.connect([
                (std_gereg_wf, asl_derivatives_wf, [
                    ("outputnode.cbf_ts_score_std", "inputnode.cbf_ts_score_std"),
                    ("outputnode.mean_cbf_score_std", "inputnode.mean_cbf_score_std"),
                    ("outputnode.mean_cbf_scrub_std", "inputnode.mean_cbf_scrub_std"),
                ]),
            ])
            # fmt:on

        if basil:
            # fmt:off
            workflow.connect([
                (std_gereg_wf, asl_derivatives_wf, [
                    ("outputnode.mean_cbf_basil_std", "inputnode.mean_cbf_basil_std"),
                    ("outputnode.mean_cbf_gm_basil_std", "inputnode.mean_cbf_gm_basil_std"),
                    ("outputnode.mean_cbf_wm_basil_std", "inputnode.mean_cbf_wm_basil_std"),
                    ("outputnode.att_std", "inputnode.att_std"),
                ]),
            ])
            # fmt:on

    # xform to 'MNI152NLin2009cAsym' is always computed, so this should always be available.
    select_xform_MNI152NLin2009cAsym_to_t1w = pe.Node(
        KeySelect(fields=["std2anat_xfm"], key="MNI152NLin2009cAsym"),
        name="carpetplot_select_std",
        run_without_submitting=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, select_xform_MNI152NLin2009cAsym_to_t1w, [
            ("std2anat_xfm", "std2anat_xfm"),
            ("template", "keys"),
        ]),
    ])
    # fmt:on

    parcellate_cbf_wf = init_parcellate_cbf_wf(
        scorescrub=scorescrub,
        basil=basil,
        name="parcellate_cbf_wf",
    )

    # fmt:off
    workflow.connect([
        (select_xform_MNI152NLin2009cAsym_to_t1w, parcellate_cbf_wf, [
            ("std2anat_xfm", "inputnode.MNI152NLin2009cAsym_to_anat_xform"),
        ]),
        (refine_mask, parcellate_cbf_wf, [("out_mask", "inputnode.asl_mask")]),
        (asl_reg_wf, parcellate_cbf_wf, [
            ("outputnode.itk_t1_to_asl", "inputnode.anat_to_asl_xform"),
        ]),
        (cbf_compt_wf, parcellate_cbf_wf, [("outputnode.mean_cbf", "inputnode.mean_cbf")]),
        (parcellate_cbf_wf, asl_derivatives_wf, [
            ("outputnode.atlas_names", "inputnode.atlas_names"),
            ("outputnode.mean_cbf_parcellated", "inputnode.mean_cbf_parcellated"),
        ]),
    ])
    # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (cbf_compt_wf, parcellate_cbf_wf, [
                ("outputnode.mean_cbf_score", "inputnode.mean_cbf_score"),
                ("outputnode.mean_cbf_scrub", "inputnode.mean_cbf_scrub"),
            ]),
            (parcellate_cbf_wf, asl_derivatives_wf, [
                ("outputnode.mean_cbf_score_parcellated", "inputnode.mean_cbf_score_parcellated"),
                ("outputnode.mean_cbf_scrub_parcellated", "inputnode.mean_cbf_scrub_parcellated"),
            ]),
        ])
        # fmt:on

    if basil:
        # fmt:off
        workflow.connect([
            (cbf_compt_wf, parcellate_cbf_wf, [
                ("outputnode.mean_cbf_basil", "inputnode.mean_cbf_basil"),
                ("outputnode.mean_cbf_gm_basil", "inputnode.mean_cbf_gm_basil"),
            ]),
            (parcellate_cbf_wf, asl_derivatives_wf, [
                ("outputnode.mean_cbf_basil_parcellated", "inputnode.mean_cbf_basil_parcellated"),
                (
                    "outputnode.mean_cbf_gm_basil_parcellated",
                    "inputnode.mean_cbf_gm_basil_parcellated",
                ),
            ]),
        ])
        # fmt:on

    # REPORTING ############################################################
    ds_report_summary = pe.Node(
        DerivativesDataSink(desc="summary", datatype="figures", dismiss_entities=("echo",)),
        name="ds_report_summary",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([(summary, ds_report_summary, [("out_report", "in_file")])])
    # fmt:on

    # Fill-in datasinks of reportlets seen so far
    for node in workflow.list_node_names():
        if node.split(".")[-1].startswith("ds_report"):
            workflow.get_node(node).inputs.base_directory = output_dir
            workflow.get_node(node).inputs.source_file = ref_file

    return workflow
