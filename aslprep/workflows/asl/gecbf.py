# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""CBF-processing workflows for GE data."""
from nipype.interfaces import utility as niu
from nipype.interfaces.fsl import Split as FSLSplit
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
from aslprep.workflows.asl.cbf import init_compute_cbf_ge_wf, init_parcellate_cbf_wf
from aslprep.workflows.asl.ge_utils import init_asl_reference_ge_wf, init_asl_reg_ge_wf
from aslprep.workflows.asl.outputs import init_asl_derivatives_wf
from aslprep.workflows.asl.plotting import init_gecbfplot_wf
from aslprep.workflows.asl.qc import init_compute_cbf_qc_wf
from aslprep.workflows.asl.registration import init_asl_t1_trans_wf
from aslprep.workflows.asl.resampling import init_asl_std_trans_wf


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
    anat_to_template_xfm
        List of transform files, collated with templates
    template_to_anat_xfm
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
    2.  Generate ASL reference image.
        -   Extract averaged, smoothed M0 image and reference image
            (which is generally the M0 image).
    3.  Register ASL to T1w.
    4.  Calculate CBF.
    5.  Apply the ASL-to-T1w transforms to get T1w-space outputs
        (passed along to derivatives workflow).
    6.  Refine the brain mask.
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
                "anat_to_template_xfm",
                "template_to_anat_xfm",
                "template",
            ],
        ),
        name="inputnode",
    )
    inputnode.inputs.asl_file = asl_file
    inputnode.inputs.m0scan = run_data["m0scan"]
    inputnode.inputs.m0scan_metadata = run_data["m0scan_metadata"]

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
    workflow.connect([(inputnode, asl_derivatives_wf, [("asl_file", "inputnode.source_file")])])
    # fmt:on

    # begin workflow
    # Extract averaged, smoothed M0 image and reference image (which is generally the M0 image).
    asl_reference_wf = init_asl_reference_ge_wf(
        metadata=metadata,
        aslcontext=run_data["aslcontext"],
        smooth_kernel=smoothkernel,
        name="asl_reference_ge_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, asl_reference_wf, [
            ("asl_file", "inputnode.asl_file"),
            ("m0scan", "inputnode.m0scan"),
            ("m0scan_metadata", "inputnode.m0scan_metadata"),
        ]),
    ])
    # fmt:on

    # Split 4D ASL file into list of 3D volumes, so that volume-wise transforms (e.g., HMC params)
    # can be applied with other transforms in single shots.
    # This will be useful for GE/non-GE integration.
    asl_split = pe.Node(FSLSplit(dimension="t"), name="asl_split", mem_gb=mem_gb["filesize"] * 3)

    workflow.connect([(inputnode, asl_split, [("asl_file", "in_file")])])

    # Set HMC xforms and fieldwarp to "identity" since neither is performed for GE data.
    # This will be useful as I swap out GE-specific resampling workflows with general ones,
    # which require HMC xforms and a fieldwarp.
    xform_buffer = pe.Node(
        niu.IdentityInterface(
            fields=[
                "hmc_xforms",
                "fieldwarp",  # from sdc
                "epi_brain",  # from sdc
                "epi_mask",  # from sdc
            ],
        ),
        name="xform_buffer",
    )
    xform_buffer.inputs.hmc_xforms = "identity"
    xform_buffer.inputs.fieldwarp = "identity"

    # fmt:off
    workflow.connect([
        (asl_reference_wf, xform_buffer, [
            ("outputnode.ref_image_brain", "epi_brain"),
            ("outputnode.asl_mask", "epi_mask"),
        ]),
    ])
    # fmt:on

    # ASL-to-T1w registration
    asl_reg_wf = init_asl_reg_ge_wf(
        use_bbr=config.workflow.use_bbr,
        asl2t1w_dof=config.workflow.asl2t1w_dof,
        asl2t1w_init=config.workflow.asl2t1w_init,
        name="asl_reg_ge_wf",
        sloppy=False,
        write_report=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, asl_reg_wf, [("t1w_dseg", "inputnode.t1w_dseg")]),
        (asl_reference_wf, asl_reg_wf, [
            ("outputnode.ref_image_brain", "inputnode.ref_asl_brain"),
        ]),
        (t1w_brain, asl_reg_wf, [("out_file", "inputnode.t1w_brain")]),
        (asl_reg_wf, asl_derivatives_wf, [
            ("outputnode.anat_to_aslref_xfm", "inputnode.anat_to_aslref_xfm"),
            ("outputnode.aslref_to_anat_xfm", "inputnode.aslref_to_anat_xfm"),
        ]),
        (asl_reg_wf, summary, [("outputnode.fallback", "fallback")]),
    ])
    # fmt:on

    # Compute CBF from the raw ASL data.
    compute_cbf_wf = init_compute_cbf_ge_wf(
        name_source=asl_file,
        aslcontext=run_data["aslcontext"],
        metadata=metadata,
        scorescrub=scorescrub,
        basil=basil,
        M0Scale=mscale,
        mem_gb=mem_gb["filesize"],
        name="compute_cbf_wf",
    )
    CBF_DERIVS = [
        "cbf_ts",
        "mean_cbf",
        "mean_cbf_basil",
        "mean_cbf_gm_basil",
        "mean_cbf_wm_basil",
        "att",
    ]
    # We don't want mean_cbf_wm_basil for this list and SCORE/SCRUB is blocked for GE data.
    MEAN_CBF_DERIVS = [
        "mean_cbf",
        "mean_cbf_basil",
        "mean_cbf_gm_basil",
    ]

    # fmt:off
    workflow.connect([
        (inputnode, compute_cbf_wf, [
            ("asl_file", "inputnode.asl_file"),
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("t1w_mask", "inputnode.t1w_mask"),
        ]),
        (asl_reference_wf, compute_cbf_wf, [
            ("outputnode.asl_mask", "inputnode.asl_mask"),
            ("outputnode.m0_file", "inputnode.m0_file"),
            ("outputnode.m0tr", "inputnode.m0tr"),
        ]),
        (asl_reg_wf, compute_cbf_wf, [
            ("outputnode.aslref_to_anat_xfm", "inputnode.aslref_to_anat_xfm"),
            ("outputnode.anat_to_aslref_xfm", "inputnode.anat_to_aslref_xfm"),
        ]),
    ])
    # fmt:on

    # Apply ASL registration to T1w
    nonstd_spaces = set(spaces.get_nonstandard())

    asl_t1_trans_wf = init_asl_t1_trans_wf(
        multiecho=False,
        output_t1space=nonstd_spaces.intersection(("T1w", "anat")),
        scorescrub=scorescrub,
        basil=basil,
        generate_reference=False,  # the GE workflow doesn't generate a new reference
        mem_gb=mem_gb["resampled"],
        omp_nthreads=omp_nthreads,
        use_compression=False,
        name="asl_t1_trans_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, asl_t1_trans_wf, [
            ("asl_file", "inputnode.name_source"),
            ("t1w_mask", "inputnode.t1w_mask"),
        ]),
        (t1w_brain, asl_t1_trans_wf, [("out_file", "inputnode.t1w_brain")]),
        (xform_buffer, asl_t1_trans_wf, [
            ("hmc_xforms", "inputnode.hmc_xforms"),
            ("fieldwarp", "inputnode.fieldwarp"),
            ("epi_brain", "inputnode.ref_asl_brain"),
            ("epi_mask", "inputnode.ref_asl_mask"),
        ]),
        # keeping this separate from the top for symmetry with non-GE workflow
        (asl_split, asl_t1_trans_wf, [("out_files", "inputnode.asl_split")]),
        # unused if multiecho, but this is safe
        (asl_reg_wf, asl_t1_trans_wf, [
            ("outputnode.aslref_to_anat_xfm", "inputnode.aslref_to_anat_xfm"),
        ]),
    ])
    # fmt:on

    for cbf_deriv in CBF_DERIVS:
        # fmt:off
        workflow.connect([
            (compute_cbf_wf, asl_t1_trans_wf, [
                (f"outputnode.{cbf_deriv}", f"inputnode.{cbf_deriv}"),
            ]),
            (asl_t1_trans_wf, asl_derivatives_wf, [
                (f"outputnode.{cbf_deriv}_t1", f"inputnode.{cbf_deriv}_t1"),
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
        (asl_reg_wf, refine_mask, [("outputnode.anat_to_aslref_xfm", "transforms")]),
        (asl_reference_wf, refine_mask, [("outputnode.asl_mask", "asl_mask")]),
    ])
    # fmt:on

    # Generate QC metrics
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
            ("outputnode.anat_to_aslref_xfm", "inputnode.anat_to_aslref_xfm"),
        ]),
        (compute_cbf_qc_wf, asl_derivatives_wf, [("outputnode.qc_file", "inputnode.qc_file")]),
        (compute_cbf_qc_wf, summary, [("outputnode.qc_file", "qc_file")]),
    ])
    # fmt:on

    for cbf_deriv in MEAN_CBF_DERIVS:
        # fmt:off
        workflow.connect([
            (compute_cbf_wf, compute_cbf_qc_wf, [
                (f"outputnode.{cbf_deriv}", f"inputnode.{cbf_deriv}"),
            ]),
        ])
        # fmt:on

    # Map native-space outputs to derivatives
    if nonstd_spaces.intersection(("func", "run", "asl", "aslref", "sbref")):
        # fmt:off
        workflow.connect([
            (inputnode, asl_derivatives_wf, [("asl_file", "inputnode.asl_native")]),
            (asl_reference_wf, asl_derivatives_wf, [
                ("outputnode.raw_ref_image", "inputnode.aslref_native"),
            ]),
            (refine_mask, asl_derivatives_wf, [("out_mask", "inputnode.asl_mask_native")]),
        ])
        # fmt:on

        for cbf_deriv in CBF_DERIVS:
            # fmt:off
            workflow.connect([
                (compute_cbf_wf, asl_derivatives_wf, [
                    (f"outputnode.{cbf_deriv}", f"inputnode.{cbf_deriv}_native"),
                ]),
            ])
            # fmt:on

    # Map T1w-space outputs to derivatives.
    # Also warp the final asl mask into T1w space.
    if nonstd_spaces.intersection(("T1w", "anat")):
        from aslprep.interfaces.ants import ApplyTransforms

        aslmask_to_t1w = pe.Node(
            ApplyTransforms(interpolation="MultiLabel"),
            name="aslmask_to_t1w",
            mem_gb=0.1,
        )

        # fmt:off
        workflow.connect([
            (asl_reg_wf, aslmask_to_t1w, [("outputnode.aslref_to_anat_xfm", "transforms")]),
            (asl_t1_trans_wf, aslmask_to_t1w, [("outputnode.asl_mask_t1", "reference_image")]),
            (refine_mask, aslmask_to_t1w, [("out_mask", "input_image")]),
            (aslmask_to_t1w, asl_derivatives_wf, [("output_image", "inputnode.asl_mask_t1")]),
        ])
        # fmt:on

    # Map standard-space outputs to derivatives.
    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # Apply transforms in 1 shot
        asl_std_trans_wf = init_asl_std_trans_wf(
            mem_gb=4,
            omp_nthreads=omp_nthreads,
            spaces=spaces,
            scorescrub=scorescrub,
            basil=basil,
            generate_reference=False,
            name="asl_std_trans_wf",
        )

        # fmt:off
        workflow.connect([
            (inputnode, asl_std_trans_wf, [
                ("template", "inputnode.templates"),
                ("anat_to_template_xfm", "inputnode.anat_to_template_xfm"),
                ("asl_file", "inputnode.name_source"),
            ]),
            (asl_reg_wf, asl_std_trans_wf, [
                ("outputnode.aslref_to_anat_xfm", "inputnode.aslref_to_anat_xfm"),
            ]),
            (refine_mask, asl_std_trans_wf, [("out_mask", "inputnode.asl_mask")]),
            (asl_split, asl_std_trans_wf, [("out_files", "inputnode.asl_split")]),
            (asl_std_trans_wf, compute_cbf_qc_wf, [
                ("outputnode.asl_mask_std", "inputnode.asl_mask_std"),
            ]),
            (asl_std_trans_wf, asl_derivatives_wf, [
                ("outputnode.template", "inputnode.template"),
                ("outputnode.spatial_reference", "inputnode.spatial_reference"),
                ("outputnode.aslref_std", "inputnode.aslref_std"),
                ("outputnode.asl_std", "inputnode.asl_std"),
                ("outputnode.asl_mask_std", "inputnode.asl_mask_std"),
            ]),
        ])
        # fmt:on

        # asl_derivatives_wf internally parametrizes over snapshotted spaces.
        for cbf_deriv in CBF_DERIVS:
            # fmt:off
            workflow.connect([
                (compute_cbf_wf, asl_std_trans_wf, [
                    (f"outputnode.{cbf_deriv}", f"inputnode.{cbf_deriv}"),
                ]),
                (asl_std_trans_wf, asl_derivatives_wf, [
                    (f"outputnode.{cbf_deriv}_std", f"inputnode.{cbf_deriv}_std"),
                ]),
            ])
            # fmt:on

    # Plot CBF outputs.
    plot_cbf_wf = init_gecbfplot_wf(
        basil=basil,
        name="plot_cbf_wf",
    )

    for cbf_deriv in MEAN_CBF_DERIVS:
        # fmt:off
        workflow.connect([
            (compute_cbf_wf, plot_cbf_wf, [
                (f"outputnode.{cbf_deriv}", f"inputnode.{cbf_deriv}"),
            ]),
        ])
        # fmt:on

    # fmt:off
    workflow.connect([
        (asl_reference_wf, plot_cbf_wf, [("outputnode.ref_image_brain", "inputnode.aslref")]),
    ])
    # fmt:on

    # xform to 'MNI152NLin2009cAsym' is always computed, so this should always be available.
    select_xform_MNI152NLin2009cAsym_to_t1w = pe.Node(
        KeySelect(fields=["template_to_anat_xfm"], key="MNI152NLin2009cAsym"),
        name="carpetplot_select_std",
        run_without_submitting=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, select_xform_MNI152NLin2009cAsym_to_t1w, [
            ("template_to_anat_xfm", "template_to_anat_xfm"),
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
            ("template_to_anat_xfm", "inputnode.MNI152NLin2009cAsym_to_anat_xfm"),
        ]),
        (refine_mask, parcellate_cbf_wf, [("out_mask", "inputnode.asl_mask")]),
        (asl_reg_wf, parcellate_cbf_wf, [
            ("outputnode.anat_to_aslref_xfm", "inputnode.anat_to_aslref_xfm"),
        ]),
        (parcellate_cbf_wf, asl_derivatives_wf, [
            ("outputnode.atlas_names", "inputnode.atlas_names"),
        ]),
    ])
    # fmt:on

    for cbf_deriv in MEAN_CBF_DERIVS:
        # fmt:off
        workflow.connect([
            (compute_cbf_wf, parcellate_cbf_wf, [
                (f"outputnode.{cbf_deriv}", f"inputnode.{cbf_deriv}"),
            ]),
            (parcellate_cbf_wf, asl_derivatives_wf, [
                (f"outputnode.{cbf_deriv}_parcellated", f"inputnode.{cbf_deriv}_parcellated"),
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
