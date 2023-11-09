# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Preprocessing workflows for ASL data."""
import nibabel as nb
from fmriprep.workflows.bold.base import get_estimator
from fmriprep.workflows.bold.registration import init_bold_reg_wf
from fmriprep.workflows.bold.resampling import (
    init_bold_fsLR_resampling_wf,
    init_bold_grayords_wf,
    init_bold_surf_wf,
)
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.utils.connections import pop_file

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.interfaces.cbf import RefineMask
from aslprep.interfaces.fsl import Split
from aslprep.interfaces.reports import FunctionalSummary
from aslprep.interfaces.utility import ReduceASLFiles
from aslprep.utils.asl import determine_multi_pld, select_processing_target
from aslprep.utils.bids import collect_run_data
from aslprep.utils.misc import _create_mem_gb, _get_wf_name
from aslprep.workflows.asl.cbf import init_compute_cbf_wf, init_parcellate_cbf_wf
from aslprep.workflows.asl.confounds import init_asl_confounds_wf, init_carpetplot_wf
from aslprep.workflows.asl.hmc import init_asl_hmc_wf
from aslprep.workflows.asl.outputs import init_asl_derivatives_wf
from aslprep.workflows.asl.plotting import init_plot_cbf_wf
from aslprep.workflows.asl.qc import init_compute_cbf_qc_wf
from aslprep.workflows.asl.registration import init_asl_t1_trans_wf
from aslprep.workflows.asl.resampling import (
    init_asl_preproc_trans_wf,
    init_asl_std_trans_wf,
)
from aslprep.workflows.asl.util import init_asl_reference_wf, init_validate_asl_wf


def init_asl_preproc_wf(asl_file, has_fieldmap=False):
    """Perform the functional preprocessing stages of ASLPrep.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep import config
            from aslprep.workflows.asl.base import init_asl_preproc_wf

            with mock_config():
                asl_file = (
                    config.execution.bids_dir / 'sub-01' / 'perf' /
                    'sub-01_asl.nii.gz'
                )
                wf = init_asl_preproc_wf(str(asl_file))

    Parameters
    ----------
    asl_file
        asl series NIfTI file
    has_fieldmap : :obj:`bool`
        Signals the workflow to use inputnode fieldmap files

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
        quality control measures

    Notes
    -----
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
    from niworkflows.interfaces.nibabel import ApplyMask
    from niworkflows.interfaces.reportlets.registration import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )
    from niworkflows.interfaces.utility import KeySelect

    # Have some options handy
    omp_nthreads = config.nipype.omp_nthreads
    spaces = config.workflow.spaces
    dummyvols = config.workflow.dummy_vols
    smooth_kernel = config.workflow.smooth_kernel
    m0_scale = config.workflow.m0_scale
    scorescrub = config.workflow.scorescrub
    basil = config.workflow.basil
    freesurfer = config.workflow.run_reconall
    freesurfer_spaces = spaces.get_fs_spaces()
    project_goodvoxels = config.workflow.project_goodvoxels and config.workflow.cifti_output
    layout = config.execution.layout

    # Take first file (only file, because we don't support multi-echo ASL) as reference
    ref_file = asl_file
    # get original image orientation
    ref_orientation = get_img_orientation(ref_file)

    asl_tlen, mem_gb = _create_mem_gb(ref_file)

    wf_name = _get_wf_name(ref_file)
    config.loggers.workflow.debug(
        (
            'Creating asl processing workflow for "%s" (%.2f GB / %d TRs). '
            "Memory resampled/largemem=%.2f/%.2f GB."
        ),
        ref_file,
        mem_gb["filesize"],
        asl_tlen,
        mem_gb["resampled"],
        mem_gb["largemem"],
    )

    # Collect associated files
    run_data = collect_run_data(layout, ref_file)
    sbref_file = run_data["sbref"]
    metadata = run_data["asl_metadata"].copy()

    is_multi_pld = determine_multi_pld(metadata=metadata)
    if scorescrub and is_multi_pld:
        config.loggers.workflow.warning(
            f"SCORE/SCRUB processing will be disabled for multi-delay {asl_file}"
        )
        scorescrub = False

    # Determine which volumes to use in the pipeline
    processing_target = select_processing_target(aslcontext=run_data["aslcontext"])
    m0type = metadata["M0Type"]

    if has_fieldmap:
        from sdcflows import fieldmaps as fm

        # We may have pruned the estimator collection due to `--ignore fieldmaps`
        estimator_key = [key for key in get_estimator(layout, asl_file) if key in fm._estimators]

        if not estimator_key:
            has_fieldmap = False
            config.loggers.workflow.critical(
                f"None of the available B0 fieldmaps are associated to <{asl_file}>"
            )
        else:
            config.loggers.workflow.info(
                f"Found usable B0-map (fieldmap) estimator(s) <{', '.join(estimator_key)}> "
                f"to correct <{asl_file}> for susceptibility-derived distortions."
            )

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
                "asl_file",
                "asl_metadata",
                "aslcontext",
                "m0scan",
                "m0scan_metadata",
                # anatomical data
                "t1w_preproc",
                "t1w_mask",
                "t1w_dseg",
                "t1w_aseg",
                "t1w_aparc",
                "t1w_tpms",
                "template",
                "anat2std_xfm",
                "std2anat_xfm",
                # Undefined if --fs-no-reconall, but this is safe
                "subjects_dir",
                "subject_id",
                "anat_ribbon",
                "t1w2fsnative_xfm",
                "fsnative2t1w_xfm",
                "surfaces",
                "morphometrics",
                "sphere_reg_fsLR",
            ],
        ),
        name="inputnode",
    )
    inputnode.inputs.asl_file = asl_file
    inputnode.inputs.asl_metadata = metadata
    inputnode.inputs.aslcontext = run_data["aslcontext"]
    inputnode.inputs.m0scan = run_data["m0scan"]
    inputnode.inputs.m0scan_metadata = run_data["m0scan_metadata"]

    if sbref_file is not None:
        from niworkflows.interfaces.header import ValidateImage

        val_sbref = pe.Node(ValidateImage(in_file=sbref_file), name="val_sbref")

    # Generate a brain-masked version of the t1w
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
            registration=("FSL", "FreeSurfer")[freesurfer],
            registration_dof=config.workflow.asl2t1w_dof,
            registration_init=config.workflow.asl2t1w_init,
            pe_direction=metadata.get("PhaseEncodingDirection"),
            tr=metadata.get("RepetitionTime", metadata["RepetitionTimePreparation"]),
            orientation=ref_orientation,
        ),
        name="summary",
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    asl_derivatives_wf = init_asl_derivatives_wf(
        bids_root=layout.root,
        cifti_output=config.workflow.cifti_output,
        freesurfer=freesurfer,
        project_goodvoxels=project_goodvoxels,
        metadata=metadata,
        spaces=spaces,
        is_multi_pld=is_multi_pld,
        scorescrub=scorescrub,
        basil=basil,
        output_confounds=True,
    )
    asl_derivatives_wf.inputs.inputnode.cifti_density = config.workflow.cifti_output
    # fmt:off
    workflow.connect([
        (inputnode, asl_derivatives_wf, [
            ("asl_file", "inputnode.source_file"),
            ("surfaces", "inputnode.surf_files"),
        ]),
    ])
    # fmt:on

    # Validate the ASL file.
    validate_asl_wf = init_validate_asl_wf(name="validate_asl_wf")
    workflow.connect([(inputnode, validate_asl_wf, [("asl_file", "inputnode.asl_file")])])

    # Generate a tentative aslref from the most appropriate available image type in the ASL file
    initial_aslref_wf = init_asl_reference_wf(
        aslcontext=run_data["aslcontext"],
        name="initial_aslref_wf",
    )
    initial_aslref_wf.inputs.inputnode.dummy_scans = config.workflow.dummy_vols

    # fmt:off
    workflow.connect([
        (validate_asl_wf, initial_aslref_wf, [("outputnode.asl_file", "inputnode.asl_file")]),
    ])
    # fmt:on

    if sbref_file is not None:
        workflow.connect([(val_sbref, initial_aslref_wf, [("out_file", "inputnode.sbref_file")])])

    # Drop volumes in the ASL file that won't be used
    # (e.g., precalculated CBF volumes if control-label pairs are available).
    reduce_asl_file = pe.Node(
        ReduceASLFiles(
            processing_target=processing_target,
            metadata=metadata,
        ),
        name="reduce_asl_file",
    )

    # fmt:off
    workflow.connect([
        (inputnode, reduce_asl_file, [("aslcontext", "aslcontext")]),
        (validate_asl_wf, reduce_asl_file, [("outputnode.asl_file", "asl_file")]),
    ])
    # fmt:on

    # Split 4D ASL file into list of 3D volumes, so that volume-wise transforms (e.g., HMC params)
    # can be applied with other transforms in single shots.
    asl_split = pe.Node(Split(dimension="t"), name="asl_split", mem_gb=mem_gb["filesize"] * 3)
    workflow.connect([(reduce_asl_file, asl_split, [("asl_file", "in_file")])])

    # Head motion correction, performed separately for each retained image type
    asl_hmc_wf = init_asl_hmc_wf(
        processing_target=processing_target,
        m0type=m0type,
        mem_gb=mem_gb["filesize"],
        omp_nthreads=omp_nthreads,
        name="asl_hmc_wf",
    )

    # fmt:off
    workflow.connect([
        (reduce_asl_file, asl_hmc_wf, [
            ("asl_file", "inputnode.asl_file"),
            ("aslcontext", "inputnode.aslcontext"),
        ]),
        (initial_aslref_wf, asl_hmc_wf, [("outputnode.raw_ref_image", "inputnode.raw_ref_image")]),
    ])
    # fmt:on

    # Calculate ASL-to-T1 registration
    asl_reg_wf = init_bold_reg_wf(
        bold2t1w_dof=config.workflow.asl2t1w_dof,
        bold2t1w_init=config.workflow.asl2t1w_init,
        freesurfer=freesurfer,
        mem_gb=mem_gb["resampled"],
        name="asl_reg_wf",
        omp_nthreads=omp_nthreads,
        sloppy=config.execution.debug,
        use_bbr=config.workflow.use_bbr,
        use_compression=False,
    )

    # fmt:off
    workflow.connect([
        (inputnode, asl_reg_wf, [
            ("t1w_dseg", "inputnode.t1w_dseg"),
            # Undefined if --fs-no-reconall, but this is safe
            ("subjects_dir", "inputnode.subjects_dir"),
            ("subject_id", "inputnode.subject_id"),
            ("fsnative2t1w_xfm", "inputnode.fsnative2t1w_xfm"),
        ]),
        (t1w_brain, asl_reg_wf, [("out_file", "inputnode.t1w_brain")]),
        (asl_reg_wf, summary, [("outputnode.fallback", "fallback")]),
    ])
    # fmt:on

    # Apply ASL to T1w registration
    asl_t1_trans_wf = init_asl_t1_trans_wf(
        freesurfer=freesurfer,
        is_multi_pld=is_multi_pld,
        scorescrub=scorescrub,
        basil=basil,
        mem_gb=mem_gb["resampled"],
        omp_nthreads=omp_nthreads,
        generate_reference=True,
        use_compression=False,
        name="asl_t1_trans_wf",
    )
    asl_t1_trans_wf.inputs.inputnode.fieldwarp = "identity"

    # fmt:off
    workflow.connect([
        (inputnode, asl_t1_trans_wf, [
            ("asl_file", "inputnode.name_source"),
            ("t1w_mask", "inputnode.t1w_mask"),
            ("t1w_aseg", "inputnode.t1w_aseg"),
            ("t1w_aparc", "inputnode.t1w_aparc"),
        ]),
        (reduce_asl_file, asl_t1_trans_wf, [("aslcontext", "inputnode.aslcontext")]),
        (t1w_brain, asl_t1_trans_wf, [("out_file", "inputnode.t1w_brain")]),
        (asl_split, asl_t1_trans_wf, [("out_files", "inputnode.asl_split")]),
        (asl_hmc_wf, asl_t1_trans_wf, [("outputnode.xforms", "inputnode.hmc_xforms")]),
        (asl_reg_wf, asl_t1_trans_wf, [("outputnode.itk_bold_to_t1", "inputnode.itk_bold_to_t1")]),
    ])
    # fmt:on

    # get confounds
    asl_confounds_wf = init_asl_confounds_wf(mem_gb=mem_gb["largemem"], name="asl_confounds_wf")
    asl_confounds_wf.get_node("inputnode").inputs.t1_transform_flags = [False]

    # fmt:off
    workflow.connect([
        (inputnode, asl_confounds_wf, [
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("t1w_mask", "inputnode.t1w_mask"),
        ]),
        (asl_hmc_wf, asl_confounds_wf, [
            ("outputnode.movpar_file", "inputnode.movpar_file"),
            ("outputnode.rmsd_file", "inputnode.rmsd_file"),
        ]),
        (asl_reg_wf, asl_confounds_wf, [
            ("outputnode.itk_t1_to_bold", "inputnode.anat_to_aslref_xfm"),
        ]),
        (initial_aslref_wf, asl_confounds_wf, [("outputnode.skip_vols", "inputnode.skip_vols")]),
        (asl_confounds_wf, asl_derivatives_wf, [
            ("outputnode.confounds_file", "inputnode.confounds"),
            ("outputnode.confounds_metadata", "inputnode.confounds_metadata"),
        ]),
        (asl_confounds_wf, summary, [("outputnode.confounds_file", "confounds_file")]),
    ])
    # fmt:on

    # Compute CBF from ASLRef-space ASL data.
    compute_cbf_wf = init_compute_cbf_wf(
        name_source=asl_file,
        processing_target=processing_target,
        dummy_vols=dummyvols,
        m0_scale=m0_scale,
        scorescrub=scorescrub,
        basil=basil,
        smooth_kernel=smooth_kernel,
        metadata=metadata,
        name="compute_cbf_wf",
    )
    cbf_derivs = ["mean_cbf"]
    mean_cbf_derivs = ["mean_cbf"]

    if is_multi_pld:
        cbf_derivs += ["att"]
    else:
        cbf_derivs += ["cbf_ts"]

    if scorescrub:
        cbf_derivs += [
            "cbf_ts_score",
            "mean_cbf_score",
            "mean_cbf_scrub",
        ]
        mean_cbf_derivs += [
            "mean_cbf_score",
            "mean_cbf_scrub",
        ]

    if basil:
        cbf_derivs += [
            "mean_cbf_basil",
            "mean_cbf_gm_basil",
            "mean_cbf_wm_basil",
            "att_basil",
        ]
        # We don't want mean_cbf_wm_basil for this list.
        mean_cbf_derivs += [
            "mean_cbf_basil",
            "mean_cbf_gm_basil",
        ]

    # fmt:off
    workflow.connect([
        (inputnode, compute_cbf_wf, [
            ("t1w_mask", "inputnode.t1w_mask"),
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("m0scan", "inputnode.m0scan"),
            ("m0scan_metadata", "inputnode.m0scan_metadata"),
        ]),
        (reduce_asl_file, compute_cbf_wf, [("aslcontext", "inputnode.aslcontext")]),
        (asl_reg_wf, compute_cbf_wf, [
            ("outputnode.itk_t1_to_bold", "inputnode.anat_to_aslref_xfm"),
            ("outputnode.itk_bold_to_t1", "inputnode.aslref_to_anat_xfm"),
        ]),
        (asl_reg_wf, asl_derivatives_wf, [
            ("outputnode.itk_t1_to_bold", "inputnode.anat_to_aslref_xfm"),
            ("outputnode.itk_bold_to_t1", "inputnode.aslref_to_anat_xfm"),
        ]),
    ])
    # fmt:on

    asl_final = pe.Node(
        niu.IdentityInterface(fields=["asl", "aslref", "mask"]),
        name="asl_final",
    )
    # fmt:off
    workflow.connect([
        (asl_final, asl_reg_wf, [("aslref", "inputnode.ref_bold_brain")]),
        (asl_final, asl_t1_trans_wf, [
            ("mask", "inputnode.ref_asl_mask"),
            ("boldref", "inputnode.ref_asl_brain"),
        ]),
    ])
    # fmt:on

    # Generate a tentative aslref from the most appropriate available image type in the ASL file
    final_aslref_wf = init_asl_reference_wf(
        aslcontext=run_data["aslcontext"],
        name="final_aslref_wf",
    )
    final_aslref_wf.__desc__ = None

    # fmt:off
    workflow.connect([
        (initial_aslref_wf, final_aslref_wf, [("outputnode.skip_vols", "inputnode.dummy_scans")]),
        (final_aslref_wf, asl_final, [
            ("outputnode.ref_image", "aslref"),
            ("outputnode.asl_mask", "mask"),
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
        (asl_reg_wf, refine_mask, [("outputnode.itk_t1_to_bold", "transforms")]),
    ])
    # fmt:on

    # Generate QC metrics
    compute_cbf_qc_wf = init_compute_cbf_qc_wf(
        is_ge=False,
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
        (asl_hmc_wf, compute_cbf_qc_wf, [("outputnode.rmsd_file", "inputnode.rmsd_file")]),
        (asl_reg_wf, compute_cbf_qc_wf, [
            ("outputnode.itk_t1_to_bold", "inputnode.anat_to_aslref_xfm"),
        ]),
        (asl_confounds_wf, compute_cbf_qc_wf, [
            ("outputnode.confounds_file", "inputnode.confounds_file"),
        ]),
        (compute_cbf_qc_wf, asl_derivatives_wf, [("outputnode.qc_file", "inputnode.qc_file")]),
        (compute_cbf_qc_wf, summary, [("outputnode.qc_file", "qc_file")]),
    ])
    # fmt:on

    for cbf_deriv in mean_cbf_derivs:
        # fmt:off
        workflow.connect([
            (compute_cbf_wf, compute_cbf_qc_wf, [
                (f"outputnode.{cbf_deriv}", f"inputnode.{cbf_deriv}"),
            ]),
        ])
        # fmt:on

    # Map final ASL mask into T1w space (if required)
    nonstd_spaces = set(spaces.get_nonstandard())
    if nonstd_spaces.intersection(("T1w", "anat")):
        from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms

        aslmask_to_t1w = pe.Node(
            ApplyTransforms(interpolation="MultiLabel"),
            name="aslmask_to_t1w",
            mem_gb=0.1,
        )
        # fmt:off
        workflow.connect([
            (asl_reg_wf, aslmask_to_t1w, [("outputnode.itk_bold_to_t1", "transforms")]),
            (asl_t1_trans_wf, aslmask_to_t1w, [("outputnode.asl_mask_t1", "reference_image")]),
            (asl_final, aslmask_to_t1w, [("mask", "input_image")]),
        ])
        # fmt:on

        # Transform CBF derivatives to T1 space.
        # These derivatives have already undergone HMC+SDC,
        # so we only need to apply the ASLRef-to-T1w transform.
        for cbf_deriv in cbf_derivs:
            kwargs = {}
            if cbf_deriv in ["cbf_ts", "cbf_ts_score"]:
                kwargs["dimension"] = 3

            cbf_to_t1w = pe.Node(
                ApplyTransforms(
                    interpolation="LanczosWindowedSinc",
                    float=True,
                    input_image_type=3,
                    **kwargs,
                ),
                name=f"warp_{cbf_deriv}_to_t1w",
                mem_gb=mem_gb * 3 * omp_nthreads,
                n_procs=omp_nthreads,
            )
            # fmt:off
            workflow.connect([
                (asl_reg_wf, cbf_to_t1w, [("outputnode.itk_bold_to_t1", "transforms")]),
                (asl_t1_trans_wf, cbf_to_t1w, [("outputnode.asl_mask_t1", "reference_image")]),
                (asl_final, cbf_to_t1w, [(cbf_deriv, "input_image")]),
                (cbf_to_t1w, asl_derivatives_wf, [("output_image", f"inputnode.{cbf_deriv}_t1")]),
            ])
            # fmt:on

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # Apply transforms in 1 shot
        asl_std_trans_wf = init_asl_std_trans_wf(
            freesurfer=freesurfer,
            mem_gb=mem_gb["resampled"],
            omp_nthreads=omp_nthreads,
            spaces=spaces,
            is_multi_pld=is_multi_pld,
            scorescrub=scorescrub,
            basil=basil,
            generate_reference=True,
            name="asl_std_trans_wf",
            use_compression=not config.execution.low_mem,
        )
        asl_std_trans_wf.inputs.inputnode.fieldwarp = "identity"

        # fmt:off
        workflow.connect([
            (inputnode, asl_std_trans_wf, [
                ("template", "inputnode.templates"),
                ("anat2std_xfm", "inputnode.anat2std_xfm"),
                ("asl_file", "inputnode.name_source"),
                ("t1w_aseg", "inputnode.asl_aseg"),
                ("t1w_aparc", "inputnode.asl_aparc"),
            ]),
            (asl_final, asl_std_trans_wf, [("mask", "inputnode.asl_mask")]),
            (asl_reg_wf, asl_std_trans_wf, [
                ("outputnode.itk_bold_to_t1", "inputnode.aslref_to_anat_xfm"),
            ]),
        ])
        # fmt:on

        if freesurfer:
            # fmt:off
            workflow.connect([
                (asl_std_trans_wf, asl_derivatives_wf, [
                    ("outputnode.asl_aseg_std", "inputnode.asl_aseg_std"),
                    ("outputnode.asl_aparc_std", "inputnode.asl_aparc_std"),
                ]),
            ])
            # fmt:on

        # fmt:off
        workflow.connect([
            (asl_split, asl_std_trans_wf, [("out_files", "inputnode.asl_split")]),
            (asl_hmc_wf, asl_std_trans_wf, [("outputnode.xforms", "inputnode.hmc_xforms")]),
        ])
        # fmt:on

        # func_derivatives_wf internally parametrizes over snapshotted spaces.
        # fmt:off
        workflow.connect([
            (asl_std_trans_wf, asl_derivatives_wf, [
                ("outputnode.template", "inputnode.template"),
                ("outputnode.spatial_reference", "inputnode.spatial_reference"),
                ("outputnode.aslref_std", "inputnode.aslref_std"),
                ("outputnode.asl_std", "inputnode.asl_std"),
                ("outputnode.asl_mask_std", "inputnode.asl_mask_std"),
            ]),
        ])
        # fmt:on

    # SURFACES ##################################################################################
    # Freesurfer
    if freesurfer and freesurfer_spaces:
        config.loggers.workflow.debug("Creating BOLD surface-sampling workflow.")
        asl_surf_wf = init_bold_surf_wf(
            mem_gb=mem_gb["resampled"],
            surface_spaces=freesurfer_spaces,
            medial_surface_nan=config.workflow.medial_surface_nan,
            name="asl_surf_wf",
        )
        # fmt:off
        workflow.connect([
            (inputnode, asl_surf_wf, [
                ("subjects_dir", "inputnode.subjects_dir"),
                ("subject_id", "inputnode.subject_id"),
                ("t1w2fsnative_xfm", "inputnode.t1w2fsnative_xfm"),
            ]),
            (asl_t1_trans_wf, asl_surf_wf, [("outputnode.asl_t1", "inputnode.source_file")]),
            (asl_surf_wf, asl_derivatives_wf, [("outputnode.target", "inputnode.surf_refs")]),
        ])
        # fmt:on

    # CIFTI output
    if config.workflow.cifti_output:
        asl_fsLR_resampling_wf = init_bold_fsLR_resampling_wf(
            estimate_goodvoxels=project_goodvoxels,
            grayord_density=config.workflow.cifti_output,
            omp_nthreads=omp_nthreads,
            mem_gb=mem_gb["resampled"],
        )
        # fmt:off
        workflow.connect([
            (inputnode, asl_fsLR_resampling_wf, [
                ("surfaces", "inputnode.surfaces"),
                ("morphometrics", "inputnode.morphometrics"),
                ("sphere_reg_fsLR", "inputnode.sphere_reg_fsLR"),
                ("anat_ribbon", "inputnode.anat_ribbon"),
            ]),
            (asl_t1_trans_wf, asl_fsLR_resampling_wf, [
                ("outputnode.asl_t1", "inputnode.bold_file"),
            ]),
            (asl_fsLR_resampling_wf, asl_derivatives_wf, [
                ("outputnode.goodvoxels_mask", "inputnode.goodvoxels_mask"),
            ]),
        ])
        # fmt:on

        asl_grayords_wf = init_bold_grayords_wf(
            grayord_density=config.workflow.cifti_output,
            mem_gb=mem_gb["resampled"],
            repetition_time=metadata["RepetitionTime"],
        )

        # fmt:off
        workflow.connect([
            (asl_std_trans_wf, asl_grayords_wf, [
                ("outputnode.asl_std", "inputnode.bold_std"),
                ("outputnode.spatial_reference", "inputnode.spatial_reference"),
            ]),
            (asl_fsLR_resampling_wf, asl_grayords_wf, [
                ("outputnode.asl_fsLR", "inputnode.bold_fsLR"),
            ]),
            (asl_grayords_wf, asl_derivatives_wf, [
                ("outputnode.cifti_bold", "inputnode.asl_cifti"),
                ("outputnode.cifti_metadata", "inputnode.cifti_metadata"),
            ]),
        ])
        # fmt:on

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        carpetplot_wf = init_carpetplot_wf(
            mem_gb=mem_gb["resampled"],
            metadata=metadata,
            cifti_output=config.workflow.cifti_output,
            name="carpetplot_wf",
        )

        # Xform to "MNI152NLin2009cAsym" is always computed.
        carpetplot_select_std = pe.Node(
            KeySelect(fields=["std2anat_xfm"], key="MNI152NLin2009cAsym"),
            name="carpetplot_select_std",
            run_without_submitting=True,
        )

        if config.workflow.cifti_output:
            # fmt:off
            workflow.connect([
                (asl_grayords_wf, carpetplot_wf, [
                    ("outputnode.cifti_bold", "inputnode.cifti_bold"),
                ]),
            ])
            # fmt:on

        def _last(inlist):
            return inlist[-1]

        # fmt:off
        workflow.connect([
            (initial_aslref_wf, carpetplot_wf, [
                ("outputnode.skip_vols", "inputnode.dummy_scans"),
            ]),
            (inputnode, carpetplot_select_std, [
                ("std2anat_xfm", "std2anat_xfm"),
                ("template", "keys"),
            ]),
            (carpetplot_select_std, carpetplot_wf, [("std2anat_xfm", "inputnode.std2anat_xfm")]),
            (asl_final, carpetplot_wf, [
                ("asl", "inputnode.bold"),
                ("mask", "inputnode.bold_mask"),
            ]),
            (asl_reg_wf, carpetplot_wf, [
                ("outputnode.itk_t1_to_bold", "inputnode.t1_bold_xform"),
            ]),
            (asl_confounds_wf, carpetplot_wf, [
                ("outputnode.confounds_file", "inputnode.confounds_file"),
                ("outputnode.crown_mask", "inputnode.crown_mask"),
                (("outputnode.acompcor_masks", _last), "inputnode.acompcor_mask"),
            ]),
        ])
        # fmt:on

    # Plot CBF outputs.
    plot_cbf_wf = init_plot_cbf_wf(
        metadata=metadata,
        plot_timeseries=not is_multi_pld,
        scorescrub=scorescrub,
        basil=basil,
        name="plot_cbf_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, plot_cbf_wf, [("t1w_dseg", "inputnode.t1w_dseg")]),
        (carpetplot_select_std, plot_cbf_wf, [("std2anat_xfm", "inputnode.std2anat_xfm")]),
        (compute_cbf_wf, plot_cbf_wf, [
            ("outputnode.score_outlier_index", "inputnode.score_outlier_index"),
        ]),
        (initial_aslref_wf, plot_cbf_wf, [("outputnode.ref_image_brain", "inputnode.aslref")]),
        (asl_reg_wf, plot_cbf_wf, [
            ("outputnode.itk_t1_to_bold", "inputnode.anat_to_aslref_xfm"),
        ]),
        (refine_mask, plot_cbf_wf, [("out_mask", "inputnode.asl_mask")]),
        (asl_confounds_wf, plot_cbf_wf, [
            ("outputnode.confounds_file", "inputnode.confounds_file"),
        ]),
    ])
    # fmt:on

    for cbf_deriv in cbf_derivs:
        # fmt:off
        workflow.connect([
            (compute_cbf_wf, plot_cbf_wf, [
                (f"outputnode.{cbf_deriv}", f"inputnode.{cbf_deriv}"),
            ]),
        ])
        # fmt:on

    parcellate_cbf_wf = init_parcellate_cbf_wf(
        basil=basil,
        scorescrub=scorescrub,
        name="parcellate_cbf_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, parcellate_cbf_wf, [("asl_file", "inputnode.source_file")]),
        (carpetplot_select_std, parcellate_cbf_wf, [
            ("std2anat_xfm", "inputnode.MNI152NLin2009cAsym_to_anat_xfm"),
        ]),
        (asl_final, parcellate_cbf_wf, [("mask", "inputnode.asl_mask")]),
        (asl_reg_wf, parcellate_cbf_wf, [
            ("outputnode.itk_t1_to_bold", "inputnode.anat_to_aslref_xfm"),
        ]),
        (parcellate_cbf_wf, asl_derivatives_wf, [
            ("outputnode.atlas_names", "inputnode.atlas_names"),
        ]),
    ])
    # fmt:on

    for cbf_deriv in mean_cbf_derivs:
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

    ds_report_validation = pe.Node(
        DerivativesDataSink(
            base_directory=config.execution.aslprep_dir,
            desc="validation",
            datatype="figures",
            dismiss_entities=("echo",),
        ),
        name="ds_report_validation",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([
        (summary, ds_report_summary, [("out_report", "in_file")]),
        (initial_aslref_wf, ds_report_validation, [("outputnode.validation_report", "in_file")]),
    ])
    # fmt:on

    # Fill-in datasinks of reportlets seen so far
    for node in workflow.list_node_names():
        if node.split(".")[-1].startswith("ds_report"):
            workflow.get_node(node).inputs.base_directory = config.execution.aslprep_dir
            workflow.get_node(node).inputs.source_file = ref_file

    if not has_fieldmap:
        # Finalize workflow without SDC connections
        summary.inputs.distortion_correction = "None"

        # Resample in native space in just one shot
        asl_asl_trans_wf = init_asl_preproc_trans_wf(
            mem_gb=mem_gb["resampled"],
            omp_nthreads=omp_nthreads,
            use_compression=not config.execution.low_mem,
            use_fieldwarp=False,
            name="asl_asl_trans_wf",
        )
        asl_asl_trans_wf.inputs.inputnode.fieldwarp = "identity"

        # fmt:off
        workflow.connect([
            # Connect asl_asl_trans_wf
            (inputnode, asl_asl_trans_wf, [("name_source", "inputnode.name_source")]),
            (asl_split, asl_asl_trans_wf, [("out_files", "inputnode.asl_file")]),
            (asl_hmc_wf, asl_asl_trans_wf, [("outputnode.xforms", "inputnode.hmc_xforms")]),
            (asl_asl_trans_wf, asl_final, [("outputnode.asl", "asl")]),
            (asl_asl_trans_wf, final_aslref_wf, [("outputnode.asl", "inputnode.asl_file")]),
        ])
        # fmt:on

        return workflow

    from niworkflows.interfaces.utility import KeySelect
    from sdcflows.workflows.apply.correction import init_unwarp_wf
    from sdcflows.workflows.apply.registration import init_coeff2epi_wf

    coeff2epi_wf = init_coeff2epi_wf(
        debug="fieldmaps" in config.execution.debug,
        omp_nthreads=config.nipype.omp_nthreads,
        sloppy=config.execution.sloppy,
        write_coeff=True,
    )
    unwarp_wf = init_unwarp_wf(
        free_mem=config.environment.free_mem,
        debug="fieldmaps" in config.execution.debug,
        omp_nthreads=config.nipype.omp_nthreads,
    )
    unwarp_wf.inputs.inputnode.metadata = metadata

    output_select = pe.Node(
        KeySelect(fields=["fmap", "fmap_ref", "fmap_coeff", "fmap_mask", "sdc_method"]),
        name="output_select",
        run_without_submitting=True,
    )
    output_select.inputs.key = estimator_key[0]
    if len(estimator_key) > 1:
        config.loggers.workflow.warning(
            f"Several fieldmaps <{', '.join(estimator_key)}> are "
            f"'IntendedFor' <{asl_file}>, using {estimator_key[0]}"
        )

    sdc_report = pe.Node(
        SimpleBeforeAfter(
            before_label="Distorted",
            after_label="Corrected",
            dismiss_affine=True,
        ),
        name="sdc_report",
        mem_gb=0.1,
    )

    ds_report_sdc = pe.Node(
        DerivativesDataSink(
            base_directory=config.execution.aslprep_dir,
            desc="sdc",
            suffix="asl",
            datatype="figures",
            dismiss_entities=("echo",),
        ),
        name="ds_report_sdc",
        run_without_submitting=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, output_select, [
            ("fmap", "fmap"),
            ("fmap_ref", "fmap_ref"),
            ("fmap_coeff", "fmap_coeff"),
            ("fmap_mask", "fmap_mask"),
            ("sdc_method", "sdc_method"),
            ("fmap_id", "keys"),
        ]),
        (output_select, coeff2epi_wf, [
            ("fmap_ref", "inputnode.fmap_ref"),
            ("fmap_coeff", "inputnode.fmap_coeff"),
            ("fmap_mask", "inputnode.fmap_mask")]),
        (output_select, summary, [("sdc_method", "distortion_correction")]),
        (initial_aslref_wf, coeff2epi_wf, [
            ("outputnode.ref_image", "inputnode.target_ref"),
            ("outputnode.asl_mask", "inputnode.target_mask")]),
        (initial_aslref_wf, unwarp_wf, [("outputnode.ref_image", "inputnode.distorted_ref")]),
        (coeff2epi_wf, unwarp_wf, [("outputnode.fmap_coeff", "inputnode.fmap_coeff")]),
        (asl_hmc_wf, unwarp_wf, [("outputnode.xforms", "inputnode.hmc_xforms")]),
        (initial_aslref_wf, sdc_report, [("outputnode.ref_image", "before")]),
        (asl_split, unwarp_wf, [("out_files", "inputnode.distorted")]),
        (final_aslref_wf, sdc_report, [
            ("outputnode.ref_image", "after"),
            ("outputnode.asl_mask", "wm_seg")]),
        (inputnode, ds_report_sdc, [("asl_file", "source_file")]),
        (sdc_report, ds_report_sdc, [("out_report", "in_file")]),

    ])
    # fmt:on

    if "fieldmaps" in config.execution.debug:
        # Generate additional reportlets to assess SDC
        from sdcflows.interfaces.reportlets import FieldmapReportlet

        # First, one for checking the co-registration between fieldmap and EPI
        sdc_coreg_report = pe.Node(
            SimpleBeforeAfter(
                before_label="Distorted target",
                after_label="Fieldmap ref.",
            ),
            name="sdc_coreg_report",
            mem_gb=0.1,
        )
        ds_report_sdc_coreg = pe.Node(
            DerivativesDataSink(
                base_directory=config.execution.aslprep_dir,
                datatype="figures",
                desc="fmapCoreg",
                dismiss_entities=("echo",),
                suffix="asl",
            ),
            name="ds_report_sdc_coreg",
            run_without_submitting=True,
        )

        # Second, showing the fieldmap reconstructed from coefficients in the EPI space
        fmap_report = pe.Node(FieldmapReportlet(), "fmap_report")

        ds_fmap_report = pe.Node(
            DerivativesDataSink(
                base_directory=config.execution.aslprep_dir,
                datatype="figures",
                desc="fieldmap",
                dismiss_entities=("echo",),
                suffix="asl",
            ),
            name="ds_fmap_report",
            run_without_submitting=True,
        )

        # fmt:off
        workflow.connect([
            (initial_aslref_wf, sdc_coreg_report, [("outputnode.ref_image", "before")]),
            (coeff2epi_wf, sdc_coreg_report, [("coregister.inverse_warped_image", "after")]),
            (final_aslref_wf, sdc_coreg_report, [("outputnode.asl_mask", "wm_seg")]),
            (inputnode, ds_report_sdc_coreg, [("asl_file", "source_file")]),
            (sdc_coreg_report, ds_report_sdc_coreg, [("out_report", "in_file")]),
            (unwarp_wf, fmap_report, [(("outputnode.fieldmap", pop_file), "fieldmap")]),
            (coeff2epi_wf, fmap_report, [("coregister.inverse_warped_image", "reference")]),
            (final_aslref_wf, fmap_report, [("outputnode.asl_mask", "mask")]),
            (fmap_report, ds_fmap_report, [("out_report", "in_file")]),
            (inputnode, ds_fmap_report, [("asl_file", "source_file")]),
        ])
        # fmt:on

    # fmt:off
    workflow.connect([
        (unwarp_wf, asl_final, [("outputnode.corrected", "asl")]),
        # remaining workflow connections
        (unwarp_wf, final_aslref_wf, [("outputnode.corrected", "inputnode.asl_file")]),
        (unwarp_wf, asl_t1_trans_wf, [
            # TEMPORARY: For the moment we can't use frame-wise fieldmaps
            (("outputnode.fieldwarp_ref", pop_file), "inputnode.fieldwarp"),
        ]),
        (unwarp_wf, asl_std_trans_wf, [
            # TEMPORARY: For the moment we can't use frame-wise fieldmaps
            (("outputnode.fieldwarp_ref", pop_file), "inputnode.fieldwarp"),
        ]),
    ])
    # fmt:on

    return workflow


def get_img_orientation(imgf):
    """Return the image orientation as a string."""
    img = nb.load(imgf)
    return "".join(nb.aff2axcodes(img.affine))
