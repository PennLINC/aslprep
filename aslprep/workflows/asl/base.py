# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Preprocessing workflows for ASL data."""
from nipype.interfaces import utility as niu
from nipype.interfaces.fsl import Split as FSLSplit
from nipype.pipeline import engine as pe

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.interfaces.cbf_computation import RefineMask
from aslprep.interfaces.reports import FunctionalSummary
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.func.util import init_asl_reference_wf
from aslprep.niworkflows.interfaces.nibabel import ApplyMask
from aslprep.niworkflows.interfaces.utility import KeySelect
from aslprep.sdcflows.workflows.base import fieldmap_wrangler, init_sdc_estimate_wf
from aslprep.utils.bids import collect_run_data
from aslprep.utils.meepi import combine_meepi_source
from aslprep.utils.misc import _create_mem_gb, _get_series_len, _get_wf_name
from aslprep.workflows.asl.cbf import init_cbf_compt_wf, init_parcellate_cbf_wf
from aslprep.workflows.asl.confounds import init_asl_confs_wf, init_carpetplot_wf
from aslprep.workflows.asl.hmc import init_asl_hmc_wf
from aslprep.workflows.asl.outputs import init_asl_derivatives_wf
from aslprep.workflows.asl.plotting import init_cbfplot_wf
from aslprep.workflows.asl.qc import init_compute_cbf_qc_wf
from aslprep.workflows.asl.registration import init_asl_reg_wf, init_asl_t1_trans_wf
from aslprep.workflows.asl.resampling import (
    init_asl_preproc_trans_wf,
    init_asl_std_trans_wf,
)
from aslprep.workflows.asl.stc import init_asl_stc_wf
from aslprep.workflows.asl.t2s import init_asl_t2s_wf


def init_asl_preproc_wf(asl_file):
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
                    'sub-01_task-restEyesOpen_asl.nii.gz'
                )
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
    anat_to_template_xfm
        List of transform files, collated with templates
    template_to_anat_xfm
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
        quality control meausres

    Notes
    -----
    1.  Brain-mask T1w.
    2.  Generate initial ASL reference image.
        I think this is just the BOLD reference workflow from niworkflows,
        without any modifications for ASL data. I could be missing something.
        -   Feeds into the aslbuffer if STC isn't run.
        -   Not in GE workflow.
    3.  Slice-timing correction (can be skipped).
        -   Feeds into the aslbuffer if STC is run.
        -   Not in GE workflow.
    4.  Motion correction.
        -   Outputs the HMC transform and the motion parameters.
        -   Not in GE workflow.
    5.  Susceptibility distortion correction.
        -   Outputs the SDC transform.
        -   Not in GE workflow.
    6.  If data are multi-echo (which they can't be), skullstrip the ASL data and perform optimal
        combination.
        -   Not in GE workflow.
    7.  Register ASL to T1w.
    8.  Apply the ASL-to-T1w transforms to get T1w-space outputs
        (passed along to derivatives workflow).
    9.  Apply the ASL-to-ASLref transforms to get native-space outputs.
    10. Calculate confounds.
        -   Not in GE workflow.
    11. Calculate CBF.
    12. Refine the brain mask.
    13. Generate distortion correction report.
        -   Not in GE workflow.
    14. Apply ASL-to-template transforms to get template-space outputs.
    15. Calculate CBF QC measures.
    16. Generate CBF plots.
    17. Generate carpet plots.
        -   Not in GE workflow.
    18. Parcellate CBF results.
    """
    mem_gb = {"filesize": 1, "resampled": 1, "largemem": 1}
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

    ref_file = asl_file
    if multiecho:
        tes = [layout.get_metadata(echo)["EchoTime"] for echo in asl_file]
        ref_file = dict(zip(tes, asl_file))[min(tes)]

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
    run_data = collect_run_data(layout, ref_file, multiecho=multiecho)
    sbref_file = run_data["sbref"]
    metadata = run_data["asl_metadata"].copy()

    # Find fieldmaps. Options: (phase1|phase2|phasediff|epi|fieldmap|syn)
    fmaps = None
    if "fieldmaps" not in config.workflow.ignore:
        fmaps = fieldmap_wrangler(
            layout,
            ref_file,
            use_syn=config.workflow.use_syn_sdc,
            force_syn=config.workflow.force_syn,
        )
    elif config.workflow.use_syn_sdc or config.workflow.force_syn:
        # If fieldmaps are not enabled, activate SyN-SDC in unforced (False) mode
        fmaps = {"syn": False}

    n_vols = _get_series_len(ref_file)
    MIN_VOLS_FOR_STC = 4
    run_stc = (
        bool(metadata.get("SliceTiming"))
        and "slicetiming" not in config.workflow.ignore
        and (n_vols > MIN_VOLS_FOR_STC)
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

    if sbref_file is not None:
        from aslprep.niworkflows.interfaces.images import ValidateImage

        val_sbref = pe.Node(ValidateImage(in_file=sbref_file), name="val_sbref")

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
            slice_timing="TooShort" if n_vols <= MIN_VOLS_FOR_STC else run_stc,
            registration=("FSL"),
            registration_dof=config.workflow.asl2t1w_dof,
            registration_init=config.workflow.asl2t1w_init,
            pe_direction=metadata.get("PhaseEncodingDirection"),
            tr=metadata.get("RepetitionTime", metadata["RepetitionTimePreparation"]),
        ),
        name="summary",
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    asl_derivatives_wf = init_asl_derivatives_wf(
        bids_root=layout.root,
        metadata=metadata,
        output_dir=output_dir,
        spaces=spaces,
        scorescrub=scorescrub,
        basil=basil,
        output_confounds=True,
    )

    # Generate a tentative aslref
    # XXX: Does the niworkflows reference workflow make sense for ASL data,
    # where volumes can be control, label, M0, deltam, or CBF?
    # It seems like the different contrasts could impact the reference image estimation.
    # For example, should just control or M0 images be used for this?
    asl_reference_wf = init_asl_reference_wf(omp_nthreads=omp_nthreads)
    asl_reference_wf.inputs.inputnode.dummy_scans = 0

    workflow.connect([(inputnode, asl_reference_wf, [("asl_file", "inputnode.asl_file")])])

    if sbref_file is not None:
        workflow.connect([(val_sbref, asl_reference_wf, [("out_file", "inputnode.sbref_file")])])

    # Split 4D ASL file into list of 3D volumes, so that volume-wise transforms (e.g., HMC params)
    # can be applied with other transforms in single shots.
    asl_split = pe.Node(FSLSplit(dimension="t"), name="asl_split", mem_gb=mem_gb["filesize"] * 3)

    # aslbuffer: an identity used as a pointer to either the original ASL file (if STD is disabled)
    # or the STC'ed ASL file for further use.
    aslbuffer = pe.Node(niu.IdentityInterface(fields=["asl_file"]), name="aslbuffer")

    workflow.connect([(aslbuffer, asl_split, [("asl_file", "in_file")])])

    # Slice timing correction
    if run_stc:
        asl_stc_wf = init_asl_stc_wf(name="asl_stc_wf", metadata=metadata)

        # fmt:off
        workflow.connect([
            (asl_reference_wf, asl_stc_wf, [("outputnode.skip_vols", "inputnode.skip_vols")]),
            (asl_stc_wf, aslbuffer, [("outputnode.stc_file", "asl_file")]),
        ])
        # fmt:on
        if not multiecho:
            # fmt:off
            workflow.connect([
                (asl_reference_wf, asl_stc_wf, [("outputnode.asl_file", "inputnode.asl_file")]),
            ])
            # fmt:on

        else:  # for meepi, iterate through stc_wf for all workflows
            meepi_echos = aslbuffer.clone(name="meepi_echos")
            meepi_echos.iterables = ("asl_file", asl_file)
            workflow.connect([(meepi_echos, asl_stc_wf, [("asl_file", "inputnode.asl_file")])])

    elif not multiecho:  # STC is too short or False
        # bypass STC from original asl to the splitter through aslbuffer
        workflow.connect([(asl_reference_wf, aslbuffer, [("outputnode.asl_file", "asl_file")])])

    else:
        # for meepi, iterate over all meepi echos to aslbuffer
        aslbuffer.iterables = ("asl_file", asl_file)

    # Head motion correction
    asl_hmc_wf = init_asl_hmc_wf(
        name="asl_hmc_wf",
        mem_gb=mem_gb["filesize"],
        omp_nthreads=omp_nthreads,
    )

    # fmt:off
    workflow.connect([
        (asl_reference_wf, asl_hmc_wf, [
            ("outputnode.raw_ref_image", "inputnode.raw_ref_image"),
            ("outputnode.asl_file", "inputnode.asl_file"),
        ]),
    ])
    # fmt:on

    # Susceptibility distortion correction
    asl_sdc_wf = init_sdc_estimate_wf(
        fmaps,
        metadata,
        omp_nthreads=omp_nthreads,
        debug=config.execution.debug,
    )

    # fmt:off
    workflow.connect([
        (t1w_brain, asl_sdc_wf, [("out_file", "inputnode.t1w_brain")]),
        (asl_reference_wf, asl_sdc_wf, [
            ("outputnode.ref_image", "inputnode.epi_file"),
            ("outputnode.ref_image_brain", "inputnode.epi_brain"),
            ("outputnode.asl_mask", "inputnode.epi_mask"),
        ]),
        (asl_sdc_wf, summary, [("outputnode.method", "distortion_correction")]),
    ])
    # fmt:on

    # EPI-T1 registration workflow
    asl_reg_wf = init_asl_reg_wf(
        asl2t1w_dof=config.workflow.asl2t1w_dof,
        asl2t1w_init=config.workflow.asl2t1w_init,
        name="asl_reg_wf",
        sloppy=config.execution.debug,
        use_bbr=config.workflow.use_bbr,
    )

    # fmt:off
    workflow.connect([
        (inputnode, asl_reg_wf, [("t1w_dseg", "inputnode.t1w_dseg")]),
        (t1w_brain, asl_reg_wf, [("out_file", "inputnode.t1w_brain")]),
        (asl_sdc_wf, asl_reg_wf, [("outputnode.epi_brain", "inputnode.ref_asl_brain")]),
        (asl_reg_wf, summary, [("outputnode.fallback", "fallback")]),
    ])
    # fmt:on

    # Apply ASL registration to T1w
    nonstd_spaces = set(spaces.get_nonstandard())

    asl_t1_trans_wf = init_asl_t1_trans_wf(
        name="asl_t1_trans_wf",
        use_fieldwarp=bool(fmaps),
        use_hmc=True,
        multiecho=multiecho,
        output_t1space=nonstd_spaces.intersection(("T1w", "anat")),
        scorescrub=scorescrub,
        basil=basil,
        mem_gb=mem_gb["resampled"],
        omp_nthreads=omp_nthreads,
        use_compression=False,
    )

    # fmt:off
    workflow.connect([
        (inputnode, asl_t1_trans_wf, [
            ("asl_file", "inputnode.name_source"),
            ("t1w_mask", "inputnode.t1w_mask"),
        ]),
        (t1w_brain, asl_t1_trans_wf, [("out_file", "inputnode.t1w_brain")]),
        (asl_sdc_wf, asl_t1_trans_wf, [
            ("outputnode.out_warp", "inputnode.fieldwarp"),
            ("outputnode.epi_mask", "inputnode.ref_asl_mask"),
            ("outputnode.epi_brain", "inputnode.ref_asl_brain"),
        ]),
        # unused if multiecho, but this is safe
        (asl_hmc_wf, asl_t1_trans_wf, [("outputnode.xforms", "inputnode.hmc_xforms")]),
        (asl_reg_wf, asl_t1_trans_wf, [
            ("outputnode.aslref_to_anat_xfm", "inputnode.aslref_to_anat_xfm"),
        ]),
        (asl_t1_trans_wf, asl_derivatives_wf, [
            ("outputnode.asl_t1", "inputnode.asl_t1"),
            ("outputnode.aslref_t1", "inputnode.aslref_t1"),
        ]),
    ])
    # fmt:on

    # Apply transforms in 1 shot
    asl_asl_trans_wf = init_asl_preproc_trans_wf(
        mem_gb=mem_gb["resampled"],
        omp_nthreads=omp_nthreads,
        use_compression=not config.execution.low_mem,
        use_fieldwarp=bool(fmaps),
        name="asl_asl_trans_wf",
    )
    asl_asl_trans_wf.inputs.inputnode.name_source = ref_file

    # fmt:off
    workflow.connect([
        (asl_split, asl_asl_trans_wf, [("out_files", "inputnode.asl_file")]),
        (asl_hmc_wf, asl_asl_trans_wf, [("outputnode.xforms", "inputnode.hmc_xforms")]),
        (asl_sdc_wf, asl_asl_trans_wf, [
            ("outputnode.out_warp", "inputnode.fieldwarp"),
            ("outputnode.epi_mask", "inputnode.asl_mask"),
        ]),
    ])
    # fmt:on

    # get confounds
    asl_confounds_wf = init_asl_confs_wf(mem_gb=mem_gb["largemem"], name="asl_confounds_wf")
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
            ("outputnode.anat_to_aslref_xfm", "inputnode.anat_to_aslref_xfm"),
        ]),
        (asl_reference_wf, asl_confounds_wf, [("outputnode.skip_vols", "inputnode.skip_vols")]),
        (asl_asl_trans_wf, asl_confounds_wf, [("outputnode.asl_mask", "inputnode.asl_mask")]),
        (asl_confounds_wf, asl_derivatives_wf, [
            ("outputnode.confounds_file", "inputnode.confounds"),
            ("outputnode.confounds_metadata", "inputnode.confounds_metadata"),
        ]),
        (asl_confounds_wf, summary, [("outputnode.confounds_file", "confounds_file")]),
    ])
    # fmt:on

    # Special steps for multi-echo data
    # XXX: This may not be reachable, since echo is not an allowed entity for ASL data.
    # XXX: Also, optimal combination via tedana may not be a good idea for multi-echo ASL data.
    if multiecho:
        from aslprep.niworkflows.func.util import init_skullstrip_asl_wf

        skullstrip_asl_wf = init_skullstrip_asl_wf(name="skullstrip_asl_wf")

        inputnode.inputs.asl_file = ref_file  # Replace reference w first echo

        join_echos = pe.JoinNode(
            niu.IdentityInterface(fields=["asl_files"]),
            joinsource=("meepi_echos" if run_stc is True else "aslbuffer"),
            joinfield=["asl_files"],
            name="join_echos",
        )

        # create optimal combination, adaptive T2* map
        asl_t2s_wf = init_asl_t2s_wf(
            echo_times=tes,
            name="asl_t2smap_wf",
        )

        # fmt:off
        workflow.connect([
            (skullstrip_asl_wf, join_echos, [("outputnode.skull_stripped_file", "asl_files")]),
            (join_echos, asl_t2s_wf, [("asl_files", "inputnode.asl_file")]),
            (inputnode, asl_derivatives_wf, [
                (("asl_file", combine_meepi_source), "inputnode.source_file"),
            ]),
            (asl_asl_trans_wf, skullstrip_asl_wf, [("outputnode.asl", "inputnode.in_file")]),
            (asl_t2s_wf, asl_confounds_wf, [("outputnode.asl", "inputnode.asl")]),
            (asl_t2s_wf, asl_t1_trans_wf, [("outputnode.asl", "inputnode.asl_split")]),
        ])
        # fmt:on

    else:
        # fmt:off
        workflow.connect([
            (inputnode, asl_derivatives_wf, [("asl_file", "inputnode.source_file")]),
            (asl_asl_trans_wf, asl_confounds_wf, [("outputnode.asl", "inputnode.asl")]),
            (asl_split, asl_t1_trans_wf, [("out_files", "inputnode.asl_split")]),
        ])
        # fmt:on

    # compute the CBF here
    compute_cbf_wf = init_cbf_compt_wf(
        name_source=asl_file,
        aslcontext=run_data["aslcontext"],
        dummy_vols=dummyvols,
        M0Scale=mscale,
        scorescrub=scorescrub,
        basil=basil,
        smooth_kernel=smoothkernel,
        metadata=metadata,
        name="compute_cbf_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, compute_cbf_wf, [
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("m0scan", "inputnode.m0scan"),
            ("m0scan_metadata", "inputnode.m0scan_metadata"),
            ("t1w_mask", "inputnode.t1w_mask"),
        ]),
        (asl_asl_trans_wf, compute_cbf_wf, [
            ("outputnode.asl", "inputnode.asl_file"),
            ("outputnode.asl_mask", "inputnode.asl_mask"),
        ]),
        (asl_reg_wf, compute_cbf_wf, [
            ("outputnode.anat_to_aslref_xfm", "inputnode.anat_to_aslref_xfm"),
            ("outputnode.aslref_to_anat_xfm", "inputnode.aslref_to_anat_xfm"),
        ]),
        (asl_reg_wf, asl_derivatives_wf, [
            ("outputnode.anat_to_aslref_xfm", "inputnode.anat_to_aslref_xfm"),
            ("outputnode.aslref_to_anat_xfm", "inputnode.aslref_to_anat_xfm"),
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
        (asl_asl_trans_wf, refine_mask, [("outputnode.asl_mask", "asl_mask")]),
    ])
    # fmt:on

    if fmaps:
        from aslprep.sdcflows.workflows.outputs import init_sdc_unwarp_report_wf

        # Report on asl correction
        fmap_unwarp_report_wf = init_sdc_unwarp_report_wf()
        # fmt:off
        workflow.connect([
            (inputnode, fmap_unwarp_report_wf, [("t1w_dseg", "inputnode.in_seg")]),
            (asl_reference_wf, fmap_unwarp_report_wf, [
                ("outputnode.ref_image", "inputnode.in_pre"),
            ]),
            (asl_reg_wf, fmap_unwarp_report_wf, [
                ("outputnode.anat_to_aslref_xfm", "inputnode.in_xfm"),
            ]),
            (asl_sdc_wf, fmap_unwarp_report_wf, [
                ("outputnode.epi_corrected", "inputnode.in_post"),
            ]),
        ])
        # fmt:on

        # Overwrite ``out_path_base`` of unwarping DataSinks
        # And ensure echo is dropped from report
        for node in fmap_unwarp_report_wf.list_node_names():
            if node.split(".")[-1].startswith("ds_"):
                fmap_unwarp_report_wf.get_node(node).interface.out_path_base = "aslprep"
                fmap_unwarp_report_wf.get_node(node).inputs.dismiss_entities = ("echo",)

        for node in asl_sdc_wf.list_node_names():
            if node.split(".")[-1].startswith("ds_"):
                asl_sdc_wf.get_node(node).interface.out_path_base = "aslprep"
                asl_sdc_wf.get_node(node).inputs.dismiss_entities = ("echo",)

        if "syn" in fmaps:
            sdc_select_std = pe.Node(
                KeySelect(fields=["template_to_anat_xfm"], key="MNI152NLin2009cAsym"),
                name="sdc_select_std",
                run_without_submitting=True,
            )

            # fmt:off
            workflow.connect([
                (inputnode, sdc_select_std, [
                    ("template_to_anat_xfm", "template_to_anat_xfm"),
                    ("template", "keys"),
                ]),
                (sdc_select_std, asl_sdc_wf, [
                    ("template_to_anat_xfm", "inputnode.template_to_anat_xfm"),
                ]),
            ])
            # fmt:on

        if fmaps.get("syn") is True:  # SyN forced
            syn_unwarp_report_wf = init_sdc_unwarp_report_wf(
                name="syn_unwarp_report_wf", forcedsyn=True
            )
            # fmt:off
            workflow.connect([
                (inputnode, syn_unwarp_report_wf, [("t1w_dseg", "inputnode.in_seg")]),
                (asl_reference_wf, syn_unwarp_report_wf, [
                    ("outputnode.ref_image", "inputnode.in_pre"),
                ]),
                (asl_reg_wf, syn_unwarp_report_wf, [
                    ("outputnode.anat_to_aslref_xfm", "inputnode.in_xfm"),
                ]),
                (asl_sdc_wf, syn_unwarp_report_wf, [
                    ("outputnode.syn_ref", "inputnode.in_post"),
                ]),
            ])
            # fmt:on

            # Overwrite ``out_path_base`` of unwarping DataSinks
            # And ensure echo is dropped from report
            for node in syn_unwarp_report_wf.list_node_names():
                if node.split(".")[-1].startswith("ds_"):
                    syn_unwarp_report_wf.get_node(node).interface.out_path_base = "aslprep"
                    syn_unwarp_report_wf.get_node(node).inputs.dismiss_entities = ("echo",)

    # Map final asl mask into T1w space (if required)
    nonstd_spaces = set(spaces.get_nonstandard())
    if nonstd_spaces.intersection(("T1w", "anat")):
        from aslprep.interfaces.ants import ApplyTransforms

        # fmt:off
        workflow.connect([
            (compute_cbf_wf, asl_t1_trans_wf, [
                ("outputnode.cbf_ts", "inputnode.cbf_ts"),
                ("outputnode.mean_cbf", "inputnode.mean_cbf"),
            ]),
            (asl_t1_trans_wf, asl_derivatives_wf, [
                ("outputnode.cbf_ts_t1", "inputnode.cbf_ts_t1"),
                ("outputnode.mean_cbf_t1", "inputnode.mean_cbf_t1"),
            ]),
        ])
        # fmt:on

        # NOTE: Can this be bundled into the ASL-T1w transform workflow?
        aslmask_to_t1w = pe.Node(
            ApplyTransforms(interpolation="MultiLabel"), name="aslmask_to_t1w", mem_gb=0.1
        )

        # fmt:off
        workflow.connect([
            (asl_reg_wf, aslmask_to_t1w, [("outputnode.aslref_to_anat_xfm", "transforms")]),
            (asl_t1_trans_wf, aslmask_to_t1w, [("outputnode.asl_mask_t1", "reference_image")]),
            (refine_mask, aslmask_to_t1w, [("out_mask", "input_image")]),
            (aslmask_to_t1w, asl_derivatives_wf, [("output_image", "inputnode.asl_mask_t1")]),
        ])
        # fmt:on

        if scorescrub:
            # fmt:off
            workflow.connect([
                (compute_cbf_wf, asl_t1_trans_wf, [
                    ("outputnode.cbf_ts_score", "inputnode.cbf_ts_score"),
                    ("outputnode.mean_cbf_score", "inputnode.mean_cbf_score"),
                    ("outputnode.mean_cbf_scrub", "inputnode.mean_cbf_scrub"),
                ]),
                (asl_t1_trans_wf, asl_derivatives_wf, [
                    ("outputnode.mean_cbf_scrub_t1", "inputnode.mean_cbf_scrub_t1"),
                    ("outputnode.cbf_ts_score_t1", "inputnode.cbf_ts_score_t1"),
                    ("outputnode.mean_cbf_score_t1", "inputnode.mean_cbf_score_t1"),
                ]),
            ])
            # fmt:on

        if basil:
            # fmt:off
            workflow.connect([
                (compute_cbf_wf, asl_t1_trans_wf, [
                    ("outputnode.mean_cbf_basil", "inputnode.mean_cbf_basil"),
                    ("outputnode.mean_cbf_gm_basil", "inputnode.mean_cbf_gm_basil"),
                    ("outputnode.mean_cbf_wm_basil", "inputnode.mean_cbf_wm_basil"),
                    ("outputnode.att", "inputnode.att"),
                ]),
                (asl_t1_trans_wf, asl_derivatives_wf, [
                    ("outputnode.mean_cbf_basil_t1", "inputnode.mean_cbf_basil_t1"),
                    ("outputnode.mean_cbf_gm_basil_t1", "inputnode.mean_cbf_gm_basil_t1"),
                    ("outputnode.mean_cbf_wm_basil_t1", "inputnode.mean_cbf_wm_basil_t1"),
                    ("outputnode.att_t1", "inputnode.att_t1"),
                ]),
            ])
            # fmt:on

    if nonstd_spaces.intersection(("func", "run", "asl", "aslref", "sbref")):
        # fmt:off
        workflow.connect([
            (asl_asl_trans_wf, asl_derivatives_wf, [
                ("outputnode.asl", "inputnode.asl_native"),
                ("outputnode.aslref", "inputnode.aslref_native"),
            ]),
            (refine_mask, asl_derivatives_wf, [
                ("out_mask", "inputnode.asl_mask_native"),
            ]),
            (compute_cbf_wf, asl_derivatives_wf, [
                ("outputnode.cbf_ts", "inputnode.cbf_ts_native"),
                ("outputnode.mean_cbf", "inputnode.mean_cbf_native"),
            ]),
        ])
        # fmt:on

        if scorescrub:
            # fmt:off
            workflow.connect([
                (compute_cbf_wf, asl_derivatives_wf, [
                    ("outputnode.cbf_ts_score", "inputnode.cbf_ts_score_native"),
                    ("outputnode.mean_cbf_score", "inputnode.mean_cbf_score_native"),
                    ("outputnode.mean_cbf_scrub", "inputnode.mean_cbf_scrub_native"),
                ]),
            ])
            # fmt:on

        if basil:
            # fmt:off
            workflow.connect([
                (compute_cbf_wf, asl_derivatives_wf, [
                    ("outputnode.mean_cbf_basil", "inputnode.mean_cbf_basil_native"),
                    ("outputnode.mean_cbf_gm_basil", "inputnode.mean_cbf_gm_basil_native"),
                    ("outputnode.mean_cbf_wm_basil", "inputnode.mean_cbf_wm_basil_native"),
                    ("outputnode.att", "inputnode.att_native"),
                ]),
            ])
            # fmt:on

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # Apply transforms in 1 shot
        # Only use uncompressed output if AROMA is to be run
        asl_std_trans_wf = init_asl_std_trans_wf(
            mem_gb=mem_gb["resampled"],
            omp_nthreads=omp_nthreads,
            spaces=spaces,
            scorescrub=scorescrub,
            basil=basil,
            name="asl_std_trans_wf",
            use_compression=not config.execution.low_mem,
            use_fieldwarp=bool(fmaps),
        )

        # fmt:off
        workflow.connect([
            (inputnode, asl_std_trans_wf, [
                ("template", "inputnode.templates"),
                ("anat_to_template_xfm", "inputnode.anat_to_template_xfm"),
                ("asl_file", "inputnode.name_source"),
            ]),
            (asl_hmc_wf, asl_std_trans_wf, [("outputnode.xforms", "inputnode.hmc_xforms")]),
            (asl_reg_wf, asl_std_trans_wf, [
                ("outputnode.aslref_to_anat_xfm", "inputnode.aslref_to_anat_xfm"),
            ]),
            (refine_mask, asl_std_trans_wf, [("out_mask", "inputnode.asl_mask")]),
            (asl_sdc_wf, asl_std_trans_wf, [("outputnode.out_warp", "inputnode.fieldwarp")]),
            (compute_cbf_wf, asl_std_trans_wf, [
                ("outputnode.cbf_ts", "inputnode.cbf_ts"),
                ("outputnode.mean_cbf", "inputnode.mean_cbf"),
            ]),
        ])
        # fmt:on

        if scorescrub:
            # fmt:off
            workflow.connect([
                (compute_cbf_wf, asl_std_trans_wf, [
                    ("outputnode.cbf_ts_score", "inputnode.cbf_ts_score"),
                    ("outputnode.mean_cbf_score", "inputnode.mean_cbf_score"),
                    ("outputnode.mean_cbf_scrub", "inputnode.mean_cbf_scrub"),
                ]),
            ])
            # fmt:on

        if basil:
            # fmt:off
            workflow.connect([
                (compute_cbf_wf, asl_std_trans_wf, [
                    ("outputnode.mean_cbf_basil", "inputnode.mean_cbf_basil"),
                    ("outputnode.mean_cbf_gm_basil", "inputnode.mean_cbf_gm_basil"),
                    ("outputnode.mean_cbf_wm_basil", "inputnode.mean_cbf_wm_basil"),
                    ("outputnode.att", "inputnode.att"),
                ]),
            ])
            # fmt:on

        if not multiecho:
            # fmt:off
            workflow.connect([
                (asl_split, asl_std_trans_wf, [("out_files", "inputnode.asl_split")]),
            ])
            # fmt:on

        else:
            split_opt_comb = asl_split.clone(name="split_opt_comb")

            # fmt:off
            workflow.connect([
                (asl_t2s_wf, split_opt_comb, [("outputnode.asl", "in_file")]),
                (split_opt_comb, asl_std_trans_wf, [("out_files", "inputnode.asl_split")]),
            ])
            # fmt:on

        # asl_derivatives_wf internally parametrizes over snapshotted spaces.
        # fmt:off
        workflow.connect([
            (asl_std_trans_wf, asl_derivatives_wf, [
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
                (asl_std_trans_wf, asl_derivatives_wf, [
                    ("outputnode.cbf_ts_score_std", "inputnode.cbf_ts_score_std"),
                    ("outputnode.mean_cbf_score_std", "inputnode.mean_cbf_score_std"),
                    ("outputnode.mean_cbf_scrub_std", "inputnode.mean_cbf_scrub_std"),
                ]),
            ])
            # fmt:on

        if basil:
            # fmt:off
            workflow.connect([
                (asl_std_trans_wf, asl_derivatives_wf, [
                    ("outputnode.mean_cbf_basil_std", "inputnode.mean_cbf_basil_std"),
                    ("outputnode.mean_cbf_gm_basil_std", "inputnode.mean_cbf_gm_basil_std"),
                    ("outputnode.mean_cbf_wm_basil_std", "inputnode.mean_cbf_wm_basil_std"),
                    ("outputnode.att_std", "inputnode.att_std"),
                ]),
            ])
            # fmt:on

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
            ("outputnode.anat_to_aslref_xfm", "inputnode.anat_to_aslref_xfm"),
        ]),
        (compute_cbf_wf, compute_cbf_qc_wf, [("outputnode.mean_cbf", "inputnode.mean_cbf")]),
        (asl_confounds_wf, compute_cbf_qc_wf, [
            ("outputnode.confounds_file", "inputnode.confounds_file"),
        ]),
        (compute_cbf_qc_wf, asl_derivatives_wf, [("outputnode.qc_file", "inputnode.qc_file")]),
        (compute_cbf_qc_wf, summary, [("outputnode.qc_file", "qc_file")]),
    ])
    # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (compute_cbf_wf, compute_cbf_qc_wf, [
                ("outputnode.mean_cbf_score", "inputnode.mean_cbf_score"),
                ("outputnode.mean_cbf_scrub", "inputnode.mean_cbf_scrub"),
            ]),
        ])
        # fmt:on

    if basil:
        # fmt:off
        workflow.connect([
            (compute_cbf_wf, compute_cbf_qc_wf, [
                ("outputnode.mean_cbf_basil", "inputnode.mean_cbf_basil"),
                ("outputnode.mean_cbf_gm_basil", "inputnode.mean_cbf_gm_basil"),
            ]),
        ])
        # fmt:on

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # fmt:off
        workflow.connect([
            (asl_std_trans_wf, compute_cbf_qc_wf, [
                ("outputnode.asl_mask_std", "inputnode.asl_mask_std"),
            ]),
        ])
        # fmt:on

    cbf_plot = init_cbfplot_wf(
        metadata=metadata,
        scorescrub=scorescrub,
        basil=basil,
        name="cbf_plot",
    )
    # fmt:off

    workflow.connect([
        (compute_cbf_wf, cbf_plot, [
            ("outputnode.mean_cbf", "inputnode.mean_cbf"),
            ("outputnode.mean_cbf_score", "inputnode.mean_cbf_score"),
            ("outputnode.mean_cbf_scrub", "inputnode.mean_cbf_scrub"),
            ("outputnode.mean_cbf_basil", "inputnode.mean_cbf_basil"),
            ("outputnode.mean_cbf_gm_basil", "inputnode.mean_cbf_gm_basil"),
            ("outputnode.cbf_ts_score", "inputnode.cbf_ts_score"),
            ("outputnode.cbf_ts", "inputnode.cbf_ts"),
            ("outputnode.score_outlier_index", "inputnode.score_outlier_index"),
        ]),
        (inputnode, cbf_plot, [("template_to_anat_xfm", "inputnode.template_to_anat_xfm")]),
        (asl_reg_wf, cbf_plot, [
            ("outputnode.anat_to_aslref_xfm", "inputnode.anat_to_aslref_xfm"),
        ]),
        (refine_mask, cbf_plot, [("out_mask", "inputnode.asl_mask")]),
        (asl_reference_wf, cbf_plot, [("outputnode.ref_image_brain", "inputnode.aslref")]),
        (asl_confounds_wf, cbf_plot, [("outputnode.confounds_file", "inputnode.confounds_file")]),
    ])
    # fmt:on

    carpetplot_wf = init_carpetplot_wf(
        mem_gb=mem_gb["resampled"],
        metadata=metadata,
        # cifti_output=config.workflow.cifti_output,
        name="carpetplot_wf",
    )

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
        (select_xform_MNI152NLin2009cAsym_to_t1w, carpetplot_wf, [
            ("template_to_anat_xfm", "inputnode.template_to_anat_xfm"),
        ]),
        (asl_asl_trans_wf if not multiecho else asl_t2s_wf, carpetplot_wf, [
            ("outputnode.asl", "inputnode.asl"),
        ]),
        (refine_mask, carpetplot_wf, [("out_mask", "inputnode.asl_mask")]),
        (asl_reg_wf, carpetplot_wf, [
            ("outputnode.anat_to_aslref_xfm", "inputnode.anat_to_aslref_xfm"),
        ]),
        (asl_confounds_wf, carpetplot_wf, [
            ("outputnode.confounds_file", "inputnode.confounds_file"),
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
        (select_xform_MNI152NLin2009cAsym_to_t1w, parcellate_cbf_wf, [
            ("template_to_anat_xfm", "inputnode.MNI152NLin2009cAsym_to_anat_xfm"),
        ]),
        (asl_asl_trans_wf, parcellate_cbf_wf, [("outputnode.asl_mask", "inputnode.asl_mask")]),
        (asl_reg_wf, parcellate_cbf_wf, [
            ("outputnode.anat_to_aslref_xfm", "inputnode.anat_to_aslref_xfm"),
        ]),
        (compute_cbf_wf, parcellate_cbf_wf, [("outputnode.mean_cbf", "inputnode.mean_cbf")]),
        (parcellate_cbf_wf, asl_derivatives_wf, [
            ("outputnode.atlas_names", "inputnode.atlas_names"),
            ("outputnode.mean_cbf_parcellated", "inputnode.mean_cbf_parcellated"),
        ]),
    ])
    # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (compute_cbf_wf, parcellate_cbf_wf, [
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
            (compute_cbf_wf, parcellate_cbf_wf, [
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

    ds_report_validation = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
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
        (asl_reference_wf, ds_report_validation, [("outputnode.validation_report", "in_file")]),
    ])
    # fmt:on

    # Fill-in datasinks of reportlets seen so far
    for node in workflow.list_node_names():
        if node.split(".")[-1].startswith("ds_report"):
            workflow.get_node(node).inputs.base_directory = output_dir
            workflow.get_node(node).inputs.source_file = ref_file

    return workflow
