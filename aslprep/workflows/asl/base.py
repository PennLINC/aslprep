# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Preprocessing workflows for ASL data."""
import os

import nibabel as nb
from nipype.interfaces import utility as niu
from nipype.interfaces.fsl import Split as FSLSplit
from nipype.pipeline import engine as pe
from niworkflows.interfaces.reportlets.registration import (
    SimpleBeforeAfterRPT as SimpleBeforeAfter,
)
from niworkflows.interfaces.utility import KeySelect
from niworkflows.utils.connections import listify, pop_file

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.interfaces.reports import FunctionalSummary
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.func.util import init_asl_reference_wf
from aslprep.niworkflows.interfaces.nibabel import ApplyMask
from aslprep.utils.meepi import combine_meepi_source
from aslprep.utils.misc import _create_mem_gb, _get_wf_name
from aslprep.workflows.asl.cbf import (
    init_cbf_compt_wf,
    init_cbfplot_wf,
    init_cbfqc_compt_wf,
    init_cbfroiquant_wf,
)
from aslprep.workflows.asl.confounds import init_asl_confs_wf, init_carpetplot_wf
from aslprep.workflows.asl.hmc import init_asl_hmc_wf
from aslprep.workflows.asl.outputs import init_asl_derivatives_wf
from aslprep.workflows.asl.registration import init_asl_reg_wf, init_asl_t1_trans_wf
from aslprep.workflows.asl.resampling import (
    init_asl_preproc_trans_wf,
    init_asl_std_trans_wf,
    init_asl_surf_wf,
)
from aslprep.workflows.asl.stc import init_asl_stc_wf
from aslprep.workflows.asl.t2s import init_asl_t2s_wf, init_t2s_reporting_wf


def init_asl_preproc_wf(asl_file, has_fieldmap=False):
    r"""Perform the functional preprocessing stages of ASLPrep.

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
        Path to NIfTI file (single echo) or list of paths to NIfTI files (multi-echo)
    has_fieldmap : :obj:`bool`
        Signals the workflow to use inputnode fieldmap files

    Inputs
    ------
    asl_file
        ASL series NIfTI file
    t1w_preproc
        Bias-corrected structural template image
    t1w_mask
        Mask of the skull-stripped template image
    t1w_dseg
        Segmentation of preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF)
    t1w_aseg
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
    asl_t1_ref
        ASL reference image, resampled to T1w space
    asl2anat_xfm
        Affine transform from ASL reference space to T1w space
    anat2asl_xfm
        Affine transform from T1w space to ASL reference space
    hmc_xforms
        Affine transforms for each ASL volume to the ASL reference
    asl_mask_t1
        ASL series mask in T1w space
    asl_aseg_t1
        FreeSurfer ``aseg`` resampled to match ``asl_t1``
    asl_aparc_t1
        FreeSurfer ``aparc+aseg`` resampled to match ``asl_t1``
    asl_std
        ASL series, resampled to template space
    asl_std_ref
        ASL reference image, resampled to template space
    asl_mask_std
        ASL series mask in template space
    asl_aseg_std
        FreeSurfer ``aseg`` resampled to match ``asl_std``
    asl_aparc_std
        FreeSurfer ``aparc+aseg`` resampled to match ``asl_std``
    asl_native
        ASL series, with distortion corrections applied (native space)
    asl_native_ref
        ASL reference image in native space
    asl_mask_native
        ASL series mask in native space
    asl_echos_native
        Per-echo ASL series, with distortion corrections applied
    asl_cifti
        ASL CIFTI image
    cifti_metadata
        Path of metadata files corresponding to ``asl_cifti``.
    surfaces
        ASL series, resampled to FreeSurfer surfaces
    t2star_asl
        Estimated T2\\* map in ASL native space
    t2star_t1
        Estimated T2\\* map in T1w space
    t2star_std
        Estimated T2\\* map in template space
    confounds
        TSV of confounds
    confounds_metadata
        Confounds metadata dictionary
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
    """
    img = nb.load(asl_file[0] if isinstance(asl_file, (list, tuple)) else asl_file)
    nvols = 1 if img.ndim < 4 else img.shape[3]
    if nvols <= 5 - config.execution.sloppy:
        config.loggers.workflow.warning(
            f"Too short ASL series (<= 5 timepoints). Skipping processing of <{asl_file}>."
        )
        return

    mem_gb = {"filesize": 1, "resampled": 1, "largemem": 1}
    asl_tlen = 10

    # Have some options handy
    layout = config.execution.layout
    omp_nthreads = config.nipype.omp_nthreads
    freesurfer = config.workflow.run_reconall
    spaces = config.workflow.spaces
    output_dir = str(config.execution.output_dir)
    freesurfer_spaces = spaces.get_fs_spaces()
    project_goodvoxels = config.workflow.project_goodvoxels
    dummy_scans = config.workflow.dummy_scans
    smoothkernel = config.workflow.smooth_kernel
    mscale = config.workflow.m0_scale
    scorescrub = config.workflow.scorescrub
    basil = config.workflow.basil

    if project_goodvoxels and freesurfer_spaces != ["fsaverage"]:
        config.loggers.workflow.critical(
            f"--project-goodvoxels only works with fsaverage (requested: {freesurfer_spaces})"
        )
        config.loggers.workflow.warn("Disabling --project-goodvoxels")
        project_goodvoxels = False

    # Extract BIDS entities and metadata from ASL file(s)
    entities = extract_entities(asl_file)
    layout = config.execution.layout

    # Extract metadata
    all_metadata = [layout.get_metadata(fname) for fname in listify(asl_file)]

    # Take first file as reference
    ref_file = pop_file(asl_file)
    metadata = all_metadata[0]
    # get original image orientation
    ref_orientation = get_img_orientation(ref_file)

    echo_idxs = listify(entities.get("echo", []))
    multiecho = len(echo_idxs) > 2
    if len(echo_idxs) == 1:
        config.loggers.workflow.warning(
            f"Running a single echo <{ref_file}> from a seemingly multi-echo dataset."
        )
        asl_file = ref_file  # Just in case - drop the list

    if len(echo_idxs) == 2:
        raise RuntimeError(
            "Multi-echo processing requires at least three different echos (found two)."
        )

    if multiecho:
        # Drop echo entity for future queries, have a boolean shorthand
        entities.pop("echo", None)
        # reorder echoes from shortest to largest
        tes, asl_file = zip(
            *sorted([(layout.get_metadata(bf)["EchoTime"], bf) for bf in asl_file])
        )
        ref_file = asl_file[0]  # Reset reference to be the shortest TE

    if os.path.isfile(ref_file):
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

    # Find associated sbref, if possible
    overrides = {
        "suffix": "sbref",
        "extension": [".nii", ".nii.gz"],
    }
    if config.execution.bids_filters:
        overrides.update(config.execution.bids_filters.get("sbref", {}))
    sb_ents = {**entities, **overrides}
    sbref_files = layout.get(return_type="file", **sb_ents)

    sbref_msg = f"No single-band-reference found for {os.path.basename(ref_file)}."
    if sbref_files and "sbref" in config.workflow.ignore:
        sbref_msg = "Single-band reference file(s) found and ignored."
        sbref_files = []
    elif sbref_files:
        sbref_msg = (
            "Using single-band reference file(s) "
            f"{','.join([os.path.basename(sbf) for sbf in sbref_files])}."
        )
    config.loggers.workflow.info(sbref_msg)

    if has_fieldmap:
        # First check if specified via B0FieldSource
        estimator_key = listify(metadata.get("B0FieldSource"))

        if not estimator_key:
            import re
            from pathlib import Path

            from sdcflows.fieldmaps import get_identifier

            # Fallback to IntendedFor
            intended_rel = re.sub(
                r"^sub-[a-zA-Z0-9]*/",
                "",
                str(Path(asl_file if not multiecho else asl_file[0]).relative_to(layout.root)),
            )
            estimator_key = get_identifier(intended_rel)

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

    # Check whether STC must/can be run
    run_stc = bool(metadata.get("SliceTiming")) and "slicetiming" not in config.workflow.ignore

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
Non-gridded (surface) resamplings were performed using `mri_vol2surf`
(FreeSurfer).
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "subjects_dir",
                "subject_id",
                "t1w_preproc",
                "t1w_mask",
                "t1w_dseg",
                "t1w_tpms",
                "t1w_aseg",
                "t1w_aparc",
                "anat2std_xfm",
                "std2anat_xfm",
                "template",
                "anat_ribbon",
                "t1w2fsnative_xfm",
                "fsnative2t1w_xfm",
                "fmap",
                "fmap_ref",
                "fmap_coeff",
                "fmap_mask",
                "fmap_id",
                "sdc_method",
            ]
        ),
        name="inputnode",
    )
    inputnode.inputs.asl_file = asl_file

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_t1",
                "asl_t1_ref",
                "asl2anat_xfm",
                "anat2asl_xfm",
                "hmc_xforms",
                "asl_mask_t1",
                "asl_aseg_t1",
                "asl_aparc_t1",
                "asl_std",
                "asl_std_ref",
                "asl_mask_std",
                "asl_aseg_std",
                "asl_aparc_std",
                "asl_native",
                "asl_native_ref",
                "asl_mask_native",
                "asl_echos_native",
                "asl_cifti",
                "cifti_metadata",
                "surfaces",
                "t2star_asl",
                "t2star_t1",
                "t2star_std",
                "confounds",
                "confounds_metadata",
                "cbf_t1",
                "cbf_std",
                "meancbf_t1",
                "meancbf_std",
                "score_t1",
                "score_std",
                "avgscore_t1",
                "avgscore_std",
                "scrub_t1",
                "scrub_std",
                "basil_t1",
                "basil_std",
                "pv_t1",
                "pv_std",
                "pv_native",
                "att",
                "att_t1",
                "att_std",
                "pvwm_t1",
                "pvwm_std",
                "qc_file",
            ]
        ),
        name="outputnode",
    )

    # Generate a brain-masked conversion of the t1w
    t1w_brain = pe.Node(ApplyMask(), name="t1w_brain")

    # Track echo index - this allows us to treat multi- and single-echo workflows
    # almost identically
    echo_index = pe.Node(niu.IdentityInterface(fields=["echoidx"]), name="echo_index")
    if multiecho:
        echo_index.iterables = [("echoidx", range(len(asl_file)))]
    else:
        echo_index.inputs.echoidx = 0

    # ASL source: track original ASL file(s)
    asl_source = pe.Node(niu.Select(inlist=asl_file), name="asl_source")

    # ASL buffer: an identity used as a pointer to either the original ASL
    # or the STC'ed one for further use.
    aslbuffer = pe.Node(niu.IdentityInterface(fields=["asl_file"]), name="aslbuffer")

    summary = pe.Node(
        FunctionalSummary(
            slice_timing=run_stc,
            registration=("FSL", "FreeSurfer")[freesurfer],
            registration_dof=config.workflow.asl2t1w_dof,
            registration_init=config.workflow.asl2t1w_init,
            pe_direction=metadata.get("PhaseEncodingDirection"),
            echo_idx=echo_idxs,
            tr=metadata["RepetitionTime"],
            orientation=ref_orientation,
        ),
        name="summary",
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )
    summary.inputs.dummy_scans = config.workflow.dummy_scans

    asl_derivatives_wf = init_asl_derivatives_wf(
        bids_root=layout.root,
        cifti_output=config.workflow.cifti_output,
        freesurfer=freesurfer,
        project_goodvoxels=project_goodvoxels,
        all_metadata=all_metadata,
        multiecho=multiecho,
        output_dir=output_dir,
        spaces=spaces,
        scorescrub=scorescrub,
        basil=basil,
    )
    asl_derivatives_wf.inputs.inputnode.all_source_files = asl_file
    asl_derivatives_wf.inputs.inputnode.cifti_density = config.workflow.cifti_output

    # fmt:off
    workflow.connect([
        (outputnode, asl_derivatives_wf, [
            ("asl_t1", "inputnode.asl_t1"),
            ("asl_t1_ref", "inputnode.asl_t1_ref"),
            ("asl2anat_xfm", "inputnode.asl2anat_xfm"),
            ("anat2asl_xfm", "inputnode.anat2asl_xfm"),
            ("hmc_xforms", "inputnode.hmc_xforms"),
            ("asl_aseg_t1", "inputnode.asl_aseg_t1"),
            ("asl_aparc_t1", "inputnode.asl_aparc_t1"),
            ("asl_mask_t1", "inputnode.asl_mask_t1"),
            ("asl_native", "inputnode.asl_native"),
            ("asl_native_ref", "inputnode.asl_native_ref"),
            ("asl_mask_native", "inputnode.asl_mask_native"),
            ("asl_echos_native", "inputnode.asl_echos_native"),
            ("confounds", "inputnode.confounds"),
            ("surfaces", "inputnode.surf_files"),
            ("asl_cifti", "inputnode.asl_cifti"),
            ("cifti_metadata", "inputnode.cifti_metadata"),
            ("t2star_asl", "inputnode.t2star_asl"),
            ("t2star_t1", "inputnode.t2star_t1"),
            ("t2star_std", "inputnode.t2star_std"),
            ("confounds_metadata", "inputnode.confounds_metadata"),
            ("acompcor_masks", "inputnode.acompcor_masks"),
            ("tcompcor_mask", "inputnode.tcompcor_mask"),
        ]),
    ])
    # fmt:on

    # Generate a tentative aslref
    initial_aslref_wf = init_asl_reference_wf(
        name="initial_aslref_wf",
        omp_nthreads=omp_nthreads,
        asl_file=asl_file,
        sbref_files=sbref_files,
        multiecho=multiecho,
    )
    initial_aslref_wf.inputs.inputnode.dummy_scans = config.workflow.dummy_scans

    # Select validated ASL files (orientations checked or corrected)
    select_asl = pe.Node(niu.Select(), name="select_asl")

    # Top-level ASL splitter
    asl_split = pe.Node(FSLSplit(dimension="t"), name="asl_split", mem_gb=mem_gb["filesize"] * 3)

    # HMC on the ASL
    asl_hmc_wf = init_asl_hmc_wf(
        name="asl_hmc_wf",
        mem_gb=mem_gb["filesize"],
        omp_nthreads=omp_nthreads,
    )

    # calculate ASL registration to T1w
    asl_reg_wf = init_asl_reg_wf(
        asl2t1w_dof=config.workflow.asl2t1w_dof,
        asl2t1w_init=config.workflow.asl2t1w_init,
        freesurfer=freesurfer,
        mem_gb=mem_gb["resampled"],
        name="asl_reg_wf",
        omp_nthreads=omp_nthreads,
        sloppy=config.execution.sloppy,
        use_bbr=config.workflow.use_bbr,
        use_compression=False,
    )

    # apply ASL registration to T1w
    asl_t1_trans_wf = init_asl_t1_trans_wf(
        name="asl_t1_trans_wf",
        scorescrub=scorescrub,
        basil=basil,
        freesurfer=freesurfer,
        mem_gb=mem_gb["resampled"],
        omp_nthreads=omp_nthreads,
        use_compression=False,
    )
    asl_t1_trans_wf.inputs.inputnode.fieldwarp = "identity"

    # get confounds
    asl_confounds_wf = init_asl_confs_wf(
        mem_gb=mem_gb["largemem"],
        metadata=metadata,
        freesurfer=freesurfer,
        regressors_all_comps=config.workflow.regressors_all_comps,
        regressors_fd_th=config.workflow.regressors_fd_th,
        regressors_dvars_th=config.workflow.regressors_dvars_th,
        name="asl_confounds_wf",
    )
    asl_confounds_wf.get_node("inputnode").inputs.t1_transform_flags = [False]

    # SLICE-TIME CORRECTION (or bypass) #############################################
    if run_stc:
        asl_stc_wf = init_asl_stc_wf(name="asl_stc_wf", metadata=metadata)
        # fmt:off
        workflow.connect([
            (initial_aslref_wf, asl_stc_wf, [("outputnode.skip_vols", "inputnode.skip_vols")]),
            (select_asl, asl_stc_wf, [("out", "inputnode.asl_file")]),
            (asl_stc_wf, aslbuffer, [("outputnode.stc_file", "asl_file")]),
        ])
        # fmt:on

    # bypass STC from original ASL in both SE and ME cases
    else:
        workflow.connect([(select_asl, aslbuffer, [("out", "asl_file")])])

    # MULTI-ECHO EPI DATA #############################################
    if multiecho:  # instantiate relevant interfaces, imports
        split_opt_comb = asl_split.clone(name="split_opt_comb")

        inputnode.inputs.asl_file = ref_file  # Replace reference w first echo

        join_echos = pe.JoinNode(
            niu.IdentityInterface(fields=["asl_files"]),
            joinsource="echo_index",
            joinfield=["asl_files"],
            name="join_echos",
        )

        # create optimal combination, adaptive T2* map
        asl_t2s_wf = init_asl_t2s_wf(
            echo_times=tes,
            mem_gb=mem_gb["filesize"],
            omp_nthreads=omp_nthreads,
            name="asl_t2smap_wf",
        )

        t2s_reporting_wf = init_t2s_reporting_wf()

        ds_report_t2scomp = pe.Node(
            DerivativesDataSink(
                desc="t2scomp",
                datatype="figures",
                dismiss_entities=("echo",),
            ),
            name="ds_report_t2scomp",
            run_without_submitting=True,
        )

        ds_report_t2star_hist = pe.Node(
            DerivativesDataSink(
                desc="t2starhist",
                datatype="figures",
                dismiss_entities=("echo",),
            ),
            name="ds_report_t2star_hist",
            run_without_submitting=True,
        )

    asl_final = pe.Node(
        niu.IdentityInterface(fields=["asl", "aslref", "mask", "asl_echos", "t2star"]),
        name="asl_final",
    )

    # Generate a final ASL reference
    # This ASL references *does not use* single-band reference images.
    final_aslref_wf = init_asl_reference_wf(
        name="final_aslref_wf",
        omp_nthreads=omp_nthreads,
        multiecho=multiecho,
    )
    final_aslref_wf.__desc__ = None  # Unset description to avoid second appearance

    # MAIN WORKFLOW STRUCTURE #######################################################
    # fmt:off
    workflow.connect([
        # Prepare masked T1w image
        (inputnode, t1w_brain, [
            ("t1w_preproc", "in_file"),
            ("t1w_mask", "in_mask"),
        ]),
        # Select validated asl files per-echo
        (initial_aslref_wf, select_asl, [("outputnode.all_asl_files", "inlist")]),
        # ASL buffer has slice-time corrected if it was run, original otherwise
        (aslbuffer, asl_split, [("asl_file", "in_file")]),
        # HMC
        (initial_aslref_wf, asl_hmc_wf, [
            ("outputnode.raw_ref_image", "inputnode.raw_ref_image"),
            ("outputnode.asl_file", "inputnode.asl_file"),
        ]),
        (asl_hmc_wf, outputnode, [
            ("outputnode.xforms", "hmc_xforms"),
        ]),
        # EPI-T1w registration workflow
        (inputnode, asl_reg_wf, [
            ("t1w_dseg", "inputnode.t1w_dseg"),
            # Undefined if --fs-no-reconall, but this is safe
            ("subjects_dir", "inputnode.subjects_dir"),
            ("subject_id", "inputnode.subject_id"),
            ("fsnative2t1w_xfm", "inputnode.fsnative2t1w_xfm"),
        ]),
        (asl_final, asl_reg_wf, [
            ("aslref", "inputnode.ref_asl_brain")]),
        (t1w_brain, asl_reg_wf, [("out_file", "inputnode.t1w_brain")]),
        (inputnode, asl_t1_trans_wf, [
            ("asl_file", "inputnode.name_source"),
            ("t1w_mask", "inputnode.t1w_mask"),
            ("t1w_aseg", "inputnode.t1w_aseg"),
            ("t1w_aparc", "inputnode.t1w_aparc"),
        ]),
        (t1w_brain, asl_t1_trans_wf, [("out_file", "inputnode.t1w_brain")]),
        (asl_reg_wf, outputnode, [
            ("outputnode.itk_asl_to_t1", "asl2anat_xfm"),
            ("outputnode.itk_t1_to_asl", "anat2asl_xfm"),
        ]),
        (asl_reg_wf, asl_t1_trans_wf, [
            ("outputnode.itk_asl_to_t1", "inputnode.itk_asl_to_t1"),
        ]),
        (asl_final, asl_t1_trans_wf, [
            ("mask", "inputnode.ref_asl_mask"),
            ("aslref", "inputnode.ref_asl_brain"),
        ]),
        (asl_t1_trans_wf, outputnode, [
            ("outputnode.asl_t1", "asl_t1"),
            ("outputnode.asl_t1_ref", "asl_t1_ref"),
            ("outputnode.asl_aseg_t1", "asl_aseg_t1"),
            ("outputnode.asl_aparc_t1", "asl_aparc_t1"),
        ]),
        # Connect asl_confounds_wf
        (inputnode, asl_confounds_wf, [
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("t1w_mask", "inputnode.t1w_mask"),
        ]),
        (asl_hmc_wf, asl_confounds_wf, [
            ("outputnode.movpar_file", "inputnode.movpar_file"),
            ("outputnode.rmsd_file", "inputnode.rmsd_file"),
        ]),
        (asl_reg_wf, asl_confounds_wf, [
            ("outputnode.itk_t1_to_asl", "inputnode.t1_asl_xform")
        ]),
        (initial_aslref_wf, asl_confounds_wf, [
            ("outputnode.skip_vols", "inputnode.skip_vols"),
        ]),
        (initial_aslref_wf, final_aslref_wf, [
            ("outputnode.skip_vols", "inputnode.dummy_scans"),
        ]),
        (final_aslref_wf, asl_final, [
            ("outputnode.ref_image", "aslref"),
            ("outputnode.asl_mask", "mask"),
        ]),
        (asl_final, asl_confounds_wf, [
            ("asl", "inputnode.asl"),
            ("mask", "inputnode.asl_mask"),
        ]),
        (asl_confounds_wf, outputnode, [
            ("outputnode.confounds_file", "confounds"),
            ("outputnode.confounds_metadata", "confounds_metadata"),
            ("outputnode.acompcor_masks", "acompcor_masks"),
            ("outputnode.tcompcor_mask", "tcompcor_mask"),
        ]),
        # Native-space ASL files (if calculated)
        (asl_final, outputnode, [
            ("asl", "asl_native"),
            ("aslref", "asl_native_ref"),
            ("mask", "asl_mask_native"),
            ("asl_echos", "asl_echos_native"),
            ("t2star", "t2star_asl"),
        ]),
        # Summary
        (initial_aslref_wf, summary, [("outputnode.algo_dummy_scans", "algo_dummy_scans")]),
        (asl_reg_wf, summary, [("outputnode.fallback", "fallback")]),
        (outputnode, summary, [("confounds", "confounds_file")]),
        # Select echo indices for original/validated ASL files
        (echo_index, asl_source, [("echoidx", "index")]),
        (echo_index, select_asl, [("echoidx", "index")]),
    ])
    # fmt:on

    # for standard EPI data, pass along correct file
    if not multiecho:
        # fmt:off
        workflow.connect([
            (inputnode, asl_derivatives_wf, [("asl_file", "inputnode.source_file")]),
            (asl_split, asl_t1_trans_wf, [("out_files", "inputnode.asl_split")]),
            (asl_hmc_wf, asl_t1_trans_wf, [("outputnode.xforms", "inputnode.hmc_xforms")]),
        ])
        # fmt:on
    else:  # for meepi, use optimal combination
        # fmt:off
        workflow.connect([
            # update name source for optimal combination
            (inputnode, asl_derivatives_wf, [
                (("asl_file", combine_meepi_source), "inputnode.source_file"),
            ]),
            (join_echos, asl_t2s_wf, [("asl_files", "inputnode.asl_file")]),
            (join_echos, asl_final, [("asl_files", "asl_echos")]),
            (asl_t2s_wf, split_opt_comb, [("outputnode.asl", "in_file")]),
            (split_opt_comb, asl_t1_trans_wf, [("out_files", "inputnode.asl_split")]),
            (asl_t2s_wf, asl_final, [
                ("outputnode.asl", "asl"),
                ("outputnode.t2star_map", "t2star"),
            ]),
            (inputnode, t2s_reporting_wf, [("t1w_dseg", "inputnode.label_file")]),
            (asl_reg_wf, t2s_reporting_wf, [
                ("outputnode.itk_t1_to_asl", "inputnode.label_asl_xform"),
            ]),
            (asl_final, t2s_reporting_wf, [
                ("t2star", "inputnode.t2star_file"),
                ("aslref", "inputnode.aslref"),
            ]),
            (t2s_reporting_wf, ds_report_t2scomp, [("outputnode.t2s_comp_report", "in_file")]),
            (t2s_reporting_wf, ds_report_t2star_hist, [("outputnode.t2star_hist", "in_file")]),
        ])
        # fmt:on

        # Already applied in asl_asl_trans_wf, which inputs to asl_t2s_wf
        asl_t1_trans_wf.inputs.inputnode.hmc_xforms = "identity"

    # compute  the CBF here
    compt_cbf_wf = init_cbf_compt_wf(
        name="compt_cbf_wf",
        dummy_vols=dummy_scans,
        M0Scale=mscale,
        bids_dir=config.execution.bids_dir / f"sub-{config.execution.participant_label[0]}",
        scorescrub=scorescrub,
        basil=basil,
        smooth_kernel=smoothkernel,
        metadata=metadata,
    )

    # cbf computation workflow
    # fmt:off
    workflow.connect([
        (asl_final, compt_cbf_wf, [
            ("asl", "inputnode.asl_file"),
            ("mask", "inputnode.asl_mask"),
        ]),
        (inputnode, compt_cbf_wf, [
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("asl_file", "inputnode.in_file"),
            ("t1w_mask", "inputnode.t1w_mask")
        ]),
        (asl_reg_wf, compt_cbf_wf, [
            ("outputnode.itk_t1_to_asl", "inputnode.t1_asl_xform"),
            ("outputnode.itk_asl_to_t1", "inputnode.itk_asl_to_t1"),
        ]),
        (asl_reg_wf, outputnode, [
            ("outputnode.itk_t1_to_asl", "itk_t1_to_asl"),
            ("outputnode.itk_asl_to_t1", "itk_asl_to_t1"),
        ]),
        (asl_reg_wf, asl_derivatives_wf, [
            ("outputnode.itk_t1_to_asl", "inputnode.itk_t1_to_asl"),
            ("outputnode.itk_asl_to_t1", "inputnode.itk_asl_to_t1"),
        ]),
    ])
    # fmt:on

    # Map final ASL mask into T1w space (if required)
    nonstd_spaces = set(spaces.get_nonstandard())
    if nonstd_spaces.intersection(("T1w", "anat")):
        from niworkflows.interfaces.fixes import (
            FixHeaderApplyTransforms as ApplyTransforms,
        )

        aslmask_to_t1w = pe.Node(
            ApplyTransforms(interpolation="MultiLabel"),
            name="aslmask_to_t1w",
            mem_gb=0.1,
        )
        # fmt:off
        workflow.connect([
            (asl_reg_wf, aslmask_to_t1w, [("outputnode.itk_asl_to_t1", "transforms")]),
            (asl_t1_trans_wf, aslmask_to_t1w, [("outputnode.asl_mask_t1", "reference_image")]),
            (asl_final, aslmask_to_t1w, [("mask", "input_image")]),
            (aslmask_to_t1w, outputnode, [("output_image", "asl_mask_t1")]),
        ])
        # fmt:on

        if multiecho:
            t2star_to_t1w = pe.Node(
                ApplyTransforms(interpolation="LanczosWindowedSinc", float=True),
                name="t2star_to_t1w",
                mem_gb=0.1,
            )
            # fmt:off
            workflow.connect([
                (asl_reg_wf, t2star_to_t1w, [("outputnode.itk_asl_to_t1", "transforms")]),
                (asl_t1_trans_wf, t2star_to_t1w, [
                    ("outputnode.asl_mask_t1", "reference_image")
                ]),
                (asl_final, t2star_to_t1w, [("t2star", "input_image")]),
                (t2star_to_t1w, outputnode, [("output_image", "t2star_t1")]),
            ])
            # fmt:on

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # Apply transforms in 1 shot
        # Only use uncompressed output if AROMA is to be run
        asl_std_trans_wf = init_asl_std_trans_wf(
            freesurfer=freesurfer,
            mem_gb=mem_gb["resampled"],
            omp_nthreads=omp_nthreads,
            spaces=spaces,
            multiecho=multiecho,
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
            (asl_final, asl_std_trans_wf, [
                ("mask", "inputnode.asl_mask"),
                ("t2star", "inputnode.t2star"),
            ]),
            (asl_reg_wf, asl_std_trans_wf, [
                ("outputnode.itk_asl_to_t1", "inputnode.itk_asl_to_t1"),
            ]),
            (asl_std_trans_wf, outputnode, [
                ("outputnode.asl_std", "asl_std"),
                ("outputnode.asl_std_ref", "asl_std_ref"),
                ("outputnode.asl_mask_std", "asl_mask_std"),
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
                (asl_std_trans_wf, outputnode, [
                    ("outputnode.asl_aseg_std", "asl_aseg_std"),
                    ("outputnode.asl_aparc_std", "asl_aparc_std"),
                ]),
            ])
            # fmt:on

        if not multiecho:
            # fmt:off
            workflow.connect([
                (asl_split, asl_std_trans_wf, [("out_files", "inputnode.asl_split")]),
                (asl_hmc_wf, asl_std_trans_wf, [("outputnode.xforms", "inputnode.hmc_xforms")]),
            ])
            # fmt:on
        else:
            # fmt:off
            workflow.connect([
                (split_opt_comb, asl_std_trans_wf, [("out_files", "inputnode.asl_split")]),
                (asl_std_trans_wf, outputnode, [("outputnode.t2star_std", "t2star_std")]),
            ])
            # fmt:on

            # Already applied in asl_asl_trans_wf, which inputs to asl_t2s_wf
            asl_std_trans_wf.inputs.inputnode.hmc_xforms = "identity"

        # fmt:off
        # asl_derivatives_wf internally parametrizes over snapshotted spaces.
        workflow.connect([
            (asl_std_trans_wf, asl_derivatives_wf, [
                ("outputnode.template", "inputnode.template"),
                ("outputnode.spatial_reference", "inputnode.spatial_reference"),
                ("outputnode.asl_std_ref", "inputnode.asl_std_ref"),
                ("outputnode.asl_std", "inputnode.asl_std"),
                ("outputnode.asl_mask_std", "inputnode.asl_mask_std"),
            ]),
        ])
        # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (asl_std_trans_wf, asl_derivatives_wf, [
                ("outputnode.score_std", "inputnode.score_std"),
                ("outputnode.avgscore_std", "inputnode.avgscore_std"),
                ("outputnode.scrub_std", "inputnode.scrub_std"),
            ]),
        ])
        # fmt:on

        if basil:
            # fmt:off
            workflow.connect([
                (asl_std_trans_wf, asl_derivatives_wf, [
                    ("outputnode.basil_std", "inputnode.basil_std"),
                    ("outputnode.pv_std", "inputnode.pv_std"),
                    ("outputnode.pvwm_std", "inputnode.pvwm_std"),
                    ("outputnode.att_std", "inputnode.att_std"),
                ]),
            ])
            # fmt:on

    compt_qccbf_wf = init_cbfqc_compt_wf(
        name="compt_qccbf_wf",
        asl_file=asl_file,
        scorescrub=scorescrub,
        basil=basil,
    )

    # fmt:off
    workflow.connect([
        (asl_final, compt_qccbf_wf, [("mask", "inputnode.asl_mask")]),
        (inputnode, compt_qccbf_wf, [
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("t1w_mask", "inputnode.t1w_mask"),
        ]),
        (asl_reg_wf, compt_qccbf_wf, [("outputnode.itk_t1_to_asl", "inputnode.t1_asl_xform")]),
        (compt_cbf_wf, compt_qccbf_wf, [("outputnode.out_mean", "inputnode.meancbf")]),
        (asl_confounds_wf, compt_qccbf_wf, [("outputnode.confounds_file", "inputnode.confmat")]),
        (compt_qccbf_wf, outputnode, [("outputnode.qc_file", "qc_file")]),
        (compt_qccbf_wf, asl_derivatives_wf, [("outputnode.qc_file", "inputnode.qc_file")]),
        (compt_qccbf_wf, summary, [("outputnode.qc_file", "qc_file")]),
        (asl_hmc_wf, compt_qccbf_wf, [("outputnode.rmsd_file", "inputnode.rmsd_file")]),
    ])
    # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (compt_cbf_wf, compt_qccbf_wf, [
                ("outputnode.out_avgscore", "inputnode.avgscore"),
                ("outputnode.out_scrub", "inputnode.scrub"),
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
        (compt_cbf_wf, cbf_plot, [
            ("outputnode.out_mean", "inputnode.cbf"),
            ("outputnode.out_avgscore", "inputnode.score"),
            ("outputnode.out_scrub", "inputnode.scrub"),
            ("outputnode.out_cbfb", "inputnode.basil"),
            ("outputnode.out_cbfpv", "inputnode.pvc"),
            ("outputnode.out_score", "inputnode.score_ts"),
            ("outputnode.out_cbf", "inputnode.cbf_ts"),
            ("outputnode.out_scoreindex", "inputnode.scoreindex"),
        ]),
        (inputnode, cbf_plot, [("std2anat_xfm", "inputnode.std2anat_xfm")]),
        (asl_reg_wf, cbf_plot, [("outputnode.itk_t1_to_asl", "inputnode.t1_asl_xform")]),
        (asl_final, cbf_plot, [("mask", "inputnode.asl_mask")]),
        (final_aslref_wf, cbf_plot, [("outputnode.ref_image_brain", "inputnode.asl_ref")]),
        (asl_confounds_wf, cbf_plot, [("outputnode.confounds_file", "inputnode.confounds_file")]),
    ])
    # fmt:on

    # SURFACES ##################################################################################
    # Freesurfer
    if freesurfer and freesurfer_spaces:
        config.loggers.workflow.debug("Creating ASL surface-sampling workflow.")
        asl_surf_wf = init_asl_surf_wf(
            mem_gb=mem_gb["resampled"],
            surface_spaces=freesurfer_spaces,
            medial_surface_nan=config.workflow.medial_surface_nan,
            project_goodvoxels=project_goodvoxels,
            name="asl_surf_wf",
        )
        # fmt:off
        workflow.connect([
            (inputnode, asl_surf_wf, [
                ("subjects_dir", "inputnode.subjects_dir"),
                ("subject_id", "inputnode.subject_id"),
                ("t1w2fsnative_xfm", "inputnode.t1w2fsnative_xfm"),
                ("anat_ribbon", "inputnode.anat_ribbon"),
                ("t1w_mask", "inputnode.t1w_mask"),
            ]),
            (asl_t1_trans_wf, asl_surf_wf, [("outputnode.asl_t1", "inputnode.source_file")]),
            (asl_surf_wf, outputnode, [("outputnode.surfaces", "surfaces")]),
            (asl_surf_wf, asl_derivatives_wf, [
                ("outputnode.target", "inputnode.surf_refs"),
                ("outputnode.goodvoxels_ribbon", "inputnode.goodvoxels_ribbon"),
            ]),
        ])
        # fmt:on

        # CIFTI output
        if config.workflow.cifti_output:
            from aslprep.workflows.asl.resampling import init_asl_grayords_wf

            asl_grayords_wf = init_asl_grayords_wf(
                grayord_density=config.workflow.cifti_output,
                mem_gb=mem_gb["resampled"],
                repetition_time=metadata["RepetitionTime"],
            )

            # fmt:off
            workflow.connect([
                (asl_std_trans_wf, asl_grayords_wf, [
                    ("outputnode.asl_std", "inputnode.asl_std"),
                    ("outputnode.spatial_reference", "inputnode.spatial_reference"),
                ]),
                (asl_surf_wf, asl_grayords_wf, [
                    ("outputnode.surfaces", "inputnode.surf_files"),
                    ("outputnode.target", "inputnode.surf_refs"),
                ]),
                (asl_grayords_wf, outputnode, [
                    ("outputnode.cifti_asl", "asl_cifti"),
                    ("outputnode.cifti_metadata", "cifti_metadata"),
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
            workflow.connect(
                asl_grayords_wf, "outputnode.cifti_asl", carpetplot_wf, "inputnode.cifti_asl",
            )
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
            (carpetplot_select_std, carpetplot_wf, [
                ("std2anat_xfm", "inputnode.std2anat_xfm"),
            ]),
            (asl_final, carpetplot_wf, [
                ("asl", "inputnode.asl"),
                ("mask", "inputnode.asl_mask"),
            ]),
            (asl_reg_wf, carpetplot_wf, [
                ("outputnode.itk_t1_to_asl", "inputnode.t1_asl_xform"),
            ]),
            (asl_confounds_wf, carpetplot_wf, [
                ("outputnode.confounds_file", "inputnode.confounds_file"),
                ("outputnode.crown_mask", "inputnode.crown_mask"),
                (("outputnode.acompcor_masks", _last), "inputnode.acompcor_mask"),
            ]),
        ])
        # fmt:on

    cbfroiqu = init_cbfroiquant_wf(
        basil=basil,
        scorescrub=scorescrub,
        name="cbf_roiquant",
    )

    # fmt:off
    workflow.connect([
        (asl_final, cbfroiqu, [("mask", "inputnode.aslmask")]),
        (inputnode, cbfroiqu, [("std2anat_xfm", "inputnode.std2anat_xfm")]),
        (asl_reg_wf, cbfroiqu, [("outputnode.itk_t1_to_asl", "inputnode.t1_asl_xform")]),
        (compt_cbf_wf, cbfroiqu, [("outputnode.out_mean", "inputnode.cbf")]),
        (cbfroiqu, asl_derivatives_wf, [
            ("outputnode.cbf_hvoxf", "inputnode.cbf_hvoxf"),
            ("outputnode.cbf_sc207", "inputnode.cbf_sc207"),
            ("outputnode.cbf_sc217", "inputnode.cbf_sc217"),
            ("outputnode.cbf_sc407", "inputnode.cbf_sc407"),
            ("outputnode.cbf_sc417", "inputnode.cbf_sc417"),
        ]),
    ])
    # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (compt_cbf_wf, cbfroiqu, [
                ("outputnode.out_avgscore", "inputnode.score"),
                ("outputnode.out_scrub", "inputnode.scrub"),
            ]),
            (cbfroiqu, asl_derivatives_wf, [
                ("outputnode.score_hvoxf", "inputnode.score_hvoxf"),
                ("outputnode.score_sc207", "inputnode.score_sc207"),
                ("outputnode.score_sc217", "inputnode.score_sc217"),
                ("outputnode.score_sc407", "inputnode.score_sc407"),
                ("outputnode.score_sc417", "inputnode.score_sc417"),
                ("outputnode.scrub_hvoxf", "inputnode.scrub_hvoxf"),
                ("outputnode.scrub_sc207", "inputnode.scrub_sc207"),
                ("outputnode.scrub_sc217", "inputnode.scrub_sc217"),
                ("outputnode.scrub_sc407", "inputnode.scrub_sc407"),
                ("outputnode.scrub_sc417", "inputnode.scrub_sc417"),
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
        DerivativesDataSink(desc="validation", datatype="figures", dismiss_entities=("echo",)),
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
            workflow.get_node(node).inputs.base_directory = output_dir
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
            (asl_source, asl_asl_trans_wf, [("out", "inputnode.name_source")]),
            (asl_split, asl_asl_trans_wf, [("out_files", "inputnode.asl_file")]),
            (asl_hmc_wf, asl_asl_trans_wf, [
                ("outputnode.xforms", "inputnode.hmc_xforms"),
            ]),
        ])

        workflow.connect([
            (asl_asl_trans_wf, asl_final, [("outputnode.asl", "asl")]),
            (asl_asl_trans_wf, final_aslref_wf, [
                ("outputnode.asl", "inputnode.asl_file"),
            ]),
        ] if not multiecho else [
            (initial_aslref_wf, asl_t2s_wf, [
                ("outputnode.asl_mask", "inputnode.asl_mask"),
            ]),
            (asl_asl_trans_wf, join_echos, [
                ("outputnode.asl", "asl_files"),
            ]),
            (join_echos, final_aslref_wf, [
                ("asl_files", "inputnode.asl_file"),
            ]),
        ])
        # fmt:on
        return workflow

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
            base_directory=output_dir,
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
        (initial_aslref_wf, unwarp_wf, [
            ("outputnode.ref_image", "inputnode.distorted_ref"),
        ]),
        (coeff2epi_wf, unwarp_wf, [
            ("outputnode.fmap_coeff", "inputnode.fmap_coeff")]),
        (asl_hmc_wf, unwarp_wf, [
            ("outputnode.xforms", "inputnode.hmc_xforms")]),
        (initial_aslref_wf, sdc_report, [
            ("outputnode.ref_image", "before")]),
        (asl_split, unwarp_wf, [
            ("out_files", "inputnode.distorted")]),
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
                base_directory=output_dir,
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
                base_directory=output_dir,
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
            (initial_aslref_wf, sdc_coreg_report, [
                ("outputnode.ref_image", "before"),
            ]),
            (coeff2epi_wf, sdc_coreg_report, [
                ("coregister.inverse_warped_image", "after"),
            ]),
            (final_aslref_wf, sdc_coreg_report, [
                ("outputnode.asl_mask", "wm_seg"),
            ]),
            (inputnode, ds_report_sdc_coreg, [("asl_file", "source_file")]),
            (sdc_coreg_report, ds_report_sdc_coreg, [("out_report", "in_file")]),
            (unwarp_wf, fmap_report, [(("outputnode.fieldmap", pop_file), "fieldmap")]),
            (coeff2epi_wf, fmap_report, [
                ("coregister.inverse_warped_image", "reference"),
            ]),
            (final_aslref_wf, fmap_report, [
                ("outputnode.asl_mask", "mask"),
            ]),
            (fmap_report, ds_fmap_report, [("out_report", "in_file")]),
            (inputnode, ds_fmap_report, [("asl_file", "source_file")]),
        ])
        # fmt:on

    if not multiecho:
        # fmt:off
        workflow.connect([
            (unwarp_wf, asl_final, [("outputnode.corrected", "asl")]),
            # remaining workflow connections
            (unwarp_wf, final_aslref_wf, [
                ("outputnode.corrected", "inputnode.asl_file"),
            ]),
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

    # Finalize connections if ME-EPI
    join_sdc_echos = pe.JoinNode(
        niu.IdentityInterface(
            fields=[
                "fieldmap",
                "fieldwarp",
                "corrected",
                "corrected_ref",
                "corrected_mask",
            ]
        ),
        joinsource="echo_index",
        joinfield=[
            "fieldmap",
            "fieldwarp",
            "corrected",
            "corrected_ref",
            "corrected_mask",
        ],
        name="join_sdc_echos",
    )

    def _dpop(list_of_lists):
        return list_of_lists[0][0]

    # fmt:off
    workflow.connect([
        (unwarp_wf, join_echos, [
            ("outputnode.corrected", "asl_files"),
        ]),
        (unwarp_wf, join_sdc_echos, [
            ("outputnode.fieldmap", "fieldmap"),
            ("outputnode.fieldwarp", "fieldwarp"),
            ("outputnode.corrected", "corrected"),
            ("outputnode.corrected_ref", "corrected_ref"),
            ("outputnode.corrected_mask", "corrected_mask"),
        ]),
        # remaining workflow connections
        (join_sdc_echos, final_aslref_wf, [
            ("corrected", "inputnode.asl_file"),
        ]),
        (join_sdc_echos, asl_t2s_wf, [
            (("corrected_mask", pop_file), "inputnode.asl_mask"),
        ]),
    ])
    # fmt:on

    return workflow


def _to_join(in_file, join_file):
    """Join two tsv files if the join_file is not ``None``."""
    from niworkflows.interfaces.utility import JoinTSVColumns

    if join_file is None:
        return in_file
    res = JoinTSVColumns(in_file=in_file, join_file=join_file).run()
    return res.outputs.out_file


def extract_entities(file_list):
    """Return a dictionary of common entities given a list of files.

    Examples
    --------
    >>> extract_entities("sub-01/anat/sub-01_T1w.nii.gz")
    {"subject": "01", "suffix": "T1w", "datatype": "anat", "extension": ".nii.gz"}
    >>> extract_entities(["sub-01/anat/sub-01_T1w.nii.gz"] * 2)
    {"subject": "01", "suffix": "T1w", "datatype": "anat", "extension": ".nii.gz"}
    >>> extract_entities(["sub-01/anat/sub-01_run-1_T1w.nii.gz",
    ...                   "sub-01/anat/sub-01_run-2_T1w.nii.gz"])
    {"subject": "01", "run": [1, 2], "suffix": "T1w", "datatype": "anat", "extension": ".nii.gz"}

    """
    from collections import defaultdict

    from bids.layout import parse_file_entities

    entities = defaultdict(list)
    for e, v in [
        ev_pair for f in listify(file_list) for ev_pair in parse_file_entities(f).items()
    ]:
        entities[e].append(v)

    def _unique(inlist):
        inlist = sorted(set(inlist))
        if len(inlist) == 1:
            return inlist[0]
        return inlist

    return {k: _unique(v) for k, v in entities.items()}


def get_img_orientation(imgf):
    """Return the image orientation as a string."""
    img = nb.load(imgf)
    return "".join(nb.aff2axcodes(img.affine))
