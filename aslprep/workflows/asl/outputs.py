# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for writing out derivative files."""
import typing as ty

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.utils.images import dseg_label
from smriprep.workflows.outputs import _bids_relative

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.interfaces.ants import ApplyTransforms

BASE_INPUT_FIELDS = {
    "asl": {
        "desc": "preproc",
        "suffix": "asl",
    },
    "aslref": {
        "suffix": "aslref",
    },
    "asl_mask": {
        "desc": "brain",
        "suffix": "mask",
    },
    # CBF outputs
    "cbf_ts": {
        "desc": "timeseries",
        "suffix": "cbf",
    },
    "mean_cbf": {
        "suffix": "cbf",
    },
    "att": {
        "suffix": "att",
    },
    # SCORE/SCRUB outputs
    "cbf_ts_score": {
        "desc": "scoreTimeseries",
        "suffix": "cbf",
    },
    "mean_cbf_score": {
        "desc": "score",
        "suffix": "cbf",
    },
    "mean_cbf_scrub": {
        "desc": "scrub",
        "suffix": "cbf",
    },
    # BASIL outputs
    "mean_cbf_basil": {
        "desc": "basil",
        "suffix": "cbf",
    },
    "mean_cbf_gm_basil": {
        "desc": "basilGM",
        "suffix": "cbf",
    },
    "mean_cbf_wm_basil": {
        "desc": "basilWM",
        "suffix": "cbf",
    },
    "att_basil": {
        "desc": "basil",
        "suffix": "att",
    },
}


def init_asl_fit_reports_wf(
    *,
    sdc_correction: bool,
    freesurfer: bool,  # noqa:U100
    output_dir: str,
    name="asl_fit_reports_wf",
) -> pe.Workflow:
    """Set up a battery of datasinks to store reports in the right location.

    Copied from fMRIPrep's init_bold_fit_reports_wf.
    Modifications include changes to fields and variable names, which aren't important,
    and DerivativesDataSink suffixes, which are.

    Parameters
    ----------
    freesurfer : :obj:`bool`
        FreeSurfer was enabled
    output_dir : :obj:`str`
        Directory in which to save derivatives
    name : :obj:`str`
        Workflow name (default: anat_reports_wf)

    Inputs
    ------
    source_file
        Input BOLD images

    std_t1w
        T1w image resampled to standard space
    std_mask
        Mask of skull-stripped template
    subject_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    t1w_conform_report
        Conformation report
    t1w_preproc
        The T1w reference map, which is calculated as the average of bias-corrected
        and preprocessed T1w images, defining the anatomical space.
    t1w_dseg
        Segmentation in T1w space
    t1w_mask
        Brain (binary) mask estimated by brain extraction.
    template
        Template space and specifications

    """
    from niworkflows.interfaces.reportlets.registration import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )
    from sdcflows.interfaces.reportlets import FieldmapReportlet

    workflow = pe.Workflow(name=name)

    inputfields = [
        "source_file",
        "sdc_aslref",
        "coreg_aslref",
        "aslref2anat_xfm",
        "aslref2fmap_xfm",
        "t1w_preproc",
        "t1w_mask",
        "t1w_dseg",
        "fieldmap",
        "fmap_ref",
        # May be missing
        "subject_id",
        "subjects_dir",
        # Report snippets
        "summary_report",
        "validation_report",
    ]
    inputnode = pe.Node(niu.IdentityInterface(fields=inputfields), name="inputnode")

    ds_summary = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc="summary",
            datatype="figures",
            dismiss_entities=("echo",),
        ),
        name="ds_report_summary",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    ds_validation = pe.Node(
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

    # Resample anatomical references into BOLD space for plotting
    t1w_aslref = pe.Node(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            float=True,
            invert_transform_flags=[True],
            interpolation="LanczosWindowedSinc",
            args="-v",
        ),
        name="t1w_aslref",
        mem_gb=1,
    )

    t1w_wm = pe.Node(
        niu.Function(function=dseg_label),
        name="t1w_wm",
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    t1w_wm.inputs.label = 2  # BIDS default is WM=2

    aslref_wm = pe.Node(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            invert_transform_flags=[True],
            interpolation="NearestNeighbor",
            args="-v",
        ),
        name="aslref_wm",
        mem_gb=1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, ds_summary, [
            ("source_file", "source_file"),
            ("summary_report", "in_file"),
        ]),
        (inputnode, ds_validation, [
            ("source_file", "source_file"),
            ("validation_report", "in_file"),
        ]),
        (inputnode, t1w_aslref, [
            ("t1w_preproc", "input_image"),
            ("coreg_aslref", "reference_image"),
            ("aslref2anat_xfm", "transforms"),
        ]),
        (inputnode, t1w_wm, [("t1w_dseg", "in_seg")]),
        (inputnode, aslref_wm, [
            ("coreg_aslref", "reference_image"),
            ("aslref2anat_xfm", "transforms"),
        ]),
        (t1w_wm, aslref_wm, [("out", "input_image")]),
    ])
    # fmt:on

    # Reportlets follow the structure of init_asl_fit_wf stages
    # - SDC1:
    #       Before: Pre-SDC aslref
    #       After: Fieldmap reference resampled on aslref
    #       Three-way: Fieldmap resampled on aslref
    # - SDC2:
    #       Before: Pre-SDC aslref with white matter mask
    #       After: Post-SDC aslref with white matter mask
    # - EPI-T1 registration:
    #       Before: T1w brain with white matter mask
    #       After: Resampled aslref with white matter mask

    if sdc_correction:
        fmapref_aslref = pe.Node(
            ApplyTransforms(
                dimension=3,
                default_value=0,
                float=True,
                invert_transform_flags=[True],
                interpolation="LanczosWindowedSinc",
                args="-v",
            ),
            name="fmapref_aslref",
            mem_gb=1,
        )

        # SDC1
        sdcreg_report = pe.Node(
            FieldmapReportlet(
                reference_label="ASL reference",
                moving_label="Fieldmap reference",
                show="both",
            ),
            name="sdecreg_report",
            mem_gb=0.1,
        )

        ds_sdcreg_report = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc="fmapCoreg",
                suffix="asl",
                datatype="figures",
                dismiss_entities=("echo",),
            ),
            name="ds_sdcreg_report",
        )

        # SDC2
        sdc_report = pe.Node(
            SimpleBeforeAfter(
                before_label="Distorted",
                after_label="Corrected",
                dismiss_affine=True,
            ),
            name="sdc_report",
            mem_gb=0.1,
        )

        ds_sdc_report = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc="sdc",
                suffix="asl",
                datatype="figures",
                dismiss_entities=("echo",),
            ),
            name="ds_sdc_report",
        )

        # fmt:off
        workflow.connect([
            (inputnode, fmapref_aslref, [
                ("fmap_ref", "input_image"),
                ("coreg_aslref", "reference_image"),
                ("aslref2fmap_xfm", "transforms"),
            ]),
            (inputnode, sdcreg_report, [
                ("sdc_aslref", "reference"),
                ("fieldmap", "fieldmap")
            ]),
            (fmapref_aslref, sdcreg_report, [("output_image", "moving")]),
            (inputnode, ds_sdcreg_report, [("source_file", "source_file")]),
            (sdcreg_report, ds_sdcreg_report, [("out_report", "in_file")]),
            (inputnode, sdc_report, [
                ("sdc_aslref", "before"),
                ("coreg_aslref", "after"),
            ]),
            (aslref_wm, sdc_report, [("output_image", "wm_seg")]),
            (inputnode, ds_sdc_report, [("source_file", "source_file")]),
            (sdc_report, ds_sdc_report, [("out_report", "in_file")]),
        ])
        # fmt:on

    # EPI-T1 registration
    # Resample T1w image onto EPI-space

    epi_t1_report = pe.Node(
        SimpleBeforeAfter(
            before_label="T1w",
            after_label="EPI",
            dismiss_affine=True,
        ),
        name="epi_t1_report",
        mem_gb=0.1,
    )

    ds_epi_t1_report = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc="coreg",
            suffix="asl",
            datatype="figures",
            dismiss_entities=("echo",),
        ),
        name="ds_epi_t1_report",
    )

    # fmt:off
    workflow.connect([
        (inputnode, epi_t1_report, [("coreg_aslref", "after")]),
        (t1w_aslref, epi_t1_report, [("output_image", "before")]),
        (aslref_wm, epi_t1_report, [("output_image", "wm_seg")]),
        (inputnode, ds_epi_t1_report, [("source_file", "source_file")]),
        (epi_t1_report, ds_epi_t1_report, [("out_report", "in_file")]),
    ])
    # fmt:on

    return workflow


def init_ds_aslref_wf(
    *,
    bids_root,
    output_dir,
    desc: str,
    name="ds_aslref_wf",
) -> pe.Workflow:
    """Write out aslref image."""
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["source_files", "aslref"]),
        name="inputnode",
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=["aslref"]), name="outputnode")

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name="raw_sources")
    raw_sources.inputs.bids_root = bids_root

    ds_aslref = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc=desc,
            suffix="aslref",
            compress=True,
            dismiss_entities=("echo",),
        ),
        name="ds_aslref",
        run_without_submitting=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, raw_sources, [("source_files", "in_files")]),
        (inputnode, ds_aslref, [
            ("aslref", "in_file"),
            ("source_files", "source_file"),
        ]),
        (raw_sources, ds_aslref, [("out", "RawSources")]),
        (ds_aslref, outputnode, [("out_file", "aslref")]),
    ])
    # fmt:on

    return workflow


def init_ds_asl_native_wf(
    *,
    bids_root: str,
    output_dir: str,
    asl_output: bool,
    metadata: ty.List[dict],
    cbf_3d: ty.List[str],
    cbf_4d: ty.List[str],
    att: ty.List[str],
    name="ds_asl_native_wf",
) -> pe.Workflow:
    """Write out aslref-space outputs."""
    workflow = pe.Workflow(name=name)

    inputnode_fields = [
        "source_files",
        "asl",
        "asl_mask",
    ]
    inputnode_fields += cbf_3d
    inputnode_fields += cbf_4d
    inputnode_fields += att
    inputnode = pe.Node(
        niu.IdentityInterface(fields=inputnode_fields),
        name="inputnode",
    )

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name="raw_sources")
    raw_sources.inputs.bids_root = bids_root
    workflow.connect([(inputnode, raw_sources, [("source_files", "in_files")])])

    # Masks should be output if any other derivatives are output
    ds_asl_mask = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc="brain",
            suffix="mask",
            compress=True,
            dismiss_entities=("echo",),
        ),
        name="ds_asl_mask",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([
        (inputnode, ds_asl_mask, [
            ("source_files", "source_file"),
            ("asl_mask", "in_file"),
        ]),
        (raw_sources, ds_asl_mask, [("out", "RawSources")]),
    ])  # fmt:skip

    if asl_output:
        ds_asl = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc="preproc",
                compress=True,
                SkullStripped=False,
                dismiss_entities=("echo",),
                **metadata,
            ),
            name="ds_asl",
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(inputnode, ds_asl, [("asl", "in_file")])])
        datasinks = [ds_asl]

        for cbf_name in cbf_4d + cbf_3d:
            # TODO: Add EstimationReference and EstimationAlgorithm
            cbf_meta = {
                "Units": "mL/100 g/min",
            }
            fields = BASE_INPUT_FIELDS[cbf_name]

            ds_cbf = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    compress=True,
                    dismiss_entities=("echo",),
                    **fields,
                    **cbf_meta,
                ),
                name=f"ds_{cbf_name}",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            datasinks.append(ds_cbf)
            workflow.connect([(inputnode, ds_cbf, [(cbf_name, "in_file")])])

        for att_name in att:
            # TODO: Add EstimationReference and EstimationAlgorithm
            att_meta = {
                "Units": "s",
            }
            fields = BASE_INPUT_FIELDS[att_name]

            ds_att = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    compress=True,
                    dismiss_entities=("echo",),
                    **fields,
                    **att_meta,
                ),
                name=f"ds_{att_name}",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            datasinks.append(ds_att)

            workflow.connect([(inputnode, ds_att, [(att_name, "in_file")])])

        workflow.connect(
            [
                (inputnode, datasink, [("source_files", "source_file")]) for datasink in datasinks
            ] + [
                (raw_sources, datasink, [("out", "RawSources")]) for datasink in datasinks
            ]
        )  # fmt:skip

    return workflow


def init_ds_volumes_wf(
    *,
    bids_root: str,
    output_dir: str,
    metadata: ty.List[dict],
    cbf_3d: ty.List[str],
    cbf_4d: ty.List[str],
    att: ty.List[str],
    name: str = "ds_volumes_wf",
) -> pe.Workflow:
    """Apply transforms from reference to anatomical/standard space and write out derivatives."""
    workflow = pe.Workflow(name=name)
    inputnode_fields = [
        "source_files",
        "ref_file",
        "asl",  # Resampled into target space
        "asl_mask",  # aslref space
        "aslref",  # aslref space
        # Anatomical
        "aslref2anat_xfm",
        # Template
        "anat2std_xfm",
        # Entities
        "space",
        "cohort",
        "resolution",
    ]
    inputnode_fields += cbf_3d
    inputnode_fields += cbf_4d
    inputnode_fields += att
    inputnode = pe.Node(
        niu.IdentityInterface(fields=inputnode_fields),
        name="inputnode",
    )

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name="raw_sources")
    raw_sources.inputs.bids_root = bids_root
    aslref2target = pe.Node(niu.Merge(2), name="aslref2target")

    # BOLD is pre-resampled
    ds_asl = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc="preproc",
            compress=True,
            SkullStripped=True,
            dismiss_entities=("echo",),
            **metadata,
        ),
        name="ds_asl",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([
        (inputnode, raw_sources, [("source_files", "in_files")]),
        (inputnode, aslref2target, [
            ("aslref2anat_xfm", "in1"),
            ("anat2std_xfm", "in2"),
        ]),
        (inputnode, ds_asl, [
            ("source_files", "source_file"),
            ("asl", "in_file"),
            ("space", "space"),
            ("cohort", "cohort"),
            ("resolution", "resolution"),
        ]),
    ])  # fmt:skip

    resample_ref = pe.Node(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            float=True,
            interpolation="LanczosWindowedSinc",
            args="-v",
        ),
        name="resample_ref",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    resample_mask = pe.Node(
        ApplyTransforms(interpolation="GenericLabel", args="-v"),
        name="resample_mask",
    )
    resamplers = [resample_ref, resample_mask]

    workflow.connect([
        (inputnode, resample_ref, [("aslref", "input_image")]),
        (inputnode, resample_mask, [("asl_mask", "input_image")]),
    ])  # fmt:skip

    ds_ref = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            suffix="aslref",
            compress=True,
            dismiss_entities=("echo",),
        ),
        name="ds_ref",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    ds_mask = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc="brain",
            suffix="mask",
            compress=True,
            dismiss_entities=("echo",),
        ),
        name="ds_mask",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    datasinks = [ds_ref, ds_mask]

    for cbf_name in cbf_4d + cbf_3d:
        # TODO: Add EstimationReference and EstimationAlgorithm
        cbf_meta = {
            "Units": "mL/100 g/min",
        }
        fields = BASE_INPUT_FIELDS[cbf_name]

        kwargs = {}
        if cbf_name in cbf_4d:
            kwargs["dimension"] = 3

        resample_cbf = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
                args="-v",
                **kwargs,
            ),
            name=f"warp_{cbf_name}_to_std",
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        ds_cbf = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                compress=True,
                dismiss_entities=("echo",),
                **fields,
                **cbf_meta,
            ),
            name=f"ds_{cbf_name}",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        resamplers.append(resample_cbf)
        datasinks.append(ds_cbf)

        workflow.connect([(inputnode, resample_cbf, [(cbf_name, "input_image")])])

    for att_name in att:
        # TODO: Add EstimationReference and EstimationAlgorithm
        att_meta = {
            "Units": "s",
        }
        fields = BASE_INPUT_FIELDS[att_name]

        resample_att = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
                args="-v",
            ),
            name=f"warp_{att_name}_to_std",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        ds_att = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                compress=True,
                dismiss_entities=("echo",),
                **fields,
                **att_meta,
            ),
            name=f"ds_{att_name}",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        resamplers.append(resample_att)
        datasinks.append(ds_att)

        workflow.connect([(inputnode, resample_att, [(att_name, "input_image")])])

    workflow.connect(
        [
            (inputnode, resampler, [("ref_file", "reference_image")])
            for resampler in resamplers
        ] + [
            (aslref2target, resampler, [("out", "transforms")])
            for resampler in resamplers
        ] + [
            (inputnode, datasink, [
                ("source_files", "source_file"),
                ("space", "space"),
                ("cohort", "cohort"),
                ("resolution", "resolution"),
            ])
            for datasink in datasinks
        ] + [
            (resampler, datasink, [("output_image", "in_file")])
            for resampler, datasink in zip(resamplers, datasinks)
        ]
    )  # fmt:skip

    return workflow


def init_ds_ciftis_wf(
    *,
    bids_root: str,
    output_dir: str,
    metadata: ty.List[dict],
    cbf_3d: ty.List[str],
    cbf_4d: ty.List[str],
    att: ty.List[str],
    name: str = "ds_ciftis_wf",
) -> pe.Workflow:
    """Apply transforms from reference to fsLR space and write out derivatives."""
    from fmriprep.workflows.bold.resampling import (
        init_bold_fsLR_resampling_wf,
        init_bold_grayords_wf,
    )

    workflow = pe.Workflow(name=name)
    inputnode_fields = [
        "asl_cifti",
        "source_files",
        # Anatomical
        "anat",
        "aslref2anat_xfm",
        # Template
        "anat2mni6_xfm",
        # Pre-computed goodvoxels mask. May be Undefined.
        "goodvoxels_mask",
        # Other inputs
        "white",
        "pial",
        "midthickness",
        "midthickness_fsLR",
        "sphere_reg_fsLR",
        "cortex_mask",
        "anat_ribbon",
    ]
    inputnode_fields += cbf_3d
    inputnode_fields += cbf_4d
    inputnode_fields += att
    inputnode = pe.Node(
        niu.IdentityInterface(fields=inputnode_fields),
        name="inputnode",
    )

    outputnode_fields = []
    outputnode_fields += cbf_3d
    outputnode_fields += cbf_4d
    outputnode_fields += att
    outputnode = pe.Node(
        niu.IdentityInterface(fields=outputnode_fields),
        name="outputnode",
    )

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name="raw_sources")
    raw_sources.inputs.bids_root = bids_root
    workflow.connect([(inputnode, raw_sources, [("source_files", "in_files")])])

    ds_asl_cifti = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            space="fsLR",
            density=config.workflow.cifti_output,
            suffix="asl",
            extension="dtseries.nii",
            compress=False,
        ),
        name="ds_asl_cifti",
        run_without_submitting=True,
    )
    workflow.connect([(inputnode, ds_asl_cifti, [("asl_cifti", "in_file")])])

    for cbf_deriv in cbf_4d + cbf_3d + att:
        kwargs = {}
        extension = "dscalar.nii"
        if cbf_deriv in cbf_4d:
            kwargs["dimension"] = 3
            extension = "dtseries.nii"

        if cbf_deriv in att:
            meta = {"Units": "s"}
        else:
            meta = {"Units": "mL/100 g/min"}

        warp_cbf_to_anat = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
                args="-v",
                **kwargs,
            ),
            name=f"warp_{cbf_deriv}_to_anat",
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([
            (inputnode, warp_cbf_to_anat, [
                (cbf_deriv, "input_image"),
                ("anat", "reference_image"),
                ("aslref2anat_xfm", "transforms"),
            ]),
        ])  # fmt:skip

        aslref2MNI6 = pe.Node(niu.Merge(2), name="aslref2MNI6")
        workflow.connect([
            (inputnode, aslref2MNI6, [
                ("aslref2anat_xfm", "in1"),
                ("anat2mni6_xfm", "in2"),
            ]),
        ])  # fmt:skip

        warp_cbf_to_MNI6 = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
                args="-v",
                **kwargs,
            ),
            name=f"warp_{cbf_deriv}_to_MNI6",
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([
            (inputnode, warp_cbf_to_MNI6, [
                ("mni6_mask", "reference_image"),
                (cbf_deriv, "input_image"),
            ]),
            (aslref2MNI6, warp_cbf_to_MNI6, [("out", "transforms")]),
        ])  # fmt:skip

        # XXX: Need to add predefined "goodvoxels_mask" support.
        cbf_fsLR_resampling_wf = init_bold_fsLR_resampling_wf(
            estimate_goodvoxels=False,  # already made from ASL data
            grayord_density=config.workflow.cifti_output,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            name=f"{cbf_deriv}_fsLR_resampling_wf",
        )
        workflow.connect([
            # Resample T1w-space CBF to fsLR surfaces
            (inputnode, cbf_fsLR_resampling_wf, [
                ("white", "inputnode.white"),
                ("pial", "inputnode.pial"),
                ("midthickness", "inputnode.midthickness"),
                ("midthickness_fsLR", "inputnode.midthickness_fsLR"),
                ("sphere_reg_fsLR", "inputnode.sphere_reg_fsLR"),
                ("cortex_mask", "inputnode.cortex_mask"),
                ("anat_ribbon", "inputnode.anat_ribbon"),
                ("goodvoxels_mask", "inputnode.goodvoxels_mask"),
            ]),
            (warp_cbf_to_anat, cbf_fsLR_resampling_wf, [("output_image", "inputnode.bold_file")])
        ])  # fmt:skip

        cbf_grayords_wf = init_bold_grayords_wf(
            grayord_density=config.workflow.cifti_output,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            repetition_time=metadata["RepetitionTime"],
            name=f"{cbf_deriv}_grayords_wf",
        )
        workflow.connect([
            (warp_cbf_to_MNI6, cbf_grayords_wf, [("output_image", "inputnode.bold_std")]),
            (cbf_fsLR_resampling_wf, cbf_grayords_wf, [
                ("outputnode.bold_fsLR", "inputnode.bold_fsLR"),
            ]),
        ])  # fmt:skip

        ds_cbf_cifti = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                space="fsLR",
                density=config.workflow.cifti_output,
                extension=extension,
                compress=False,
                **BASE_INPUT_FIELDS[cbf_deriv],
                **meta,
            ),
            name=f"ds_{cbf_deriv}_cifti",
            run_without_submitting=True,
        )
        workflow.connect([
            (raw_sources, ds_cbf_cifti, [("out", "RawSources")]),
            (cbf_grayords_wf, ds_cbf_cifti, [
                ("outputnode.cifti_bold", "in_file"),
                (("outputnode.cifti_metadata", _read_json), "meta_dict"),
            ]),
            (ds_cbf_cifti, outputnode, [("out_file", cbf_deriv)])
        ])  # fmt:skip

    return workflow


def _read_json(in_file):
    from json import loads
    from pathlib import Path

    if not isinstance(in_file, str):
        raise ValueError(f"_read_json: input is not str ({in_file})")

    return loads(Path(in_file).read_text())
