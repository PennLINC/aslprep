# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for writing out derivative files."""
import typing as ty

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from niworkflows.interfaces.utility import KeySelect
from niworkflows.utils.images import dseg_label

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.utils.spaces import SpatialReferences

# from smriprep.workflows.outputs import _bids_relative  # fixed in new release, which reqs 3.10

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
        "desc": "pvGM",
        "suffix": "cbf",
    },
    "mean_cbf_wm_basil": {
        "desc": "pvWM",
        "suffix": "cbf",
    },
    "att_basil": {
        "desc": "basil",
        "suffix": "att",
    },
}


def _bids_relative(in_files, bids_root):
    from pathlib import Path

    if not isinstance(in_files, (list, tuple)):
        in_files = [in_files]
    ret = []
    for file in in_files:
        try:
            ret.append(str(Path(file).relative_to(bids_root)))
        except ValueError:
            ret.append(file)
    return in_files


def init_asl_fit_reports_wf(
    *,
    sdc_correction: bool,
    freesurfer: bool,
    output_dir: str,
    name="asl_fit_reports_wf",
) -> pe.Workflow:
    """Set up a battery of datasinks to store reports in the right location.

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
                reference_label="BOLD reference",
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


def init_ds_registration_wf(
    *,
    bids_root: str,
    output_dir: str,
    source: str,
    dest: str,
    name: str,
) -> pe.Workflow:
    """Write out registration transform.

    Copied from fMRIPrep next.
    """
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["source_files", "xform"]),
        name="inputnode",
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=["xform"]), name="outputnode")

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name="raw_sources")
    raw_sources.inputs.bids_root = bids_root

    ds_xform = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            mode="image",
            suffix="xfm",
            extension=".txt",
            dismiss_entities=("echo",),
            **{"from": source, "to": dest},
        ),
        name="ds_xform",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([
        (inputnode, raw_sources, [("source_files", "in_files")]),
        (inputnode, ds_xform, [
            ("xform", "in_file"),
            ("source_files", "source_file"),
        ]),
        (raw_sources, ds_xform, [("out", "RawSources")]),
        (ds_xform, outputnode, [("out_file", "xform")]),
    ])
    # fmt:on

    return workflow


def init_ds_hmc_wf(
    *,
    bids_root,
    output_dir,
    name="ds_hmc_wf",
) -> pe.Workflow:
    """Write out motion correction derivatives."""
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["source_files", "xforms"]),
        name="inputnode",
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=["xforms"]), name="outputnode")

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name="raw_sources")
    raw_sources.inputs.bids_root = bids_root

    ds_xforms = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc="hmc",
            suffix="xfm",
            extension=".txt",
            compress=True,
            dismiss_entities=("echo",),
            **{"from": "orig", "to": "aslref"},
        ),
        name="ds_xforms",
        run_without_submitting=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, raw_sources, [("source_files", "in_files")]),
        (inputnode, ds_xforms, [
            ("xforms", "in_file"),
            ("source_files", "source_file"),
        ]),
        (raw_sources, ds_xforms, [("out", "RawSources")]),
        (ds_xforms, outputnode, [("out_file", "xforms")]),
    ])
    # fmt:on

    return workflow


def init_ds_asl_native_wf(
    *,
    bids_root: str,
    output_dir: str,
    asl_output: bool,
    name="ds_asl_native_wf",
) -> pe.Workflow:
    """Write out aslref-space outputs."""
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "source_files",
                "asl",
                "asl_mask",
            ],
        ),
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
            ),
            name="ds_asl",
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([
            (inputnode, ds_asl, [
                ("source_files", "source_file"),
                ("asl", "in_file"),
            ]),
        ])  # fmt:skip

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
        ApplyTransforms(interpolation="MultiLabel", args="-v"),
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
                suffix="att",
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


def init_asl_derivatives_wf(
    bids_root: str,
    metadata: dict,
    spaces: SpatialReferences,
    cifti_output: bool,
    freesurfer: bool,
    project_goodvoxels: bool,
    is_multi_pld: bool,
    output_confounds=True,
    scorescrub=False,
    basil=False,
    name="asl_derivatives_wf",
):
    """Set up a battery of datasinks to store derivatives in the right location.

    Parameters
    ----------
    bids_root : :obj:`str`
        Original BIDS dataset path.
    cifti_output : :obj:`bool`
        Whether the ``--cifti-output`` flag was set.
    freesurfer : :obj:`bool`
        Whether FreeSurfer anatomical processing was run.
    project_goodvoxels : :obj:`bool`
        Whether the option was used to exclude voxels with
        locally high coefficient of variation, or that lie outside the
        cortical surfaces, from the surface projection.
    metadata : :obj:`dict`
        Metadata dictionary associated to the ASL run.
    spaces : :py:class:`~aslprep.utils.spaces.SpatialReferences`
        A container for storing, organizing, and parsing spatial normalizations. Composed of
        :py:class:`~aslprep.utils.spaces.Reference` objects representing spatial references.
        Each ``Reference`` contains a space, which is a string of either TemplateFlow template IDs
        (e.g., ``MNI152Lin``, ``MNI152NLin6Asym``, ``MNIPediatricAsym``), nonstandard references
        (e.g., ``T1w`` or ``anat``, ``sbref``, ``run``, etc.), or a custom template located in
        the TemplateFlow root directory. Each ``Reference`` may also contain a spec, which is a
        dictionary with template specifications (e.g., a specification of ``{'resolution': 2}``
        would lead to resampling on a 2mm resolution of the space).
    is_multi_pld : :obj:`bool`
        True if data are multi-delay, False otherwise.
    name : :obj:`str`
        This workflow's identifier (default: ``asl_derivatives_wf``).
    """
    nonstd_spaces = set(spaces.get_nonstandard())
    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "source_file",
                "template",
                "spatial_reference",
                "qc_file",
                # Preprocessed ASL files
                "asl_native",
                "asl_cifti",
                "aslref_native",
                "asl_mask_native",
                # Transforms
                "hmc_xforms",
                "aslref2anat_xfm",
                # Standard CBF outputs
                "cbf_ts_native",
                "cbf_ts_cifti",
                "mean_cbf_native",
                "mean_cbf_cifti",
                "att_native",
                "att_cifti",
                # SCORE/SCRUB outputs
                "cbf_ts_score_native",
                "cbf_ts_score_cifti",
                "mean_cbf_score_native",
                "mean_cbf_score_cifti",
                "mean_cbf_scrub_native",
                "mean_cbf_scrub_cifti",
                # BASIL outputs
                "mean_cbf_basil_native",
                "mean_cbf_basil_cifti",
                "mean_cbf_gm_basil_native",
                "mean_cbf_gm_basil_cifti",
                "mean_cbf_wm_basil_native",
                "mean_cbf_wm_basil_cifti",
                "att_basil_native",
                "att_basil_cifti",
                # Parcellated CBF outputs
                "atlas_names",
                "mean_cbf_parcellated",
                "mean_cbf_score_parcellated",
                "mean_cbf_scrub_parcellated",
                "mean_cbf_basil_parcellated",
                "mean_cbf_gm_basil_parcellated",
                # Non-GE outputs
                "confounds",
                "confounds_metadata",
                # CIFTI-related outputs
                "cifti_metadata",
                "cifti_density",
                "goodvoxels_mask",
                "surf_files",
                "surf_refs",
                # Freesurfer outputs
                "asl_aseg_t1",
                "asl_aseg_std",
                "asl_aparc_t1",
                "asl_aparc_std",
            ],
        ),
        name="inputnode",
    )

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name="raw_sources")
    raw_sources.inputs.bids_root = bids_root

    # fmt:off
    workflow.connect([(inputnode, raw_sources, [("source_file", "in_files")])])
    # fmt:on

    if output_confounds:
        ds_confounds = pe.Node(
            DerivativesDataSink(
                base_directory=config.execution.aslprep_dir,
                desc="confounds",
                suffix="regressors",
            ),
            name="ds_confounds",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        # fmt:off
        workflow.connect([
            (inputnode, ds_confounds, [
                ("source_file", "source_file"),
                ("confounds", "in_file"),
                ("confounds_metadata", "meta_dict"),
            ]),
        ])
        # fmt:on

    ds_qcfile = pe.Node(
        DerivativesDataSink(
            base_directory=config.execution.aslprep_dir,
            desc="qualitycontrol",
            suffix="cbf",
            compress=False,
        ),
        name="ds_qcfile",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([
        (inputnode, ds_qcfile, [
            ("source_file", "source_file"),
            ("qc_file", "in_file"),
        ]),
    ])
    # fmt:on

    ds_asl_to_t1w_xform = pe.Node(
        DerivativesDataSink(
            base_directory=config.execution.aslprep_dir,
            dismiss_entities=("echo",),
            to="T1w",
            mode="image",
            suffix="xfm",
            extension=".txt",
            **{"from": "scanner"},
        ),
        name="ds_asl_to_t1w_xform",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([
        (inputnode, ds_asl_to_t1w_xform, [
            ("source_file", "source_file"),
            ("aslref2anat_xfm", "in_file"),
        ]),
    ])
    # fmt:on

    # Parcellated data
    PARCELLATED_INPUT_FIELDS = {
        "mean_cbf_parcellated": {},
        "mean_cbf_score_parcellated": {
            "desc": "score",
        },
        "mean_cbf_scrub_parcellated": {
            "desc": "scrub",
        },
        "mean_cbf_basil_parcellated": {
            "desc": "basil",
        },
        "mean_cbf_gm_basil_parcellated": {
            "desc": "basilGM",
        },
    }
    parcellated_inputs = ["mean_cbf_parcellated"]
    if scorescrub:
        parcellated_inputs += ["mean_cbf_score_parcellated", "mean_cbf_scrub_parcellated"]

    if basil:
        parcellated_inputs += ["mean_cbf_basil_parcellated", "mean_cbf_gm_basil_parcellated"]

    for parcellated_input in parcellated_inputs:
        ds_parcellated_input = pe.MapNode(
            DerivativesDataSink(
                base_directory=config.execution.aslprep_dir,
                suffix="cbf",
                compress=False,
                **PARCELLATED_INPUT_FIELDS[parcellated_input],
            ),
            name=f"ds_{parcellated_input}",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            iterfield=["atlas", "in_file"],
        )

        # fmt:off
        workflow.connect([
            (inputnode, ds_parcellated_input, [
                ("source_file", "source_file"),
                ("atlas_names", "atlas"),
                (parcellated_input, "in_file"),
            ]),
        ])
        # fmt:on

    base_inputs = ["asl", "aslref", "asl_mask", "mean_cbf"]
    if is_multi_pld:
        # ATT is only calculated for multi-delay data
        base_inputs += ["att"]
    else:
        # CBF time series is only calculated for single-delay data
        base_inputs += ["cbf_ts"]

    if scorescrub:
        base_inputs += ["cbf_ts_score", "mean_cbf_score", "mean_cbf_scrub"]

    if basil:
        base_inputs += ["mean_cbf_basil", "mean_cbf_gm_basil", "mean_cbf_wm_basil", "att_basil"]

    # Native-space derivatives
    if nonstd_spaces.intersection(("run", "asl", "aslref", "sbref")):
        for base_input in base_inputs:
            base_input_native = f"{base_input}_native"
            ds_base_input_native = pe.Node(
                DerivativesDataSink(
                    base_directory=config.execution.aslprep_dir,
                    compress=True,
                    **BASE_INPUT_FIELDS[base_input],
                ),
                name=f"ds_{base_input_native}",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )

            # fmt:off
            workflow.connect([
                (inputnode, ds_base_input_native, [
                    ("source_file", "source_file"),
                    (base_input_native, "in_file"),
                ]),
                (raw_sources, ds_base_input_native, [("out", "RawSources")]),
            ])
            # fmt:on

    # T1w-space derivatives
    if nonstd_spaces.intersection(("T1w", "anat")):
        for base_input in base_inputs:
            base_input_t1 = f"{base_input}_t1"
            ds_base_input_t1 = pe.Node(
                DerivativesDataSink(
                    base_directory=config.execution.aslprep_dir,
                    compress=True,
                    space="T1w",
                    **BASE_INPUT_FIELDS[base_input],
                ),
                name=f"ds_{base_input_t1}",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )

            # fmt:off
            workflow.connect([
                (inputnode, ds_base_input_t1, [
                    ("source_file", "source_file"),
                    (base_input_t1, "in_file"),
                ]),
                (raw_sources, ds_base_input_t1, [("out", "RawSources")]),
            ])
            # fmt:on

        if freesurfer:
            ds_asl_aseg_t1 = pe.Node(
                DerivativesDataSink(
                    base_directory=config.execution.aslprep_dir,
                    space="T1w",
                    desc="aseg",
                    suffix="dseg",
                    compress=True,
                    dismiss_entities=("echo",),
                ),
                name="ds_asl_aseg_t1",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            ds_asl_aparc_t1 = pe.Node(
                DerivativesDataSink(
                    base_directory=config.execution.aslprep_dir,
                    space="T1w",
                    desc="aparcaseg",
                    suffix="dseg",
                    compress=True,
                    dismiss_entities=("echo",),
                ),
                name="ds_asl_aparc_t1",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            # fmt:off
            workflow.connect([
                (inputnode, ds_asl_aseg_t1, [
                    ("source_file", "source_file"),
                    ("asl_aseg_t1", "in_file"),
                ]),
                (inputnode, ds_asl_aparc_t1, [
                    ("source_file", "source_file"),
                    ("asl_aparc_t1", "in_file"),
                ]),
            ])
            # fmt:on

    if getattr(spaces, "_cached") is None:
        return workflow

    # Store resamplings in standard spaces when listed in --output-spaces
    if spaces.cached.references:
        from niworkflows.interfaces.space import SpaceDataSource

        spacesource = pe.Node(SpaceDataSource(), name="spacesource", run_without_submitting=True)
        spacesource.iterables = (
            "in_tuple",
            [(s.fullname, s.spec) for s in spaces.cached.get_standard(dim=(3,))],
        )

        select_std = pe.Node(
            KeySelect(fields=[f"{base_input}_std" for base_input in base_inputs] + ["template"]),
            name="select_std",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        # fmt:off
        workflow.connect([
            (inputnode, select_std, [
                ("template", "template"),
                ("spatial_reference", "keys"),
            ]),
            (spacesource, select_std, [("uid", "key")]),
        ])
        # fmt:on

        for base_input in base_inputs:
            base_input_std = f"{base_input}_std"
            ds_base_input_std = pe.Node(
                DerivativesDataSink(
                    base_directory=config.execution.aslprep_dir,
                    compress=True,
                    **BASE_INPUT_FIELDS[base_input],
                ),
                name=f"ds_{base_input_std}",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )

            # fmt:off
            workflow.connect([
                (inputnode, ds_base_input_std, [("source_file", "source_file")]),
                (raw_sources, ds_base_input_std, [("out", "RawSources")]),
                (inputnode, select_std, [(base_input_std, base_input_std)]),
                (select_std, ds_base_input_std, [(base_input_std, "in_file")]),
                (spacesource, ds_base_input_std, [
                    ("space", "space"),
                    ("cohort", "cohort"),
                    ("resolution", "resolution"),
                    ("density", "density"),
                ]),
            ])
            # fmt:on

        if freesurfer:
            select_fs_std = pe.Node(
                KeySelect(fields=["asl_aseg_std", "asl_aparc_std", "template"]),
                name="select_fs_std",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            ds_asl_aseg_std = pe.Node(
                DerivativesDataSink(
                    base_directory=config.execution.aslprep_dir,
                    desc="aseg",
                    suffix="dseg",
                    compress=True,
                    dismiss_entities=("echo",),
                ),
                name="ds_asl_aseg_std",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            ds_asl_aparc_std = pe.Node(
                DerivativesDataSink(
                    base_directory=config.execution.aslprep_dir,
                    desc="aparcaseg",
                    suffix="dseg",
                    compress=True,
                    dismiss_entities=("echo",),
                ),
                name="ds_asl_aparc_std",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            # fmt:off
            workflow.connect([
                (spacesource, select_fs_std, [("uid", "key")]),
                (inputnode, select_fs_std, [
                    ("asl_aseg_std", "asl_aseg_std"),
                    ("asl_aparc_std", "asl_aparc_std"),
                    ("template", "template"),
                    ("spatial_reference", "keys"),
                ]),
                (select_fs_std, ds_asl_aseg_std, [("asl_aseg_std", "in_file")]),
                (spacesource, ds_asl_aseg_std, [
                    ("space", "space"),
                    ("cohort", "cohort"),
                    ("resolution", "resolution"),
                    ("density", "density"),
                ]),
                (select_fs_std, ds_asl_aparc_std, [("asl_aparc_std", "in_file")]),
                (spacesource, ds_asl_aparc_std, [
                    ("space", "space"),
                    ("cohort", "cohort"),
                    ("resolution", "resolution"),
                    ("density", "density"),
                ]),
                (inputnode, ds_asl_aseg_std, [("source_file", "source_file")]),
                (inputnode, ds_asl_aparc_std, [("source_file", "source_file")])
            ])
            # fmt:on

    fs_outputs = spaces.cached.get_fs_spaces()
    if freesurfer and fs_outputs:
        from niworkflows.interfaces.surf import Path2BIDS

        select_fs_surf = pe.Node(
            KeySelect(fields=["surfaces", "surf_kwargs"]),
            name="select_fs_surf",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        select_fs_surf.iterables = [("key", fs_outputs)]
        select_fs_surf.inputs.surf_kwargs = [{"space": s} for s in fs_outputs]

        name_surfs = pe.MapNode(
            Path2BIDS(pattern=r"(?P<hemi>[lr])h.\w+"),
            iterfield="in_file",
            name="name_surfs",
            run_without_submitting=True,
        )

        ds_asl_surfs = pe.MapNode(
            DerivativesDataSink(
                base_directory=config.execution.aslprep_dir,
                extension=".asl.gii",
                TaskName=metadata.get("TaskName"),
            ),
            iterfield=["in_file", "hemi"],
            name="ds_asl_surfs",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        # fmt:off
        workflow.connect([
            (inputnode, select_fs_surf, [
                ("surf_files", "surfaces"),
                ("surf_refs", "keys"),
            ]),
            (select_fs_surf, name_surfs, [("surfaces", "in_file")]),
            (inputnode, ds_asl_surfs, [("source_file", "source_file")]),
            (select_fs_surf, ds_asl_surfs, [
                ("surfaces", "in_file"),
                ("key", "space"),
            ]),
            (name_surfs, ds_asl_surfs, [("hemi", "hemi")]),
        ])
        # fmt:on

    if freesurfer and project_goodvoxels:
        ds_goodvoxels_mask = pe.Node(
            DerivativesDataSink(
                base_directory=config.execution.aslprep_dir,
                space="T1w",
                desc="goodvoxels",
                suffix="mask",
                Type="ROI",  # Metadata
                compress=True,
                dismiss_entities=("echo",),
            ),
            name="ds_goodvoxels_mask",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        # fmt:off
        workflow.connect([
            (inputnode, ds_goodvoxels_mask, [
                ("source_file", "source_file"),
                ("goodvoxels_mask", "in_file"),
            ]),
        ])
        # fmt:on

    # CIFTI output
    if cifti_output:
        ds_asl_cifti = pe.Node(
            DerivativesDataSink(
                base_directory=config.execution.aslprep_dir,
                suffix="asl",
                compress=False,
                TaskName=metadata.get("TaskName"),
                space="fsLR",
            ),
            name="ds_asl_cifti",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        # fmt:off
        workflow.connect([
            (inputnode, ds_asl_cifti, [
                ("asl_cifti", "in_file"),
                ("source_file", "source_file"),
                ("cifti_density", "density"),
                (("cifti_metadata", _read_json), "meta_dict"),
            ]),
        ])
        # fmt:on

    return workflow


def _read_json(in_file):
    from json import loads
    from pathlib import Path

    return loads(Path(in_file).read_text())
