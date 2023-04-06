# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for writing out derivative files."""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from smriprep.workflows.outputs import _bids_relative

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.interfaces.utility import KeySelect


def init_asl_derivatives_wf(
    bids_root,
    metadata,
    output_dir,
    spaces,
    scorescrub=False,
    basil=False,
    name="asl_derivatives_wf",
):
    """Set up a battery of datasinks to store derivatives in the right location.

    Parameters
    ----------
    bids_root : :obj:`str`
        Original BIDS dataset path.
    metadata : :obj:`dict`
        Metadata dictionary associated to the ASL run.
    output_dir : :obj:`str`
        Where derivatives should be written out to.
    spaces : :py:class:`~niworkflows.utils.spaces.SpatialReferences`
        A container for storing, organizing, and parsing spatial normalizations. Composed of
        :py:class:`~niworkflows.utils.spaces.Reference` objects representing spatial references.
        Each ``Reference`` contains a space, which is a string of either TemplateFlow template IDs
        (e.g., ``MNI152Lin``, ``MNI152NLin6Asym``, ``MNIPediatricAsym``), nonstandard references
        (e.g., ``T1w`` or ``anat``, ``sbref``, ``run``, etc.), or a custom template located in
        the TemplateFlow root directory. Each ``Reference`` may also contain a spec, which is a
        dictionary with template specifications (e.g., a specification of ``{'resolution': 2}``
        would lead to resampling on a 2mm resolution of the space).
    name : :obj:`str`
        This workflow's identifier (default: ``func_derivatives_wf``).
    """
    nonstd_spaces = set(spaces.get_nonstandard())
    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "confounds",
                "confounds_metadata",
                "source_file",
                "template",
                "spatial_reference",
                "qc_file",
                # Preprocessed ASL files
                "asl_native",
                "asl_t1",
                "asl_std",
                "asl_native_ref",
                "asl_t1_ref",
                "asl_std_ref",
                "asl_mask_native",
                "asl_mask_t1",
                "asl_mask_std",
                # Transforms
                "itk_asl_to_t1",
                "itk_t1_to_asl",
                # Standard CBF outputs
                "cbf",
                "cbf_t1",
                "cbf_std",
                "meancbf",
                "meancbf_t1",
                "meancbf_std",
                # SCORE/SCRUB outputs
                "score",
                "score_t1",
                "score_std",
                "avgscore",
                "avgscore_t1",
                "avgscore_std",
                "scrub",
                "scrub_t1",
                "scrub_std",
                # BASIL outputs
                "basil",
                "basil_t1",
                "basil_std",
                "pv",
                "pv_t1",
                "pv_std",
                "pvwm",
                "pvwm_t1",
                "pvwm_std",
                "att",
                "att_t1",
                "att_std",
                # Parcellated CBF outputs
                "atlas_names",
                "mean_cbf_parcellated",
                "mean_cbf_score_parcellated",
                "mean_cbf_scrub_parcellated",
                "mean_cbf_basil_parcellated",
                "mean_cbf_gm_basil_parcellated",
            ]
        ),
        name="inputnode",
    )

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name="raw_sources")
    raw_sources.inputs.bids_root = bids_root

    ds_confounds = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc="confounds",
            suffix="regressors",
        ),
        name="ds_confounds",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    # fmt:off
    workflow.connect([
        (inputnode, raw_sources, [("source_file", "in_files")]),
        (inputnode, ds_confounds, [
            ("source_file", "source_file"),
            ("confounds", "in_file"),
            ("confounds_metadata", "meta_dict"),
        ]),
    ])
    # fmt:on

    qcfile = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc="qualitycontrol",
            suffix="cbf",
            compress=False,
        ),
        name="qcfile",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    # fmt:off
    workflow.connect([
        (inputnode, qcfile, [
            ("source_file", "source_file"),
            ("qc_file", "in_file"),
        ]),
    ])
    # fmt:on

    # write transform matrix file between asl native space and T1w
    ds_t1w_to_asl_xform = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            to="T1w",
            mode="image",
            suffix="xfm",
            extension=".txt",
            dismiss_entities=("echo",),
            **{"from": "scanner"},
        ),
        name="ds_t1w_to_asl_xform",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([
        (inputnode, ds_t1w_to_asl_xform, [
            ("source_file", "source_file"),
            ("itk_asl_to_t1", "in_file"),
        ]),
    ])
    # fmt:on

    ds_asl_to_t1w_xform = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            dismiss_entities=("echo",),
            to="scanner",
            mode="image",
            suffix="xfm",
            extension=".txt",
            **{"from": "T1w"},
        ),
        name="ds_asl_to_t1w_xform",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([
        (inputnode, ds_asl_to_t1w_xform, [
            ("source_file", "source_file"),
            ("itk_t1_to_asl", "in_file"),
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
                base_directory=output_dir,
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

    timeseries_metadata = {
        "SkullStripped": False,
        "RepetitionTime": metadata.get("RepetitionTime"),
        "RepetitionTimePreparation": metadata.get("RepetitionTimePreparation"),
        "TaskName": metadata.get("TaskName"),
    }

    BASE_INPUT_FIELDS = {
        "asl": {
            "desc": "preproc",
            "suffix": "asl",
            **timeseries_metadata,
        },
        "aslref": {
            "suffix": "aslref",
            "dismiss_entities": ("echo",),
        },
        "asl_mask": {
            "desc": "brain",
            "suffix": "mask",
            "dismiss_entities": ("echo",),
        },
        "cbf": {
            "suffix": "cbf",
        },
        "meancbf": {
            "desc": "mean",
            "suffix": "cbf",
        },
        "score": {
            "desc": "score",
            "suffix": "cbf",
        },
        "avgscore": {
            "desc": "meanScore",
            "suffix": "cbf",
        },
        "scrub": {
            "desc": "scrub",
            "suffix": "cbf",
        },
        "basil": {
            "desc": "basil",
            "suffix": "cbf",
        },
        "pv": {
            "desc": "pvGM",
            "suffix": "cbf",
        },
        "pvwm": {
            "desc": "pvWM",
            "suffix": "cbf",
        },
        "att": {
            "desc": "bat",
            "suffix": "cbf",
        },
    }

    base_inputs = ["asl", "aslref", "asl_mask", "meancbf"]
    if scorescrub:
        base_inputs += ["score", "avgscore", "scrub"]

    if basil:
        base_inputs += ["basil", "pv", "pvwm", "att"]

    if nonstd_spaces.intersection(("func", "run", "asl", "sbref")):
        for base_input in base_inputs:
            ds_base_input_native = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    compress=True,
                    **BASE_INPUT_FIELDS[base_input],
                ),
                name=f"ds_{base_input}_native",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )

            # fmt:off
            workflow.connect([
                (inputnode, ds_base_input_native, [
                    ("source_file", "source_file"),
                    (f"{base_input}_native", "in_file"),
                ]),
            ])
            # fmt:on

            # The mask file gets one extra bit of metadata.
            if base_input == "asl_mask":
                # fmt:off
                workflow.connect([(raw_sources, ds_base_input_native, [("out", "RawSources")])])
                # fmt:on

    # Resample to T1w space
    if nonstd_spaces.intersection(("T1w", "anat")):
        for base_input in base_inputs:
            ds_base_input_t1 = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    compress=True,
                    space="T1w",
                    **BASE_INPUT_FIELDS[base_input],
                ),
                name=f"ds_{base_input}_t1",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )

            # fmt:off
            workflow.connect([
                (inputnode, ds_base_input_t1, [
                    ("source_file", "source_file"),
                    (f"{base_input}_t1", "in_file"),
                ]),
            ])
            # fmt:on

            # The mask file gets one extra bit of metadata.
            if base_input == "asl_mask":
                # fmt:off
                workflow.connect([(raw_sources, ds_base_input_t1, [("out", "RawSources")])])
                # fmt:on

    if getattr(spaces, "_cached") is None:
        return workflow

    # Store resamplings in standard spaces when listed in --output-spaces
    if spaces.cached.references:
        from niworkflows.interfaces.space import SpaceDataSource

        spacesource = pe.Node(SpaceDataSource(), name="spacesource", run_without_submitting=True)
        config.loggers.interface.warning(spaces.cached.get_standard(dim=(3,)))
        spacesource.iterables = (
            "in_tuple",
            [(s.fullname, s.spec) for s in spaces.cached.get_standard(dim=(3,))],
        )
        out_names = [
            "template",
            "asl_std",
            "asl_std_ref",
            "asl_mask_std",
            "cbf_std",
            "meancbf_std",
        ]
        if scorescrub:
            out_names += ["score_std", "avgscore_std", "scrub_std"]

        if basil:
            out_names += ["basil_std", "pv_std", "pvwm_std", "att_std"]

        select_std = pe.Node(
            KeySelect(fields=out_names),
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
            ds_base_input_std = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    compress=True,
                    **BASE_INPUT_FIELDS[base_input],
                ),
                name=f"ds_{base_input}_std",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )

            # fmt:off
            workflow.connect([
                (inputnode, select_std, [(f"{base_input}_std", f"{base_input}_std")]),
                (inputnode, ds_base_input_std, [
                    ("source_file", "source_file"),
                    (f"{base_input}_std", "in_file"),
                ]),
                (select_std, ds_base_input_std, [(f"{base_input}_std", "in_file")]),
                (spacesource, ds_base_input_std, [
                    ("space", "space"),
                    ("cohort", "cohort"),
                    ("resolution", "resolution"),
                    ("density", "density"),
                ]),
            ])
            # fmt:on

            # The mask file gets one extra bit of metadata.
            if base_input == "asl_mask":
                # fmt:off
                workflow.connect([(raw_sources, ds_base_input_std, [("out", "RawSources")])])
                # fmt:on

    return workflow


def init_geasl_derivatives_wf(
    bids_root,
    metadata,
    output_dir,
    spaces,
    scorescrub=False,
    basil=False,
    name="asl_derivatives_wf",
):
    """Set up a battery of datasinks to store derivatives in the right location, for GE data.

    Parameters
    ----------
    bids_root : :obj:`str`
        Original BIDS dataset path.
    cifti_output : :obj:`bool`
        Whether the ``--cifti-output`` flag was set.
    freesurfer : :obj:`bool`
        Whether FreeSurfer anatomical processing was run.
    metadata : :obj:`dict`
        Metadata dictionary associated to the ASL run.
    output_dir : :obj:`str`
        Where derivatives should be written out to.
    spaces : :py:class:`~niworkflows.utils.spaces.SpatialReferences`
        A container for storing, organizing, and parsing spatial normalizations. Composed of
        :py:class:`~niworkflows.utils.spaces.Reference` objects representing spatial references.
        Each ``Reference`` contains a space, which is a string of either TemplateFlow template IDs
        (e.g., ``MNI152Lin``, ``MNI152NLin6Asym``, ``MNIPediatricAsym``), nonstandard references
        (e.g., ``T1w`` or ``anat``, ``sbref``, ``run``, etc.), or a custom template located in
        the TemplateFlow root directory. Each ``Reference`` may also contain a spec, which is a
        dictionary with template specifications (e.g., a specification of ``{'resolution': 2}``
        would lead to resampling on a 2mm resolution of the space).
    name : :obj:`str`
        This workflow's identifier (default: ``func_derivatives_wf``).
    """
    nonstd_spaces = set(spaces.get_nonstandard())
    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_mask_std",
                "asl_mask_t1",
                "asl_std",
                "asl_std_ref",
                "asl_t1",
                "asl_t1_ref",
                "asl_native",
                "asl_native_ref",
                "asl_mask_native",
                "source_file",
                "itk_asl_to_t1",
                "itk_t1_to_asl",
                "template",
                "spatial_reference",
                "cbf",
                "meancbf",
                "score",
                "avgscore",
                "scrub",
                "basil",
                "pv",
                "pvwm",
                "cbf_t1",
                "meancbf_t1",
                "att_t1",
                "score_t1",
                "avgscore_t1",
                "scrub_t1",
                "basil_t1",
                "pv_t1",
                "pvwm_t1",
                "cbf_std",
                "meancbf_std",
                "score_std",
                "avgscore_std",
                "scrub_std",
                "basil_std",
                "pv_std",
                "pvwm_std",
                "att",
                "att_std",
                "qc_file",
                # Parcellated CBF outputs
                "atlas_names",
                "mean_cbf_parcellated",
                "mean_cbf_score_parcellated",
                "mean_cbf_scrub_parcellated",
                "mean_cbf_basil_parcellated",
                "mean_cbf_gm_basil_parcellated",
            ],
        ),
        name="inputnode",
    )

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name="raw_sources")
    raw_sources.inputs.bids_root = bids_root

    # fmt:off
    workflow.connect([(inputnode, raw_sources, [("source_file", "in_files")])])
    # fmt:on

    qcfile = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc="qualitycontrol",
            suffix="cbf",
            compress=False,
        ),
        name="qcfile",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    # fmt:off
    workflow.connect([
        (inputnode, qcfile, [
            ("source_file", "source_file"),
            ("qc_file", "in_file"),
        ]),
    ])
    # fmt:on

    # write transform matrix file between asl native space and T1w
    ds_t1w_to_asl_xform = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            to="scanner",
            mode="image",
            suffix="xfm",
            extension=".txt",
            dismiss_entities=("echo",),
            **{"from": "T1w"},
        ),
        name="ds_t1w_to_asl_xform",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([
        (inputnode, ds_t1w_to_asl_xform, [
            ("source_file", "source_file"),
            ("itk_asl_to_t1", "in_file"),
        ]),
    ])
    # fmt:on

    ds_asl_to_t1w_xform = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            to="T1w",
            mode="image",
            suffix="xfm",
            extension=".txt",
            dismiss_entities=("echo",),
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
            ("itk_t1_to_asl", "in_file"),
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
                base_directory=output_dir,
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

    if nonstd_spaces.intersection(("func", "run", "asl", "sbref")):
        ds_asl_native = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc="preproc",
                compress=True,
                SkullStripped=False,
                RepetitionTime=metadata.get("RepetitionTime"),
                RepetitionTimePreparation=metadata.get("RepetitionTimePreparation"),
                TaskName=metadata.get("TaskName"),
            ),
            name="ds_asl_native",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        ds_asl_native_ref = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                suffix="aslref",
                compress=True,
                dismiss_entities=("echo",),
            ),
            name="ds_asl_native_ref",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        ds_asl_mask_native = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc="brain",
                suffix="mask",
                compress=True,
                dismiss_entities=("echo",),
            ),
            name="ds_asl_mask_native",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        cbfnative = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                suffix="cbf",
                compress=True,
            ),
            name="cbfnative",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        meancbfnative = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc="mean",
                suffix="cbf",
                compress=True,
            ),
            name="meancbfnative",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        # fmt:off
        workflow.connect([
            (inputnode, ds_asl_native, [
                ("source_file", "source_file"),
                ("asl_native", "in_file"),
            ]),
            (inputnode, ds_asl_native_ref, [
                ("source_file", "source_file"),
                ("asl_native_ref", "in_file"),
            ]),
            (inputnode, ds_asl_mask_native, [
                ("source_file", "source_file"),
                ("asl_mask_native", "in_file"),
            ]),
            (inputnode, cbfnative, [
                ("source_file", "source_file"),
                ("cbf", "in_file"),
            ]),
            (inputnode, meancbfnative, [
                ("source_file", "source_file"),
                ("meancbf", "in_file"),
            ]),
        ])
        # fmt:on

        if scorescrub:
            scorenative = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="score",
                    suffix="cbf",
                    compress=True,
                ),
                name="scorenative",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            meanscorenative = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="meanScore",
                    suffix="cbf",
                    compress=True,
                ),
                name="meanscorenative",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            scrubnative = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="scrub",
                    suffix="cbf",
                    compress=True,
                ),
                name="scrubnative",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )

            # fmt:off
            workflow.connect([
                (inputnode, scorenative, [
                    ("source_file", "source_file"),
                    ("score", "in_file"),
                ]),
                (inputnode, meanscorenative, [
                    ("source_file", "source_file"),
                    ("avgscore", "in_file"),
                ]),
                (inputnode, scrubnative, [
                    ("source_file", "source_file"),
                    ("scrub", "in_file"),
                ]),
            ])
            # fmt:on

        if basil:
            basilnative = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="basil",
                    suffix="cbf",
                    compress=True,
                ),
                name="basilnative",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            pvnative = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="pvGM",
                    suffix="cbf",
                    compress=True,
                ),
                name="pvcnative",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            pvnativewm = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="pvWM",
                    suffix="cbf",
                    compress=True,
                ),
                name="pvcnativewm",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            attnative = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="bat",
                    suffix="cbf",
                    compress=True,
                ),
                name="attcnative",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )

            # fmt:off
            workflow.connect([
                (inputnode, basilnative, [
                    ("source_file", "source_file"),
                    ("basil", "in_file"),
                ]),
                (inputnode, pvnative, [
                    ("source_file", "source_file"),
                    ("pv", "in_file"),
                ]),
                (inputnode, pvnativewm, [
                    ("source_file", "source_file"),
                    ("pvwm", "in_file"),
                ]),
                (inputnode, attnative, [
                    ("source_file", "source_file"),
                    ("att", "in_file"),
                ]),
                (raw_sources, ds_asl_mask_native, [("out", "RawSources")]),
            ])
            # fmt:on

    # Resample to T1w space
    if nonstd_spaces.intersection(("T1w", "anat")):
        ds_asl_t1 = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                space="T1w",
                desc="preproc",
                compress=True,
                SkullStripped=False,
                RepetitionTime=metadata.get("RepetitionTime"),
                RepetitionTimePreparation=metadata.get("RepetitionTimePreparation"),
                TaskName=metadata.get("TaskName"),
                dismiss_entities=("echo",),
            ),
            name="ds_asl_t1",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        ds_asl_t1_ref = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                space="T1w",
                suffix="aslref",
                compress=True,
                dismiss_entities=("echo",),
            ),
            name="ds_asl_t1_ref",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        ds_asl_mask_t1 = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                space="T1w",
                desc="brain",
                suffix="mask",
                compress=True,
                dismiss_entities=("echo",),
            ),
            name="ds_asl_mask_t1",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        cbfnativet1 = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                suffix="cbf",
                space="T1w",
                compress=True,
            ),
            name="cbfnativet1",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        meancbfnativet1 = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc="mean",
                suffix="cbf",
                space="T1w",
                compress=True,
            ),
            name="meancbfnativet1",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        # fmt:off
        workflow.connect([
            (inputnode, ds_asl_t1, [
                ("source_file", "source_file"),
                ("asl_t1", "in_file"),
            ]),
            (inputnode, ds_asl_t1_ref, [
                ("source_file", "source_file"),
                ("asl_t1_ref", "in_file"),
            ]),
            (inputnode, ds_asl_mask_t1, [
                ("source_file", "source_file"),
                ("asl_mask_t1", "in_file"),
            ]),
            (inputnode, cbfnativet1, [
                ("source_file", "source_file"),
                ("cbf_t1", "in_file"),
            ]),
            (inputnode, meancbfnativet1, [
                ("source_file", "source_file"),
                ("meancbf_t1", "in_file"),
            ]),
        ])
        # fmt:on

        if scorescrub:
            scorenativet1 = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="score",
                    suffix="cbf",
                    space="T1w",
                    compress=True,
                ),
                name="scorenativet1",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            meanscorenativet1 = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    suffix="cbf",
                    desc="meanScore",
                    space="T1w",
                    compress=True,
                ),
                name="meanscorenativet1",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            scrubnativet1 = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="scrub",
                    suffix="cbf",
                    space="T1w",
                    compress=True,
                ),
                name="scrubnativet1",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )

            # fmt:off
            workflow.connect([
                (inputnode, scorenativet1, [
                    ("source_file", "source_file"),
                    ("score_t1", "in_file"),
                ]),
                (inputnode, meanscorenativet1, [
                    ("source_file", "source_file"),
                    ("avgscore_t1", "in_file"),
                ]),
                (inputnode, scrubnativet1, [
                    ("source_file", "source_file"),
                    ("scrub_t1", "in_file"),
                ]),
            ])
            # fmt:on

        if basil:
            basilnativet1 = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="basil",
                    suffix="cbf",
                    space="T1w",
                    compress=True,
                ),
                name="basilnativet1",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            pvnativet1 = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="pvGM",
                    suffix="cbf",
                    space="T1w",
                    compress=True,
                ),
                name="pvcnativet1",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            pvnativet1wm = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="pvWM",
                    suffix="cbf",
                    space="T1w",
                    compress=True,
                ),
                name="pvcnativet1wm",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            attnativet1 = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="bat",
                    suffix="cbf",
                    space="T1w",
                    compress=True,
                ),
                name="attnativet1",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )

            # fmt:off
            workflow.connect([
                (inputnode, basilnativet1, [
                    ("source_file", "source_file"),
                    ("basil_t1", "in_file"),
                ]),
                (inputnode, pvnativet1, [
                    ("source_file", "source_file"),
                    ("pv_t1", "in_file"),
                ]),
                (inputnode, pvnativet1wm, [
                    ("source_file", "source_file"),
                    ("pvwm_t1", "in_file"),
                ]),
                (inputnode, attnativet1, [
                    ("source_file", "source_file"),
                    ("att_t1", "in_file"),
                ]),
            ])
            # fmt:on

        # fmt:off
        workflow.connect([(raw_sources, ds_asl_mask_t1, [("out", "RawSources")])])
        # fmt:on

    if getattr(spaces, "_cached") is None:
        return workflow

    # Store resamplings in standard spaces when listed in --output-spaces
    if spaces.cached.references:
        from aslprep.niworkflows.interfaces.space import SpaceDataSource

        spacesource = pe.Node(SpaceDataSource(), name="spacesource", run_without_submitting=True)
        spacesource.iterables = (
            "in_tuple",
            [(s.fullname, s.spec) for s in spaces.cached.get_standard(dim=(3,))],
        )

        out_names = [
            "template",
            "asl_std",
            "asl_std_ref",
            "asl_mask_std",
            "cbf_std",
            "meancbf_std",
        ]
        if scorescrub:
            out_names += ["score_std", "avgscore_std", "scrub_std"]

        if basil:
            out_names += ["basil_std", "pv_std", "pvwm_std", "att_std"]

        select_std = pe.Node(
            KeySelect(fields=out_names),
            name="select_std",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        ds_asl_std = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc="preproc",
                compress=True,
                SkullStripped=False,
                RepetitionTime=metadata.get("RepetitionTime"),
                RepetitionTimePreparation=metadata.get("RepetitionTimePreparation"),
                TaskName=metadata.get("TaskName"),
                dismiss_entities=("echo",),
            ),
            name="ds_asl_std",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        ds_asl_std_ref = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                suffix="aslref",
                compress=True,
                dismiss_entities=("echo",),
            ),
            name="ds_asl_std_ref",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        ds_asl_mask_std = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc="brain",
                suffix="mask",
                compress=True,
                dismiss_entities=("echo",),
            ),
            name="ds_asl_mask_std",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        cbfstd = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                suffix="cbf",
                compress=True,
            ),
            name="cbfstd",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        meancbfstd = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc="mean",
                suffix="cbf",
                compress=True,
            ),
            name="meancbfstd",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        # fmt:off
        workflow.connect([
            (inputnode, ds_asl_std, [("source_file", "source_file")]),
            (inputnode, ds_asl_std_ref, [("source_file", "source_file")]),
            (inputnode, ds_asl_mask_std, [("source_file", "source_file")]),
            (inputnode, cbfstd, [("source_file", "source_file")]),
            (inputnode, meancbfstd, [("source_file", "source_file")]),
            (inputnode, select_std, [
                ("asl_std", "asl_std"),
                ("asl_std_ref", "asl_std_ref"),
                ("asl_mask_std", "asl_mask_std"),
                ("cbf_std", "cbf_std"),
                ("meancbf_std", "meancbf_std"),
                ("template", "template"),
                ("spatial_reference", "keys"),
            ]),
            (spacesource, select_std, [("uid", "key")]),
            (select_std, ds_asl_std, [("asl_std", "in_file")]),
            (spacesource, ds_asl_std, [
                ("space", "space"),
                ("cohort", "cohort"),
                ("resolution", "resolution"),
                ("density", "density"),
            ]),
            (select_std, ds_asl_std_ref, [("asl_std_ref", "in_file")]),
            (spacesource, ds_asl_std_ref, [
                ("space", "space"),
                ("cohort", "cohort"),
                ("resolution", "resolution"),
                ("density", "density"),
            ]),
            (select_std, ds_asl_mask_std, [("asl_mask_std", "in_file")]),
            (spacesource, ds_asl_mask_std, [
                ("space", "space"),
                ("cohort", "cohort"),
                ("resolution", "resolution"),
                ("density", "density"),
            ]),
            (select_std, cbfstd, [("cbf_std", "in_file")]),
            (spacesource, cbfstd, [
                ("space", "space"),
                ("cohort", "cohort"),
                ("resolution", "resolution"),
                ("density", "density"),
            ]),
            (select_std, meancbfstd, [("meancbf_std", "in_file")]),
            (spacesource, meancbfstd, [
                ("space", "space"),
                ("cohort", "cohort"),
                ("resolution", "resolution"),
                ("density", "density"),
            ]),
            (raw_sources, ds_asl_mask_std, [("out", "RawSources")]),
        ])
        # fmt:on

        if scorescrub:
            scorestd = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="score",
                    suffix="cbf",
                    compress=True,
                ),
                name="scorestd",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            meanscorestd = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="meanScore",
                    suffix="cbf",
                    compress=True,
                ),
                name="meanscorestd",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            scrubstd = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="scrub",
                    suffix="cbf",
                    compress=True,
                ),
                name="scrubstd",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )

            # fmt:off
            workflow.connect([
                (inputnode, scorestd, [("source_file", "source_file")]),
                (inputnode, meanscorestd, [("source_file", "source_file")]),
                (inputnode, scrubstd, [("source_file", "source_file")]),
                (inputnode, select_std, [
                    ("score_std", "score_std"),
                    ("avgscore_std", "avgscore_std"),
                    ("scrub_std", "scrub_std"),
                ]),
                (select_std, scorestd, [("score_std", "in_file")]),
                (spacesource, scorestd, [
                    ("space", "space"),
                    ("cohort", "cohort"),
                    ("resolution", "resolution"),
                    ("density", "density"),
                ]),
                (select_std, meanscorestd, [("avgscore_std", "in_file")]),
                (spacesource, meanscorestd, [
                    ("space", "space"),
                    ("cohort", "cohort"),
                    ("resolution", "resolution"),
                    ("density", "density"),
                ]),
                (select_std, scrubstd, [("scrub_std", "in_file")]),
                (spacesource, scrubstd, [
                    ("space", "space"),
                    ("cohort", "cohort"),
                    ("resolution", "resolution"),
                    ("density", "density"),
                ]),
            ])
            # fmt:on

        if basil:
            basilstd = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="basil",
                    suffix="cbf",
                    compress=True,
                ),
                name="basilstd",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            pvstd = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="pvGM",
                    suffix="cbf",
                    compress=True,
                ),
                name="pvcstd",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            pvstdwm = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="pvWM",
                    suffix="cbf",
                    compress=True,
                ),
                name="pvcstdwm",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            attstd = pe.Node(
                DerivativesDataSink(
                    base_directory=output_dir,
                    desc="bat",
                    suffix="cbf",
                    compress=True,
                ),
                name="attstd",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )

            # fmt:off
            workflow.connect([
                (inputnode, basilstd, [("source_file", "source_file")]),
                (inputnode, pvstd, [("source_file", "source_file")]),
                (inputnode, pvstdwm, [("source_file", "source_file")]),
                (inputnode, attstd, [("source_file", "source_file")]),
                (inputnode, select_std, [
                    ("basil_std", "basil_std"),
                    ("pv_std", "pv_std"),
                    ("pvwm_std", "pvwm_std"),
                    ("att_std", "att_std"),
                ]),
                (select_std, basilstd, [("basil_std", "in_file")]),
                (spacesource, basilstd, [
                    ("space", "space"),
                    ("cohort", "cohort"),
                    ("resolution", "resolution"),
                    ("density", "density"),
                ]),
                (select_std, pvstd, [("pv_std", "in_file")]),
                (spacesource, pvstd, [
                    ("space", "space"),
                    ("cohort", "cohort"),
                    ("resolution", "resolution"),
                    ("density", "density"),
                ]),
                (select_std, pvstdwm, [("pvwm_std", "in_file")]),
                (spacesource, pvstdwm, [
                    ("space", "space"),
                    ("cohort", "cohort"),
                    ("resolution", "resolution"),
                    ("density", "density"),
                ]),
                (select_std, attstd, [("att_std", "in_file")]),
                (spacesource, attstd, [
                    ("space", "space"),
                    ("cohort", "cohort"),
                    ("resolution", "resolution"),
                    ("density", "density"),
                ]),
            ])
            # fmt:on

    return workflow
