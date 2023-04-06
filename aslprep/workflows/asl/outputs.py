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
                "source_file",
                "template",
                "spatial_reference",
                "qc_file",
                # Preprocessed ASL files
                "asl_native",
                "asl_t1",
                "asl_std",
                "aslref_native",
                "aslref_t1",
                "aslref_std",
                "asl_mask_native",
                "asl_mask_t1",
                "asl_mask_std",
                # Transforms
                "itk_asl_to_t1",
                "itk_t1_to_asl",
                # Standard CBF outputs
                "cbf_native",
                "cbf_t1",
                "cbf_std",
                "meancbf_native",
                "meancbf_t1",
                "meancbf_std",
                # SCORE/SCRUB outputs
                "score_native",
                "score_t1",
                "score_std",
                "avgscore_native",
                "avgscore_t1",
                "avgscore_std",
                "scrub_native",
                "scrub_t1",
                "scrub_std",
                # BASIL outputs
                "basil_native",
                "basil_t1",
                "basil_std",
                "pv_native",
                "pv_t1",
                "pv_std",
                "pvwm_native",
                "pvwm_t1",
                "pvwm_std",
                "att_native",
                "att_t1",
                "att_std",
                # Parcellated CBF outputs
                "atlas_names",
                "mean_cbf_parcellated",
                "mean_cbf_score_parcellated",
                "mean_cbf_scrub_parcellated",
                "mean_cbf_basil_parcellated",
                "mean_cbf_gm_basil_parcellated",
                # non-GE outputs
                "confounds",
                "confounds_metadata",
            ],
        ),
        name="inputnode",
    )

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name="raw_sources")
    raw_sources.inputs.bids_root = bids_root

    if output_confounds:
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

    ds_qcfile = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
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
            ("itk_t1_to_asl", "in_file"),
        ]),
    ])
    # fmt:on

    ds_asl_to_t1w_xform = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
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
            ("itk_asl_to_t1", "in_file"),
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

    # Now prepare to write out primary imaging derivatives
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

    # Native-space derivatives
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

    # T1w-space derivatives
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

    # Standard-space derivatives
    from niworkflows.interfaces.space import SpaceDataSource

    spacesource = pe.Node(SpaceDataSource(), name="spacesource", run_without_submitting=True)
    spacesource.iterables = (
        "in_tuple",
        [(s.fullname, s.spec) for s in spaces.cached.get_standard(dim=(3,))],
    )

    select_std = pe.Node(
        KeySelect(fields=[f"{bi}_std" for bi in base_inputs] + ["template"]),
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
            (inputnode, ds_base_input_std, [("source_file", "source_file")]),
            (inputnode, select_std, [(f"{base_input}_std", f"{base_input}_std")]),
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
