# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for writing out derivative files."""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.utility import KeySelect
from smriprep.workflows.outputs import _bids_relative

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.utils.spaces import SpatialReferences


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
                "asl_cifti",
                "aslref_native",
                "aslref_t1",
                "aslref_std",
                "asl_mask_native",
                "asl_mask_t1",
                "asl_mask_std",
                # Transforms
                "hmc_xforms",
                "aslref_to_anat_xfm",
                "anat_to_aslref_xfm",
                # Standard CBF outputs
                "cbf_ts_native",
                "cbf_ts_t1",
                "cbf_ts_std",
                "cbf_ts_cifti",
                "mean_cbf_native",
                "mean_cbf_t1",
                "mean_cbf_std",
                "mean_cbf_cifti",
                "att_native",
                "att_t1",
                "att_std",
                "att_cifti",
                # SCORE/SCRUB outputs
                "cbf_ts_score_native",
                "cbf_ts_score_t1",
                "cbf_ts_score_std",
                "cbf_ts_score_cifti",
                "mean_cbf_score_native",
                "mean_cbf_score_t1",
                "mean_cbf_score_std",
                "mean_cbf_score_cifti",
                "mean_cbf_scrub_native",
                "mean_cbf_scrub_t1",
                "mean_cbf_scrub_std",
                "mean_cbf_scrub_cifti",
                # BASIL outputs
                "mean_cbf_basil_native",
                "mean_cbf_basil_t1",
                "mean_cbf_basil_std",
                "mean_cbf_basil_cifti",
                "mean_cbf_gm_basil_native",
                "mean_cbf_gm_basil_t1",
                "mean_cbf_gm_basil_std",
                "mean_cbf_gm_basil_cifti",
                "mean_cbf_wm_basil_native",
                "mean_cbf_wm_basil_t1",
                "mean_cbf_wm_basil_std",
                "mean_cbf_wm_basil_cifti",
                "att_basil_native",
                "att_basil_t1",
                "att_basil_std",
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

    # write transform matrix file between asl native space and T1w
    ds_t1w_to_asl_xform = pe.Node(
        DerivativesDataSink(
            base_directory=config.execution.aslprep_dir,
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
            ("anat_to_aslref_xfm", "in_file"),
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
            ("aslref_to_anat_xfm", "in_file"),
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

    # Now prepare to write out primary imaging derivatives
    asl_metadata = {
        "SkullStripped": False,
        "RepetitionTime": metadata.get("RepetitionTime"),
        "RepetitionTimePreparation": metadata.get("RepetitionTimePreparation"),
    }
    cbf_metadata = {
        "Units": "mL/100 g/min",
    }

    BASE_INPUT_FIELDS = {
        "asl": {
            "desc": "preproc",
            "suffix": "asl",
            **asl_metadata,
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
        # CBF outputs
        "cbf_ts": {
            "desc": "timeseries",
            "suffix": "cbf",
            **cbf_metadata,
        },
        "mean_cbf": {
            "suffix": "cbf",
            **cbf_metadata,
        },
        "att": {
            "suffix": "att",
            "Units": "s",
        },
        # SCORE/SCRUB outputs
        "cbf_ts_score": {
            "desc": "scoreTimeseries",
            "suffix": "cbf",
            **cbf_metadata,
        },
        "mean_cbf_score": {
            "desc": "score",
            "suffix": "cbf",
            **cbf_metadata,
        },
        "mean_cbf_scrub": {
            "desc": "scrub",
            "suffix": "cbf",
            **cbf_metadata,
        },
        # BASIL outputs
        "mean_cbf_basil": {
            "desc": "basil",
            "suffix": "cbf",
            **cbf_metadata,
        },
        "mean_cbf_gm_basil": {
            "desc": "pvGM",
            "suffix": "cbf",
            **cbf_metadata,
        },
        "mean_cbf_wm_basil": {
            "desc": "pvWM",
            "suffix": "cbf",
            **cbf_metadata,
        },
        "att_basil": {
            "desc": "basil",
            "suffix": "att",
            "Units": "s",
        },
    }

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
            ds_bold_aseg_t1 = pe.Node(
                DerivativesDataSink(
                    base_directory=config.execution.aslprep_dir,
                    space="T1w",
                    desc="aseg",
                    suffix="dseg",
                    compress=True,
                    dismiss_entities=("echo",),
                ),
                name="ds_bold_aseg_t1",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            ds_bold_aparc_t1 = pe.Node(
                DerivativesDataSink(
                    base_directory=config.execution.aslprep_dir,
                    space="T1w",
                    desc="aparcaseg",
                    suffix="dseg",
                    compress=True,
                    dismiss_entities=("echo",),
                ),
                name="ds_bold_aparc_t1",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            # fmt:off
            workflow.connect([
                (inputnode, ds_bold_aseg_t1, [
                    ("source_file", "source_file"),
                    ("bold_aseg_t1", "in_file"),
                ]),
                (inputnode, ds_bold_aparc_t1, [
                    ("source_file", "source_file"),
                    ("bold_aparc_t1", "in_file"),
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
                KeySelect(fields=["bold_aseg_std", "bold_aparc_std", "template"]),
                name="select_fs_std",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            ds_bold_aseg_std = pe.Node(
                DerivativesDataSink(
                    base_directory=config.execution.aslprep_dir,
                    desc="aseg",
                    suffix="dseg",
                    compress=True,
                    dismiss_entities=("echo",),
                ),
                name="ds_bold_aseg_std",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            ds_bold_aparc_std = pe.Node(
                DerivativesDataSink(
                    base_directory=config.execution.aslprep_dir,
                    desc="aparcaseg",
                    suffix="dseg",
                    compress=True,
                    dismiss_entities=("echo",),
                ),
                name="ds_bold_aparc_std",
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            # fmt:off
            workflow.connect([
                (spacesource, select_fs_std, [("uid", "key")]),
                (inputnode, select_fs_std, [
                    ("bold_aseg_std", "bold_aseg_std"),
                    ("bold_aparc_std", "bold_aparc_std"),
                    ("template", "template"),
                    ("spatial_reference", "keys"),
                ]),
                (select_fs_std, ds_bold_aseg_std, [("bold_aseg_std", "in_file")]),
                (spacesource, ds_bold_aseg_std, [
                    ("space", "space"),
                    ("cohort", "cohort"),
                    ("resolution", "resolution"),
                    ("density", "density"),
                ]),
                (select_fs_std, ds_bold_aparc_std, [("bold_aparc_std", "in_file")]),
                (spacesource, ds_bold_aparc_std, [
                    ("space", "space"),
                    ("cohort", "cohort"),
                    ("resolution", "resolution"),
                    ("density", "density"),
                ]),
                (inputnode, ds_bold_aseg_std, [("source_file", "source_file")]),
                (inputnode, ds_bold_aparc_std, [("source_file", "source_file")])
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

        ds_bold_surfs = pe.MapNode(
            DerivativesDataSink(
                base_directory=config.execution.aslprep_dir,
                extension=".func.gii",
                TaskName=metadata.get("TaskName"),
            ),
            iterfield=["in_file", "hemi"],
            name="ds_bold_surfs",
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
            (inputnode, ds_bold_surfs, [("source_file", "source_file")]),
            (select_fs_surf, ds_bold_surfs, [
                ("surfaces", "in_file"),
                ("key", "space"),
            ]),
            (name_surfs, ds_bold_surfs, [("hemi", "hemi")]),
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
        ds_bold_cifti = pe.Node(
            DerivativesDataSink(
                base_directory=config.execution.aslprep_dir,
                suffix="bold",
                compress=False,
                TaskName=metadata.get("TaskName"),
                space="fsLR",
            ),
            name="ds_bold_cifti",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        # fmt:off
        workflow.connect([
            (inputnode, ds_bold_cifti, [
                ("bold_cifti", "in_file"),
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
