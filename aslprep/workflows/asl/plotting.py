# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for plotting ASLPrep derivatives."""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.confounds import GatherCBFConfounds
from aslprep.interfaces.plotting import CBFByTissueTypePlot, CBFSummaryPlot
from aslprep.interfaces.reports import CBFSummary
from aslprep.utils.misc import _select_last_in_list
from aslprep.workflows.asl.confounds import init_carpetplot_wf


def init_plot_cbf_wf(
    metadata,
    plot_timeseries=True,
    scorescrub=False,
    basil=False,
    name="plot_cbf_wf",
):
    """Plot CBF results.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.plotting import init_plot_cbf_wf

            wf = init_plot_cbf_wf(
                metadata={
                    "RepetitionTime": 4,
                    "RepetitionTimePreparation": 4,
                },
            )
    """
    from niworkflows.interfaces.images import SignalExtraction

    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "aslref",
                "asl_mask",
                "t1w_dseg",
                "aslref2anat_xfm",
                "std2anat_xfm",
                "confounds_file",
                "qc_file",
                # If plot_timeseries is True
                "crown_mask",
                "acompcor_masks",
                # CBF outputs
                "mean_cbf",
                # Single-delay outputs
                "cbf_ts",  # only for non-GE
                # Multi-delay outputs
                "att",
                # SCORE/SCRUB outputs
                "cbf_ts_score",  # unused
                "mean_cbf_score",
                "mean_cbf_scrub",
                "score_outlier_index",
                # BASIL outputs
                "mean_cbf_basil",
                "mean_cbf_gm_basil",
                "mean_cbf_wm_basil",  # unused
                "att_basil",  # unused
            ],
        ),
        name="inputnode",
    )

    summary = pe.Node(
        CBFSummary(),
        name="summary",
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )
    workflow.connect([
        (inputnode, summary, [
            ("confounds_file", "confounds_file"),
            ("qc_file", "qc_file"),
        ])
    ])  # fmt:skip

    # Warp dseg file from T1w space to ASL reference space
    warp_t1w_dseg_to_aslref = pe.Node(
        ApplyTransforms(
            float=True,
            dimension=3,
            default_value=0,
            interpolation="GenericLabel",
            invert_transform_flags=[True],
            args="-v",
        ),
        name="warp_t1w_dseg_to_aslref",
    )

    # fmt:off
    workflow.connect([
        (inputnode, warp_t1w_dseg_to_aslref, [
            ("asl_mask", "reference_image"),
            ("t1w_dseg", "input_image"),
            ("aslref2anat_xfm", "transforms"),
        ]),
    ])
    # fmt:on

    if plot_timeseries:
        # Global and segment regressors
        signals_class_labels = [
            "global_signal",
            "csf",
            "white_matter",
            "csf_wm",
        ]
        merge_rois = pe.Node(
            niu.Merge(2, ravel_inputs=True),
            name="merge_rois",
            run_without_submitting=True,
        )
        signals = pe.Node(
            SignalExtraction(class_labels=signals_class_labels),
            name="signals",
            mem_gb=2,
        )
        # fmt:off
        workflow.connect([
            (inputnode, merge_rois, [
                ("asl_mask", "in1"),
                ("acompcor_masks", "in2"),
            ]),
            (inputnode, signals, [("cbf_ts", "in_file")]),
            (merge_rois, signals, [("out", "label_files")]),
        ])
        # fmt:on

        # Time series are only available for non-GE data.
        # Create confounds file with SCORE index
        create_cbf_confounds = pe.Node(
            GatherCBFConfounds(),
            name="create_cbf_confounds",
        )
        # fmt:off
        workflow.connect([
            (inputnode, create_cbf_confounds, [("score_outlier_index", "score")]),
            (signals, create_cbf_confounds, [("out_file", "signals")]),
        ])
        # fmt:on

        carpetplot_wf = init_carpetplot_wf(
            mem_gb=2,
            confounds_list=[
                ("global_signal", None, "GS"),
                ("csf", None, "GSCSF"),
                ("white_matter", None, "GSWM"),
            ]
            + ([("score_outlier_index", None, "SCORE Index")] if scorescrub else []),
            metadata=metadata,
            cifti_output=False,
            suffix="cbf",
            name="cbf_carpetplot_wf",
        )
        carpetplot_wf.inputs.inputnode.dummy_scans = 0

        # fmt:off
        workflow.connect([
            (inputnode, carpetplot_wf, [
                ("std2anat_xfm", "inputnode.std2anat_xfm"),
                ("cbf_ts", "inputnode.asl"),
                ("asl_mask", "inputnode.asl_mask"),
                ("aslref2anat_xfm", "inputnode.aslref2anat_xfm"),
                ("crown_mask", "inputnode.crown_mask"),
                (("acompcor_masks", _select_last_in_list), "inputnode.acompcor_mask"),
            ]),
            (create_cbf_confounds, carpetplot_wf, [
                ("confounds_file", "inputnode.confounds_file"),
            ]),
        ])
        # fmt:on

    cbf_summary = pe.Node(CBFSummaryPlot(label="cbf", vmax=100), name="cbf_summary", mem_gb=1)

    # fmt:off
    workflow.connect([
        (inputnode, cbf_summary, [
            ("mean_cbf", "cbf"),
            ("aslref", "ref_vol"),
        ]),
    ])
    # fmt:on

    ds_report_cbf_summary = pe.Node(
        DerivativesDataSink(desc="cbfplot", datatype="figures", keep_dtype=True),
        name="ds_report_cbf_summary",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    workflow.connect([(cbf_summary, ds_report_cbf_summary, [("out_file", "in_file")])])

    cbf_by_tt_plot = pe.Node(
        CBFByTissueTypePlot(),
        name="cbf_by_tt_plot",
    )

    # fmt:off
    workflow.connect([
        (inputnode, cbf_by_tt_plot, [("mean_cbf", "cbf")]),
        (warp_t1w_dseg_to_aslref, cbf_by_tt_plot, [("output_image", "seg_file")]),
    ])
    # fmt:on

    ds_report_cbf_to_tt_plot = pe.Node(
        DerivativesDataSink(desc="cbfByTissueType", datatype="figures", keep_dtype=True),
        name="ds_report_cbf_to_tt_plot",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    workflow.connect([(cbf_by_tt_plot, ds_report_cbf_to_tt_plot, [("out_file", "in_file")])])

    if scorescrub:
        score_summary = pe.Node(
            CBFSummaryPlot(label="score", vmax=100),
            name="score_summary",
            mem_gb=1,
        )

        # fmt:off
        workflow.connect([
            (inputnode, score_summary, [
                ("mean_cbf_score", "cbf"),
                ("aslref", "ref_vol"),
            ]),
        ])
        # fmt:on

        ds_report_score_summary = pe.Node(
            DerivativesDataSink(desc="scoreplot", datatype="figures", keep_dtype=True),
            name="ds_report_score_summary",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        workflow.connect([(score_summary, ds_report_score_summary, [("out_file", "in_file")])])

        score_by_tt_plot = pe.Node(
            CBFByTissueTypePlot(),
            name="score_by_tt_plot",
        )

        # fmt:off
        workflow.connect([
            (inputnode, score_by_tt_plot, [("mean_cbf_score", "cbf")]),
            (warp_t1w_dseg_to_aslref, score_by_tt_plot, [("output_image", "seg_file")]),
        ])
        # fmt:on

        ds_report_score_by_tt_plot = pe.Node(
            DerivativesDataSink(desc="scoreByTissueType", datatype="figures", keep_dtype=True),
            name="ds_report_score_by_tt_plot",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        # fmt:off
        workflow.connect([
            (score_by_tt_plot, ds_report_score_by_tt_plot, [("out_file", "in_file")]),
        ])
        # fmt:on

        scrub_summary = pe.Node(
            CBFSummaryPlot(label="scrub", vmax=100),
            name="scrub_summary",
            mem_gb=1,
        )

        # fmt:off
        workflow.connect([
            (inputnode, scrub_summary, [
                ("mean_cbf_scrub", "cbf"),
                ("aslref", "ref_vol"),
            ]),
        ])
        # fmt:on

        ds_report_scrub_summary = pe.Node(
            DerivativesDataSink(desc="scrubplot", datatype="figures", keep_dtype=True),
            name="ds_report_scrub_summary",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        workflow.connect([(scrub_summary, ds_report_scrub_summary, [("out_file", "in_file")])])

        scrub_by_tt_plot = pe.Node(
            CBFByTissueTypePlot(),
            name="scrub_by_tt_plot",
        )

        # fmt:off
        workflow.connect([
            (inputnode, scrub_by_tt_plot, [("mean_cbf_scrub", "cbf")]),
            (warp_t1w_dseg_to_aslref, scrub_by_tt_plot, [("output_image", "seg_file")]),
        ])
        # fmt:on

        ds_report_scrub_by_tt_plot = pe.Node(
            DerivativesDataSink(desc="scrubByTissueType", datatype="figures", keep_dtype=True),
            name="ds_report_scrub_by_tt_plot",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        # fmt:off
        workflow.connect([
            (scrub_by_tt_plot, ds_report_scrub_by_tt_plot, [("out_file", "in_file")]),
        ])
        # fmt:on

    if basil:
        basil_summary = pe.Node(
            CBFSummaryPlot(label="basil", vmax=100),
            name="basil_summary",
            mem_gb=1,
        )

        # fmt:off
        workflow.connect([
            (inputnode, basil_summary, [
                ("mean_cbf_basil", "cbf"),
                ("aslref", "ref_vol"),
            ]),
        ])
        # fmt:on

        ds_report_basil_summary = pe.Node(
            DerivativesDataSink(desc="basilplot", datatype="figures", keep_dtype=True),
            name="ds_report_basil_summary",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        workflow.connect([(basil_summary, ds_report_basil_summary, [("out_file", "in_file")])])

        basil_by_tt_plot = pe.Node(
            CBFByTissueTypePlot(),
            name="basil_by_tt_plot",
        )

        # fmt:off
        workflow.connect([
            (inputnode, basil_by_tt_plot, [("mean_cbf_basil", "cbf")]),
            (warp_t1w_dseg_to_aslref, basil_by_tt_plot, [("output_image", "seg_file")]),
        ])
        # fmt:on

        ds_report_basil_by_tt_plot = pe.Node(
            DerivativesDataSink(desc="basilByTissueType", datatype="figures", keep_dtype=True),
            name="ds_report_basil_by_tt_plot",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        # fmt:off
        workflow.connect([
            (basil_by_tt_plot, ds_report_basil_by_tt_plot, [("out_file", "in_file")]),
        ])
        # fmt:on

        pvc_summary = pe.Node(
            CBFSummaryPlot(label="pvc", vmax=120),
            name="pvc_summary",
            mem_gb=1,
        )

        # fmt:off
        workflow.connect([
            (inputnode, pvc_summary, [
                ("mean_cbf_gm_basil", "cbf"),
                ("aslref", "ref_vol"),
            ]),
        ])
        # fmt:on

        ds_report_pvc_summary = pe.Node(
            DerivativesDataSink(desc="pvcplot", datatype="figures", keep_dtype=True),
            name="ds_report_pvc_summary",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        workflow.connect([(pvc_summary, ds_report_pvc_summary, [("out_file", "in_file")])])

        pvc_by_tt_plot = pe.Node(
            CBFByTissueTypePlot(),
            name="pvc_by_tt_plot",
        )

        # fmt:off
        workflow.connect([
            (inputnode, pvc_by_tt_plot, [("mean_cbf_gm_basil", "cbf")]),
            (warp_t1w_dseg_to_aslref, pvc_by_tt_plot, [("output_image", "seg_file")]),
        ])
        # fmt:on

        ds_report_pvc_by_tt_plot = pe.Node(
            DerivativesDataSink(desc="pvcByTissueType", datatype="figures", keep_dtype=True),
            name="ds_report_pvc_by_tt_plot",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        workflow.connect([(pvc_by_tt_plot, ds_report_pvc_by_tt_plot, [("out_file", "in_file")])])

    return workflow
