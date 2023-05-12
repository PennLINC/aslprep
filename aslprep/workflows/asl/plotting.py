# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for plotting ASLPrep derivatives."""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.plotting import CBFByTissueTypePlot, CBFSummary, CBFtsSummary
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.utils.misc import get_template_str


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
                metadata={"RepetitionTimePreparation": 4},
            )
    """
    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "aslref",
                "asl_mask",
                "t1w_dseg",
                "anat_to_aslref_xfm",
                "template_to_anat_xfm",
                "confounds_file",  # only for non-GE
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
            ]
        ),
        name="inputnode",
    )

    # String together transforms from MNI152NLin2009cAsym to ASL reference
    mrg_xfms = pe.Node(niu.Merge(2), name="mrg_xfms")

    # fmt:off
    workflow.connect([
        (inputnode, mrg_xfms, [
            ("template_to_anat_xfm", "in1"),
            ("anat_to_aslref_xfm", "in2"),
        ]),
    ])
    # fmt:on

    # Warp dseg file from T1w space to ASL reference space
    warp_t1w_dseg_to_aslref = pe.Node(
        ApplyTransforms(
            float=True,
            dimension=3,
            default_value=0,
            interpolation="MultiLabel",
        ),
        name="warp_t1w_dseg_to_aslref",
    )

    # fmt:off
    workflow.connect([
        (inputnode, warp_t1w_dseg_to_aslref, [
            ("asl_mask", "reference_image"),
            ("t1w_dseg", "input_image"),
            ("anat_to_aslref_xfm", "transforms"),
        ]),
    ])
    # fmt:on

    grab_carpet_dseg = pe.Node(
        niu.Function(
            input_names=["template", "kwargs"],
            output_names=["carpet_dseg"],
            function=get_template_str,
        ),
        name="grab_carpet_dseg",
    )
    grab_carpet_dseg.inputs.template = "MNI152NLin2009cAsym"
    grab_carpet_dseg.inputs.kwargs = {
        "resolution": 1,
        "desc": "carpet",
        "suffix": "dseg",
    }

    warp_carpet_dseg_to_aslref = pe.Node(
        ApplyTransforms(
            float=True,
            dimension=3,
            default_value=0,
            interpolation="MultiLabel",
        ),
        name="warp_carpet_dseg_to_aslref",
    )

    # fmt:off
    workflow.connect([
        (inputnode, warp_carpet_dseg_to_aslref, [("asl_mask", "reference_image")]),
        (mrg_xfms, warp_carpet_dseg_to_aslref, [("out", "transforms")]),
        (grab_carpet_dseg, warp_carpet_dseg_to_aslref, [("carpet_dseg", "input_image")]),
    ])
    # fmt:on

    if plot_timeseries:
        # Time series are only available for non-GE data.
        cbf_ts_summary = pe.Node(
            CBFtsSummary(tr=metadata.get("RepetitionTime", metadata["RepetitionTimePreparation"])),
            name="cbf_ts_summary",
            mem_gb=2,
        )

        # fmt:off
        workflow.connect([
            (inputnode, cbf_ts_summary, [
                ("cbf_ts", "cbf_ts"),
                ("confounds_file", "confounds_file"),
                ("score_outlier_index", "score_outlier_index"),
            ]),
            (warp_carpet_dseg_to_aslref, cbf_ts_summary, [("output_image", "seg_file")]),
        ])
        # fmt:on

        ds_report_cbf_ts_summary = pe.Node(
            DerivativesDataSink(desc="cbftsplot", datatype="figures", keep_dtype=True),
            name="ds_report_cbf_ts_summary",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        workflow.connect([(cbf_ts_summary, ds_report_cbf_ts_summary, [("out_file", "in_file")])])

    cbf_summary = pe.Node(CBFSummary(label="cbf", vmax=100), name="cbf_summary", mem_gb=1)

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
            CBFSummary(label="score", vmax=100),
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
            CBFSummary(label="scrub", vmax=100),
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
            CBFSummary(label="basil", vmax=100),
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
            CBFSummary(label="pvc", vmax=120),
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
