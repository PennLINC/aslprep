# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for plotting ASLPrep derivatives."""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from templateflow.api import get as get_template

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.plotting import CBFSummary, CBFtsSummary
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow


def init_cbfplot_wf(
    metadata,
    scorescrub=False,
    basil=False,
    name="cbf_plot",
):
    """Plot CBF results.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.plotting import init_cbfplot_wf

            wf = init_cbfplot_wf(
                metadata={"RepetitionTimePreparation": 4},
            )
    """
    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "aslref",
                "asl_mask",
                "anat_to_aslref_xfm",
                "template_to_anat_xfm",
                "confounds_file",
                # CBF outputs
                "cbf_ts",
                "mean_cbf",
                # SCORE/SCRUB outputs
                "cbf_ts_score",
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
    mrg_xfms = pe.Node(niu.Merge(2), name="mrg_xfms")

    seg = get_template("MNI152NLin2009cAsym", resolution=1, desc="carpet", suffix="dseg")

    resample_parc = pe.Node(
        ApplyTransforms(
            float=True,
            input_image=str(seg),
            dimension=3,
            default_value=0,
            interpolation="MultiLabel",
        ),
        name="resample_parc",
    )

    cbftssummary = pe.Node(
        CBFtsSummary(tr=metadata.get("RepetitionTime", metadata["RepetitionTimePreparation"])),
        name="cbf_ts_summary",
        mem_gb=2,
    )
    cbfsummary = pe.Node(CBFSummary(label="cbf", vmax=90), name="cbf_summary", mem_gb=1)
    ds_report_cbftsplot = pe.Node(
        DerivativesDataSink(desc="cbftsplot", datatype="figures", keep_dtype=True),
        name="ds_report_cbftsplot",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    ds_report_cbfplot = pe.Node(
        DerivativesDataSink(desc="cbfplot", datatype="figures", keep_dtype=True),
        name="ds_report_cbfplot",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    # fmt:off
    workflow.connect([
        (inputnode, mrg_xfms, [
            ("template_to_anat_xfm", "in1"),
            ("anat_to_aslref_xfm", "in2"),
        ]),
        (inputnode, resample_parc, [("asl_mask", "reference_image")]),
        (mrg_xfms, resample_parc, [("out", "transforms")]),
        (resample_parc, cbftssummary, [("output_image", "seg_file")]),
        (inputnode, cbftssummary, [
            ("cbf_ts", "cbf_ts"),
            ("confounds_file", "confounds_file"),
            ("score_outlier_index", "score_outlier_index"),
        ]),
        (cbftssummary, ds_report_cbftsplot, [("out_file", "in_file")]),
        (inputnode, cbfsummary, [
            ("mean_cbf", "cbf"),
            ("aslref", "ref_vol"),
        ]),
        (cbfsummary, ds_report_cbfplot, [("out_file", "in_file")]),
    ])
    # fmt:on

    if scorescrub:
        scoresummary = pe.Node(CBFSummary(label="score", vmax=90), name="score_summary", mem_gb=1)
        scrubsummary = pe.Node(CBFSummary(label="scrub", vmax=90), name="scrub_summary", mem_gb=1)
        ds_report_scoreplot = pe.Node(
            DerivativesDataSink(desc="scoreplot", datatype="figures", keep_dtype=True),
            name="ds_report_scoreplot",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        ds_report_scrubplot = pe.Node(
            DerivativesDataSink(desc="scrubplot", datatype="figures", keep_dtype=True),
            name="ds_report_scrubplot",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        # fmt:off
        workflow.connect([
            (inputnode, scoresummary, [
                ("mean_cbf_score", "cbf"),
                ("aslref", "ref_vol"),
            ]),
            (scoresummary, ds_report_scoreplot, [("out_file", "in_file")]),
            (inputnode, scrubsummary, [
                ("mean_cbf_scrub", "cbf"),
                ("aslref", "ref_vol"),
            ]),
            (scrubsummary, ds_report_scrubplot, [("out_file", "in_file")]),
        ])
        # fmt:on

    if basil:
        basilsummary = pe.Node(CBFSummary(label="basil", vmax=100), name="basil_summary", mem_gb=1)
        pvcsummary = pe.Node(CBFSummary(label="pvc", vmax=120), name="pvc_summary", mem_gb=1)
        ds_report_basilplot = pe.Node(
            DerivativesDataSink(desc="basilplot", datatype="figures", keep_dtype=True),
            name="ds_report_basilplot",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        ds_report_pvcplot = pe.Node(
            DerivativesDataSink(desc="pvcplot", datatype="figures", keep_dtype=True),
            name="ds_report_pvcplot",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        # fmt:off
        workflow.connect([
            (inputnode, basilsummary, [
                ("mean_cbf_basil", "cbf"),
                ("aslref", "ref_vol"),
            ]),
            (basilsummary, ds_report_basilplot, [("out_file", "in_file")]),
            (inputnode, pvcsummary, [
                ("mean_cbf_gm_basil", "cbf"),
                ("aslref", "ref_vol"),
            ]),
            (pvcsummary, ds_report_pvcplot, [("out_file", "in_file")]),
        ])
        # fmt:on

    return workflow


def init_gecbfplot_wf(basil=False, name="cbf_plot"):
    """Plot CBF results for GE data.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.plotting import init_gecbfplot_wf

            wf = init_gecbfplot_wf(basil=True)
    """
    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "mean_cbf",
                "aslref",
                "mean_cbf_basil",
                "mean_cbf_gm_basil",
                "mean_cbf_wm_basil",  # unused
            ],
        ),
        name="inputnode",
    )

    cbfsummary = pe.Node(CBFSummary(label="cbf", vmax=90), name="cbf_summary", mem_gb=1)
    ds_report_cbfplot = pe.Node(
        DerivativesDataSink(desc="cbfplot", datatype="figures", keep_dtype=True),
        name="ds_report_cbfplot",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    # fmt:off
    workflow.connect([
        (inputnode, cbfsummary, [
            ("mean_cbf", "cbf"),
            ("aslref", "ref_vol"),
        ]),
        (cbfsummary, ds_report_cbfplot, [("out_file", "in_file")]),
    ])
    # fmt:on

    if basil:
        basilsummary = pe.Node(CBFSummary(label="basil", vmax=100), name="basil_summary", mem_gb=1)
        pvcsummary = pe.Node(CBFSummary(label="pvc", vmax=120), name="pvc_summary", mem_gb=1)
        ds_report_basilplot = pe.Node(
            DerivativesDataSink(desc="basilplot", datatype="figures", keep_dtype=True),
            name="ds_report_basilplot",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        ds_report_pvcplot = pe.Node(
            DerivativesDataSink(desc="pvcplot", datatype="figures", keep_dtype=True),
            name="ds_report_pvcplot",
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        # fmt:off
        workflow.connect([
            (inputnode, basilsummary, [
                ("mean_cbf_basil", "cbf"),
                ("aslref", "ref_vol"),
            ]),
            (basilsummary, ds_report_basilplot, [("out_file", "in_file")]),
            (inputnode, pvcsummary, [
                ("mean_cbf_gm_basil", "cbf"),
                ("aslref", "ref_vol"),
            ]),
            (pvcsummary, ds_report_pvcplot, [("out_file", "in_file")]),
        ])
        # fmt:on

    return workflow
