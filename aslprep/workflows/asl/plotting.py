# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for plotting ASLPrep derivatives."""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from templateflow.api import get as get_template

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.interfaces.plotting import CBFSummary, CBFtsSummary
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.interfaces.fixes import (
    FixHeaderApplyTransforms as ApplyTransforms,
)


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
                metadata={},
            )
    """
    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "cbf",
                "cbf_ts",
                "score_ts",
                "score",
                "scrub",
                "asl_ref",
                "basil",
                "pvc",
                "asl_mask",
                "t1_asl_xform",
                "std2anat_xfm",
                "confounds_file",
                "scoreindex",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "cbf_carpetplot",
                "score_carpetplot",
                "cbf_summary_plot",
                "cbf_summary_plot",
                "score_summary_plot",
                "scrub_summary_plot",
                "basil_summary_plot",
                "pvc_summary_plot",
            ]
        ),
        name="outputnode",
    )
    mrg_xfms = pe.Node(niu.Merge(2), name="mrg_xfms")

    seg = get_template("MNI152NLin2009cAsym", resolution=1, desc="carpet", suffix="dseg")
    print(seg)
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
        CBFtsSummary(tr=metadata["RepetitionTimePreparation"]),
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
        (inputnode, mrg_xfms, [("t1_asl_xform", "in2"), ("std2anat_xfm", "in1")]),
        (inputnode, resample_parc, [("asl_mask", "reference_image")]),
        (mrg_xfms, resample_parc, [("out", "transforms")]),
        (resample_parc, cbftssummary, [("output_image", "seg_file")]),
        (inputnode, cbftssummary, [
            ("cbf_ts", "cbf_ts"),
            ("confounds_file", "conf_file"),
            ("scoreindex", "score_file"),
        ]),
        (cbftssummary, ds_report_cbftsplot, [("out_file", "in_file")]),
        (cbftssummary, outputnode, [("out_file", "cbf_carpetplot")]),
        (inputnode, cbfsummary, [
            ("cbf", "cbf"),
            ("asl_ref", "ref_vol"),
        ]),
        (cbfsummary, ds_report_cbfplot, [("out_file", "in_file")]),
        (cbfsummary, outputnode, [("out_file", "cbf_summary_plot")]),
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
                ("score", "cbf"),
                ("asl_ref", "ref_vol"),
            ]),
            (scoresummary, ds_report_scoreplot, [("out_file", "in_file")]),
            (scoresummary, outputnode, [("out_file", "score_summary_plot")]),
            (inputnode, scrubsummary, [
                ("scrub", "cbf"),
                ("asl_ref", "ref_vol"),
            ]),
            (scrubsummary, ds_report_scrubplot, [("out_file", "in_file")]),
            (scrubsummary, outputnode, [("out_file", "scrub_summary_plot")]),
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
                ("basil", "cbf"),
                ("asl_ref", "ref_vol"),
            ]),
            (basilsummary, ds_report_basilplot, [("out_file", "in_file")]),
            (basilsummary, outputnode, [("out_file", "basil_summary_plot")]),
            (inputnode, pvcsummary, [
                ("pvc", "cbf"),
                ("asl_ref", "ref_vol"),
            ]),
            (pvcsummary, ds_report_pvcplot, [("out_file", "in_file")]),
            (pvcsummary, outputnode, [("out_file", "pvc_summary_plot")]),
        ])
        # fmt:on

    return workflow


def init_gecbfplot_wf(scorescrub=False, basil=False, name="cbf_plot"):
    """Plot CBF results for GE data.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.plotting import init_gecbfplot_wf

            wf = init_gecbfplot_wf()
    """
    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["cbf", "score", "scrub", "asl_ref", "basil", "pvc"]),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "cbf_summary_plot",
                "score_summary_plot",
                "scrub_summary_plot",
                "basil_summary_plot",
                "pvc_summary_plot",
            ]
        ),
        name="outputnode",
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
            ("cbf", "cbf"),
            ("asl_ref", "ref_vol"),
        ]),
        (cbfsummary, ds_report_cbfplot, [("out_file", "in_file")]),
        (cbfsummary, outputnode, [("out_file", "cbf_summary_plot")]),
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
                ("score", "cbf"),
                ("asl_ref", "ref_vol"),
            ]),
            (scoresummary, ds_report_scoreplot, [("out_file", "in_file")]),
            (scoresummary, outputnode, [("out_file", "score_summary_plot")]),
            (inputnode, scrubsummary, [
                ("scrub", "cbf"),
                ("asl_ref", "ref_vol"),
            ]),
            (scrubsummary, ds_report_scrubplot, [("out_file", "in_file")]),
            (scrubsummary, outputnode, [("out_file", "scrub_summary_plot")]),
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
                ("basil", "cbf"),
                ("asl_ref", "ref_vol"),
            ]),
            (basilsummary, ds_report_basilplot, [("out_file", "in_file")]),
            (basilsummary, outputnode, [("out_file", "basil_summary_plot")]),
            (inputnode, pvcsummary, [
                ("pvc", "cbf"),
                ("asl_ref", "ref_vol"),
            ]),
            (pvcsummary, ds_report_pvcplot, [("out_file", "in_file")]),
            (pvcsummary, outputnode, [("out_file", "pvc_summary_plot")]),
        ])
        # fmt:on

    return workflow