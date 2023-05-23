# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows to process GE ASL data."""
from nipype.interfaces import fsl
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.interfaces.ge import GeReferenceFile
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.interfaces.masks import SimpleShowMaskRPT
from aslprep.workflows.asl.registration import init_fsl_bbr_wf

DEFAULT_MEMORY_MIN_GB = config.DEFAULT_MEMORY_MIN_GB
LOGGER = config.loggers.workflow


def init_asl_reference_ge_wf(
    metadata,
    aslcontext,
    smooth_kernel=5,
    name="asl_reference_ge_wf",
):
    """Generate a reference volume and its skull-stripped version.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.ge_utils import init_asl_reference_ge_wf

            wf = init_asl_reference_ge_wf(
                metadata={},
                aslcontext="sub-01_aslcontext.tsv",
            )
    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
First, a reference volume and its skull-stripped version were generated.
    """

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "m0scan",
                "m0scan_metadata",
            ],
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "raw_ref_image",
                "ref_image_brain",
                "asl_mask",
                "m0_file",
                "m0tr",
                "mask_report",
            ]
        ),
        name="outputnode",
    )

    gen_ref = pe.Node(
        GeReferenceFile(fwhm=smooth_kernel, metadata=metadata, aslcontext=aslcontext),
        omp_nthreads=1,
        mem_gb=1,
        name="gen_ge_ref",
    )

    # fmt:off
    workflow.connect([
        (inputnode, gen_ref, [
            ("asl_file", "in_file"),
            ("m0scan", "m0scan"),
            ("m0scan_metadata", "m0scan_metadata"),
        ]),
        (gen_ref, outputnode, [
            ("ref_file", "raw_ref_image"),
            ("m0_file", "m0_file"),
            ("m0tr", "m0tr"),
        ]),
    ])
    # fmt:on

    skull_strip_wf = pe.Node(fsl.BET(frac=0.5, mask=True), name="fslbet")

    # fmt:off
    workflow.connect([
        (gen_ref, skull_strip_wf, [("ref_file", "in_file")]),
        (skull_strip_wf, outputnode, [("mask_file", "asl_mask")]),
    ])
    # fmt:on

    apply_mask = pe.Node(fsl.ApplyMask(), name="apply_mask")

    # fmt:off
    workflow.connect([
        (gen_ref, apply_mask, [("ref_file", "in_file")]),
        (skull_strip_wf, apply_mask, [("mask_file", "mask_file")]),
        (apply_mask, outputnode, [("out_file", "ref_image_brain")]),
    ])
    # fmt:on

    mask_reportlet = pe.Node(SimpleShowMaskRPT(), name="mask_reportlet")

    # fmt:off
    workflow.connect([
        (gen_ref, mask_reportlet, [("ref_file", "background_file")]),
        (skull_strip_wf, mask_reportlet, [("mask_file", "mask_file")]),
    ])
    # fmt:on

    return workflow


def init_asl_reg_ge_wf(
    use_bbr,
    asl2t1w_dof,
    asl2t1w_init,
    sloppy=False,
    write_report=True,
    name="asl_reg_ge_wf",
):
    """Calculate registration transforms from ASL reference volume to T1w space.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.ge_utils import init_asl_reg_ge_wf

            wf = init_asl_reg_ge_wf(
                use_bbr=True,
                asl2t1w_dof=9,
                asl2t1w_init="register",
            )
    """
    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=["ref_asl_brain", "t1w_brain", "t1w_dseg"]),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["aslref_to_anat_xfm", "anat_to_aslref_xfm", "fallback"]),
        name="outputnode",
    )

    bbr_wf = init_fsl_bbr_wf(
        use_bbr=use_bbr,
        asl2t1w_dof=asl2t1w_dof,
        asl2t1w_init=asl2t1w_init,
        sloppy=sloppy,
    )

    # fmt:off
    workflow.connect([
        (inputnode, bbr_wf, [
            ("ref_asl_brain", "inputnode.in_file"),
            ("t1w_dseg", "inputnode.t1w_dseg"),
            ("t1w_brain", "inputnode.t1w_brain"),
        ]),
        (bbr_wf, outputnode, [
            ("outputnode.aslref_to_anat_xfm", "aslref_to_anat_xfm"),
            ("outputnode.anat_to_aslref_xfm", "anat_to_aslref_xfm"),
            ("outputnode.fallback", "fallback"),
        ]),
    ])
    # fmt:on

    if write_report:
        ds_report_reg = pe.Node(
            DerivativesDataSink(datatype="figures", dismiss_entities=("echo",)),
            name="ds_report_reg",
            run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB,
        )

        def _asl_reg_suffix(fallback):  # noqa: U100
            return "flirtbbr"

        # fmt:off
        workflow.connect([
            (bbr_wf, ds_report_reg, [
                ("outputnode.out_report", "in_file"),
                (("outputnode.fallback", _asl_reg_suffix), "desc"),
            ]),
        ])
        # fmt:on

    return workflow
