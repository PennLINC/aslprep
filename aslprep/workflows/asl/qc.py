# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for calculating CBF QC metrics."""
from nipype.interfaces import utility as niu
from nipype.interfaces.afni import Resample
from nipype.pipeline import engine as pe
from templateflow.api import get as get_template

from aslprep.interfaces.cbf_computation import ComputeCBFQC, ComputeCBFQCforGE
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.interfaces.fixes import (
    FixHeaderApplyTransforms as ApplyTransforms,
)


def init_cbfqc_compt_wf(
    asl_file,
    scorescrub=False,
    basil=False,
    name="cbfqc_compt_wf",
):
    """Create a workflow for :abbr:`dolui2017automated (compute cbf)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.qc import init_cbfqc_compt_wf

            wf = init_cbfqc_compt_wf(
                asl_file="",
            )

    Parameters
    ----------
    metadata : :obj:`dict`
        BIDS metadata for asl file
    name : :obj:`str`
        Name of workflow (default: ``cbfqc_compt_wf'``)

    Inputs
    ------
    *cbf
        all cbf
    asl_mask
        asl mask NIFTI file
    t1w_tpms
        t1w probability maps
    t1_asl_xform
        t1w to asl transfromation file

    Outputs
    -------
    qc_file
       qc measures in tsv
    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """
The Quality evaluation index (QEI) was computed for each CBF map [@dolui2017automated].
QEI is based on the similarity between the CBF and the structural images, the spatial
variability of the CBF image, and the percentage of grey matter voxels containing
negative CBF values.
"""
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "meancbf",
                "avgscore",
                "scrub",
                "basil",
                "asl_mask",
                "t1w_tpms",
                "confmat",
                "asl_mask_std",
                "t1_asl_xform",
                "pv",
                "t1w_mask",
                "rmsd_file",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=["qc_file"]), name="outputnode")

    def _pick_gm(files):
        return files[0]

    def _pick_wm(files):
        return files[1]

    def _pick_csf(files):
        return files[-1]

    gm_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="gm_tfm",
        mem_gb=0.1,
    )
    wm_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="wm_tfm",
        mem_gb=0.1,
    )
    csf_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="csf_tfm",
        mem_gb=0.1,
    )

    mask_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="masktonative",
        mem_gb=0.1,
    )

    brain_mask = str(
        get_template("MNI152NLin2009cAsym", resolution=2, desc="brain", suffix="mask")
    )

    resample = pe.Node(
        Resample(in_file=brain_mask, outputtype="NIFTI_GZ"),
        name="resample",
        mem_gb=0.1,
    )

    qccompute = pe.Node(
        ComputeCBFQC(in_file=asl_file),
        name="qccompute",
        run_without_submitting=True,
        mem_gb=0.2,
    )

    # fmt:off
    workflow.connect([
        (inputnode, csf_tfm, [
            ("asl_mask", "reference_image"),
            ("t1_asl_xform", "transforms"),
        ]),
        (inputnode, csf_tfm, [(("t1w_tpms", _pick_csf), "input_image")]),
        (inputnode, wm_tfm, [
            ("asl_mask", "reference_image"),
            ("t1_asl_xform", "transforms"),
        ]),
        (inputnode, wm_tfm, [(("t1w_tpms", _pick_wm), "input_image")]),
        (inputnode, gm_tfm, [
            ("asl_mask", "reference_image"),
            ("t1_asl_xform", "transforms"),
        ]),
        (inputnode, gm_tfm, [(("t1w_tpms", _pick_gm), "input_image")]),
        (inputnode, mask_tfm, [
            ("asl_mask", "reference_image"),
            ("t1_asl_xform", "transforms"),
            ("t1w_mask", "input_image"),
        ]),
        (mask_tfm, qccompute, [("output_image", "in_t1mask")]),
        (inputnode, qccompute, [
            ("asl_mask", "in_aslmask"),
            ("confmat", "in_confmat"),
        ]),
        (inputnode, qccompute, [(("asl_mask_std", _pick_csf), "in_aslmaskstd")]),
        (inputnode, qccompute, [("rmsd_file", "rmsd_file")]),
        (inputnode, resample, [(("asl_mask_std", _pick_csf), "master")]),
        (resample, qccompute, [("out_file", "in_templatemask")]),
        (gm_tfm, qccompute, [("output_image", "in_greyM")]),
        (wm_tfm, qccompute, [("output_image", "in_whiteM")]),
        (csf_tfm, qccompute, [("output_image", "in_csf")]),
        (inputnode, qccompute, [("meancbf", "in_meancbf")]),
    ])
    # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (inputnode, qccompute, [
                ("scrub", "in_scrub"),
                ("avgscore", "in_avgscore"),
            ]),
        ])
        # fmt:on

    if basil:
        # fmt:off
        workflow.connect([
            (inputnode, qccompute, [
                ("basil", "in_basil"),
                ("pv", "in_pvc"),
            ]),
        ])
        # fmt:on

    # fmt:off
    workflow.connect([(qccompute, outputnode, [("qc_file", "qc_file")])])
    # fmt:on

    return workflow


def init_cbfgeqc_compt_wf(
    asl_file,
    scorescrub=False,
    basil=False,
    name="cbfqc_compt_wf",
):
    """Create a workflow for :abbr:`dolui2017automated (compute cbf)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.qc import init_cbfgeqc_compt_wf

            wf = init_cbfgeqc_compt_wf(
                asl_file="",
            )

    Parameters
    ----------
    metadata : :obj:`dict`
        BIDS metadata for asl file
    name : :obj:`str`
        Name of workflow (default: ``cbfqc_compt_wf'``)

    Inputs
    ------
    *cbf
        all cbf
    asl_mask
        asl mask NIFTI file
    t1w_tpms
        t1w probability maps
    t1_asl_xform
        t1w to asl transfromation file

    Outputs
    -------
    qc_file
       qc measures in tsv
    """
    workflow = Workflow(name=name)
    # workflow.__desc__ = """\

    # """
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "meancbf",
                "avgscore",
                "scrub",
                "basil",
                "asl_mask",
                "t1w_tpms",
                "asl_mask_std",
                "t1_asl_xform",
                "pv",
                "t1w_mask",
                "rmsd_file",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=["qc_file"]), name="outputnode")

    def _pick_csf(files):
        return files[-1]

    def _pick_gm(files):
        return files[0]

    def _pick_wm(files):
        return files[1]

    csf_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="csf_tfm",
        mem_gb=0.1,
    )
    wm_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="wm_tfm",
        mem_gb=0.1,
    )
    gm_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="gm_tfm",
        mem_gb=0.1,
    )

    mask_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="masktonative",
        mem_gb=0.1,
    )

    brain_mask = str(
        get_template("MNI152NLin2009cAsym", resolution=2, desc="brain", suffix="mask")
    )

    resample = pe.Node(
        Resample(in_file=brain_mask, outputtype="NIFTI_GZ"),
        name="resample",
        mem_gb=0.1,
    )

    qccompute = pe.Node(
        ComputeCBFQCforGE(in_file=asl_file),
        name="qccompute",
        run_without_submitting=True,
        mem_gb=0.2,
    )

    # fmt:off
    workflow.connect([
        (inputnode, csf_tfm, [
            ("asl_mask", "reference_image"),
            ("t1_asl_xform", "transforms"),
        ]),
        (inputnode, csf_tfm, [(("t1w_tpms", _pick_csf), "input_image")]),
        (inputnode, wm_tfm, [
            ("asl_mask", "reference_image"),
            ("t1_asl_xform", "transforms"),
        ]),
        (inputnode, wm_tfm, [(("t1w_tpms", _pick_wm), "input_image")]),
        (inputnode, gm_tfm, [
            ("asl_mask", "reference_image"),
            ("t1_asl_xform", "transforms"),
        ]),
        (inputnode, gm_tfm, [(("t1w_tpms", _pick_gm), "input_image")]),
        (inputnode, mask_tfm, [
            ("asl_mask", "reference_image"),
            ("t1_asl_xform", "transforms"),
            ("t1w_mask", "input_image"),
        ]),
        (mask_tfm, qccompute, [("output_image", "in_t1mask")]),
        (inputnode, qccompute, [("asl_mask", "in_aslmask")]),
        (inputnode, qccompute, [(("asl_mask_std", _pick_csf), "in_aslmaskstd")]),
        (inputnode, resample, [(("asl_mask_std", _pick_csf), "master")]),
        (resample, qccompute, [("out_file", "in_templatemask")]),
        (gm_tfm, qccompute, [("output_image", "in_greyM")]),
        (wm_tfm, qccompute, [("output_image", "in_whiteM")]),
        (csf_tfm, qccompute, [("output_image", "in_csf")]),
        (inputnode, qccompute, [("meancbf", "in_meancbf")]),
    ])
    # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (inputnode, qccompute, [
                ("scrub", "in_scrub"),
                ("avgscore", "in_avgscore"),
            ]),
        ])
        # fmt:on

    if basil:
        # fmt:off
        workflow.connect([(inputnode, qccompute, [("basil", "in_basil"), ("pv", "in_pvc")])])
        # fmt:on

    # fmt:off
    workflow.connect([(qccompute, outputnode, [("qc_file", "qc_file")])])
    # fmt:on
    return workflow
