# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for calculating CBF QC metrics."""
from nipype.interfaces import utility as niu
from nipype.interfaces.afni import Resample
from nipype.pipeline import engine as pe
from templateflow.api import get as get_template

from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.qc import ComputeCBFQC
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.utils.misc import _select_last_in_list


def init_compute_cbf_qc_wf(
    is_ge,
    scorescrub=False,
    basil=False,
    name="compute_cbf_qc_wf",
):
    """Create a workflow for :abbr:`dolui2017automated (compute cbf)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.qc import init_compute_cbf_qc_wf

            wf = init_compute_cbf_qc_wf(
                is_ge=False,
                scorescrub=True,
                basil=True,
                name="compute_cbf_qc_wf",
            )

    Parameters
    ----------
    is_ge : bool
    scorescrub : bool
    basil : bool
    name : :obj:`str`
        Name of workflow (default: "compute_cbf_qc_wf")

    Inputs
    ------
    *cbf
        all cbf
    asl_mask
        asl mask NIFTI file
    t1w_tpms
        t1w probability maps
    anat_to_aslref_xfm
        t1w to asl transformation file

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
                "name_source",
                "asl_mask",
                "t1w_mask",
                "t1w_tpms",
                "anat_to_aslref_xfm",
                "asl_mask_std",
                # CBF inputs
                "mean_cbf",
                # SCORE/SCRUB inputs
                "mean_cbf_score",
                "mean_cbf_scrub",
                # BASIL inputs
                "mean_cbf_basil",
                "mean_cbf_gm_basil",
                # non-GE inputs
                "confounds_file",
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
        return files[2]

    gm_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="gm_tfm",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, gm_tfm, [
            ("asl_mask", "reference_image"),
            ("anat_to_aslref_xfm", "transforms"),
            (("t1w_tpms", _pick_gm), "input_image"),
        ]),
    ])
    # fmt:on

    wm_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="wm_tfm",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, wm_tfm, [
            ("asl_mask", "reference_image"),
            ("anat_to_aslref_xfm", "transforms"),
            (("t1w_tpms", _pick_wm), "input_image"),
        ]),
    ])
    # fmt:on

    csf_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="csf_tfm",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, csf_tfm, [
            ("asl_mask", "reference_image"),
            ("anat_to_aslref_xfm", "transforms"),
            (("t1w_tpms", _pick_csf), "input_image"),
        ]),
    ])
    # fmt:on

    mask_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="masktonative",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, mask_tfm, [
            ("asl_mask", "reference_image"),
            ("anat_to_aslref_xfm", "transforms"),
            ("t1w_mask", "input_image"),
        ]),
    ])
    # fmt:on

    brain_mask = str(
        get_template("MNI152NLin2009cAsym", resolution=2, desc="brain", suffix="mask")
    )

    resample = pe.Node(
        Resample(in_file=brain_mask, outputtype="NIFTI_GZ"),
        name="resample",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([(inputnode, resample, [(("asl_mask_std", _select_last_in_list), "master")])])
    # fmt:on

    compute_qc_metrics = pe.Node(
        ComputeCBFQC(tpm_threshold=0.7),
        name="compute_qc_metrics",
        run_without_submitting=True,
        mem_gb=0.2,
    )

    # fmt:off
    workflow.connect([
        (mask_tfm, compute_qc_metrics, [("output_image", "t1w_mask")]),
        (inputnode, compute_qc_metrics, [
            ("name_source", "name_source"),
            ("asl_mask", "asl_mask"),
            (("asl_mask_std", _select_last_in_list), "asl_mask_std"),
            ("mean_cbf", "mean_cbf"),
        ]),
        (resample, compute_qc_metrics, [("out_file", "template_mask")]),
        (gm_tfm, compute_qc_metrics, [("output_image", "gm_tpm")]),
        (wm_tfm, compute_qc_metrics, [("output_image", "wm_tpm")]),
        (csf_tfm, compute_qc_metrics, [("output_image", "csf_tpm")]),
        (compute_qc_metrics, outputnode, [("qc_file", "qc_file")]),
    ])
    # fmt:on

    if not is_ge:
        # The QC node only expects a confounds file and RMSD file for non-GE data.
        # fmt:off
        workflow.connect([
            (inputnode, compute_qc_metrics, [
                ("confounds_file", "confounds_file"),
                ("rmsd_file", "rmsd_file"),
            ]),
        ])
        # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (inputnode, compute_qc_metrics, [
                ("mean_cbf_scrub", "mean_cbf_scrub"),
                ("mean_cbf_score", "mean_cbf_score"),
            ]),
        ])
        # fmt:on

    if basil:
        # fmt:off
        workflow.connect([
            (inputnode, compute_qc_metrics, [
                ("mean_cbf_basil", "mean_cbf_basil"),
                ("mean_cbf_gm_basil", "mean_cbf_gm_basil"),
            ]),
        ])
        # fmt:on

    return workflow
