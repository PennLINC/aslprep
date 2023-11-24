# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for calculating CBF QC metrics."""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from templateflow.api import get as get_template

from aslprep import config
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.bids import DerivativesDataSink
from aslprep.interfaces.qc import ComputeCBFQC


def init_cbf_qc_wf(
    is_ge,
    scorescrub=False,
    basil=False,
    name="cbf_qc_wf",
):
    """Create a workflow for :abbr:`dolui2017automated (compute cbf)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep.workflows.asl.qc import init_cbf_qc_wf

            with mock_config():
                wf = init_cbf_qc_wf(
                    is_ge=False,
                    scorescrub=True,
                    basil=True,
                    name="cbf_qc_wf",
                )

    Parameters
    ----------
    is_ge : bool
    scorescrub : bool
    basil : bool
    name : :obj:`str`
        Name of workflow (default: "cbf_qc_wf")

    Inputs
    ------
    *cbf
        all cbf
    asl_mask
        asl mask NIFTI file
    t1w_tpms
        t1w probability maps
    aslref2anat_xfm
        aslref to t1w transformation file
    asl_mask_std : list
        Since ASLPrep always includes MNI152NLin2009cAsym as a standard space,
        this should always be provided.

    Outputs
    -------
    qc_file
        qc measures in tsv
    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """
The quality evaluation index (QEI) was computed for each CBF map [@dolui2017automated].
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
                "aslref2anat_xfm",
                "mni2009c2anat_xfm",
                # CBF inputs
                "mean_cbf",
                # SCORE/SCRUB inputs
                "mean_cbf_score",
                "mean_cbf_scrub",
                # BASIL inputs
                "mean_cbf_basil",
                "mean_cbf_gm_basil",
                "mean_cbf_wm_basil",
                # non-GE inputs
                "confounds_file",
                "rmsd_file",
            ],
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
        ApplyTransforms(
            interpolation="NearestNeighbor",
            float=True,
            invert_transform_flags=[True],
        ),
        name="gm_tfm",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, gm_tfm, [
            ("asl_mask", "reference_image"),
            ("aslref2anat_xfm", "transforms"),
            (("t1w_tpms", _pick_gm), "input_image"),
        ]),
    ])
    # fmt:on

    wm_tfm = pe.Node(
        ApplyTransforms(
            interpolation="NearestNeighbor",
            float=True,
            invert_transform_flags=[True],
        ),
        name="wm_tfm",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, wm_tfm, [
            ("asl_mask", "reference_image"),
            ("aslref2anat_xfm", "transforms"),
            (("t1w_tpms", _pick_wm), "input_image"),
        ]),
    ])
    # fmt:on

    csf_tfm = pe.Node(
        ApplyTransforms(
            interpolation="NearestNeighbor",
            float=True,
            invert_transform_flags=[True],
        ),
        name="csf_tfm",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, csf_tfm, [
            ("asl_mask", "reference_image"),
            ("aslref2anat_xfm", "transforms"),
            (("t1w_tpms", _pick_csf), "input_image"),
        ]),
    ])
    # fmt:on

    warp_t1w_mask_to_aslref = pe.Node(
        ApplyTransforms(
            interpolation="NearestNeighbor",
            float=True,
            invert_transform_flags=[True],
        ),
        name="warp_t1w_mask_to_aslref",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, warp_t1w_mask_to_aslref, [
            ("asl_mask", "reference_image"),
            ("aslref2anat_xfm", "transforms"),
            ("t1w_mask", "input_image"),
        ]),
    ])
    # fmt:on

    template_brain_mask = str(
        get_template("MNI152NLin2009cAsym", resolution=2, desc="brain", suffix="mask")
    )

    aslref2mni152nlin2009casym = pe.Node(niu.Merge(2), name="aslref2mni152nlin2009casym")
    workflow.connect([
        (inputnode, aslref2mni152nlin2009casym, [
            ("aslref2anat_xfm", "in1"),
            ("mni2009c2anat_xfm", "in2"),
        ]),
    ])  # fmt:skip

    warp_asl_mask_to_mni152nlin2009casym = pe.Node(
        ApplyTransforms(
            interpolation="NearestNeighbor",
            float=True,
            invert_transform_flags=[False, True],
            reference_image=template_brain_mask,
        ),
        name="warp_asl_mask_to_mni152nlin2009casym",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, warp_asl_mask_to_mni152nlin2009casym, [
            ("asl_mask", "reference_image"),
            ("asl_mask", "input_image"),
        ]),
        (aslref2mni152nlin2009casym, warp_asl_mask_to_mni152nlin2009casym, [
            ("out", "transforms"),
        ])
    ])
    # fmt:on

    compute_qc_metrics = pe.Node(
        ComputeCBFQC(
            tpm_threshold=0.7,
            template_mask=template_brain_mask,
        ),
        name="compute_qc_metrics",
        run_without_submitting=True,
        mem_gb=0.2,
    )

    # fmt:off
    workflow.connect([
        (warp_t1w_mask_to_aslref, compute_qc_metrics, [("output_image", "t1w_mask")]),
        (inputnode, compute_qc_metrics, [
            ("name_source", "name_source"),
            ("asl_mask", "asl_mask"),
            ("mean_cbf", "mean_cbf"),
        ]),
        (warp_asl_mask_to_mni152nlin2009casym, compute_qc_metrics, [
            ("output_image", "asl_mask_std"),
        ]),
        (gm_tfm, compute_qc_metrics, [("output_image", "gm_tpm")]),
        (wm_tfm, compute_qc_metrics, [("output_image", "wm_tpm")]),
        (csf_tfm, compute_qc_metrics, [("output_image", "csf_tpm")]),
        (compute_qc_metrics, outputnode, [("qc_file", "qc_file")]),
    ])
    # fmt:on

    ds_qc_metadata = pe.Node(
        DerivativesDataSink(
            base_directory=config.execution.aslprep_dir,
            dismiss_entities=list(DerivativesDataSink._allowed_entities),
            allowed_entities=[],
            suffix="qc",
            extension=".json",
        ),
        name="ds_qc_metadata",
        run_without_submitting=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, ds_qc_metadata, [("name_source", "source_file")]),
        (compute_qc_metrics, ds_qc_metadata, [("qc_metadata", "in_file")]),
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
