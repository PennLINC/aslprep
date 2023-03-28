# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for calculating CBF."""
import pandas as pd
from nipype.interfaces import utility as niu
from nipype.interfaces.afni import Resample
from nipype.interfaces.base import Undefined
from nipype.interfaces.fsl import Info, MultiImageMaths
from nipype.pipeline import engine as pe
from templateflow.api import get as get_template

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.interfaces.cbf_computation import (
    BASILCBF,
    ComputeCBF,
    ComputeCBFQC,
    ComputeCBFQCforGE,
    ExtractCBF,
    ExtractCBForDeltaM,
    ParcellateCBF,
    RefineMask,
    ScoreAndScrubCBF,
)
from aslprep.interfaces.plotting import CBFSummary, CBFtsSummary
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.interfaces.fixes import (
    FixHeaderApplyTransforms as ApplyTransforms,
)
from aslprep.utils.misc import get_atlas, get_tis, pcasl_or_pasl


def init_cbf_compt_wf(
    name_source,
    metadata,
    bids_dir,
    dummy_vols,
    scorescrub=False,
    basil=False,
    M0Scale=1,
    smooth_kernel=5,
    name="cbf_compt_wf",
):
    """Create a workflow for :abbr:`CCBF (compute cbf)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.cbf import init_cbf_compt_wf

            wf = init_cbf_compt_wf(
                name_source="",
                metadata={},
                bids_dir="",
                dummy_vols=0,
            )

    Parameters
    ----------
    name_source : :obj:`str`
        Path to the raw ASL file.
    metadata : :obj:`dict`
        BIDS metadata for asl file
    name : :obj:`str`
        Name of workflow (default: ``cbf_compt_wf``)

    Inputs
    ------
    in_file
        Raw asl file. Equivalent to name_source.
    bids_dir
        bids dir of the subject
    asl_file
        asl series NIfTI file
    asl_mask
        asl mask NIFTI file
    t1w_tpms
        t1w probability maps
    t1w_mask
        t1w mask Nifti
    t1_asl_xform
        t1w to asl transfromation file
    itk_asl_to_t1
        asl to t1q transfromation file
    scorescrub
        compute score and scrub
    basil
        compute basil and PVC

    Outputs
    -------
    *cbf
       all cbf outputs
       cbf,score, scrub, pv, and basil
    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """
### Cerebral blood flow computation and denoising

*ASLPrep* was configured to calculate cerebral blood flow (CBF) using the following methods.

The cerebral blood flow (CBF) was quantified from preprocessed ASL data using a general kinetic
model [@buxton1998general].
"""
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "in_file",
                "asl_mask",
                "t1w_tpms",
                "t1w_mask",
                "t1_asl_xform",
                "itk_asl_to_t1",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "out_cbf",
                "out_mean",
                "out_score",
                "out_cbfpvwm",
                "out_avgscore",
                "out_scrub",
                "out_cbfb",
                "out_scoreindex",
                "out_cbfpv",
                "out_att",
            ]
        ),
        name="outputnode",
    )

    refine_mask = pe.Node(
        RefineMask(),
        mem_gb=0.2,
        run_without_submitting=True,
        name="refine_mask",
    )

    # fmt:off
    workflow.connect([
        (inputnode, refine_mask, [
            ("t1w_mask", "in_t1mask"),
            ("asl_mask", "in_aslmask"),
            ("t1_asl_xform", "transforms"),
        ]),
    ])
    # fmt:on

    # convert tmps to asl_space
    def _pick_csf(files):
        return files[-1]

    def _pick_gm(files):
        return files[0]

    def _pick_wm(files):
        return files[1]

    def _getfiledir(file):
        import os

        return os.path.dirname(file)

    csf_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="csf_tfm",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, csf_tfm, [
            ("asl_mask", "reference_image"),
            ("t1_asl_xform", "transforms"),
            (("t1w_tpms", _pick_csf), "input_image"),
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
            ("t1_asl_xform", "transforms"),
            (("t1w_tpms", _pick_wm), "input_image"),
        ]),
    ])
    # fmt:on

    gm_tfm = pe.Node(
        ApplyTransforms(interpolation="NearestNeighbor", float=True),
        name="gm_tfm",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, gm_tfm, [
            ("asl_mask", "reference_image"),
            ("t1_asl_xform", "transforms"),
            (("t1w_tpms", _pick_gm), "input_image"),
        ]),
    ])
    # fmt:on

    tiscbf = get_tis(metadata)

    is_casl = pcasl_or_pasl(metadata=metadata)
    # XXX: This is a bad way to find this file.
    aslcontext = name_source.replace("_asl.nii.gz", "_aslcontext.tsv")
    aslcontext_df = pd.read_table(aslcontext)
    cbf_only = all(aslcontext_df["volume_type"].isin(("m0scan", "cbf")))
    if cbf_only and not basil:
        config.loggers.workflow.info(f"Only CBF volumes are detected in {name_source}.")
    elif cbf_only:
        config.loggers.workflow.warning(
            f"Only CBF volumes are detected in {name_source}. "
            "BASIL will automatically be disabled."
        )
        basil = False

    extract_deltam = pe.Node(
        ExtractCBF(
            aslcontext=aslcontext,
            dummy_vols=dummy_vols,
            bids_dir=bids_dir,
            fwhm=smooth_kernel,
            in_metadata=metadata,
        ),
        mem_gb=0.2,
        run_without_submitting=True,
        name="extract_deltam",
    )

    # fmt:off
    workflow.connect([
        (inputnode, extract_deltam, [("asl_file", "asl_file")]),
        (refine_mask, extract_deltam, [("out_mask", "in_mask")]),
    ])
    # fmt:on

    compute_cbf = pe.Node(
        ComputeCBF(
            cbf_only=cbf_only,
            m0scale=M0Scale,
        ),
        mem_gb=0.2,
        run_without_submitting=True,
        name="compute_cbf",
    )

    # fmt:off
    workflow.connect([
        (refine_mask, compute_cbf, [("out_mask", "mask")]),
        (extract_deltam, compute_cbf, [
            ("out_file", "deltam"),
            ("m0_file", "m0_file"),
            ("metadata", "metadata"),
        ]),
        (compute_cbf, outputnode, [
            ("cbf", "out_cbf"),
            ("mean_cbf", "out_mean"),
        ]),
    ])
    # fmt:on

    if scorescrub:
        score_and_scrub_cbf = pe.Node(
            ScoreAndScrubCBF(in_thresh=0.7, in_wfun="huber"),
            mem_gb=0.2,
            name="score_and_scrub_cbf",
            run_without_submitting=True,
        )

        workflow.__desc__ += """

Structural Correlation based Outlier Rejection (SCORE) algorithm was applied to the CBF timeseries
to discard CBF volumes with outlying values [@dolui2017structural] before computing the mean CBF.
Following SCORE, the Structural Correlation with RobUst Bayesian (SCRUB) algorithm was applied to
the CBF maps using structural tissue probability maps to reweight the mean CBF
[@dolui2017structural;@dolui2016scrub].
"""
        # fmt:off
        workflow.connect([
            (refine_mask, score_and_scrub_cbf, [("out_mask", "in_mask")]),
            (compute_cbf, score_and_scrub_cbf, [("cbf", "in_file")]),
            (gm_tfm, score_and_scrub_cbf, [("output_image", "in_greyM")]),
            (wm_tfm, score_and_scrub_cbf, [("output_image", "in_whiteM")]),
            (csf_tfm, score_and_scrub_cbf, [("output_image", "in_csf")]),
            (score_and_scrub_cbf, outputnode, [
                ("out_score", "out_score"),
                ("out_scoreindex", "out_scoreindex"),
                ("out_avgscore", "out_avgscore"),
                ("out_scrub", "out_scrub"),
            ]),
        ])
        # fmt:on

    if basil:
        workflow.__desc__ += f"""

CBF was also computed with Bayesian Inference for Arterial Spin Labeling (BASIL)
[@chappell2008variational], as implemented in *FSL* {Info.version().split(':')[0]}.
BASIL computes CBF using a spatial regularization of the estimated perfusion image and
additionally calculates a partial-volume corrected CBF image [@chappell_pvc].
"""
        if is_casl:
            bolus = metadata["LabelingDuration"]
        else:  # pasl
            bolus = metadata["BolusCutOffDelayTime"]

        # NOTE: Do we need the TR of just the M0 volumes?
        # TODO: Extract M0 TR in extract_deltam.
        basilcbf = pe.Node(
            BASILCBF(
                m0scale=M0Scale,
                bolus=bolus,
                m0tr=metadata["RepetitionTime"],
                alpha=metadata.get("LabelingEfficiency", Undefined),
                pvc=True,
                tis=tiscbf,
                pcasl=is_casl,
            ),
            name="basilcbf",
            run_without_submitting=True,
            mem_gb=0.2,
        )

        # fmt:off
        workflow.connect([
            (refine_mask, basilcbf, [("out_mask", "mask")]),
            (extract_deltam, basilcbf, [
                (("m0_file", _getfiledir), "out_basename"),
                ("out_file", "in_file"),
                ("m0_file", "mzero"),
            ]),
            (gm_tfm, basilcbf, [("output_image", "pvgm")]),
            (wm_tfm, basilcbf, [("output_image", "pvwm")]),
            (basilcbf, outputnode, [
                ("out_cbfb", "out_cbfb"),
                ("out_cbfpv", "out_cbfpv"),
                ("out_cbfpvwm", "out_cbfpvwm"),
                ("out_att", "out_att"),
            ]),
        ])
        # fmt:on

    return workflow


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

            from aslprep.workflows.asl.cbf import init_cbfqc_compt_wf

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

            from aslprep.workflows.asl.cbf import init_cbfgeqc_compt_wf

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

            from aslprep.workflows.asl.cbf import init_cbfplot_wf

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
        CBFtsSummary(tr=metadata["RepetitionTime"]), name="cbf_ts_summary", mem_gb=2
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

            from aslprep.workflows.asl.cbf import init_gecbfplot_wf

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


def init_cbfroiquant_wf(scorescrub=False, basil=False, name="cbf_roiquant"):
    """Parcellate CBF results using a set of atlases.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.cbf import init_cbfroiquant_wf

            wf = init_cbfroiquant_wf()
    """
    workflow = Workflow(name=name)

    workflow.__desc__ = """\
For each CBF map, the ROIs for the following atlases were extracted:
the Harvard-Oxford and the Schaefer 200 and 400-parcel resolution atlases.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "cbf",
                "score",
                "scrub",
                "basil",
                "pvc",
                "aslmask",
                "t1_asl_xform",
                "std2anat_xfm",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "cbf_hvoxf",
                "score_hvoxf",
                "scrub_hvoxf",
                "basil_hvoxf",
                "pvc_hvoxf",
                "cbf_sc207",
                "score_sc207",
                "scrub_sc207",
                "basil_sc207",
                "pvc_sc207",
                "cbf_sc217",
                "score_sc217",
                "scrub_sc217",
                "basil_sc217",
                "pvc_sc217",
                "cbf_sc407",
                "score_sc407",
                "scrub_sc407",
                "basil_sc407",
                "pvc_sc407",
                "cbf_sc417",
                "score_sc417",
                "scrub_sc417",
                "basil_sc417",
                "pvc_sc417",
            ]
        ),
        name="outputnode",
    )

    hvoxfile, hvoxdata, hvoxlabel = get_atlas(atlasname="HarvardOxford")
    sc207file, sc207data, sc207label = get_atlas(atlasname="schaefer200x7")
    sc217file, sc217data, sc217label = get_atlas(atlasname="schaefer200x17")
    sc407file, sc407data, sc407label = get_atlas(atlasname="schaefer400x7")
    sc417file, sc417data, sc417label = get_atlas(atlasname="schaefer400x17")

    mrg_xfms = pe.Node(niu.Merge(2), name="mrg_xfms")
    hvoftrans = pe.Node(
        ApplyTransforms(
            float=True,
            input_image=hvoxfile,
            dimension=3,
            default_value=0,
            interpolation="NearestNeighbor",
        ),
        name="hvoftrans",
    )
    sc207trans = pe.Node(
        ApplyTransforms(
            float=True,
            input_image=sc207file,
            dimension=3,
            default_value=0,
            interpolation="NearestNeighbor",
        ),
        name="sc207trans",
    )
    sc217trans = pe.Node(
        ApplyTransforms(
            float=True,
            input_image=sc217file,
            dimension=3,
            default_value=0,
            interpolation="NearestNeighbor",
        ),
        name="sc217trans",
    )
    sc407trans = pe.Node(
        ApplyTransforms(
            float=True,
            input_image=sc407file,
            dimension=3,
            default_value=0,
            interpolation="NearestNeighbor",
        ),
        name="sc407trans",
    )
    sc417trans = pe.Node(
        ApplyTransforms(
            float=True,
            input_image=sc417file,
            dimension=3,
            default_value=0,
            interpolation="NearestNeighbor",
        ),
        name="sc417trans",
    )

    cbfroihv = pe.Node(ParcellateCBF(atlaslabel=hvoxlabel, atlasdata=hvoxdata), name="cbfroihv")
    cbfroi207 = pe.Node(ParcellateCBF(atlaslabel=sc207label, atlasdata=sc207data), name="cbf207")
    cbfroi217 = pe.Node(ParcellateCBF(atlaslabel=sc217label, atlasdata=sc217data), name="cbf217")
    cbfroi407 = pe.Node(ParcellateCBF(atlaslabel=sc407label, atlasdata=sc407data), name="cbf407")
    cbfroi417 = pe.Node(ParcellateCBF(atlaslabel=sc417label, atlasdata=sc417data), name="cbf417")
    if scorescrub:
        scorehv = pe.Node(ParcellateCBF(atlaslabel=hvoxlabel, atlasdata=hvoxdata), name="scorehv")
        score207 = pe.Node(
            ParcellateCBF(atlaslabel=sc207label, atlasdata=sc207data),
            name="score207",
        )
        score217 = pe.Node(
            ParcellateCBF(atlaslabel=sc217label, atlasdata=sc217data),
            name="score217",
        )
        score407 = pe.Node(
            ParcellateCBF(atlaslabel=sc407label, atlasdata=sc407data),
            name="score407",
        )
        score417 = pe.Node(
            ParcellateCBF(atlaslabel=sc417label, atlasdata=sc417data),
            name="score417",
        )
        scrubhv = pe.Node(ParcellateCBF(atlaslabel=hvoxlabel, atlasdata=hvoxdata), name="scrubhv")
        scrub207 = pe.Node(
            ParcellateCBF(atlaslabel=sc207label, atlasdata=sc207data),
            name="scrub207",
        )
        scrub217 = pe.Node(
            ParcellateCBF(atlaslabel=sc217label, atlasdata=sc217data),
            name="scrub217",
        )
        scrub407 = pe.Node(
            ParcellateCBF(atlaslabel=sc407label, atlasdata=sc407data),
            name="scrub407",
        )
        scrub417 = pe.Node(
            ParcellateCBF(atlaslabel=sc417label, atlasdata=sc417data),
            name="scrub417",
        )

    if basil:
        basilhv = pe.Node(ParcellateCBF(atlaslabel=hvoxlabel, atlasdata=hvoxdata), name="basilhv")
        basil207 = pe.Node(
            ParcellateCBF(atlaslabel=sc207label, atlasdata=sc207data),
            name="basil207",
        )
        basil217 = pe.Node(
            ParcellateCBF(atlaslabel=sc217label, atlasdata=sc217data),
            name="basil217",
        )
        basil407 = pe.Node(
            ParcellateCBF(atlaslabel=sc407label, atlasdata=sc407data),
            name="basil407",
        )
        basil417 = pe.Node(
            ParcellateCBF(atlaslabel=sc417label, atlasdata=sc417data),
            name="basil417",
        )
        pvchv = pe.Node(ParcellateCBF(atlaslabel=hvoxlabel, atlasdata=hvoxdata), name="pvchv")
        pvc207 = pe.Node(ParcellateCBF(atlaslabel=sc207label, atlasdata=sc207data), name="pvc207")
        pvc217 = pe.Node(ParcellateCBF(atlaslabel=sc217label, atlasdata=sc217data), name="pvc217")
        pvc407 = pe.Node(ParcellateCBF(atlaslabel=sc407label, atlasdata=sc407data), name="pvc407")
        pvc417 = pe.Node(ParcellateCBF(atlaslabel=sc417label, atlasdata=sc417data), name="pvc417")

    # fmt:off
    workflow.connect([
        (inputnode, mrg_xfms, [
            ("t1_asl_xform", "in2"),
            ("std2anat_xfm", "in1"),
        ]),
        (inputnode, hvoftrans, [("aslmask", "reference_image")]),
        (mrg_xfms, hvoftrans, [("out", "transforms")]),
        (inputnode, sc207trans, [("aslmask", "reference_image")]),
        (mrg_xfms, sc207trans, [("out", "transforms")]),
        (inputnode, sc217trans, [("aslmask", "reference_image")]),
        (mrg_xfms, sc217trans, [("out", "transforms")]),
        (inputnode, sc407trans, [("aslmask", "reference_image")]),
        (mrg_xfms, sc407trans, [("out", "transforms")]),
        (inputnode, sc417trans, [("aslmask", "reference_image")]),
        (mrg_xfms, sc417trans, [("out", "transforms")]),
        (hvoftrans, cbfroihv, [("output_image", "atlasfile")]),
        (sc207trans, cbfroi207, [("output_image", "atlasfile")]),
        (sc217trans, cbfroi217, [("output_image", "atlasfile")]),
        (sc407trans, cbfroi407, [("output_image", "atlasfile")]),
        (sc417trans, cbfroi417, [("output_image", "atlasfile")]),
        (inputnode, cbfroihv, [("cbf", "in_cbf")]),
        (inputnode, cbfroi207, [("cbf", "in_cbf")]),
        (inputnode, cbfroi217, [("cbf", "in_cbf")]),
        (inputnode, cbfroi407, [("cbf", "in_cbf")]),
        (inputnode, cbfroi417, [("cbf", "in_cbf")]),
        (cbfroihv, outputnode, [("atlascsv", "cbf_hvoxf")]),
        (cbfroi207, outputnode, [("atlascsv", "cbf_sc207")]),
        (cbfroi217, outputnode, [("atlascsv", "cbf_sc217")]),
        (cbfroi407, outputnode, [("atlascsv", "cbf_sc407")]),
        (cbfroi417, outputnode, [("atlascsv", "cbf_sc417")]),
    ])
    # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (hvoftrans, scorehv, [("output_image", "atlasfile")]),
            (hvoftrans, scrubhv, [("output_image", "atlasfile")]),
            (sc207trans, score207, [("output_image", "atlasfile")]),
            (sc207trans, scrub207, [("output_image", "atlasfile")]),
            (sc217trans, score217, [("output_image", "atlasfile")]),
            (sc217trans, scrub217, [("output_image", "atlasfile")]),
            (sc407trans, score407, [("output_image", "atlasfile")]),
            (sc407trans, scrub407, [("output_image", "atlasfile")]),
            (sc417trans, score417, [("output_image", "atlasfile")]),
            (sc417trans, scrub417, [("output_image", "atlasfile")]),
            (inputnode, scorehv, [("score", "in_cbf")]),
            (inputnode, scrubhv, [("scrub", "in_cbf")]),
            (inputnode, score207, [("score", "in_cbf")]),
            (inputnode, scrub207, [("scrub", "in_cbf")]),
            (inputnode, score217, [("score", "in_cbf")]),
            (inputnode, scrub217, [("scrub", "in_cbf")]),
            (inputnode, score407, [("score", "in_cbf")]),
            (inputnode, scrub407, [("scrub", "in_cbf")]),
            (inputnode, score417, [("score", "in_cbf")]),
            (inputnode, scrub417, [("scrub", "in_cbf")]),
            (scorehv, outputnode, [("atlascsv", "score_hvoxf")]),
            (score207, outputnode, [("atlascsv", "score_sc207")]),
            (score217, outputnode, [("atlascsv", "score_sc217")]),
            (score407, outputnode, [("atlascsv", "score_sc407")]),
            (score417, outputnode, [("atlascsv", "score_sc417")]),
            (scrubhv, outputnode, [("atlascsv", "scrub_hvoxf")]),
            (scrub207, outputnode, [("atlascsv", "scrub_sc207")]),
            (scrub217, outputnode, [("atlascsv", "scrub_sc217")]),
            (scrub407, outputnode, [("atlascsv", "scrub_sc407")]),
            (scrub417, outputnode, [("atlascsv", "scrub_sc417")]),
        ])
        # fmt:on

    if basil:
        # fmt:off
        workflow.connect([
            (hvoftrans, basilhv, [("output_image", "atlasfile")]),
            (hvoftrans, pvchv, [("output_image", "atlasfile")]),
            (sc207trans, pvc207, [("output_image", "atlasfile")]),
            (sc207trans, basil207, [("output_image", "atlasfile")]),
            (sc217trans, pvc217, [("output_image", "atlasfile")]),
            (sc217trans, basil217, [("output_image", "atlasfile")]),
            (sc407trans, pvc407, [("output_image", "atlasfile")]),
            (sc407trans, basil407, [("output_image", "atlasfile")]),
            (sc417trans, pvc417, [("output_image", "atlasfile")]),
            (sc417trans, basil417, [("output_image", "atlasfile")]),
            (inputnode, basilhv, [("basil", "in_cbf")]),
            (inputnode, pvchv, [("pvc", "in_cbf")]),
            (inputnode, basil207, [("basil", "in_cbf")]),
            (inputnode, pvc207, [("pvc", "in_cbf")]),
            (inputnode, basil217, [("basil", "in_cbf")]),
            (inputnode, pvc217, [("pvc", "in_cbf")]),
            (inputnode, basil407, [("basil", "in_cbf")]),
            (inputnode, pvc407, [("pvc", "in_cbf")]),
            (inputnode, basil417, [("basil", "in_cbf")]),
            (inputnode, pvc417, [("pvc", "in_cbf")]),
            (basilhv, outputnode, [("atlascsv", "basil_hvoxf")]),
            (basil207, outputnode, [("atlascsv", "basil_sc207")]),
            (basil217, outputnode, [("atlascsv", "basil_sc217")]),
            (basil407, outputnode, [("atlascsv", "basil_sc407")]),
            (basil417, outputnode, [("atlascsv", "basil_sc417")]),
            (pvchv, outputnode, [("atlascsv", "pvc_hvoxf")]),
            (pvc207, outputnode, [("atlascsv", "pvc_sc207")]),
            (pvc217, outputnode, [("atlascsv", "pvc_sc217")]),
            (pvc407, outputnode, [("atlascsv", "pvc_sc407")]),
            (pvc417, outputnode, [("atlascsv", "pvc_sc417")]),
        ])
        # fmt:on

    return workflow


def init_gecbf_compt_wf(
    name_source,
    metadata,
    mem_gb,
    M0Scale=1,
    scorescrub=False,
    basil=False,
    name="cbf_compt_wf",
):
    """Calculate CBF for GE data.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.cbf import init_gecbf_compt_wf

            wf = init_gecbf_compt_wf(
                name_source="",
                metadata={},
                mem_gb=0.1,
            )
    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The CBF was quantified from *preproccessed* ASL data using a standard
model [@detre_perfusion_1992;@alsop_recommended_2015].
"""
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "in_file",
                "asl_mask",
                "t1w_tpms",
                "t1w_mask",
                "t1_asl_xform",
                "itk_asl_to_t1",
                "m0_file",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "out_cbf",
                "out_mean",
                "out_score",
                "out_avgscore",
                "out_scrub",
                "out_cbfb",
                "out_scoreindex",
                "out_cbfpv",
                "out_att",
                "out_cbfpvwm",
            ]
        ),
        name="outputnode",
    )
    # convert tmps to asl_space
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

    # XXX: This is a bad way to find this file.
    aslcontext = name_source.replace("_asl.nii.gz", "_aslcontext.tsv")
    aslcontext_df = pd.read_table(aslcontext)
    cbf_only = all(aslcontext_df["volume_type"].isin(("m0scan", "cbf")))
    if cbf_only and not basil:
        config.loggers.workflow.info(f"Only CBF volumes are detected in {name_source}.")
    elif cbf_only:
        config.loggers.workflow.warning(
            f"Only CBF volumes are detected in {name_source}. "
            "BASIL will automatically be disabled."
        )
        basil = False

    vol_types = aslcontext_df["volume_type"].tolist()
    deltam_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "deltam"]
    cbf_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "CBF"]
    control_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "control"]

    tiscbf = get_tis(metadata)

    refine_mask = pe.Node(RefineMask(), mem_gb=1, run_without_submitting=True, name="refine_mask")

    # fmt:off
    workflow.connect([
        (inputnode, refine_mask, [
            ("t1w_mask", "in_t1mask"),
            ("asl_mask", "in_aslmask"),
            ("t1_asl_xform", "transforms"),
        ]),
    ])
    # fmt:on

    def _pick_gm(files):
        return files[0]

    def _pick_wm(files):
        return files[1]

    def _pick_csf(files):
        return files[-1]

    def _getfiledir(file):
        import os

        return os.path.dirname(file)

    collect_cbf = niu.IdentityInterface(
        fields=["deltam", "cbf"],
        name="collect_cbf",
    )

    is_casl = pcasl_or_pasl(metadata=metadata)

    if deltam_volume_idx or control_volume_idx:
        # If deltaM or label-control pairs are available, then calculate CBF.
        extract_deltam = pe.Node(
            ExtractCBForDeltaM(file_type="d"),
            mem_gb=1,
            run_without_submitting=True,
            name="extract_deltam",
        )

        # fmt:off
        workflow.connect([
            # extract deltaM data
            (inputnode, extract_deltam, [
                ("asl_file", "in_asl"),
                ("asl_mask", "in_aslmask"),
            ]),
        ])
        # fmt:on

        compute_cbf = pe.Node(
            ComputeCBF(
                metadata=metadata,
                m0scale=M0Scale,
                cbf_only=cbf_only,
            ),
            mem_gb=mem_gb,
            run_without_submitting=True,
            name="compute_cbf",
        )

        # fmt:off
        workflow.connect([
            (inputnode, compute_cbf, [("m0_file", "m0_file")]),
            (extract_deltam, compute_cbf, [("out_file", "deltam")]),
            (extract_deltam, collect_cbf, [("out_file", "deltam")]),
            (refine_mask, compute_cbf, [("out_mask", "mask")]),
            (compute_cbf, collect_cbf, [("cbf", "cbf")]),
            (compute_cbf, outputnode, [
                ("cbf", "out_cbf"),
                ("mean_cbf", "out_mean"),
            ]),
        ])
        # fmt:on

    elif cbf_volume_idx:
        extract_cbf = pe.Node(
            ExtractCBForDeltaM(file_type="c"),
            mem_gb=1,
            run_without_submitting=True,
            name="extract_cbf",
        )

        # fmt:off
        workflow.connect([
            (inputnode, extract_cbf, [
                ("asl_file", "in_asl"),
                ("asl_mask", "in_aslmask"),
            ]),
        ])
        # fmt:on

        mask_cbf = pe.Node(
            MultiImageMaths(op_string=" -mul  %s "),
            mem_gb=1,
            run_without_submitting=True,
            name="mask_cbf",
        )

        # fmt:off
        workflow.connect([
            (refine_mask, mask_cbf, [("out_mask", "in_file")]),
            (extract_cbf, mask_cbf, [("out_file", "operand_files")]),
            (mask_cbf, collect_cbf, [("out_file", "cbf")]),
            (mask_cbf, outputnode, [
                ("out_file", "out_cbf"),
                ("out_file", "out_mean"),
            ]),
        ])
        # fmt:on

    if scorescrub:
        workflow.__desc__ += """\
CBF is susceptible to artifacts due to low signal to noise ratio and high sensitivity to  motion.
Therefore, the Structural Correlation with RobUst Bayesian (SCRUB) algorithm was applied to the CBF
timeseries to discard few extreme outlier volumes (if present) and iteratively reweight the mean
CBF with structural tissues probability maps [@dolui2017structural;@dolui2016scrub].
"""

        score_and_scrub_cbf = pe.Node(
            ScoreAndScrubCBF(in_thresh=0.7, in_wfun="huber"),
            mem_gb=mem_gb,
            name="scorescrub",
            run_without_submitting=True,
        )

        # fmt:off
        workflow.connect([
            # extract probability maps
            (inputnode, csf_tfm, [
                ("asl_mask", "reference_image"),
                ("t1_asl_xform", "transforms"),
                (("t1w_tpms", _pick_csf), "input_image"),
            ]),
            (inputnode, wm_tfm, [
                ("asl_mask", "reference_image"),
                ("t1_asl_xform", "transforms"),
                (("t1w_tpms", _pick_wm), "input_image"),
            ]),
            (inputnode, gm_tfm, [
                ("asl_mask", "reference_image"),
                ("t1_asl_xform", "transforms"),
                (("t1w_tpms", _pick_gm), "input_image"),
            ]),
            (refine_mask, score_and_scrub_cbf, [("out_mask", "in_mask")]),
            (gm_tfm, score_and_scrub_cbf, [("output_image", "in_greyM")]),
            (wm_tfm, score_and_scrub_cbf, [("output_image", "in_whiteM")]),
            (csf_tfm, score_and_scrub_cbf, [("output_image", "in_csf")]),
            (collect_cbf, score_and_scrub_cbf, [("cbf", "in_file")]),
            (score_and_scrub_cbf, outputnode, [
                ("out_score", "out_score"),
                ("out_scoreindex", "out_scoreindex"),
                ("out_avgscore", "out_avgscore"),
                ("out_scrub", "out_scrub"),
            ]),
        ])
        # fmt:on

    if basil:
        workflow.__desc__ += """\
In addition, CBF was also computed by Bayesian Inference for Arterial Spin Labeling
(BASIL) as implemented in FSL. BASIL is based on Bayesian inference principles
[@chappell2008variational], and computed CBF from ASL by incorporating natural
variability of other model parameters and spatial regularization of the estimated
perfusion image, including correction of partial volume effects [@chappell_pvc].
"""

        if is_casl:
            bolus = metadata["LabelingDuration"]
        else:
            assert metadata.get("BolusCutOffFlag"), "BolusCutOffFlag must be True for PASL"
            bolus = metadata["BolusCutOffDelayTime"]
            if metadata["BolusCutOffTechnique"] == "Q2TIPS":
                # BolusCutOffDelayTime is a list, and the first entry should be used.
                bolus = bolus[0]

        basilcbf = pe.Node(
            BASILCBF(
                m0scale=M0Scale,
                bolus=bolus,
                m0tr=metadata["RepetitionTime"],
                alpha=metadata.get("LabelingEfficiency", Undefined),
                pvc=True,
                tis=tiscbf,
                pcasl=is_casl,
            ),
            name="basilcbf",
            run_without_submitting=True,
            mem_gb=mem_gb,
        )

        # fmt:off
        workflow.connect([
            (inputnode, basilcbf, [
                (("asl_mask", _getfiledir), "out_basename"),
                ("m0_file", "mzero"),
            ]),
            (collect_cbf, basilcbf, [("deltam", "in_file")]),
            (gm_tfm, basilcbf, [("output_image", "pvgm")]),
            (wm_tfm, basilcbf, [("output_image", "pvwm")]),
            (refine_mask, basilcbf, [("out_mask", "mask")]),
            (basilcbf, outputnode, [
                ("out_cbfb", "out_cbfb"),
                ("out_cbfpv", "out_cbfpv"),
                ("out_cbfpvwm", "out_cbfpvwm"),
                ("out_att", "out_att"),
            ]),
        ])
        # fmt:on

    return workflow
