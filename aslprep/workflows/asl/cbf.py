# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for calculating CBF."""
import pandas as pd
from nipype.interfaces import utility as niu
from nipype.interfaces.fsl import Info, MultiImageMaths
from nipype.pipeline import engine as pe
from templateflow.api import get as get_template

from aslprep import config
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.cbf_computation import (
    BASILCBF,
    ComputeCBF,
    ExtractCBF,
    ExtractCBForDeltaM,
    RefineMask,
    ScoreAndScrubCBF,
)
from aslprep.interfaces.parcellation import ParcellateCBF
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.utils.atlas import get_atlas_names, get_atlas_nifti
from aslprep.utils.misc import estimate_labeling_efficiency, get_tis, pcasl_or_pasl


def init_compute_cbf_wf(
    name_source,
    aslcontext,
    metadata,
    dummy_vols,
    scorescrub=False,
    basil=False,
    m0_scale=1,
    smooth_kernel=5,
    name="compute_cbf_wf",
):
    """Create a workflow for :abbr:`CCBF (compute cbf)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.cbf import init_compute_cbf_wf

            wf = init_compute_cbf_wf(
                name_source="",
                metadata={},
                dummy_vols=0,
            )

    Parameters
    ----------
    name_source : :obj:`str`
        Path to the raw ASL file.
    aslcontext : :obj:`str`
        Path to the aslcontext file associated with the ASL file being processed.
        Used to set the aslcontext input.
    metadata : :obj:`dict`
        BIDS metadata for asl file
    scorescrub
    basil
    m0_scale
    smooth_kernel
    name : :obj:`str`
        Name of workflow (default: ``compute_cbf_wf``)

    Inputs
    ------
    asl_file
        asl series NIfTI file, after preprocessing
    aslcontext : :obj:`str`
        Defined from the parameter.
    m0scan : :obj:`str` or None
    m0scan_metadata : :obj:`dict` or None
    asl_mask
        asl mask NIFTI file
    t1w_tpms
        t1w probability maps
    t1w_mask
        t1w mask Nifti
    anat_to_aslref_xfm
        t1w to asl transformation file
    aslref_to_anat_xfm
        asl to t1w transformation file

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

    if m0_scale != 1:
        workflow.__desc__ += (
            f"Prior to calculating CBF, the M0 volumes were scaled by a factor of {m0_scale}."
        )

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "m0scan",
                "m0scan_metadata",
                "asl_mask",
                "t1w_tpms",
                "t1w_mask",
                "anat_to_aslref_xfm",
                "aslref_to_anat_xfm",
            ]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
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
                "mean_cbf_wm_basil",
                "att",
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
            ("t1w_mask", "t1w_mask"),
            ("asl_mask", "asl_mask"),
            ("anat_to_aslref_xfm", "transforms"),
        ]),
    ])
    # fmt:on

    # Warp tissue probability maps to ASL space
    def _pick_gm(files):
        return files[0]

    def _pick_wm(files):
        return files[1]

    def _pick_csf(files):
        return files[2]

    def _getfiledir(file):
        import os

        return os.path.dirname(file)

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

    tiscbf = get_tis(metadata)
    is_casl = pcasl_or_pasl(metadata=metadata)

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
            name_source=name_source,
            aslcontext=aslcontext,
            dummy_vols=dummy_vols,
            fwhm=smooth_kernel,
            metadata=metadata,
        ),
        mem_gb=0.2,
        run_without_submitting=True,
        name="extract_deltam",
    )

    # fmt:off
    workflow.connect([
        (inputnode, extract_deltam, [
            ("asl_file", "asl_file"),
            ("m0scan", "m0scan"),
            ("m0scan_metadata", "m0scan_metadata"),
        ]),
        (refine_mask, extract_deltam, [("out_mask", "in_mask")]),
    ])
    # fmt:on

    compute_cbf = pe.Node(
        ComputeCBF(
            cbf_only=cbf_only,
            m0_scale=m0_scale,
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
            ("cbf", "cbf_ts"),
            ("mean_cbf", "mean_cbf"),
        ]),
    ])
    # fmt:on

    if scorescrub:
        score_and_scrub_cbf = pe.Node(
            ScoreAndScrubCBF(tpm_threshold=0.7, wavelet_function="huber"),
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
            (refine_mask, score_and_scrub_cbf, [("out_mask", "mask")]),
            (compute_cbf, score_and_scrub_cbf, [("cbf", "cbf_ts")]),
            (gm_tfm, score_and_scrub_cbf, [("output_image", "gm_tpm")]),
            (wm_tfm, score_and_scrub_cbf, [("output_image", "wm_tpm")]),
            (csf_tfm, score_and_scrub_cbf, [("output_image", "csf_tpm")]),
            (score_and_scrub_cbf, outputnode, [
                ("cbf_ts_score", "cbf_ts_score"),
                ("score_outlier_index", "score_outlier_index"),
                ("mean_cbf_score", "mean_cbf_score"),
                ("mean_cbf_scrub", "mean_cbf_scrub"),
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
            if metadata["BolusCutOffTechnique"] == "Q2TIPS":
                # BolusCutOffDelayTime is a list, and the first entry should be used.
                bolus = bolus[0]

        basilcbf = pe.Node(
            BASILCBF(
                m0_scale=m0_scale,
                bolus=bolus,
                alpha=estimate_labeling_efficiency(metadata),
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
                ("out_file", "deltam"),
                ("m0_file", "mzero"),
                ("m0tr", "m0tr"),
            ]),
            (gm_tfm, basilcbf, [("output_image", "gm_tpm")]),
            (wm_tfm, basilcbf, [("output_image", "wm_tpm")]),
            (basilcbf, outputnode, [
                ("mean_cbf_basil", "mean_cbf_basil"),
                ("mean_cbf_gm_basil", "mean_cbf_gm_basil"),
                ("mean_cbf_wm_basil", "mean_cbf_wm_basil"),
                ("att", "att"),
            ]),
        ])
        # fmt:on

    return workflow


def init_compute_cbf_ge_wf(
    name_source,
    aslcontext,
    metadata,
    mem_gb,
    m0_scale=1,
    scorescrub=False,
    basil=False,
    name="compute_cbf_wf",
):
    """Calculate CBF for GE data.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.cbf import init_compute_cbf_ge_wf

            wf = init_compute_cbf_ge_wf(
                name_source="",
                metadata={},
                mem_gb=0.1,
            )
    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The CBF was quantified from *preprocessed* ASL data using a standard
model [@detre_perfusion_1992;@alsop_recommended_2015].
"""

    if m0_scale != 1:
        workflow.__desc__ += (
            f"Prior to calculating CBF, the M0 volumes were scaled by a factor of {m0_scale}."
        )

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "asl_mask",
                "t1w_tpms",
                "t1w_mask",
                "anat_to_aslref_xfm",
                "aslref_to_anat_xfm",
                "m0_file",
                "m0tr",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
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
                "mean_cbf_wm_basil",
                "att",
            ]
        ),
        name="outputnode",
    )

    def _pick_gm(files):
        return files[0]

    def _pick_wm(files):
        return files[1]

    def _pick_csf(files):
        return files[2]

    def _getfiledir(file):
        import os

        return os.path.dirname(file)

    # convert tmps to asl_space
    # extract probability maps
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
    cbf_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "cbf"]
    control_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "control"]

    tiscbf = get_tis(metadata)

    refine_mask = pe.Node(RefineMask(), mem_gb=1, run_without_submitting=True, name="refine_mask")

    # fmt:off
    workflow.connect([
        (inputnode, refine_mask, [
            ("t1w_mask", "t1w_mask"),
            ("asl_mask", "asl_mask"),
            ("anat_to_aslref_xfm", "transforms"),
        ]),
    ])
    # fmt:on

    collect_cbf = pe.Node(
        niu.IdentityInterface(
            fields=["deltam", "cbf"],
        ),
        name="collect_cbf",
    )

    is_casl = pcasl_or_pasl(metadata=metadata)

    if deltam_volume_idx or control_volume_idx:
        # If deltaM or label-control pairs are available, then calculate CBF.
        extract_deltam = pe.Node(
            ExtractCBForDeltaM(file_type="d", aslcontext=aslcontext),
            mem_gb=1,
            run_without_submitting=True,
            name="extract_deltam",
        )

        # fmt:off
        workflow.connect([
            # extract deltaM data
            (inputnode, extract_deltam, [
                ("asl_file", "asl_file"),
                ("asl_mask", "asl_mask"),
            ]),
        ])
        # fmt:on

        compute_cbf = pe.Node(
            ComputeCBF(
                metadata=metadata,
                m0_scale=m0_scale,
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
                ("cbf", "cbf_ts"),
                ("mean_cbf", "mean_cbf"),
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
                ("asl_file", "asl_file"),
                ("asl_mask", "asl_mask"),
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
                ("out_file", "cbf_ts"),
                ("out_file", "mean_cbf"),
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
            ScoreAndScrubCBF(tpm_threshold=0.7, wavelet_function="huber"),
            mem_gb=mem_gb,
            name="score_and_scrub_cbf",
            run_without_submitting=True,
        )

        # fmt:off
        workflow.connect([
            (refine_mask, score_and_scrub_cbf, [("out_mask", "mask")]),
            (gm_tfm, score_and_scrub_cbf, [("output_image", "gm_tpm")]),
            (wm_tfm, score_and_scrub_cbf, [("output_image", "wm_tpm")]),
            (csf_tfm, score_and_scrub_cbf, [("output_image", "csf_tpm")]),
            (collect_cbf, score_and_scrub_cbf, [("cbf", "cbf_ts")]),
            (score_and_scrub_cbf, outputnode, [
                ("cbf_ts_score", "cbf_ts_score"),
                ("score_outlier_index", "score_outlier_index"),
                ("mean_cbf_score", "mean_cbf_score"),
                ("mean_cbf_scrub", "mean_cbf_scrub"),
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
                m0_scale=m0_scale,
                bolus=bolus,
                alpha=estimate_labeling_efficiency(metadata),
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
                ("m0tr", "m0tr"),
            ]),
            (collect_cbf, basilcbf, [("deltam", "deltam")]),
            (gm_tfm, basilcbf, [("output_image", "gm_tpm")]),
            (wm_tfm, basilcbf, [("output_image", "wm_tpm")]),
            (refine_mask, basilcbf, [("out_mask", "mask")]),
            (basilcbf, outputnode, [
                ("mean_cbf_basil", "mean_cbf_basil"),
                ("mean_cbf_gm_basil", "mean_cbf_gm_basil"),
                ("mean_cbf_wm_basil", "mean_cbf_wm_basil"),
                ("att", "att"),
            ]),
        ])
        # fmt:on

    return workflow


def init_parcellate_cbf_wf(
    min_coverage=0.5,
    scorescrub=False,
    basil=False,
    mem_gb=0.1,
    omp_nthreads=1,
    name="parcellate_cbf_wf",
):
    """Parcellate CBF results using a set of atlases.

    Atlases are in MNI152NLin6Asym space, but ASLPrep is not guaranteed to generate
    MNI152NLin6Asym-space derivatives.
    However, it *is* guaranteed to generate MNI152NLin2009cAsym derivatives,
    so this workflow warps the atlases to asl reference space via the path:
    MNI152NLin6Asym --> MNI152NLin2009cAsym --> anat --> asl.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.cbf import init_parcellate_cbf_wf

            wf = init_parcellate_cbf_wf()

    Parameters
    ----------
    min_coverage
    scorescrub
    basil
    mem_gb
    omp_nthreads
    name
        Default is "parcellate_cbf_wf".

    Inputs
    ------
    mean_cbf : str
    mean_cbf_score : Undefined or str
    mean_cbf_scrub : Undefined or str
    mean_cbf_basil : Undefined or str
    mean_cbf_gm_basil : Undefined or str
    asl_mask : str
    anat_to_aslref_xfm : str
    MNI152NLin2009cAsym_to_anat_xfm : str
        The transform from MNI152NLin2009cAsym to the subject's anatomical space.

    Outputs
    -------
    atlas_names : list of str
        A list of atlases used for parcellating the CBF results.
        The list of atlas names is generated by :func:`aslprep.utils.atlas.get_atlas_names`.
        The atlases include: "Schaefer117", "Schaefer217", "Schaefer317", "Schaefer417",
        "Schaefer517", "Schaefer617", "Schaefer717", "Schaefer817", "Schaefer917",
        "Schaefer1017", "Glasser", "Gordon", and "subcortical" (Tian).
    mean_cbf_parcellated : list of str
    mean_cbf_score_parcellated : Undefined or list of str
        Only defined if ``scorescrub`` is True.
    mean_cbf_scrub_parcellated : Undefined or list of str
        Only defined if ``scorescrub`` is True.
    mean_cbf_basil_parcellated : Undefined or list of str
        Only defined if ``basil`` is True.
    mean_cbf_gm_basil_parcellated : Undefined or list of str
        Only defined if ``basil`` is True.
    """
    workflow = Workflow(name=name)

    workflow.__desc__ = """\
For each CBF map, the ROIs for the following atlases were extracted:
the Harvard-Oxford and the Schaefer 200 and 400-parcel resolution atlases.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "mean_cbf",
                "mean_cbf_score",
                "mean_cbf_scrub",
                "mean_cbf_basil",
                "mean_cbf_gm_basil",
                "asl_mask",
                "anat_to_aslref_xfm",
                "MNI152NLin2009cAsym_to_anat_xfm",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "atlas_names",
                "mean_cbf_parcellated",
                "mean_cbf_score_parcellated",
                "mean_cbf_scrub_parcellated",
                "mean_cbf_basil_parcellated",
                "mean_cbf_gm_basil_parcellated",
            ],
        ),
        name="outputnode",
    )

    atlas_name_grabber = pe.Node(
        niu.Function(output_names=["atlas_names"], function=get_atlas_names),
        name="atlas_name_grabber",
    )

    # fmt:off
    workflow.connect([(atlas_name_grabber, outputnode, [("atlas_names", "atlas_names")])])
    # fmt:on

    # get atlases via pkgrf
    atlas_file_grabber = pe.MapNode(
        niu.Function(
            input_names=["atlas_name"],
            output_names=["atlas_file", "atlas_labels_file"],
            function=get_atlas_nifti,
        ),
        name="atlas_file_grabber",
        iterfield=["atlas_name"],
    )

    # fmt:off
    workflow.connect([(atlas_name_grabber, atlas_file_grabber, [("atlas_names", "atlas_name")])])
    # fmt:on

    # Atlases are in MNI152NLin6Asym
    MNI152NLin6Asym_to_MNI152NLin2009cAsym = str(
        get_template(
            template="MNI152NLin2009cAsym",
            mode="image",
            suffix="xfm",
            extension=".h5",
            **{"from": "MNI152NLin6Asym"},
        ),
    )

    # Prepare to warp atlases from MNI152NLin6Asym to ASL reference space.
    # One of the output spaces selected by the user *may* be MNI152NLin6Asym,
    # but MNI152NLin2009cAsym is always used, so it's safer to go:
    # MNI152NLin6Asym --> MNI152NLin2009cAsym --> anat --> asl
    merge_xforms = pe.Node(niu.Merge(3), name="merge_xforms")
    merge_xforms.inputs.in1 = MNI152NLin6Asym_to_MNI152NLin2009cAsym

    # fmt:off
    workflow.connect([
        (inputnode, merge_xforms, [
            ("MNI152NLin2009cAsym_to_anat_xfm", "in2"),
            ("anat_to_aslref_xfm", "in3"),
        ]),
    ])
    # fmt:on

    # Using the generated transforms, apply them to get everything in the correct MNI form
    warp_atlases_to_asl_space = pe.MapNode(
        ApplyTransforms(
            interpolation="GenericLabel",
            input_image_type=3,
            dimension=3,
        ),
        name="warp_atlases_to_asl_space",
        iterfield=["input_image"],
        mem_gb=mem_gb,
        n_procs=omp_nthreads,
    )

    # fmt:off
    workflow.connect([
        (inputnode, warp_atlases_to_asl_space, [("asl_mask", "reference_image")]),
        (atlas_file_grabber, warp_atlases_to_asl_space, [("atlas_file", "input_image")]),
        (merge_xforms, warp_atlases_to_asl_space, [("out", "transforms")]),
    ])
    # fmt:on

    cbf_types = ["mean_cbf"]
    if scorescrub:
        cbf_types += ["mean_cbf_score", "mean_cbf_scrub"]

    if basil:
        cbf_types += ["mean_cbf_basil", "mean_cbf_gm_basil"]

    for cbf_type in cbf_types:
        parcellate_cbf = pe.MapNode(
            ParcellateCBF(min_coverage=min_coverage),
            name=f"parcellate_{cbf_type}",
            iterfield=["atlas", "atlas_labels"],
            mem_gb=mem_gb,
        )

        # fmt:off
        workflow.connect([
            (inputnode, parcellate_cbf, [
                (cbf_type, "in_file"),
                ("asl_mask", "mask"),
            ]),
            (atlas_file_grabber, parcellate_cbf, [("atlas_labels_file", "atlas_labels")]),
            (warp_atlases_to_asl_space, parcellate_cbf, [("output_image", "atlas")]),
            (parcellate_cbf, outputnode, [
                ("timeseries", f"{cbf_type}_parcellated"),
                ("coverage", f"{cbf_type}_coverage"),
            ]),
        ])
        # fmt:on

    return workflow
