# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for calculating CBF."""
import numpy as np
from nipype.interfaces import utility as niu
from nipype.interfaces.fsl import Info
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.func.util import init_enhance_and_skullstrip_bold_wf
from niworkflows.interfaces.images import RobustAverage
from templateflow.api import get as get_template

from aslprep import config
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.bids import DerivativesDataSink
from aslprep.interfaces.cbf import (
    BASILCBF,
    ComputeCBF,
    ExtractCBF,
    RefineMask,
    ScoreAndScrubCBF,
)
from aslprep.interfaces.parcellation import ParcellateCBF
from aslprep.utils.asl import (
    determine_multi_pld,
    estimate_labeling_efficiency,
    get_bolus_duration,
    get_inflow_times,
    pcasl_or_pasl,
)
from aslprep.utils.atlas import get_atlas_names, get_atlas_nifti
from aslprep.utils.bids import find_atlas_entities


def init_cbf_wf(
    name_source,
    processing_target,
    metadata,
    dummy_scans,
    scorescrub=False,
    basil=False,
    m0_scale=1,
    smooth_kernel=5,
    name="cbf_wf",
):
    """Create a workflow for :abbr:`CCBF (compute cbf)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            import json

            from aslprep.tests.tests import mock_config
            from aslprep import config
            from aslprep.workflows.asl.cbf import init_cbf_wf

            with mock_config():
                perf_dir = config.execution.bids_dir / "sub-01" / "perf"
                with open(perf_dir / "sub-01_asl.json", "r") as fo:
                    metadata = json.load(fo)

                wf = init_cbf_wf(
                    name_source=str(perf_dir / "sub-01_asl.nii.gz"),
                    processing_target="control",
                    metadata=metadata,
                    dummy_scans=0,
                    scorescrub=True,
                    basil=True,
                )

    Parameters
    ----------
    name_source : :obj:`str`
        Path to the raw ASL file.
    metadata : :obj:`dict`
        BIDS metadata for asl file
    scorescrub
    basil
    m0_scale
    smooth_kernel
    name : :obj:`str`
        Name of workflow (default: ``cbf_wf``)

    Inputs
    ------
    asl_file
        asl series NIfTI file, after preprocessing
    aslcontext : :obj:`str`
    m0scan : :obj:`str` or None
        An M0 scan (if available as a separate file) in aslref space.
    m0scan_metadata : :obj:`dict` or None
    asl_mask
        asl mask NIFTI file
    t1w_tpms
        t1w probability maps
    t1w_mask
        t1w mask Nifti
    aslref2anat_xfm
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

"""

    m0type = metadata["M0Type"]
    is_casl = pcasl_or_pasl(metadata=metadata)
    is_multi_pld = determine_multi_pld(metadata=metadata)
    if (processing_target == "cbf") and not basil:
        config.loggers.workflow.info(f"Only CBF volumes are detected in {name_source}.")
    elif processing_target == "cbf":
        config.loggers.workflow.warning(
            f"Only CBF volumes are detected in {name_source}. "
            "BASIL will automatically be disabled."
        )
        basil = False

    if m0type in ("Included", "Separate"):
        m0_str = (
            "Calibration (M0) volumes associated with the ASL scan were smoothed with a "
            f"Gaussian kernel (FWHM={smooth_kernel} mm) and the average calibration image was "
            f"calculated and scaled by {m0_scale}."
        )
    elif m0type == "Estimate":
        m0_str = (
            f"A single M0 estimate of {metadata['M0Estimate']} was used to produce a calibration "
            f"'image' and was scaled by {m0_scale}."
        )
    else:
        m0_str = (
            f"As no calibration images or provided M0 estimate was available for the ASL scan, "
            "the control volumes used as a substitute. "
            "The control volumes in the ASL scans were smoothed with a "
            f"Gaussian kernel (FWHM={smooth_kernel} mm) and the average control image was "
            f"calculated and scaled by {m0_scale}."
        )

    if processing_target == "cbf":
        workflow.__desc__ += """\
*ASLPrep* loaded pre-calculated cerebral blood flow (CBF) data from the ASL file.
"""

    elif is_casl:
        if is_multi_pld:
            workflow.__desc__ += f"""\
*ASLPrep* calculated cerebral blood flow (CBF) from the multi-delay
{metadata['ArterialSpinLabelingType']} data using the following method.

First, delta-M values were averaged over time for each post-labeling delay (PLD).
{m0_str}

Next, arterial transit time (ATT) was estimated on a voxel-wise basis according to
@dai2012reduced.

CBF was then calculated for each delay using the mean delta-M values and the estimated ATT,
according to the formula from @fan2017long.

CBF was then averaged over delays according to @juttukonda2021characterizing,
in which an unweighted average is calculated for each voxel across all delays in which
PLD + labeling duration > ATT.
"""
        else:
            # Single-delay (P)CASL data
            workflow.__desc__ += f"""\
*ASLPrep* calculated cerebral blood flow (CBF) from the single-delay
{metadata['ArterialSpinLabelingType']} using a single-compartment general kinetic model
[@buxton1998general].
{m0_str}
"""

    else:
        bcut = metadata.get("BolusCutOffTechnique")

        if is_multi_pld:
            workflow.__desc__ += """HAHAHAHA!"""

        # Single-delay PASL data, with different bolus cut-off techniques
        if bcut == "QUIPSS":
            workflow.__desc__ += f"""\
*ASLPrep* calculated cerebral blood flow (CBF) from the single-delay PASL
using a single-compartment general kinetic model [@buxton1998general]
using the QUIPSS modification, as described in @wong1998quantitative.
{m0_str}
"""
        elif bcut == "QUIPSSII":
            workflow.__desc__ += f"""\
*ASLPrep* calculated cerebral blood flow (CBF) from the single-delay PASL
using a single-compartment general kinetic model [@buxton1998general]
using the QUIPSS II modification, as described in @alsop_recommended_2015.
{m0_str}
"""
        elif bcut == "Q2TIPS":
            workflow.__desc__ += f"""\
*ASLPrep* calculated cerebral blood flow (CBF) from the single-delay PASL
using a single-compartment general kinetic model [@buxton1998general]
using the Q2TIPS modification, as described in @noguchi2015technical.
{m0_str}
"""
        else:
            # No bolus cutoff delay technique
            raise ValueError("PASL without a bolus cut-off technique is not supported in ASLPrep.")

    if "SliceTiming" in metadata:
        workflow.__desc__ += (
            "Prior to calculating CBF, post-labeling delay values were shifted on a slice-wise "
            "basis based on the slice timing."
        )

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "aslcontext",
                "m0scan",
                "m0scan_metadata",
                "asl_mask",
                "t1w_tpms",
                "t1w_mask",
                "aslref2anat_xfm",
            ],
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "mean_cbf",
                "cbf_ts",  # Only calculated for single-delay data
                "att",  # Only calculated for multi-delay data
                "abat",  # Only calculated for multi-delay data
                "abv",  # Only calculated for multi-delay data
                "plds",
                # SCORE/SCRUB outputs
                "cbf_ts_score",
                "mean_cbf_score",
                "mean_cbf_scrub",
                "score_outlier_index",
                # BASIL outputs
                "mean_cbf_basil",
                "mean_cbf_gm_basil",
                "mean_cbf_wm_basil",
                "att_basil",
            ]
        ),
        name="outputnode",
    )

    warp_t1w_mask_to_asl = pe.Node(
        ApplyTransforms(
            dimension=3,
            float=True,
            interpolation="NearestNeighbor",
            invert_transform_flags=[True],
            input_image_type=3,
            args="-v",
        ),
        name="warp_t1w_mask_to_asl",
    )
    workflow.connect([
        (inputnode, warp_t1w_mask_to_asl, [
            ("asl_mask", "reference_image"),
            ("t1w_mask", "input_image"),
            ("aslref2anat_xfm", "transforms"),
        ]),
    ])  # fmt:skip

    reduce_mask = pe.Node(
        RefineMask(),
        mem_gb=0.2,
        run_without_submitting=True,
        name="reduce_mask",
    )

    workflow.connect([
        (inputnode, reduce_mask, [("asl_mask", "asl_mask")]),
        (warp_t1w_mask_to_asl, reduce_mask, [("output_image", "t1w_mask")]),
    ])  # fmt:skip

    # Warp tissue probability maps to ASL space
    def _pick_gm(files):
        if not isinstance(files, list):
            raise ValueError(f"_pick_gm: input is not list ({files})")
        return files[0]

    def _pick_wm(files):
        if not isinstance(files, list):
            raise ValueError(f"_pick_wm: input is not list ({files})")
        return files[1]

    def _pick_csf(files):
        if not isinstance(files, list):
            raise ValueError(f"_pick_csf: input is not list ({files})")
        return files[2]

    def _getfiledir(file):
        import os

        return os.path.dirname(file)

    gm_tfm = pe.Node(
        ApplyTransforms(
            interpolation="NearestNeighbor",
            float=True,
            invert_transform_flags=[True],
            args="-v",
        ),
        name="gm_tfm",
        mem_gb=0.1,
    )

    workflow.connect([
        (inputnode, gm_tfm, [
            ("asl_mask", "reference_image"),
            ("aslref2anat_xfm", "transforms"),
            (("t1w_tpms", _pick_gm), "input_image"),
        ]),
    ])  # fmt:skip

    wm_tfm = pe.Node(
        ApplyTransforms(
            interpolation="NearestNeighbor",
            float=True,
            invert_transform_flags=[True],
            args="-v",
        ),
        name="wm_tfm",
        mem_gb=0.1,
    )

    workflow.connect([
        (inputnode, wm_tfm, [
            ("asl_mask", "reference_image"),
            ("aslref2anat_xfm", "transforms"),
            (("t1w_tpms", _pick_wm), "input_image"),
        ]),
    ])  # fmt:skip

    csf_tfm = pe.Node(
        ApplyTransforms(
            interpolation="NearestNeighbor",
            float=True,
            invert_transform_flags=[True],
            args="-v",
        ),
        name="csf_tfm",
        mem_gb=0.1,
    )

    workflow.connect([
        (inputnode, csf_tfm, [
            ("asl_mask", "reference_image"),
            ("aslref2anat_xfm", "transforms"),
            (("t1w_tpms", _pick_csf), "input_image"),
        ]),
    ])  # fmt:skip

    extract_deltam = pe.Node(
        ExtractCBF(
            name_source=name_source,
            dummy_scans=dummy_scans,
            fwhm=smooth_kernel,
            metadata=metadata,
        ),
        mem_gb=0.2,
        run_without_submitting=True,
        name="extract_deltam",
    )

    workflow.connect([
        (inputnode, extract_deltam, [
            ("asl_file", "asl_file"),
            ("aslcontext", "aslcontext"),
            ("m0scan_metadata", "m0scan_metadata"),
        ]),
        (reduce_mask, extract_deltam, [("out_mask", "in_mask")]),
    ])  # fmt:skip

    if metadata["M0Type"] == "Separate":
        mean_m0 = pe.Node(RobustAverage(), name="mean_m0", mem_gb=1)
        workflow.connect([
            (inputnode, mean_m0, [("m0scan", "in_file")]),
            (mean_m0, extract_deltam, [("out_file", "m0scan")]),
        ])  # fmt:skip

        enhance_and_skullstrip_m0scan_wf = init_enhance_and_skullstrip_bold_wf(
            pre_mask=False,
            omp_nthreads=1,
            name="enhance_and_skullstrip_m0scan_wf",
        )
        workflow.connect([
            (mean_m0, enhance_and_skullstrip_m0scan_wf, [("out_file", "inputnode.in_file")]),
            (enhance_and_skullstrip_m0scan_wf, reduce_mask, [("outputnode.mask_file", "m0_mask")]),
        ])  # fmt:skip

    compute_cbf = pe.Node(
        ComputeCBF(
            cbf_only=processing_target == "cbf",
            m0_scale=m0_scale,
        ),
        mem_gb=0.2,
        run_without_submitting=True,
        name="compute_cbf",
    )

    workflow.connect([
        (reduce_mask, compute_cbf, [("out_mask", "mask")]),
        (extract_deltam, compute_cbf, [
            ("out_file", "deltam"),
            ("m0_file", "m0_file"),
            ("metadata", "metadata"),
        ]),
        (compute_cbf, outputnode, [
            ("cbf_ts", "cbf_ts"),
            ("mean_cbf", "mean_cbf"),
            ("att", "att"),
            ("plds", "plds"),
            ("abat", "abat"),
            ("abv", "abv"),
        ]),
    ])  # fmt:skip

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
        workflow.connect([
            (reduce_mask, score_and_scrub_cbf, [("out_mask", "mask")]),
            (compute_cbf, score_and_scrub_cbf, [("cbf_ts", "cbf_ts")]),
            (gm_tfm, score_and_scrub_cbf, [("output_image", "gm_tpm")]),
            (wm_tfm, score_and_scrub_cbf, [("output_image", "wm_tpm")]),
            (csf_tfm, score_and_scrub_cbf, [("output_image", "csf_tpm")]),
            (score_and_scrub_cbf, outputnode, [
                ("cbf_ts_score", "cbf_ts_score"),
                ("score_outlier_index", "score_outlier_index"),
                ("mean_cbf_score", "mean_cbf_score"),
                ("mean_cbf_scrub", "mean_cbf_scrub"),
            ]),
        ])  # fmt:skip

    if basil:
        workflow.__desc__ += f"""

CBF was also computed with Bayesian Inference for Arterial Spin Labeling (BASIL)
[@chappell2008variational], as implemented in *FSL* {Info.version()}.
BASIL computes CBF using a spatial regularization of the estimated perfusion image and
additionally calculates a partial-volume corrected CBF image [@chappell_pvc].
"""
        # Node to define bolus
        determine_bolus_duration = pe.Node(
            niu.Function(
                function=get_bolus_duration,
                input_names=["metadata", "is_casl"],
                output_names=["bolus"],
            ),
            name="determine_bolus_duration",
        )
        determine_bolus_duration.inputs.is_casl = is_casl
        workflow.connect([(extract_deltam, determine_bolus_duration, [("metadata", "metadata")])])

        # Node to define tis
        determine_inflow_times = pe.Node(
            niu.Function(
                function=get_inflow_times,
                input_names=["metadata", "is_casl"],
                output_names=["tis"],
            ),
            name="determine_inflow_times",
        )
        determine_inflow_times.inputs.is_casl = is_casl

        workflow.connect([(extract_deltam, determine_inflow_times, [("metadata", "metadata")])])

        # Node to estimate labeling efficiency
        estimate_alpha = pe.Node(
            niu.Function(
                function=estimate_labeling_efficiency,
                input_names=["metadata"],
                output_names=["labeling_efficiency"],
            ),
            name="estimate_alpha",
        )

        workflow.connect([(extract_deltam, estimate_alpha, [("metadata", "metadata")])])

        basil_kwargs = {}
        if "SliceTiming" in metadata.keys():
            slicetime_diffs = np.unique(np.diff(metadata["SliceTiming"]))
            # Check if slice times are monotonic
            monotonic_slicetimes = slicetime_diffs.size == 1
            # Check if slice times are ascending
            ascending_slicetimes = np.all(slicetime_diffs > 0)
            # Only set slicedt for ascending slice orders.
            if monotonic_slicetimes and ascending_slicetimes:
                basil_kwargs["slice_spacing"] = slicetime_diffs[0]
            else:
                config.loggers.interface.warning(
                    "Slice times are not ascending. They will be ignored in the BASIL call."
                )

        basilcbf = pe.Node(
            BASILCBF(
                m0_scale=m0_scale,
                pvc=True,
                pcasl=is_casl,
                **basil_kwargs,
            ),
            name="basilcbf",
            run_without_submitting=True,
            mem_gb=0.2,
        )

        workflow.connect([
            (reduce_mask, basilcbf, [("out_mask", "mask")]),
            (extract_deltam, basilcbf, [
                (("m0_file", _getfiledir), "out_basename"),
                ("out_file", "deltam"),
                ("m0_file", "mzero"),
            ]),
            (determine_bolus_duration, basilcbf, [("bolus", "bolus")]),
            (determine_inflow_times, basilcbf, [("tis", "tis")]),
            (estimate_alpha, basilcbf, [("labeling_efficiency", "alpha")]),
            (gm_tfm, basilcbf, [("output_image", "gm_tpm")]),
            (wm_tfm, basilcbf, [("output_image", "wm_tpm")]),
            (basilcbf, outputnode, [
                ("mean_cbf_basil", "mean_cbf_basil"),
                ("mean_cbf_gm_basil", "mean_cbf_gm_basil"),
                ("mean_cbf_wm_basil", "mean_cbf_wm_basil"),
                ("att_basil", "att_basil"),
            ]),
        ])  # fmt:skip

        if metadata["M0Type"] != "Estimate":
            workflow.connect([(extract_deltam, basilcbf, [("m0tr", "m0tr")])])

    return workflow


def init_parcellate_cbf_wf(
    cbf_3d,
    min_coverage=0.5,
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

            wf = init_parcellate_cbf_wf(cbf_3d=["mean_cbf", "mean_cbf_score"])

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
    source_file : str
    mean_cbf : str
    mean_cbf_score : Undefined or str
    mean_cbf_scrub : Undefined or str
    mean_cbf_basil : Undefined or str
    mean_cbf_gm_basil : Undefined or str
    asl_mask : str
    aslref2anat_xfm : str
    MNI152NLin2009cAsym_to_anat_xfm : str
        The transform from MNI152NLin2009cAsym to the subject's anatomical space.

    Outputs
    -------
    atlas_names : list of str
        A list of atlases used for parcellating the CBF results.
        The list of atlas names is generated by :func:`aslprep.utils.atlas.get_atlas_names`.
        The atlases include: "4S156Parcels", "4S256Parcels", "4S356Parcels", "4S456Parcels",
        "4S556Parcels", "4S656Parcels", "4S756Parcels", "4S856Parcels", "4S956Parcels",
        "4S1056Parcels", "Glasser", "Gordon", "Tian", and "HCP".
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
    CBF_ENTITIES = {
        "mean_cbf": {},
        "mean_cbf_score": {
            "desc": "score",
        },
        "mean_cbf_scrub": {
            "desc": "scrub",
        },
        "mean_cbf_basil": {
            "desc": "basil",
        },
        "mean_cbf_gm_basil": {
            "desc": "basilGM",
        },
        "mean_cbf_wm_basil": {
            "desc": "basilWM",
        },
    }
    workflow = Workflow(name=name)

    workflow.__desc__ = f"""
Parcellated CBF estimates were extracted for the following atlases:
the Schaefer Supplemented with Subcortical Structures (4S) atlas
[@Schaefer_2017; @pauli2018high; @king2019functional; @najdenovska2018vivo; @glasser2013minimal] at
10 different resolutions (156, 256, 356, 456, 556, 656, 756, 856, 956, and 1056 parcels),
the Glasser atlas [@Glasser_2016], the Gordon atlas [@Gordon_2014],
the Tian subcortical atlas [@tian2020topographic], and the HCP CIFTI subcortical atlas
[@glasser2013minimal].
In cases of partial coverage, either uncovered voxels (values of all zeros or NaNs) were
ignored (when the parcel had >{min_coverage * 100}% coverage)
or the whole parcel was set to zero (when the parcel had <{min_coverage * 100}% coverage).
"""

    input_fields = [
        "source_file",
        "asl_mask",
        "aslref2anat_xfm",
        "MNI152NLin2009cAsym_to_anat_xfm",
    ]
    input_fields += cbf_3d
    inputnode = pe.Node(
        niu.IdentityInterface(fields=input_fields),
        name="inputnode",
    )
    output_fields = ["atlas_names"] + [f"{field}_parcellated" for field in cbf_3d]
    outputnode = pe.Node(
        niu.IdentityInterface(fields=output_fields),
        name="outputnode",
    )

    atlas_name_grabber = pe.Node(
        niu.Function(
            input_names=["subset"],
            output_names=["atlas_names"],
            function=get_atlas_names,
        ),
        name="atlas_name_grabber",
    )
    atlas_name_grabber.inputs.subset = "all"
    workflow.connect([(atlas_name_grabber, outputnode, [("atlas_names", "atlas_names")])])

    # get atlases via aslprep.data.load
    atlas_file_grabber = pe.MapNode(
        niu.Function(
            input_names=["atlas_name"],
            output_names=["atlas_file", "atlas_labels_file", "atlas_metadata_file"],
            function=get_atlas_nifti,
        ),
        name="atlas_file_grabber",
        iterfield=["atlas_name"],
    )
    workflow.connect([(atlas_name_grabber, atlas_file_grabber, [("atlas_names", "atlas_name")])])

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

    workflow.connect([
        (inputnode, merge_xforms, [
            ("MNI152NLin2009cAsym_to_anat_xfm", "in2"),
            ("aslref2anat_xfm", "in3"),
        ]),
    ])  # fmt:skip

    # Using the generated transforms, apply them to get everything in the correct MNI form
    warp_atlases_to_asl_space = pe.MapNode(
        ApplyTransforms(
            interpolation="GenericLabel",
            input_image_type=3,
            dimension=3,
            invert_transform_flags=[False, False, True],
            args="-v",
        ),
        name="warp_atlases_to_asl_space",
        iterfield=["input_image"],
        mem_gb=mem_gb,
        n_procs=omp_nthreads,
    )

    workflow.connect([
        (inputnode, warp_atlases_to_asl_space, [("asl_mask", "reference_image")]),
        (atlas_file_grabber, warp_atlases_to_asl_space, [("atlas_file", "input_image")]),
        (merge_xforms, warp_atlases_to_asl_space, [("out", "transforms")]),
    ])  # fmt:skip

    for cbf_type in cbf_3d:
        parcellate_cbf = pe.MapNode(
            ParcellateCBF(min_coverage=min_coverage),
            name=f"parcellate_{cbf_type}",
            iterfield=["atlas", "atlas_labels"],
            mem_gb=mem_gb,
        )

        workflow.connect([
            (inputnode, parcellate_cbf, [
                (cbf_type, "in_file"),
                ("asl_mask", "mask"),
            ]),
            (atlas_file_grabber, parcellate_cbf, [("atlas_labels_file", "atlas_labels")]),
            (warp_atlases_to_asl_space, parcellate_cbf, [("output_image", "atlas")]),
        ])  # fmt:skip

        ds_cbf = pe.MapNode(
            DerivativesDataSink(
                base_directory=config.execution.aslprep_dir,
                check_hdr=False,
                suffix="cbf",
                **CBF_ENTITIES[cbf_type],
            ),
            name=f"ds_{cbf_type}",
            iterfield=["atlas", "in_file"],
            run_without_submitting=True,
        )
        workflow.connect([
            (inputnode, ds_cbf, [("source_file", "source_file")]),
            (atlas_name_grabber, ds_cbf, [("atlas_names", "atlas")]),
            (parcellate_cbf, ds_cbf, [("timeseries", "in_file")]),
        ])  # fmt:skip

        if cbf_type in ("mean_cbf", "mean_cbf_basil"):
            # I think it is easier to only retain the coverage file for the regular CBF estimates.
            # SCORE/SCRUB CBF should have the same coverage as the regular CBF.
            # BASIL might have different coverage because it drops any voxels with negative CBF,
            # but the different PVC types will have the same coverage as the main BASIL CBF.
            ds_coverage = pe.MapNode(
                DerivativesDataSink(
                    base_directory=config.execution.aslprep_dir,
                    check_hdr=False,
                    suffix="coverage",
                    **CBF_ENTITIES[cbf_type],
                ),
                name=f"ds_coverage_{cbf_type}",
                iterfield=["atlas", "in_file"],
                run_without_submitting=True,
            )

            workflow.connect([
                (inputnode, ds_coverage, [("source_file", "source_file")]),
                (atlas_name_grabber, ds_coverage, [("atlas_names", "atlas")]),
                (parcellate_cbf, ds_coverage, [("coverage", "in_file")]),
            ])  # fmt:skip

    # Get entities from atlas for datasinks
    get_atlas_entities = pe.MapNode(
        niu.Function(
            input_names=["filename"],
            output_names=["tpl", "atlas", "res", "suffix", "extension"],
            function=find_atlas_entities,
        ),
        name="get_atlas_entities",
        iterfield=["filename"],
    )
    workflow.connect([(atlas_file_grabber, get_atlas_entities, [("atlas_file", "filename")])])

    # Write out standard-space atlas file.
    # This won't be in the same space that the data were parcellated in,
    # but it's useful as a reference.
    ds_atlas = pe.MapNode(
        DerivativesDataSink(
            base_directory=config.execution.aslprep_dir,
            check_hdr=False,
            dismiss_entities=["datatype", "subject", "session", "task", "run", "desc"],
            allowed_entities=["space", "res", "den", "atlas", "desc", "cohort"],
        ),
        name="ds_atlas",
        iterfield=["space", "atlas", "resolution", "suffix", "extension", "in_file"],
        run_without_submitting=True,
    )

    workflow.connect([
        (inputnode, ds_atlas, [("source_file", "source_file")]),
        (atlas_file_grabber, ds_atlas, [("atlas_file", "in_file")]),
        (get_atlas_entities, ds_atlas, [
            ("tpl", "space"),
            ("atlas", "atlas"),
            ("res", "resolution"),
            ("suffix", "suffix"),
            ("extension", "extension"),
        ]),
    ])  # fmt:skip

    ds_atlas_labels_file = pe.MapNode(
        DerivativesDataSink(
            base_directory=config.execution.aslprep_dir,
            check_hdr=False,
            dismiss_entities=[
                "datatype",
                "subject",
                "session",
                "task",
                "run",
                "desc",
                "space",
                "res",
                "den",
                "cohort",
            ],
            allowed_entities=["atlas"],
            extension=".tsv",
        ),
        name="ds_atlas_labels_file",
        iterfield=["atlas", "suffix", "in_file"],
        run_without_submitting=True,
    )

    workflow.connect([
        (inputnode, ds_atlas_labels_file, [("source_file", "source_file")]),
        (atlas_file_grabber, ds_atlas_labels_file, [("atlas_labels_file", "in_file")]),
        (get_atlas_entities, ds_atlas_labels_file, [
            ("atlas", "atlas"),
            ("suffix", "suffix"),
        ]),
    ])  # fmt:skip

    ds_atlas_metadata = pe.MapNode(
        DerivativesDataSink(
            base_directory=config.execution.aslprep_dir,
            check_hdr=False,
            dismiss_entities=[
                "datatype",
                "subject",
                "session",
                "task",
                "run",
                "desc",
                "space",
                "res",
                "den",
                "cohort",
            ],
            allowed_entities=["atlas"],
            extension=".json",
        ),
        name="ds_atlas_metadata",
        iterfield=["atlas", "suffix", "in_file"],
        run_without_submitting=True,
    )

    workflow.connect([
        (inputnode, ds_atlas_metadata, [("source_file", "source_file")]),
        (atlas_file_grabber, ds_atlas_metadata, [("atlas_metadata_file", "in_file")]),
        (get_atlas_entities, ds_atlas_metadata, [
            ("atlas", "atlas"),
            ("suffix", "suffix"),
        ]),
    ])  # fmt:skip

    return workflow
