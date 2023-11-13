# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for calculating confounds for ASL data."""
from nipype.algorithms import confounds as nac
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.utility import AddTSVHeader
from templateflow.api import get as get_template

from aslprep.config import DEFAULT_MEMORY_MIN_GB
from aslprep.interfaces import ASLSummary, DerivativesDataSink, GatherConfounds
from aslprep.interfaces.ants import ApplyTransforms


def init_asl_confounds_wf(
    mem_gb: float,
    freesurfer: bool = False,
    name: str = "asl_confounds_wf",
):
    """Build a workflow to generate and write out confounding signals.

    This workflow calculates confounds for a asl series, and aggregates them
    into a :abbr:`TSV (tab-separated value)` file, for use as nuisance
    regressors in a :abbr:`GLM (general linear model)`.
    The following confounds are calculated, with column headings in parentheses:

    #. DVARS - original and standardized variants (``dvars``, ``std_dvars``)
    #. Framewise displacement, based on head-motion parameters
       (``framewise_displacement``)
    #. Estimated head-motion parameters, in mm and rad
       (``trans_x``, ``trans_y``, ``trans_z``, ``rot_x``, ``rot_y``, ``rot_z``)

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.confounds import init_asl_confounds_wf

            wf = init_asl_confounds_wf(
                mem_gb=1,
            )

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of asl file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    freesurfer : :obj:`bool`
        True if FreeSurfer derivatives were used.
    name : :obj:`str`
        Name of workflow (default: ``asl_confounds_wf``)

    Inputs
    ------
    asl
        asl image, after the prescribed corrections (HMC and SDC)
        when available.
    asl_mask
        asl series mask
    movpar_file
        SPM-formatted motion parameters file
    rmsd_file
        Framewise displacement as measured by ``fsl_motion_outliers``.
    skip_vols
        number of non steady state volumes
    t1w_mask
        Mask of the skull-stripped template image
    t1w_tpms
        List of tissue probability maps in T1w space
    anat_to_aslref_xfm
        Affine matrix that maps the T1w space into alignment with
        the native asl space

    Outputs
    -------
    confounds_file
        TSV of all aggregated confounds
    confounds_metadata
        Confounds metadata dictionary.
    crown_mask
        Mask of brain edge voxels
    acompcor_masks
    """
    from fmriprep.interfaces.confounds import aCompCorMasks
    from niworkflows.interfaces.morphology import BinaryDilation, BinarySubtraction
    from niworkflows.interfaces.nibabel import ApplyMask, Binarize

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
Several confounding timeseries were calculated, including both framewise displacement
(FD) and DVARS. FD and DVARS are calculated using the implementations in in *Nipype*
(following the definition by [@power_fd_dvars]) for each ASL run.  ASLPrep summarizes
in-scanner motion as the mean framewise displacement and relative root-mean square displacement.
"""
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl",
                "asl_mask",
                "movpar_file",
                "rmsd_file",
                "skip_vols",
                "t1w_mask",
                "t1w_tpms",
                "anat_to_aslref_xfm",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "confounds_file",
                "confounds_metadata",
                "acompcor_masks",
                "crown_mask",
            ],
        ),
        name="outputnode",
    )

    # DVARS
    dvars = pe.Node(
        nac.ComputeDVARS(save_nstd=True, save_std=True, remove_zerovariance=True),
        name="dvars",
        mem_gb=mem_gb,
    )
    # fmt:off
    workflow.connect([
        (inputnode, dvars, [
            ("asl", "in_file"),
            ("asl_mask", "in_mask"),
        ]),
    ])
    # fmt:on

    # Frame displacement
    fdisp = pe.Node(nac.FramewiseDisplacement(parameter_source="SPM"), name="fdisp", mem_gb=mem_gb)
    workflow.connect([(inputnode, fdisp, [("movpar_file", "in_file")])])

    # Global and segment regressors
    # signals_class_labels = ["csf", "white_matter", "global_signal"]

    # Arrange confounds
    add_dvars_header = pe.Node(
        AddTSVHeader(columns=["dvars"]),
        name="add_dvars_header",
        mem_gb=0.01,
        run_without_submitting=True,
    )
    add_std_dvars_header = pe.Node(
        AddTSVHeader(columns=["std_dvars"]),
        name="add_std_dvars_header",
        mem_gb=0.01,
        run_without_submitting=True,
    )
    add_motion_headers = pe.Node(
        AddTSVHeader(columns=["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"]),
        name="add_motion_headers",
        mem_gb=0.01,
        run_without_submitting=True,
    )
    add_rmsd_header = pe.Node(
        AddTSVHeader(columns=["rmsd"]),
        name="add_rmsd_header",
        mem_gb=0.01,
        run_without_submitting=True,
    )
    concat = pe.Node(GatherConfounds(), name="concat", mem_gb=0.01, run_without_submitting=True)

    # fmt:off
    workflow.connect([
        # Collate computed confounds together
        (inputnode, add_motion_headers, [("movpar_file", "in_file")]),
        (inputnode, add_rmsd_header, [("rmsd_file", "in_file")]),
        (dvars, add_dvars_header, [("out_nstd", "in_file")]),
        (dvars, add_std_dvars_header, [("out_std", "in_file")]),
        (fdisp, concat, [("out_file", "fd")]),
        (add_motion_headers, concat, [("out_file", "motion")]),
        (add_rmsd_header, concat, [("out_file", "rmsd")]),
        (add_dvars_header, concat, [("out_file", "dvars")]),
        (add_std_dvars_header, concat, [("out_file", "std_dvars")]),
        # Set outputs
        (concat, outputnode, [("confounds_file", "confounds_file")]),
    ])
    # fmt:on

    # Project T1w mask into BOLD space and merge with BOLD brainmask
    t1w_mask_tfm = pe.Node(
        ApplyTransforms(interpolation="MultiLabel"),
        name="t1w_mask_tfm",
    )
    union_mask = pe.Node(niu.Function(function=_binary_union), name="union_mask")

    # Create the crown mask
    dilated_mask = pe.Node(BinaryDilation(), name="dilated_mask")
    subtract_mask = pe.Node(BinarySubtraction(), name="subtract_mask")

    # fmt:off
    workflow.connect([
        # Brain mask
        (inputnode, t1w_mask_tfm, [
            ("t1w_mask", "input_image"),
            ("bold_mask", "reference_image"),
            ("t1_bold_xform", "transforms"),
        ]),
        (inputnode, union_mask, [("bold_mask", "mask1")]),
        (t1w_mask_tfm, union_mask, [("output_image", "mask2")]),
        (union_mask, dilated_mask, [("out", "in_mask")]),
        (union_mask, subtract_mask, [("out", "in_subtract")]),
        (dilated_mask, subtract_mask, [("out_mask", "in_base")]),
        (subtract_mask, outputnode, [("out_mask", "crown_mask")]),
    ])
    # fmt:on

    # Generate aCompCor probseg maps
    acc_masks = pe.Node(aCompCorMasks(is_aseg=freesurfer), name="acc_masks")
    # fmt:off
    workflow.connect([
        (inputnode, acc_masks, [
            ("t1w_tpms", "in_vfs"),
            (("bold", _get_zooms), "bold_zooms"),
        ]),
    ])
    # fmt:on

    # Resample probseg maps in BOLD space via T1w-to-BOLD transform
    acc_msk_tfm = pe.MapNode(
        ApplyTransforms(interpolation="Gaussian"),
        iterfield=["input_image"],
        name="acc_msk_tfm",
        mem_gb=0.1,
    )
    # fmt:off
    workflow.connect([
        (inputnode, acc_msk_tfm, [
            ("t1_bold_xform", "transforms"),
            ("bold_mask", "reference_image"),
        ]),
        (acc_masks, acc_msk_tfm, [("out_masks", "input_image")]),
    ])
    # fmt:on

    acc_msk_brain = pe.MapNode(ApplyMask(), name="acc_msk_brain", iterfield=["in_file"])
    # fmt:off
    workflow.connect([
        (inputnode, acc_msk_brain, [("asl_mask", "in_mask")]),
        (acc_msk_tfm, acc_msk_brain, [("output_image", "in_file")]),
    ])
    # fmt:on

    acc_msk_bin = pe.MapNode(Binarize(thresh_low=0.99), name="acc_msk_bin", iterfield=["in_file"])
    # fmt:off
    workflow.connect([
        (acc_msk_brain, acc_msk_bin, [("out_file", "in_file")]),
        (acc_msk_bin, outputnode, [("out_file", "acompcor_masks")]),
    ])
    # fmt:on

    return workflow


def init_carpetplot_wf(mem_gb, metadata, name="carpetplot_wf"):
    """Build a workflow to generate carpet plots.

    Resamples the MNI parcellation (ad-hoc parcellation derived from the
    Harvard-Oxford template and others).

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of ASL file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    metadata : :obj:`dict`
        BIDS metadata for ASL file
    name : :obj:`str`
        Name of workflow (default: ``carpetplot_wf``)

    Inputs
    ------
    asl
        asl image, after the prescribed corrections (HMC and SDC)
        when available.
    asl_mask
        ASL series mask
    confounds_file
        TSV of all aggregated confounds
    anat_to_aslref_xfm
        Affine matrix that maps the T1w space into alignment with
        the native ASL space
    std2anat_xfm
        ANTs-compatible affine-and-warp transform file
    """
    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl",
                "asl_mask",
                "confounds_file",
                "anat_to_aslref_xfm",
                "std2anat_xfm",
            ]
        ),
        name="inputnode",
    )

    # List transforms
    mrg_xfms = pe.Node(niu.Merge(2), name="mrg_xfms")

    # Warp segmentation into EPI space
    resample_parc = pe.Node(
        ApplyTransforms(
            float=True,
            input_image=str(
                get_template(
                    "MNI152NLin2009cAsym",
                    resolution=1,
                    desc="carpet",
                    suffix="dseg",
                    extension=[".nii", ".nii.gz"],
                )
            ),
            dimension=3,
            default_value=0,
            interpolation="MultiLabel",
        ),
        name="resample_parc",
    )

    # Carpetplot and confounds plot
    conf_plot = pe.Node(
        ASLSummary(
            tr=metadata.get("RepetitionTime", metadata["RepetitionTimePreparation"]),
            confounds_list=[("std_dvars", None, "DVARS"), ("framewise_displacement", "mm", "FD")],
        ),
        name="conf_plot",
        mem_gb=mem_gb,
    )
    ds_report_asl_conf = pe.Node(
        DerivativesDataSink(desc="carpetplot", datatype="figures", keep_dtype=True),
        name="ds_report_asl_conf",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([
        (inputnode, mrg_xfms, [
            ("anat_to_aslref_xfm", "in1"),
            ("std2anat_xfm", "in2"),
        ]),
        (inputnode, resample_parc, [("asl_mask", "reference_image")]),
        (mrg_xfms, resample_parc, [("out", "transforms")]),
        # Carpetplot
        (inputnode, conf_plot, [
            ("asl", "in_func"),
            ("asl_mask", "in_mask"),
            ("confounds_file", "confounds_file"),
        ]),
        (resample_parc, conf_plot, [("output_image", "in_segm")]),
        (conf_plot, ds_report_asl_conf, [("out_file", "in_file")]),
    ])
    # fmt:on

    return workflow


def _binary_union(mask1, mask2):
    """Generate the union of two masks."""
    from pathlib import Path

    import nibabel as nb
    import numpy as np

    img = nb.load(mask1)
    mskarr1 = np.asanyarray(img.dataobj, dtype=int) > 0
    mskarr2 = np.asanyarray(nb.load(mask2).dataobj, dtype=int) > 0
    out = img.__class__(mskarr1 | mskarr2, img.affine, img.header)
    out.set_data_dtype("uint8")
    out_name = Path("mask_union.nii.gz").absolute()
    out.to_filename(out_name)
    return str(out_name)


def _get_zooms(in_file):
    import nibabel as nb

    return tuple(nb.load(in_file).header.get_zooms()[:3])
