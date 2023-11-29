# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for calculating confounds for ASL data."""
from fmriprep.workflows.bold.confounds import _carpet_parcellation
from nipype.algorithms import confounds as nac
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.utility import AddTSVHeader
from templateflow.api import get as get_template

from aslprep.config import DEFAULT_MEMORY_MIN_GB
from aslprep.interfaces import ASLCarpetPlot, DerivativesDataSink, GatherConfounds
from aslprep.interfaces.ants import ApplyTransforms


def init_asl_confounds_wf(
    n_volumes: int,
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

            from aslprep.tests.tests import mock_config
            from aslprep.workflows.asl.confounds import init_asl_confounds_wf

            with mock_config():
                wf = init_asl_confounds_wf(
                    n_volumes=50,
                    mem_gb=1,
                    freesurfer=True,
                )

    Parameters
    ----------
    n_volumes : :obj:`int`
        Number of volumes in the ASL file.
        Some relative measures (e.g., FD, DVARS, RMSD) will produce empty arrays for single-volume
        datasets, so we must skip those.
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
    aslref2anat_xfm
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
    from niworkflows.interfaces.images import SignalExtraction
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
                "aslref2anat_xfm",
            ],
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

    add_motion_headers = pe.Node(
        AddTSVHeader(columns=["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"]),
        name="add_motion_headers",
        mem_gb=0.01,
        run_without_submitting=True,
    )
    workflow.connect([(inputnode, add_motion_headers, [("movpar_file", "in_file")])])

    if n_volumes > 2:  # set to 2 bc relative arrays will be 1D instead of 2D for 1-volume data
        # DVARS
        dvars = pe.Node(
            nac.ComputeDVARS(save_nstd=True, save_std=True, remove_zerovariance=True),
            name="dvars",
            mem_gb=mem_gb,
        )
        workflow.connect([
            (inputnode, dvars, [
                ("asl", "in_file"),
                ("asl_mask", "in_mask"),
            ]),
        ])  # fmt:skip

        # Frame displacement
        fdisp = pe.Node(
            nac.FramewiseDisplacement(parameter_source="SPM"),
            name="fdisp",
            mem_gb=mem_gb,
        )
        workflow.connect([(inputnode, fdisp, [("movpar_file", "in_file")])])

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
        add_rmsd_header = pe.Node(
            AddTSVHeader(columns=["rmsd"]),
            name="add_rmsd_header",
            mem_gb=0.01,
            run_without_submitting=True,
        )

        workflow.connect([
            # Collate computed confounds together
            (inputnode, add_rmsd_header, [("rmsd_file", "in_file")]),
            (dvars, add_dvars_header, [("out_nstd", "in_file")]),
            (dvars, add_std_dvars_header, [("out_std", "in_file")]),
        ])  # fmt:skip

    # Project T1w mask into BOLD space and merge with BOLD brainmask
    t1w_mask_tfm = pe.Node(
        ApplyTransforms(interpolation="MultiLabel", invert_transform_flags=[True], args="-v"),
        name="t1w_mask_tfm",
    )
    union_mask = pe.Node(niu.Function(function=_binary_union), name="union_mask")

    # Create the crown mask
    dilated_mask = pe.Node(BinaryDilation(), name="dilated_mask")
    subtract_mask = pe.Node(BinarySubtraction(), name="subtract_mask")

    workflow.connect([
        # Brain mask
        (inputnode, t1w_mask_tfm, [
            ("t1w_mask", "input_image"),
            ("asl_mask", "reference_image"),
            ("aslref2anat_xfm", "transforms"),
        ]),
        (inputnode, union_mask, [("asl_mask", "mask1")]),
        (t1w_mask_tfm, union_mask, [("output_image", "mask2")]),
        (union_mask, dilated_mask, [("out", "in_mask")]),
        (union_mask, subtract_mask, [("out", "in_subtract")]),
        (dilated_mask, subtract_mask, [("out_mask", "in_base")]),
        (subtract_mask, outputnode, [("out_mask", "crown_mask")]),
    ])  # fmt:skip

    # Generate aCompCor probseg maps
    acc_masks = pe.Node(aCompCorMasks(is_aseg=freesurfer), name="acc_masks")
    workflow.connect([
        (inputnode, acc_masks, [
            ("t1w_tpms", "in_vfs"),
            (("asl", _get_zooms), "bold_zooms"),
        ]),
    ])  # fmt:skip

    # Resample probseg maps in BOLD space via T1w-to-BOLD transform
    acc_msk_tfm = pe.MapNode(
        ApplyTransforms(interpolation="Gaussian", invert_transform_flags=[True], args="-v"),
        iterfield=["input_image"],
        name="acc_msk_tfm",
        mem_gb=0.1,
    )
    workflow.connect([
        (inputnode, acc_msk_tfm, [
            ("aslref2anat_xfm", "transforms"),
            ("asl_mask", "reference_image"),
        ]),
        (acc_masks, acc_msk_tfm, [("out_masks", "input_image")]),
    ])  # fmt:skip

    acc_msk_brain = pe.MapNode(ApplyMask(), name="acc_msk_brain", iterfield=["in_file"])
    workflow.connect([
        (inputnode, acc_msk_brain, [("asl_mask", "in_mask")]),
        (acc_msk_tfm, acc_msk_brain, [("output_image", "in_file")]),
    ])  # fmt:skip

    acc_msk_bin = pe.MapNode(Binarize(thresh_low=0.99), name="acc_msk_bin", iterfield=["in_file"])
    workflow.connect([
        (acc_msk_brain, acc_msk_bin, [("out_file", "in_file")]),
        (acc_msk_bin, outputnode, [("out_file", "acompcor_masks")]),
    ])  # fmt:skip

    # Global and segment regressors
    signals_class_labels = [
        "global_signal",
        "csf",
        "white_matter",
        "csf_wm",
    ]
    merge_rois = pe.Node(
        niu.Merge(2, ravel_inputs=True),
        name="merge_rois",
        run_without_submitting=True,
    )
    signals = pe.Node(
        SignalExtraction(class_labels=signals_class_labels),
        name="signals",
        mem_gb=mem_gb,
    )
    workflow.connect([
        (inputnode, merge_rois, [("asl_mask", "in1")]),
        (acc_msk_bin, merge_rois, [("out_file", "in2")]),
        (inputnode, signals, [("asl", "in_file")]),
        (merge_rois, signals, [("out", "label_files")]),
    ])  # fmt:skip

    concat = pe.Node(
        GatherConfounds(
            fd=None,
            rmsd=None,
            dvars=None,
            std_dvars=None,
        ),
        name="concat",
        mem_gb=0.01,
        run_without_submitting=True,
    )
    workflow.connect([
        (add_motion_headers, concat, [("out_file", "motion")]),
        (signals, concat, [("out_file", "signals")]),
        (concat, outputnode, [("confounds_file", "confounds_file")]),
    ])  # fmt:skip
    if n_volumes > 2:
        workflow.connect([
            (fdisp, concat, [("out_file", "fd")]),
            (add_rmsd_header, concat, [("out_file", "rmsd")]),
            (add_dvars_header, concat, [("out_file", "dvars")]),
            (add_std_dvars_header, concat, [("out_file", "std_dvars")]),
        ])  # fmt:skip

    return workflow


def init_carpetplot_wf(
    mem_gb: float,
    confounds_list: list,
    metadata: dict,
    cifti_output: bool,
    suffix: str = "asl",
    name: str = "asl_carpet_wf",
):
    """
    Build a workflow to generate *carpet* plots.

    Resamples the MNI parcellation (ad-hoc parcellation derived from the
    Harvard-Oxford template and others).

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    confounds_list : :obj:`list` of length-3 :obj:`tuple` of :obj:`str`
    metadata : :obj:`dict`
        BIDS metadata for BOLD file
    name : :obj:`str`
        Name of workflow (default: ``asl_carpet_wf``)

    Inputs
    ------
    asl
        BOLD image, after the prescribed corrections (STC, HMC and SDC)
        when available.
    asl_mask
        BOLD series mask
    confounds_file
        TSV of all aggregated confounds
    aslref2anat_xfm
        Affine matrix that maps the T1w space into alignment with
        the native BOLD space
    std2anat_xfm
        ANTs-compatible affine-and-warp transform file
    cifti_asl
        BOLD image in CIFTI format, to be used in place of volumetric BOLD
    crown_mask
        Mask of brain edge voxels
    acompcor_mask
        Mask of deep WM+CSF
    dummy_scans
        Number of nonsteady states to be dropped at the beginning of the timeseries.

    Outputs
    -------
    out_carpetplot
        Path of the generated SVG file

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl",
                "asl_mask",
                "confounds_file",
                "aslref2anat_xfm",
                "std2anat_xfm",
                "cifti_asl",
                "crown_mask",
                "acompcor_mask",
                "dummy_scans",
            ]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(niu.IdentityInterface(fields=["out_carpetplot"]), name="outputnode")

    # Carpetplot and confounds plot
    conf_plot = pe.Node(
        ASLCarpetPlot(
            tr=metadata["RepetitionTime"],
            confounds_list=confounds_list,
        ),
        name="conf_plot",
        mem_gb=mem_gb,
    )
    ds_report_asl_conf = pe.Node(
        DerivativesDataSink(
            desc="carpetplot",
            datatype="figures",
            suffix=suffix,
            extension="svg",
            dismiss_entities=("echo",),
        ),
        name="ds_report_asl_conf",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    parcels = pe.Node(niu.Function(function=_carpet_parcellation), name="parcels")
    parcels.inputs.nifti = not cifti_output
    # List transforms
    mrg_xfms = pe.Node(niu.Merge(2), name="mrg_xfms")

    # Warp segmentation into EPI space
    resample_parc = pe.Node(
        ApplyTransforms(
            dimension=3,
            input_image=str(
                get_template(
                    "MNI152NLin2009cAsym",
                    resolution=1,
                    desc="carpet",
                    suffix="dseg",
                    extension=[".nii", ".nii.gz"],
                ),
            ),
            interpolation="MultiLabel",
            invert_transform_flags=[True, False],
            args="-u int -v",
        ),
        name="resample_parc",
    )

    workflow = Workflow(name=name)
    if cifti_output:
        workflow.connect(inputnode, "cifti_asl", conf_plot, "in_cifti")

    # fmt:off
    workflow.connect([
        (inputnode, mrg_xfms, [
            ("aslref2anat_xfm", "in1"),
            ("std2anat_xfm", "in2"),
        ]),
        (inputnode, resample_parc, [("asl_mask", "reference_image")]),
        (inputnode, parcels, [
            ("crown_mask", "crown_mask"),
            ("acompcor_mask", "acompcor_mask"),
        ]),
        (inputnode, conf_plot, [
            ("asl", "in_nifti"),
            ("confounds_file", "confounds_file"),
            ("dummy_scans", "drop_trs"),
        ]),
        (mrg_xfms, resample_parc, [("out", "transforms")]),
        (resample_parc, parcels, [("output_image", "segmentation")]),
        (parcels, conf_plot, [("out", "in_segm")]),
        (conf_plot, ds_report_asl_conf, [("out_file", "in_file")]),
        (conf_plot, outputnode, [("out_file", "out_carpetplot")]),
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
