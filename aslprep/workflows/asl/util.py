# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Utility workflows."""
from nipype.interfaces import afni, fsl
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.fixes import (
    FixN4BiasFieldCorrection as N4BiasFieldCorrection,
)
from niworkflows.interfaces.header import CopyXForm

DEFAULT_MEMORY_MIN_GB = 0.01


def init_enhance_and_skullstrip_asl_wf(pre_mask=False, name="enhance_and_skullstrip_asl_wf"):
    """Enhance and run brain extraction on an ASL image.

    This workflow takes in an :abbr:`ASL (Arterial Spin Labeling)` average/summary
    (e.g., a reference image averaging non-steady-state timepoints),
    and sharpens the histogram with the application of the N4 algorithm for removing the
    :abbr:`INU (intensity non-uniformity)` bias field and calculates a signal mask.

    Steps of this workflow are:

        1.  Calculate a tentative mask by registering (9-parameters) to *fMRIPrep*'s
            :abbr:`EPI (echo-planar imaging)` -*aslref* template, which
            is in MNI space.
            The tentative mask is obtained by resampling the MNI template's
            brainmask into *aslref*-space.
        2.  Binary dilation of the tentative mask with a sphere of 3mm diameter.
        3.  Run ANTs' ``N4BiasFieldCorrection`` on the input
            :abbr:`ASL (arterial spin labeling)` average, using the
            mask generated in 1) instead of the internal Otsu thresholding.
        4.  Calculate a loose mask using FSL's ``bet``, with one mathematical morphology
            dilation of one iteration and a sphere of 6mm as structuring element.
        5.  Mask the :abbr:`INU (intensity non-uniformity)`-corrected image
            with the latest mask calculated in 3), then use AFNI's ``3dUnifize``
            to *standardize* the T2* contrast distribution.
        6.  Calculate a mask using AFNI's ``3dAutomask`` after the contrast
            enhancement of 4).
        7.  Calculate a final mask as the intersection of 4) and 6).
        8.  Apply final mask on the enhanced reference.

    Step 1 can be skipped if the ``pre_mask`` argument is set to ``True`` and
    a tentative mask is passed in to the workflow throught the ``pre_mask``
    Nipype input.

    Workflow graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.util import init_enhance_and_skullstrip_asl_wf

            wf = init_enhance_and_skullstrip_asl_wf()

    .. _N4BiasFieldCorrection: https://hdl.handle.net/10380/3053

    Parameters
    ----------
    pre_mask : bool
        Indicates whether the ``pre_mask`` input will be set (and thus, step 1 should be skipped).
    name : str
        Name of workflow (default: ``enhance_and_skullstrip_asl_wf``)

    Inputs
    ------
    in_file : str
        ASL image (single volume)
    pre_mask : bool
        A tentative brain mask to initialize the workflow (requires ``pre_mask``
        parameter set ``True``).

    Outputs
    -------
    bias_corrected_file : str
        the ``in_file`` after `N4BiasFieldCorrection`_
    skull_stripped_file : str
        the ``bias_corrected_file`` after skull-stripping
    mask_file : str
        mask of the skull-stripped input file
    out_report : str
        reportlet for the skull-stripping
    """
    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=["in_file", "pre_mask"]), name="inputnode")
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["mask_file", "skull_stripped_file", "bias_corrected_file"]),
        name="outputnode",
    )

    pre_mask = pre_mask

    # Ensure mask's header matches reference's
    # check_hdr = pe.Node(MatchHeader(), name="check_hdr", run_without_submitting=True)

    # Run N4 normally, force num_threads=1 for stability (images are small, no need for >1)
    n4_correct = pe.Node(
        N4BiasFieldCorrection(dimension=3, copy_header=True, bspline_fitting_distance=200),
        shrink_factor=2,
        name="n4_correct",
        n_procs=1,
    )
    n4_correct.inputs.rescale_intensities = True

    # Create a generous BET mask out of the bias-corrected EPI
    skullstrip_first_pass = pe.Node(fsl.BET(frac=0.2, mask=True), name="skullstrip_first_pass")
    bet_dilate = pe.Node(
        fsl.DilateImage(
            operation="max",
            kernel_shape="sphere",
            kernel_size=6.0,
            internal_datatype="char",
        ),
        name="skullstrip_first_dilate",
    )
    bet_mask = pe.Node(fsl.ApplyMask(), name="skullstrip_first_mask")

    # Use AFNI's unifize for T2 constrast & fix header
    unifize = pe.Node(
        afni.Unifize(
            t2=True,
            outputtype="NIFTI_GZ",
            # Default -clfrac is 0.1, 0.4 was too conservative
            # -rbt because I'm a Jedi AFNI Master (see 3dUnifize's documentation)
            args="-clfrac 0.2 -rbt 18.3 65.0 90.0",
            out_file="uni.nii.gz",
        ),
        name="unifize",
    )
    fixhdr_unifize = pe.Node(CopyXForm(), name="fixhdr_unifize", mem_gb=0.1)

    # Run AFNI's 3dAutomask to extract a refined brain mask
    skullstrip_second_pass = pe.Node(
        afni.Automask(dilate=1, outputtype="NIFTI_GZ"),
        name="skullstrip_second_pass",
    )
    fixhdr_skullstrip2 = pe.Node(CopyXForm(), name="fixhdr_skullstrip2", mem_gb=0.1)

    # Take intersection of both masks
    combine_masks = pe.Node(fsl.BinaryMaths(operation="mul"), name="combine_masks")

    # Compute masked brain
    apply_mask = pe.Node(fsl.ApplyMask(), name="apply_mask")

    # binarize_mask = pe.Node(Binarize(thresh_low=brainmask_thresh), name="binarize_mask")

    # fmt: off
    workflow.connect([
        (inputnode, n4_correct, [("in_file", "mask_image")]),
        (inputnode, n4_correct, [("in_file", "input_image")]),
        (inputnode, fixhdr_unifize, [("in_file", "hdr_file")]),
        (inputnode, fixhdr_skullstrip2, [("in_file", "hdr_file")]),
        (n4_correct, skullstrip_first_pass, [("output_image", "in_file")]),
        (skullstrip_first_pass, bet_dilate, [("mask_file", "in_file")]),
        (bet_dilate, bet_mask, [("out_file", "mask_file")]),
        (skullstrip_first_pass, bet_mask, [("out_file", "in_file")]),
        (bet_mask, unifize, [("out_file", "in_file")]),
        (unifize, fixhdr_unifize, [("out_file", "in_file")]),
        (fixhdr_unifize, skullstrip_second_pass, [("out_file", "in_file")]),
        (skullstrip_first_pass, combine_masks, [("mask_file", "in_file")]),
        (skullstrip_second_pass, fixhdr_skullstrip2, [("out_file", "in_file")]),
        (fixhdr_skullstrip2, combine_masks, [("out_file", "operand_file")]),
        (fixhdr_unifize, apply_mask, [("out_file", "in_file")]),
        (combine_masks, apply_mask, [("out_file", "mask_file")]),
        (combine_masks, outputnode, [("out_file", "mask_file")]),
        (apply_mask, outputnode, [("out_file", "skull_stripped_file")]),
        (n4_correct, outputnode, [("output_image", "bias_corrected_file")]),
    ])
    # fmt: on

    return workflow
