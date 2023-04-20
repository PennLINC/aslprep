# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for registering ASL data."""
import os
import os.path as op

import pkg_resources as pkgr
from nipype.interfaces import c3, fsl
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.interfaces.registration import FLIRTRPT
from aslprep.niworkflows.utils.images import dseg_label
from aslprep.utils.misc import _conditional_downsampling

DEFAULT_MEMORY_MIN_GB = config.DEFAULT_MEMORY_MIN_GB
LOGGER = config.loggers.workflow


def init_asl_reg_wf(
    use_bbr,
    asl2t1w_dof,
    asl2t1w_init,
    sloppy=False,
    write_report=True,
    name="asl_reg_wf",
):
    """Build a workflow to run same-subject, ASL-to-T1w image-registration.

    Calculates the registration between a reference ASL image and T1w-space
    using a boundary-based registration (BBR) cost function.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.registration import init_asl_reg_wf

            wf = init_asl_reg_wf(
                use_bbr=True,
                asl2t1w_dof=9,
                asl2t1w_init="register",
            )

    Parameters
    ----------
    use_bbr : :obj:`bool` or None
        Enable/disable boundary-based registration refinement.
        If ``None``, test BBR result for distortion before accepting.
    asl2t1w_dof : 6, 9 or 12
        Degrees-of-freedom for ASL-T1w registration
    asl2t1w_init : str, 'header' or 'register'
        If ``'header'``, use header information for initialization of ASL and T1 images.
        If ``'register'``, align volumes by their centers.
    mem_gb : :obj:`float`
        Size of ASL file in GB
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    name : :obj:`str`
        Name of workflow (default: ``asl_reg_wf``)
    use_compression : :obj:`bool`
        Save registered ASL series as ``.nii.gz``
    write_report : :obj:`bool`
        Whether a reportlet should be stored

    Inputs
    ------
    ref_asl_brain
        Reference image to which ASL series is aligned
        If ``fieldwarp == True``, ``ref_asl_brain`` should be unwarped
    t1w_brain
        Skull-stripped ``t1w_preproc``
    t1w_dseg
        Segmentation of preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF)

    Outputs
    -------
    aslref_to_anat_xfm
        Affine transform from ``ref_asl_brain`` to T1 space (ITK format)
    anat_to_aslref_xfm
        Affine transform from T1 space to ASL space (ITK format)
    fallback
        Boolean indicating whether BBR was rejected (mri_coreg registration returned)

    See Also
    --------
      * :py:func:`~aslprep.workflows.asl.registration.init_fsl_bbr_wf`

    """
    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "ref_asl_brain",
                "t1w_brain",
                "t1w_dseg",
            ]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "aslref_to_anat_xfm",
                "anat_to_aslref_xfm",
                "fallback",
            ]
        ),
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


def init_fsl_bbr_wf(use_bbr, asl2t1w_dof, asl2t1w_init, sloppy=False, name="fsl_bbr_wf"):
    """Build a workflow to run FSL's ``flirt``.

    This workflow uses FSL FLIRT to register a ASL image to a T1-weighted
    structural image, using a boundary-based registration (BBR) cost function.
    It is a counterpart to :py:func:`~aslprep.workflows.asl.registration.init_bbreg_wf`,
    which performs the same task using FreeSurfer's ``bbregister``.

    The ``use_bbr`` option permits a high degree of control over registration.
    If ``False``, standard, rigid coregistration will be performed by FLIRT.
    If ``True``, FLIRT-BBR will be seeded with the initial transform found by
    the rigid coregistration.
    If ``None``, after FLIRT-BBR is run, the resulting affine transform
    will be compared to the initial transform found by FLIRT.
    Excessive deviation will result in rejecting the BBR refinement and
    accepting the original, affine registration.

    Workflow Graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.registration import init_fsl_bbr_wf

            wf = init_fsl_bbr_wf(
                use_bbr=True,
                asl2t1w_dof=9,
                asl2t1w_init="register",
            )


    Parameters
    ----------
    use_bbr : :obj:`bool` or None
        Enable/disable boundary-based registration refinement.
        If ``None``, test BBR result for distortion before accepting.
    asl2t1w_dof : 6, 9 or 12
        Degrees-of-freedom for ASL-T1w registration
    asl2t1w_init : str, 'header' or 'register'
        If ``'header'``, use header information for initialization of ASL and T1 images.
        If ``'register'``, align volumes by their centers.
    name : :obj:`str`, optional
        Workflow name (default: fsl_bbr_wf)

    Inputs
    ------
    in_file
        Reference ASL image to be registered
    t1w_brain
        Skull-stripped T1-weighted structural image
    t1w_dseg
        FAST segmentation of ``t1w_brain``


    Outputs
    -------
    aslref_to_anat_xfm
        Affine transform from ``ref_asl_brain`` to T1w space (ITK format)
    anat_to_aslref_xfm
        Affine transform from T1 space to ASL space (ITK format)
    out_report
        Reportlet for assessing registration quality
    fallback
        Boolean indicating whether BBR was rejected (rigid FLIRT registration returned)

    """
    workflow = Workflow(name=name)
    workflow.__desc__ = f"""\
ASLPrep co-registered the ASL reference to the T1w reference using *FSL*’s `flirt` [@flirt], which
implemented the boundary-based registration cost-function [@bbr]. Co-registration used
{asl2t1w_dof} degrees of freedom. The quality of co-registration and normalization to template was
quantified using the Dice and Jaccard indices, the cross-correlation with the reference image,
and the overlap between the ASL and reference images (e.g., image coverage).
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            [
                "in_file",
                "t1w_dseg",
                "t1w_brain",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            [
                "aslref_to_anat_xfm",
                "anat_to_aslref_xfm",
                "out_report",
                "fallback",
            ]
        ),
        name="outputnode",
    )

    wm_mask = pe.Node(niu.Function(function=dseg_label), name="wm_mask")
    wm_mask.inputs.label = 2  # BIDS default is WM=2
    flt_bbr_init = pe.Node(
        FLIRTRPT(dof=6, generate_report=not use_bbr, uses_qform=True),
        name="flt_bbr_init",
    )

    if asl2t1w_init not in ("register", "header"):
        raise ValueError(f"Unknown ASL-T1w initialization option: {asl2t1w_init}")

    if asl2t1w_init == "header":
        raise NotImplementedError("Header-based registration initialization not supported for FSL")

    invt_bbr = pe.Node(
        fsl.ConvertXFM(invert_xfm=True),
        name="invt_bbr",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    # ASL to T1 transform matrix is from fsl, using c3 tools to convert to
    # something ANTs will like.
    fsl2itk_fwd = pe.Node(
        c3.C3dAffineTool(fsl2ras=True, itk_transform=True),
        name="fsl2itk_fwd",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    fsl2itk_inv = pe.Node(
        c3.C3dAffineTool(fsl2ras=True, itk_transform=True),
        name="fsl2itk_inv",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([
        (inputnode, flt_bbr_init, [
            ("in_file", "in_file"),
            ("t1w_brain", "reference"),
        ]),
        (inputnode, fsl2itk_fwd, [
            ("t1w_brain", "reference_file"),
            ("in_file", "source_file"),
        ]),
        (inputnode, fsl2itk_inv, [
            ("in_file", "reference_file"),
            ("t1w_brain", "source_file"),
        ]),
        (invt_bbr, fsl2itk_inv, [("out_file", "transform_file")]),
        (fsl2itk_fwd, outputnode, [("itk_transform", "aslref_to_anat_xfm")]),
        (fsl2itk_inv, outputnode, [("itk_transform", "anat_to_aslref_xfm")]),
    ])
    # fmt:on

    # Short-circuit workflow building, use rigid registration
    if use_bbr is False:
        # fmt:off
        workflow.connect([
            (flt_bbr_init, invt_bbr, [("out_matrix_file", "in_file")]),
            (flt_bbr_init, fsl2itk_fwd, [("out_matrix_file", "transform_file")]),
            (flt_bbr_init, outputnode, [("out_report", "out_report")]),
        ])
        # fmt:on
        outputnode.inputs.fallback = True

        return workflow

    flt_bbr = pe.Node(
        FLIRTRPT(cost_func="bbr", dof=asl2t1w_dof, generate_report=True),
        name="flt_bbr",
    )

    FSLDIR = os.getenv("FSLDIR")
    if FSLDIR:
        flt_bbr.inputs.schedule = op.join(FSLDIR, "etc/flirtsch/bbr.sch")
    else:
        # Should mostly be hit while building docs
        LOGGER.warning("FSLDIR unset - using packaged BBR schedule")
        flt_bbr.inputs.schedule = pkgr.resource_filename("aslprep", "data/flirtsch/bbr.sch")

    # fmt:off
    workflow.connect([
        (inputnode, wm_mask, [("t1w_dseg", "in_seg")]),
        (inputnode, flt_bbr, [("in_file", "in_file")]),
        (flt_bbr_init, flt_bbr, [("out_matrix_file", "in_matrix_file")]),
    ])
    # fmt:on

    if sloppy is True:
        downsample = pe.Node(
            niu.Function(
                function=_conditional_downsampling,
                output_names=["out_file", "out_mask"],
            ),
            name="downsample",
        )
        # fmt:off
        workflow.connect([
            (inputnode, downsample, [("t1w_brain", "in_file")]),
            (wm_mask, downsample, [("out", "in_mask")]),
            (downsample, flt_bbr, [
                ("out_file", "reference"),
                ("out_mask", "wm_seg"),
            ]),
        ])
        # fmt:on
    else:
        # fmt:off
        workflow.connect([
            (inputnode, flt_bbr, [("t1w_brain", "reference")]),
            (wm_mask, flt_bbr, [("out", "wm_seg")]),
        ])
        # fmt:on

    # Short-circuit workflow building, use boundary-based registration
    if use_bbr is True:
        # fmt:off
        workflow.connect([
            (flt_bbr, invt_bbr, [("out_matrix_file", "in_file")]),
            (flt_bbr, fsl2itk_fwd, [("out_matrix_file", "transform_file")]),
            (flt_bbr, outputnode, [("out_report", "out_report")]),
        ])
        # fmt:on
        outputnode.inputs.fallback = False

        return workflow

    # fmt:off
    workflow.connect([
        (flt_bbr, invt_bbr, [("out_matrix_file", "in_file")]),
        (flt_bbr, fsl2itk_fwd, [("out_matrix_file", "transform_file")]),
        (flt_bbr, outputnode, [("out_report", "out_report")]),
    ])
    # fmt:on

    return workflow
