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
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.func.util import init_asl_reference_wf
from aslprep.niworkflows.interfaces.itk import MultiApplyTransforms
from aslprep.niworkflows.interfaces.nilearn import Merge
from aslprep.niworkflows.interfaces.registration import FLIRTRPT
from aslprep.niworkflows.interfaces.utils import GenerateSamplingReference
from aslprep.niworkflows.utils.images import dseg_label
from aslprep.utils.misc import _conditional_downsampling

DEFAULT_MEMORY_MIN_GB = config.DEFAULT_MEMORY_MIN_GB
LOGGER = config.loggers.workflow


def init_asl_reg_wf(
    use_bbr,
    asl2t1w_dof,
    asl2t1w_init,
    name="asl_reg_wf",
    sloppy=False,
    write_report=True,
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
                asl2t1w_init='register',
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
    use_fieldwarp : :obj:`bool`
        Include SDC warp in single-shot transform from ASL to T1
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
    aslref_to_t1w_xfm
        Affine transform from ``ref_asl_brain`` to T1 space (ITK format)
    itk_t1_to_asl
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
                "aslref_to_t1w_xfm",
                "itk_t1_to_asl",
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
            ("outputnode.aslref_to_t1w_xfm", "aslref_to_t1w_xfm"),
            ("outputnode.itk_t1_to_asl", "itk_t1_to_asl"),
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


def init_asl_t1_trans_wf(
    mem_gb,
    omp_nthreads,
    scorescrub=False,
    basil=False,
    cbft1space=False,
    multiecho=False,
    use_fieldwarp=False,
    use_compression=True,
    name="asl_t1_trans_wf",
):
    """
    Co-register the reference ASL image to T1w-space.

    The workflow uses :abbr:`BBR (boundary-based registration)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.registration import init_asl_t1_trans_wf

            wf = init_asl_t1_trans_wf(
                mem_gb=3,
                omp_nthreads=1,
            )

    Parameters
    ----------
    use_fieldwarp : :obj:`bool`
        Include SDC warp in single-shot transform from ASLto T1
    multiecho : :obj:`bool`
        If multiecho data was supplied, HMC already performed
    mem_gb : :obj:`float`
        Size of ASL file in GB
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    use_compression : :obj:`bool`
        Save registered ASL series as ``.nii.gz``
    name : :obj:`str`
        Name of workflow (default: ``asl_reg_wf``)

    Inputs
    ------
    name_source
        ASL series NIfTI file
        Used to recover original information lost during processing
    ref_asl_brain
        Reference image to which ASL series is aligned
        If ``fieldwarp == True``, ``ref_asl_brain`` should be unwarped
    ref_asl_mask
        Skull-stripping mask of reference image
    t1w_brain
        Skull-stripped bias-corrected structural template image
    t1w_mask
        Mask of the skull-stripped template image
    asl_split
        Individual 3D ASL volumes, not motion corrected
    hmc_xforms
        List of affine transforms aligning each volume to ``ref_image`` in ITK format
    aslref_to_t1w_xfm
        Affine transform from ``ref_asl_brain`` to T1 space (ITK format)
    fieldwarp
        a :abbr:`DFM (displacements field map)` in ITK format

    Outputs
    -------
    asl_t1
        Motion-corrected ASL series in T1 space
    aslref_t1
        Reference, contrast-enhanced summary of the motion-corrected ASL series in T1w space
    asl_mask_t1
        ASL mask in T1 space

    See also
    --------
      * :py:func:`~aslprep.workflows.asl.registration.init_fsl_bbr_wf`

    """
    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "name_source",
                "ref_asl_brain",
                "ref_asl_mask",
                "t1w_brain",
                "t1w_mask",
                "asl_split",
                "fieldwarp",
                "hmc_xforms",
                "aslref_to_t1w_xfm",
                "cbf",
                "mean_cbf",
                # SCORE/SCRUB outputs
                "cbf_ts_score",
                "mean_cbf_score",
                "mean_cbf_scrub",
                # BASIL outputs
                "mean_cbf_basil",
                "mean_cbf_gm_basil",
                "mean_cbf_wm_basil",
                "att",
            ]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_t1",
                "aslref_t1",
                "asl_mask_t1",
                "cbf_ts_t1",
                "mean_cbf_t1",
                # SCORE/SCRUB outputs
                "cbf_ts_score_t1",
                "mean_cbf_score_t1",
                "mean_cbf_scrub_t1",
                # BASIL outputs
                "mean_cbf_basil_t1",
                "mean_cbf_gm_basil_t1",
                "mean_cbf_gm_basil_t1",
                "att_t1",
            ]
        ),
        name="outputnode",
    )

    gen_ref = pe.Node(
        GenerateSamplingReference(),
        name="gen_ref",
        mem_gb=0.3,
    )  # 256x256x256 * 64 / 8 ~ 150MB

    mask_t1w_tfm = pe.Node(
        ApplyTransforms(interpolation="MultiLabel"),
        name="mask_t1w_tfm",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, gen_ref, [
            ("ref_asl_brain", "moving_image"),
            ("t1w_brain", "fixed_image"),
            ("t1w_mask", "fov_mask"),
        ]),
        (inputnode, mask_t1w_tfm, [("ref_asl_mask", "input_image")]),
        (gen_ref, mask_t1w_tfm, [("out_file", "reference_image")]),
        (inputnode, mask_t1w_tfm, [("aslref_to_t1w_xfm", "transforms")]),
        (mask_t1w_tfm, outputnode, [("output_image", "asl_mask_t1")]),
    ])
    # fmt:on

    asl_to_t1w_transform = pe.Node(
        MultiApplyTransforms(interpolation="LanczosWindowedSinc", float=True, copy_dtype=True),
        name="asl_to_t1w_transform",
        mem_gb=mem_gb * 3 * omp_nthreads,
        n_procs=omp_nthreads,
    )
    # merge 3D volumes into 4D timeseries
    merge = pe.Node(Merge(compress=use_compression), name="merge", mem_gb=mem_gb)

    # Generate a reference on the target T1w space
    gen_final_ref = init_asl_reference_wf(omp_nthreads, pre_mask=True)

    if not multiecho:
        # Merge transforms placing the head motion correction last
        nforms = 2 + int(use_fieldwarp)
        merge_xforms = pe.Node(
            niu.Merge(nforms),
            name="merge_xforms",
            run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB,
        )
        if use_fieldwarp:
            # fmt:off
            workflow.connect([(inputnode, merge_xforms, [("fieldwarp", "in2")])])
            # fmt:on

        # fmt:off
        workflow.connect([
            # merge transforms
            (inputnode, merge_xforms, [
                ("hmc_xforms", f"in{nforms}"),
                ("aslref_to_t1w_xfm", "in1"),
            ]),
            (merge_xforms, asl_to_t1w_transform, [("out", "transforms")]),
            (inputnode, asl_to_t1w_transform, [("asl_split", "input_image")]),
            (inputnode, merge, [("name_source", "header_source")]),
            (gen_ref, asl_to_t1w_transform, [("out_file", "reference_image")]),
            (asl_to_t1w_transform, merge, [("out_files", "in_files")]),
            (merge, gen_final_ref, [("out_file", "inputnode.asl_file")]),
            (mask_t1w_tfm, gen_final_ref, [("output_image", "inputnode.asl_mask")]),
            (merge, outputnode, [("out_file", "asl_t1")]),
        ])
        # fmt:on

    else:
        from nipype.interfaces.fsl import Split as FSLSplit

        asl_split = pe.Node(
            FSLSplit(dimension="t"),
            name="asl_split",
            mem_gb=DEFAULT_MEMORY_MIN_GB,
        )

        # fmt:off
        workflow.connect([
            (inputnode, asl_split, [("asl_split", "in_file")]),
            (asl_split, asl_to_t1w_transform, [("out_files", "input_image")]),
            (inputnode, asl_to_t1w_transform, [("aslref_to_t1w_xfm", "transforms")]),
            (inputnode, merge, [("name_source", "header_source")]),
            (gen_ref, asl_to_t1w_transform, [("out_file", "reference_image")]),
            (asl_to_t1w_transform, merge, [("out_files", "in_files")]),
            (merge, gen_final_ref, [("out_file", "inputnode.asl_file")]),
            (mask_t1w_tfm, gen_final_ref, [("output_image", "inputnode.asl_mask")]),
            (merge, outputnode, [("out_file", "asl_t1")]),
        ])
        # fmt:on

    if cbft1space:
        cbf_to_t1w_transform = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
                dimension=3,
            ),
            name="cbf_to_t1w_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )
        meancbf_to_t1w_transform = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
            ),
            name="meancbf_to_t1w_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )

        # fmt:off
        workflow.connect([
            (gen_final_ref, outputnode, [("outputnode.ref_image", "aslref_t1")]),
            (inputnode, cbf_to_t1w_transform, [("cbf", "input_image")]),
            (cbf_to_t1w_transform, outputnode, [("output_image", "cbf_ts_t1")]),
            (inputnode, cbf_to_t1w_transform, [("aslref_to_t1w_xfm", "transforms")]),
            (gen_ref, cbf_to_t1w_transform, [("out_file", "reference_image")]),
            (inputnode, meancbf_to_t1w_transform, [("mean_cbf", "input_image")]),
            (meancbf_to_t1w_transform, outputnode, [("output_image", "mean_cbf_t1")]),
            (inputnode, meancbf_to_t1w_transform, [("aslref_to_t1w_xfm", "transforms")]),
            (gen_ref, meancbf_to_t1w_transform, [("out_file", "reference_image")]),
        ])
        # fmt:on

    if cbft1space and scorescrub:
        score_to_t1w_transform = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
                dimension=3,
            ),
            name="score_to_t1w_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )
        avgscore_to_t1w_transform = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
            ),
            name="avgscore_to_t1w_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )
        scrub_to_t1w_transform = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
            ),
            name="scrub_to_t1w_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )

        # fmt:off
        workflow.connect([
            (inputnode, score_to_t1w_transform, [("score", "input_image")]),
            (score_to_t1w_transform, outputnode, [("output_image", "cbf_ts_score_t1")]),
            (inputnode, score_to_t1w_transform, [("aslref_to_t1w_xfm", "transforms")]),
            (gen_ref, score_to_t1w_transform, [("out_file", "reference_image")]),
            (inputnode, avgscore_to_t1w_transform, [("mean_cbf_score", "input_image")]),
            (avgscore_to_t1w_transform, outputnode, [("output_image", "mean_cbf_score_t1")]),
            (inputnode, avgscore_to_t1w_transform, [("aslref_to_t1w_xfm", "transforms")]),
            (gen_ref, avgscore_to_t1w_transform, [("out_file", "reference_image")]),
            (inputnode, scrub_to_t1w_transform, [("scrub", "input_image")]),
            (scrub_to_t1w_transform, outputnode, [("output_image", "mean_cbf_scrub_t1")]),
            (inputnode, scrub_to_t1w_transform, [("aslref_to_t1w_xfm", "transforms")]),
            (gen_ref, scrub_to_t1w_transform, [("out_file", "reference_image")]),
        ])
        # fmt:on

    if cbft1space and basil:
        basil_to_t1w_transform = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
            ),
            name="basil_to_t1w_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )
        pv_to_t1w_transform = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
            ),
            name="pv_to_t1w_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )
        pvwm_to_t1w_transform = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
            ),
            name="pv_to_t1w_transformwm",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )
        att_to_t1w_transform = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
            ),
            name="att_to_t1w_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )

        # fmt:off
        workflow.connect([
            (inputnode, basil_to_t1w_transform, [("mean_cbf_basil", "input_image")]),
            (basil_to_t1w_transform, outputnode, [("output_image", "mean_cbf_basil_t1")]),
            (inputnode, basil_to_t1w_transform, [("aslref_to_t1w_xfm", "transforms")]),
            (gen_ref, basil_to_t1w_transform, [("out_file", "reference_image")]),
            (inputnode, pv_to_t1w_transform, [("mean_cbf_gm_basil", "input_image")]),
            (pv_to_t1w_transform, outputnode, [("output_image", "mean_cbf_gm_basil_t1")]),
            (inputnode, pv_to_t1w_transform, [("aslref_to_t1w_xfm", "transforms")]),
            (gen_ref, pv_to_t1w_transform, [("out_file", "reference_image")]),
            (inputnode, pvwm_to_t1w_transform, [("mean_cbf_wm_basil", "input_image")]),
            (pvwm_to_t1w_transform, outputnode, [("output_image", "mean_cbf_wm_basil_t1")]),
            (inputnode, pvwm_to_t1w_transform, [("aslref_to_t1w_xfm", "transforms")]),
            (gen_ref, pvwm_to_t1w_transform, [("out_file", "reference_image")]),
            (inputnode, att_to_t1w_transform, [("att", "input_image")]),
            (att_to_t1w_transform, outputnode, [("output_image", "att_t1")]),
            (inputnode, att_to_t1w_transform, [("aslref_to_t1w_xfm", "transforms")]),
            (gen_ref, att_to_t1w_transform, [("out_file", "reference_image")]),
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
                asl2t1w_init='register',
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
    aslref_to_t1w_xfm
        Affine transform from ``ref_asl_brain`` to T1w space (ITK format)
    itk_t1_to_asl
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
                "aslref_to_t1w_xfm",
                "itk_t1_to_asl",
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
        (fsl2itk_fwd, outputnode, [("itk_transform", "aslref_to_t1w_xfm")]),
        (fsl2itk_inv, outputnode, [("itk_transform", "itk_t1_to_asl")]),
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
