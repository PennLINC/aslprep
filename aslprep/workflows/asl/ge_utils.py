# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows to process GE ASL data."""
import os

from nipype.interfaces import fsl
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from niworkflows.interfaces.nibabel import GenerateSamplingReference
from niworkflows.interfaces.reportlets.masks import SimpleShowMaskRPT
from niworkflows.interfaces.utility import KeySelect

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.interfaces.ge import GeReferenceFile
from aslprep.utils.misc import _aslist, _is_native, _select_template, _split_spec
from aslprep.utils.niworkflows import format_reference
from aslprep.workflows.asl.registration import init_fsl_bbr_wf

DEFAULT_MEMORY_MIN_GB = config.DEFAULT_MEMORY_MIN_GB
LOGGER = config.loggers.workflow


def init_asl_geref_wf(
    metadata,
    bids_dir,
    smooth_kernel=5,
    name="asl_gereference_wf",
):
    """Generate a reference volume and its skull-stripped version."""
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
First, a reference volume and its skull-stripped version were generated.
        """

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["asl_file"]),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "raw_ref_image",
                "ref_image_brain",
                "asl_mask",
                "m0_file",
                "mask_report",
            ]
        ),
        name="outputnode",
    )

    gen_ref = pe.Node(
        GeReferenceFile(bids_dir=bids_dir, fwhm=smooth_kernel, in_metadata=metadata),
        omp_nthreads=1,
        mem_gb=1,
        name="gen_ge_ref",
    )
    gen_ref.base_dir = os.getcwd()
    skull_strip_wf = pe.Node(fsl.BET(frac=0.5, mask=True), name="fslbet")
    apply_mask = pe.Node(fsl.ApplyMask(), name="apply_mask")
    mask_reportlet = pe.Node(SimpleShowMaskRPT(), name="mask_reportlet")

    # fmt:off
    workflow.connect([
        (inputnode, gen_ref, [("asl_file", "in_file")]),
        (gen_ref, skull_strip_wf, [("ref_file", "in_file")]),
        (gen_ref, outputnode, [("ref_file", "raw_ref_image")]),
        (gen_ref, apply_mask, [("ref_file", "in_file")]),
        (skull_strip_wf, outputnode, [("mask_file", "asl_mask")]),
        (skull_strip_wf, apply_mask, [("mask_file", "mask_file")]),
        (apply_mask, outputnode, [("out_file", "ref_image_brain")]),
        (gen_ref, mask_reportlet, [("ref_file", "background_file")]),
        (skull_strip_wf, mask_reportlet, [("mask_file", "mask_file")]),
        (gen_ref, outputnode, [("m0_file", "m0_file")]),
    ])
    # fmt:on

    return workflow


def init_asl_gereg_wf(
    use_bbr,
    asl2t1w_dof,
    asl2t1w_init,
    name="asl_reg_wf",
    sloppy=False,
    write_report=True,
):
    """Calculate registration transforms from ASL reference volume to T1w space."""
    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=["ref_asl_brain", "t1w_brain", "t1w_dseg"]),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["itk_asl_to_t1", "itk_t1_to_asl", "fallback"]),
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
            ("outputnode.itk_asl_to_t1", "itk_asl_to_t1"),
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


def init_asl_t1_getrans_wf(
    mem_gb,
    omp_nthreads,
    cbft1space=False,
    scorescrub=False,
    basil=False,
    name="asl_t1_trans_wf",
):
    """Co-register the reference ASL image to T1w-space.

    The workflow uses :abbr:`BBR (boundary-based registration)`.
    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms

    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "name_source",
                "ref_asl_brain",
                "ref_asl_mask",
                "asl_file",
                "t1w_brain",
                "t1w_mask",
                "cbf",
                "meancbf",
                "att",
                "score",
                "avgscore",
                "scrub",
                "basil",
                "pv",
                "pvwm",
                "itk_asl_to_t1",
            ]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_t1",
                "asl_t1_ref",
                "asl_mask_t1",
                "att_t1",
                "cbf_t1",
                "meancbf_t1",
                "score_t1",
                "avgscore_t1",
                "scrub_t1",
                "basil_t1",
                "pv_t1",
                "pvwm_t1",
            ]
        ),
        name="outputnode",
    )

    # gen_ref = pe.Node(GenerateSamplingReference(), name='gen_ref',
    # mem_gb=0.3)  # 256x256x256 * 64 / 8 ~ 150MB

    mask_t1w_tfm = pe.Node(
        ApplyTransforms(interpolation="MultiLabel"),
        name="mask_t1w_tfm",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, mask_t1w_tfm, [("ref_asl_mask", "input_image")]),
        (inputnode, mask_t1w_tfm, [("t1w_brain", "reference_image")]),
        (inputnode, mask_t1w_tfm, [("itk_asl_to_t1", "transforms")]),
        (mask_t1w_tfm, outputnode, [("output_image", "asl_mask_t1")]),
    ])
    # fmt:on

    asl_to_t1w_transform = pe.Node(
        ApplyTransforms(
            interpolation="LanczosWindowedSinc",
            float=True,
            input_image_type=3,
            dimension=3,
        ),
        name="asl_to_t1w_transform",
        mem_gb=mem_gb,
    )

    # Generate a reference on the target T1w space
    # fmt:off
    workflow.connect([
        (inputnode, asl_to_t1w_transform, [("ref_asl_brain", "input_image")]),
        (inputnode, asl_to_t1w_transform, [("itk_asl_to_t1", "transforms")]),
        (inputnode, asl_to_t1w_transform, [("t1w_brain", "reference_image")]),
        (asl_to_t1w_transform, outputnode, [("output_image", "asl_t1")]),
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
            (asl_to_t1w_transform, outputnode, [("output_image", "asl_t1_ref")]),
            (inputnode, cbf_to_t1w_transform, [("cbf", "input_image")]),
            (cbf_to_t1w_transform, outputnode, [("output_image", "cbf_t1")]),
            (inputnode, cbf_to_t1w_transform, [("itk_asl_to_t1", "transforms")]),
            (inputnode, cbf_to_t1w_transform, [("t1w_brain", "reference_image")]),
            (inputnode, meancbf_to_t1w_transform, [("meancbf", "input_image")]),
            (meancbf_to_t1w_transform, outputnode, [("output_image", "meancbf_t1")]),
            (inputnode, meancbf_to_t1w_transform, [("itk_asl_to_t1", "transforms")]),
            (inputnode, meancbf_to_t1w_transform, [("t1w_brain", "reference_image")]),
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
            name="pvwm_to_t1w_transform",
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
            (inputnode, score_to_t1w_transform, [
                ("score", "input_image"),
                ("itk_asl_to_t1", "transforms"),
                ("t1w_brain", "reference_image"),
            ]),
            (score_to_t1w_transform, outputnode, [("output_image", "score_t1")]),
            (inputnode, avgscore_to_t1w_transform, [
                ("avgscore", "input_image"),
                ("itk_asl_to_t1", "transforms"),
                ("t1w_brain", "reference_image"),
            ]),
            (avgscore_to_t1w_transform, outputnode, [("output_image", "avgscore_t1")]),
            (inputnode, scrub_to_t1w_transform, [
                ("scrub", "input_image"),
                ("itk_asl_to_t1", "transforms"),
                ("t1w_brain", "reference_image"),
            ]),
            (scrub_to_t1w_transform, outputnode, [("output_image", "scrub_t1")]),
            (inputnode, basil_to_t1w_transform, [
                ("basil", "input_image"),
                ("itk_asl_to_t1", "transforms"),
                ("t1w_brain", "reference_image"),
            ]),
            (basil_to_t1w_transform, outputnode, [("output_image", "basil_t1")]),
            (inputnode, pv_to_t1w_transform, [
                ("pv", "input_image"),
                ("itk_asl_to_t1", "transforms"),
                ("t1w_brain", "reference_image"),
            ]),
            (pv_to_t1w_transform, outputnode, [("output_image", "pv_t1")]),
            (inputnode, pvwm_to_t1w_transform, [
                ("pvwm", "input_image"),
                ("itk_asl_to_t1", "transforms"),
                ("t1w_brain", "reference_image"),
            ]),
            (pvwm_to_t1w_transform, outputnode, [("output_image", "pvwm_t1")]),
            (inputnode, att_to_t1w_transform, [
                ("att", "input_image"),
                ("itk_asl_to_t1", "transforms"),
                ("t1w_brain", "reference_image"),
            ]),
            (att_to_t1w_transform, outputnode, [("output_image", "att_t1")]),
        ])
        # fmt:on

    return workflow


def init_asl_gestd_trans_wf(
    mem_gb,
    omp_nthreads,
    spaces,
    scorescrub=False,
    basil=False,
    name="asl_gestd_trans_wf",
):
    """Resample ASL and CBF derivatives into target spaces."""
    workflow = Workflow(name=name)
    output_references = spaces.cached.get_spaces(nonstandard=False, dim=(3,))
    std_vol_references = [
        (s.fullname, s.spec) for s in spaces.references if s.standard and s.dim == 3
    ]

    if len(output_references) == 1:
        workflow.__desc__ = f"""\
The ASL and CBF dreivatives  were resampled into standard space,
generating a *preprocessed ASL and computed CBF in {output_references[0]} space*.
"""
    elif len(output_references) > 1:
        workflow.__desc__ = f"""\
The ASL and CBF dreivatives were resampled into several standard spaces,
correspondingly generating the following *spatially-normalized,
preprocessed ASL runs*: {', '.join(output_references)}.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "anat2std_xfm",
                "cbf",
                "meancbf",
                "att",
                "asl_file",
                "score",
                "avgscore",
                "scrub",
                "basil",
                "pv",
                "pvwm",
                "asl_mask",
                "itk_asl_to_t1",
                "name_source",
                "templates",
            ]
        ),
        name="inputnode",
    )

    iterablesource = pe.Node(niu.IdentityInterface(fields=["std_target"]), name="iterablesource")
    # Generate conversions for every template+spec at the input
    iterablesource.iterables = [("std_target", std_vol_references)]

    split_target = pe.Node(
        niu.Function(
            function=_split_spec,
            input_names=["in_target"],
            output_names=["space", "template", "spec"],
        ),
        run_without_submitting=True,
        name="split_target",
    )

    select_std = pe.Node(
        KeySelect(fields=["anat2std_xfm"]), name="select_std", run_without_submitting=True
    )

    select_tpl = pe.Node(
        niu.Function(function=_select_template), name="select_tpl", run_without_submitting=True
    )

    mask_std_tfm = pe.Node(
        ApplyTransforms(interpolation="MultiLabel"), name="mask_std_tfm", mem_gb=1
    )

    # Write corrected file in the designated output dir
    mask_merge_tfms = pe.Node(
        niu.Merge(2),
        name="mask_merge_tfms",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    nxforms = 3
    merge_xforms = pe.Node(
        niu.Merge(nxforms),
        name="merge_xforms",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    asl_to_std_transform = pe.Node(
        ApplyTransforms(
            interpolation="LanczosWindowedSinc", float=True, input_image_type=3, dimension=3
        ),
        name="asl_to_std_transform",
        mem_gb=mem_gb * 3 * omp_nthreads,
        n_procs=omp_nthreads,
    )
    cbf_to_std_transform = pe.Node(
        ApplyTransforms(
            interpolation="LanczosWindowedSinc", float=True, input_image_type=3, dimension=3
        ),
        name="cbf_to_std_transform",
        mem_gb=mem_gb * 3 * omp_nthreads,
        n_procs=omp_nthreads,
    )
    meancbf_to_std_transform = pe.Node(
        ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
        name="meancbf_to_std_transform",
        mem_gb=mem_gb * 3 * omp_nthreads,
        n_procs=omp_nthreads,
    )

    if scorescrub:
        score_to_std_transform = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc", float=True, input_image_type=3, dimension=3
            ),
            name="score_to_std_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )

        avgscore_to_std_transform = pe.Node(
            ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
            name="avgscore_to_std_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )

        scrub_to_std_transform = pe.Node(
            ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
            name="scrub_to_std_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )

    if basil:
        basil_to_std_transform = pe.Node(
            ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
            name="basil_to_std_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )

        pv_to_std_transform = pe.Node(
            ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
            name="pv_to_std_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )
        pvwm_to_std_transform = pe.Node(
            ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
            name="pvwm_to_std_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )
        att_to_std_transform = pe.Node(
            ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
            name="att_to_std_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )
    # merge = pe.Node(Merge(compress=use_compression), name='merge',
    # mem_gb=mem_gb * 3)
    mask_merge_tfms = pe.Node(
        niu.Merge(2),
        name="mask_merge_tfms",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    # Generate a reference on the target standard space
    gen_ref = pe.Node(GenerateSamplingReference(), name="gen_ref", mem_gb=0.3)

    # fmt:off
    workflow.connect([
        (iterablesource, split_target, [("std_target", "in_target")]),
        (iterablesource, select_tpl, [("std_target", "template")]),
        (inputnode, select_std, [
            ("anat2std_xfm", "anat2std_xfm"),
            ("templates", "keys"),
        ]),
        (inputnode, mask_std_tfm, [("asl_mask", "input_image")]),
        (inputnode, gen_ref, [("asl_file", "moving_image")]),
        (inputnode, merge_xforms, [(("itk_asl_to_t1", _aslist), "in2")]),
        (inputnode, mask_merge_tfms, [(("itk_asl_to_t1", _aslist), "in2")]),
        (inputnode, asl_to_std_transform, [("asl_file", "input_image")]),
        (split_target, select_std, [("space", "key")]),
        (select_std, merge_xforms, [("anat2std_xfm", "in1")]),
        (select_std, mask_merge_tfms, [("anat2std_xfm", "in1")]),
        (split_target, gen_ref, [(("spec", _is_native), "keep_native")]),
        (select_tpl, gen_ref, [("out", "fixed_image")]),
        (merge_xforms, asl_to_std_transform, [("out", "transforms")]),
        (gen_ref, asl_to_std_transform, [("out_file", "reference_image")]),
        (gen_ref, mask_std_tfm, [("out_file", "reference_image")]),
        (mask_merge_tfms, mask_std_tfm, [("out", "transforms")]),
    ])
    # fmt:on

    output_names = [
        "asl_mask_std",
        "asl_std",
        "asl_std_ref",
        "spatial_reference",
        "template",
        "cbf_std",
        "meancbf_std",
    ]

    if scorescrub:
        output_names = output_names + ["score_std", "avgscore_std", "scrub_std"]

    if basil:
        output_names = output_names + ["basil_std", "pv_std", "pvwm_std", "att_std"]

    poutputnode = pe.Node(niu.IdentityInterface(fields=output_names), name="poutputnode")

    # fmt:off
    workflow.connect([
        # Connecting outputnode
        (iterablesource, poutputnode, [(("std_target", format_reference), "spatial_reference")]),
        (asl_to_std_transform, poutputnode, [("output_image", "asl_std")]),
        (asl_to_std_transform, poutputnode, [("output_image", "asl_std_ref")]),
        (mask_std_tfm, poutputnode, [("output_image", "asl_mask_std")]),
        (select_std, poutputnode, [("key", "template")]),
        (mask_merge_tfms, cbf_to_std_transform, [("out", "transforms")]),
        (gen_ref, cbf_to_std_transform, [("out_file", "reference_image")]),
        (inputnode, cbf_to_std_transform, [("cbf", "input_image")]),
        (cbf_to_std_transform, poutputnode, [("output_image", "cbf_std")]),
        (mask_merge_tfms, meancbf_to_std_transform, [("out", "transforms")]),
        (gen_ref, meancbf_to_std_transform, [("out_file", "reference_image")]),
        (inputnode, meancbf_to_std_transform, [("cbf", "input_image")]),
        (meancbf_to_std_transform, poutputnode, [("output_image", "meancbf_std")]),
    ])
    # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (mask_merge_tfms, avgscore_to_std_transform, [("out", "transforms")]),
            (gen_ref, avgscore_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, avgscore_to_std_transform, [("avgscore", "input_image")]),
            (avgscore_to_std_transform, poutputnode, [("output_image", "avgscore_std")]),
            (mask_merge_tfms, score_to_std_transform, [("out", "transforms")]),
            (gen_ref, score_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, score_to_std_transform, [("score", "input_image")]),
            (score_to_std_transform, poutputnode, [("output_image", "score_std")]),
            (mask_merge_tfms, scrub_to_std_transform, [("out", "transforms")]),
            (gen_ref, scrub_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, scrub_to_std_transform, [("scrub", "input_image")]),
            (scrub_to_std_transform, poutputnode, [("output_image", "scrub_std")]),
        ])
        # fmt:on

    if basil:
        # fmt:off
        workflow.connect([
            (mask_merge_tfms, basil_to_std_transform, [("out", "transforms")]),
            (gen_ref, basil_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, basil_to_std_transform, [("basil", "input_image")]),
            (basil_to_std_transform, poutputnode, [("output_image", "basil_std")]),
            (mask_merge_tfms, pv_to_std_transform, [("out", "transforms")]),
            (gen_ref, pv_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, pv_to_std_transform, [("pv", "input_image")]),
            (pv_to_std_transform, poutputnode, [("output_image", "pv_std")]),
            (mask_merge_tfms, pvwm_to_std_transform, [("out", "transforms")]),
            (gen_ref, pvwm_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, pvwm_to_std_transform, [("pvwm", "input_image")]),
            (pvwm_to_std_transform, poutputnode, [("output_image", "pvwm_std")]),
            (mask_merge_tfms, att_to_std_transform, [("out", "transforms")]),
            (gen_ref, att_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, att_to_std_transform, [("att", "input_image")]),
            (att_to_std_transform, poutputnode, [("output_image", "att_std")]),
        ])
        # fmt:on

    # Connect parametric outputs to a Join outputnode
    outputnode = pe.JoinNode(
        niu.IdentityInterface(fields=output_names),
        name="outputnode",
        joinsource="iterablesource",
    )
    # fmt:off
    workflow.connect([(poutputnode, outputnode, [(f, f) for f in output_names])])
    # fmt:on

    return workflow
