# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows to process GE ASL data."""
from nipype.interfaces import fsl
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from aslprep import config
from aslprep.interfaces import DerivativesDataSink
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.ge import GeReferenceFile
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.interfaces.masks import SimpleShowMaskRPT
from aslprep.niworkflows.interfaces.utility import KeySelect
from aslprep.niworkflows.interfaces.utils import GenerateSamplingReference
from aslprep.niworkflows.utils.spaces import format_reference
from aslprep.utils.misc import _aslist, _is_native, _select_template, _split_spec
from aslprep.workflows.asl.registration import init_fsl_bbr_wf

DEFAULT_MEMORY_MIN_GB = config.DEFAULT_MEMORY_MIN_GB
LOGGER = config.loggers.workflow


def init_asl_reference_ge_wf(
    metadata,
    aslcontext,
    smooth_kernel=5,
    name="asl_reference_ge_wf",
):
    """Generate a reference volume and its skull-stripped version.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.ge_utils import init_asl_reference_ge_wf

            wf = init_asl_reference_ge_wf(
                metadata={},
                bids_dir=".",
            )
    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
First, a reference volume and its skull-stripped version were generated.
        """

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "m0scan",
                "m0scan_metadata",
            ],
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "raw_ref_image",
                "ref_image_brain",
                "asl_mask",
                "m0_file",
                "m0tr",
                "mask_report",
            ]
        ),
        name="outputnode",
    )

    gen_ref = pe.Node(
        GeReferenceFile(fwhm=smooth_kernel, metadata=metadata, aslcontext=aslcontext),
        omp_nthreads=1,
        mem_gb=1,
        name="gen_ge_ref",
    )

    # fmt:off
    workflow.connect([
        (inputnode, gen_ref, [
            ("asl_file", "in_file"),
            ("m0scan", "m0scan"),
            ("m0scan_metadata", "m0scan_metadata"),
        ]),
        (gen_ref, outputnode, [
            ("ref_file", "raw_ref_image"),
            ("m0_file", "m0_file"),
            ("m0tr", "m0tr"),
        ]),
    ])
    # fmt:on

    skull_strip_wf = pe.Node(fsl.BET(frac=0.5, mask=True), name="fslbet")

    # fmt:off
    workflow.connect([
        (gen_ref, skull_strip_wf, [("ref_file", "in_file")]),
        (skull_strip_wf, outputnode, [("mask_file", "asl_mask")]),
    ])
    # fmt:on

    apply_mask = pe.Node(fsl.ApplyMask(), name="apply_mask")

    # fmt:off
    workflow.connect([
        (gen_ref, apply_mask, [("ref_file", "in_file")]),
        (skull_strip_wf, apply_mask, [("mask_file", "mask_file")]),
        (apply_mask, outputnode, [("out_file", "ref_image_brain")]),
    ])
    # fmt:on

    mask_reportlet = pe.Node(SimpleShowMaskRPT(), name="mask_reportlet")

    # fmt:off
    workflow.connect([
        (gen_ref, mask_reportlet, [("ref_file", "background_file")]),
        (skull_strip_wf, mask_reportlet, [("mask_file", "mask_file")]),
    ])
    # fmt:on

    return workflow


def init_asl_reg_ge_wf(
    use_bbr,
    asl2t1w_dof,
    asl2t1w_init,
    sloppy=False,
    write_report=True,
    name="asl_reg_ge_wf",
):
    """Calculate registration transforms from ASL reference volume to T1w space.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.ge_utils import init_asl_reg_ge_wf

            wf = init_asl_reg_ge_wf(
                use_bbr=True,
                asl2t1w_dof=9,
                asl2t1w_init="register",
            )
    """
    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=["ref_asl_brain", "t1w_brain", "t1w_dseg"]),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["aslref_to_anat_xfm", "anat_to_aslref_xfm", "fallback"]),
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


def init_asl_t1_trans_ge_wf(
    mem_gb,
    omp_nthreads,
    output_t1space=False,
    scorescrub=False,
    basil=False,
    name="asl_t1_trans_ge_wf",
):
    """Co-register the reference ASL image to T1w-space.

    The workflow uses :abbr:`BBR (boundary-based registration)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.ge_utils import init_asl_t1_trans_ge_wf

            wf = init_asl_t1_trans_ge_wf(
                mem_gb=0.1,
                omp_nthreads=1,
            )
    """
    from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow

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
                "aslref_to_anat_xfm",
                # CBF outputs
                "cbf_ts",
                "mean_cbf",
                # SCORE/SCRUB outputs
                "cbf_ts_score",
                "mean_cbf_score",
                "mean_cbf_scrub",
                # BASIL outputs
                "mean_cbf_basil",
                "mean_cbf_gm_basil",
                "mean_cbf_wm_basil",
                "att_basil",
            ],
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_t1",
                "aslref_t1",
                "asl_mask_t1",
                # CBF outputs
                "cbf_ts_t1",
                "mean_cbf_t1",
                # SCORE/SCRUB outputs
                "cbf_ts_score_t1",
                "mean_cbf_score_t1",
                "mean_cbf_scrub_t1",
                # BASIL outputs
                "mean_cbf_basil_t1",
                "mean_cbf_gm_basil_t1",
                "mean_cbf_wm_basil_t1",
                "att_t1",
            ]
        ),
        name="outputnode",
    )

    # gen_ref = pe.Node(GenerateSamplingReference(), name='gen_ref',
    # mem_gb=0.3)  # 256x256x256 * 64 / 8 ~ 150MB

    mask_to_t1w_transform = pe.Node(
        ApplyTransforms(interpolation="MultiLabel"),
        name="mask_to_t1w_transform",
        mem_gb=0.1,
    )

    # fmt:off
    workflow.connect([
        (inputnode, mask_to_t1w_transform, [("ref_asl_mask", "input_image")]),
        (inputnode, mask_to_t1w_transform, [("t1w_brain", "reference_image")]),
        (inputnode, mask_to_t1w_transform, [("aslref_to_anat_xfm", "transforms")]),
        (mask_to_t1w_transform, outputnode, [("output_image", "asl_mask_t1")]),
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
        (inputnode, asl_to_t1w_transform, [("aslref_to_anat_xfm", "transforms")]),
        (inputnode, asl_to_t1w_transform, [("t1w_brain", "reference_image")]),
        (asl_to_t1w_transform, outputnode, [("output_image", "asl_t1")]),
    ])
    # fmt:on

    if not output_t1space:
        return workflow

    input_names = ["cbf_ts", "mean_cbf"]
    if scorescrub:
        input_names += ["cbf_ts_score", "mean_cbf_score", "mean_cbf_scrub"]

    if basil:
        input_names += ["mean_cbf_basil", "mean_cbf_gm_basil", "mean_cbf_wm_basil", "att_basil"]

    for input_name in input_names:
        kwargs = {}
        if input_name in ["cbf_ts", "cbf_ts_score"]:
            kwargs["dimension"] = 3

        warp_input_to_t1w = pe.Node(
            ApplyTransforms(
                interpolation="LanczosWindowedSinc",
                float=True,
                input_image_type=3,
                **kwargs,
            ),
            name=f"warp_{input_name}_to_t1w",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )

        # fmt:off
        workflow.connect([
            (inputnode, warp_input_to_t1w, [
                (input_name, "input_image"),
                ("aslref_to_anat_xfm", "transforms"),
                ("t1w_brain", "reference_image"),
            ]),
            (warp_input_to_t1w, outputnode, [("output_image", f"{input_name}_t1")]),
        ])
        # fmt:on

    return workflow


def init_asl_std_trans_ge_wf(
    mem_gb,
    omp_nthreads,
    spaces,
    scorescrub=False,
    basil=False,
    name="asl_std_trans_ge_wf",
):
    """Resample ASL and CBF derivatives into target spaces.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.ge_utils import init_asl_std_trans_ge_wf

            wf = init_asl_std_trans_ge_wf(
                mem_gb=0.1,
                omp_nthreads=1,
                spaces="",
            )

    Parameters
    ----------
    mem_gb
    omp_nthreads
    spaces
    scorescrub
    basil
    name

    Outputs
    -------
    asl_mask_std
    asl_std
    aslref_std
    spatial_reference
    template
    cbf_ts_std
    mean_cbf_std
    cbf_ts_score_std
    mean_cbf_score_std
    mean_cbf_scrub_std
    mean_cbf_basil_std
    mean_cbf_gm_basil_std
    mean_cbf_wm_basil_std
    att_std
    """
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
                "aslref_to_anat_xfm",
                "name_source",
                "templates",
                "anat_to_template_xfm",
                "asl_file",
                "asl_mask",
                # CBF outputs
                "cbf_ts",
                "mean_cbf",
                # SCORE/SCRUB outputs
                "cbf_ts_score",
                "mean_cbf_score",
                "mean_cbf_scrub",
                # BASIL outputs
                "mean_cbf_basil",
                "mean_cbf_gm_basil",
                "mean_cbf_wm_basil",
                "att_basil",
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
        KeySelect(fields=["anat_to_template_xfm"]),
        name="select_std",
        run_without_submitting=True,
    )

    select_tpl = pe.Node(
        niu.Function(function=_select_template),
        name="select_tpl",
        run_without_submitting=True,
    )

    mask_std_tfm = pe.Node(
        ApplyTransforms(interpolation="MultiLabel"),
        name="mask_std_tfm",
        mem_gb=1,
    )

    # Write corrected file in the designated output dir
    mask_merge_tfms = pe.Node(
        niu.Merge(2),
        name="mask_merge_tfms",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    merge_xforms = pe.Node(
        niu.Merge(3),
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
            ("anat_to_template_xfm", "anat_to_template_xfm"),
            ("templates", "keys"),
        ]),
        (inputnode, mask_std_tfm, [("asl_mask", "input_image")]),
        (inputnode, gen_ref, [("asl_file", "moving_image")]),
        (inputnode, merge_xforms, [(("aslref_to_anat_xfm", _aslist), "in2")]),
        (inputnode, mask_merge_tfms, [(("aslref_to_anat_xfm", _aslist), "in2")]),
        (inputnode, asl_to_std_transform, [("asl_file", "input_image")]),
        (split_target, select_std, [("space", "key")]),
        (select_std, merge_xforms, [("anat_to_template_xfm", "in1")]),
        (select_std, mask_merge_tfms, [("anat_to_template_xfm", "in1")]),
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
        "aslref_std",
        "spatial_reference",
        "template",
        "cbf_ts_std",
        "mean_cbf_std",
    ]

    if scorescrub:
        output_names += [
            "cbf_ts_score_std",
            "mean_cbf_score_std",
            "mean_cbf_scrub_std",
        ]

    if basil:
        output_names += [
            "mean_cbf_basil_std",
            "mean_cbf_gm_basil_std",
            "mean_cbf_wm_basil_std",
            "att_std",
        ]

    poutputnode = pe.Node(niu.IdentityInterface(fields=output_names), name="poutputnode")

    # fmt:off
    workflow.connect([
        # Connecting outputnode
        (iterablesource, poutputnode, [(("std_target", format_reference), "spatial_reference")]),
        (asl_to_std_transform, poutputnode, [("output_image", "asl_std")]),
        (asl_to_std_transform, poutputnode, [("output_image", "aslref_std")]),
        (mask_std_tfm, poutputnode, [("output_image", "asl_mask_std")]),
        (select_std, poutputnode, [("key", "template")]),
        (mask_merge_tfms, cbf_to_std_transform, [("out", "transforms")]),
        (gen_ref, cbf_to_std_transform, [("out_file", "reference_image")]),
        (inputnode, cbf_to_std_transform, [("cbf_ts", "input_image")]),
        (cbf_to_std_transform, poutputnode, [("output_image", "cbf_ts_std")]),
        (mask_merge_tfms, meancbf_to_std_transform, [("out", "transforms")]),
        (gen_ref, meancbf_to_std_transform, [("out_file", "reference_image")]),
        (inputnode, meancbf_to_std_transform, [("cbf_ts", "input_image")]),
        (meancbf_to_std_transform, poutputnode, [("output_image", "mean_cbf_std")]),
    ])
    # fmt:on

    if scorescrub:
        # fmt:off
        workflow.connect([
            (mask_merge_tfms, avgscore_to_std_transform, [("out", "transforms")]),
            (gen_ref, avgscore_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, avgscore_to_std_transform, [("mean_cbf_score", "input_image")]),
            (avgscore_to_std_transform, poutputnode, [("output_image", "mean_cbf_score_std")]),
            (mask_merge_tfms, score_to_std_transform, [("out", "transforms")]),
            (gen_ref, score_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, score_to_std_transform, [("cbf_ts_score", "input_image")]),
            (score_to_std_transform, poutputnode, [("output_image", "cbf_ts_score_std")]),
            (mask_merge_tfms, scrub_to_std_transform, [("out", "transforms")]),
            (gen_ref, scrub_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, scrub_to_std_transform, [("mean_cbf_scrub", "input_image")]),
            (scrub_to_std_transform, poutputnode, [("output_image", "mean_cbf_scrub_std")]),
        ])
        # fmt:on

    if basil:
        # fmt:off
        workflow.connect([
            (mask_merge_tfms, basil_to_std_transform, [("out", "transforms")]),
            (gen_ref, basil_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, basil_to_std_transform, [("mean_cbf_basil", "input_image")]),
            (basil_to_std_transform, poutputnode, [("output_image", "mean_cbf_basil_std")]),
            (mask_merge_tfms, pv_to_std_transform, [("out", "transforms")]),
            (gen_ref, pv_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, pv_to_std_transform, [("mean_cbf_gm_basil", "input_image")]),
            (pv_to_std_transform, poutputnode, [("output_image", "mean_cbf_gm_basil_std")]),
            (mask_merge_tfms, pvwm_to_std_transform, [("out", "transforms")]),
            (gen_ref, pvwm_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, pvwm_to_std_transform, [("mean_cbf_wm_basil", "input_image")]),
            (pvwm_to_std_transform, poutputnode, [("output_image", "mean_cbf_wm_basil_std")]),
            (mask_merge_tfms, att_to_std_transform, [("out", "transforms")]),
            (gen_ref, att_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, att_to_std_transform, [("att_basil", "input_image")]),
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
