# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for resampling data."""
from nipype.interfaces import freesurfer as fs
from nipype.interfaces import utility as niu
from nipype.interfaces.fsl import Split as FSLSplit
from nipype.interfaces.io import FreeSurferSource
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.func.util import init_bold_reference_wf
from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from niworkflows.interfaces.itk import MultiApplyTransforms
from niworkflows.interfaces.nibabel import GenerateSamplingReference
from niworkflows.interfaces.nilearn import Merge
from niworkflows.interfaces.surf import GiftiSetAnatomicalStructure
from niworkflows.interfaces.utility import KeySelect

from aslprep.config import DEFAULT_MEMORY_MIN_GB
from aslprep.utils.misc import (
    _aslist,
    _is_native,
    _itk2lta,
    _select_first_in_list,
    _select_template,
    _split_spec,
    select_target,
)
from aslprep.utils.niworkflows import format_reference


def init_asl_surf_wf(mem_gb, surface_spaces, medial_surface_nan, name="asl_surf_wf"):
    """Sample functional images to FreeSurfer surfaces.

    For each vertex, the cortical ribbon is sampled at six points (spaced 20% of thickness apart)
    and averaged.
    Outputs are in GIFTI format.

    Workflow Graph
        .. workflow::
            :graph2use: colored
            :simple_form: yes

            from aslprep.workflows.asl import init_asl_surf_wf
            wf = init_asl_surf_wf(mem_gb=0.1,
                                   surface_spaces=['fsnative', 'fsaverage5'],
                                   medial_surface_nan=False)

    Parameters
    ----------
    surface_spaces : :obj:`list`
        List of FreeSurfer surface-spaces (either ``fsaverage{3,4,5,6,}`` or ``fsnative``)
        the functional images are to be resampled to.
        For ``fsnative``, images will be resampled to the individual subject's
        native surface.
    medial_surface_nan : :obj:`bool`
        Replace medial wall values with NaNs on functional GIFTI files

    Inputs
    ------
    source_file
        Motion-corrected ASL series in T1 space
    t1w_preproc
        Bias-corrected structural template image
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    t1w2fsnative_xfm
        LTA-style affine matrix translating from T1w to FreeSurfer-conformed subject space

    Outputs
    -------
    surfaces
        ASL series, resampled to FreeSurfer surfaces

    """
    workflow = Workflow(name=name)
    workflow.__desc__ = f"""\
The ASL time-series were resampled onto the following surfaces
(FreeSurfer reconstruction nomenclature):
{', '.join(['*' + s + '*' for s in surface_spaces])}.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=["source_file", "subject_id", "subjects_dir", "t1w2fsnative_xfm"]
        ),
        name="inputnode",
    )
    itersource = pe.Node(niu.IdentityInterface(fields=["target"]), name="itersource")
    itersource.iterables = [("target", surface_spaces)]

    get_fsnative = pe.Node(FreeSurferSource(), name="get_fsnative", run_without_submitting=True)

    targets = pe.Node(
        niu.Function(function=select_target),
        name="targets",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    # Rename the source file to the output space to simplify naming later
    rename_src = pe.Node(
        niu.Rename(format_string="%(subject)s", keep_ext=True),
        name="rename_src",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    itk2lta = pe.Node(niu.Function(function=_itk2lta), name="itk2lta", run_without_submitting=True)
    sampler = pe.MapNode(
        fs.SampleToSurface(
            cortex_mask=True,
            interp_method="trilinear",
            out_type="gii",
            override_reg_subj=True,
            sampling_method="average",
            sampling_range=(0, 1, 0.2),
            sampling_units="frac",
        ),
        iterfield=["hemi"],
        name="sampler",
        mem_gb=mem_gb * 3,
    )
    sampler.inputs.hemi = ["lh", "rh"]
    update_metadata = pe.MapNode(
        GiftiSetAnatomicalStructure(),
        iterfield=["in_file"],
        name="update_metadata",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    outputnode = pe.JoinNode(
        niu.IdentityInterface(fields=["surfaces", "target"]),
        joinsource="itersource",
        name="outputnode",
    )

    workflow.connect(
        [
            (
                inputnode,
                get_fsnative,
                [("subject_id", "subject_id"), ("subjects_dir", "subjects_dir")],
            ),
            (inputnode, targets, [("subject_id", "subject_id")]),
            (inputnode, rename_src, [("source_file", "in_file")]),
            (inputnode, itk2lta, [("source_file", "src_file"), ("t1w2fsnative_xfm", "in_file")]),
            (get_fsnative, itk2lta, [("T1", "dst_file")]),
            (inputnode, sampler, [("subjects_dir", "subjects_dir"), ("subject_id", "subject_id")]),
            (itersource, targets, [("target", "space")]),
            (itersource, rename_src, [("target", "subject")]),
            (itk2lta, sampler, [("out", "reg_file")]),
            (targets, sampler, [("out", "target_subject")]),
            (rename_src, sampler, [("out_file", "source_file")]),
            (update_metadata, outputnode, [("out_file", "surfaces")]),
            (itersource, outputnode, [("target", "target")]),
        ]
    )

    if not medial_surface_nan:
        workflow.connect(sampler, "out_file", update_metadata, "in_file")
        return workflow

    from niworkflows.interfaces.freesurfer import MedialNaNs

    # Refine if medial vertices should be NaNs
    medial_nans = pe.MapNode(
        MedialNaNs(), iterfield=["in_file"], name="medial_nans", mem_gb=DEFAULT_MEMORY_MIN_GB
    )

    workflow.connect(
        [
            (inputnode, medial_nans, [("subjects_dir", "subjects_dir")]),
            (sampler, medial_nans, [("out_file", "in_file")]),
            (medial_nans, update_metadata, [("out_file", "in_file")]),
        ]
    )
    return workflow


def init_asl_std_trans_wf(
    mem_gb,
    omp_nthreads,
    spaces,
    scorescrub=False,
    basil=False,
    name="asl_std_trans_wf",
    use_compression=True,
    use_fieldwarp=False,
):
    """Sample ASL into standard space with a single-step resampling of the original ASL series.

    .. important::
        This workflow provides two outputnodes.
        One output node (with name ``poutputnode``) will be parameterized in a Nipype sense
        (see `Nipype iterables
        <https://miykael.github.io/nipype_tutorial/notebooks/basic_iteration.html>`__), and a
        second node (``outputnode``) will collapse the parameterized outputs into synchronous
        lists of the output fields listed below.

    Workflow Graph
        .. workflow::
            :graph2use: colored
            :simple_form: yes

            from aslprep.utils.niworkflows import SpatialReferences
            from aslprep.workflows.asl import init_asl_std_trans_wf
            wf = init_asl_std_trans_wf(
                mem_gb=3,
                omp_nthreads=1,
                spaces=SpatialReferences(
                    spaces=['MNI152Lin',
                            ('MNIPediatricAsym', {'cohort': '6'})],
                    checkpoint=True),
            )

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of ASL file in GB
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    spaces : :py:class:`~aslprep.utils.niworkflows.SpatialReferences`
        A container for storing, organizing, and parsing spatial normalizations. Composed of
        :py:class:`~aslprep.utils.niworkflows.Reference` objects representing spatial references.
        Each ``Reference`` contains a space, which is a string of either TemplateFlow template IDs.
    name : :obj:`str`
        Name of workflow (default: ``asl_std_trans_wf``)
    use_compression : :obj:`bool`
        Save registered ASL series as ``.nii.gz``
    use_fieldwarp : :obj:`bool`
        Include SDC warp in single-shot transform from ASL to MNI

    Inputs
    ------
    anat2std_xfm
        List of anatomical-to-standard space transforms generated during
        spatial normalization.
    asl_mask
        Skull-stripping mask of reference image
    asl_split
        Individual 3D volumes, not motion corrected
    fieldwarp
        a :abbr:`DFM (displacements field map)` in ITK format
    hmc_xforms
        List of affine transforms aligning each volume to ``ref_image`` in ITK format
    itk_asl_to_t1
        Affine transform from ``ref_asl_brain`` to T1 space (ITK format)
    name_source
        ASL series NIfTI file
        Used to recover original information lost during processing
    templates
        List of templates that were applied as targets during
        spatial normalization.

    Outputs
    -------
    asl_std
        ASL series, resampled to template space
    cbf_std, *cbf
        cbf series, resampled to template space
    asl_std_ref
        Reference, contrast-enhanced summary of the ASL series, resampled to template space
    asl_mask_std
        ASL series mask in template space
    template
        Template identifiers synchronized correspondingly to previously
        described outputs.

    """
    workflow = Workflow(name=name)
    std_vol_references = [
        (s.fullname, s.spec) for s in spaces.references if s.standard and s.dim == 3
    ]

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "anat2std_xfm",
                "cbf",
                "meancbf",
                "att",
                "score",
                "avgscore",
                "scrub",
                "basil",
                "pv",
                "pvwm",
                "asl_mask",
                "asl_split",
                "fieldwarp",
                "hmc_xforms",
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

    gen_ref = pe.Node(
        GenerateSamplingReference(), name="gen_ref", mem_gb=0.3
    )  # 256x256x256 * 64 / 8 ~ 150MB)

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

    nxforms = 3 + use_fieldwarp
    merge_xforms = pe.Node(
        niu.Merge(nxforms),
        name="merge_xforms",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([(inputnode, merge_xforms, [("hmc_xforms", f"in{nxforms}")])])

    if use_fieldwarp:
        workflow.connect([(inputnode, merge_xforms, [("fieldwarp", "in3")])])

    asl_to_std_transform = pe.Node(
        MultiApplyTransforms(interpolation="LanczosWindowedSinc", float=True, copy_dtype=True),
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
            name="pv_to_std_transformwm",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )
        att_to_std_transform = pe.Node(
            ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
            name="att_to_std_transform",
            mem_gb=mem_gb * 3 * omp_nthreads,
            n_procs=omp_nthreads,
        )

    merge = pe.Node(Merge(compress=use_compression), name="merge", mem_gb=mem_gb * 3)

    # Generate a reference on the target standard space
    gen_final_ref = init_bold_reference_wf(pre_mask=True, brainmask_thresh=0.1)

    workflow.connect(
        [
            (iterablesource, split_target, [("std_target", "in_target")]),
            (iterablesource, select_tpl, [("std_target", "template")]),
            (inputnode, select_std, [("anat2std_xfm", "anat2std_xfm"), ("templates", "keys")]),
            (inputnode, mask_std_tfm, [("asl_mask", "input_image")]),
            (inputnode, gen_ref, [(("asl_split", _select_first_in_list), "moving_image")]),
            (inputnode, merge_xforms, [(("itk_asl_to_t1", _aslist), "in2")]),
            (inputnode, merge, [("name_source", "header_source")]),
            (inputnode, mask_merge_tfms, [(("itk_asl_to_t1", _aslist), "in2")]),
            (inputnode, asl_to_std_transform, [("asl_split", "input_image")]),
            (split_target, select_std, [("space", "key")]),
            (select_std, merge_xforms, [("anat2std_xfm", "in1")]),
            (select_std, mask_merge_tfms, [("anat2std_xfm", "in1")]),
            (split_target, gen_ref, [(("spec", _is_native), "keep_native")]),
            (select_tpl, gen_ref, [("out", "fixed_image")]),
            (merge_xforms, asl_to_std_transform, [("out", "transforms")]),
            (gen_ref, asl_to_std_transform, [("out_file", "reference_image")]),
            (gen_ref, mask_std_tfm, [("out_file", "reference_image")]),
            (mask_merge_tfms, mask_std_tfm, [("out", "transforms")]),
            (mask_std_tfm, gen_final_ref, [("output_image", "inputnode.asl_mask")]),
            (asl_to_std_transform, merge, [("out_files", "in_files")]),
            (merge, gen_final_ref, [("out_file", "inputnode.asl_file")]),
        ]
    )

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
        output_names = output_names + ["basil_std", "pv_std", "att_std", "pvwm_std"]

    poutputnode = pe.Node(niu.IdentityInterface(fields=output_names), name="poutputnode")
    workflow.connect(
        [
            # Connecting outputnode
            (
                iterablesource,
                poutputnode,
                [(("std_target", format_reference), "spatial_reference")],
            ),
            (merge, poutputnode, [("out_file", "asl_std")]),
            (gen_final_ref, poutputnode, [("outputnode.ref_image", "asl_std_ref")]),
            (mask_std_tfm, poutputnode, [("output_image", "asl_mask_std")]),
            (select_std, poutputnode, [("key", "template")]),
            (mask_merge_tfms, cbf_to_std_transform, [("out", "transforms")]),
            (gen_ref, cbf_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, cbf_to_std_transform, [("cbf", "input_image")]),
            (cbf_to_std_transform, poutputnode, [("output_image", "cbf_std")]),
            (mask_merge_tfms, meancbf_to_std_transform, [("out", "transforms")]),
            (gen_ref, meancbf_to_std_transform, [("out_file", "reference_image")]),
            (inputnode, meancbf_to_std_transform, [("meancbf", "input_image")]),
            (meancbf_to_std_transform, poutputnode, [("output_image", "meancbf_std")]),
        ]
    )

    if scorescrub:
        workflow.connect(
            [
                (mask_merge_tfms, score_to_std_transform, [("out", "transforms")]),
                (gen_ref, score_to_std_transform, [("out_file", "reference_image")]),
                (inputnode, score_to_std_transform, [("score", "input_image")]),
                (score_to_std_transform, poutputnode, [("output_image", "score_std")]),
                (mask_merge_tfms, avgscore_to_std_transform, [("out", "transforms")]),
                (gen_ref, avgscore_to_std_transform, [("out_file", "reference_image")]),
                (inputnode, avgscore_to_std_transform, [("avgscore", "input_image")]),
                (avgscore_to_std_transform, poutputnode, [("output_image", "avgscore_std")]),
                (mask_merge_tfms, scrub_to_std_transform, [("out", "transforms")]),
                (gen_ref, scrub_to_std_transform, [("out_file", "reference_image")]),
                (inputnode, scrub_to_std_transform, [("scrub", "input_image")]),
                (scrub_to_std_transform, poutputnode, [("output_image", "scrub_std")]),
            ]
        )

    if basil:
        workflow.connect(
            [
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
            ]
        )

    # Connect parametric outputs to a Join outputnode
    outputnode = pe.JoinNode(
        niu.IdentityInterface(fields=output_names), name="outputnode", joinsource="iterablesource"
    )
    workflow.connect(
        [
            (poutputnode, outputnode, [(f, f) for f in output_names]),
        ]
    )
    return workflow


def init_asl_preproc_trans_wf(
    mem_gb,
    omp_nthreads,
    name="asl_preproc_trans_wf",
    use_compression=True,
    use_fieldwarp=False,
    split_file=False,
    interpolation="LanczosWindowedSinc",
):
    """Resample in native (original) space.

    This workflow resamples the input fMRI in its native (original)
    space in a "single shot" from the original asl series.

    Workflow Graph
        .. workflow::
            :graph2use: colored
            :simple_form: yes

            from aslprep.workflows.asl import init_asl_preproc_trans_wf
            wf = init_asl_preproc_trans_wf(mem_gb=3, omp_nthreads=1)

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of asl file in GB
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    name : :obj:`str`
        Name of workflow (default: ``asl_std_trans_wf``)
    use_compression : :obj:`bool`
        Save registered asl series as ``.nii.gz``
    use_fieldwarp : :obj:`bool`
        Include SDC warp in single-shot transform from asl to MNI
    split_file : :obj:`bool`
        Whether the input file should be splitted (it is a 4D file)
        or it is a list of 3D files (default ``False``, do not split)
    interpolation : :obj:`str`
        Interpolation type to be used by ANTs' ``applyTransforms``
        (default ``'LanczosWindowedSinc'``)

    Inputs
    ------
    asl_file
        Individual 3D volumes, not motion corrected
    asl_mask
        Skull-stripping mask of reference image
    name_source
        asl series NIfTI file
        Used to recover original information lost during processing
    hmc_xforms
        List of affine transforms aligning each volume to ``ref_image`` in ITK format
    fieldwarp
        a :abbr:`DFM (displacements field map)` in ITK format

    Outputs
    -------
    asl
        asl series, resampled in native space, including all preprocessing
    asl_mask
        asl series mask calculated with the new time-series
    asl_ref
        asl reference image: an average-like 3D image of the time-series
    asl_ref_brain
        Same as ``asl_ref``, but once the brain mask has been applied

    """
    workflow = Workflow(name=name)
    # workflow.__desc__ = """\
    # The ASL timeseries were resampled onto their original,
    # native space by applying the transforms to correct for head-motion.
    # These resampled ASL timeseries are referred to as preprocessed ASL
    # """

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=["name_source", "asl_file", "asl_mask", "hmc_xforms", "fieldwarp"]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=["asl", "asl_mask", "asl_ref", "asl_ref_brain"]),
        name="outputnode",
    )

    asl_transform = pe.Node(
        MultiApplyTransforms(interpolation=interpolation, float=True, copy_dtype=True),
        name="asl_transform",
        mem_gb=mem_gb * 3 * omp_nthreads,
        n_procs=omp_nthreads,
    )

    merge = pe.Node(Merge(compress=use_compression), name="merge", mem_gb=mem_gb * 3)

    # Generate a new asl reference
    asl_reference_wf = init_bold_reference_wf(brainmask_thresh=0.1)
    asl_reference_wf.__desc__ = None  # Unset description to avoid second appearance

    workflow.connect(
        [
            (inputnode, merge, [("name_source", "header_source")]),
            (asl_transform, merge, [("out_files", "in_files")]),
            (merge, asl_reference_wf, [("out_file", "inputnode.asl_file")]),
            (merge, outputnode, [("out_file", "asl")]),
            (
                asl_reference_wf,
                outputnode,
                [
                    ("outputnode.ref_image", "asl_ref"),
                    ("outputnode.ref_image_brain", "asl_ref_brain"),
                    ("outputnode.asl_mask", "asl_mask"),
                ],
            ),
        ]
    )

    # Input file is not splitted
    if split_file:
        asl_split = pe.Node(FSLSplit(dimension="t"), name="asl_split", mem_gb=mem_gb * 3)
        workflow.connect(
            [
                (inputnode, asl_split, [("asl_file", "in_file")]),
                (
                    asl_split,
                    asl_transform,
                    [
                        ("out_files", "input_image"),
                        (("out_files", _select_first_in_list), "reference_image"),
                    ],
                ),
            ]
        )
    else:
        workflow.connect(
            [
                (
                    inputnode,
                    asl_transform,
                    [
                        ("asl_file", "input_image"),
                        (("asl_file", _select_first_in_list), "reference_image"),
                    ],
                ),
            ]
        )

    if use_fieldwarp:
        merge_xforms = pe.Node(
            niu.Merge(2),
            name="merge_xforms",
            run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect(
            [
                (inputnode, merge_xforms, [("fieldwarp", "in1"), ("hmc_xforms", "in2")]),
                (merge_xforms, asl_transform, [("out", "transforms")]),
            ]
        )
    else:

        def _aslist(val):
            return [val]

        workflow.connect(
            [
                (inputnode, asl_transform, [(("hmc_xforms", _aslist), "transforms")]),
            ]
        )
    return workflow
