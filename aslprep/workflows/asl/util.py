# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Utility workflows."""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.header import ValidateImage
from niworkflows.interfaces.reportlets.masks import SimpleShowMaskRPT
from niworkflows.utils.connections import listify
from niworkflows.utils.misc import pass_dummy_scans

from aslprep.interfaces.utility import SplitReferenceTarget
from aslprep.niworkflows.func.util import init_enhance_and_skullstrip_asl_wf
from aslprep.interfaces.niworkflows import EstimateReferenceImage

DEFAULT_MEMORY_MIN_GB = 0.01


def init_validate_asl_wf(asl_file=None, name="validate_asl_wf"):
    """Build a workflow to validate an ASL file.

    This was split out of :func:`~aslprep.workflows.asl.util.init_asl_reference_wf`,
    as that workflow only takes in a subset of the ASL run (typically the M0 scans).

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.util import init_validate_asl_wf

            wf = init_validate_asl_wf()

    Parameters
    ----------
    asl_file : :obj:`str`
        ASL series NIfTI file
    name : :obj:`str`
        Name of workflow (default: ``validate_asl_wf``)

    Inputs
    ------
    asl_file : :obj:`str`
        ASL series NIfTI file

    Outputs
    -------
    asl_file : :obj:`str`
        Validated ASL series NIfTI file
    validation_report : :obj:`str`
        HTML reportlet indicating whether ``asl_file`` had a valid affine
    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
First, the middle volume of the ASL timeseries was selected as the refernce volume and
brain extracted using *Nipype*'s custom brain extraction workflow.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["asl_file"]),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "validation_report",
            ],
        ),
        name="outputnode",
    )

    # Simplify manually setting input image
    if asl_file is not None:
        inputnode.inputs.asl_file = asl_file

    val_asl = pe.MapNode(
        ValidateImage(),
        name="val_asl",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        iterfield=["in_file"],
    )

    workflow.connect([(inputnode, val_asl, [(("asl_file", listify), "in_file")])])

    asl_1st = pe.Node(niu.Select(index=[0]), name="asl_1st", run_without_submitting=True)

    # fmt:off
    workflow.connect([
        (val_asl, asl_1st, [(("out_file", listify), "inlist")]),
        (asl_1st, outputnode, [("out", "asl_file")]),
    ])
    # fmt:on

    validate_1st = pe.Node(niu.Select(index=[0]), name="validate_1st", run_without_submitting=True)

    # fmt:off
    workflow.connect([
        (val_asl, validate_1st, [(("out_report", listify), "inlist")]),
        (validate_1st, outputnode, [("out", "validation_report")]),
    ])
    # fmt:on

    return workflow


def init_asl_reference_wf(
    omp_nthreads,
    asl_file=None,
    sbref_files=None,
    brainmask_thresh=0.1,
    pre_mask=False,
    name="asl_reference_wf",
    gen_report=False,
):
    """
    Build a workflow that generates reference ASL images for a series.

    The raw reference image is the target of :abbr:`HMC (head motion correction)`, and a
    contrast-enhanced reference is the subject of distortion correction, as well as
    boundary-based registration to T1w and template spaces.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.util import init_asl_reference_wf

            wf = init_asl_reference_wf(omp_nthreads=1)

    Parameters
    ----------
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    asl_file : :obj:`str`
        ASL series NIfTI file
    sbref_files : :obj:`list` or :obj:`bool`
        Single band (as opposed to multi band) reference NIfTI file.
        If ``True`` is passed, the workflow is built to accommodate SBRefs,
        but the input is left undefined (i.e., it is left open for connection)
    brainmask_thresh: :obj:`float`
        Lower threshold for the probabilistic brainmask to obtain
        the final binary mask (default: 0.5).
    pre_mask : :obj:`bool`
        Indicates whether the ``pre_mask`` input will be set (and thus, step 1
        should be skipped).
    name : :obj:`str`
        Name of workflow (default: ``asl_reference_wf``)
    gen_report : :obj:`bool`
        Whether a mask report node should be appended in the end

    Inputs
    ------
    asl_file : str
        ASL series NIfTI file
    asl_mask : bool
        A tentative brain mask to initialize the workflow (requires ``pre_mask``
        parameter set ``True``).
    dummy_scans : int or None
        Number of non-steady-state volumes specified by user at beginning of ``asl_file``
    sbref_file : str
        single band (as opposed to multi band) reference NIfTI file

    Outputs
    -------
    asl_file : str
        Validated ASL series NIfTI file
    raw_ref_image : str
        Reference image to which ASL series is motion corrected
    skip_vols : int
        Number of non-steady-state volumes selected at beginning of ``asl_file``
    algo_dummy_scans : int
        Number of non-steady-state volumes agorithmically detected at
        beginning of ``asl_file``
    ref_image : str
        Contrast-enhanced reference image
    ref_image_brain : str
        Skull-stripped reference image
    asl_mask : str
        Skull-stripping mask of reference image
    validation_report : str
        HTML reportlet indicating whether ``asl_file`` had a valid affine


    Subworkflows
        * :py:func:`~niworkflows.func.util.init_enhance_and_skullstrip_wf`

    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
First, the middle volume of the ASL timeseries was selected as the refernce volume and
brain extracted using *Nipype*'s custom brain extraction workflow.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "aslcontext",
                "asl_mask",
                "dummy_scans",
                "sbref_file",
            ],
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "raw_ref_image",
                "skip_vols",
                "algo_dummy_scans",
                "ref_image",
                "ref_image_brain",
                "asl_mask",
                "validation_report",
                "mask_report",
            ]
        ),
        name="outputnode",
    )

    # Simplify manually setting input image
    if asl_file is not None:
        inputnode.inputs.asl_file = asl_file

    # Split up the ASL data into image-type-specific files,
    # and reduce aslcontext, metadata, and asl_file if necessary.
    select_reference_volumes = pe.Node(
        SplitReferenceTarget(),
        name="select_reference_volumes",
    )

    # fmt:off
    workflow.connect([
        (inputnode, select_reference_volumes, [
            ("asl_file", "asl_file"),
            ("aslcontext", "aslcontext"),
        ]),
    ])
    # fmt:on

    val_asl = pe.MapNode(
        ValidateImage(),
        name="val_asl",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        iterfield=["in_file"],
    )

    workflow.connect([(select_reference_volumes, val_asl, [(("out_file", listify), "in_file")])])

    gen_ref = pe.Node(EstimateReferenceImage(), name="gen_ref", mem_gb=1)

    # fmt:off
    workflow.connect([
        (val_asl, gen_ref, [("out_file", "in_file")]),
        (gen_ref, outputnode, [
            ("ref_image", "raw_ref_image"),
            ("n_volumes_to_discard", "algo_dummy_scans"),
        ]),
    ])
    # fmt:on

    enhance_and_skullstrip_asl_wf = init_enhance_and_skullstrip_asl_wf(
        brainmask_thresh=brainmask_thresh,
        omp_nthreads=omp_nthreads,
        pre_mask=pre_mask,
    )

    # fmt:off
    workflow.connect([
        (inputnode, enhance_and_skullstrip_asl_wf, [("asl_mask", "inputnode.pre_mask")]),
        (gen_ref, enhance_and_skullstrip_asl_wf, [("ref_image", "inputnode.in_file")]),
        (enhance_and_skullstrip_asl_wf, outputnode, [
            ("outputnode.bias_corrected_file", "ref_image"),
            ("outputnode.mask_file", "asl_mask"),
            ("outputnode.skull_stripped_file", "ref_image_brain"),
        ]),
    ])
    # fmt:on

    calc_dummy_scans = pe.Node(
        niu.Function(function=pass_dummy_scans, output_names=["skip_vols_num"]),
        name="calc_dummy_scans",
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    # fmt:off
    workflow.connect([
        (inputnode, calc_dummy_scans, [("dummy_scans", "dummy_scans")]),
        (gen_ref, calc_dummy_scans, [("n_volumes_to_discard", "algo_dummy_scans")]),
        (calc_dummy_scans, outputnode, [("skip_vols_num", "skip_vols")]),
    ])
    # fmt:on

    asl_1st = pe.Node(niu.Select(index=[0]), name="asl_1st", run_without_submitting=True)

    # fmt:off
    workflow.connect([
        (val_asl, asl_1st, [(("out_file", listify), "inlist")]),
        (asl_1st, outputnode, [("out", "asl_file")]),
    ])
    # fmt:on

    validate_1st = pe.Node(niu.Select(index=[0]), name="validate_1st", run_without_submitting=True)

    # fmt: off
    workflow.connect([
        (val_asl, validate_1st, [(("out_report", listify), "inlist")]),
        (validate_1st, outputnode, [("out", "validation_report")]),
    ])
    # fmt: on

    if sbref_files:
        nsbrefs = 0
        if sbref_files is not True:
            # If not boolean, then it is a list-of or pathlike.
            inputnode.inputs.sbref_file = sbref_files
            nsbrefs = 1 if isinstance(sbref_files, str) else len(sbref_files)

        val_sbref = pe.MapNode(
            ValidateImage(),
            name="val_sbref",
            mem_gb=DEFAULT_MEMORY_MIN_GB,
            iterfield=["in_file"],
        )
        # fmt: off
        workflow.connect([
            (inputnode, val_sbref, [(("sbref_file", listify), "in_file")]),
            (val_sbref, gen_ref, [("out_file", "sbref_file")]),
        ])
        # fmt: on

        # Edit the boilerplate as the SBRef will be the reference
        workflow.__desc__ = f"""\
First, a reference volume and its skull-stripped version were generated
by aligning and averaging
{nsbrefs or ''} single-band references (SBRefs).
"""

    if gen_report:
        mask_reportlet = pe.Node(SimpleShowMaskRPT(), name="mask_reportlet")
        # fmt: off
        workflow.connect([
            (enhance_and_skullstrip_asl_wf, mask_reportlet, [
                ("outputnode.bias_corrected_file", "background_file"),
                ("outputnode.mask_file", "mask_file"),
            ]),
        ])
        # fmt: on

    return workflow
