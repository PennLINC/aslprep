# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# Copyright 2021 The NiPreps Developers <nipreps@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# We support and encourage derived works from this project, please read
# about our expectations at
#
#     https://www.nipreps.org/community/licensing/
#
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.header import ValidateImage
from niworkflows.utils.misc import pass_dummy_scans

from aslprep import config
from aslprep.interfaces.reference import SelectHighestContrastVolumes


def init_raw_aslref_wf(
    *,
    asl_file=None,
    m0scan=False,
    name="raw_aslref_wf",
):
    """Build a workflow that generates reference BOLD images for a series.

    CBF images have the best GM-WM contrast, but raw ASL files rarely have precomputed CBF volumes.
    As such, ASLPrep will select volumes in the following order: CBF, M0, deltam, control.

    The raw reference image is the target of :abbr:`HMC (head motion correction)`, and a
    contrast-enhanced reference is the subject of distortion correction, as well as
    boundary-based registration to T1w and template spaces.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.reference import init_raw_aslref_wf

            wf = init_raw_aslref_wf()

    Parameters
    ----------
    asl_file : :obj:`str`
        ASL seris NIfTI file
    m0scan : :obj:`bool`
        True if a separate M0 file is available. False if not.
    name : :obj:`str`
        Name of workflow (default: ``asl_reference_wf``)

    Inputs
    ------
    asl_file : str
        BOLD series NIfTI file
    dummy_scans : int or None
        Number of non-steady-state volumes specified by user at beginning of ``asl_file``

    Outputs
    -------
    asl_file : str
        Validated BOLD series NIfTI file
    aslref : str
        Reference image to which BOLD series is motion corrected
    skip_vols : int
        Number of non-steady-state volumes selected at beginning of ``asl_file``
    algo_dummy_scans : int
        Number of non-steady-state volumes algorithmically detected at
        beginning of ``asl_file``
    """
    from niworkflows.interfaces.images import RobustAverage

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
First, a reference volume was generated using a custom methodology of *ASLPrep*,
for use in head motion correction.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["asl_file", "dummy_scans"]),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "aslref",
                "validation_report",
            ],
        ),
        name="outputnode",
    )

    # Simplify manually setting input image
    if asl_file is not None:
        inputnode.inputs.asl_file = asl_file

    val_asl = pe.Node(
        ValidateImage(),
        name="val_asl",
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    if m0scan:
        val_m0scan = pe.Node(
            ValidateImage(),
            name="val_asl",
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(inputnode, val_m0scan, [("m0scan", "in_file")])])

    select_highest_contrast_volumes = pe.Node(
        SelectHighestContrastVolumes(),
        name="select_highest_contrast_volumes",
        mem_gb=1,
    )
    gen_avg = pe.Node(RobustAverage(), name="gen_avg", mem_gb=1)

    # fmt:off
    workflow.connect([
        (inputnode, val_asl, [("asl_file", "in_file")]),
        (val_asl, select_highest_contrast_volumes, [("out_file", "asl_file")]),
        (inputnode, select_highest_contrast_volumes, [("aslcontext", "aslcontext")]),
        (select_highest_contrast_volumes, gen_avg, [("selected_volumes_file", "in_file")]),
        (val_asl, outputnode, [
            ("out_file", "asl_file"),
            ("out_report", "validation_report"),
        ]),
        (gen_avg, outputnode, [("out_file", "aslref")]),
    ])
    # fmt:on

    if m0scan:
        workflow.connect([(val_m0scan, select_highest_contrast_volumes, [("m0scan", "m0scan")])])

    return workflow


def init_asl_reference_wf(
    asl_file=None,
    aslcontext=None,
    sbref_files=None,
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

            wf = init_asl_reference_wf()

    Parameters
    ----------
    asl_file : :obj:`str`
        ASL series NIfTI file
    sbref_files : :obj:`list` or :obj:`bool`
        Single band (as opposed to multi band) reference NIfTI file.
        If ``True`` is passed, the workflow is built to accommodate SBRefs,
        but the input is left undefined (i.e., it is left open for connection)
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

    reference_target = ""
    if aslcontext is not None:
        reference_target = select_processing_target(aslcontext)
        aslcontext_df = pd.read_table(aslcontext)
        if "m0scan" in aslcontext_df["volume_type"].values:
            reference_target = "M0"

        reference_target += " "

    workflow.__desc__ = f"""\
First, the middle {reference_target}volume of the ASL timeseries was selected as the
reference volume and brain extracted using *Nipype*'s custom brain extraction workflow.
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

    enhance_and_skullstrip_asl_wf = init_enhance_and_skullstrip_asl_wf(pre_mask=pre_mask)

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