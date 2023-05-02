# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Utility workflows."""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.images import ValidateImage
from niworkflows.utils.connections import listify

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
