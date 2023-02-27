# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Slice-Timing Correction (STC) of ASL images."""
from nipype.interfaces import afni
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from aslprep import config
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.interfaces.utils import CopyXForm

LOGGER = config.loggers.workflow


def init_asl_stc_wf(metadata, name="asl_stc_wf"):
    """Create a workflow for :abbr:`STC (slice-timing correction)`.

    This workflow performs :abbr:`STC (slice-timing correction)` over the input
    :abbr:`ASL (arterial spin labeling)` image.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl import init_asl_stc_wf
            wf = init_asl_stc_wf(
                metadata={"RepetitionTime": 2.0,
                          "SliceTiming": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]},
                )

    Parameters
    ----------
    metadata : :obj:`dict`
        BIDS metadata for ASL file
    name : :obj:`str`
        Name of workflow (default: ``asl_stc_wf``)

    Inputs
    ------
    asl_file
        ASL series NIfTI file
    skip_vols
        Number of non-steady-state volumes detected at beginning of ``asl_file``

    Outputs
    -------
    stc_file
        Slice-timing corrected ASL series NIfTI file

    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
ASL runs were slice-time corrected using `3dTshift` from AFNI [@afni].
"""
    inputnode = pe.Node(niu.IdentityInterface(fields=["asl_file", "skip_vols"]), name="inputnode")
    outputnode = pe.Node(niu.IdentityInterface(fields=["stc_file"]), name="outputnode")

    LOGGER.log(25, "Slice-timing correction will be included.")

    # It would be good to fingerprint memory use of afni.TShift
    slice_timing_correction = pe.Node(
        afni.TShift(
            outputtype="NIFTI_GZ",
            tr=f"{metadata['RepetitionTime']}s",
            slice_timing=metadata["SliceTiming"],
            slice_encoding_direction=metadata.get("SliceEncodingDirection", "k"),
        ),
        name="slice_timing_correction",
    )

    copy_xform = pe.Node(CopyXForm(), name="copy_xform", mem_gb=0.1)

    # fmt:off
    workflow.connect([
        (inputnode, slice_timing_correction, [
            ("asl_file", "in_file"), ("skip_vols", "ignore"),
        ]),
        (slice_timing_correction, copy_xform, [("out_file", "in_file")]),
        (inputnode, copy_xform, [("asl_file", "hdr_file")]),
        (copy_xform, outputnode, [("out_file", "stc_file")]),
    ])
    # fmt:on

    return workflow
