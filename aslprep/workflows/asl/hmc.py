# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for estimating and correcting head motion in ASL images."""
from nipype.interfaces import fsl
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from aslprep.config import DEFAULT_MEMORY_MIN_GB
from aslprep.interfaces.confounds import ZigZagCorrection
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.interfaces import NormalizeMotionParams
from aslprep.niworkflows.interfaces.itk import MCFLIRT2ITK
from aslprep.utils.misc import _select_last_in_list


def init_asl_hmc_wf(
    processing_target,
    m0type,
    mem_gb,
    omp_nthreads,
    name="control_label_hmc_wf",
):
    """Estimate head-motion parameters and optionally correct them for intensity differences.

    This workflow first estimates motion parameters using MCFLIRT,
    then applies the zig-zag regressor method :footcite:p:`wang2012improving`
    to correct the parameters for intensity differences between control and label volumes.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.hmc import init_asl_hmc_wf

            wf = init_asl_hmc_wf(
                mem_gb=3,
                omp_nthreads=1,
            )

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of ASL file in GB
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    name : :obj:`str`
        Name of workflow (default: ``asl_hmc_wf``)

    Inputs
    ------
    asl_file
        Control-label pair series NIfTI file.
        If an ASL run contains M0 volumes, deltaM volumes, or CBF volumes,
        those volumes should be removed before running this workflow.
    aslcontext
        ASL context TSV file.
    raw_ref_image
        Reference image to which ASL series is motion corrected

    Outputs
    -------
    xforms
        ITKTransform file aligning each volume to ``ref_image``
    movpar_file
        MCFLIRT motion parameters, normalized to SPM format (X, Y, Z, Rx, Ry, Rz)
    rms_file
        Framewise displacement as measured by ``fsl_motion_outliers``
    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
Head-motion parameters were estimated to the control-label time series using
*FSL*'s `mcflirt` [@mcflirt].
Next, ASLPrep applied zig-zag regression [@wang2012improving] to correct the motion parameters
for intensity differences between the control and label volumes.
ASLPrep wrote the modified head-motion parameters to the ASL run's confound file.

"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "control_file",
                "label_file",
                "deltam_file",
                "cbf_file",
                "m0_file",
                "aslcontext",
                "raw_ref_image",
            ],
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "movpar_file",
                "xforms",
                "rmsd_file",
            ],
        ),
        name="outputnode",
    )

    files_to_mcflirt = []
    if m0type in ("Included", "Separate"):
        files_to_mcflirt.append("m0_file")

    if processing_target == "controllabel":
        files_to_mcflirt += ["control_file", "label_file"]
    else:
        files_to_mcflirt.append(f"{processing_target}_file")

    # Head motion correction (hmc)
    mcflirt = pe.Node(
        fsl.MCFLIRT(save_mats=True, save_plots=True, save_rms=True),
        name="mcflirt",
        mem_gb=mem_gb * 3,
    )

    # fmt:off
    workflow.connect([
        (inputnode, mcflirt, [
            ("raw_ref_image", "ref_file"),
            ("asl_file", "in_file"),
        ]),
    ])
    # fmt:on

    fsl2itk = pe.Node(MCFLIRT2ITK(), name="fsl2itk", mem_gb=0.05, n_procs=omp_nthreads)

    # fmt:off
    workflow.connect([
        (inputnode, fsl2itk, [
            ("raw_ref_image", "in_source"),
            ("raw_ref_image", "in_reference"),
        ]),
        (fsl2itk, outputnode, [("out_file", "xforms")]),
    ])
    # fmt:on

    normalize_motion = pe.Node(
        NormalizeMotionParams(format="FSL"),
        name="normalize_motion",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    workflow.connect([(normalize_motion, outputnode, [("out_file", "movpar_file")])])

    if use_zigzag:
        correct_motion = pe.Node(
            ZigZagCorrection(),
            name="correct_motion",
            mem_gb=DEFAULT_MEMORY_MIN_GB,
        )

        # fmt:off
        workflow.connect([
            (inputnode, correct_motion, [("aslcontext", "aslcontext")]),
            (mcflirt, correct_motion, [("par_file", "movpar_file")]),
            # XXX: This probably won't work since it's a TSV rather than a MAT file
            # Is there a trick for going from TSV files to ITK files?
            (correct_motion, fsl2itk, [("movpar_file", "in_files")]),
            (correct_motion, normalize_motion, [("movpar_file", "in_file")]),
            (correct_motion, outputnode, [("rmsd_file", "rmsd_file")]),
        ])
        # fmt:on
    else:
        # fmt:off
        workflow.connect([
            (mcflirt, fsl2itk, [("mat_file", "in_files")]),
            (mcflirt, normalize_motion, [("par_file", "in_file")]),
            (mcflirt, outputnode, [
                (("rms_files", _select_last_in_list), "rmsd_file"),
            ]),
        ])
        # fmt:on

    return workflow
