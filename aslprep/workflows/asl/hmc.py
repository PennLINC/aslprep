# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for estimating and correcting head motion in ASL images."""
from nipype.interfaces import fsl
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from aslprep.config import DEFAULT_MEMORY_MIN_GB
from aslprep.interfaces.utility import CombineMotionParameters
from aslprep.niworkflows.engine.workflows import LiterateWorkflow as Workflow
from aslprep.niworkflows.interfaces import NormalizeMotionParams
from aslprep.niworkflows.interfaces.itk import MCFLIRT2ITK
from aslprep.utils.misc import _select_first_in_list


def init_asl_hmc_wf(
    processing_target,
    m0type,
    mem_gb,
    omp_nthreads,
    name="asl_hmc_wf",
):
    """Estimate head-motion parameters and optionally correct them for intensity differences.

    This workflow estimates motion parameters for each unique type of volume
    (e.g., control, label, deltam, M0, CBF).

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.hmc import init_asl_hmc_wf

            wf = init_asl_hmc_wf(
                processing_target="controllabel",
                m0type="Separate",
                mem_gb=3,
                omp_nthreads=1,
                name="asl_hmc_wf",
            )

    Parameters
    ----------
    processing_target : {"controllabel", "deltam", "cbf"}
    m0type
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
                "m0scan_file",
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

    # Combine the motpars files, mat files, and rms files across the different MCFLIRTed files,
    # based on the aslcontext file.
    # XXX: What about a separate M0 scan?
    combine_motpars = pe.Node(
        CombineMotionParameters(m0type=m0type, processing_target=processing_target),
        name="combine_motpars",
    )

    workflow.connect([(inputnode, combine_motpars, [("aslcontext", "aslcontext")])])

    files_to_mcflirt = []
    if m0type == "Included":
        files_to_mcflirt.append("m0scan")

    if processing_target == "controllabel":
        files_to_mcflirt += ["control", "label"]
    else:
        files_to_mcflirt.append(processing_target)

    for file_to_mcflirt in files_to_mcflirt:
        # Head motion correction (hmc)
        mcflirt = pe.Node(
            fsl.MCFLIRT(save_mats=True, save_plots=True, save_rms=True),
            name=f"mcflirt_{file_to_mcflirt}",
            mem_gb=mem_gb * 3,
        )

        # fmt:off
        workflow.connect([
            (inputnode, mcflirt, [
                ("raw_ref_image", "ref_file"),
                (f"{file_to_mcflirt}_file", "in_file"),
            ]),
            (mcflirt, combine_motpars, [
                ("mat_file", f"{file_to_mcflirt}_mat_file"),
                ("par_file", f"{file_to_mcflirt}_par_file"),
                # XXX: Typically grabs relative RMS, but I switched to grab absolute.
                (("rms_files", _select_first_in_list), f"{file_to_mcflirt}_rms_file"),
            ]),
        ])
        # fmt:on

    # TODO: Use rmsdiff to calculate relative rms from transform files.
    fsl2itk = pe.Node(MCFLIRT2ITK(), name="fsl2itk", mem_gb=0.05, n_procs=omp_nthreads)

    # fmt:off
    workflow.connect([
        (inputnode, fsl2itk, [
            ("raw_ref_image", "in_source"),
            ("raw_ref_image", "in_reference"),
        ]),
        (combine_motpars, fsl2itk, [("combined_mat_file", "in_files")]),
        (fsl2itk, outputnode, [("out_file", "xforms")]),
    ])
    # fmt:on

    normalize_motion = pe.Node(
        NormalizeMotionParams(format="FSL"),
        name="normalize_motion",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    workflow.connect([(normalize_motion, outputnode, [("out_file", "movpar_file")])])

    # fmt:off
    workflow.connect([
        (combine_motpars, normalize_motion, [("combined_par_file", "in_file")]),
        (combine_motpars, outputnode, [("combined_rms_file", "rmsd_file")]),
    ])
    # fmt:on

    return workflow
