# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for estimating and correcting head motion in ASL images."""

from nipype.interfaces import fsl
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.itk import MCFLIRT2ITK
from niworkflows.utils.connections import listify

from aslprep.config import DEFAULT_MEMORY_MIN_GB
from aslprep.interfaces.confounds import CreateFakeMotionOutputs, NormalizeMotionParams
from aslprep.interfaces.utility import (
    CombineMotionParameters,
    PairwiseRMSDiff,
    SplitByVolumeType,
)


def init_asl_hmc_wf(
    use_ge,
    mem_gb,
    omp_nthreads,
    name='asl_hmc_wf',
):
    """Estimate head-motion parameters and optionally correct them for intensity differences.

    This workflow separately estimates motion parameters for each unique type of volume
    (e.g., control, label, deltam, M0, CBF), and then stitches the resulting parameters
    back together according to the aslcontext file.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.hmc import init_asl_hmc_wf

            wf = init_asl_hmc_wf(
                mem_gb=3,
                omp_nthreads=1,
                name="asl_hmc_wf",
            )

    Parameters
    ----------
    use_ge : :obj:`bool`
        Disable HMC for GE data.
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
    rmsd_file
        Framewise displacement as measured by ``fsl_motion_outliers``

    Notes
    -----
    ASLPrep uses volume type-wise motion correction :footcite:p:`wang2008empirical` instead of the
    zig-zag regression approach :footcite:p:`wang2012improving` because it is unclear how
    M0 volumes should be treated in the zig-zag method.

    References
    ----------
    .. footbibliography::
    """
    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'asl_file',
                'aslcontext',
                'processing_target',
                'raw_ref_image',
            ],
        ),
        name='inputnode',
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'movpar_file',
                'xforms',
                'rmsd_file',
            ],
        ),
        name='outputnode',
    )

    if use_ge:
        # Skip HMC, just mock up movpar and xform files with identity transforms
        workflow.__desc__ = """\
Head-motion correction was not performed on the ASL data.

"""
        create_fake_motion_outputs = pe.Node(
            CreateFakeMotionOutputs(),
            name='create_fake_motion_outputs',
        )
        workflow.connect([
            (inputnode, create_fake_motion_outputs, [('asl_file', 'asl_file')]),
            (create_fake_motion_outputs, [
                ('movpar_file', 'movpar_file'),
                ('xforms', 'xforms'),
                ('rmsd_file', 'rmsd_file'),
            ])
        ])  # fmt:skip
        return workflow

    # Otherwise, perform ASL-specific HMC
    workflow.__desc__ = """\
Head-motion parameters were estimated for the ASL data using *FSL*'s `mcflirt` [@mcflirt].
Motion correction was performed separately for each of the volume types
in order to account for intensity differences between different contrasts,
which, when motion corrected together, can conflate intensity differences with
head motions [@wang2008empirical].
Next, ASLPrep concatenated the motion parameters across volume types and
re-calculated relative root mean-squared deviation.

"""

    split_by_volume_type = pe.Node(
        SplitByVolumeType(),
        name='split_by_volume_type',
    )
    workflow.connect([
        (inputnode, split_by_volume_type, [
            ('aslcontext', 'aslcontext'),
            ('asl_file', 'asl_file'),
        ]),
    ])  # fmt:skip

    mcflirt = pe.MapNode(
        fsl.MCFLIRT(save_mats=True, save_plots=True, save_rms=False),
        name='mcflirt',
        mem_gb=mem_gb * 3,
        iterfield=['in_file'],
    )
    workflow.connect([
        (inputnode, mcflirt, [('raw_ref_image', 'ref_file')]),
        (split_by_volume_type, mcflirt, [('out_files', 'in_file')]),
    ])  # fmt:skip

    listify_mat_files = pe.MapNode(
        niu.Function(
            function=listify,
            input_names=['value'],
            output_names=['lst'],
        ),
        name='listify_mat_files',
        iterfield=['value'],
    )
    workflow.connect([(mcflirt, listify_mat_files, [('mat_file', 'value')])])

    # Combine the motpars files, mat files, and rms files across the different MCFLIRTed files,
    # based on the aslcontext file.
    combine_motpars = pe.Node(
        CombineMotionParameters(),
        name='combine_motpars',
    )
    workflow.connect([
        (inputnode, combine_motpars, [('aslcontext', 'aslcontext')]),
        (split_by_volume_type, combine_motpars, [('volume_types', 'volume_types')]),
        (mcflirt, combine_motpars, [('par_file', 'par_files')]),
        (listify_mat_files, combine_motpars, [('lst', 'mat_files')]),
    ])  # fmt:skip

    # Use rmsdiff to calculate relative rms from transform files.
    rmsdiff = pe.Node(PairwiseRMSDiff(), name='rmsdiff')

    workflow.connect([
        (inputnode, rmsdiff, [('raw_ref_image', 'ref_file')]),
        (combine_motpars, rmsdiff, [('mat_file_list', 'in_files')]),
        (rmsdiff, outputnode, [('out_file', 'rmsd_file')]),
    ])  # fmt:skip

    fsl2itk = pe.Node(MCFLIRT2ITK(), name='fsl2itk', mem_gb=0.05, n_procs=omp_nthreads)

    workflow.connect([
        (inputnode, fsl2itk, [
            ('raw_ref_image', 'in_source'),
            ('raw_ref_image', 'in_reference'),
        ]),
        (combine_motpars, fsl2itk, [('mat_file_list', 'in_files')]),
        (fsl2itk, outputnode, [('out_file', 'xforms')]),
    ])  # fmt:skip

    normalize_motion = pe.Node(
        NormalizeMotionParams(format='FSL'),
        name='normalize_motion',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([
        (combine_motpars, normalize_motion, [('combined_par_file', 'in_file')]),
        (normalize_motion, outputnode, [('out_file', 'movpar_file')]),
    ])  # fmt:skip

    return workflow
