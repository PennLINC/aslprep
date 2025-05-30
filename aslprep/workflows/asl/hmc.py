# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for estimating and correcting head motion in ASL images."""


def init_asl_hmc_wf(
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
    processing_target : {"control", "deltam", "cbf"}
    m0type : {"Separate", "Included", "Absent", "Estimate"}
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

    Notes
    -----
    ASLPrep uses volume type-wise motion correction :footcite:p:`wang2008empirical` instead of the
    zig-zag regression approach :footcite:p:`wang2012improving` because it is unclear how
    M0 volumes should be treated in the zig-zag method.

    References
    ----------
    .. footbibliography::
    """
    from nipype.interfaces import fsl
    from nipype.interfaces import utility as niu
    from nipype.pipeline import engine as pe
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.itk import MCFLIRT2ITK
    from niworkflows.utils.connections import listify

    from aslprep.config import DEFAULT_MEMORY_MIN_GB
    from aslprep.interfaces.confounds import NormalizeMotionParams
    from aslprep.interfaces.utility import (
        CombineMotionParameters,
        PairwiseRMSDiff,
        SplitByVolumeType,
    )

    workflow = Workflow(name=name)

    workflow.__desc__ = """\
Head-motion parameters were estimated for the ASL data using *FSL*'s `mcflirt` [@mcflirt].
Motion correction was performed separately for each of the volume types
in order to account for intensity differences between different contrasts,
which, when motion corrected together, can conflate intensity differences with
head motions [@wang2008empirical].
Next, ASLPrep concatenated the motion parameters across volume types and
re-calculated relative root mean-squared deviation.

"""

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
        fsl.MCFLIRT(save_mats=True, save_plots=True, save_rms=False, cost='mutualinfo'),
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


def init_m0scan_hmc_wf(
    output_dir,
    mem_gb,
    omp_nthreads,
    name='m0scan_hmc_wf',
):
    """Estimate head-motion parameters for the M0 scan.

    This workflow resamples the M0 scan to the ASL reference resolution/orientation,
    and then applies MCFLIRT to the M0 scan.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.hmc import init_m0scan_hmc_wf

            wf = init_m0scan_hmc_wf(
                mem_gb=3,
                omp_nthreads=1,
                name="asl_hmc_wf",
            )

    Parameters
    ----------
    output_dir : :obj:`str`
        Output directory for the M0 scan motion correction transforms
    mem_gb : :obj:`float`
        Memory usage for the workflow in GB.
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    name : :obj:`str`
        Name of workflow (default: ``m0scan_hmc_wf``)

    Inputs
    ------
    m0scan
        M0 scan NIfTI file.
    aslref
        ASL reference image NIfTI file.

    Outputs
    -------
    m0scan2aslref_xfm
        ITKTransform file aligning the M0 scan to the ASL reference.
    m0scan_aslref
        Mean M0 scan image in ASL reference space. Used for fit report workflow.
    """
    from nipype.interfaces import fsl
    from nipype.interfaces import utility as niu
    from nipype.pipeline import engine as pe
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.itk import MCFLIRT2ITK

    from aslprep.interfaces.ants import ApplyTransforms
    from aslprep.interfaces.bids import DerivativesDataSink
    from aslprep.interfaces.utility import Ensure4D, MeanImage

    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'm0scan',
                'aslref',
            ],
        ),
        name='inputnode',
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'm0scan2aslref_xfm',
                'm0scan_aslref',
            ],
        ),
        name='outputnode',
    )

    # Ensure M0 scan is 4D
    ensure_4d = pe.Node(
        Ensure4D(),
        name='ensure_4d',
    )
    workflow.connect([(inputnode, ensure_4d, [('m0scan', 'in_file')])])

    # Resample M0 scan to ASL reference resolution/orientation
    resample_m0scan_to_asl = pe.Node(
        ApplyTransforms(
            input_image_type=3,
            interpolation='Gaussian',
            transforms=['identity'],
            args='--verbose',
        ),
        name='resample_m0scan_to_asl',
    )
    workflow.connect([
        (inputnode, resample_m0scan_to_asl, [('aslref', 'reference_image')]),
        (ensure_4d, resample_m0scan_to_asl, [('out_file', 'input_image')]),
    ])  # fmt:skip

    # Register the M0 scan to the ASL reference.
    # By using MCFLIRT, we can support 4D M0 scans.
    # Register the M0 scan to the ASL reference.
    mcflirt = pe.Node(
        fsl.MCFLIRT(save_mats=True, cost='mutualinfo'),
        name='mcflirt',
        mem_gb=mem_gb,
    )
    workflow.connect([
        (inputnode, mcflirt, [('aslref', 'ref_file')]),
        (resample_m0scan_to_asl, mcflirt, [('output_image', 'in_file')]),
    ])  # fmt:skip

    # Calculate mean image of M0 scan
    mean_m0scan = pe.Node(
        MeanImage(),
        name='mean_m0scan',
        mem_gb=mem_gb,
    )
    workflow.connect([
        (mcflirt, mean_m0scan, [('out_file', 'in_file')]),
        (mean_m0scan, outputnode, [('out_file', 'm0scan_aslref')]),
    ])  # fmt:skip

    fsl2itk = pe.Node(MCFLIRT2ITK(), name='fsl2itk', mem_gb=0.05, n_procs=omp_nthreads)
    workflow.connect([
        (inputnode, fsl2itk, [
            ('aslref', 'in_source'),
            ('aslref', 'in_reference'),
        ]),
        (mcflirt, fsl2itk, [('mat_file', 'in_files')]),
        (fsl2itk, outputnode, [('out_file', 'm0scan2aslref_xfm')]),
    ])  # fmt:skip

    ds_m0scan2aslref_xfm = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            datatype='perf',
            suffix='xfm',
            extension='.txt',
            **{
                'from': 'm0scan',
                'to': 'aslref',
            },
        ),
        name='ds_m0scan2aslref_xfm',
    )
    workflow.connect([
        (inputnode, ds_m0scan2aslref_xfm, [('m0scan', 'source_file')]),
        (outputnode, ds_m0scan2aslref_xfm, [('m0scan2aslref_xfm', 'in_file')]),
    ])  # fmt:skip

    return workflow
