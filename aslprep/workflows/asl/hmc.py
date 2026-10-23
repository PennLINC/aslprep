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

    ``config.workflow.hmc_method`` selects the method:

    - ``typewise``: MCFLIRT registers each volume type to the ASL reference.
    - ``typewise_offset`` (default): as ``typewise``, then the offset left between the
      corrected label and control volumes is removed (per post-labeling delay of at least 1 s,
      or inversion time of at least 0.8 s for PASL, the corrected label mean is registered to
      the corrected control mean with ``3dvolreg``; shorter delays are not corrected).
    - ``typewise_volreg`` (**experimental**): each volume type is registered by ``3dvolreg``
      (least squares) to its own mean, which FLIRT aligns to the ASL reference once, and the
      label offset correction always follows.

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
    metadata
        Metadata of ``asl_file``, reduced to its volumes (post-labeling delays and labeling
        type); used by the label offset correction (``typewise_offset`` and
        ``typewise_volreg``).

    Outputs
    -------
    xforms
        ITKTransform file aligning each volume to ``ref_image``. With ``--hmc-method
        typewise_offset`` (the default) or ``typewise_volreg``, the transforms of label volumes
        at delays of at least 1 s (0.8 s for PASL) that also have control volumes remove the
        offset between the corrected labels and controls
        (:class:`~aslprep.interfaces.hmc.CorrectLabelOffset`); other labels are unchanged.
    label_offset_report
        JSON record of the estimated label offsets (``typewise_offset`` and
        ``typewise_volreg`` only).
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

    from aslprep import config
    from aslprep.config import DEFAULT_MEMORY_MIN_GB
    from aslprep.interfaces.confounds import NormalizeMotionParams
    from aslprep.interfaces.hmc import CorrectLabelOffset, VolregSeries
    from aslprep.interfaces.utility import (
        CombineMotionParameters,
        MeanImage,
        PairwiseRMSDiff,
        SplitByVolumeType,
    )

    hmc_method = config.workflow.hmc_method
    workflow = Workflow(name=name)

    if hmc_method == 'typewise_volreg':
        workflow.__desc__ = """\
Head-motion parameters were estimated for the ASL data separately for each of the volume types,
in order to account for intensity differences between different contrasts,
which, when motion corrected together, can conflate intensity differences with
head motions [@wang2008empirical].
The volumes of each type were registered rigidly to the type's middle volume with *AFNI*'s
`3dvolreg` [@afni] and averaged; the mean was aligned to the ASL reference with *FSL*'s `flirt`
(six degrees of freedom, mutual information) [@flirt], and the type's volumes were registered
to that aligned mean with `3dvolreg`.
Next, ASLPrep concatenated the motion parameters across volume types and
re-calculated relative root mean-squared deviation.

"""
    else:
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
                'metadata',
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
                'label_offset_report',
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

    if hmc_method == 'typewise_volreg':
        # Experimental. Each volume type is registered, by least squares, to a reference of its
        # own contrast: its volumes to their middle volume, then their mean, aligned once to the
        # ASL reference by FLIRT (mutual information), and the volumes again to that aligned
        # mean, so the transforms still map to the ASL reference (PennLINC/aslprep#767). The
        # types' separate cross-contrast alignments leave an offset between labels and
        # controls, so the label offset correction below always follows.
        type_volreg = pe.MapNode(
            VolregSeries(),
            name='type_volreg',
            mem_gb=mem_gb * 3,
            iterfield=['in_file'],
        )
        type_mean = pe.MapNode(MeanImage(), name='type_mean', mem_gb=mem_gb, iterfield=['in_file'])
        type_mean_to_ref = pe.MapNode(
            fsl.FLIRT(dof=6, cost='mutualinfo', interp='spline'),
            name='type_mean_to_ref',
            mem_gb=DEFAULT_MEMORY_MIN_GB,
            iterfield=['in_file'],
        )
        registration = pe.MapNode(
            VolregSeries(),
            name='volreg',
            mem_gb=mem_gb * 3,
            iterfield=['in_file', 'ref_file'],
        )
        workflow.connect([
            (split_by_volume_type, type_volreg, [('out_files', 'in_file')]),
            (type_volreg, type_mean, [('out_file', 'in_file')]),
            (type_mean, type_mean_to_ref, [('out_file', 'in_file')]),
            (inputnode, type_mean_to_ref, [('raw_ref_image', 'reference')]),
            (split_by_volume_type, registration, [('out_files', 'in_file')]),
            (type_mean_to_ref, registration, [('out_file', 'ref_file')]),
        ])  # fmt:skip
    else:
        registration = pe.MapNode(
            fsl.MCFLIRT(save_mats=True, save_plots=True, save_rms=False, cost='mutualinfo'),
            name='mcflirt',
            mem_gb=mem_gb * 3,
            iterfield=['in_file'],
        )
        workflow.connect([
            (inputnode, registration, [('raw_ref_image', 'ref_file')]),
            (split_by_volume_type, registration, [('out_files', 'in_file')]),
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
    workflow.connect([(registration, listify_mat_files, [('mat_file', 'value')])])

    # Combine the motpars files, mat files, and rms files across the different MCFLIRTed files,
    # based on the aslcontext file.
    combine_motpars = pe.Node(
        CombineMotionParameters(),
        name='combine_motpars',
    )
    workflow.connect([
        (inputnode, combine_motpars, [('aslcontext', 'aslcontext')]),
        (split_by_volume_type, combine_motpars, [('volume_types', 'volume_types')]),
        (registration, combine_motpars, [('par_file', 'par_files')]),
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
    ])  # fmt:skip

    if hmc_method in ('typewise_offset', 'typewise_volreg'):
        workflow.__desc__ += """\
The systematic offset that type-wise motion correction leaves between the corrected label and
control volumes was then removed: for each post-labeling delay of at least 1.0 s (inversion
time of at least 0.8 s for PASL), the mean of the corrected label volumes was registered rigidly
to the mean of the corrected control volumes with *AFNI*'s `3dvolreg` [@afni], and that
correction was composed into the transforms of the delay's label volumes. Labels at shorter
delays, which may still hold arterial signal, were not corrected.

"""
        correct_label_offset = pe.Node(
            CorrectLabelOffset(),
            name='correct_label_offset',
            mem_gb=mem_gb * 3,
        )
        workflow.connect([
            (inputnode, correct_label_offset, [
                ('asl_file', 'asl_file'),
                ('aslcontext', 'aslcontext'),
                ('metadata', 'metadata'),
                ('raw_ref_image', 'ref_file'),
            ]),
            (fsl2itk, correct_label_offset, [('out_file', 'xforms')]),
            (correct_label_offset, outputnode, [
                ('out_file', 'xforms'),
                ('out_report', 'label_offset_report'),
            ]),
        ])  # fmt:skip
    else:
        workflow.connect([(fsl2itk, outputnode, [('out_file', 'xforms')])])

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

    from aslprep.interfaces.ants import RegridImage
    from aslprep.interfaces.bids import DerivativesDataSink
    from aslprep.interfaces.utility import Ensure4D, GetImageType, MeanImage

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

    # In some cases, the separate M0 scan may have a different resolution than the ASL scan.
    # We resample the M0 scan to match the ASL scan at this point.
    get_image_type = pe.Node(GetImageType(), name='get_image_type')
    workflow.connect([(inputnode, get_image_type, [('m0scan', 'image')])])

    # Resample M0 scan to ASL reference resolution/orientation.
    # An M0 scan already on that grid is passed through, as Gaussian interpolation would only
    # smooth it.
    resample_m0scan_to_asl = pe.Node(
        RegridImage(
            interpolation='Gaussian',
            transforms=['identity'],
            args='--verbose',
        ),
        name='resample_m0scan_to_asl',
    )
    workflow.connect([
        (inputnode, resample_m0scan_to_asl, [
            ('aslref', 'reference_image'),
            ('m0scan', 'input_image'),
        ]),
        (get_image_type, resample_m0scan_to_asl, [('image_type', 'input_image_type')]),
    ])  # fmt:skip

    ensure_4d = pe.Node(Ensure4D(), name='ensure_4d')
    workflow.connect([(resample_m0scan_to_asl, ensure_4d, [('output_image', 'in_file')])])

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
        (ensure_4d, mcflirt, [('out_file', 'in_file')]),
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
