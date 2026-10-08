"""Workflows to apply changes to ASL data."""

from __future__ import annotations

import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe

from aslprep import config


def init_asl_volumetric_resample_wf(
    *,
    metadata: dict,
    mem_gb: dict[str, float],
    jacobian: bool,
    gradwarp: bool = False,
    gradwarp_jacobian: bool = True,
    fallback_total_readout_time: str | float | None = None,
    fieldmap_id: str | None = None,
    omp_nthreads: int = 1,
    name: str = 'asl_volumetric_resample_wf',
) -> pe.Workflow:
    """Resample an ASL series to a volumetric target space.

    This workflow collates a sequence of transforms to resample an ASL series in a single shot,
    including motion correction, fieldmap correction, and gradient nonlinearity correction,
    if requested.

    This is fMRIPrep's ``init_bold_volumetric_resample_wf`` (fMRIPrep 25.2),
    with gradient nonlinearity correction added.
    Its input and output names are kept for compatibility.

    .. workflow::

        from aslprep.workflows.asl.apply import init_asl_volumetric_resample_wf
        wf = init_asl_volumetric_resample_wf(
            metadata={
                'RepetitionTime': 2.0,
                'PhaseEncodingDirection': 'j-',
                'TotalReadoutTime': 0.03
            },
            mem_gb={'resampled': 1},
            jacobian=True,
            gradwarp=True,
            fieldmap_id='my_fieldmap',
        )

    Parameters
    ----------
    metadata
        BIDS metadata for the ASL file.
    mem_gb
        Memory estimates for the ASL series, in GB. The ``'resampled'`` key is used.
    jacobian
        Whether to apply the Jacobian determinant of the fieldmap.
    gradwarp
        Whether to apply a gradient nonlinearity displacement field (``gradwarp_field``).
    gradwarp_jacobian
        Whether to modulate intensities by the Jacobian determinant of ``gradwarp_field``.
    fallback_total_readout_time
        Total readout time to use if it cannot be determined from the metadata
        (a number, or ``'estimated'``).
    fieldmap_id
        Fieldmap identifier, if fieldmap correction is to be applied.
    omp_nthreads
        Maximum number of threads an individual process may use.
    name
        Name of workflow (default: ``asl_volumetric_resample_wf``)

    Inputs
    ------
    bold_file
        ASL series to resample.
    bold_ref_file
        Reference image to which the ASL series is aligned.
    target_ref_file
        Reference image defining the target space.
    target_mask
        Brain mask corresponding to ``target_ref_file``.
        This is used to define the field of view for the resampled ASL series.
    motion_xfm
        List of affine transforms aligning each volume to ``bold_ref_file``.
        If undefined, no motion correction is performed.
    boldref2fmap_xfm
        Affine transform from ``bold_ref_file`` to the fieldmap reference image.
    fmap_ref
        Fieldmap reference image defining the valid field of view for the fieldmap.
    fmap_coeff
        B-Spline coefficients for the fieldmap.
    fmap_id
        Fieldmap identifier, to select correct fieldmap in case there are multiple.
    gradwarp_field
        Gradient nonlinearity displacement field, on the ASL reference grid.
    boldref2anat_xfm
        Affine transform from ``bold_ref_file`` to the anatomical reference image.
    anat2std_xfm
        Affine transform from the anatomical reference image to standard space.
        Leave undefined to resample to anatomical reference space.

    Outputs
    -------
    bold_file
        The ``bold_file`` input, resampled to ``target_ref_file`` space.
    resampling_reference
        An empty reference image with the correct affine and header for resampling
        further images into the ASL series' space.

    Returns
    -------
    workflow : nipype.pipeline.engine.Workflow
        The workflow.
    """
    from fmriprep.interfaces.resampling import DistortionParameters, ReconstructFieldmap
    from niworkflows.interfaces.nibabel import GenerateSamplingReference
    from niworkflows.interfaces.utility import KeySelect

    from aslprep.interfaces.resampling import ResampleSeries

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'bold_file',
                'bold_ref_file',
                'target_ref_file',
                'target_mask',
                # HMC
                'motion_xfm',
                # SDC
                'boldref2fmap_xfm',
                'fmap_ref',
                'fmap_coeff',
                'fmap_id',
                # Gradient nonlinearity
                'gradwarp_field',
                # Anatomical
                'boldref2anat_xfm',
                # Template
                'anat2std_xfm',
                # Entity for selecting target resolution
                'resolution',
            ],
        ),
        name='inputnode',
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['bold_file', 'resampling_reference']),
        name='outputnode',
    )

    gen_ref = pe.Node(GenerateSamplingReference(), name='gen_ref', mem_gb=0.3)

    boldref2target = pe.Node(niu.Merge(2), name='boldref2target', run_without_submitting=True)
    bold2target = pe.Node(niu.Merge(2), name='bold2target', run_without_submitting=True)
    resample = pe.Node(
        ResampleSeries(
            jacobian=jacobian,
            gradwarp_jacobian=gradwarp_jacobian,
        ),
        name='resample',
        n_procs=omp_nthreads,
        mem_gb=mem_gb['resampled'],
    )

    workflow.connect([
        (inputnode, gen_ref, [
            ('bold_ref_file', 'moving_image'),
            ('target_ref_file', 'fixed_image'),
            ('target_mask', 'fov_mask'),
            (('resolution', _is_native), 'keep_native'),
        ]),
        (inputnode, boldref2target, [
            ('boldref2anat_xfm', 'in1'),
            ('anat2std_xfm', 'in2'),
        ]),
        (inputnode, bold2target, [('motion_xfm', 'in1')]),
        (inputnode, resample, [('bold_file', 'in_file')]),
        (gen_ref, resample, [('out_file', 'ref_file')]),
        (boldref2target, bold2target, [('out', 'in2')]),
        (bold2target, resample, [('out', 'transforms')]),
        (gen_ref, outputnode, [('out_file', 'resampling_reference')]),
        (resample, outputnode, [('out_file', 'bold_file')]),
    ])  # fmt:skip

    if gradwarp:
        workflow.connect([(inputnode, resample, [('gradwarp_field', 'gradwarp_field')])])

    if not fieldmap_id:
        return workflow

    fmap_select = pe.Node(
        KeySelect(fields=['fmap_ref', 'fmap_coeff'], key=fieldmap_id),
        name='fmap_select',
        run_without_submitting=True,
    )
    distortion_params = pe.Node(
        DistortionParameters(
            metadata=metadata,
            fallback=fallback_total_readout_time,
        ),
        name='distortion_params',
        run_without_submitting=True,
    )
    fmap2target = pe.Node(niu.Merge(2), name='fmap2target', run_without_submitting=True)
    inverses = pe.Node(
        niu.Function(function=_gen_inverses),
        name='inverses',
        run_without_submitting=True,
    )

    fmap_recon = pe.Node(ReconstructFieldmap(), name='fmap_recon', mem_gb=1)

    workflow.connect([
        (inputnode, fmap_select, [
            ('fmap_ref', 'fmap_ref'),
            ('fmap_coeff', 'fmap_coeff'),
            ('fmap_id', 'keys'),
        ]),
        (inputnode, distortion_params, [('bold_file', 'in_file')]),
        (inputnode, fmap2target, [('boldref2fmap_xfm', 'in1')]),
        (gen_ref, fmap_recon, [('out_file', 'target_ref_file')]),
        (boldref2target, fmap2target, [('out', 'in2')]),
        (boldref2target, inverses, [('out', 'inlist')]),
        (fmap_select, fmap_recon, [
            ('fmap_coeff', 'in_coeffs'),
            ('fmap_ref', 'fmap_ref_file'),
        ]),
        (fmap2target, fmap_recon, [('out', 'transforms')]),
        (inverses, fmap_recon, [('out', 'inverse')]),
        # Inject fieldmap correction into resample node
        (distortion_params, resample, [
            ('readout_time', 'ro_time'),
            ('pe_direction', 'pe_dir'),
        ]),
        (fmap_recon, resample, [('out_file', 'fieldmap')]),
    ])  # fmt:skip

    return workflow


def _gen_inverses(inlist: list) -> list[bool]:
    """Create a list indicating the first transform should be inverted.

    Parameters
    ----------
    inlist : list or str or None
        The transforms that follow the inverted one.

    Returns
    -------
    list of bool
        True for the first transform, and False for each transform in ``inlist``.
    """
    from niworkflows.utils.connections import listify

    if not inlist:
        return [True]
    return [True] + [False] * len(listify(inlist))


def _is_native(value):
    """Check whether a resolution entity requests the native resolution.

    Parameters
    ----------
    value : str or None
        The ``resolution`` input.

    Returns
    -------
    bool
        True if ``value`` is ``'native'``.
    """
    return value == 'native'


def init_asl_cifti_resample_wf(
    *,
    asl_file: str,
    metadata: dict,
    mem_gb: dict,
    fieldmap_id: str | None = None,
    jacobian: bool = False,
    gradwarp: bool = False,
    gradwarp_jacobian: bool = True,
    omp_nthreads: int = 1,
    name: str = 'asl_cifti_resample_wf',
) -> pe.Workflow:
    """Resample an ASL series to a CIFTI target space.

    This workflow collates a sequence of transforms to resample an ASL series in a single shot,
    including motion correction, fieldmap correction, and gradient nonlinearity correction,
    if requested.

    This is an ASLPrep-specific workflow collecting steps from
    fmriprep.workflows.bold.base.init_bold_wf.

    .. workflow::

        from aslprep.workflows.asl.apply import init_asl_cifti_resample_wf

        wf = init_asl_cifti_resample_wf(
            metadata={
                "RepetitionTime": 2.0,
                "PhaseEncodingDirection": "j-",
                "TotalReadoutTime": 0.03
            },
            mem_gb={
                "resampled": 1,
            },
            fieldmap_id="my_fieldmap",
        )

    Parameters
    ----------
    metadata
        BIDS metadata for ASL file.
    mem_gb
    fieldmap_id
        Fieldmap identifier, if fieldmap correction is to be applied.
    jacobian
        Whether to apply the Jacobian determinant of the fieldmap.
    gradwarp
        Whether to apply a gradient nonlinearity displacement field (``gradwarp_field``).
    gradwarp_jacobian
        Whether to modulate intensities by the Jacobian determinant of ``gradwarp_field``.
    omp_nthreads
        Maximum number of threads an individual process may use.
    name
        Name of workflow (default: ``asl_cifti_resample_wf``)

    Inputs
    ------
    asl_file
        ASL series to resample.
    bold_ref_file
        Reference image to which ASL series is aligned.
    target_ref_file
        Reference image defining the target space.
    target_mask
        Brain mask corresponding to ``target_ref_file``.
        This is used to define the field of view for the resampled ASL series.
    motion_xfm
        List of affine transforms aligning each volume to ``bold_ref_file``.
        If undefined, no motion correction is performed.
    boldref2fmap_xfm
        Affine transform from ``bold_ref_file`` to the fieldmap reference image.
    fmap_ref
        Fieldmap reference image defining the valid field of view for the fieldmap.
    fmap_coeff
        B-Spline coefficients for the fieldmap.
    fmap_id
        Fieldmap identifier, to select correct fieldmap in case there are multiple.
    boldref2anat_xfm
        Affine transform from ``bold_ref_file`` to the anatomical reference image.
    anat2std_xfm
        Affine transform from the anatomical reference image to standard space.
        Leave undefined to resample to anatomical reference space.

    Outputs
    -------
    bold_file
        The ``bold_file`` input, resampled to ``target_ref_file`` space.
    resampling_reference
        An empty reference image with the correct affine and header for resampling
        further images into the ASL series' space.
    """
    from fmriprep.workflows.bold.resampling import (
        init_bold_fsLR_resampling_wf,
        init_bold_grayords_wf,
    )
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                # Raw ASL file (asl_minimal)
                'asl_file',
                # ASL file in T1w space
                'asl_anat',
                # Other inputs
                'mni6_mask',
                'aslref2fmap_xfm',
                'aslref2anat_xfm',
                'anat2mni6_xfm',
                'fmap_ref',
                'fmap_coeff',
                'fmap_id',
                'motion_xfm',
                'gradwarp_field',
                'coreg_aslref',
                'white',
                'pial',
                'midthickness',
                'midthickness_fsLR',
                'sphere_reg_fsLR',
                'cortex_mask',
                # Pre-computed goodvoxels mask in T1w space. May be Undefined.
                'goodvoxels_mask',
            ],
        ),
        name='inputnode',
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['asl_cifti', 'cifti_metadata']),
        name='outputnode',
    )

    asl_MNI6_wf = init_asl_volumetric_resample_wf(
        metadata=metadata,
        fieldmap_id=fieldmap_id,
        jacobian=jacobian,
        gradwarp=gradwarp,
        gradwarp_jacobian=gradwarp_jacobian,
        fallback_total_readout_time=config.workflow.fallback_total_readout_time,
        omp_nthreads=omp_nthreads,
        mem_gb=mem_gb,
        name='asl_MNI6_wf',
    )

    asl_fsLR_resampling_wf = init_bold_fsLR_resampling_wf(
        grayord_density=config.workflow.cifti_output,
        omp_nthreads=omp_nthreads,
        mem_gb=mem_gb['resampled'],
        name='asl_fsLR_resampling_wf',
    )

    if config.workflow.project_goodvoxels:
        workflow.connect([
            (inputnode, asl_fsLR_resampling_wf, [('goodvoxels_mask', 'inputnode.volume_roi')]),
        ])  # fmt:skip

    asl_grayords_wf = init_bold_grayords_wf(
        grayord_density=config.workflow.cifti_output,
        mem_gb=mem_gb['resampled'],
        repetition_time=metadata['RepetitionTime'],
        name='asl_grayords_wf',
    )

    workflow.connect([
        # Resample ASL to MNI152NLin6Asym, may duplicate asl_std_wf above
        (inputnode, asl_MNI6_wf, [
            ('mni6_mask', 'inputnode.target_ref_file'),
            ('mni6_mask', 'inputnode.target_mask'),
            ('anat2mni6_xfm', 'inputnode.anat2std_xfm'),
            ('fmap_ref', 'inputnode.fmap_ref'),
            ('fmap_coeff', 'inputnode.fmap_coeff'),
            ('fmap_id', 'inputnode.fmap_id'),
            ('asl_file', 'inputnode.bold_file'),
            ('motion_xfm', 'inputnode.motion_xfm'),
            ('coreg_aslref', 'inputnode.bold_ref_file'),
            ('aslref2fmap_xfm', 'inputnode.boldref2fmap_xfm'),
            ('aslref2anat_xfm', 'inputnode.boldref2anat_xfm'),
            ('gradwarp_field', 'inputnode.gradwarp_field'),
        ]),
        # Resample T1w-space ASL to fsLR surfaces
        (inputnode, asl_fsLR_resampling_wf, [
            ('asl_anat', 'inputnode.bold_file'),
            ('white', 'inputnode.white'),
            ('pial', 'inputnode.pial'),
            ('midthickness', 'inputnode.midthickness'),
            ('midthickness_fsLR', 'inputnode.midthickness_fsLR'),
            ('sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR'),
            ('cortex_mask', 'inputnode.cortex_mask'),
        ]),
        (asl_MNI6_wf, asl_grayords_wf, [('outputnode.bold_file', 'inputnode.bold_std')]),
        (asl_fsLR_resampling_wf, asl_grayords_wf, [
            ('outputnode.bold_fsLR', 'inputnode.bold_fsLR'),
        ]),
        (asl_grayords_wf, outputnode, [
            ('outputnode.cifti_bold', 'asl_cifti'),
            ('outputnode.cifti_metadata', 'cifti_metadata'),
        ]),
    ])  # fmt:skip

    return workflow
