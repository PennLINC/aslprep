"""Workflows to apply changes to ASL data."""

from __future__ import annotations

import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe

from aslprep import config


def init_asl_cifti_resample_wf(
    *,
    metadata: dict,
    mem_gb: dict,
    fieldmap_id: str | None = None,
    omp_nthreads: int = 1,
    name: str = 'asl_cifti_resample_wf',
) -> pe.Workflow:
    """Resample an ASL series to a CIFTI target space.

    This workflow collates a sequence of transforms to resample an ASL series in a single shot,
    including motion correction and fieldmap correction, if requested.

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
        BIDS metadata for BOLD file.
    mem_gb
    fieldmap_id
        Fieldmap identifier, if fieldmap correction is to be applied.
    omp_nthreads
        Maximum number of threads an individual process may use.
    name
        Name of workflow (default: ``bold_volumetric_resample_wf``)

    Inputs
    ------
    asl_file
        ASL series to resample.
    bold_ref_file
        Reference image to which BOLD series is aligned.
    target_ref_file
        Reference image defining the target space.
    target_mask
        Brain mask corresponding to ``target_ref_file``.
        This is used to define the field of view for the resampled BOLD series.
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
        further images into the BOLD series' space.
    """
    from fmriprep.workflows.bold.apply import init_bold_volumetric_resample_wf
    from fmriprep.workflows.bold.resampling import (
        init_bold_grayords_wf,
        init_goodvoxels_bold_mask_wf,
    )

    from aslprep.workflows.asl.resampling import init_bold_fsLR_resampling_wf

    workflow = pe.Workflow(name=name)

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
                'coreg_aslref',
                'white',
                'pial',
                'midthickness',
                'midthickness_fsLR',
                'sphere_reg_fsLR',
                'cortex_mask',
                'anat_ribbon',
            ],
        ),
        name='inputnode',
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['asl_cifti', 'cifti_metadata', 'goodvoxels_mask']),
        name='outputnode',
    )

    asl_MNI6_wf = init_bold_volumetric_resample_wf(
        metadata=metadata,
        fieldmap_id=fieldmap_id,
        omp_nthreads=omp_nthreads,
        mem_gb=mem_gb,
        jacobian='fmap-jacobian' not in config.workflow.ignore,
        name='asl_MNI6_wf',
    )

    asl_fsLR_resampling_wf = init_bold_fsLR_resampling_wf(
        grayord_density=config.workflow.cifti_output,
        omp_nthreads=omp_nthreads,
        mem_gb=mem_gb['resampled'],
        name='asl_fsLR_resampling_wf',
    )

    if config.workflow.project_goodvoxels:
        goodvoxels_bold_mask_wf = init_goodvoxels_bold_mask_wf(mem_gb['resampled'])

        workflow.connect([
            (inputnode, goodvoxels_bold_mask_wf, [
                ('asl_anat', 'inputnode.bold_file'),
                ('anat_ribbon', 'inputnode.anat_ribbon'),
            ]),
            (goodvoxels_bold_mask_wf, asl_fsLR_resampling_wf, [
                ('outputnode.goodvoxels_mask', 'inputnode.volume_roi'),
            ]),
            (goodvoxels_bold_mask_wf, outputnode, [
                ('outputnode.goodvoxels_mask', 'goodvoxels_mask'),
            ]),
        ])  # fmt:skip

        asl_fsLR_resampling_wf.__desc__ += """\
            A "goodvoxels" mask was applied during volume-to-surface sampling in fsLR space,
            excluding voxels whose time-series have a locally high coefficient of variation.
            """

    asl_grayords_wf = init_bold_grayords_wf(
        grayord_density=config.workflow.cifti_output,
        mem_gb=mem_gb['resampled'],
        repetition_time=metadata['RepetitionTime'],
        name='asl_grayords_wf',
    )

    workflow.connect([
        # Resample BOLD to MNI152NLin6Asym, may duplicate asl_std_wf above
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
        ]),
        # Resample T1w-space BOLD to fsLR surfaces
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
