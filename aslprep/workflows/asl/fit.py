# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# Copyright 2023 The NiPreps Developers <nipreps@gmail.com>
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
"""Fit workflows for ASLPrep."""

import os
import re
import typing as ty

import bids
import nibabel as nb
import numpy as np
from fmriprep.interfaces.resampling import DistortionParameters, ReconstructFieldmap
from fmriprep.utils.bids import extract_entities
from fmriprep.workflows.bold import outputs as output_workflows
from fmriprep.workflows.bold.reference import init_validation_and_dummies_wf
from fmriprep.workflows.bold.registration import init_bold_reg_wf
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.interfaces.header import ValidateImage
from niworkflows.interfaces.nitransforms import ConcatenateXFMs
from niworkflows.interfaces.utility import KeySelect
from sdcflows.workflows.apply.registration import init_coeff2epi_wf

# ASL workflows
from aslprep import config
from aslprep.interfaces.bids import OverrideDerivativesDataSink
from aslprep.interfaces.hmc import IdentityTransforms
from aslprep.interfaces.reports import FunctionalSummary
from aslprep.interfaces.resampling import ResampleSeries
from aslprep.interfaces.utility import ReduceASLFiles
from aslprep.utils.asl import select_processing_target
from aslprep.utils.gradwarp import GradwarpPlan
from aslprep.workflows.asl.hmc import init_asl_hmc_wf
from aslprep.workflows.asl.outputs import init_asl_fit_reports_wf, init_ds_aslref_wf
from aslprep.workflows.asl.reference import (
    init_anat_aslmask_wf,
    init_enhance_and_skullstrip_asl_wf,
    init_raw_aslref_wf,
    init_skullstrip_aslref_wf,
    init_synthstrip_aslref_wf,
)


def get_sbrefs(
    asl_file: str,
    entity_overrides: dict[str, ty.Any],
    layout: bids.BIDSLayout,
) -> list[str]:
    """Find single-band reference(s) associated with ASL file.

    Parameters
    ----------
    asl_file
        Absolute path to the ASL file
    entity_overrides
        Query parameters to override defaults
    layout
        :class:`~bids.layout.BIDSLayout` to query

    Returns
    -------
    sbref_files
        List of absolute paths to sbref files associated with the input ASL file
    """
    entities = extract_entities(asl_file)
    entities.update(suffix='sbref', extension=['.nii', '.nii.gz'])
    entities.update(entity_overrides)

    return layout.get(return_type='file', **entities)


def init_asl_fit_wf(
    *,
    asl_file: str,
    aslcontext: str,
    m0scan: str | None,
    use_ge: bool,
    precomputed: dict | None = None,
    fieldmap_id: str | None = None,
    jacobian: bool = False,
    gradwarp_plan: GradwarpPlan | None = None,
    m0scan_gradwarp_plan: GradwarpPlan | None = None,
    asl2anat_init: str | None = None,
    coreg_r1map: dict | None = None,
    omp_nthreads: int = 1,
    name: str = 'asl_fit_wf',
) -> pe.Workflow:
    """Control the minimal estimation steps for functional preprocessing.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep import config
            from aslprep.workflows.asl.fit import init_asl_fit_wf

            with mock_config():
                asl_file = (
                    config.execution.bids_dir / "sub-01" / "perf" /
                    "sub-01_asl.nii.gz"
                )
                wf = init_asl_fit_wf(
                    asl_file=str(asl_file),
                    aslcontext=str(asl_file).replace('.nii.gz', 'context.tsv'),
                    m0scan=None,
                    use_ge=False,
                )

    Parameters
    ----------
    asl_file
        Path to ASL NIfTI file.
    aslcontext
        Path to ASL context file.
    m0scan
        Path to M0 NIfTI file, if available.
    use_ge
        If True, the M0 scan (when available) will be prioritized as the volume type for the
        reference image, as GE deltam volumes exhibit extreme background noise, and head-motion
        correction is skipped if the series has no label-control pairs.
    precomputed
        Dictionary containing precomputed derivatives to reuse, if possible.
    fieldmap_id
        ID of the fieldmap to use to correct this ASL series. If :obj:`None`,
        no correction will be applied.
    jacobian
        Whether to apply the Jacobian determinant of the fieldmap during resampling.
    gradwarp_plan
        Gradient nonlinearity correction to apply, if any.
        If it calls for a correction, a displacement field is computed on the grid of
        ``hmc_aslref`` and applied to the coregistration reference and the fieldmap
        registration target.
    m0scan_gradwarp_plan
        Gradient nonlinearity correction to apply to the separate M0 scan, if any,
        resolved from the M0 scan's own metadata.
    asl2anat_init
        Anatomical image to use as the initial target for ASL-to-anatomical
        coregistration (``'t1w'`` or ``'t2w'``).
        If :obj:`None`, ``config.workflow.asl2anat_init`` is used,
        with ``'auto'`` resolved to ``'t1w'``.
    coreg_r1map
        A T1 or R1 map on the grid of the separate M0 scan, as returned by
        :func:`~aslprep.utils.bids.collect_m0scan_t1map`, or :obj:`None`.
        If the separate M0 scan is the reference volume type, and has one volume, the map's R1,
        moved into the frame of the coregistration reference with the M0 scan's transform and
        masked with the ASL brain mask, is registered to the anatomical image in place of the
        coregistration reference, with the boundary-based registration's contrast reversed.
        The coregistration reference itself, and the transform's name, do not change.

    Inputs
    ------
    asl_file
        ASL series NIfTI file
    t1w_preproc
        Bias-corrected structural template image
    t1w_mask
        Mask of the skull-stripped template image
    t1w_dseg
        Segmentation of preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF)
    anat2std_xfm
        List of transform files, collated with templates
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    fsnative2t1w_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed subject space to T1w
    fmap_id
        Unique identifiers to select fieldmap files
    fmap
        List of estimated fieldmaps (collated with fmap_id)
    fmap_ref
        List of fieldmap reference files (collated with fmap_id)
    fmap_coeff
        List of lists of spline coefficient files (collated with fmap_id)
    fmap_mask
        List of fieldmap masks (collated with fmap_id)
    sdc_method
        List of fieldmap correction method names (collated with fmap_id)

    Outputs
    -------
    hmc_aslref
        ASL reference image used for head motion correction.
        Minimally processed to ensure consistent contrast with ASL series.
    coreg_aslref
        ASL reference image used for coregistration. Contrast-enhanced
        and fieldmap-corrected for greater anatomical fidelity, and aligned
        with ``hmc_aslref``.
    asl_mask
        Mask of ``coreg_aslref``.
    motion_xfm
        Affine transforms from each ASL volume to ``hmc_aslref``, written
        as concatenated ITK affine transforms.
    aslref2anat_xfm
        Affine transform mapping from ASL reference space to the anatomical
        space.
    aslref2fmap_xfm
        Affine transform mapping from ASL reference space to the fieldmap
        space, if applicable.
    gradwarp_field
        Gradient nonlinearity displacement field on the ``hmc_aslref`` grid, if applicable.
    m0scan_gradwarp_field
        Gradient nonlinearity displacement field for the separate M0 scan,
        on the ``hmc_aslref`` grid, if applicable.
    skullstrip_qc
        With ``--skull-strip-asl anatomical`` (and no precomputed mask), the quality-control
        measures of the anatomical brain mask's check against SynthStrip.
    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    from aslprep.utils.misc import estimate_asl_mem_usage
    from aslprep.workflows.asl.gradwarp import gradwarp_boilerplate, init_gradwarp_wf

    if precomputed is None:
        precomputed = {}
    if asl2anat_init is None:
        asl2anat_init = config.workflow.asl2anat_init
    if asl2anat_init == 'auto':
        # init_single_subject_wf resolves 'auto' based on T2w availability.
        # Without subject data, fall back to the T1w.
        asl2anat_init = 't1w'
    layout = config.execution.layout
    bids_filters = config.execution.get().get('bids_filters', {})

    # Collect asl and sbref files.
    # sbrefs aren't supported for ASL data, but I think that might change in the future.
    sbref_files = get_sbrefs(
        asl_file,
        entity_overrides=bids_filters.get('sbref', {}),
        layout=layout,
    )
    basename = os.path.basename(asl_file)
    sbref_msg = f'No single-band-reference found for {basename}.'
    if sbref_files and 'sbref' in config.workflow.ignore:
        sbref_msg = f'Single-band reference file(s) found for {basename} and ignored.'
        sbref_files = []
    elif sbref_files:
        sbref_msg = (
            'Using single-band reference file(s) '
            f'{", ".join(os.path.basename(f) for f in sbref_files)}.'
        )
    config.loggers.workflow.info(sbref_msg)

    # Get metadata from ASL file(s)
    metadata = layout.get_metadata(asl_file)
    # Patch RepetitionTimePreparation into RepetitionTime,
    # for the sake of ASL-based interfaces and workflows.
    # This value shouldn't be used for anything except figures and reportlets.
    metadata['RepetitionTime'] = metadata.get(
        'RepetitionTime',
        np.mean(metadata['RepetitionTimePreparation']),
    )

    orientation = ''.join(nb.aff2axcodes(nb.load(asl_file).affine))

    _, mem_gb = estimate_asl_mem_usage(asl_file)

    hmc_aslref = precomputed.get('hmc_aslref')
    coreg_aslref = precomputed.get('coreg_aslref')
    # Can contain
    #  1) aslref2fmap
    #  2) aslref2anat
    #  3) hmc
    transforms = precomputed.get('transforms', {})
    hmc_xforms = transforms.get('hmc')
    aslref2fmap_xform = transforms.get('aslref2fmap')
    aslref2anat_xform = transforms.get('aslref2anat')
    m0scan2aslref_xform = transforms.get('m0scan2aslref')
    aslref_mask = precomputed.get('aslref_mask')

    gradwarp = gradwarp_plan is not None and gradwarp_plan.warp_dim is not None
    m0scan_gradwarp = bool(
        m0scan and m0scan_gradwarp_plan is not None and m0scan_gradwarp_plan.warp_dim is not None
    )
    gradwarp_jacobian = 'gradwarp-jacobian' not in (config.workflow.ignore or [])
    if gradwarp:
        # These depend on the corrected geometry, and derivatives do not record whether (or how)
        # they were corrected, so they are recomputed rather than risk mixing geometries.
        # The HMC reference and transforms are estimated on raw data, so they remain reusable.
        reused = {
            'coregistration reference': coreg_aslref,
            'ASL reference mask': aslref_mask,
            'ASL-to-anatomical transform': aslref2anat_xform,
            'ASL-to-fieldmap transform': aslref2fmap_xform,
        }
        reused = [name for name, value in reused.items() if value]
        if reused:
            config.loggers.workflow.warning(
                'Recomputing the precomputed %s, as gradient nonlinearity correction is enabled.',
                ', '.join(reused),
            )
        coreg_aslref = aslref_mask = aslref2anat_xform = aslref2fmap_xform = None

    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'asl_file',
                'aslcontext',
                # Fieldmap registration
                'fmap',
                'fmap_ref',
                'fmap_coeff',
                'fmap_mask',
                'fmap_id',
                'sdc_method',
                # Anatomical coregistration
                't1w_preproc',
                't1w_mask',
                't1w_dseg',
                'subjects_dir',
                'subject_id',
                'fsnative2t1w_xfm',
            ],
        ),
        name='inputnode',
    )
    inputnode.inputs.asl_file = asl_file
    inputnode.inputs.aslcontext = aslcontext
    inputnode.inputs.m0scan = m0scan

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'hmc_aslref',
                'coreg_aslref',
                'asl_mask',
                'motion_xfm',
                'aslref2anat_xfm',
                'aslref2fmap_xfm',
                'm0scan2aslref_xfm',
                'gradwarp_field',
                'm0scan_gradwarp_field',
                'skullstrip_qc',
            ],
        ),
        name='outputnode',
    )

    # If all derivatives exist, inputnode could go unconnected, so add explicitly
    workflow.add_nodes([inputnode])

    reference_volume_type = select_reference_volume_type(
        aslcontext=aslcontext,
        metadata=metadata,
        prioritize_m0=use_ge,
    )
    skip_hmc = skips_motion_correction(aslcontext, use_ge)
    use_r1 = use_coregistration_r1(
        coreg_r1map,
        m0scan=m0scan,
        reference_volume_type=reference_volume_type,
        aslref2anat_xform=aslref2anat_xform,
    )

    hmcref_buffer = pe.Node(
        niu.IdentityInterface(fields=['aslref', 'asl_file']),
        name='hmcref_buffer',
    )
    fmapref_buffer = pe.Node(niu.Function(function=_select_ref), name='fmapref_buffer')
    hmc_buffer = pe.Node(niu.IdentityInterface(fields=['hmc_xforms']), name='hmc_buffer')
    m0scanreg_buffer = pe.Node(
        niu.IdentityInterface(fields=['m0scan2aslref_xfm']),
        name='m0scanreg_buffer',
    )
    fmapreg_buffer = pe.Node(
        niu.IdentityInterface(fields=['aslref2fmap_xfm']),
        name='fmapreg_buffer',
    )
    regref_buffer = pe.Node(
        niu.IdentityInterface(fields=['aslref']),
        name='regref_buffer',
    )
    aslmask_buffer = pe.Node(
        niu.IdentityInterface(fields=['aslmask']),
        name='aslmask_buffer',
    )

    if hmc_aslref:
        hmcref_buffer.inputs.aslref = hmc_aslref
        config.loggers.workflow.debug('Reusing motion correction reference: %s', hmc_aslref)
    if hmc_xforms:
        hmc_buffer.inputs.hmc_xforms = hmc_xforms
        # Only transforms whose sidecar records the requested --hmc-method reach this point
        # (merge_functional_derivatives drops the others).
        config.loggers.workflow.debug('Reusing motion correction transforms: %s', hmc_xforms)
    if m0scan2aslref_xform:
        m0scanreg_buffer.inputs.m0scan2aslref_xfm = m0scan2aslref_xform
        config.loggers.workflow.debug(
            'Reusing M0 scan to ASL reference transform(s): %s',
            m0scan2aslref_xform,
        )
    if aslref2fmap_xform:
        fmapreg_buffer.inputs.aslref2fmap_xfm = aslref2fmap_xform
        config.loggers.workflow.debug('Reusing ASL-to-fieldmap transform: %s', aslref2fmap_xform)
    if coreg_aslref:
        regref_buffer.inputs.aslref = coreg_aslref
        config.loggers.workflow.debug('Reusing coregistration reference: %s', coreg_aslref)
    if aslref_mask:
        aslmask_buffer.inputs.aslmask = aslref_mask
        config.loggers.workflow.debug('Reusing ASL reference mask: %s', aslref_mask)
    fmapref_buffer.inputs.sbref_files = sbref_files

    summary = pe.Node(
        FunctionalSummary(
            distortion_correction='None',  # Can override with connection
            registration=(
                'Precomputed'
                if aslref2anat_xform
                else 'FreeSurfer'
                if config.workflow.run_reconall
                else 'FSL'
            ),
            registration_dof=config.workflow.asl2anat_dof,
            registration_init=asl2anat_init,
            pe_direction=metadata.get('PhaseEncodingDirection'),
            tr=metadata['RepetitionTime'],
            orientation=orientation,
        ),
        name='summary',
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    asl_fit_reports_wf = init_asl_fit_reports_wf(
        source_file=asl_file,
        # TODO: Enable sdc report even if we find coregref
        sdc_correction=fieldmap_id is not None,
        # The M0 report shows the M0 registered here; with a reused M0 transform there is none
        # (the report node failed for want of its 'before' image and stalled the run).
        separate_m0scan=m0scan is not None and not m0scan2aslref_xform,
        freesurfer=config.workflow.run_reconall,
        output_dir=config.execution.aslprep_dir,
    )

    workflow.connect([
        (hmcref_buffer, fmapref_buffer, [('aslref', 'aslref_files')]),
        (hmcref_buffer, outputnode, [('aslref', 'hmc_aslref')]),
        (regref_buffer, outputnode, [('aslref', 'coreg_aslref')]),
        (aslmask_buffer, outputnode, [('aslmask', 'asl_mask')]),
        (fmapreg_buffer, outputnode, [('aslref2fmap_xfm', 'aslref2fmap_xfm')]),
        (hmc_buffer, outputnode, [('hmc_xforms', 'motion_xfm')]),
        (m0scanreg_buffer, outputnode, [('m0scan2aslref_xfm', 'm0scan2aslref_xfm')]),
        (inputnode, asl_fit_reports_wf, [
            ('t1w_preproc', 'inputnode.t1w_preproc'),
            # May not need all of these
            ('t1w_mask', 'inputnode.t1w_mask'),
            ('t1w_dseg', 'inputnode.t1w_dseg'),
            ('subjects_dir', 'inputnode.subjects_dir'),
            ('subject_id', 'inputnode.subject_id'),
        ]),
        (outputnode, asl_fit_reports_wf, [
            ('coreg_aslref', 'inputnode.coreg_aslref'),
            ('asl_mask', 'inputnode.asl_mask'),
            ('aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
        ]),
        (summary, asl_fit_reports_wf, [('out_report', 'inputnode.summary_report')]),
    ])  # fmt:skip

    # Stage 1: Generate motion correction aslref
    hmc_aslref_source_buffer = pe.Node(
        niu.IdentityInterface(fields=['in_file']),
        name='hmc_aslref_source_buffer',
    )
    if not hmc_aslref:
        config.loggers.workflow.info('Stage 1: Adding HMC aslref workflow')
        hmc_aslref_wf = init_raw_aslref_wf(
            name='hmc_aslref_wf',
            asl_file=asl_file,
            reference_volume_type=reference_volume_type,
            m0scan=(metadata['M0Type'] == 'Separate'),
            use_ge=use_ge,
        )
        hmc_aslref_wf.inputs.inputnode.m0scan = m0scan

        workflow.connect([(inputnode, hmc_aslref_wf, [('aslcontext', 'inputnode.aslcontext')])])

        ds_hmc_aslref_wf = init_ds_aslref_wf(
            source_file=asl_file,
            bids_root=layout.root,
            output_dir=config.execution.aslprep_dir,
            desc='hmc',
            name='ds_hmc_aslref_wf',
        )
        ds_hmc_aslref_wf.inputs.inputnode.source_files = [asl_file]

        workflow.connect([
            (hmc_aslref_wf, ds_hmc_aslref_wf, [('outputnode.aslref', 'inputnode.aslref')]),
            (ds_hmc_aslref_wf, hmcref_buffer, [('outputnode.aslref', 'aslref')]),
            (hmc_aslref_wf, hmcref_buffer, [('outputnode.asl_file', 'asl_file')]),
            (hmc_aslref_wf, asl_fit_reports_wf, [
                ('outputnode.validation_report', 'inputnode.validation_report'),
            ]),
            (ds_hmc_aslref_wf, hmc_aslref_source_buffer, [('outputnode.aslref', 'in_file')]),
        ])  # fmt:skip
    else:
        config.loggers.workflow.info('Found HMC aslref - skipping Stage 1')

        validation_and_dummies_wf = init_validation_and_dummies_wf(bold_file=asl_file)

        workflow.connect([
            (validation_and_dummies_wf, hmcref_buffer, [('outputnode.bold_file', 'asl_file')]),
            (validation_and_dummies_wf, asl_fit_reports_wf, [
                ('outputnode.validation_report', 'inputnode.validation_report'),
            ]),
            (hmcref_buffer, hmc_aslref_source_buffer, [('aslref', 'in_file')]),
        ])  # fmt:skip

    # Stage 1b: Gradient nonlinearity displacement field, on the HMC reference's (raw) grid.
    # It is cheap and deterministic, so it is always recomputed rather than reused.
    if gradwarp:
        config.loggers.workflow.info(
            'Stage 1b: Adding gradient nonlinearity correction (%s, from %s)',
            gradwarp_plan.warp_dim,
            'command line' if gradwarp_plan.basis == 'forced' else 'ImageType',
        )
        gradwarp_wf = init_gradwarp_wf(plan=gradwarp_plan, asl_file=asl_file)
        gradwarp_wf.__desc__ = gradwarp_boilerplate(gradwarp_plan, jacobian=gradwarp_jacobian)
        workflow.connect([
            (hmcref_buffer, gradwarp_wf, [('aslref', 'inputnode.ref_image')]),
            (gradwarp_wf, outputnode, [('outputnode.gradwarp_field', 'gradwarp_field')]),
        ])  # fmt:skip
    elif gradwarp_plan is not None:
        config.loggers.workflow.info(
            'Skipping gradient nonlinearity correction: '
            'the ASL data were corrected on the scanner (ImageType: DIS3D).'
        )
        workflow.__desc__ = (workflow.__desc__ or '') + gradwarp_boilerplate(
            gradwarp_plan, jacobian=False
        )

    # The separate M0 scan may have been corrected differently on the scanner, and may have
    # a different slice orientation, so it gets its own field. As it is resampled onto the ASL
    # reference grid before correction, the field is generated on that grid.
    if m0scan_gradwarp:
        m0scan_gradwarp_wf = init_gradwarp_wf(
            plan=m0scan_gradwarp_plan,
            asl_file=m0scan,
            report=False,
            name='m0scan_gradwarp_wf',
        )
        m0scan_gradwarp_wf.inputs.inputnode.slice_ref_image = m0scan
        workflow.connect([
            (hmcref_buffer, m0scan_gradwarp_wf, [('aslref', 'inputnode.ref_image')]),
            (m0scan_gradwarp_wf, outputnode, [
                ('outputnode.gradwarp_field', 'm0scan_gradwarp_field'),
            ]),
        ])  # fmt:skip
    elif m0scan and m0scan_gradwarp_plan is not None:
        config.loggers.workflow.info(
            'Skipping gradient nonlinearity correction of the M0 scan: '
            'it was corrected on the scanner (ImageType: DIS3D).'
        )

    # Reduce the ASL series to only include volumes that need to be processed.
    processing_target = pe.Node(
        niu.Function(
            function=select_processing_target,
            input_names=['aslcontext'],
            output_names=['processing_target'],
        ),
        name='processing_target',
    )

    reduce_asl_file = pe.Node(
        ReduceASLFiles(metadata=metadata),
        name='reduce_asl_file',
    )

    workflow.connect([
        (inputnode, processing_target, [('aslcontext', 'aslcontext')]),
        (inputnode, reduce_asl_file, [('aslcontext', 'aslcontext')]),
        (processing_target, reduce_asl_file, [('processing_target', 'processing_target')]),
        (hmcref_buffer, reduce_asl_file, [('asl_file', 'asl_file')]),
    ])  # fmt:skip

    # Stage 2: Estimate head motion
    if not hmc_xforms:
        with OverrideDerivativesDataSink(output_workflows):
            ds_hmc_wf = output_workflows.init_ds_hmc_wf(
                source_file=asl_file,
                bids_root=layout.root,
                output_dir=config.execution.aslprep_dir,
            )

        ds_hmc_wf.get_node('inputnode').inputs.source_files = [asl_file]
        # fMRIPrep will write out an orig-to-boldref transform to anat, so we need to overwrite
        # some fields.
        ds_hmc_wf.get_node('ds_xforms').inputs.datatype = 'perf'
        ds_hmc_wf.get_node('ds_xforms').inputs.to = 'aslref'
        # Provenance: reuse of these transforms is limited to runs with the same method.
        ds_hmc_wf.get_node('ds_xforms').inputs.meta_dict = {
            'HMCMethod': 'none' if skip_hmc else config.workflow.hmc_method
        }

    if not hmc_xforms and skip_hmc:
        config.loggers.workflow.info(
            'Stage 2: Skipping motion correction (GE-specific path, no label-control pairs)'
        )
        identity_xforms = pe.Node(IdentityTransforms(), name='identity_xforms')
        workflow.__desc__ = (workflow.__desc__ or '') + (
            '\nHead-motion correction was not performed, as the series has too few volumes '
            'and no label-control pairs.\n'
        )
        workflow.connect([
            (hmcref_buffer, identity_xforms, [('aslref', 'ref_file')]),
            (reduce_asl_file, identity_xforms, [('aslcontext', 'aslcontext')]),
            (identity_xforms, ds_hmc_wf, [('out_file', 'inputnode.xforms')]),
            (identity_xforms, hmc_buffer, [('out_file', 'hmc_xforms')]),
        ])  # fmt:skip
    elif not hmc_xforms:
        config.loggers.workflow.info('Stage 2: Adding motion correction workflow')
        asl_hmc_wf = init_asl_hmc_wf(
            name='asl_hmc_wf',
            mem_gb=mem_gb['filesize'],
            omp_nthreads=omp_nthreads,
        )
        workflow.connect([
            (hmcref_buffer, asl_hmc_wf, [('aslref', 'inputnode.raw_ref_image')]),
            (reduce_asl_file, asl_hmc_wf, [
                ('asl_file', 'inputnode.asl_file'),
                ('aslcontext', 'inputnode.aslcontext'),
                ('metadata', 'inputnode.metadata'),
            ]),
            (asl_hmc_wf, ds_hmc_wf, [('outputnode.xforms', 'inputnode.xforms')]),
            (asl_hmc_wf, hmc_buffer, [('outputnode.xforms', 'hmc_xforms')]),
        ])  # fmt:skip
    else:
        config.loggers.workflow.info('Found motion correction transforms - skipping Stage 2')

    # Stage 2b: Register M0 scan to aslref
    if m0scan and not m0scan2aslref_xform:
        from aslprep.workflows.asl.hmc import init_m0scan_hmc_wf

        config.loggers.workflow.info('Stage 2b: Adding M0 scan registration workflow')

        m0scan_hmc_wf = init_m0scan_hmc_wf(
            output_dir=config.execution.aslprep_dir,
            mem_gb=mem_gb['filesize'],
            omp_nthreads=omp_nthreads,
            name='m0scan_hmc_wf',
        )
        workflow.connect([
            (inputnode, m0scan_hmc_wf, [('m0scan', 'inputnode.m0scan')]),
            (hmcref_buffer, m0scan_hmc_wf, [('aslref', 'inputnode.aslref')]),
            (m0scan_hmc_wf, m0scanreg_buffer, [
                ('outputnode.m0scan2aslref_xfm', 'm0scan2aslref_xfm'),
            ]),
            (m0scan_hmc_wf, asl_fit_reports_wf, [
                ('outputnode.m0scan_aslref', 'inputnode.m0scan_aslref'),
            ]),
        ])  # fmt:skip
    elif m0scan:
        config.loggers.workflow.info('Found M0 motion correction transforms - skipping Stage 2b')

    # Stage 3: Register fieldmap to aslref and reconstruct in ASL space
    if fieldmap_id:
        config.loggers.workflow.info('Stage 3: Adding fieldmap reconstruction workflow')
        fmap_select = pe.Node(
            KeySelect(
                fields=['fmap_ref', 'fmap_coeff', 'fmap_mask', 'sdc_method'],
                key=fieldmap_id,
            ),
            name='fmap_select',
            run_without_submitting=True,
        )

        aslref_fmap = pe.Node(ReconstructFieldmap(inverse=[True]), name='aslref_fmap', mem_gb=1)

        workflow.connect([
            (inputnode, fmap_select, [
                ('fmap_ref', 'fmap_ref'),
                ('fmap_coeff', 'fmap_coeff'),
                ('fmap_mask', 'fmap_mask'),
                ('sdc_method', 'sdc_method'),
                ('fmap_id', 'keys'),
            ]),
            (fmapref_buffer, aslref_fmap, [('out', 'target_ref_file')]),
            (fmapreg_buffer, aslref_fmap, [('aslref2fmap_xfm', 'transforms')]),
            (fmap_select, aslref_fmap, [
                ('fmap_coeff', 'in_coeffs'),
                ('fmap_ref', 'fmap_ref_file'),
            ]),
            (fmap_select, asl_fit_reports_wf, [('fmap_ref', 'inputnode.fmap_ref')]),
            (fmap_select, summary, [('sdc_method', 'distortion_correction')]),
            (fmapref_buffer, asl_fit_reports_wf, [('out', 'inputnode.sdc_aslref')]),
            (fmapreg_buffer, asl_fit_reports_wf, [
                ('aslref2fmap_xfm', 'inputnode.aslref2fmap_xfm'),
            ]),
            (aslref_fmap, asl_fit_reports_wf, [('out_file', 'inputnode.fieldmap')]),
        ])  # fmt:skip

        if not aslref2fmap_xform:
            config.loggers.workflow.info('Stage 3: Registering fieldmap to aslref')
            fmapreg_wf = init_coeff2epi_wf(
                debug='fieldmaps' in config.execution.debug,
                omp_nthreads=config.nipype.omp_nthreads,
                sloppy=config.execution.sloppy,
                name='fmapreg_wf',
            )

            itk_mat2txt = pe.Node(ConcatenateXFMs(out_fmt='itk'), name='itk_mat2txt')
            fmapreg_source_files = pe.Node(
                niu.Merge(2),
                name='fmapreg_source_files',
                run_without_submitting=True,
            )

            # fMRIPrep's init_ds_registration_wf will write out the ASL xfms to `anat` for
            # some reason, so we must override it.
            with OverrideDerivativesDataSink(output_workflows):
                ds_fmapreg_wf = output_workflows.init_ds_registration_wf(
                    source_file=asl_file,
                    bids_root=layout.root,
                    output_dir=config.execution.aslprep_dir,
                    source='aslref',
                    dest=re.sub(r'[^a-zA-Z0-9]', '', fieldmap_id),
                    desc='fmap',
                    name='ds_fmapreg_wf',
                )
            ds_fmapreg_wf.get_node('ds_xform').inputs.datatype = 'perf'

            workflow.connect([
                (fmap_select, fmapreg_wf, [
                    ('fmap_ref', 'inputnode.fmap_ref'),
                    ('fmap_mask', 'inputnode.fmap_mask'),
                ]),
                (fmapref_buffer, fmapreg_source_files, [('out', 'in1')]),
                (fmap_select, fmapreg_source_files, [('fmap_ref', 'in2')]),
                (fmapreg_wf, itk_mat2txt, [('outputnode.target2fmap_xfm', 'in_xfms')]),
                (itk_mat2txt, ds_fmapreg_wf, [('out_xfm', 'inputnode.xform')]),
                (fmapreg_source_files, ds_fmapreg_wf, [('out', 'inputnode.source_files')]),
                (ds_fmapreg_wf, fmapreg_buffer, [('outputnode.xform', 'aslref2fmap_xfm')]),
            ])  # fmt:skip
        else:
            config.loggers.workflow.info(
                'Stage 3: Found fieldmap transform - skipping registration'
            )
    else:
        config.loggers.workflow.info('No fieldmap correction - skipping Stage 3')

    # The brain mask estimated from the coregistration reference (Stage 4).
    # With --skull-strip-asl anatomical, the mask used is the anatomical one instead, checked
    # against this one once the coregistration (Stage 5) exists.
    aslstrip_buffer = pe.Node(niu.IdentityInterface(fields=['aslmask']), name='aslstrip_buffer')
    final_mask = (aslstrip_buffer, 'aslmask')
    if config.workflow.skull_strip_asl == 'anatomical' and not aslref_mask:
        anat_aslmask_wf = init_anat_aslmask_wf()
        final_mask = (anat_aslmask_wf, 'outputnode.mask_file')
        workflow.connect([
            (inputnode, anat_aslmask_wf, [('t1w_mask', 'inputnode.t1w_mask')]),
            (regref_buffer, anat_aslmask_wf, [('aslref', 'inputnode.aslref')]),
            (aslstrip_buffer, anat_aslmask_wf, [('aslmask', 'inputnode.asl_mask')]),
            (anat_aslmask_wf, outputnode, [('outputnode.qc', 'skullstrip_qc')]),
        ])  # fmt:skip

    # Stage 4: Create coregistration reference
    # Fieldmap correction only happens during fit if this stage is needed
    if not coreg_aslref:
        config.loggers.workflow.info('Stage 4: Adding coregistration aslref workflow')

        # Select initial boldref, enhance contrast, and generate mask
        # XXX: I'm not sure if this is reachable
        if sbref_files and nb.load(sbref_files[0]).ndim > 3:
            raw_sbref_wf = init_raw_aslref_wf(
                name='raw_sbref_wf',
                asl_file=sbref_files[0],
            )
            workflow.connect(raw_sbref_wf, 'outputnode.aslref', fmapref_buffer, 'sbref_files')

        if config.workflow.skull_strip_asl == 'legacy':
            enhance_aslref_wf = init_enhance_and_skullstrip_asl_wf(
                disable_n4=config.workflow.disable_n4,
                omp_nthreads=omp_nthreads,
            )
        else:
            enhance_aslref_wf = init_synthstrip_aslref_wf(
                disable_n4=config.workflow.disable_n4,
            )
        coreg_ref_source_files = pe.Node(
            niu.Merge(3),
            name='coreg_ref_source_files',
            run_without_submitting=True,
        )

        ds_coreg_aslref_wf = init_ds_aslref_wf(
            source_file=asl_file,
            bids_root=layout.root,
            output_dir=config.execution.aslprep_dir,
            desc='coreg',
            name='ds_coreg_aslref_wf',
        )
        # Provenance: it depends on the brain extraction method, and is reused only with it.
        ds_coreg_aslref_wf.get_node('ds_aslref').inputs.meta_dict = {
            'SkullStripASL': config.workflow.skull_strip_asl
        }

        workflow.connect([
            (fmapref_buffer, enhance_aslref_wf, [('out', 'inputnode.in_file')]),
            (fmapref_buffer, coreg_ref_source_files, [('out', 'in1')]),
            (coreg_ref_source_files, ds_coreg_aslref_wf, [('out', 'inputnode.source_files')]),
            (ds_coreg_aslref_wf, regref_buffer, [('outputnode.aslref', 'aslref')]),
        ])  # fmt:skip

        # A precomputed mask is used as it is (a connection would override it)
        if not aslref_mask:
            ds_aslmask_wf = output_workflows.init_ds_boldmask_wf(
                source_file=asl_file,
                output_dir=config.execution.aslprep_dir,
                desc='brain',
                name='ds_aslmask_wf',
            )
            ds_aslmask_wf.get_node('ds_boldmask').inputs.meta_dict = {
                'SkullStripASL': config.workflow.skull_strip_asl
            }
            workflow.connect([
                (ds_coreg_aslref_wf, ds_aslmask_wf, [
                    ('outputnode.aslref', 'inputnode.source_files'),
                ]),
                (final_mask[0], ds_aslmask_wf, [(final_mask[1], 'inputnode.boldmask')]),
                (ds_aslmask_wf, aslmask_buffer, [('outputnode.boldmask', 'aslmask')]),
            ])  # fmt:skip
            if final_mask[0] is not aslstrip_buffer:
                # Records which mask was used, as well as the method
                workflow.connect([
                    (anat_aslmask_wf, ds_aslmask_wf, [
                        ('outputnode.mask_metadata', 'ds_boldmask.meta_dict'),
                    ]),
                ])  # fmt:skip

        if fieldmap_id or gradwarp:
            unwarp_aslref = pe.Node(
                ResampleSeries(jacobian=jacobian, gradwarp_jacobian=gradwarp_jacobian),
                name='unwarp_aslref',
                n_procs=omp_nthreads,
                mem_gb=mem_gb['resampled'],
            )

            workflow.connect([
                (fmapref_buffer, unwarp_aslref, [('out', 'ref_file')]),
                (enhance_aslref_wf, unwarp_aslref, [
                    ('outputnode.bias_corrected_file', 'in_file'),
                ]),
                (unwarp_aslref, ds_coreg_aslref_wf, [('out_file', 'inputnode.aslref')]),
            ])  # fmt:skip

            if not aslref_mask:
                skullstrip_asl_wf = init_skullstrip_aslref_wf(name='skullstrip_asl_wf')
                # enhance_aslref_wf already describes the brain extraction method
                skullstrip_asl_wf.__desc__ = None
                workflow.connect([
                    (ds_coreg_aslref_wf, skullstrip_asl_wf, [
                        ('outputnode.aslref', 'inputnode.in_file'),
                    ]),
                    (skullstrip_asl_wf, aslstrip_buffer, [('outputnode.mask_file', 'aslmask')]),
                ])  # fmt:skip

            if gradwarp:
                workflow.connect([
                    (gradwarp_wf, unwarp_aslref, [('outputnode.gradwarp_field', 'gradwarp_field')]),
                ])  # fmt:skip

        if fieldmap_id:
            distortion_params = pe.Node(
                DistortionParameters(
                    metadata=metadata,
                    in_file=asl_file,
                    fallback=config.workflow.fallback_total_readout_time,
                ),
                name='distortion_params',
                run_without_submitting=True,
            )

            workflow.connect([
                (aslref_fmap, unwarp_aslref, [('out_file', 'fieldmap')]),
                (distortion_params, unwarp_aslref, [
                    ('readout_time', 'ro_time'),
                    ('pe_direction', 'pe_dir'),
                ]),
                (fmapreg_buffer, coreg_ref_source_files, [('aslref2fmap_xfm', 'in2')]),
                (fmap_select, coreg_ref_source_files, [('fmap_coeff', 'in3')]),
            ])  # fmt:skip

            if not aslref2fmap_xform and gradwarp:
                # Register the fieldmap to the gradient-corrected reference: the frame in which
                # the fieldmap is later evaluated, as the gradient correction is applied first.
                gradwarp_fmapreg_ref = pe.Node(
                    ResampleSeries(jacobian=False, gradwarp_jacobian=gradwarp_jacobian),
                    name='gradwarp_fmapreg_ref',
                    mem_gb=mem_gb['resampled'],
                )
                gradwarp_fmapreg_mask = pe.Node(
                    ResampleSeries(jacobian=False, gradwarp_jacobian=False, order=0),
                    name='gradwarp_fmapreg_mask',
                    mem_gb=mem_gb['resampled'],
                )
                workflow.connect([
                    (enhance_aslref_wf, gradwarp_fmapreg_ref, [
                        ('outputnode.bias_corrected_file', 'in_file'),
                        ('outputnode.bias_corrected_file', 'ref_file'),
                    ]),
                    (enhance_aslref_wf, gradwarp_fmapreg_mask, [
                        ('outputnode.mask_file', 'in_file'),
                        ('outputnode.mask_file', 'ref_file'),
                    ]),
                    (gradwarp_wf, gradwarp_fmapreg_ref, [
                        ('outputnode.gradwarp_field', 'gradwarp_field'),
                    ]),
                    (gradwarp_wf, gradwarp_fmapreg_mask, [
                        ('outputnode.gradwarp_field', 'gradwarp_field'),
                    ]),
                    (gradwarp_fmapreg_ref, fmapreg_wf, [('out_file', 'inputnode.target_ref')]),
                    (gradwarp_fmapreg_mask, fmapreg_wf, [('out_file', 'inputnode.target_mask')]),
                ])  # fmt:skip
            elif not aslref2fmap_xform:
                workflow.connect([
                    (enhance_aslref_wf, fmapreg_wf, [
                        ('outputnode.bias_corrected_file', 'inputnode.target_ref'),
                        ('outputnode.mask_file', 'inputnode.target_mask'),
                    ]),
                ])  # fmt:skip
        elif not gradwarp:
            workflow.connect([
                (enhance_aslref_wf, ds_coreg_aslref_wf, [
                    ('outputnode.bias_corrected_file', 'inputnode.aslref'),
                ]),
                (enhance_aslref_wf, aslstrip_buffer, [('outputnode.mask_file', 'aslmask')]),
            ])  # fmt:skip
    else:
        config.loggers.workflow.info('Found coregistration reference - skipping Stage 4')

        if not aslref_mask:
            skullstrip_precomp_ref_wf = init_skullstrip_aslref_wf(name='skullstrip_precomp_ref_wf')
            skullstrip_precomp_ref_wf.inputs.inputnode.in_file = coreg_aslref
            workflow.connect([
                (skullstrip_precomp_ref_wf, aslstrip_buffer, [('outputnode.mask_file', 'aslmask')]),
                (final_mask[0], aslmask_buffer, [(final_mask[1], 'aslmask')]),
            ])  # fmt:skip

    # Stage 5: Register ASL to anatomical space
    if not aslref2anat_xform:
        config.loggers.workflow.info('Stage 5: Adding coregistration workflow')
        use_bbr = (
            True
            if 'bbr' in config.workflow.force
            else False
            if 'no-bbr' in config.workflow.force
            else None
        )
        # calculate ASL registration to T1w
        asl_reg_wf = init_bold_reg_wf(
            bold2anat_dof=config.workflow.asl2anat_dof,
            bold2anat_init=asl2anat_init,
            use_bbr=use_bbr,
            freesurfer=config.workflow.run_reconall,
            omp_nthreads=omp_nthreads,
            mem_gb=mem_gb['resampled'],
            sloppy=config.execution.sloppy,
        )

        if config.execution.sloppy and config.workflow.run_reconall:
            asl_reg_wf.inputs.bbreg_wf.bbregister.args = '--nmax 5 --subsamp 4'

        # fMRIPrep's init_ds_registration_wf will write out the ASL xfms to `anat` for some reason,
        # so we must override it.
        with OverrideDerivativesDataSink(output_workflows):
            ds_aslreg_wf = output_workflows.init_ds_registration_wf(
                source_file=asl_file,
                bids_root=layout.root,
                output_dir=config.execution.aslprep_dir,
                source='aslref',
                dest='T1w',
                desc='coreg',
                name='ds_aslreg_wf',
            )

        # The image registered to the anatomical image: the coregistration reference, or R1
        regsource = (regref_buffer, 'aslref')
        regmetadata = (asl_reg_wf, 'outputnode.metadata')
        if use_r1:
            config.loggers.workflow.info(
                'Stage 5: Registering the R1 of %s to the anatomical image',
                coreg_r1map['uri'],
            )
            asl_r1_wf = init_asl_r1_wf(
                coreg_r1map=coreg_r1map,
                fieldmap_id=fieldmap_id,
                metadata=metadata,
                asl_file=asl_file,
                m0scan_gradwarp=m0scan_gradwarp,
            )
            workflow.connect([
                (fmapref_buffer, asl_r1_wf, [('out', 'inputnode.aslref')]),
                (m0scanreg_buffer, asl_r1_wf, [('m0scan2aslref_xfm', 'inputnode.m0scan2aslref_xfm')]),
                # The mask estimated from the ASL data, as the anatomical one needs this registration
                (aslmask_buffer if aslref_mask else aslstrip_buffer, asl_r1_wf, [
                    ('aslmask', 'inputnode.asl_mask'),
                ]),
            ])  # fmt:skip
            if fieldmap_id:
                workflow.connect([(aslref_fmap, asl_r1_wf, [('out_file', 'inputnode.fieldmap')])])

            if m0scan_gradwarp:
                workflow.connect([
                    (m0scan_gradwarp_wf, asl_r1_wf, [
                        ('outputnode.gradwarp_field', 'inputnode.m0scan_gradwarp_field'),
                    ]),
                ])  # fmt:skip

            reverse_bbr_contrast(asl_reg_wf)
            bbr_wf = asl_reg_wf.get_node('bbreg_wf') or asl_reg_wf.get_node('fsl_bbr_wf')
            bbr_wf.__desc__ = (bbr_wf.__desc__ or '') + R1_COREGISTRATION_DESC.format(
                uri=coreg_r1map['uri'],
            )
            regsource = (asl_r1_wf, 'outputnode.r1_brain')

            # The transform's sidecar records the image that was registered
            r1_metadata = pe.Node(
                niu.Function(function=_update_dict, output_names=['out_dict']),
                name='r1_metadata',
                run_without_submitting=True,
            )
            r1_metadata.inputs.update = {
                'RegistrationImage': 'R1map',
                'RegistrationImageSources': [coreg_r1map['uri']],
            }
            workflow.connect([(asl_reg_wf, r1_metadata, [('outputnode.metadata', 'in_dict')])])
            regmetadata = (r1_metadata, 'out_dict')

        workflow.connect([
            (inputnode, asl_reg_wf, [
                ('t1w_preproc', 'inputnode.t1w_preproc'),
                ('t1w_mask', 'inputnode.t1w_mask'),
                ('t1w_dseg', 'inputnode.t1w_dseg'),
                # Undefined if --fs-no-reconall, but this is safe
                ('subjects_dir', 'inputnode.subjects_dir'),
                ('subject_id', 'inputnode.subject_id'),
                ('fsnative2t1w_xfm', 'inputnode.fsnative2t1w_xfm'),
            ]),
            (regsource[0], asl_reg_wf, [(regsource[1], 'inputnode.ref_bold_brain')]),
            (regref_buffer, ds_aslreg_wf, [('aslref', 'inputnode.source_files')]),
            (asl_reg_wf, ds_aslreg_wf, [('outputnode.itk_bold_to_t1', 'inputnode.xform')]),
            (regmetadata[0], ds_aslreg_wf, [(regmetadata[1], 'inputnode.metadata')]),
            (ds_aslreg_wf, outputnode, [('outputnode.xform', 'aslref2anat_xfm')]),
            (asl_reg_wf, summary, [('outputnode.fallback', 'fallback')]),
        ])  # fmt:skip
        if final_mask[0] is not aslstrip_buffer:
            workflow.connect([
                (ds_aslreg_wf, anat_aslmask_wf, [('outputnode.xform', 'inputnode.aslref2anat_xfm')]),
            ])  # fmt:skip
    else:
        config.loggers.workflow.info('Found coregistration transform - skipping Stage 5')
        outputnode.inputs.aslref2anat_xfm = aslref2anat_xform
        if final_mask[0] is not aslstrip_buffer:
            anat_aslmask_wf.inputs.inputnode.aslref2anat_xfm = aslref2anat_xform

    return workflow


def init_asl_native_wf(
    *,
    asl_file: str,
    m0scan: str | None = None,
    fieldmap_id: str | None = None,
    jacobian: bool = False,
    gradwarp: bool = False,
    m0scan_gradwarp: bool = False,
    gradwarp_jacobian: bool = True,
    omp_nthreads: int = 1,
    name: str = 'asl_native_wf',
) -> pe.Workflow:
    r"""Apply minimal resampling workflow.

    This workflow resamples to aslref space with head motion and susceptibility distortion
    correction. It also selects the transforms needed to perform further resampling.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep import config
            from aslprep.workflows.asl.fit import init_asl_native_wf

            with mock_config():
                asl_file = (
                    config.execution.bids_dir / "sub-01" / "perf" /
                    "sub-01_asl.nii.gz"
                )
                wf = init_asl_native_wf(asl_file=str(asl_file))

    Parameters
    ----------
    asl_file
        List of paths to NIfTI files.
    m0scan
        Path to M0 NIfTI file, if available.
    fieldmap_id
        ID of the fieldmap to use to correct this ASL series. If :obj:`None`,
        no correction will be applied.
    jacobian
        Whether to apply the Jacobian determinant of the fieldmap during resampling.
    gradwarp
        Whether to apply a gradient nonlinearity displacement field (``gradwarp_field``).
    m0scan_gradwarp
        Whether to apply a gradient nonlinearity displacement field to the separate M0 scan
        (``m0scan_gradwarp_field``).
    gradwarp_jacobian
        Whether to modulate intensities by the Jacobian determinant of the fields.

    Inputs
    ------
    aslref
        ASL reference file
    asl_mask
        Mask of ASL reference file
    m0scan
        If M0Type is 'Separate', then the M0 file will resampled into ASL reference space.
        Otherwise, this field will be undefined.
    motion_xfm
        Affine transforms from each ASL volume to ``hmc_aslref``, written
        as concatenated ITK affine transforms.
    m0scan2aslref_xfm
        Affine transform mapping from M0 scan space to ASL reference space, if applicable.
    aslref2fmap_xfm
        Affine transform mapping from ASL reference space to the fieldmap
        space, if applicable.
    fmap_id
        Unique identifiers to select fieldmap files
    fmap_ref
        List of fieldmap reference files (collated with fmap_id)
    fmap_coeff
        List of lists of spline coefficient files (collated with fmap_id)
    gradwarp_field
        Gradient nonlinearity displacement field on the ASL reference grid, if applicable.
    m0scan_gradwarp_field
        Gradient nonlinearity displacement field for the M0 scan, if applicable.

    Outputs
    -------
    asl_minimal
        ASL series ready for further resampling.
        This ASL file will only contain the volumes needed for processing.
    asl_native
        ASL series resampled into ASL reference space.
        Head motion and susceptibility distortion correction will be applied to each file.
    m0scan_native
        If M0Type is 'Separate', then the M0 file will resampled into ASL reference space.
        Otherwise, this field will be undefined.
    metadata
        Metadata dictionary of ASL series with the shortest echo
    motion_xfm
        Motion correction transforms for further correcting asl_minimal.
        For multi-echo data, motion correction has already been applied, so
        this will be undefined.
    """
    from aslprep.interfaces.ants import RegridImage
    from aslprep.interfaces.utility import GetImageType
    from aslprep.utils.misc import estimate_asl_mem_usage

    layout = config.execution.layout
    metadata = layout.get_metadata(asl_file)

    _, mem_gb = estimate_asl_mem_usage(asl_file)

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'aslcontext',
                # ASL fit
                'aslref',
                'asl_mask',
                'm0scan',
                'motion_xfm',
                'm0scan2aslref_xfm',
                'aslref2fmap_xfm',
                # Fieldmap fit
                'fmap_ref',
                'fmap_coeff',
                'fmap_id',
                # Gradient nonlinearity
                'gradwarp_field',
                'm0scan_gradwarp_field',
            ],
        ),
        name='inputnode',
    )
    inputnode.inputs.m0scan = m0scan

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'asl_minimal',
                'asl_native',
                'm0scan_native',
                'aslcontext',
                'metadata',
                # Transforms
                'motion_xfm',
                'm0scan2aslref_xfm',
            ],
        ),
        name='outputnode',
    )

    aslbuffer = pe.Node(
        niu.IdentityInterface(fields=['asl_file', 'ro_time', 'pe_dir']),
        name='aslbuffer',
    )

    # Validate the ASL file
    validate_asl = pe.Node(ValidateImage(in_file=asl_file), name='validate_asl')

    # Drop volumes in the ASL file that won't be used
    # (e.g., precalculated CBF volumes if control-label pairs are available).
    processing_target = pe.Node(
        niu.Function(
            function=select_processing_target,
            input_names=['aslcontext'],
            output_names=['processing_target'],
        ),
        name='processing_target',
    )
    reduce_asl_file = pe.Node(
        ReduceASLFiles(metadata=metadata),
        name='reduce_asl_file',
    )

    workflow.connect([
        (inputnode, processing_target, [('aslcontext', 'aslcontext')]),
        (processing_target, reduce_asl_file, [('processing_target', 'processing_target')]),
        (inputnode, reduce_asl_file, [('aslcontext', 'aslcontext')]),
        (validate_asl, reduce_asl_file, [('out_file', 'asl_file')]),
        (reduce_asl_file, aslbuffer, [('asl_file', 'asl_file')]),
        (reduce_asl_file, outputnode, [
            ('aslcontext', 'aslcontext'),
            ('metadata', 'metadata'),
        ]),
    ])  # fmt:skip

    # Prepare fieldmap metadata
    if fieldmap_id:
        fmap_select = pe.Node(
            KeySelect(fields=['fmap_ref', 'fmap_coeff'], key=fieldmap_id),
            name='fmap_select',
            run_without_submitting=True,
        )

        distortion_params = pe.Node(
            DistortionParameters(
                metadata=metadata,
                in_file=asl_file,
                fallback=config.workflow.fallback_total_readout_time,
            ),
            name='distortion_params',
            run_without_submitting=True,
        )
        workflow.connect([
            (inputnode, fmap_select, [
                ('fmap_ref', 'fmap_ref'),
                ('fmap_coeff', 'fmap_coeff'),
                ('fmap_id', 'keys'),
            ]),
            (reduce_asl_file, distortion_params, [
                ('metadata', 'metadata'),
                ('asl_file', 'in_file'),
            ]),
            (distortion_params, aslbuffer, [
                ('readout_time', 'ro_time'),
                ('pe_direction', 'pe_dir'),
            ]),
        ])  # fmt:skip

    # Resample ASL to aslref
    aslref_asl = pe.Node(
        ResampleSeries(jacobian=jacobian, gradwarp_jacobian=gradwarp_jacobian),
        name='aslref_asl',
        n_procs=omp_nthreads,
        mem_gb=mem_gb['resampled'],
    )

    workflow.connect([
        (inputnode, aslref_asl, [
            ('aslref', 'ref_file'),
            ('motion_xfm', 'transforms'),
        ]),
        (aslbuffer, aslref_asl, [
            ('asl_file', 'in_file'),
            ('ro_time', 'ro_time'),
            ('pe_dir', 'pe_dir'),
        ]),
    ])  # fmt:skip

    if m0scan:
        # In some cases, the separate M0 scan may have a different resolution than the ASL scan.
        # We resample the M0 scan to match the ASL scan at this point.
        get_image_type = pe.Node(GetImageType(), name='get_image_type')
        workflow.connect([(inputnode, get_image_type, [('m0scan', 'image')])])

        # Resample separate M0 file to aslref resolution/orientation.
        # Gaussian interpolation smooths the image, which reduces aliasing when a finer M0 scan
        # is moved to the ASL grid. An M0 scan that is already on that grid is passed through,
        # as the smoothing would only dilute M0 at tissue boundaries (and so inflate CBF there).
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

        # Apply transforms to M0 scan to aslref space
        aslref_m0scan = pe.Node(
            ResampleSeries(
                jacobian=jacobian,
                gradwarp_jacobian=gradwarp_jacobian,
            ),
            name='aslref_m0scan',
            n_procs=omp_nthreads,
        )
        workflow.connect([
            (inputnode, aslref_m0scan, [
                ('aslref', 'ref_file'),
                ('m0scan2aslref_xfm', 'transforms'),
            ]),
            (resample_m0scan_to_asl, aslref_m0scan, [('output_image', 'in_file')]),
            (aslbuffer, aslref_m0scan, [
                ('ro_time', 'ro_time'),
                ('pe_dir', 'pe_dir'),
            ]),
            (aslref_m0scan, outputnode, [('out_file', 'm0scan_native')]),
        ])  # fmt:skip

        if m0scan_gradwarp:
            # The M0 scan was moved onto the ASL reference grid above without changing its
            # world coordinates, and the field is defined in scanner coordinates.
            workflow.connect([
                (inputnode, aslref_m0scan, [('m0scan_gradwarp_field', 'gradwarp_field')]),
            ])  # fmt:skip

    if fieldmap_id:
        aslref_fmap = pe.Node(ReconstructFieldmap(inverse=[True]), name='aslref_fmap', mem_gb=1)
        workflow.connect([
            (inputnode, aslref_fmap, [
                ('aslref', 'target_ref_file'),
                ('aslref2fmap_xfm', 'transforms'),
            ]),
            (fmap_select, aslref_fmap, [
                ('fmap_coeff', 'in_coeffs'),
                ('fmap_ref', 'fmap_ref_file'),
            ]),
            (aslref_fmap, aslref_asl, [('out_file', 'fieldmap')]),
        ])  # fmt:skip

        if m0scan:
            # This assumes that the M0 scan and the ASL file have the same distortion.
            # In some edge cases, that might not be true.
            workflow.connect([(aslref_fmap, aslref_m0scan, [('out_file', 'fieldmap')])])

    if gradwarp:
        workflow.connect([(inputnode, aslref_asl, [('gradwarp_field', 'gradwarp_field')])])

    workflow.connect([
        (inputnode, outputnode, [('motion_xfm', 'motion_xfm')]),
        (aslbuffer, outputnode, [('asl_file', 'asl_minimal')]),
        (aslref_asl, outputnode, [('out_file', 'asl_native')]),
    ])  # fmt:skip

    return workflow


#: Boilerplate for the registration of R1 to the anatomical image, formatted with the map's URI.
R1_COREGISTRATION_DESC = (
    ' The image registered to the anatomical image was not the ASL reference, but the R1 '
    '(1 / T1) of a quantitative map on the grid of the M0 scan ({uri}), moved into the frame '
    "of the ASL reference with the M0 scan's transform and masked with the ASL brain mask, "
    'as its gray-white contrast is much stronger than that of the M0 scan. As white matter is '
    'brighter than gray matter in R1, the contrast of the boundary-based registration was '
    'reversed.'
)


def use_coregistration_r1(coreg_r1map, *, m0scan, reference_volume_type, aslref2anat_xform):
    """Decide whether R1, from a map on the grid of the M0 scan, is registered to the anatomy.

    The M0 scan's transform to the ASL reference moves the map into the coregistration
    reference's frame. It is used only when the reference was made from the separate M0 scan,
    so that the transform is between images of the same contrast, and only when the M0 scan
    has one volume, so that there is one transform.

    Parameters
    ----------
    coreg_r1map : :obj:`dict` or None
        The map, from :func:`~aslprep.utils.bids.collect_m0scan_t1map`.
    m0scan : :obj:`str` or None
        The separate M0 scan.
    reference_volume_type : :obj:`str`
        The volume type of the ASL reference (see :func:`select_reference_volume_type`).
    aslref2anat_xform : :obj:`str` or None
        A precomputed ASL-to-anatomical transform, in which case there is no registration.

    Returns
    -------
    :obj:`bool`
    """
    if coreg_r1map is None or not m0scan or aslref2anat_xform:
        return False

    reason = None
    if reference_volume_type != 'separate_m0scan':
        reason = f'the ASL reference was made from the {reference_volume_type} volumes'
    elif nb.load(m0scan).ndim == 4 and nb.load(m0scan).shape[3] > 1:
        reason = 'the M0 scan has more than one volume'

    if reason is not None:
        config.loggers.workflow.info(
            'Not registering the R1 of %s to the anatomical image, as %s; '
            'the ASL reference is registered instead.',
            coreg_r1map['uri'],
            reason,
        )
        return False

    return True


def reverse_bbr_contrast(asl_reg_wf):
    """Set a registration workflow's boundary-based registration for white matter brighter.

    fMRIPrep's boundary-based registrations expect gray matter brighter than white matter, as in
    a T2*-weighted image. FLIRT's default slope (-0.5) is reversed, and FreeSurfer's
    ``bbregister`` is given the T1 contrast. Without boundary-based registration, nothing changes.

    Parameters
    ----------
    asl_reg_wf : :obj:`~nipype.pipeline.engine.Workflow`
        The workflow from :func:`~fmriprep.workflows.bold.registration.init_bold_reg_wf`.
    """
    fsl_bbr_wf = asl_reg_wf.get_node('fsl_bbr_wf')
    flt_bbr = fsl_bbr_wf.get_node('flt_bbr') if fsl_bbr_wf is not None else None
    if flt_bbr is not None:
        flt_bbr.inputs.bbrslope = 0.5

    bbreg_wf = asl_reg_wf.get_node('bbreg_wf')
    bbregister = bbreg_wf.get_node('bbregister') if bbreg_wf is not None else None
    if bbregister is not None:
        bbregister.inputs.contrast_type = 't1'


def init_asl_r1_wf(
    *,
    coreg_r1map: dict,
    fieldmap_id: str | None,
    metadata: dict,
    asl_file: str,
    m0scan_gradwarp: bool = False,
    name: str = 'asl_r1_wf',
) -> pe.Workflow:
    """Bring R1, from a map on the grid of the M0 scan, into the coregistration reference's frame.

    The map is converted to R1 on its own grid, put on the ASL reference grid, moved with the
    M0 scan's transform to the ASL reference, corrected for distortion as the M0 scan is, and
    masked with the ASL brain mask.
    (It has no workflow graph here, as it needs an existing map to be built.)

    Parameters
    ----------
    coreg_r1map
        The map, from :func:`~aslprep.utils.bids.collect_m0scan_t1map`.
    fieldmap_id
        ID of the fieldmap that corrects the ASL series, if any.
    metadata
        The ASL file's metadata, for its distortion parameters.
    asl_file
        The ASL file.
    m0scan_gradwarp
        Whether the M0 scan has a gradient nonlinearity field (``m0scan_gradwarp_field``).

    Inputs
    ------
    aslref
        The image on whose grid the coregistration reference lies.
    m0scan2aslref_xfm
        The M0 scan's transform to the ASL reference.
    asl_mask
        The ASL brain mask.
    fieldmap
        The fieldmap, in the frame of the ASL reference, if any.
    m0scan_gradwarp_field
        The M0 scan's gradient nonlinearity field, if any.

    Outputs
    -------
    r1_brain
        R1, in the frame of the coregistration reference, masked with the ASL brain mask.
    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.nibabel import ApplyMask

    from aslprep.interfaces.ants import RegridImage
    from aslprep.interfaces.t1map import CoregistrationR1

    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'aslref',
                'm0scan2aslref_xfm',
                'asl_mask',
                'fieldmap',
                'm0scan_gradwarp_field',
            ],
        ),
        name='inputnode',
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=['r1_brain']), name='outputnode')

    coreg_r1 = pe.Node(
        CoregistrationR1(
            in_file=coreg_r1map['t1map'],
            metadata=coreg_r1map['metadata'],
            suffix=coreg_r1map['suffix'],
        ),
        name='coreg_r1',
        mem_gb=1,
    )
    # As the M0 scan is, before its transform is applied (an M0 scan on the grid is passed through)
    regrid_r1 = pe.Node(
        RegridImage(interpolation='Linear', transforms=['identity'], input_image_type=0),
        name='regrid_r1',
    )
    # R1 is a quantitative map, so it is not modulated by any Jacobian
    aslref_r1 = pe.Node(
        ResampleSeries(jacobian=False, gradwarp_jacobian=False),
        name='aslref_r1',
    )
    mask_r1 = pe.Node(ApplyMask(), name='mask_r1')

    workflow.connect([
        (inputnode, regrid_r1, [('aslref', 'reference_image')]),
        (coreg_r1, regrid_r1, [('out_file', 'input_image')]),
        (inputnode, aslref_r1, [
            ('aslref', 'ref_file'),
            ('m0scan2aslref_xfm', 'transforms'),
        ]),
        (regrid_r1, aslref_r1, [('output_image', 'in_file')]),
        (aslref_r1, mask_r1, [('out_file', 'in_file')]),
        (inputnode, mask_r1, [('asl_mask', 'in_mask')]),
        (mask_r1, outputnode, [('out_file', 'r1_brain')]),
    ])  # fmt:skip

    if fieldmap_id:
        # This assumes that the M0 scan and the ASL file have the same distortion, as the
        # resampling of the M0 scan does.
        distortion_params = pe.Node(
            DistortionParameters(
                metadata=metadata,
                in_file=asl_file,
                fallback=config.workflow.fallback_total_readout_time,
            ),
            name='distortion_params',
            run_without_submitting=True,
        )
        workflow.connect([
            (inputnode, aslref_r1, [('fieldmap', 'fieldmap')]),
            (distortion_params, aslref_r1, [
                ('readout_time', 'ro_time'),
                ('pe_direction', 'pe_dir'),
            ]),
        ])  # fmt:skip

    if m0scan_gradwarp:
        workflow.connect([(inputnode, aslref_r1, [('m0scan_gradwarp_field', 'gradwarp_field')])])

    return workflow


def _update_dict(in_dict, update):
    """Return a copy of a dictionary, updated with another."""
    return {**(in_dict or {}), **update}


def _select_ref(sbref_files, aslref_files):
    """Select first sbref or aslref file, preferring sbref if available."""
    from niworkflows.utils.connections import listify

    refs = sbref_files or aslref_files
    return listify(refs)[0]


def select_reference_volume_type(aslcontext: str, metadata: dict, prioritize_m0: bool) -> str:
    """Select the reference volume type for the ASL series.

    This function selects the reference volume type for the ASL series based on the ASL context
    file and the metadata.
    """
    import pandas as pd

    aslcontext_df = pd.read_table(aslcontext)
    separate_m0scan = metadata['M0Type'] == 'Separate'
    if 'm0scan' in aslcontext_df['volume_type'].tolist() and prioritize_m0:
        target_type = 'm0scan'
    elif separate_m0scan and prioritize_m0:
        target_type = 'separate_m0scan'
    elif 'cbf' in aslcontext_df['volume_type'].tolist():
        target_type = 'cbf'
    elif 'deltam' in aslcontext_df['volume_type'].tolist():
        target_type = 'deltam'
    elif 'm0scan' in aslcontext_df['volume_type'].tolist():
        target_type = 'm0scan'
    elif separate_m0scan:
        target_type = 'separate_m0scan'
    else:
        target_type = 'control'

    return target_type


def skips_motion_correction(aslcontext: str, use_ge: bool) -> bool:
    """Whether head-motion correction is skipped for the ASL series.

    The GE-specific path skips it for series without label-control pairs (an M0 with delta-M or
    CBF images), whose few volumes share too little contrast for registration to beat none
    (PennLINC/aslprep#769). Short label-control series are still motion corrected.
    """
    import pandas as pd

    volume_types = set(pd.read_table(aslcontext)['volume_type'])
    return use_ge and not volume_types & {'control', 'label'}
