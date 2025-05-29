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
import typing as ty

import bids
import nibabel as nb
import numpy as np
from fmriprep.interfaces.resampling import (
    DistortionParameters,
    ReconstructFieldmap,
    ResampleSeries,
)
from fmriprep.utils.bids import extract_entities
from fmriprep.workflows.bold import outputs as output_workflows
from fmriprep.workflows.bold.reference import init_validation_and_dummies_wf
from fmriprep.workflows.bold.registration import init_bold_reg_wf
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.func.util import init_enhance_and_skullstrip_bold_wf, init_skullstrip_bold_wf
from niworkflows.interfaces.header import ValidateImage
from niworkflows.interfaces.nitransforms import ConcatenateXFMs
from niworkflows.interfaces.utility import KeySelect
from sdcflows.workflows.apply.registration import init_coeff2epi_wf

# ASL workflows
from aslprep import config
from aslprep.interfaces.bids import OverrideDerivativesDataSink
from aslprep.interfaces.reports import FunctionalSummary
from aslprep.interfaces.utility import ReduceASLFiles
from aslprep.utils.asl import select_processing_target
from aslprep.workflows.asl.hmc import init_asl_hmc_wf
from aslprep.workflows.asl.outputs import init_asl_fit_reports_wf, init_ds_aslref_wf
from aslprep.workflows.asl.reference import init_raw_aslref_wf


def get_sbrefs(
    asl_file: str,
    entity_overrides: dict[str, ty.Any],
    layout: bids.BIDSLayout,
) -> list[str]:
    """Find single-band reference(s) associated with ASL file.

    Parameters
    ----------
    asl_file
        List of absolute paths to ASL files
    entity_overrides
        Query parameters to override defaults
    layout
        :class:`~bids.layout.BIDSLayout` to query

    Returns
    -------
    sbref_files
        List of absolute paths to sbref files associated with input ASL files,
        sorted by EchoTime
    """
    entities = extract_entities(asl_file)
    entities.update(suffix='sbref', extension=['.nii', '.nii.gz'], **entity_overrides)

    return layout.get(return_type='file', **entities)


def init_asl_fit_wf(
    *,
    asl_file: str,
    aslcontext: str,
    m0scan: str | None,
    use_ge: bool,
    precomputed: dict = None,
    fieldmap_id: str | None = None,
    jacobian: bool = False,
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
                    aslcontext=str(asl_file.replace('.nii.gz', 'context.tsv')),
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
        reference image, as GE deltam volumes exhibit extreme background noise.
    precomputed
        Dictionary containing precomputed derivatives to reuse, if possible.
    fieldmap_id
        ID of the fieldmap to use to correct this ASL series. If :obj:`None`,
        no correction will be applied.

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
    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    from aslprep.utils.misc import estimate_asl_mem_usage

    if precomputed is None:
        precomputed = {}
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
        niu.IdentityInterface(fields=['aslref', 'aslmask']),
        name='regref_buffer',
    )

    if hmc_aslref:
        hmcref_buffer.inputs.aslref = hmc_aslref
        config.loggers.workflow.debug('Reusing motion correction reference: %s', hmc_aslref)
    if hmc_xforms:
        hmc_buffer.inputs.hmc_xforms = hmc_xforms
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
            registration_init=config.workflow.asl2anat_init,
            pe_direction=metadata.get('PhaseEncodingDirection'),
            tr=metadata['RepetitionTime'],
            orientation=orientation,
        ),
        name='summary',
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    asl_fit_reports_wf = init_asl_fit_reports_wf(
        # TODO: Enable sdc report even if we find coregref
        sdc_correction=fieldmap_id is not None,
        separate_m0scan=m0scan,
        freesurfer=config.workflow.run_reconall,
        output_dir=config.execution.aslprep_dir,
    )

    workflow.connect([
        (hmcref_buffer, fmapref_buffer, [('aslref', 'aslref_files')]),
        (hmcref_buffer, outputnode, [('aslref', 'hmc_aslref')]),
        (regref_buffer, outputnode, [
            ('aslref', 'coreg_aslref'),
            ('aslmask', 'asl_mask'),
        ]),
        (fmapreg_buffer, outputnode, [('aslref2fmap_xfm', 'aslref2fmap_xfm')]),
        (hmc_buffer, outputnode, [('hmc_xforms', 'motion_xfm')]),
        (m0scanreg_buffer, outputnode, [('m0scan2aslref_xfm', 'm0scan2aslref_xfm')]),
        (inputnode, asl_fit_reports_wf, [
            ('asl_file', 'inputnode.source_file'),
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
            bids_root=layout.root,
            output_dir=config.execution.aslprep_dir,
            desc='hmc',
            name='ds_hmc_aslref_wf',
        )
        ds_hmc_aslref_wf.inputs.inputnode.source_files = [asl_file]

        workflow.connect([
            (hmc_aslref_wf, hmcref_buffer, [
                ('outputnode.asl_file', 'asl_file'),
                ('outputnode.aslref', 'aslref'),
            ]),
            (hmcref_buffer, ds_hmc_aslref_wf, [('aslref', 'inputnode.aslref')]),
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
        config.loggers.workflow.info('Stage 2: Adding motion correction workflow')
        asl_hmc_wf = init_asl_hmc_wf(
            name='asl_hmc_wf',
            mem_gb=mem_gb['filesize'],
            omp_nthreads=omp_nthreads,
        )

        with OverrideDerivativesDataSink(output_workflows):
            ds_hmc_wf = output_workflows.init_ds_hmc_wf(
                bids_root=layout.root,
                output_dir=config.execution.aslprep_dir,
            )

        ds_hmc_wf.get_node('inputnode').inputs.source_files = [asl_file]
        # fMRIPrep will write out an orig-to-boldref transform to anat, so we need to overwrite
        # some fields.
        ds_hmc_wf.get_node('ds_xforms').inputs.datatype = 'perf'
        ds_hmc_wf.get_node('ds_xforms').inputs.to = 'aslref'

        workflow.connect([
            (hmcref_buffer, asl_hmc_wf, [('aslref', 'inputnode.raw_ref_image')]),
            (reduce_asl_file, asl_hmc_wf, [
                ('asl_file', 'inputnode.asl_file'),
                ('aslcontext', 'inputnode.aslcontext'),
            ]),
            (asl_hmc_wf, ds_hmc_wf, [('outputnode.xforms', 'inputnode.xforms')]),
            (asl_hmc_wf, hmc_buffer, [('outputnode.xforms', 'hmc_xforms')]),
        ])  # fmt:skip
    else:
        config.loggers.workflow.info('Found motion correction transforms - skipping Stage 2')

    # Stage 2b: Register M0 scan to aslref
    if m0scan:
        from nipype.interfaces import fsl
        from niworkflows.interfaces.itk import MCFLIRT2ITK

        from aslprep.interfaces.bids import DerivativesDataSink
        from aslprep.interfaces.utility import MeanImage

        config.loggers.workflow.info('Stage 2b: Adding M0 scan registration workflow')

        # Register the M0 scan to the ASL reference.
        # By using MCFLIRT, we can support 4D M0 scans.
        # Register the M0 scan to the ASL reference.
        mcflirt = pe.Node(
            fsl.MCFLIRT(save_mats=True, cost='mutualinfo'),
            name='mcflirt',
            mem_gb=mem_gb['filesize'],
        )
        workflow.connect([
            (inputnode, mcflirt, [('m0scan', 'in_file')]),
            (hmcref_buffer, mcflirt, [('aslref', 'ref_file')]),
        ])  # fmt:skip

        # Calculate mean image of M0 scan
        mean_m0scan = pe.Node(
            MeanImage(),
            name='mean_m0scan',
            mem_gb=mem_gb['filesize'],
        )
        workflow.connect([
            (mcflirt, mean_m0scan, [('out_file', 'in_file')]),
            (mean_m0scan, asl_fit_reports_wf, [('out_file', 'inputnode.m0scan_aslref')]),
        ])  # fmt:skip

        fsl2itk = pe.Node(MCFLIRT2ITK(), name='fsl2itk', mem_gb=0.05, n_procs=omp_nthreads)
        workflow.connect([
            (mcflirt, fsl2itk, [('mat_file', 'in_files')]),
            (fsl2itk, m0scanreg_buffer, [('out_file', 'm0scan2aslref_xfm')]),
        ])  # fmt:skip

        ds_m0scan2aslref_xfm = pe.Node(
            DerivativesDataSink(
                base_directory=config.execution.aslprep_dir,
                source_file=m0scan,
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
            (m0scanreg_buffer, ds_m0scan2aslref_xfm, [('m0scan2aslref_xfm', 'in_file')]),
        ])  # fmt:skip

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

            # fMRIPrep's init_ds_registration_wf will write out the ASL xfms to `anat` for
            # some reason, so we must override it.
            with OverrideDerivativesDataSink(output_workflows):
                ds_fmapreg_wf = output_workflows.init_ds_registration_wf(
                    bids_root=layout.root,
                    output_dir=config.execution.aslprep_dir,
                    source='aslref',
                    dest=fieldmap_id.replace('_', ''),
                    name='ds_fmapreg_wf',
                )
            ds_fmapreg_wf.get_node('inputnode').inputs.source_files = [asl_file]
            ds_fmapreg_wf.get_node('ds_xform').inputs.datatype = 'perf'

            workflow.connect([
                (fmap_select, fmapreg_wf, [
                    ('fmap_ref', 'inputnode.fmap_ref'),
                    ('fmap_mask', 'inputnode.fmap_mask'),
                ]),
                (fmapreg_wf, itk_mat2txt, [('outputnode.target2fmap_xfm', 'in_xfms')]),
                (itk_mat2txt, ds_fmapreg_wf, [('out_xfm', 'inputnode.xform')]),
                (ds_fmapreg_wf, fmapreg_buffer, [('outputnode.xform', 'aslref2fmap_xfm')]),
            ])  # fmt:skip
        else:
            config.loggers.workflow.info(
                'Stage 3: Found fieldmap transform - skipping registration'
            )
    else:
        config.loggers.workflow.info('No fieldmap correction - skipping Stage 3')

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

        enhance_aslref_wf = init_enhance_and_skullstrip_bold_wf(omp_nthreads=omp_nthreads)

        ds_coreg_aslref_wf = init_ds_aslref_wf(
            bids_root=layout.root,
            output_dir=config.execution.aslprep_dir,
            desc='coreg',
            name='ds_coreg_aslref_wf',
        )
        ds_aslmask_wf = output_workflows.init_ds_boldmask_wf(
            output_dir=config.execution.aslprep_dir,
            desc='brain',
            name='ds_aslmask_wf',
        )
        ds_aslmask_wf.inputs.inputnode.source_files = [asl_file]

        workflow.connect([
            (fmapref_buffer, enhance_aslref_wf, [('out', 'inputnode.in_file')]),
            (hmc_aslref_source_buffer, ds_coreg_aslref_wf, [
                ('in_file', 'inputnode.source_files'),
            ]),
            (ds_coreg_aslref_wf, regref_buffer, [('outputnode.aslref', 'aslref')]),
            (ds_aslmask_wf, regref_buffer, [('outputnode.boldmask', 'aslmask')]),
        ])  # fmt:skip

        if fieldmap_id:
            distortion_params = pe.Node(
                DistortionParameters(
                    metadata=metadata,
                    in_file=asl_file,
                ),
                name='distortion_params',
                run_without_submitting=True,
            )

            unwarp_aslref = pe.Node(
                ResampleSeries(jacobian=jacobian),
                name='unwarp_aslref',
                n_procs=omp_nthreads,
                mem_gb=mem_gb['resampled'],
            )

            skullstrip_asl_wf = init_skullstrip_bold_wf()

            workflow.connect([
                (fmapref_buffer, unwarp_aslref, [('out', 'ref_file')]),
                (enhance_aslref_wf, unwarp_aslref, [
                    ('outputnode.bias_corrected_file', 'in_file'),
                ]),
                (aslref_fmap, unwarp_aslref, [('out_file', 'fieldmap')]),
                (distortion_params, unwarp_aslref, [
                    ('readout_time', 'ro_time'),
                    ('pe_direction', 'pe_dir'),
                ]),
                (unwarp_aslref, ds_coreg_aslref_wf, [('out_file', 'inputnode.aslref')]),
                (ds_coreg_aslref_wf, skullstrip_asl_wf, [
                    ('outputnode.aslref', 'inputnode.in_file'),
                ]),
                (skullstrip_asl_wf, ds_aslmask_wf, [
                    ('outputnode.mask_file', 'inputnode.boldmask'),
                ]),
            ])  # fmt:skip

            if not aslref2fmap_xform:
                workflow.connect([
                    (enhance_aslref_wf, fmapreg_wf, [
                        ('outputnode.bias_corrected_file', 'inputnode.target_ref'),
                        ('outputnode.mask_file', 'inputnode.target_mask'),
                    ]),
                ])  # fmt:skip
        else:
            workflow.connect([
                (enhance_aslref_wf, ds_coreg_aslref_wf, [
                    ('outputnode.bias_corrected_file', 'inputnode.aslref'),
                ]),
                (enhance_aslref_wf, ds_aslmask_wf, [
                    ('outputnode.mask_file', 'inputnode.boldmask'),
                ]),
            ])  # fmt:skip
    else:
        config.loggers.workflow.info('Found coregistration reference - skipping Stage 4')

        # TODO: Allow precomputed bold masks to be passed
        # Also needs consideration for how it interacts above
        skullstrip_precomp_ref_wf = init_skullstrip_bold_wf(name='skullstrip_precomp_ref_wf')
        skullstrip_precomp_ref_wf.inputs.inputnode.in_file = coreg_aslref
        workflow.connect([
            (skullstrip_precomp_ref_wf, regref_buffer, [('outputnode.mask_file', 'aslmask')])
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
            bold2anat_init=config.workflow.asl2anat_init,
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
                bids_root=layout.root,
                output_dir=config.execution.aslprep_dir,
                source='aslref',
                dest='T1w',
                name='ds_aslreg_wf',
            )
        ds_aslreg_wf.get_node('inputnode').inputs.source_files = [asl_file]

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
            (regref_buffer, asl_reg_wf, [('aslref', 'inputnode.ref_bold_brain')]),
            # Incomplete sources
            (regref_buffer, ds_aslreg_wf, [('aslref', 'inputnode.source_files')]),
            (asl_reg_wf, ds_aslreg_wf, [('outputnode.itk_bold_to_t1', 'inputnode.xform')]),
            (ds_aslreg_wf, outputnode, [('outputnode.xform', 'aslref2anat_xfm')]),
            (asl_reg_wf, summary, [('outputnode.fallback', 'fallback')]),
        ])  # fmt:skip
    else:
        config.loggers.workflow.info('Found coregistration transform - skipping Stage 5')
        outputnode.inputs.aslref2anat_xfm = aslref2anat_xform

    return workflow


def init_asl_native_wf(
    *,
    asl_file: str,
    m0scan: str | None = None,
    fieldmap_id: str | None = None,
    jacobian: bool = False,
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
            ],
        ),
        name='inputnode',
    )

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
            DistortionParameters(metadata=metadata, in_file=asl_file),
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
        ResampleSeries(jacobian=jacobian),
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
        # Resample separate M0 file to aslref
        aslref_m0scan = pe.Node(
            ResampleSeries(
                jacobian=jacobian,
                in_file=m0scan,
            ),
            name='aslref_m0scan',
            n_procs=omp_nthreads,
        )

        workflow.connect([
            (inputnode, aslref_m0scan, [
                ('aslref', 'ref_file'),
                ('m0scan2aslref_xfm', 'transforms'),
            ]),
            (aslbuffer, aslref_m0scan, [
                ('ro_time', 'ro_time'),
                ('pe_dir', 'pe_dir'),
            ]),
            (aslref_m0scan, outputnode, [('out_file', 'm0scan_native')]),
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
            workflow.connect([(aslref_fmap, aslref_m0scan, [('out_file', 'fieldmap')])])

    workflow.connect([
        (inputnode, outputnode, [('motion_xfm', 'motion_xfm')]),
        (aslbuffer, outputnode, [('asl_file', 'asl_minimal')]),
        (aslref_asl, outputnode, [('out_file', 'asl_native')]),
    ])  # fmt:skip

    return workflow


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
