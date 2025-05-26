# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for writing out derivative files."""

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.utils.images import dseg_label
from smriprep.workflows.outputs import _bids_relative

from aslprep import config
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.bids import DerivativesDataSink

BASE_INPUT_FIELDS = {
    'asl': {
        'desc': 'preproc',
        'suffix': 'asl',
    },
    'aslref': {
        'suffix': 'aslref',
    },
    'asl_mask': {
        'desc': 'brain',
        'suffix': 'mask',
    },
    # CBF outputs
    'cbf_ts': {
        'desc': 'timeseries',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
    },
    'mean_cbf': {
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
    },
    'att': {
        'suffix': 'att',
        'Units': 's',
    },
    'abat': {
        'suffix': 'abat',
        'Units': 's',
    },
    'abv': {
        'suffix': 'abv',
        'Units': 'fraction',
    },
    # SCORE/SCRUB outputs
    'cbf_ts_score': {
        'desc': 'scoreTimeseries',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
    },
    'mean_cbf_score': {
        'desc': 'score',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
    },
    'mean_cbf_scrub': {
        'desc': 'scrub',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
    },
    # BASIL outputs
    'mean_cbf_basil': {
        'desc': 'basil',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
    },
    'mean_cbf_gm_basil': {
        'desc': 'basilGM',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
    },
    'mean_cbf_wm_basil': {
        'desc': 'basilWM',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
    },
    'att_basil': {
        'desc': 'basil',
        'suffix': 'att',
        'Units': 's',
    },
}


def prepare_timing_parameters(metadata: dict):
    """Convert initial timing metadata to post-realignment timing metadata.

    In particular, SliceTiming metadata is invalid once STC or any realignment is applied,
    as a matrix of voxels no longer corresponds to an acquisition slice.
    Therefore, if SliceTiming is present in the metadata dictionary, and a sparse
    acquisition paradigm is detected, DelayTime or AcquisitionDuration must be derived to
    preserve the timing interpretation.

    Examples
    --------
    .. testsetup::

        >>> from unittest import mock

    If SliceTiming metadata is absent, then the only change is to note that
    STC has not been applied:

    >>> prepare_timing_parameters(dict(RepetitionTime=2))
    {"RepetitionTime": 2, "SliceTimingCorrected": False}
    >>> prepare_timing_parameters(dict(RepetitionTime=2, DelayTime=0.5))
    {"RepetitionTime": 2, "DelayTime": 0.5, "SliceTimingCorrected": False}
    >>> prepare_timing_parameters(dict(VolumeTiming=[0.0, 1.0, 2.0, 5.0, 6.0, 7.0],
    ...                                AcquisitionDuration=1.0))  #doctest: +NORMALIZE_WHITESPACE
    {"VolumeTiming": [0.0, 1.0, 2.0, 5.0, 6.0, 7.0], "AcquisitionDuration": 1.0,
     "SliceTimingCorrected": False}

    When SliceTiming is available and used, then ``SliceTimingCorrected`` is ``True``
    and the ``StartTime`` indicates a series offset.

    >>> with mock.patch("fmriprep.config.workflow.ignore", []):
    ...     prepare_timing_parameters(dict(RepetitionTime=2, SliceTiming=[0.0, 0.2, 0.4, 0.6]))
    {"RepetitionTime": 2, "SliceTimingCorrected": True, "DelayTime": 1.2, "StartTime": 0.3}
    >>> with mock.patch("fmriprep.config.workflow.ignore", []):
    ...     prepare_timing_parameters(
    ...         dict(VolumeTiming=[0.0, 1.0, 2.0, 5.0, 6.0, 7.0],
    ...              SliceTiming=[0.0, 0.2, 0.4, 0.6, 0.8]))  #doctest: +NORMALIZE_WHITESPACE
    {"VolumeTiming": [0.0, 1.0, 2.0, 5.0, 6.0, 7.0], "SliceTimingCorrected": True,
     "AcquisitionDuration": 1.0, "StartTime": 0.4}

    When SliceTiming is available and not used, then ``SliceTimingCorrected`` is ``False``
    and TA is indicated with ``DelayTime`` or ``AcquisitionDuration``.

    >>> with mock.patch("fmriprep.config.workflow.ignore", ["slicetiming"]):
    ...     prepare_timing_parameters(dict(RepetitionTime=2, SliceTiming=[0.0, 0.2, 0.4, 0.6]))
    {"RepetitionTime": 2, "SliceTimingCorrected": False, "DelayTime": 1.2}
    >>> with mock.patch("fmriprep.config.workflow.ignore", ["slicetiming"]):
    ...     prepare_timing_parameters(
    ...         dict(VolumeTiming=[0.0, 1.0, 2.0, 5.0, 6.0, 7.0],
    ...              SliceTiming=[0.0, 0.2, 0.4, 0.6, 0.8]))  #doctest: +NORMALIZE_WHITESPACE
    {"VolumeTiming": [0.0, 1.0, 2.0, 5.0, 6.0, 7.0], "SliceTimingCorrected": False,
     "AcquisitionDuration": 1.0}

    If SliceTiming metadata is present but empty, then treat it as missing:

    >>> with mock.patch("fmriprep.config.workflow.ignore", []):
    ...     prepare_timing_parameters(dict(RepetitionTime=2, SliceTiming=[]))
    {"RepetitionTime": 2, "SliceTimingCorrected": False}
    >>> with mock.patch("fmriprep.config.workflow.ignore", []):
    ...     prepare_timing_parameters(dict(RepetitionTime=2, SliceTiming=[0.0]))
    {"RepetitionTime": 2, "SliceTimingCorrected": False}
    """
    timing_parameters = {
        key: metadata[key]
        for key in (
            'RepetitionTimePreparation',
            'VolumeTiming',
            'DelayTime',
            'AcquisitionDuration',
            'SliceTiming',
        )
        if key in metadata
    }

    # Treat SliceTiming of [] or length 1 as equivalent to missing and remove it in any case
    slice_timing = timing_parameters.pop('SliceTiming', [])

    run_stc = len(slice_timing) > 1 and 'slicetiming' not in config.workflow.ignore
    timing_parameters['SliceTimingCorrected'] = run_stc

    return timing_parameters


def init_asl_fit_reports_wf(
    *,
    sdc_correction: bool,
    freesurfer: bool,  # noqa:U100
    output_dir: str,
    name='asl_fit_reports_wf',
) -> pe.Workflow:
    """Set up a battery of datasinks to store reports in the right location.

    Copied from fMRIPrep's init_bold_fit_reports_wf.
    Modifications include changes to fields and variable names, which aren't important,
    and DerivativesDataSink suffixes, which are.

    Parameters
    ----------
    freesurfer : :obj:`bool`
        FreeSurfer was enabled
    output_dir : :obj:`str`
        Directory in which to save derivatives
    name : :obj:`str`
        Workflow name (default: anat_reports_wf)

    Inputs
    ------
    source_file
        Input BOLD images
    std_t1w
        T1w image resampled to standard space
    std_mask
        Mask of skull-stripped template
    subject_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    t1w_conform_report
        Conformation report
    t1w_preproc
        The T1w reference map, which is calculated as the average of bias-corrected
        and preprocessed T1w images, defining the anatomical space.
    t1w_dseg
        Segmentation in T1w space
    t1w_mask
        Brain (binary) mask estimated by brain extraction.
    template
        Template space and specifications

    """
    from nireports.interfaces.reporting.base import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )
    from sdcflows.interfaces.reportlets import FieldmapReportlet

    workflow = pe.Workflow(name=name)

    inputfields = [
        'source_file',
        'sdc_aslref',
        'coreg_aslref',
        'aslref2anat_xfm',
        'aslref2fmap_xfm',
        't1w_preproc',
        't1w_mask',
        't1w_dseg',
        'asl_mask',
        'fieldmap',
        'fmap_ref',
        # May be missing
        'subject_id',
        'subjects_dir',
        # Report snippets
        'summary_report',
        'validation_report',
    ]
    inputnode = pe.Node(niu.IdentityInterface(fields=inputfields), name='inputnode')

    ds_summary = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='summary',
            datatype='figures',
            dismiss_entities=('echo',),
        ),
        name='ds_report_summary',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    ds_validation = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='validation',
            datatype='figures',
            dismiss_entities=('echo',),
        ),
        name='ds_report_validation',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # Resample anatomical references into BOLD space for plotting
    t1w_aslref = pe.Node(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            float=True,
            invert_transform_flags=[True],
            interpolation='LanczosWindowedSinc',
            args='-v',
        ),
        name='t1w_aslref',
        mem_gb=1,
    )

    t1w_wm = pe.Node(
        niu.Function(function=dseg_label),
        name='t1w_wm',
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    t1w_wm.inputs.label = 2  # BIDS default is WM=2

    aslref_wm = pe.Node(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            invert_transform_flags=[True],
            interpolation='NearestNeighbor',
            args='-v',
        ),
        name='aslref_wm',
        mem_gb=1,
    )

    workflow.connect([
        (inputnode, ds_summary, [
            ('source_file', 'source_file'),
            ('summary_report', 'in_file'),
        ]),
        (inputnode, ds_validation, [
            ('source_file', 'source_file'),
            ('validation_report', 'in_file'),
        ]),
        (inputnode, t1w_aslref, [
            ('t1w_preproc', 'input_image'),
            ('coreg_aslref', 'reference_image'),
            ('aslref2anat_xfm', 'transforms'),
        ]),
        (inputnode, t1w_wm, [('t1w_dseg', 'in_seg')]),
        (inputnode, aslref_wm, [
            ('coreg_aslref', 'reference_image'),
            ('aslref2anat_xfm', 'transforms'),
        ]),
        (t1w_wm, aslref_wm, [('out', 'input_image')]),
    ])  # fmt:skip

    # Reportlets follow the structure of init_asl_fit_wf stages
    # - SDC1:
    #       Before: Pre-SDC aslref
    #       After: Fieldmap reference resampled on aslref
    #       Three-way: Fieldmap resampled on aslref
    # - SDC2:
    #       Before: Pre-SDC aslref with white matter mask
    #       After: Post-SDC aslref with white matter mask
    # - EPI-T1 registration:
    #       Before: T1w brain with white matter mask
    #       After: Resampled aslref with white matter mask

    if sdc_correction:
        fmapref_aslref = pe.Node(
            ApplyTransforms(
                dimension=3,
                default_value=0,
                float=True,
                invert_transform_flags=[True],
                interpolation='LanczosWindowedSinc',
                args='-v',
            ),
            name='fmapref_aslref',
            mem_gb=1,
        )

        fmap_aslref = pe.Node(
            ApplyTransforms(
                dimension=3,
                default_value=0,
                float=True,
                invert_transform_flags=[True],
                interpolation='LanczosWindowedSinc',
            ),
            name='fmap_aslref',
            mem_gb=1,
        )

        # SDC1
        sdcreg_report = pe.Node(
            FieldmapReportlet(
                reference_label='ASL reference',
                moving_label='Fieldmap reference',
                show='both',
            ),
            name='sdcreg_report',
            mem_gb=0.1,
        )

        ds_sdcreg_report = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc='fmapCoreg',
                suffix='asl',
                datatype='figures',
                dismiss_entities=('echo',),
            ),
            name='ds_sdcreg_report',
        )

        # SDC2
        sdc_report = pe.Node(
            SimpleBeforeAfter(
                before_label='Distorted',
                after_label='Corrected',
                dismiss_affine=True,
            ),
            name='sdc_report',
            mem_gb=0.1,
        )

        ds_sdc_report = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc='sdc',
                suffix='asl',
                datatype='figures',
                dismiss_entities=('echo',),
            ),
            name='ds_sdc_report',
        )

        workflow.connect([
            (inputnode, fmapref_aslref, [
                ('fmap_ref', 'input_image'),
                ('coreg_aslref', 'reference_image'),
                ('aslref2fmap_xfm', 'transforms'),
            ]),
            (inputnode, fmap_aslref, [
                ('fieldmap', 'input_image'),
                ('coreg_aslref', 'reference_image'),
                ('aslref2fmap_xfm', 'transforms'),
            ]),
            (inputnode, sdcreg_report, [
                ('sdc_aslref', 'reference'),
                ('asl_mask', 'mask'),
            ]),
            (fmapref_aslref, sdcreg_report, [('output_image', 'moving')]),
            (fmap_aslref, sdcreg_report, [('output_image', 'fieldmap')]),
            (inputnode, ds_sdcreg_report, [('source_file', 'source_file')]),
            (sdcreg_report, ds_sdcreg_report, [('out_report', 'in_file')]),
            (inputnode, sdc_report, [
                ('sdc_aslref', 'before'),
                ('coreg_aslref', 'after'),
            ]),
            (aslref_wm, sdc_report, [('output_image', 'wm_seg')]),
            (inputnode, ds_sdc_report, [('source_file', 'source_file')]),
            (sdc_report, ds_sdc_report, [('out_report', 'in_file')]),
        ])  # fmt:skip

    # EPI-T1 registration
    # Resample T1w image onto EPI-space

    epi_t1_report = pe.Node(
        SimpleBeforeAfter(
            before_label='T1w',
            after_label='EPI',
            dismiss_affine=True,
        ),
        name='epi_t1_report',
        mem_gb=0.1,
    )

    ds_epi_t1_report = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='coreg',
            suffix='asl',
            datatype='figures',
            dismiss_entities=('echo',),
        ),
        name='ds_epi_t1_report',
    )

    workflow.connect([
        (inputnode, epi_t1_report, [('coreg_aslref', 'after')]),
        (t1w_aslref, epi_t1_report, [('output_image', 'before')]),
        (aslref_wm, epi_t1_report, [('output_image', 'wm_seg')]),
        (inputnode, ds_epi_t1_report, [('source_file', 'source_file')]),
        (epi_t1_report, ds_epi_t1_report, [('out_report', 'in_file')]),
    ])  # fmt:skip

    return workflow


def init_ds_aslref_wf(
    *,
    bids_root,
    output_dir,
    desc: str,
    name='ds_aslref_wf',
) -> pe.Workflow:
    """Write out aslref image."""
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['source_files', 'aslref']),
        name='inputnode',
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=['aslref']), name='outputnode')

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name='raw_sources')
    raw_sources.inputs.bids_root = bids_root

    ds_aslref = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc=desc,
            suffix='aslref',
            compress=True,
            dismiss_entities=('echo',),
        ),
        name='ds_aslref',
        run_without_submitting=True,
    )

    workflow.connect([
        (inputnode, raw_sources, [('source_files', 'in_files')]),
        (inputnode, ds_aslref, [
            ('aslref', 'in_file'),
            ('source_files', 'source_file'),
        ]),
        (raw_sources, ds_aslref, [('out', 'RawSources')]),
        (ds_aslref, outputnode, [('out_file', 'aslref')]),
    ])  # fmt:skip

    return workflow


def init_ds_asl_native_wf(
    *,
    bids_root: str,
    output_dir: str,
    asl_output: bool,
    metadata: list[dict],
    cbf_3d: list[str],
    cbf_4d: list[str],
    att: list[str],
    name='ds_asl_native_wf',
) -> pe.Workflow:
    """Write out aslref-space outputs."""
    workflow = pe.Workflow(name=name)

    inputnode_fields = [
        'source_files',
        'asl',
    ]
    inputnode_fields += cbf_3d
    inputnode_fields += cbf_4d
    inputnode_fields += att
    inputnode = pe.Node(
        niu.IdentityInterface(fields=inputnode_fields),
        name='inputnode',
    )

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name='raw_sources')
    raw_sources.inputs.bids_root = bids_root
    workflow.connect([(inputnode, raw_sources, [('source_files', 'in_files')])])

    datasinks = []
    # Write out CBF and ATT maps in aslref space
    for cbf_name in cbf_4d + cbf_3d:
        # TODO: Add EstimationReference and EstimationAlgorithm
        fields = BASE_INPUT_FIELDS[cbf_name]

        ds_cbf = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                compress=True,
                dismiss_entities=('echo',),
                **fields,
            ),
            name=f'ds_{cbf_name}',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        datasinks.append(ds_cbf)
        workflow.connect([(inputnode, ds_cbf, [(cbf_name, 'in_file')])])

    for att_name in att:
        # TODO: Add EstimationReference and EstimationAlgorithm
        fields = BASE_INPUT_FIELDS[att_name]

        ds_att = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                compress=True,
                dismiss_entities=('echo',),
                **fields,
            ),
            name=f'ds_{att_name}',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        datasinks.append(ds_att)

        workflow.connect([(inputnode, ds_att, [(att_name, 'in_file')])])

    if asl_output:
        # Write out the preprocessed ASL time series
        ds_asl = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc='preproc',
                compress=True,
                SkullStripped=False,
                dismiss_entities=('echo',),
                **metadata,
            ),
            name='ds_asl',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(inputnode, ds_asl, [('asl', 'in_file')])])
        datasinks.append(ds_asl)

    workflow.connect(
        [(inputnode, datasink, [('source_files', 'source_file')]) for datasink in datasinks] +
        [(raw_sources, datasink, [('out', 'RawSources')]) for datasink in datasinks]
    )  # fmt:skip

    return workflow


def init_ds_volumes_wf(
    *,
    bids_root: str,
    output_dir: str,
    metadata: list[dict],
    cbf_3d: list[str],
    cbf_4d: list[str],
    att: list[str],
    name: str = 'ds_volumes_wf',
) -> pe.Workflow:
    """Apply transforms from reference to anatomical/standard space and write out derivatives."""
    workflow = pe.Workflow(name=name)
    inputnode_fields = [
        'source_files',
        'ref_file',
        'asl',  # Resampled into target space
        'asl_mask',  # aslref space
        'aslref',  # aslref space
        # Anatomical
        'aslref2anat_xfm',
        # Template
        'anat2std_xfm',
        # Entities
        'space',
        'cohort',
        'resolution',
    ]
    inputnode_fields += cbf_3d
    inputnode_fields += cbf_4d
    inputnode_fields += att
    inputnode = pe.Node(
        niu.IdentityInterface(fields=inputnode_fields),
        name='inputnode',
    )

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name='raw_sources')
    raw_sources.inputs.bids_root = bids_root
    aslref2target = pe.Node(niu.Merge(2), name='aslref2target')

    # BOLD is pre-resampled
    ds_asl = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='preproc',
            compress=True,
            SkullStripped=True,
            dismiss_entities=('echo',),
            **metadata,
        ),
        name='ds_asl',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([
        (inputnode, raw_sources, [('source_files', 'in_files')]),
        # Note that ANTs expects transforms in target-to-source order
        # Reverse this for nitransforms-based resamplers
        (inputnode, aslref2target, [
            ('anat2std_xfm', 'in1'),
            ('aslref2anat_xfm', 'in2'),
        ]),
        (inputnode, ds_asl, [
            ('source_files', 'source_file'),
            ('asl', 'in_file'),
            ('space', 'space'),
            ('cohort', 'cohort'),
            ('resolution', 'resolution'),
        ]),
    ])  # fmt:skip

    resample_ref = pe.Node(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            float=True,
            interpolation='LanczosWindowedSinc',
            args='-v',
        ),
        name='resample_ref',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    resample_mask = pe.Node(
        ApplyTransforms(interpolation='GenericLabel', args='-v'),
        name='resample_mask',
    )
    resamplers = [resample_ref, resample_mask]

    workflow.connect([
        (inputnode, resample_ref, [('aslref', 'input_image')]),
        (inputnode, resample_mask, [('asl_mask', 'input_image')]),
    ])  # fmt:skip

    ds_ref = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            suffix='aslref',
            compress=True,
            dismiss_entities=('echo',),
        ),
        name='ds_ref',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    ds_mask = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='brain',
            suffix='mask',
            compress=True,
            dismiss_entities=('echo',),
        ),
        name='ds_mask',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    datasinks = [ds_ref, ds_mask]

    for cbf_name in cbf_4d + cbf_3d:
        # TODO: Add EstimationReference and EstimationAlgorithm
        fields = BASE_INPUT_FIELDS[cbf_name]

        kwargs = {}
        if cbf_name in cbf_4d:
            kwargs['dimension'] = 3

        resample_cbf = pe.Node(
            ApplyTransforms(
                interpolation='LanczosWindowedSinc',
                float=True,
                input_image_type=3,
                args='-v',
                **kwargs,
            ),
            name=f'warp_{cbf_name}_to_std',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        ds_cbf = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                compress=True,
                dismiss_entities=('echo',),
                **fields,
            ),
            name=f'ds_{cbf_name}',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        resamplers.append(resample_cbf)
        datasinks.append(ds_cbf)

        workflow.connect([(inputnode, resample_cbf, [(cbf_name, 'input_image')])])

    for att_name in att:
        # TODO: Add EstimationReference and EstimationAlgorithm
        fields = BASE_INPUT_FIELDS[att_name]

        resample_att = pe.Node(
            ApplyTransforms(
                interpolation='LanczosWindowedSinc',
                float=True,
                input_image_type=3,
                args='-v',
            ),
            name=f'warp_{att_name}_to_std',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        ds_att = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                compress=True,
                dismiss_entities=('echo',),
                **fields,
            ),
            name=f'ds_{att_name}',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        resamplers.append(resample_att)
        datasinks.append(ds_att)

        workflow.connect([(inputnode, resample_att, [(att_name, 'input_image')])])

    workflow.connect(
        [
            (inputnode, resampler, [('ref_file', 'reference_image')])
            for resampler in resamplers
        ] + [
            (aslref2target, resampler, [('out', 'transforms')])
            for resampler in resamplers
        ] + [
            (inputnode, datasink, [
                ('source_files', 'source_file'),
                ('space', 'space'),
                ('cohort', 'cohort'),
                ('resolution', 'resolution'),
            ])
            for datasink in datasinks
        ] + [
            (resampler, datasink, [('output_image', 'in_file')])
            for resampler, datasink in zip(resamplers, datasinks, strict=False)
        ]
    )  # fmt:skip

    return workflow


def init_ds_ciftis_wf(
    *,
    bids_root: str,
    output_dir: str,
    metadata: list[dict],
    cbf_3d: list[str],
    cbf_4d: list[str],
    att: list[str],
    omp_nthreads: int,
    name: str = 'ds_ciftis_wf',
) -> pe.Workflow:
    """Apply transforms from reference to fsLR space and write out derivatives."""
    from fmriprep.workflows.bold.resampling import (
        init_bold_fsLR_resampling_wf,
        init_bold_grayords_wf,
    )

    workflow = pe.Workflow(name=name)
    inputnode_fields = [
        'asl_cifti',
        'source_files',
        # ASL-resolution, anatomical-space reference image
        'anat_ref_file',
        'aslref2anat_xfm',
        # Template
        'anat2mni6_xfm',
        # Template reference image. Resolution depends on CIFTI output resolution
        'mni6_mask',
        # Pre-computed goodvoxels mask. May be Undefined.
        'goodvoxels_mask',
        # Other inputs
        'white',
        'pial',
        'midthickness',
        'midthickness_fsLR',
        'sphere_reg_fsLR',
        'cortex_mask',
    ]
    inputnode_fields += cbf_3d
    inputnode_fields += cbf_4d
    inputnode_fields += att
    inputnode = pe.Node(
        niu.IdentityInterface(fields=inputnode_fields),
        name='inputnode',
    )

    outputnode_fields = []
    outputnode_fields += cbf_3d
    outputnode_fields += cbf_4d
    outputnode_fields += att
    outputnode = pe.Node(
        niu.IdentityInterface(fields=outputnode_fields),
        name='outputnode',
    )

    raw_sources = pe.Node(niu.Function(function=_bids_relative), name='raw_sources')
    raw_sources.inputs.bids_root = bids_root
    workflow.connect([(inputnode, raw_sources, [('source_files', 'in_files')])])

    ds_asl_cifti = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            space='fsLR',
            density=config.workflow.cifti_output,
            suffix='asl',
            extension='dtseries.nii',
            compress=False,
        ),
        name='ds_asl_cifti',
        run_without_submitting=True,
    )
    workflow.connect([
        (inputnode, ds_asl_cifti, [
            ('asl_cifti', 'in_file'),
            ('source_files', 'source_file'),
        ]),
    ])  # fmt:skip

    aslref2MNI6 = pe.Node(niu.Merge(2), name='aslref2MNI6')
    workflow.connect([
        (inputnode, aslref2MNI6, [
            ('aslref2anat_xfm', 'in1'),
            ('anat2mni6_xfm', 'in2'),
        ]),
    ])  # fmt:skip

    for cbf_deriv in cbf_4d + cbf_3d + att:
        kwargs = {}
        extension = 'dscalar.nii'
        if cbf_deriv in cbf_4d:
            kwargs['dimension'] = 3
            extension = 'dtseries.nii'

        warp_cbf_to_anat = pe.Node(
            ApplyTransforms(
                interpolation='LanczosWindowedSinc',
                float=True,
                input_image_type=3,
                args='-v',
                **kwargs,
            ),
            name=f'warp_{cbf_deriv}_to_anat',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([
            (inputnode, warp_cbf_to_anat, [
                (cbf_deriv, 'input_image'),
                ('anat_ref_file', 'reference_image'),
                ('aslref2anat_xfm', 'transforms'),
            ]),
        ])  # fmt:skip

        warp_cbf_to_MNI6 = pe.Node(
            ApplyTransforms(
                interpolation='LanczosWindowedSinc',
                float=True,
                input_image_type=3,
                args='-v',
                **kwargs,
            ),
            name=f'warp_{cbf_deriv}_to_MNI6',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([
            (inputnode, warp_cbf_to_MNI6, [
                ('mni6_mask', 'reference_image'),
                (cbf_deriv, 'input_image'),
            ]),
            (aslref2MNI6, warp_cbf_to_MNI6, [('out', 'transforms')]),
        ])  # fmt:skip

        cbf_fsLR_resampling_wf = init_bold_fsLR_resampling_wf(
            grayord_density=config.workflow.cifti_output,
            omp_nthreads=omp_nthreads,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            name=f'{cbf_deriv}_fsLR_resampling_wf',
        )
        workflow.connect([
            # Resample T1w-space CBF to fsLR surfaces
            (inputnode, cbf_fsLR_resampling_wf, [
                ('white', 'inputnode.white'),
                ('pial', 'inputnode.pial'),
                ('midthickness', 'inputnode.midthickness'),
                ('midthickness_fsLR', 'inputnode.midthickness_fsLR'),
                ('sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR'),
                ('cortex_mask', 'inputnode.cortex_mask'),
                ('goodvoxels_mask', 'inputnode.volume_roi'),
            ]),
            (warp_cbf_to_anat, cbf_fsLR_resampling_wf, [('output_image', 'inputnode.bold_file')])
        ])  # fmt:skip

        cbf_grayords_wf = init_bold_grayords_wf(
            grayord_density=config.workflow.cifti_output,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            repetition_time=metadata['RepetitionTime'],
            name=f'{cbf_deriv}_grayords_wf',
        )
        workflow.connect([
            (warp_cbf_to_MNI6, cbf_grayords_wf, [('output_image', 'inputnode.bold_std')]),
            (cbf_fsLR_resampling_wf, cbf_grayords_wf, [
                ('outputnode.bold_fsLR', 'inputnode.bold_fsLR'),
            ]),
        ])  # fmt:skip

        ds_cbf_cifti = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                space='fsLR',
                density=config.workflow.cifti_output,
                extension=extension,
                compress=False,
                **BASE_INPUT_FIELDS[cbf_deriv],
            ),
            name=f'ds_{cbf_deriv}_cifti',
            run_without_submitting=True,
        )
        workflow.connect([
            (inputnode, ds_cbf_cifti, [('source_files', 'source_file')]),
            (raw_sources, ds_cbf_cifti, [('out', 'RawSources')]),
            (cbf_grayords_wf, ds_cbf_cifti, [
                ('outputnode.cifti_bold', 'in_file'),
                (('outputnode.cifti_metadata', _read_json), 'meta_dict'),
            ]),
            (ds_cbf_cifti, outputnode, [('out_file', cbf_deriv)])
        ])  # fmt:skip

    return workflow


def _read_json(in_file):
    from json import loads
    from pathlib import Path

    if not isinstance(in_file, str):
        raise ValueError(f'_read_json: input is not str ({in_file})')

    return loads(Path(in_file).read_text())
