# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for writing out derivative files."""

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.utils.images import dseg_label

from aslprep import config
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.bids import BIDSURI, DerivativesDataSink
from aslprep.interfaces.pvc import DivideByWeight

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
        'Description': 'Volume-wise cerebral blood flow, estimated by ASLPrep.',
    },
    'mean_cbf': {
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
        'Description': 'Cerebral blood flow, estimated by ASLPrep.',
    },
    'att': {
        'suffix': 'att',
        'Units': 's',
        'Description': (
            'Arterial transit time. The time from labeling until the labeled blood water first '
            'arrives in the voxel. '
            'It depends on where the labeling plane or slab was placed. '
            'Where a range of values fits the data equally well, this is one value from the range, '
            'and the ends of the range are in the desc-lower and desc-upper maps.'
        ),
    },
    'abat': {
        'suffix': 'abat',
        'Units': 's',
        'Description': (
            'Arterial bolus arrival time. '
            'The time from labeling until the labeled blood water first arrives in the arteries '
            'within the voxel. This is not the tissue arterial transit time (att). '
            'The data only locate it within an interval set by the sampled delays, '
            'so this is the midpoint of that interval, '
            'and the ends of the interval are in the desc-lower and desc-upper maps.'
        ),
    },
    'abv': {
        'suffix': 'abv',
        'Units': 'fraction',
        'Description': (
            'Arterial blood volume. The fraction of the voxel occupied by the arterial vessel. '
            'For (P)CASL data, this value depends on the reported arterial bolus arrival time.'
        ),
    },
    'att_lower': {
        'desc': 'lower',
        'suffix': 'att',
        'Units': 's',
        'Description': (
            'Lower bound of the arterial transit time. '
            'The shortest arterial transit time that fits the data as well as the estimate.'
        ),
    },
    'att_upper': {
        'desc': 'upper',
        'suffix': 'att',
        'Units': 's',
        'Description': (
            'Upper bound of the arterial transit time. '
            'The longest arterial transit time that fits the data as well as the estimate.'
        ),
    },
    'abat_lower': {
        'desc': 'lower',
        'suffix': 'abat',
        'Units': 's',
        'Description': (
            'Lower bound of the arterial bolus arrival time. '
            'The shortest arterial bolus arrival time that fits the data as well as the estimate.'
        ),
    },
    'abat_upper': {
        'desc': 'upper',
        'suffix': 'abat',
        'Units': 's',
        'Description': (
            'Upper bound of the arterial bolus arrival time. '
            'The longest arterial bolus arrival time that fits the data as well as the estimate.'
        ),
    },
    'abv_decayed': {
        'desc': 'decayed',
        'suffix': 'abv',
        'Units': 'fraction',
        'Description': (
            'Arterial blood volume multiplied by exp(-aBAT / T1blood), '
            'the decay of the label over the arterial bolus arrival time. '
            'Unlike the arterial blood volume itself, '
            'this product does not depend on where the arterial bolus arrival time falls '
            'between the two sampled delays that enclose it.'
        ),
    },
    # SCORE/SCRUB outputs
    'cbf_ts_score': {
        'desc': 'scoreTimeseries',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
        'Description': (
            'Volume-wise cerebral blood flow after SCORE denoising, estimated by ASLPrep.'
        ),
    },
    'mean_cbf_score': {
        'desc': 'score',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
        'Description': 'Cerebral blood flow after SCORE denoising, estimated by ASLPrep.',
    },
    'mean_cbf_scrub': {
        'desc': 'scrub',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
        'Description': (
            'Cerebral blood flow after SCORE and SCRUB denoising, estimated by ASLPrep.'
        ),
    },
    # BASIL outputs
    'mean_cbf_basil': {
        'desc': 'basil',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
        'Description': 'Cerebral blood flow, estimated by BASIL.',
    },
    'mean_cbf_gm_basil': {
        'desc': 'basilGM',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
        'Description': (
            'Cerebral blood flow after gray matter partial volume correction, estimated by BASIL.'
        ),
    },
    'mean_cbf_wm_basil': {
        'desc': 'basilWM',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
        'Description': (
            'Cerebral blood flow after white matter partial volume correction, estimated by BASIL.'
        ),
    },
    'att_basil': {
        'desc': 'basil',
        'suffix': 'att',
        'Units': 's',
        'Description': (
            'Arterial transit time, estimated by BASIL. '
            'The time from labeling until the labeled blood water first arrives in the voxel. '
            'It depends on where the labeling plane or slab was placed.'
        ),
    },
    'abv_basil': {
        'desc': 'basil',
        'suffix': 'abv',
        'Units': 'fraction',
        'Description': (
            'Arterial (cerebral) blood volume, estimated by BASIL. '
            'The fraction of the voxel occupied by the arterial vessel.'
        ),
    },
    # BASIL posterior standard deviations
    'cbf_sd_basil': {
        'desc': 'basilSD',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
        'Description': (
            'Posterior standard deviation of the cerebral blood flow estimated by BASIL '
            '(desc-basil). The square root of the posterior variance written by BASIL.'
        ),
    },
    'cbf_gm_sd_basil': {
        'desc': 'basilGMSD',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
        'Description': (
            'Posterior standard deviation of the gray matter partial-volume corrected '
            'cerebral blood flow estimated by BASIL (desc-basilGM).'
        ),
    },
    'cbf_wm_sd_basil': {
        'desc': 'basilWMSD',
        'suffix': 'cbf',
        'Units': 'mL/100 g/min',
        'Description': (
            'Posterior standard deviation of the white matter partial-volume corrected '
            'cerebral blood flow estimated by BASIL (desc-basilWM).'
        ),
    },
    'att_sd_basil': {
        'desc': 'basilSD',
        'suffix': 'att',
        'Units': 's',
        'Description': (
            'Posterior standard deviation of the arterial transit time estimated by BASIL '
            '(desc-basil).'
        ),
    },
    # Partial volume maps used by local-regression partial volume correction
    'gm_pv': {
        'label': 'GM',
        'desc': 'psf',
        'suffix': 'probseg',
        'Description': (
            'Gray matter fraction on the ASL reference grid, modelling the effective '
            'resolution of the ASL acquisition.'
        ),
    },
    'wm_pv': {
        'label': 'WM',
        'desc': 'psf',
        'suffix': 'probseg',
        'Description': (
            'White matter fraction on the ASL reference grid, modelling the effective '
            'resolution of the ASL acquisition.'
        ),
    },
    'csf_pv': {
        'label': 'CSF',
        'desc': 'psf',
        'suffix': 'probseg',
        'Description': (
            'CSF fraction on the ASL reference grid, modelling the effective resolution of the '
            'ASL acquisition.'
        ),
    },
}

#: Local-regression partial volume correction outputs, per corrected CBF map.
#: Keys are the corrected map's name in the CBF workflow; values are the desc prefix.
LR_PVC_SOURCES = {'mean_cbf': 'lr', 'mean_cbf_score': 'scoreLR'}


def lr_pvc_names(source):
    """Return the output names of local-regression PVC applied to one CBF map.

    Parameters
    ----------
    source : :obj:`str`
        ``'mean_cbf'`` or ``'mean_cbf_score'``.

    Returns
    -------
    :obj:`dict`
        ``cbf`` (corrected CBF maps), ``se`` (their standard errors), ``valid`` (where they
        were estimated, by corrected map), ``diagnostics`` (native-only fit diagnostics), and
        ``csf`` (the native-only CSF quality-control map, its standard error, fit status and
        noise amplification).
    """
    suffix = '' if source == 'mean_cbf' else '_score'
    return {
        'csf': [
            f'mean_cbf_csf_lr{suffix}',
            f'cbf_csf_se_lr{suffix}',
            f'status_csf_lr{suffix}',
            f'csf_noise_lr{suffix}',
        ],
        'cbf': [f'mean_cbf_gm_lr{suffix}', f'mean_cbf_wm_lr{suffix}'],
        'se': [f'cbf_gm_se_lr{suffix}', f'cbf_wm_se_lr{suffix}'],
        'valid': {
            f'mean_cbf_gm_lr{suffix}': f'gm_valid_lr{suffix}',
            f'mean_cbf_wm_lr{suffix}': f'wm_valid_lr{suffix}',
        },
        'diagnostics': [
            f'status_lr{suffix}',
            f'condition_lr{suffix}',
            f'gm_noise_lr{suffix}',
            f'wm_noise_lr{suffix}',
            f'rmse_lr{suffix}',
        ],
    }


def _add_lr_pvc_fields():
    for source, prefix in LR_PVC_SOURCES.items():
        names = lr_pvc_names(source)
        source_desc = BASE_INPUT_FIELDS[source].get('desc')
        source_str = f' (desc-{source_desc})' if source_desc else ''
        tissue_names = {'GM': 'gray matter', 'WM': 'white matter'}
        for name, se_name, tissue in zip(names['cbf'], names['se'], ('GM', 'WM'), strict=True):
            BASE_INPUT_FIELDS[name] = {
                'desc': f'{prefix}{tissue}',
                'suffix': 'cbf',
                'Units': 'mL/100 g/min',
                'Description': (
                    f'{tissue_names[tissue].capitalize()} cerebral blood flow after partial '
                    f'volume correction by local linear regression of the CBF map{source_str}. '
                    'Zero where it was not estimated.'
                ),
            }
            BASE_INPUT_FIELDS[se_name] = {
                'desc': f'{prefix}{tissue}SE',
                'suffix': 'cbf',
                'Units': 'mL/100 g/min',
                'Description': (
                    f'Standard error of the {tissue_names[tissue]} CBF estimated by local linear '
                    f'regression (desc-{prefix}{tissue}), assuming independent, identically '
                    'distributed noise and conditional on M0. '
                    'Valid in this (native) space only.'
                ),
            }
        status, condition, gm_noise, wm_noise, rmse = names['diagnostics']
        BASE_INPUT_FIELDS[status] = {
            'desc': prefix,
            'suffix': 'pvstatus',
            'Description': 'Outcome of the local-regression partial volume correction fit.',
            'Levels': {
                '0': 'Outside the brain mask.',
                '1': 'Gray and white matter fit.',
                '2': 'Gray matter only fit (white matter contamination below the limit).',
                '3': 'White matter only fit (gray matter contamination below the limit).',
                '4': 'Not identifiable: too little tissue support, ill-conditioned, or no M0.',
            },
        }
        BASE_INPUT_FIELDS[condition] = {
            'desc': prefix,
            'suffix': 'pvcondition',
            'Description': (
                'Condition number of the normalized gray/white matter Gram matrix of each '
                'local regression. Zero where undefined.'
            ),
        }
        for name, tissue in ((gm_noise, 'GM'), (wm_noise, 'WM')):
            BASE_INPUT_FIELDS[name] = {
                'desc': f'{prefix}{tissue}',
                'suffix': 'pvnoise',
                'Description': (
                    f'Noise amplification of the {tissue_names[tissue]} estimate: its standard '
                    'error per unit of voxel noise, assuming independent noise.'
                ),
            }
        BASE_INPUT_FIELDS[rmse] = {
            'desc': prefix,
            'suffix': 'pvrmse',
            'Description': (
                'Root mean square residual of each local regression, in units of CBF times M0.'
            ),
        }
        csf_cbf, csf_se, csf_status, csf_noise = names['csf']
        csf_fit_str = (
            'from a separate local regression on gray matter, white matter and CSF partial '
            f'volume maps of the CBF map{source_str}, for quality control'
        )
        BASE_INPUT_FIELDS[csf_cbf] = {
            'desc': f'{prefix}CSF',
            'suffix': 'cbf',
            'Units': 'mL/100 g/min',
            'Description': (
                f'Apparent CSF cerebral blood flow, {csf_fit_str}. CSF has no perfusion, so '
                'values far from zero indicate partial volume map errors, macrovascular or '
                'choroid plexus signal, or model error. The gray and white matter maps '
                f'(desc-{prefix}GM, desc-{prefix}WM) do not come from this fit. '
                'Zero where it was not estimated. Valid in this (native) space only.'
            ),
        }
        BASE_INPUT_FIELDS[csf_se] = {
            'desc': f'{prefix}CSFSE',
            'suffix': 'cbf',
            'Units': 'mL/100 g/min',
            'Description': (
                f'Standard error of desc-{prefix}CSF, assuming independent, identically '
                'distributed noise and conditional on M0. Valid in this (native) space only.'
            ),
        }
        BASE_INPUT_FIELDS[csf_status] = {
            'desc': f'{prefix}CSF',
            'suffix': 'pvstatus',
            'Description': f'Outcome of the CSF quality-control fit (desc-{prefix}CSF).',
            'Levels': {
                '0': 'Outside the brain mask.',
                '1': 'Gray matter, white matter and CSF fit.',
                '2': 'Gray matter and CSF fit (white matter contamination below the limit).',
                '3': 'White matter and CSF fit (gray matter contamination below the limit).',
                '4': 'Not identifiable.',
                '5': 'CSF only fit (gray and white matter contamination below the limit).',
            },
        }
        BASE_INPUT_FIELDS[csf_noise] = {
            'desc': f'{prefix}CSF',
            'suffix': 'pvnoise',
            'Description': (
                'Noise amplification of the CSF estimate: its standard error per unit of voxel '
                'noise, assuming independent noise.'
            ),
        }


_add_lr_pvc_fields()


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
    source_file: str,
    sdc_correction: bool,
    separate_m0scan: bool,
    freesurfer: bool,
    output_dir: str,
    name='asl_fit_reports_wf',
) -> pe.Workflow:
    """Set up a battery of datasinks to store reports in the right location.

    Copied from fMRIPrep's init_bold_fit_reports_wf.
    Modifications include changes to fields and variable names, which aren't important,
    and DerivativesDataSink suffixes, which are.

    Parameters
    ----------
    source_file : :obj:`str`
        Input ASL file, used to name the reports.
    sdc_correction : :obj:`bool`
        Whether SDC correction was performed or not.
    separate_m0scan : :obj:`bool`
        Whether a separate M0 scan was used or not.
    freesurfer : :obj:`bool`
        FreeSurfer was enabled
    output_dir : :obj:`str`
        Directory in which to save derivatives
    name : :obj:`str`
        Workflow name (default: anat_reports_wf)

    Inputs
    ------
    sdc_aslref
        ASL reference before SDC, in ASL space
    coreg_aslref
        ASL reference after SDC, in ASL space
    m0scan_aslref
        Mean of the separate M0 scan, in ASL reference space
    aslref2anat_xfm
        Affine transform from ASL reference to anatomical space
    aslref2fmap_xfm
        Affine transform from ASL reference to fieldmap space
    t1w_preproc
        The T1w reference map, which is calculated as the average of bias-corrected
        and preprocessed T1w images, defining the anatomical space.
    t1w_mask
        Brain (binary) mask estimated by brain extraction.
    t1w_dseg
        Segmentation in T1w space
    asl_mask
        Brain mask in ASL reference space
    fieldmap
        Reconstructed fieldmap, in ASL reference space
    fmap_ref
        Fieldmap reference image, in fieldmap space
    subject_id
        FreeSurfer subject ID
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    summary_report
        HTML snippet summarizing ASL processing decisions
    validation_report
        HTML snippet indicating whether affine metadata were overwritten

    """
    from nireports.interfaces.reporting.base import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )
    from sdcflows.interfaces.reportlets import FieldmapReportlet

    workflow = pe.Workflow(name=name)

    inputfields = [
        'sdc_aslref',
        'coreg_aslref',
        'm0scan_aslref',
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
            source_file=source_file,
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
            source_file=source_file,
            base_directory=output_dir,
            desc='validation',
            datatype='figures',
            dismiss_entities=('echo',),
        ),
        name='ds_report_validation',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # Resample anatomical references into ASL space for plotting
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
        (inputnode, ds_summary, [('summary_report', 'in_file')]),
        (inputnode, ds_validation, [('validation_report', 'in_file')]),
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
    # - M0-ASL registration (if applicable):
    #       Before: M0 scan
    #       After: ASL reference
    # - ASL-T1 registration:
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
                source_file=source_file,
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
                source_file=source_file,
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
            (inputnode, sdcreg_report, [
                ('sdc_aslref', 'reference'),
                ('fieldmap', 'fieldmap'),
                ('asl_mask', 'mask'),
            ]),
            (fmapref_aslref, sdcreg_report, [('output_image', 'moving')]),
            (sdcreg_report, ds_sdcreg_report, [('out_report', 'in_file')]),
            (inputnode, sdc_report, [
                ('sdc_aslref', 'before'),
                ('coreg_aslref', 'after'),
            ]),
            (aslref_wm, sdc_report, [('output_image', 'wm_seg')]),
            (sdc_report, ds_sdc_report, [('out_report', 'in_file')]),
        ])  # fmt:skip

    # M0-ASL registration
    #       Before: M0 scan
    #       After: ASL reference

    if separate_m0scan:
        m0_asl_report = pe.Node(
            SimpleBeforeAfter(
                before_label='M0',
                after_label='ASL',
                dismiss_affine=True,
            ),
            name='m0_asl_report',
            mem_gb=0.1,
        )

        ds_m0_asl_report = pe.Node(
            DerivativesDataSink(
                source_file=source_file,
                base_directory=output_dir,
                desc='m0reg',
                suffix='asl',
                datatype='figures',
                dismiss_entities=('echo',),
            ),
            name='ds_m0_asl_report',
        )

        workflow.connect([
            (inputnode, m0_asl_report, [
                ('m0scan_aslref', 'before'),
                ('coreg_aslref', 'after'),
            ]),
            (aslref_wm, m0_asl_report, [('output_image', 'wm_seg')]),
            (m0_asl_report, ds_m0_asl_report, [('out_report', 'in_file')]),
        ])  # fmt:skip

    # ASL-T1 registration
    # Resample T1w image onto ASL-space

    epi_t1_report = pe.Node(
        SimpleBeforeAfter(
            before_label='T1w',
            after_label='ASL',
            dismiss_affine=True,
        ),
        name='epi_t1_report',
        mem_gb=0.1,
    )

    ds_epi_t1_report = pe.Node(
        DerivativesDataSink(
            source_file=source_file,
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
        (epi_t1_report, ds_epi_t1_report, [('out_report', 'in_file')]),
    ])  # fmt:skip

    return workflow


def init_ds_aslref_wf(
    *,
    source_file: str,
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

    sources = pe.Node(
        BIDSURI(
            numinputs=1,
            dataset_links=config.execution.dataset_links,
            out_dir=str(output_dir),
        ),
        name='sources',
    )

    ds_aslref = pe.Node(
        DerivativesDataSink(
            source_file=source_file,
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
        (inputnode, sources, [('source_files', 'in1')]),
        (inputnode, ds_aslref, [('aslref', 'in_file')]),
        (sources, ds_aslref, [('out', 'Sources')]),
        (ds_aslref, outputnode, [('out_file', 'aslref')]),
    ])  # fmt:skip

    return workflow


def init_ds_asl_native_wf(
    *,
    source_file: str,
    bids_root: str,
    output_dir: str,
    asl_output: bool,
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
        # Metadata for the preprocessed ASL series
        'asl_metadata',
        # Transforms previously used to generate the outputs
        'motion_xfm',
        'aslref2fmap_xfm',
    ]
    inputnode_fields += cbf_3d
    inputnode_fields += cbf_4d
    inputnode_fields += att
    inputnode = pe.Node(
        niu.IdentityInterface(fields=inputnode_fields),
        name='inputnode',
    )

    sources = pe.Node(
        BIDSURI(
            numinputs=3,
            dataset_links=config.execution.dataset_links,
            out_dir=str(output_dir),
        ),
        name='sources',
    )
    workflow.connect([
        (inputnode, sources, [
            ('source_files', 'in1'),
            ('motion_xfm', 'in2'),
            ('aslref2fmap_xfm', 'in3'),
        ]),
    ])  # fmt:skip

    datasinks = []
    # Write out CBF and ATT maps in aslref space
    for cbf_name in cbf_4d + cbf_3d:
        # TODO: Add EstimationReference and EstimationAlgorithm
        fields = BASE_INPUT_FIELDS[cbf_name]

        ds_cbf = pe.Node(
            DerivativesDataSink(
                source_file=source_file,
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
                source_file=source_file,
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
                source_file=source_file,
                base_directory=output_dir,
                desc='preproc',
                compress=True,
                SkullStripped=False,
                dismiss_entities=('echo',),
            ),
            name='ds_asl',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([
            (inputnode, ds_asl, [
                ('asl', 'in_file'),
                ('asl_metadata', 'meta_dict'),
            ]),
        ])  # fmt:skip
        datasinks.append(ds_asl)

    workflow.connect(
        [(sources, datasink, [('out', 'Sources')]) for datasink in datasinks]
    )  # fmt:skip

    return workflow


def init_ds_pvc_native_wf(
    *,
    source_file: str,
    output_dir: str,
    sources: list[str],
    name: str = 'ds_pvc_native_wf',
) -> pe.Workflow:
    """Write out local-regression partial volume correction outputs in aslref space.

    The corrected CBF maps, their standard errors, the fit diagnostics and the PV maps are
    always written in aslref space, because the standard errors, diagnostics and PV maps only
    exist there.

    Parameters
    ----------
    source_file
    output_dir
    sources : list of str
        The corrected CBF maps (``mean_cbf``, and ``mean_cbf_score`` with SCORE).
    name
    """
    workflow = pe.Workflow(name=name)

    names = []
    for source in sources:
        lr_names = lr_pvc_names(source)
        names += lr_names['cbf'] + lr_names['se'] + lr_names['diagnostics'] + lr_names['csf']

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['source_files', 'gm_pv', 'wm_pv', 'csf_pv', 'pv_metadata', 'pvc_metadata']
            + names
        ),
        name='inputnode',
    )

    sources_node = pe.Node(
        BIDSURI(
            numinputs=1,
            dataset_links=config.execution.dataset_links,
            out_dir=str(output_dir),
        ),
        name='sources',
    )
    workflow.connect([(inputnode, sources_node, [('source_files', 'in1')])])

    for out_name in ['gm_pv', 'wm_pv', 'csf_pv'] + names:
        fields = BASE_INPUT_FIELDS[out_name]
        meta_field = 'pv_metadata' if out_name in ('gm_pv', 'wm_pv', 'csf_pv') else 'pvc_metadata'
        ds = pe.Node(
            DerivativesDataSink(
                source_file=source_file,
                base_directory=output_dir,
                compress=True,
                dismiss_entities=('echo',),
                **fields,
            ),
            name=f'ds_{out_name}',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([
            (inputnode, ds, [(out_name, 'in_file'), (meta_field, 'meta_dict')]),
            (sources_node, ds, [('out', 'Sources')]),
        ])  # fmt:skip

    return workflow


def init_ds_volumes_wf(
    *,
    source_file: str,
    bids_root: str,
    output_dir: str,
    cbf_3d: list[str],
    cbf_4d: list[str],
    att: list[str],
    sparse_3d: dict[str, str] | None = None,
    name: str = 'ds_volumes_wf',
) -> pe.Workflow:
    """Apply transforms from reference to anatomical/standard space and write out derivatives.

    Parameters
    ----------
    sparse_3d : dict or None
        Maps that are only estimated in some voxels (zero elsewhere), keyed by name, with the
        name of their validity mask as the value. They are resampled with validity-normalized
        linear interpolation, so the zeros that stand for "not estimated" do not blend into the
        estimates.
    """
    sparse_3d = sparse_3d or {}
    workflow = pe.Workflow(name=name)
    inputnode_fields = [
        'source_files',
        'ref_file',
        'asl',  # Resampled into target space
        'asl_metadata',  # Metadata for the preprocessed ASL series
        'asl_mask',  # aslref space
        'aslref',  # aslref space
        # Anatomical
        'aslref2anat_xfm',
        # Template
        'anat2std_xfm',
        'template',  # target reference image from original transform
        # Transforms previously used to generate the outputs
        'motion_xfm',
        'aslref2fmap_xfm',
        # Entities
        'space',
        'cohort',
        'resolution',
    ]
    inputnode_fields += cbf_3d
    inputnode_fields += cbf_4d
    inputnode_fields += att
    inputnode_fields += list(sparse_3d) + sorted(set(sparse_3d.values()))
    inputnode = pe.Node(
        niu.IdentityInterface(fields=inputnode_fields),
        name='inputnode',
    )

    sources = pe.Node(
        BIDSURI(
            numinputs=6,
            dataset_links=config.execution.dataset_links,
            out_dir=str(output_dir),
        ),
        name='sources',
    )
    aslref2target = pe.Node(niu.Merge(2), name='aslref2target')

    # ASL is pre-resampled
    ds_asl = pe.Node(
        DerivativesDataSink(
            source_file=source_file,
            base_directory=output_dir,
            desc='preproc',
            compress=True,
            SkullStripped=True,
            dismiss_entities=('echo',),
        ),
        name='ds_asl',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([
        (inputnode, sources, [
            ('source_files', 'in1'),
            ('motion_xfm', 'in2'),
            ('aslref2fmap_xfm', 'in3'),
            ('aslref2anat_xfm', 'in4'),
            ('anat2std_xfm', 'in5'),
            ('template', 'in6'),
        ]),
        # Note that ANTs expects transforms in target-to-source order
        # Reverse this for nitransforms-based resamplers
        (inputnode, aslref2target, [
            ('anat2std_xfm', 'in1'),
            ('aslref2anat_xfm', 'in2'),
        ]),
        (sources, ds_asl, [('out', 'Sources')]),
        (inputnode, ds_asl, [
            ('asl', 'in_file'),
            ('asl_metadata', 'meta_dict'),
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
            source_file=source_file,
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
            source_file=source_file,
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
                source_file=source_file,
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
                source_file=source_file,
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

    # Validity-normalized resampling: warp(value) / warp(valid), with linear interpolation.
    # The values are already zero where they are not valid.
    sparse_resamplers = []
    for sparse_name, valid_name in sparse_3d.items():
        fields = BASE_INPUT_FIELDS[sparse_name]
        resample_value = pe.Node(
            ApplyTransforms(interpolation='Linear', float=True, input_image_type=3, args='-v'),
            name=f'warp_{sparse_name}_to_std',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        resample_valid = pe.Node(
            ApplyTransforms(interpolation='Linear', float=True, input_image_type=3, args='-v'),
            name=f'warp_{sparse_name}_valid_to_std',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        normalize = pe.Node(
            DivideByWeight(threshold=0.5),
            name=f'normalize_{sparse_name}',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        ds_sparse = pe.Node(
            DerivativesDataSink(
                source_file=source_file,
                base_directory=output_dir,
                compress=True,
                dismiss_entities=('echo',),
                **fields,
            ),
            name=f'ds_{sparse_name}',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        sparse_resamplers += [resample_value, resample_valid]
        datasinks.append(ds_sparse)
        workflow.connect([
            (inputnode, resample_value, [(sparse_name, 'input_image')]),
            (inputnode, resample_valid, [(valid_name, 'input_image')]),
            (resample_value, normalize, [('output_image', 'in_file')]),
            (resample_valid, normalize, [('output_image', 'weight_file')]),
            (normalize, ds_sparse, [('out_file', 'in_file')]),
        ])  # fmt:skip

    workflow.connect(
        [
            (inputnode, resampler, [('ref_file', 'reference_image')])
            for resampler in sparse_resamplers
        ] + [
            (aslref2target, resampler, [('out', 'transforms')])
            for resampler in sparse_resamplers
        ]
    )  # fmt:skip

    workflow.connect(
        [
            (inputnode, resampler, [('ref_file', 'reference_image')])
            for resampler in resamplers
        ] + [
            (aslref2target, resampler, [('out', 'transforms')])
            for resampler in resamplers
        ] + [
            (sources, datasink, [('out', 'Sources')])
            for datasink in datasinks
        ] + [
            (inputnode, datasink, [
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
    source_file: str,
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
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    workflow = Workflow(name=name)
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

    sources = pe.Node(
        BIDSURI(
            numinputs=1,
            dataset_links=config.execution.dataset_links,
            out_dir=str(output_dir),
        ),
        name='sources',
    )
    workflow.connect([(inputnode, sources, [('source_files', 'in1')])])

    ds_asl_cifti = pe.Node(
        DerivativesDataSink(
            source_file=source_file,
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
        (inputnode, ds_asl_cifti, [('asl_cifti', 'in_file')]),
        (sources, ds_asl_cifti, [('out', 'Sources')]),
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
                source_file=source_file,
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
            (sources, ds_cbf_cifti, [('out', 'Sources')]),
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


def _remove_keys(metadata, keys, updates=None):
    """Return a copy of a metadata dictionary without the given keys, plus any updates.

    Parameters
    ----------
    metadata : dict
        The metadata dictionary.
    keys : list of str
        Keys to remove.
    updates : dict or None, optional
        Keys and values to add or replace after removal.

    Returns
    -------
    dict
        The new metadata dictionary. ``metadata`` is not modified.
    """
    out = {k: v for k, v in metadata.items() if k not in keys}
    out.update(updates or {})
    return out
