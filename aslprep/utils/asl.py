"""Functions for working with ASL data."""

from __future__ import annotations

import numpy as np

# Volume types in the aslcontext file that users may ask ASLPrep to ignore with --ignore.
IGNORABLE_VOLUME_TYPES = ('cbf', 'deltam')


def pcasl_or_pasl(metadata):
    """Determine if metadata indicates a PCASL or ASL scan."""
    aslt = metadata['ArterialSpinLabelingType']

    if aslt in ['CASL', 'PCASL']:
        is_casl = True
    elif aslt == 'PASL':
        is_casl = False
    else:
        raise ValueError(
            'Labeling type cannot be classified as (P)CASL or PASL based on '
            f"ArterialSpinLabelingType: '{aslt}'."
        )

    return is_casl


def determine_multi_pld(metadata):
    """Determine if a run is multi-delay or not.

    Parameters
    ----------
    metadata : :obj:`dict`
        Dictionary of metadata from the ASL file.

    Returns
    -------
    :obj:`bool`
        True if the data are multi-delay/TI. False if not.
    """
    plds = np.array(metadata['PostLabelingDelay'])
    return np.unique(plds).size > 1


def ignore_volume_types(aslcontext, ignore, out_dir):
    """Hide volume types the user asked to ignore from the rest of the workflow.

    The ignored volumes are relabeled in a copy of the aslcontext file,
    so that ASLPrep neither builds the reference image from them nor processes them.
    The ASL file itself is not changed, which keeps the aslcontext file in line with it.

    Parameters
    ----------
    aslcontext : :obj:`str`
        Path to the aslcontext file.
    ignore : :obj:`list` of :obj:`str`
        The ``--ignore`` values. Only the volume types (``cbf`` and ``deltam``) are used.
    out_dir : :obj:`str` or :obj:`pathlib.Path`
        Directory in which to write the modified aslcontext file.

    Returns
    -------
    aslcontext : :obj:`str`
        Path to the modified aslcontext file,
        or the original path if the file has no volumes to ignore.

    Raises
    ------
    ValueError
        If ignoring the volumes leaves nothing that ASLPrep can process.
    """
    from pathlib import Path

    import pandas as pd

    ignored_types = sorted(set(ignore or []) & set(IGNORABLE_VOLUME_TYPES))
    aslcontext_df = pd.read_table(aslcontext)
    to_ignore = aslcontext_df['volume_type'].isin(ignored_types)
    if not to_ignore.any():
        return aslcontext

    aslcontext_df.loc[to_ignore, 'volume_type'] = 'ignored'
    voltypes = aslcontext_df['volume_type'].tolist()
    if not (('control' in voltypes and 'label' in voltypes) or {'deltam', 'cbf'} & set(voltypes)):
        raise ValueError(
            f'Ignoring the {" and ".join(ignored_types)} volumes (--ignore) leaves no '
            f'control-label pairs, deltam volumes or cbf volumes to process in {aslcontext}.'
        )

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / Path(aslcontext).name
    aslcontext_df.to_csv(out_file, sep='\t', index=False)

    return str(out_file)


def select_processing_target(aslcontext):
    """Determine how to handle ASL and M0 data based on dataset configuration."""
    import pandas as pd

    try:
        aslcontext_df = pd.read_table(aslcontext)
    except Exception:  # noqa: BLE001
        raise FileNotFoundError(aslcontext) from None

    voltypes = aslcontext_df['volume_type'].tolist()

    if 'control' in voltypes and 'label' in voltypes:
        processing_target = 'control'
    elif 'deltam' in voltypes:
        processing_target = 'deltam'
    elif 'cbf' in voltypes:
        processing_target = 'cbf'
    else:
        raise ValueError("aslcontext doesn't have control, label, deltam, or cbf volumes.")

    return processing_target


def get_kept_volume_types(processing_target):
    """List the volume types that are kept when the ASL series is reduced.

    Parameters
    ----------
    processing_target : {"control", "deltam", "cbf"}
        The processing target, from :func:`select_processing_target`.

    Returns
    -------
    :obj:`list` of :obj:`str`
        Volume types, as written in the aslcontext file, that are kept.
    """
    if processing_target == 'control':
        return ['control', 'label', 'm0scan']
    elif processing_target == 'deltam':
        return ['deltam', 'm0scan']
    else:
        return ['cbf', 'm0scan']


def count_reduced_volumes(aslcontext, processing_target):
    """Count the volumes in the ASL series after it is reduced to the processing target.

    Workflows downstream of the reduction (e.g., the confounds workflow) see this series,
    not the full ASL file, so they must use this count.

    Parameters
    ----------
    aslcontext : :obj:`str`
        Path to the aslcontext file.
    processing_target : {"control", "deltam", "cbf"}
        The processing target, from :func:`select_processing_target`.

    Returns
    -------
    :obj:`int`
        Number of volumes in the reduced series.
    """
    import pandas as pd

    volume_types = pd.read_table(aslcontext)['volume_type']
    return int(volume_types.isin(get_kept_volume_types(processing_target)).sum())


#: Typical labeling efficiency of each labeling type, used when the metadata do not give one.
#: PCASL and PASL values are from Alsop et al. (2015); the CASL value is from Wang et al. (2005).
TYPICAL_LABELING_EFFICIENCY = {
    'CASL': 0.68,
    'PCASL': 0.85,
    'PASL': 0.98,
}

#: Fraction of the label left by each background suppression pulse (Alsop et al., 2015).
BACKGROUND_SUPPRESSION_PULSE_EFFICIENCY = 0.95


def background_suppression_pulses(metadata):
    """Determine how many background suppression pulses reduced the label.

    Parameters
    ----------
    metadata : :obj:`dict`
        Dictionary of metadata from the ASL file.

    Returns
    -------
    n_pulses : :obj:`int`
        The number of pulses. Zero if background suppression was not applied.
    source : {"metadata", "inferred", "default"}
        Where the number came from:
        ``BackgroundSuppressionNumberPulses`` (or ``BackgroundSuppression`` being false or
        absent), the number of ``BackgroundSuppressionPulseTime`` values,
        or an assumed single pulse when neither field is defined.
    """
    if not metadata.get('BackgroundSuppression', False):
        return 0, 'metadata'

    if 'BackgroundSuppressionNumberPulses' in metadata:
        n_pulses = metadata['BackgroundSuppressionNumberPulses']
        if isinstance(n_pulses, bool) or n_pulses != int(n_pulses) or n_pulses < 0:
            raise ValueError(
                'BackgroundSuppressionNumberPulses must be a non-negative integer, '
                f'but it is {n_pulses}.'
            )
        return int(n_pulses), 'metadata'

    pulse_times = metadata.get('BackgroundSuppressionPulseTime')
    if pulse_times is not None and np.size(pulse_times) > 0:
        return int(np.size(pulse_times)), 'inferred'

    return 1, 'default'


def estimate_labeling_efficiency(metadata, warn=True):
    """Estimate labeling efficiency based on the available metadata.

    Parameters
    ----------
    metadata : :obj:`dict`
        Dictionary of metadata from the ASL file.
    warn : :obj:`bool`, optional
        Whether to warn when the number of background suppression pulses must be assumed.
        Default is True.

    Returns
    -------
    labeleff : :obj:`float`
        Labeling efficiency.

    Notes
    -----
    If LabelingEfficiency is defined in the metadata, then this value will be used as the
    efficiency of the labeling itself.
    Otherwise, that efficiency will be chosen based on the ASL type.
    PCASL and PASL values come from :footcite:t:`alsop_recommended_2015`.
    The CASL value comes from :footcite:t:`wang2005amplitude`.

    If background suppression was applied, the labeling efficiency is then multiplied by 0.95
    for each background suppression pulse, whether or not LabelingEfficiency is defined,
    since each inversion pulse reduces the ASL signal by approximately 5%
    :footcite:p:`alsop_recommended_2015`.
    The BIDS LabelingEfficiency field describes the labeling alone, so it does not include
    this loss.
    If BackgroundSuppressionNumberPulses is not defined, the number of pulses is taken from
    BackgroundSuppressionPulseTime, or one pulse is assumed if that is not defined either
    (see :func:`background_suppression_pulses`).

    References
    ----------
    .. footbibliography::
    """
    if 'LabelingEfficiency' in metadata.keys():
        labeleff = metadata['LabelingEfficiency']
    else:
        labeleff = TYPICAL_LABELING_EFFICIENCY[metadata['ArterialSpinLabelingType']]

    n_pulses, source = background_suppression_pulses(metadata)
    if warn and source == 'default':
        from aslprep import config

        config.loggers.utils.warning(
            'BackgroundSuppression is true, but BackgroundSuppressionNumberPulses is not '
            'defined. Assuming one background suppression pulse to estimate the labeling '
            'efficiency.'
        )

    return labeleff * BACKGROUND_SUPPRESSION_PULSE_EFFICIENCY**n_pulses


def get_inflow_times(metadata, is_casl):
    """Determine the appropriate inflow times for BASIL.

    For PASL data, the inflow time is the inversion time (TI),
    which BIDS stores in the ``PostLabelingDelay`` field.
    For (P)CASL data, the inflow time is the post-labeling delay plus the labeling duration.

    Parameters
    ----------
    metadata : :obj:`dict`
        Dictionary of metadata associated with ASL file.
    is_casl : :obj:`bool`
        True if the data are (P)CASL. False if the data are PASL.

    Returns
    -------
    :obj:`numpy.ndarray`
        1D array of PostLabelingDelay values.
    """
    import numpy as np

    if is_casl:
        return np.add(metadata['PostLabelingDelay'], metadata['LabelingDuration']).tolist()
    else:
        return np.array(metadata['PostLabelingDelay']).tolist()


def get_bolus_duration(metadata, is_casl):
    """Determine the appropriate bolus duration for BASIL.

    For PASL data, the bolus cutoff delay is the first BolusCutOffDelayTime.
    For (P)CASL data, it is the labeling duration.

    Parameters
    ----------
    metadata : :obj:`dict`
        Dictionary of metadata associated with ASL file.
    is_casl : :obj:`bool`
        True if the data are (P)CASL. False if the data are PASL.

    Returns
    -------
    bolus : :obj:`float`
        The bolus value.
    """
    if is_casl:
        return metadata['LabelingDuration']
    elif not metadata['BolusCutOffFlag']:
        raise ValueError('PASL without a bolus cutoff technique is not supported.')
    elif metadata['BolusCutOffTechnique'] == 'Q2TIPS':
        # BolusCutOffDelayTime is a list, and the first entry should be used.
        return metadata['BolusCutOffDelayTime'][0]
    else:  # QUIPSS or QUIPSSII
        return metadata['BolusCutOffDelayTime']


def reduce_metadata_lists(metadata, n_volumes, keep_idx):
    """Reduce any volume-wise metadata fields to only contain values for selected volumes."""
    # A hardcoded list of fields that may have one value for each volume.
    VOLUME_WISE_FIELDS = [
        'PostLabelingDelay',
        'VascularCrushingVENC',
        'LabelingDuration',
        'EchoTime',
        'FlipAngle',
        'RepetitionTimePreparation',
    ]

    for field in VOLUME_WISE_FIELDS:
        if field not in metadata:
            continue

        value = metadata[field]
        if isinstance(value, list):
            if len(value) != n_volumes:
                raise ValueError(
                    f'Number of elements in list-type metadata field {field} ({len(value)}) '
                    f"doesn't equal the number of volumes in the ASL file ({n_volumes})."
                )
            # Reduce to only the selected volumes
            metadata[field] = [value[i] for i in keep_idx]

    return metadata


def check_m0tr(m0tr):
    """Raise an error if the repetition time of the M0 volumes is not positive.

    A repetition time of zero leaves no longitudinal recovery,
    so the M0 scaling in :footcite:t:`alsop_recommended_2015` divides by zero and
    every voxel's CBF becomes zero.

    Parameters
    ----------
    m0tr : float, array-like, or None
        The TR of the M0 volumes, in seconds. None means that no M0 TR is needed.

    Raises
    ------
    ValueError
        If any M0 TR is zero or negative.
    """
    if m0tr is not None and np.any(np.asarray(m0tr, dtype=float) <= 0):
        raise ValueError(
            f'RepetitionTimePreparation for the M0 volumes is {m0tr}, but it must be positive. '
            'ASLPrep uses this value to scale the M0 data, so CBF cannot be calculated. '
            'Please set RepetitionTimePreparation in the sidecar to the time, in seconds, '
            'between successive excitations of the same slice.'
        )


def infer_m0tr(
    *,
    aslcontext,
    metadata,
    m0scan_metadata,
):
    """Infer the repetition time of the M0 volumes based on metadata.

    Parameters
    ----------
    aslcontext : str
        Path to aslcontext file.
    metadata : dict
        Metadata for ASL file.
    m0scan_metadata : dict or None
        Metadata for M0 file, if one exists. Otherwise None.

    Returns
    -------
    m0tr : float or None
        The TR of the M0 scan, if available.
    """
    import pandas as pd

    if metadata['M0Type'] == 'Separate':
        m0tr = m0scan_metadata['RepetitionTimePreparation']
        if np.array(m0tr).size > 1 and np.std(m0tr) > 0:
            raise ValueError('M0 scans have variable TR. ASLPrep does not support this.')

    elif metadata['M0Type'] == 'Included':
        aslcontext = pd.read_table(aslcontext)
        vol_types = aslcontext['volume_type'].tolist()
        m0_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == 'm0scan']
        if np.array(metadata['RepetitionTimePreparation']).size > 1:
            m0tr = np.array(metadata['RepetitionTimePreparation'])[m0_volume_idx]
        else:
            m0tr = metadata['RepetitionTimePreparation']

        if np.array(m0tr).size > 1 and np.std(m0tr) > 0:
            raise ValueError('M0 scans have variable TR. ASLPrep does not support this.')

    elif metadata['M0Type'] == 'Estimate':
        m0tr = None

    elif metadata['M0Type'] == 'Absent':
        aslcontext = pd.read_table(aslcontext)
        vol_types = aslcontext['volume_type'].tolist()
        control_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == 'control']
        cbf_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == 'cbf']
        if control_volume_idx and not cbf_volume_idx:
            # BackgroundSuppression is required, so no need to use get().
            if metadata['BackgroundSuppression']:
                raise ValueError(
                    'Background-suppressed control volumes cannot be used for calibration.'
                )

        if control_volume_idx:
            # Use the control volumes' TR as the M0 TR.
            if np.array(metadata['RepetitionTimePreparation']).size > 1:
                m0tr = np.array(metadata['RepetitionTimePreparation'])[control_volume_idx[0]]
            else:
                m0tr = metadata['RepetitionTimePreparation']

        elif cbf_volume_idx:
            m0tr = None

        else:
            raise RuntimeError(
                'm0scan is absent, '
                'and there are no control volumes that can be used as a substitute'
            )

    else:
        raise RuntimeError('no pathway to m0scan')

    check_m0tr(m0tr)

    return m0tr


def prepare_basil_kwargs(metadata):
    """Prepare keyword arguments for BASIL based on slice timing and multiband metadata.

    Parameters
    ----------
    metadata : dict
        Dictionary of metadata from the ASL file.

    Returns
    -------
    basil_kwargs : dict
        Dictionary of keyword arguments to pass to BASILCBF interface.

    Notes
    -----
    This function handles slice timing correction for both single-band and multiband acquisitions.
    For multiband data, it checks that slice times within each band are monotonic and ascending
    before setting the sliceband and slice_spacing parameters.
    """
    from aslprep import config

    basil_kwargs = {}

    if 'SliceTiming' not in metadata:
        return basil_kwargs

    mb_factor = metadata.get('MultibandAccelerationFactor', 1)
    slice_times = metadata['SliceTiming']

    if mb_factor > 1:
        # Multiband acquisition
        n_slices = len(slice_times)
        n_slices_in_band = n_slices // mb_factor
        ascending_slicetimes = True
        monotonic_slicetimes = True

        for i_band in range(mb_factor):
            band_start = n_slices_in_band * i_band
            band_end = band_start + n_slices_in_band
            band_slice_times = slice_times[band_start:band_end]
            # Round to handle floating point precision issues
            slicetime_diffs = np.unique(np.round(np.diff(band_slice_times), 10))

            # Check if slice times are monotonic for this band
            if slicetime_diffs.size != 1:
                monotonic_slicetimes = False

            # Check if slice times are ascending for this band
            if not np.all(slicetime_diffs > 0):
                ascending_slicetimes = False

        if ascending_slicetimes and monotonic_slicetimes:
            basil_kwargs['sliceband'] = n_slices_in_band
            basil_kwargs['slice_spacing'] = slicetime_diffs[0]
        else:
            config.loggers.utils.warning(
                'Slice times are not ascending. They will be ignored in the BASIL call.'
            )
    else:
        # Single-band acquisition
        # Round to handle floating point precision issues
        slicetime_diffs = np.unique(np.round(np.diff(slice_times), 10))
        # Check if slice times are monotonic
        monotonic_slicetimes = slicetime_diffs.size == 1
        # Check if slice times are ascending
        ascending_slicetimes = np.all(slicetime_diffs > 0)

        if monotonic_slicetimes and ascending_slicetimes:
            basil_kwargs['slice_spacing'] = slicetime_diffs[0]
        else:
            config.loggers.utils.warning(
                'Slice times are not ascending. They will be ignored in the BASIL call.'
            )

    return basil_kwargs


#: Names of the ASL labeling types, as written in the methods boilerplate.
ASL_TYPE_NAMES = {
    'CASL': 'continuous ASL (CASL)',
    'PCASL': 'pseudo-continuous ASL (PCASL)',
    'PASL': 'pulsed ASL (PASL)',
}

#: Names and citations of the bolus cut-off techniques, keyed by BIDS ``BolusCutOffTechnique``.
BOLUS_CUTOFF_TECHNIQUES = {
    'QUIPSS': ('QUIPSS', '@wong1998quantitative'),
    'QUIPSSII': ('QUIPSS II', '@wong1998quantitative'),
    'Q2TIPS': ('Q2TIPS', '@luh1999quipss'),
}

# Volume types whose delay times and labeling durations describe the ASL acquisition.
LABELED_VOLUME_TYPES = ('control', 'label', 'deltam', 'cbf')


def _join_words(words):
    """Join a list of strings with commas and a final "and"."""
    if len(words) <= 2:
        return ' and '.join(words)
    return f'{", ".join(words[:-1])}, and {words[-1]}'


def _plural(count, noun):
    """Write a count and a noun, adding an "s" to the noun unless the count is one."""
    return f'{count} {noun}' + ('' if count == 1 else 's')


def _labeled_values(values, volume_types):
    """Get the sorted, unique values of a metadata field, ignoring the values for M0 volumes."""
    values = np.atleast_1d(np.asarray(values, dtype=float))
    if values.size > 1 and values.size == len(volume_types):
        keep = np.isin(volume_types, LABELED_VOLUME_TYPES)
        if keep.any():
            values = values[keep]
    return np.unique(np.round(values, 6)).tolist()


def _describe_timing(name, values):
    """Describe one or more timing values, in seconds."""
    if len(values) == 1:
        return f'the {name} was {values[0]:g} s'
    return f'the {name}s were {_join_words([f"{value:g}" for value in values])} s'


def describe_asl_acquisition(metadata, aslcontext):
    """Describe the ASL acquisition for the methods boilerplate.

    The description covers the acquisition parameters that the OSIPI ASL lexicon and reporting
    recommendations :footcite:p:`suzuki2024asl` ask authors to report,
    for those that are available in the BIDS metadata.
    Parameters that are not in the metadata are left out.

    Parameters
    ----------
    metadata : :obj:`dict`
        BIDS metadata for the ASL file, before any volumes were removed.
    aslcontext : :obj:`str` or :obj:`pandas.DataFrame`
        The aslcontext file for the ASL file, or its contents.

    Returns
    -------
    :obj:`str`
        The description, with Markdown citations.

    References
    ----------
    .. footbibliography::
    """
    import pandas as pd

    if not isinstance(aslcontext, pd.DataFrame):
        aslcontext = pd.read_table(aslcontext)
    volume_types = aslcontext['volume_type'].to_numpy()

    asl_type = metadata.get('ArterialSpinLabelingType')
    is_pasl = asl_type == 'PASL'
    sentences = []

    # Labeling type and background suppression
    sequence = ASL_TYPE_NAMES.get(asl_type, f'{asl_type} ASL')
    if asl_type == 'PCASL' and metadata.get('PCASLType'):
        sequence = f'{metadata["PCASLType"]} {sequence}'
    elif is_pasl and metadata.get('PASLType'):
        sequence = f'{sequence[:-1]}, {metadata["PASLType"]})'
    article = 'an' if sequence[0] in 'aeiou' else 'a'
    sentence = f'The ASL data were acquired with {article} {sequence} sequence'
    if asl_type == 'CASL' and metadata.get('CASLType') == 'double-coil':
        sentence += ' with a separate labeling coil'
    if 'BackgroundSuppression' in metadata:
        if metadata['BackgroundSuppression']:
            sentence += ', with background suppression'
            n_pulses, pulse_source = background_suppression_pulses(metadata)
            if pulse_source != 'default':
                sentence += f' ({_plural(n_pulses, "pulse")})'
        else:
            sentence += ', without background suppression'
    sentences.append(f'{sentence}.')

    pulse_times = metadata.get('BackgroundSuppressionPulseTime')
    if metadata.get('BackgroundSuppression') and pulse_times:
        pulse_times = _join_words([f'{time:g}' for time in np.atleast_1d(pulse_times)])
        sentences.append(
            f'The background suppression pulses were applied {pulse_times} s after the start of '
            'labeling.'
        )

    # Labeling duration and delay times
    timing = []
    if not is_pasl and 'LabelingDuration' in metadata:
        timing.append(
            _describe_timing(
                'labeling duration',
                _labeled_values(metadata['LabelingDuration'], volume_types),
            )
        )
    if 'PostLabelingDelay' in metadata:
        timing.append(
            _describe_timing(
                'inversion time' if is_pasl else 'post-labeling delay',
                _labeled_values(metadata['PostLabelingDelay'], volume_types),
            )
        )
    if timing:
        timing = ' and '.join(timing)
        sentences.append(f'{timing[0].upper()}{timing[1:]}.')

    # Bolus cut-off
    if is_pasl and 'BolusCutOffFlag' in metadata:
        technique = metadata.get('BolusCutOffTechnique')
        if not metadata['BolusCutOffFlag']:
            sentences.append('No bolus cut-off technique was used.')
        elif technique is None:
            sentences.append('A bolus cut-off technique was used.')
        else:
            name, citation = BOLUS_CUTOFF_TECHNIQUES.get(technique, (technique, None))
            sentence = f'The {name} bolus cut-off technique was used'
            if citation:
                sentence += f' [{citation}]'
            delay_times = np.atleast_1d(metadata.get('BolusCutOffDelayTime', [])).tolist()
            if delay_times:
                sentence += f', with a bolus cut-off delay time (TI1) of {delay_times[0]:g} s'
                if technique == 'Q2TIPS' and len(delay_times) > 1:
                    sentence += f' and a last saturation pulse (TI1s) at {delay_times[-1]:g} s'
            sentences.append(f'{sentence}.')

    # Number of label-control pairs and scanner-calculated volumes
    if 'TotalAcquiredPairs' in metadata:
        n_pairs = metadata['TotalAcquiredPairs']
        verb = 'was' if n_pairs == 1 else 'were'
        sentences.append(f'A total of {_plural(n_pairs, "control-label pair")} {verb} acquired.')

    contents = []
    if 'TotalAcquiredPairs' not in metadata:
        n_pairs = int(min((volume_types == 'control').sum(), (volume_types == 'label').sum()))
        if n_pairs:
            contents.append(_plural(n_pairs, 'control-label pair'))
    for volume_type, noun in (
        ('deltam', 'delta-M volume'),
        ('cbf', 'CBF volume'),
    ):
        n_volumes = int((volume_types == volume_type).sum())
        if n_volumes:
            contents.append(f'{_plural(n_volumes, noun)} calculated on the scanner')
    if contents:
        sentences.append(f'The ASL file contained {_join_words(contents)}.')

    # Voxel size
    if 'AcquisitionVoxelSize' in metadata:
        voxel_size = ' x '.join(f'{size:g}' for size in metadata['AcquisitionVoxelSize'])
        sentences.append(f'The acquisition voxel size was {voxel_size} mm.')

    # Vascular crushing
    if 'VascularCrushing' in metadata:
        if metadata['VascularCrushing']:
            sentence = 'Vascular crushing was applied'
            venc = np.atleast_1d(metadata.get('VascularCrushingVENC', []))
            venc = [f'{value:g}' for value in np.unique(venc) if value]
            if venc:
                sentence += f' (VENC = {_join_words(venc)} cm/s)'
            sentences.append(f'{sentence}.')
        else:
            sentences.append('No vascular crushing was applied.')

    # Labeling location
    distance = metadata.get('LabelingDistance')
    if distance is not None:
        target = 'leading edge of the labeling slab' if is_pasl else 'center of the labeling plane'
        sentences.append(
            f'The {target} was {abs(distance):g} mm from the center of the imaging slab.'
        )
    if is_pasl and metadata.get('LabelingSlabThickness'):
        sentences.append(f'The labeling slab was {metadata["LabelingSlabThickness"]:g} mm thick.')
    if metadata.get('LabelingLocationDescription'):
        description = metadata['LabelingLocationDescription'].strip().rstrip('.')
        sentences.append(f'The labeling location was described as follows: {description}.')

    # (P)CASL labeling pulses
    pulse_fields = (
        ('LabelingPulseAverageGradient', 'an average gradient', 'mT/m'),
        ('LabelingPulseMaximumGradient', 'a maximum gradient', 'mT/m'),
        ('LabelingPulseAverageB1', 'an average B1', 'uT'),
        ('LabelingPulseFlipAngle', 'a flip angle', 'degrees'),
        ('LabelingPulseInterval', 'an interval', 'ms'),
        ('LabelingPulseDuration', 'a duration', 'ms'),
    )
    pulses = [
        f'{label} of {metadata[field]:g} {unit}'
        for field, label, unit in pulse_fields
        if not is_pasl and field in metadata
    ]
    if pulses:
        sentences.append(f'The labeling pulses had {_join_words(pulses)}.')

    return ' '.join(sentences)


def describe_labeling_efficiency(metadata):
    """Describe the labeling efficiency used to calculate CBF, for the methods boilerplate.

    This describes the value from :func:`estimate_labeling_efficiency`.

    Parameters
    ----------
    metadata : :obj:`dict`
        Dictionary of metadata from the ASL file.

    Returns
    -------
    :obj:`str`
        The description, with Markdown citations.
    """
    labeleff = estimate_labeling_efficiency(metadata, warn=False)
    asl_type = metadata['ArterialSpinLabelingType']
    if 'LabelingEfficiency' in metadata:
        source = 'the labeling efficiency in the BIDS metadata'
    else:
        citation = '@wang2005amplitude' if asl_type == 'CASL' else '@alsop_recommended_2015'
        source = f'the typical labeling efficiency for {asl_type} [{citation}]'

    text = f'a labeling efficiency of {labeleff:.3g} ({source}'
    n_pulses, pulse_source = background_suppression_pulses(metadata)
    if n_pulses:
        if pulse_source == 'default':
            text += (
                ', multiplied by 0.95 for the loss of label from one background suppression '
                'pulse, which was assumed because the number of pulses was not specified '
                '[@alsop_recommended_2015]'
            )
        elif pulse_source == 'inferred':
            text += (
                ', multiplied by 0.95 for the loss of label from each of the '
                f'{_plural(n_pulses, "background suppression pulse")}, counted from the '
                'background suppression pulse times [@alsop_recommended_2015]'
            )
        else:
            text += (
                ', multiplied by 0.95 for the loss of label from each of the '
                f'{_plural(n_pulses, "background suppression pulse")} [@alsop_recommended_2015]'
            )
    return f'{text})'
