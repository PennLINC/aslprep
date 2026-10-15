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


def estimate_labeling_efficiency(metadata):
    """Estimate labeling efficiency based on the available metadata.

    Parameters
    ----------
    metadata : :obj:`dict`
        Dictionary of metadata from the ASL file.

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
    If BackgroundSuppressionNumberPulses is not defined, one pulse is assumed.

    References
    ----------
    .. footbibliography::
    """
    BASE_LABELEFF = {
        'CASL': 0.68,
        'PCASL': 0.85,
        'PASL': 0.98,
    }
    BS_PULSE_EFF = 0.95  # hardcoded BackgroundSuppressionPulse efficiency

    if 'LabelingEfficiency' in metadata.keys():
        labeleff = metadata['LabelingEfficiency']
    else:
        labeleff = BASE_LABELEFF[metadata['ArterialSpinLabelingType']]

    if metadata.get('BackgroundSuppression', False):
        if 'BackgroundSuppressionNumberPulses' not in metadata.keys():
            from aslprep import config

            config.loggers.utils.warning(
                'BackgroundSuppression is true, but BackgroundSuppressionNumberPulses is not '
                'defined. Assuming one background suppression pulse to estimate the labeling '
                'efficiency.'
            )

        labeleff *= BS_PULSE_EFF ** metadata.get('BackgroundSuppressionNumberPulses', 1)

    return labeleff


def get_inflow_times(metadata, is_casl):
    """Determine the appropriate inflow times for BASIL.

    For PASL data, the inflow time (TI) is just the post-labeling delay (PLD).
    For (P)CASL data, TI is PLD plus the labeling duration.

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
