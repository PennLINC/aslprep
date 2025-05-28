"""Functions for working with ASL data."""

from __future__ import annotations

import numpy as np


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
    If LabelingEfficiency is defined in the metadata, then this value will be used.
    Otherwise, efficiency will be estimated based on the ASL type and number of background
    suppression pulses (if any).
    PCASL and PASL values come from :footcite:t:`alsop_recommended_2015`.
    The CASL value comes from :footcite:t:`wang2005amplitude`.
    The adjustment based on number of background suppression pulses is not described in any papers,
    as far as we know, but is apparently widely used.

    References
    ----------
    .. footbibliography::
    """
    if 'LabelingEfficiency' in metadata.keys():
        labeleff = metadata['LabelingEfficiency']
    else:
        BASE_LABELEFF = {
            'CASL': 0.68,
            'PCASL': 0.85,
            'PASL': 0.98,
        }
        labeleff = BASE_LABELEFF[metadata['ArterialSpinLabelingType']]

        if metadata.get('BackgroundSuppression', False):
            BS_PULSE_EFF = 0.95  # hardcoded BackgroundSuppressionPulse efficiency
            # We assume there was one pulse if suppression was applied,
            # but the number of pulses isn't defined.
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
            basil_kwargs['sliceband'] = mb_factor
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
            config.loggers.interface.warning(
                'Slice times are not ascending. They will be ignored in the BASIL call.'
            )

    return basil_kwargs
