"""Functions for working with ASL data."""

from __future__ import annotations

from typing import Any

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
        True if the data are multi-delay/TI. Fale if not.
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


def get_inflow_times(metadata: dict[str, Any], is_casl: bool) -> list:
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


def get_bolus_duration(metadata: dict[str, Any], is_casl: bool) -> float:
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
