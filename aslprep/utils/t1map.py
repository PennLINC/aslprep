"""Measure tissue T1s from a quantitative T1 map, for the tissue-T1 kinetic model."""

from __future__ import annotations

import numpy as np

#: Tissue probability above which a voxel counts as pure gray or white matter.
PURE_TISSUE_PROBABILITY = 0.95

#: CSF probability below which a voxel counts as free of CSF.
MAX_CSF_PROBABILITY = 0.05

#: Fewest pure-tissue voxels from which a tissue's T1 is measured.
MIN_PURE_VOXELS = 100

#: A measured T1 outside this range, as a multiple of the default for the field strength,
#: is reported as implausible (usually wrong units, or a map that is not of T1).
PLAUSIBLE_T1_RATIO = (0.6, 1.5)

#: Brain median above which a T1 map without ``Units`` is read as milliseconds.
MILLISECONDS_THRESHOLD = 100

_UNITS = {
    's': 1.0,
    'second': 1.0,
    'seconds': 1.0,
    'ms': 1e-3,
    'millisecond': 1e-3,
    'milliseconds': 1e-3,
}


def t1_scale(values, metadata):
    """Find the factor that converts a T1 map to seconds, and say how it was decided.

    Parameters
    ----------
    values : :obj:`numpy.ndarray`
        The T1 map's values.
    metadata : :obj:`dict`
        The T1 map's sidecar. ``Units`` (``s`` or ``ms``) decides the scale when present.

    Returns
    -------
    scale : :obj:`float`
        1 for seconds, 0.001 for milliseconds.
    how : :obj:`str`
        How the units were decided.

    Raises
    ------
    ValueError
        If ``Units`` is not a unit of time ASLPrep knows, or the map has no positive values.
    """
    units = metadata.get('Units')
    if units is not None:
        scale = _UNITS.get(str(units).strip().lower())
        if scale is None:
            raise ValueError(f'The T1 map has Units {units!r}; expected seconds or milliseconds.')

        return scale, f'Units {units!r} in its sidecar'

    positive = values[np.isfinite(values) & (values > 0)]
    if not positive.size:
        raise ValueError('The T1 map has no positive values.')

    median = float(np.median(positive))
    if median > MILLISECONDS_THRESHOLD:
        return 1e-3, f'read as milliseconds, as its median ({median:g}) exceeds 100'
    return 1.0, f'read as seconds, as its median ({median:g}) does not exceed 100'


def pure_tissue_mask(tissue_probability, csf_probability):
    """Select voxels of one tissue that partial volume effects leave unmixed.

    Voxels whose tissue probability is at least :data:`PURE_TISSUE_PROBABILITY` and whose CSF
    probability is below :data:`MAX_CSF_PROBABILITY` are selected, then eroded by one voxel,
    because a T1 map's own blurring mixes tissues across boundaries the segmentation does not see.
    """
    from scipy.ndimage import binary_erosion

    pure = (tissue_probability >= PURE_TISSUE_PROBABILITY) & (
        csf_probability < MAX_CSF_PROBABILITY
    )
    return binary_erosion(pure)


def tissue_t1_medians(t1, gm_probability, wm_probability, csf_probability):
    """Measure the T1 of gray and white matter as medians over pure-tissue voxels.

    Parameters
    ----------
    t1 : :obj:`numpy.ndarray`
        T1, in seconds, on the grid of the tissue probability maps.
    gm_probability, wm_probability, csf_probability : :obj:`numpy.ndarray`
        Tissue probability maps.

    Returns
    -------
    :obj:`dict`
        For each of ``'gray'`` and ``'white'``, a dictionary with the median T1 (``t1``),
        its interquartile range (``iqr``) and the number of voxels (``n``).

    Raises
    ------
    ValueError
        If a tissue has fewer than :data:`MIN_PURE_VOXELS` pure voxels with a finite, positive T1.
    """
    results = {}
    for tissue, probability in (('gray', gm_probability), ('white', wm_probability)):
        mask = pure_tissue_mask(probability, csf_probability)
        values = t1[mask]
        values = values[np.isfinite(values) & (values > 0)]
        if values.size < MIN_PURE_VOXELS:
            raise ValueError(
                f'Only {values.size} pure {tissue} matter voxels have a finite, positive T1 '
                f'(at least {MIN_PURE_VOXELS} are needed). Check that the T1 map covers the '
                'brain and is registered to the anatomical image.'
            )

        q1, median, q3 = np.percentile(values, [25, 50, 75])
        results[tissue] = {'t1': float(median), 'iqr': float(q3 - q1), 'n': int(values.size)}
    return results


def check_t1map_field_strength(t1map, field_strength, asl_file):
    """Refuse a T1 map acquired at another field strength than an ASL run, as T1 depends on it.

    Parameters
    ----------
    t1map : :obj:`dict`
        The T1 map, from :func:`~aslprep.utils.bids.collect_t1map`. Without a
        ``MagneticFieldStrength`` in its sidecar, nothing is checked.
    field_strength : :obj:`float`
        The ASL run's field strength, in tesla.
    asl_file : :obj:`str`
        The ASL run, for the error message.

    Raises
    ------
    ValueError
        If the field strengths differ.
    """
    t1map_field = t1map['metadata'].get('MagneticFieldStrength')
    if t1map_field is not None and not np.isclose(t1map_field, field_strength):
        raise ValueError(
            f'The T1 map {t1map["t1map"]} was acquired at {t1map_field} T, but {asl_file} at '
            f'{field_strength} T, and T1 depends on the field strength. '
            'Use --ignore t1map to use the default tissue T1s instead.'
        )


def implausible_t1_warnings(medians, field_strength):
    """Warn about measured T1s far from the defaults for the field strength.

    Parameters
    ----------
    medians : :obj:`dict`
        From :func:`tissue_t1_medians`.
    field_strength : :obj:`float` or None
        Magnetic field strength, in tesla. Without it, or for field strengths without defaults,
        no warnings are made.

    Returns
    -------
    :obj:`list` of :obj:`str`
    """
    from aslprep.utils.cbf import T1TISSUE_BY_FIELD_STRENGTH, T1WM_BY_FIELD_STRENGTH

    defaults = {
        'gray': T1TISSUE_BY_FIELD_STRENGTH.get(field_strength),
        'white': T1WM_BY_FIELD_STRENGTH.get(field_strength),
    }
    low, high = PLAUSIBLE_T1_RATIO
    warnings = []
    for tissue, default in defaults.items():
        if default is None:
            continue

        t1 = medians[tissue]['t1']
        if not low * default <= t1 <= high * default:
            warnings.append(
                f'The {tissue} matter T1 measured from the T1 map ({t1:.3g} s) is far from the '
                f'default for {field_strength:g} T ({default:g} s). Check the units of the map, '
                'and that it is a T1 map registered to the anatomical image.'
            )
    return warnings
