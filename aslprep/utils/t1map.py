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

#: Default edge, in mm, of the cube in which the local-regression method unmixes R1.
DEFAULT_KERNEL_MM = 7.0

#: Fewest voxels along each axis of that cube.
MIN_KERNEL_VOXELS = 3

#: Least support (sum of probabilities) of each tissue in a local-regression fit, as a fraction
#: of the cube's voxels.
MIN_SUPPORT_FRACTION = 0.05

#: Largest condition number of a local-regression fit's normalized Gram matrix.
MAX_CONDITION = 100.0

#: Largest contamination of a gray or white matter R1 by a tissue left out of the fit.
MAX_CONTAMINATION = 0.05

#: A local gray or white matter R1 further than this factor from the tissue's median R1 is not
#: used, and the median is used instead.
MAX_R1_RATIO = 2.0

#: Local-regression fit status codes of :func:`local_regression_r1`.
#: Each says which tissues were fit; a tissue was left out only when its contamination of the
#: gray and white matter estimates was below the limit.
R1_FIT_STATUS = {
    0: 'outside the brain mask',
    1: 'gray matter, white matter and CSF fit',
    2: 'gray and white matter fit',
    3: 'gray matter and CSF fit',
    4: 'white matter and CSF fit',
    5: 'gray matter only fit',
    6: 'white matter only fit',
    7: 'not identifiable',
}

# Tissue indices (GM, WM, CSF) of each fit, in order of precedence, with its status code.
_R1_SUBSETS = ((1, (0, 1, 2)), (2, (0, 1)), (3, (0, 2)), (4, (1, 2)), (5, (0,)), (6, (1,)))

# Anatomical slices unmixed at once, which bounds the memory used.
_SLAB = 24

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


def kernel_voxels(kernel_mm, affine):
    """Find the odd number of voxels along each axis that best spans ``kernel_mm``.

    At least :data:`MIN_KERNEL_VOXELS` voxels are used along each axis, so that every fit has
    neighbors, even on a grid whose voxels are close to ``kernel_mm`` or larger.

    Parameters
    ----------
    kernel_mm : :obj:`float`
        The cube's edge, in mm.
    affine : (4, 4) :obj:`numpy.ndarray`
        The grid's voxel-to-RAS affine.

    Returns
    -------
    :obj:`tuple` of :obj:`int`
    """
    if not np.isfinite(kernel_mm) or kernel_mm <= 0:
        raise ValueError(f'The kernel size must be finite and positive, got {kernel_mm}.')

    voxel_size = np.sqrt((np.asarray(affine)[:3, :3] ** 2).sum(axis=0))
    # The odd number of voxels nearest the kernel's span along each axis
    nearest = [int(2 * np.round((kernel_mm / v - 1) / 2) + 1) for v in voxel_size]
    return tuple(max(n, MIN_KERNEL_VOXELS) for n in nearest)


def local_regression_r1(
    t1,
    gm_probability,
    wm_probability,
    csf_probability,
    mask,
    kernel_size,
    medians,
):
    """Estimate the R1 (1 / T1) of gray and white matter in each voxel by local regression.

    Under fast exchange, a voxel's R1 is approximately the sum of its tissues' R1s, weighted by
    their fractions. So in a box around each voxel of ``mask``, the R1 of the voxels with a
    finite, positive T1 is regressed on the tissue probability maps, without an intercept, and the
    coefficients are the tissues' R1s there (as :func:`~aslprep.utils.pvc.local_regression_pvc`
    does for perfusion, but with CSF, whose R1 is far from zero, as a regressor).
    The fits are tried in the order of :data:`R1_FIT_STATUS`
    (see :func:`~aslprep.utils.pvc.local_regression_subsets`), with :data:`MIN_SUPPORT_FRACTION`,
    :data:`MAX_CONDITION` and :data:`MAX_CONTAMINATION`.

    A tissue's local R1 is used where that tissue was in the fit and the R1 is within a factor of
    :data:`MAX_R1_RATIO` of the tissue's median R1. Elsewhere, including outside ``mask``, the
    median is used.

    Parameters
    ----------
    t1 : :obj:`numpy.ndarray`
        T1, in seconds, on the grid of the tissue probability maps.
    gm_probability, wm_probability, csf_probability : :obj:`numpy.ndarray`
        Tissue probability maps.
    mask : :obj:`numpy.ndarray` of :obj:`bool`
        The brain mask.
    kernel_size : :obj:`tuple` of :obj:`int`
        The box size along each array axis, in voxels. All values must be odd.
    medians : :obj:`dict`
        The tissues' median T1s, from :func:`tissue_t1_medians`.

    Returns
    -------
    :obj:`dict`
        ``gm_r1`` and ``wm_r1`` (1/s), ``status`` (see :data:`R1_FIT_STATUS`), the boolean maps
        ``gm_estimated`` and ``wm_estimated`` of where the local R1s were used, and
        ``estimated_fraction``, the fraction of each tissue (``'gray'``, ``'white'``), by
        probability within ``mask``, whose R1 was estimated locally.
    """
    from aslprep.utils.pvc import local_regression_subsets

    t1 = np.asarray(t1, dtype=np.float64)
    probabilities = (gm_probability, wm_probability, csf_probability)
    tpms = [np.asarray(p, dtype=np.float64) for p in probabilities]
    mask = np.asarray(mask, dtype=bool)
    kernel_size = tuple(int(k) for k in kernel_size)
    shape = t1.shape
    with np.errstate(divide='ignore', invalid='ignore'):
        r1 = np.where(np.isfinite(t1) & (t1 > 0), 1 / t1, np.nan)

    median_r1 = {tissue: 1 / medians[tissue]['t1'] for tissue in ('gray', 'white')}
    out = {
        'gm_r1': np.full(shape, median_r1['gray']),
        'wm_r1': np.full(shape, median_r1['white']),
        'status': np.zeros(shape, dtype=np.uint8),
        'gm_estimated': np.zeros(shape, dtype=bool),
        'wm_estimated': np.zeros(shape, dtype=bool),
    }
    tau = MIN_SUPPORT_FRACTION * float(np.prod(kernel_size))

    # Unmix the mask's bounding box, in slabs along the first axis, each with a halo of half the
    # kernel, so that every box sum is the same as over the whole image.
    if mask.any():
        lo = [int(idx.min()) for idx in np.nonzero(mask)]
        hi = [int(idx.max()) + 1 for idx in np.nonzero(mask)]
        half = [k // 2 for k in kernel_size]
        for start in range(lo[0], hi[0], _SLAB):
            stop = min(start + _SLAB, hi[0])
            outer = (
                slice(max(start - half[0], 0), min(stop + half[0], shape[0])),
                *(
                    slice(max(lo[ax] - half[ax], 0), min(hi[ax] + half[ax], shape[ax]))
                    for ax in (1, 2)
                ),
            )
            fit = local_regression_subsets(
                r1[outer],
                [tpm[outer] for tpm in tpms],
                mask[outer],
                np.isfinite(r1[outer]),
                subsets=_R1_SUBSETS,
                targets=(0, 1),
                unidentified=7,
                kernel_size=kernel_size,
                tau=tau,
                kappa_max=MAX_CONDITION,
                c_max=MAX_CONTAMINATION,
                residuals=False,
            )
            inner = (
                slice(start - outer[0].start, stop - outer[0].start),
                slice(None),
                slice(None),
            )
            target = (slice(start, stop), outer[1], outer[2])
            status = fit['status'][inner]
            out['status'][target] = status
            for i_tissue, (key, tissue) in enumerate((('gm', 'gray'), ('wm', 'white'))):
                codes = [code for code, subset in _R1_SUBSETS if i_tissue in subset]
                included = np.isin(status, codes)
                estimate = fit['beta'][inner][..., i_tissue]
                estimated = (
                    included
                    & np.isfinite(estimate)
                    & (estimate >= median_r1[tissue] / MAX_R1_RATIO)
                    & (estimate <= median_r1[tissue] * MAX_R1_RATIO)
                )
                out[f'{key}_r1'][target] = np.where(estimated, estimate, median_r1[tissue])
                out[f'{key}_estimated'][target] = estimated

    out['estimated_fraction'] = {}
    for key, tissue, tpm in (('gm', 'gray', tpms[0]), ('wm', 'white', tpms[1])):
        weight = np.where(mask, np.clip(np.nan_to_num(tpm), 0, 1), 0.0)
        total = weight.sum()
        out['estimated_fraction'][tissue] = (
            float((weight * out[f'{key}_estimated']).sum() / total) if total > 0 else 0.0
        )
    return out


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
