"""Use an arterial transit time map from a multi-delay run in single-delay tissue-T1 CBF."""

from __future__ import annotations

import numpy as np

#: Status codes of the external ATT on the single-delay grid (``desc-external_attstatus``).
EXTERNAL_ATT_STATUS = {
    0: (
        "no ATT: outside the multi-delay run's coverage or the brain mask, or where the "
        'multi-delay fit failed'
    ),
    1: 'measured: located by the multi-delay fit',
    2: 'filled with the tissue median, clipped into the interval of the multi-delay fit',
    3: 'filled with the tissue median (no interval from the multi-delay fit)',
    4: 'not located, and not filled (--att-fill none)',
    5: 'not located, with too little gray and white matter to fill',
    6: "no tissue signal at this run's delay with this ATT, so no CBF",
    7: (
        'not located, and not filled because no gray or no white matter median could be '
        'measured from the located voxels'
    ),
}

#: Source run status values (:data:`aslprep.utils.cbf.ATT_STATUS`) with an ATT interval.
_WITH_INTERVAL = (1, 2, 4)

#: Fraction of a voxel's interpolation weight that must come from located voxels for its ATT
#: to count as measured.
MIN_LOCATED_WEIGHT = 0.5

#: Gray plus white matter partial volume below which an unlocated voxel is not filled.
MIN_FILL_TISSUE = 0.01

#: A tissue basis below this fraction of the voxel's all-arrived basis (ATT of zero) has too
#: little signal to quantify CBF: it would amplify the noise more than tenfold.
MIN_RELATIVE_BASIS = 0.1


def check_external_att_metadata(att_metadata, source_metadata, asl_metadata, att_name, asl_file):
    """Check that an ATT map can be used for another run's tissue-T1 CBF.

    Parameters
    ----------
    att_metadata : :obj:`dict`
        The ATT map's sidecar, with the ``Model`` ASLPrep records.
    source_metadata : :obj:`dict` or None
        The multi-delay run's raw ASL metadata, or None if it was not found.
    asl_metadata : :obj:`dict`
        The single-delay run's metadata.
    att_name, asl_file : :obj:`str`
        The ATT map and the ASL run, for the messages.

    Returns
    -------
    :obj:`list` of :obj:`str`
        Warnings: labeling that differs between the runs, which changes ATT, or that could not
        be compared.

    Raises
    ------
    ValueError
        If the map was not estimated with the tissue-T1 model, whose ATT is the one this model
        needs (the consensus model biases white matter ATT low).
    """
    from aslprep.utils.calibration import TISSUE_T1_PARAMETERS

    parameters = (att_metadata.get('Model') or {}).get('Parameters') or {}
    if not all(name in parameters for name in TISSUE_T1_PARAMETERS):
        raise ValueError(
            f'The ATT map {att_name} was not estimated with the tissue-T1 kinetic model '
            '(its sidecar has no tissue T1s in Model.Parameters). The consensus model, in which '
            'the label decays with the T1 of blood throughout, biases ATT, most of all in white '
            'matter, so its map would bias the single-delay CBF of '
            f'{asl_file}. Process the multi-delay run with --kinetic-model tissue-t1, or use '
            '--ignore att.'
        )

    if source_metadata is None:
        warning = (
            f'The multi-delay run of the ATT map {att_name} is not in the raw dataset, so its '
            f'labeling could not be compared with that of {asl_file}. ATT depends on where the '
            'label is applied, so check that the runs label in the same way.'
        )
        return [warning]

    warnings = []
    for field in ('ArterialSpinLabelingType', 'LabelingDistance'):
        source = source_metadata.get(field)
        target = asl_metadata.get(field)
        if source is not None and target is not None and source != target:
            warnings.append(
                f'The ATT map {att_name} comes from a run with {field} {source}, but '
                f'{asl_file} has {target}. ATT depends on where the label is applied, so it may '
                'not be the same in both runs.'
            )
    return warnings


def _affine_matrix(xfm_file):
    """Load a single affine ITK transform as the RAS point mapping it applies when resampling."""
    import nitransforms as nt

    matrix = np.asarray(nt.linear.load(xfm_file).matrix)
    if matrix.ndim == 3:
        if matrix.shape[0] != 1:
            raise ValueError(f'Expected a single affine transform in {xfm_file}.')
        matrix = matrix[0]
    return matrix


def target_to_source_voxels(target_affine, target_shape, source_affine, target_xfm, source_xfm):
    """Map the voxel centers of one run's ASL grid to voxel coordinates of another's.

    Each run's aslref-to-anatomical transform (as ASLPrep writes it) maps anatomical points to
    points of that run's ASL reference, so a target point goes to the anatomical image through
    the inverse of the target run's transform, and from there to the source run's reference
    through the source run's transform.

    Parameters
    ----------
    target_affine : (4, 4) :obj:`numpy.ndarray`
        The voxel-to-RAS affine of the target (single-delay) grid.
    target_shape : :obj:`tuple` of :obj:`int`
        Its shape.
    source_affine : (4, 4) :obj:`numpy.ndarray`
        The voxel-to-RAS affine of the source (multi-delay) grid.
    target_xfm, source_xfm : :obj:`str`
        The runs' aslref-to-anatomical transforms (ITK text files).

    Returns
    -------
    :obj:`numpy.ndarray` of shape (3, n_voxels)
        Source voxel coordinates of every target voxel, in C order.
    """
    target_to_anat = np.linalg.inv(_affine_matrix(target_xfm))
    anat_to_source = _affine_matrix(source_xfm)
    vox2vox = np.linalg.inv(source_affine) @ anat_to_source @ target_to_anat @ target_affine
    grid = np.indices(target_shape[:3]).reshape(3, -1)
    points = np.vstack([grid, np.ones(grid.shape[1])])
    return (vox2vox @ points)[:3]


def normalized_sample(values, weights, coords):
    """Interpolate ``values`` linearly at ``coords``, using only voxels where ``weights`` is 1.

    Returns the weighted values divided by the interpolated weights, and the interpolated
    weights, so that voxels without a value do not blend into their neighbors' values.
    """
    from scipy.ndimage import map_coordinates

    weights = np.asarray(weights, dtype=np.float64)
    values = np.where(weights > 0, np.nan_to_num(np.asarray(values, dtype=np.float64)), 0.0)
    kwargs = {'order': 1, 'mode': 'constant', 'cval': 0.0, 'prefilter': False}
    weight = map_coordinates(weights, coords, **kwargs)
    weighted = map_coordinates(values * weights, coords, **kwargs)
    with np.errstate(divide='ignore', invalid='ignore'):
        sampled = np.where(weight > 0, weighted / weight, np.nan)
    return sampled, weight


def weighted_median(values, weights):
    """The median of ``values`` weighted by ``weights``, or NaN without positive weights."""
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    keep = np.isfinite(values) & (weights > 0)
    if not keep.any():
        return np.nan

    order = np.argsort(values[keep])
    cumulative = np.cumsum(weights[keep][order])
    index = np.searchsorted(cumulative, cumulative[-1] / 2)
    return float(values[keep][order][index])


def resample_external_att(
    att,
    status,
    ci_lower,
    ci_upper,
    coords,
    gm_pv,
    wm_pv,
    mask,
    nominal_cbf,
    fill='median',
):
    """Bring an ATT map onto another run's grid, and fill the voxels where it was not located.

    Parameters
    ----------
    att, status, ci_lower, ci_upper : :obj:`numpy.ndarray`
        The source run's ATT, its :data:`~aslprep.utils.cbf.ATT_STATUS`, and its noise-aware
        interval, on the source grid.
    coords : :obj:`numpy.ndarray` of shape (3, n_voxels)
        The source voxel coordinates of the target voxels
        (see :func:`target_to_source_voxels`).
    gm_pv, wm_pv : :obj:`numpy.ndarray` of shape (n_voxels,)
        The target run's gray and white matter partial volumes.
    mask : :obj:`numpy.ndarray` of shape (n_voxels,)
        The target run's brain mask.
    nominal_cbf : :obj:`tuple` of :obj:`float`
        The gray and white matter nominal CBFs of the tissue-T1 model, which weight the tissue
        medians in mixed voxels as the model weights the tissues.
    fill : {"median", "none"}
        How to fill target voxels where the source run did not locate ATT.

    Returns
    -------
    target_att : :obj:`numpy.ndarray` of shape (n_voxels,)
        ATT on the target grid, NaN where there is none.
    target_status : :obj:`numpy.ndarray` of shape (n_voxels,)
        :data:`EXTERNAL_ATT_STATUS` codes (all but 6, which quantification sets where the ATT
        gives no signal).
    summary : :obj:`dict`
        The tissue medians and the fraction of brain voxels with each status.
    """
    from scipy.ndimage import map_coordinates

    from aslprep.utils.cbf import ATT_STATUS, tissue_mixing_weights

    if fill not in ('median', 'none'):
        raise ValueError(f'Unknown ATT fill {fill!r}; expected "median" or "none".')

    status = np.rint(np.asarray(status)).astype(int)
    located = (status == ATT_STATUS['located']).astype(np.float64)
    with_interval = np.isin(status, _WITH_INTERVAL).astype(np.float64)
    sampled_att, located_weight = normalized_sample(att, located, coords)
    lower, _ = normalized_sample(ci_lower, with_interval, coords)
    upper, interval_weight = normalized_sample(ci_upper, with_interval, coords)
    # Nearest neighbor, so that the categories are not mixed; outside the source grid, zero
    nearest = map_coordinates(status, coords, order=0, mode='constant', cval=0)
    covered = map_coordinates(
        np.ones(status.shape),
        coords,
        order=1,
        mode='constant',
        cval=0.0,
        prefilter=False,
    )

    mask = np.asarray(mask, dtype=bool)
    gm_pv = np.nan_to_num(np.asarray(gm_pv, dtype=np.float64))
    wm_pv = np.nan_to_num(np.asarray(wm_pv, dtype=np.float64))
    inside = mask & (covered >= 0.5) & (nearest != ATT_STATUS['not_fitted'])
    measured = (
        inside
        & (located_weight >= MIN_LOCATED_WEIGHT)
        & np.isfinite(sampled_att)
        & (sampled_att >= 0)
    )
    target_att = np.where(measured, sampled_att, np.nan)
    target_status = np.zeros(gm_pv.shape, dtype=np.uint8)
    target_status[measured] = 1

    medians = {
        'gray': weighted_median(sampled_att[measured], gm_pv[measured]),
        'white': weighted_median(sampled_att[measured], wm_pv[measured]),
    }
    unlocated = inside & ~measured
    tissue = (gm_pv + wm_pv) >= MIN_FILL_TISSUE
    target_status[unlocated & ~tissue] = 5
    to_fill = unlocated & tissue
    if fill == 'none':
        target_status[to_fill] = 4
    elif not np.all(np.isfinite(list(medians.values()))):
        target_status[to_fill] = 7
    else:
        weights = tissue_mixing_weights(gm_pv[to_fill], wm_pv[to_fill], *nominal_cbf)
        filled = weights[:, 0] * medians['gray'] + weights[:, 1] * medians['white']
        # A fill never contradicts what the multi-delay data showed about the voxel
        has_interval = (
            (interval_weight[to_fill] >= MIN_LOCATED_WEIGHT)
            & np.isfinite(lower[to_fill])
            & np.isfinite(upper[to_fill])
        )
        clipped = np.where(
            has_interval,
            np.clip(filled, lower[to_fill], np.maximum(upper[to_fill], lower[to_fill])),
            filled,
        )
        target_att[to_fill] = clipped
        target_status[to_fill] = np.where(has_interval, 2, 3)

    n_mask = max(int(mask.sum()), 1)
    summary = {
        # None (null in the sidecars, as JSON has no NaN) where no median could be measured
        'GrayMatterMedianATT': medians['gray'] if np.isfinite(medians['gray']) else None,
        'WhiteMatterMedianATT': medians['white'] if np.isfinite(medians['white']) else None,
        # Before quantification, which marks the voxels without tissue signal (status 6)
        'ResamplingStatusFractions': {
            str(code): float((target_status[mask] == code).sum() / n_mask)
            for code in EXTERNAL_ATT_STATUS
            if code != 6
        },
    }
    return target_att, target_status, summary
