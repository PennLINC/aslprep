"""Model-based detection of corrupted control-label pairs in multi-delay ASL data.

A corrupted pair (from motion between its control and label images, a slice dropout,
a failed labeling pulse, or a noise spike) can bias the kinetic model fit of every voxel,
because each delay is represented by only a few pairs.
These functions predict each pair from the general kinetic model,
optionally corrected toward the other pairs acquired at the same delay,
score every pair by its residual, and refit without the pairs whose scores are outliers.
"""

import numpy as np

from aslprep.utils.cbf import _gkm_arterial_basis, _gkm_tissue_basis, fit_deltam_multipld

#: Median of a chi-squared variable with one degree of freedom, which converts the median
#: squared residual into a variance.
CHI2_1_MEDIAN = 0.454936

#: Scale factor that makes the median absolute deviation a consistent estimator of the
#: standard deviation of normal data.
MAD_SCALE = 1.4826

#: Fewest pairs for which outliers are detected. The robust z-scores need several pairs.
MIN_PAIRS = 4

#: Fewest other pairs at a delay for the same-delay correction. The median of fewer pairs is
#: not robust: with one or two, a single corrupted pair moves it, and the correction would
#: then spread that pair's artifact to the clean pairs at its delay.
MIN_SAME_DELAY = 3

#: Range, relative to its median over voxels, to which each voxel's noise variance is clipped,
#: so that a few voxels with very small estimated noise do not dominate the scores.
NOISE_CLIP = (0.25, 4.0)


def predict_deltam_multipld(
    cbf,
    att,
    abat,
    abv,
    scaled_m0data,
    plds,
    labeleff,
    t1blood,
    partition_coefficient,
    is_casl,
    tau=None,
    ti1=None,
):
    """Predict delta-M from fitted general kinetic model parameters.

    This is the forward model fit by :func:`~aslprep.utils.cbf.fit_deltam_multipld`,
    with the same arguments.

    Parameters
    ----------
    cbf, att, abat, abv : :obj:`numpy.ndarray` of shape (n_voxels,)
        Fitted parameters.
    scaled_m0data : :obj:`numpy.ndarray` of shape (n_voxels,)
        M0 data, after any scaling requested by the user.
    plds : :obj:`numpy.ndarray` of shape (n_voxels, n_volumes)
        Delay times: post-labeling delays for (P)CASL data, inversion times for PASL data.
    labeleff, t1blood, partition_coefficient : :obj:`float`
        Labeling efficiency, T1 of arterial blood, and blood-brain partition coefficient.
    is_casl : :obj:`bool`
        True for (P)CASL data.
    tau : :obj:`float` or array of shape (n_volumes,) or None, optional
        Labeling duration. Only defined for (P)CASL data.
    ti1 : :obj:`float` or None, optional
        TI1. Only defined for PASL data.

    Returns
    -------
    deltam : :obj:`numpy.ndarray` of shape (n_voxels, n_volumes)
        Predicted delta-M. Voxels with non-finite parameters are NaN.
    """
    m0_a = np.asarray(scaled_m0data, dtype=float) / partition_coefficient
    if is_casl:
        bolus = np.broadcast_to(np.asarray(tau, dtype=float), (plds.shape[1],))[None, :]
    else:
        bolus = float(ti1)

    kwargs = {'is_casl': is_casl, 'labeleff': labeleff, 't1blood': t1blood, 'xp': np}
    tissue = _gkm_tissue_basis(np.asarray(att)[:, None], plds, bolus, **kwargs)
    arterial = _gkm_arterial_basis(np.asarray(abat)[:, None], plds, bolus, **kwargs)
    return m0_a[:, None] * (
        np.asarray(cbf)[:, None] * tissue + np.asarray(abv)[:, None] * arterial
    )


def _delay_groups(plds, tau, is_casl):
    """Assign each volume to a group of volumes with the same delay (and labeling duration).

    Volumes are grouped as :func:`~aslprep.utils.cbf.fit_deltam_multipld` averages them:
    by their delays in every voxel, so slice-shifted delays group by acquisition delay.
    """
    keys = np.vstack([plds, np.asarray(tau, dtype=float)[None, :]]) if is_casl else plds
    _, groups = np.unique(keys, axis=1, return_inverse=True)
    return groups.ravel()


def _median_variance(n):
    """Return the variance of the median of ``n`` standard normal values, approximately."""
    return 1.0 / n if n < 3 else np.pi / (2 * n)


def _robust_z(values):
    """Standardize values by their median and median absolute deviation."""
    center = np.median(values)
    spread = MAD_SCALE * np.median(np.abs(values - center))
    spread = max(spread, np.finfo(float).tiny)
    return (values - center) / spread


def score_pairs(deltam, predicted, groups, keep, voxels=None, shrinkage=2.0, smooth=None):
    """Score each pair by its residual from a prediction.

    Each pair is predicted by the kinetic model, plus a correction toward the other kept pairs
    with the same delay: the voxelwise median of their model residuals, weighted by
    ``m / (m + shrinkage)``, where ``m`` is the number of those pairs.
    A ``shrinkage`` of zero predicts each pair from the others at its delay alone,
    and an infinite one from the model alone.
    The correction absorbs systematic misfit of the model at a delay, which is common to all of
    its pairs. It is only applied when at least :data:`MIN_SAME_DELAY` other kept pairs share
    the delay, because the median of fewer is not robust to a corrupted pair. Otherwise the pair
    is predicted from the model, and so from the pairs at the other delays.

    Parameters
    ----------
    deltam : :obj:`numpy.ndarray` of shape (n_voxels, n_pairs)
        Observed delta-M.
    predicted : :obj:`numpy.ndarray` of shape (n_voxels, n_pairs)
        Delta-M predicted by the kinetic model.
    groups : :obj:`numpy.ndarray` of shape (n_pairs,)
        Delay group of each pair.
    keep : :obj:`numpy.ndarray` of shape (n_pairs,)
        Pairs that inform the noise estimate and the correction.
    voxels : :obj:`numpy.ndarray` of shape (n_voxels,) or None, optional
        Voxels to score. Defaults to every voxel with finite data and predictions.
    shrinkage : :obj:`float`, optional
        Number of pairs at which the correction gets half weight.
    smooth : callable or None, optional
        Spatial smoothing. It takes and returns an array of shape (n_voxels, n_columns),
        in the voxel order of ``deltam``. If given, a second residual energy is computed from
        smoothed residuals: artifacts such as motion and intensity changes leave spatially
        coherent residuals, while thermal noise is independent across voxels, so smoothing
        raises the contrast between them. The unsmoothed energy is kept, because smoothing
        hides artifacts that are not coherent, such as excess thermal noise.
        The predicted signal is smoothed for the amplitude regression.
        Smoothing is normalized by the smoothed voxel mask, so voxels outside it do not leak in.

    Returns
    -------
    scores : :obj:`dict`
        ``energy``: mean squared standardized residual of each pair (about 1 for clean pairs).
        ``energy_smooth``: the same, for smoothed residuals (only if ``smooth`` is given).
        ``amplitude``: change in each pair's perfusion signal, relative to its prediction,
        from a regression of its residual on its predicted signal and a constant
        (-1 for a pair with no perfusion signal at all).
        The constant absorbs global offsets, which the energies catch,
        so that labeling failures are not masked by them.
        ``amplitude_t``: the amplitude in units of its standard error, ignoring spatial
        correlation.
        ``noise``: estimated noise standard deviation of each scored voxel.
    """
    finite = np.all(np.isfinite(deltam), axis=1) & np.all(np.isfinite(predicted), axis=1)
    voxels = finite if voxels is None else (np.asarray(voxels, dtype=bool) & finite)
    if not voxels.any():
        raise ValueError('No voxels to score: the mask is empty or the model fit failed in it.')

    n_pairs = deltam.shape[1]
    residuals = np.where(voxels[:, None], deltam - predicted, 0.0)
    expected = np.where(voxels[:, None], predicted, 0.0)
    # Variance of each pair's residual, relative to the noise variance of one pair
    variance_factor = np.ones(n_pairs)
    for group in np.unique(groups):
        members = np.flatnonzero(groups == group)
        kept = members[keep[members]]
        for pair in members:
            others = kept[kept != pair]
            if others.size < MIN_SAME_DELAY:
                continue

            weight = 0.0 if np.isinf(shrinkage) else others.size / (others.size + shrinkage)
            correction = weight * np.median(residuals[:, others], axis=1)
            residuals[:, pair] -= correction
            expected[:, pair] += correction
            variance_factor[pair] += weight**2 * _median_variance(others.size)

    def energy(values):
        values = values[voxels]
        noise_var = np.median(values[:, keep] ** 2 / variance_factor[keep], axis=1)
        noise_var /= CHI2_1_MEDIAN
        typical = np.median(noise_var)
        if typical <= 0:
            typical = max(np.mean(noise_var), np.finfo(float).tiny)
        noise_var = np.clip(noise_var, NOISE_CLIP[0] * typical, NOISE_CLIP[1] * typical)
        pair_var = noise_var[:, None] * variance_factor[None, :]
        return np.mean(values**2 / pair_var, axis=0), noise_var

    scores = {}
    scores['energy'], noise_var = energy(residuals)
    if smooth is not None:
        coverage = smooth(voxels[:, None].astype(float))[:, 0]
        coverage = np.maximum(coverage, np.finfo(float).tiny)[:, None]
        scores['energy_smooth'], _ = energy(smooth(residuals) / coverage)
        expected = smooth(expected) / coverage

    # Weighted regression of each pair's residual on a constant and its predicted signal
    weights = 1 / noise_var
    resid = residuals[voxels]
    signal = expected[voxels]
    amplitude = np.zeros(n_pairs)
    amplitude_t = np.zeros(n_pairs)
    for pair in range(n_pairs):
        design = np.column_stack([np.ones(signal.shape[0]), signal[:, pair]])
        gram = design.T @ (design * weights[:, None])
        if np.linalg.cond(gram) > 1e12:
            continue

        inverse = np.linalg.inv(gram)
        beta = inverse @ (design.T @ (weights * resid[:, pair]))
        amplitude[pair] = beta[1]
        amplitude_t[pair] = beta[1] / np.sqrt(inverse[1, 1] * variance_factor[pair])

    scores['amplitude'] = amplitude
    scores['amplitude_t'] = amplitude_t
    scores['noise'] = np.sqrt(noise_var)
    return scores


def flag_pairs(
    scores,
    groups,
    threshold=3.5,
    min_energy_ratio=1.25,
    min_amplitude=0.25,
    min_pairs_per_delay=1,
    use_amplitude=True,
):
    """Flag pairs whose scores are outliers among all pairs.

    A pair is flagged if the robust z-score of its log residual energy (raw or smoothed)
    exceeds ``threshold`` and that energy exceeds its median by the factor ``min_energy_ratio``,
    or if the robust z-score of its amplitude exceeds ``threshold`` in magnitude
    and its amplitude differs from the median by more than ``min_amplitude``.
    The effect size floors keep the robust z-scores, which can be very large when the pairs are
    much alike, from flagging pairs whose differences do not matter.

    Returns
    -------
    flags : :obj:`numpy.ndarray` of shape (n_pairs,)
        True for the pairs to exclude.
    severity : :obj:`numpy.ndarray` of shape (n_pairs,)
        The largest of the energy and amplitude z-scores, used to rank pairs.
    """
    n_pairs = scores['energy'].size
    flags = np.zeros(n_pairs, dtype=bool)
    severity = np.full(n_pairs, -np.inf)
    for key in ('energy', 'energy_smooth'):
        if key not in scores:
            continue

        log_energy = np.log(np.maximum(scores[key], np.finfo(float).tiny))
        energy_z = _robust_z(log_energy)
        energy_ratio = np.exp(log_energy - np.median(log_energy))
        flags |= (energy_z > threshold) & (energy_ratio > min_energy_ratio)
        severity = np.maximum(severity, energy_z)

    if use_amplitude:
        amplitude_z = _robust_z(scores['amplitude_t'])
        amplitude_shift = np.abs(scores['amplitude'] - np.median(scores['amplitude']))
        flags |= (np.abs(amplitude_z) > threshold) & (amplitude_shift > min_amplitude)
        severity = np.maximum(severity, np.abs(amplitude_z))

    # Keep the least severe pairs at each delay, so that every delay stays in the fit.
    for group in np.unique(groups):
        members = np.flatnonzero(groups == group)
        n_missing = min_pairs_per_delay - np.sum(~flags[members])
        if n_missing > 0:
            flagged = members[flags[members]]
            flags[flagged[np.argsort(severity[flagged])[:n_missing]]] = False

    return flags, severity


def detect_outlier_pairs_multipld(
    deltam_arr,
    scaled_m0data,
    plds,
    labeleff,
    t1blood,
    partition_coefficient,
    is_casl,
    tau=None,
    ti1=None,
    mask=None,
    shrinkage=2.0,
    smooth=None,
    threshold=3.5,
    min_energy_ratio=1.25,
    min_amplitude=0.25,
    use_amplitude=True,
    robust_start=True,
    min_pairs_per_delay=1,
    max_iter=5,
    fit_kwargs=None,
):
    """Find corrupted control-label pairs in multi-delay data and fit the model without them.

    The model is first fit to the voxelwise median of the pairs at each delay
    (``robust_start``), so that a corrupted pair does not distort the first prediction.
    Every pair is then scored against its prediction (see :func:`score_pairs`),
    outliers are flagged (see :func:`flag_pairs`), and the model is refit to the pairs that
    were not flagged. Scoring and refitting repeat until the flagged pairs do not change,
    so a pair flagged early can be restored once the fit it was scored against improves.

    Parameters
    ----------
    deltam_arr : :obj:`numpy.ndarray` of shape (n_voxels, n_pairs)
        Subtracted (control - label) signal, one column per control-label pair.
    scaled_m0data, plds, labeleff, t1blood, partition_coefficient, is_casl, tau, ti1
        As in :func:`~aslprep.utils.cbf.fit_deltam_multipld`.
    mask : :obj:`numpy.ndarray` of shape (n_voxels,) or None, optional
        Voxels used to score pairs, such as a brain mask. Voxels where the fit fails are
        always left out. Defaults to every voxel.
    shrinkage, smooth
        See :func:`score_pairs`.
    threshold, min_energy_ratio, min_amplitude, use_amplitude, min_pairs_per_delay
        See :func:`flag_pairs`.
    robust_start : :obj:`bool`, optional
        Whether to start from a fit to the median of the pairs at each delay.
    max_iter : :obj:`int`, optional
        Maximum number of times pairs are scored. At least 1.
    fit_kwargs : :obj:`dict` or None, optional
        Extra arguments to :func:`~aslprep.utils.cbf.fit_deltam_multipld`,
        such as ``n_jobs``.

    Returns
    -------
    result : :obj:`dict`
        ``keep``: boolean array of shape (n_pairs,), the pairs used in the final fit.
        ``fit``: the output of :func:`~aslprep.utils.cbf.fit_deltam_multipld` on those pairs.
        ``scores``: the scores of the last scoring pass (see :func:`score_pairs`), with
        ``severity`` added. If the flagged pairs did not settle within ``max_iter`` passes,
        these are the scores that selected ``keep``, computed before the final fit.
        None if there were too few pairs to score, or no delay with a pair to spare.
        ``history``: the pairs kept for each scoring pass.
        ``n_fits``: number of model fits.
    """
    deltam_arr = np.asarray(deltam_arr, dtype=float)
    n_voxels, n_pairs = deltam_arr.shape
    fit_kwargs = dict(fit_kwargs or {})
    if is_casl:
        tau = np.broadcast_to(np.asarray(tau, dtype=float), (n_pairs,)).copy()

    groups = _delay_groups(plds, tau, is_casl)
    mask = np.ones(n_voxels, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)

    def fit(data, keep):
        return fit_deltam_multipld(
            deltam_arr=data[:, keep],
            scaled_m0data=scaled_m0data,
            plds=plds[:, keep],
            labeleff=labeleff,
            t1blood=t1blood,
            partition_coefficient=partition_coefficient,
            is_casl=is_casl,
            tau=tau[keep] if is_casl else None,
            ti1=ti1,
            **fit_kwargs,
        )

    if max_iter < 1:
        raise ValueError(f'max_iter must be at least 1, but got {max_iter}')

    keep = np.ones(n_pairs, dtype=bool)
    _, group_sizes = np.unique(groups, return_counts=True)
    if n_pairs < MIN_PAIRS or np.all(group_sizes <= min_pairs_per_delay):
        # Too few pairs to score, or no delay has a pair to spare.
        result = fit(deltam_arr, keep)
        return {'keep': keep, 'fit': result, 'scores': None, 'history': [], 'n_fits': 1}

    if robust_start:
        start = deltam_arr.copy()
        for group in np.unique(groups):
            members = groups == group
            start[:, members] = np.median(deltam_arr[:, members], axis=1, keepdims=True)
        result = fit(start, keep)
        fit_is_final = False
    else:
        result = fit(deltam_arr, keep)
        fit_is_final = True
    n_fits = 1

    history = []
    scores = None
    for _ in range(max_iter):
        cbf, att, abat, abv, failures, _ = result
        predicted = predict_deltam_multipld(
            cbf,
            att,
            abat,
            abv,
            scaled_m0data,
            plds,
            labeleff,
            t1blood,
            partition_coefficient,
            is_casl,
            tau=tau,
            ti1=ti1,
        )
        scores = score_pairs(
            deltam_arr,
            predicted,
            groups,
            keep,
            voxels=mask & ~failures,
            shrinkage=shrinkage,
            smooth=smooth,
        )
        flags, severity = flag_pairs(
            scores,
            groups,
            threshold=threshold,
            min_energy_ratio=min_energy_ratio,
            min_amplitude=min_amplitude,
            min_pairs_per_delay=min_pairs_per_delay,
            use_amplitude=use_amplitude,
        )
        scores['severity'] = severity
        history.append(keep.copy())
        new_keep = ~flags
        if np.array_equal(new_keep, keep) and fit_is_final:
            break

        keep = new_keep
        result = fit(deltam_arr, keep)
        fit_is_final = True
        n_fits += 1

    return {'keep': keep, 'fit': result, 'scores': scores, 'history': history, 'n_fits': n_fits}
