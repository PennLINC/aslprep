"""Observation-level quality control for delta-M (or CBF) series.

Each observation (a control-label pair, or a precomputed delta-M or CBF volume) is compared with
the other observations acquired at the same delay. The main measure is the fraction of brain
voxels whose value lies more than ``z`` robust standard deviations from that voxel's median over
the same-delay observations.

This is the rule of :footcite:t:`harteveld2022perfusion`, who reject a subtraction image when more
than 20% of a region of interest lies beyond 1.5 SD of the voxelwise mean, with two changes that
an evaluation on simulated brain ASL (aslscan, with within-volume head motion and spin history)
and on real runs supported:

- Median and MAD replace mean and SD. Under realistic head motion several of ten pairs can be
  corrupted, and they inflate the mean and SD they are compared against.
- A stricter threshold (3.5 robust SD, more than 5% of the brain) flags pairs as outliers. It
  matched or beat the published thresholds on held-out simulations and flagged no pairs of
  motion-free runs; the 1.5 SD fraction is reported as well.

Nothing here excludes data: the measures are written to a table for review.

References
----------
.. footbibliography::
"""

import numpy as np
import pandas as pd

#: Scale factor from the median absolute deviation to the SD of a normal distribution.
MAD_TO_SD = 1.4826
#: Robust z-score above which a voxel counts as deviating, for the ``outlier`` flag.
Z_OUTLIER = 3.5
#: Fraction of brain voxels beyond ``Z_OUTLIER`` above which an observation is flagged.
FRACTION_OUTLIER = 0.05
#: The published threshold (as a robust z-score), reported for comparison.
Z_PUBLISHED = 1.5
#: Fewer same-delay observations than this give no stable median and MAD.
MIN_OBSERVATIONS = 4
#: Slice groups with fewer brain voxels than this are left out of the slice-level measure.
MIN_SLICE_VOXELS = 50


def robust_z(data):
    """Absolute robust z-scores of ``data`` (voxels x observations) along the observations.

    Voxels with a MAD of zero get z-scores of zero (they cannot deviate).
    """
    median = np.median(data, axis=1, keepdims=True)
    scale = MAD_TO_SD * np.median(np.abs(data - median), axis=1, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        z = np.abs(data - median) / scale
    return np.where(scale > 0, z, 0.0)


def delay_groups(delays, labeling_durations=None):
    """Integer group labels for observations that share a delay (and labeling duration)."""
    delays = np.round(np.asarray(delays, dtype=float), 6)
    if labeling_durations is None:
        keys = [(d,) for d in delays]
    else:
        durations = np.round(np.asarray(labeling_durations, dtype=float), 6)
        keys = list(zip(delays, durations, strict=True))
    unique = sorted(set(keys))
    return np.array([unique.index(k) for k in keys])


def slice_groups(slice_timing, n_slices):
    """Slice-group label of each slice: slices read at the same time, else one per slice.

    Under multiband acquisition (or any 2D readout), the slices that share a slice time were
    excited together, so a movement during the readout corrupts them together.
    """
    if slice_timing is None or len(slice_timing) != n_slices:
        return np.arange(n_slices)
    times = np.round(np.asarray(slice_timing, dtype=float), 6)
    unique = sorted(set(times))
    return np.array([unique.index(t) for t in times])


def observation_measures(deltam, groups, slices=None):
    """Per-observation outlier measures.

    Parameters
    ----------
    deltam : (V, N) :obj:`numpy.ndarray`
        Delta-M (or CBF) of the V brain voxels for the N observations.
    groups : (N,) array
        Delay group of each observation (see :func:`delay_groups`).
    slices : (V,) array or None
        Slice-group label of each voxel (see :func:`slice_groups`).

    Returns
    -------
    :obj:`pandas.DataFrame`
        One row per observation, with ``n_same_delay``, ``deviant_fraction``,
        ``deviant_fraction_z1p5``, ``max_slice_deviant_fraction``, ``mean_deltam_z`` and
        ``outlier``. Measures are NaN (and ``outlier`` False) for observations with fewer than
        :data:`MIN_OBSERVATIONS` at their delay.
    """
    deltam = np.asarray(deltam, dtype=float)
    groups = np.asarray(groups)
    n_obs = deltam.shape[1]
    if groups.shape != (n_obs,):
        raise ValueError(f'{groups.size} delay groups for {n_obs} observations')
    if slices is not None:
        slices = np.asarray(slices)
        if slices.shape != (deltam.shape[0],):
            raise ValueError(f'{slices.size} slice labels for {deltam.shape[0]} voxels')
        labels = [s for s in np.unique(slices) if np.sum(slices == s) >= MIN_SLICE_VOXELS]

    out = pd.DataFrame(
        {
            'n_same_delay': np.zeros(n_obs, dtype=int),
            'deviant_fraction': np.full(n_obs, np.nan),
            'deviant_fraction_z1p5': np.full(n_obs, np.nan),
            'max_slice_deviant_fraction': np.full(n_obs, np.nan),
            'mean_deltam_z': np.full(n_obs, np.nan),
        }
    )
    for group in np.unique(groups):
        idx = np.flatnonzero(groups == group)
        out.loc[idx, 'n_same_delay'] = len(idx)
        if len(idx) < MIN_OBSERVATIONS:
            continue

        z = robust_z(deltam[:, idx])
        deviant = z > Z_OUTLIER
        out.loc[idx, 'deviant_fraction'] = deviant.mean(axis=0)
        out.loc[idx, 'deviant_fraction_z1p5'] = (z > Z_PUBLISHED).mean(axis=0)
        if slices is not None and labels:
            per_slice = np.stack([deviant[slices == s].mean(axis=0) for s in labels])
            out.loc[idx, 'max_slice_deviant_fraction'] = per_slice.max(axis=0)

        # A shift of the whole image (labeling failure, global signal change) moves the
        # brain-mean delta-M even when few voxels pass the voxelwise threshold
        means = deltam[:, idx].mean(axis=0)
        out.loc[idx, 'mean_deltam_z'] = robust_z(means[None, :])[0] * np.sign(
            means - np.median(means)
        )

    out['outlier'] = out['deviant_fraction'] > FRACTION_OUTLIER
    return out


def observation_volumes(volume_types):
    """Which ASL volumes each delta-M observation came from, as ASLPrep pairs them.

    ASLPrep subtracts the k-th label volume from the k-th control volume. Precomputed delta-M
    or CBF volumes are observations by themselves.

    Returns
    -------
    observation_type : :obj:`str`
        ``'control-label'``, ``'deltam'`` or ``'cbf'``.
    volumes : :obj:`list` of :obj:`tuple`
        For each observation, ``(control, label)`` volume indices, or ``(volume,)``.
    """
    volume_types = list(volume_types)

    def where(kind):
        return [i for i, t in enumerate(volume_types) if t == kind]

    if where('deltam'):
        return 'deltam', [(i,) for i in where('deltam')]
    if where('label'):
        return 'control-label', list(zip(where('control'), where('label'), strict=True))
    if where('cbf'):
        return 'cbf', [(i,) for i in where('cbf')]
    raise ValueError('No control-label pairs, delta-M or CBF volumes in the ASL series.')
