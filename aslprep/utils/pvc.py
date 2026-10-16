"""Partial volume correction utilities.

Tissue partial-volume (PV) maps that model the ASL acquisition's effective resolution, and
local linear regression (LR) partial volume correction (Asllani et al., 2008).

See ``docs/specs/2026-10-15-partial-volume-correction-design.md`` (Sections 4 and 5).
"""

from __future__ import annotations

import re

import numpy as np
from scipy import ndimage, optimize, special

#: Relative effective-resolution factors (FWHM / voxel size), ordered (in-plane, in-plane, slice).
#: Single-shot 3D GRASE: Vidorreta et al. (2013).
#: 3D spiral: one GE protocol measured by Petr et al. (2018)
#: (4.6 x 4.3 x 11.8 mm FWHM at 3.8 x 3.8 x 4 mm), with the in-plane factors averaged.
EFFECTIVE_RESOLUTION_FACTORS = {
    'auto:grase-single-shot': (1.1, 1.1, 1.38),
    'auto:spiral-petr2018': (
        (4.6 / 3.8 + 4.3 / 3.8) / 2,
        (4.6 / 3.8 + 4.3 / 3.8) / 2,
        11.8 / 4,
    ),
}

#: LR status codes (spec Section 5.3).
PVC_STATUS = {
    0: 'outside the brain mask',
    1: 'two-tissue fit',
    2: 'gray matter only fit',
    3: 'white matter only fit',
    4: 'not identifiable',
}

_FWHM_TO_SIGMA = 1 / (2 * np.sqrt(2 * np.log(2)))


def get_slice_axis(metadata):
    """Return the array axis that slices are encoded along.

    Parameters
    ----------
    metadata : :obj:`dict`
        ASL metadata.

    Returns
    -------
    :obj:`int`
        0, 1 or 2, from ``SliceEncodingDirection`` (``k`` when it is absent).
    """
    direction = metadata.get('SliceEncodingDirection', 'k')
    return 'ijk'.index(direction[0])


def kernel_size_to_axes(kernel_size, slice_axis):
    """Map a size given as (in-plane, in-plane, slice) to array axes (i, j, k).

    The first in-plane value goes to the lower-numbered in-plane axis.
    """
    inplane = [ax for ax in range(3) if ax != slice_axis]
    out = [0, 0, 0]
    out[inplane[0]] = kernel_size[0]
    out[inplane[1]] = kernel_size[1]
    out[slice_axis] = kernel_size[2]
    return tuple(out)


def has_varying_slice_timing(metadata):
    """Whether slices are acquired at different times, so their delays differ."""
    slice_timing = metadata.get('SliceTiming')
    if slice_timing is None:
        return False
    return np.ptp(np.asarray(slice_timing, dtype=float)) > 0


def check_lr_pvc_supported(metadata, processing_target, is_multi_pld, kernel_size):
    """Raise an error if local-regression PVC cannot be applied to these data.

    Parameters
    ----------
    metadata : :obj:`dict`
        ASL metadata.
    processing_target : :obj:`str`
        The volumes CBF is computed from (``'cbf'`` for scanner-provided CBF).
    is_multi_pld : :obj:`bool`
        Whether the data have more than one delay.
    kernel_size : :obj:`tuple` of :obj:`int`
        The kernel size, ordered (in-plane, in-plane, slice).

    Raises
    ------
    ValueError
        For multi-delay data, scanner-provided CBF, or a kernel spanning several slices of a
        slice-timed acquisition.
    """
    if is_multi_pld:
        raise ValueError(
            '--pvc-method local-regression only supports single-delay data. A kinetic fit of a '
            'mixed signal is not a linear mixture of the tissue fits, because gray and white '
            'matter arrival times differ. Use --basil, whose partial volume correction models '
            'both tissues.'
        )
    if processing_target == 'cbf':
        raise ValueError(
            '--pvc-method local-regression cannot be applied to scanner-provided CBF, whose '
            'calibration is unknown. Use --basil if delta-M volumes are available.'
        )
    if any(int(k) < 1 or int(k) % 2 == 0 for k in kernel_size):
        raise ValueError(f'--pvc-kernel-size values must be odd and positive, got {kernel_size}.')
    if int(kernel_size[2]) > 1 and has_varying_slice_timing(metadata):
        raise ValueError(
            '--pvc-kernel-size must span a single slice for slice-timed acquisitions: each '
            "slice's tissue signal depends on its own delay, through tissue T1 and arrival "
            f'time, so tissue amplitudes are not constant across slices. Got {kernel_size}.'
        )


def resolve_effective_resolution(option, metadata, voxel_size, slice_axis):
    """Resolve ``--asl-effective-resolution`` to a total response FWHM per array axis.

    Parameters
    ----------
    option : :obj:`str` or :obj:`list` of :obj:`float`
        ``'auto'``, ``'nominal'``, or three FWHM values in mm, ordered
        (in-plane, in-plane, slice).
    metadata : :obj:`dict`
        ASL metadata.
    voxel_size : :obj:`tuple` of :obj:`float`
        The aslref voxel size along array axes (i, j, k), in mm.
    slice_axis : :obj:`int`
        The slice-encoding array axis.

    Returns
    -------
    :obj:`dict`
        ``fwhm`` (mm, per array axis), ``source`` (a label), ``factors`` (FWHM / voxel size,
        ordered in-plane, in-plane, slice, or None), ``metadata_used`` and ``warnings``.

    Raises
    ------
    ValueError
        If user-supplied values are not positive, or are smaller than the voxel size.
    """
    voxel_size = tuple(float(v) for v in voxel_size)
    warnings = []
    metadata_used = {
        key: metadata[key]
        for key in (
            'MRAcquisitionType',
            'PulseSequenceType',
            'NumberShots',
            'SliceEncodingDirection',
            'AcquisitionVoxelSize',
        )
        if key in metadata
    }
    if 'AcquisitionVoxelSize' in metadata:
        acq = np.asarray(metadata['AcquisitionVoxelSize'], dtype=float)
        if acq.size == 3 and not np.allclose(np.sort(acq), np.sort(voxel_size), rtol=0.01):
            warnings.append(
                f'AcquisitionVoxelSize ({acq.tolist()}) differs from the reconstructed voxel size '
                f'({list(voxel_size)}). The effective-resolution defaults assume the '
                'reconstruction did not interpolate.'
            )

    if isinstance(option, str):
        option = option.lower()
    if option == 'nominal':
        return {
            'fwhm': voxel_size,
            'source': 'nominal',
            'factors': None,
            'metadata_used': metadata_used,
            'warnings': warnings,
        }

    if option == 'auto':
        source = _auto_source(metadata)
        if source in ('auto:2d', 'nominal-fallback'):
            if source == 'nominal-fallback':
                warnings.append(
                    'The effective resolution of this acquisition could not be inferred from '
                    'its metadata (MRAcquisitionType, PulseSequenceType, NumberShots), so the '
                    'tissue maps model only the voxel itself. '
                    'Set --asl-effective-resolution to model readout blurring.'
                )
            return {
                'fwhm': voxel_size,
                'source': source,
                'factors': None,
                'metadata_used': metadata_used,
                'warnings': warnings,
            }

        factors = EFFECTIVE_RESOLUTION_FACTORS[source]
        if source == 'auto:spiral-petr2018':
            warnings.append(
                'The 3D spiral effective-resolution factors were measured for one GE protocol '
                '(Petr et al., 2018) and are an approximation for other protocols. '
                'Set --asl-effective-resolution if the effective resolution is known.'
            )
        factors_axes = kernel_size_to_axes(factors, slice_axis)
        fwhm = tuple(f * v for f, v in zip(factors_axes, voxel_size, strict=True))
        return {
            'fwhm': fwhm,
            'source': source,
            'factors': list(factors),
            'metadata_used': metadata_used,
            'warnings': warnings,
        }

    values = [float(v) for v in option]
    if len(values) != 3:
        raise ValueError(f'Expected three effective-resolution values, got {values}.')
    if not all(np.isfinite(v) and v > 0 for v in values):
        raise ValueError(f'Effective-resolution values must be finite and positive, got {values}.')
    fwhm = kernel_size_to_axes(values, slice_axis)
    too_small = [
        (ax, f, v) for ax, (f, v) in enumerate(zip(fwhm, voxel_size, strict=True)) if f < v - 1e-6
    ]
    if too_small:
        details = ', '.join(f'axis {ax}: {f} mm < {v} mm' for ax, f, v in too_small)
        raise ValueError(
            'The effective resolution cannot be finer than the voxel size '
            f'({details}). Use "nominal" for a response of exactly one voxel.'
        )
    return {
        'fwhm': tuple(fwhm),
        'source': 'user',
        'factors': None,
        'metadata_used': metadata_used,
        'warnings': warnings,
    }


def _auto_source(metadata):
    acq_type = str(metadata.get('MRAcquisitionType', '')).upper()
    sequence = str(metadata.get('PulseSequenceType', ''))
    if acq_type == '2D':
        return 'auto:2d'
    if acq_type == '3D':
        if re.search('grase', sequence, flags=re.IGNORECASE):
            n_shots = metadata.get('NumberShots')
            if n_shots is None or np.all(np.asarray(n_shots) == 1):
                return 'auto:grase-single-shot'
            return 'nominal-fallback'
        if re.search('spiral', sequence, flags=re.IGNORECASE):
            return 'auto:spiral-petr2018'
    return 'nominal-fallback'


def box_gaussian_profile(x, width, sigma):
    """Evaluate the box (``width``) convolved with a Gaussian (``sigma``), up to a constant.

    The peak is at zero. A ``sigma`` of zero gives the box itself.
    """
    x = np.asarray(x, dtype=float)
    half = width / 2
    if sigma <= 0:
        return (np.abs(x) < half).astype(float) + 0.5 * (np.abs(x) == half)
    scale = sigma * np.sqrt(2)
    return special.erf((x + half) / scale) - special.erf((x - half) / scale)


def box_gaussian_fwhm(width, sigma):
    """Return the FWHM of the box convolved with a Gaussian."""
    if sigma <= 0:
        return float(width)
    peak = box_gaussian_profile(0.0, width, sigma)

    def _half(x):
        return box_gaussian_profile(x, width, sigma) - peak / 2

    upper = width / 2 + 10 * sigma
    return 2 * optimize.brentq(_half, 0.0, upper, xtol=1e-12 * upper)


def solve_box_gaussian_sigma(fwhm, width):
    """Find the Gaussian sigma that gives the box-Gaussian kernel the requested FWHM.

    Parameters
    ----------
    fwhm : :obj:`float`
        The requested total FWHM. Must be at least ``width``.
    width : :obj:`float`
        The box width (the voxel size).

    Returns
    -------
    :obj:`float`
        The Gaussian sigma, in the units of ``width``. Zero when ``fwhm`` equals ``width``.
    """
    if fwhm < width - 1e-9 * width:
        raise ValueError(f'FWHM ({fwhm}) must be at least the box width ({width}).')
    if fwhm <= width * (1 + 1e-9):
        return 0.0

    # The FWHM increases monotonically with sigma, from the box width up.
    return optimize.brentq(
        lambda sigma: box_gaussian_fwhm(width, sigma) - fwhm,
        1e-9 * width,
        fwhm * _FWHM_TO_SIGMA * 2,
        xtol=1e-12 * width,
    )


def discrete_box_gaussian_kernel(width, sigma, step):
    """Integrate the box-Gaussian kernel over bins of size ``step``, centred on zero.

    The kernel is truncated at ``width / 2 + 4 * sigma`` and renormalized to sum to one.

    Parameters
    ----------
    width : :obj:`float`
        The box width.
    sigma : :obj:`float`
        The Gaussian sigma (zero for the box alone).
    step : :obj:`float`
        The bin size (the subvoxel size).

    Returns
    -------
    kernel : :obj:`numpy.ndarray`
        Odd-length, symmetric weights that sum to one.
    """
    half = width / 2
    half_width = half + 4 * sigma
    n_half = int(np.ceil(half_width / step - 0.5 - 1e-9))
    centres = np.arange(-n_half, n_half + 1) * step
    lo = centres - step / 2
    hi = centres + step / 2

    if sigma <= 0:
        weights = np.clip(np.minimum(hi, half) - np.maximum(lo, -half), 0, None)
    else:
        # Antiderivative of Phi(t / sigma): t * Phi(t / sigma) + sigma * phi(t / sigma)
        def _cdf_integral(t):
            z = t / sigma
            return t * special.ndtr(z) + sigma * np.exp(-0.5 * z**2) / np.sqrt(2 * np.pi)

        # The kernel is Phi((x + half) / sigma) - Phi((x - half) / sigma).
        weights = (_cdf_integral(hi + half) - _cdf_integral(lo + half)) - (
            _cdf_integral(hi - half) - _cdf_integral(lo - half)
        )

    weights = np.clip(weights, 0, None)
    return weights / weights.sum()


def choose_subdivision(voxel_size, max_subvoxel):
    """Return the smallest odd subdivision whose subvoxel is no larger than ``max_subvoxel``."""
    n = int(np.ceil(voxel_size / max_subvoxel - 1e-9))
    return n if n % 2 == 1 else n + 1


def build_pv_maps(
    tpms,
    tpm_affine,
    ref_shape,
    ref_affine,
    ref2anat,
    fwhm,
    max_subvoxel=None,
    max_points=8_000_000,
    slice_axis=2,
):
    """Build tissue PV maps on the ASL reference grid that model its effective resolution.

    Each tissue probability map is treated as a continuous field sampled at the anatomical voxel
    centres, reconstructed by linear interpolation, and zero outside the anatomical field of
    view. It is sampled on a supersampled grid aligned with the ASL reference, convolved with
    the effective response (a voxel-wide box, convolved with a Gaussian when the FWHM exceeds
    the voxel size), and read out at the ASL voxel centres.

    Parameters
    ----------
    tpms : :obj:`list` of :obj:`numpy.ndarray`
        3D tissue probability maps on the anatomical grid.
    tpm_affine : (4, 4) :obj:`numpy.ndarray`
        The anatomical grid's voxel-to-RAS affine.
    ref_shape : :obj:`tuple` of :obj:`int`
        The ASL reference grid's shape.
    ref_affine : (4, 4) :obj:`numpy.ndarray`
        The ASL reference grid's voxel-to-RAS affine.
    ref2anat : (4, 4) :obj:`numpy.ndarray`
        The RAS point mapping from the ASL reference to the anatomical image.
    fwhm : :obj:`tuple` of :obj:`float`
        The total response FWHM along each ASL reference axis, in mm.
    max_subvoxel : :obj:`float` or None
        The largest subvoxel, in mm. Default: half the smallest anatomical voxel, at most 1 mm.
    max_points : :obj:`int`
        The number of supersampled points processed at once, which bounds the memory used.
    slice_axis : :obj:`int`
        The ASL reference axis along which slabs are processed.

    Returns
    -------
    pv_maps : :obj:`list` of :obj:`numpy.ndarray`
        One PV map per tissue, on the ASL reference grid.
    info : :obj:`dict`
        The subdivision, the subvoxel size, the Gaussian sigmas and the kernels used.
    """
    ref_shape = tuple(int(s) for s in ref_shape[:3])
    voxel_size = np.sqrt((ref_affine[:3, :3] ** 2).sum(axis=0))
    tpm_voxel = np.sqrt((tpm_affine[:3, :3] ** 2).sum(axis=0))
    if max_subvoxel is None:
        max_subvoxel = min(0.5 * tpm_voxel.min(), 1.0)

    subdiv = [choose_subdivision(v, max_subvoxel) for v in voxel_size]
    step = voxel_size / subdiv
    sigmas = [solve_box_gaussian_sigma(f, v) for f, v in zip(fwhm, voxel_size, strict=True)]
    kernels = [
        discrete_box_gaussian_kernel(v, s, st)
        for v, s, st in zip(voxel_size, sigmas, step, strict=True)
    ]
    # The halo, in ASL voxels, covers the kernel's half-width beyond the voxel.
    kernel_half = [(len(k) - 1) // 2 for k in kernels]
    halo = [
        int(np.ceil(max(h - (s - 1) // 2, 0) / s))
        for h, s in zip(kernel_half, subdiv, strict=True)
    ]

    # Supersampled index coordinates (in ASL voxel units) -> anatomical voxel coordinates
    ref2anat_vox = np.linalg.inv(tpm_affine) @ ref2anat @ ref_affine

    # How many ASL slices to process at once, along slice_axis
    other_axes = [ax for ax in range(3) if ax != slice_axis]
    n_other = np.prod(
        [(ref_shape[ax] + 2 * halo[ax]) * subdiv[ax] for ax in other_axes], dtype=float
    )
    per_slice = n_other * subdiv[slice_axis]
    slab = max(1, int(max_points // max(per_slice, 1)) - 2 * halo[slice_axis])
    slab = min(slab, ref_shape[slice_axis])

    tpms = [np.asarray(tpm, dtype=np.float64) for tpm in tpms]
    pv_maps = [np.zeros(ref_shape, dtype=np.float64) for _ in tpms]

    for start in range(0, ref_shape[slice_axis], slab):
        stop = min(start + slab, ref_shape[slice_axis])
        # ASL voxel index ranges, with halos
        ranges = []
        for ax in range(3):
            if ax == slice_axis:
                lo, hi = start - halo[ax], stop + halo[ax]
            else:
                lo, hi = -halo[ax], ref_shape[ax] + halo[ax]
            ranges.append((lo, hi))

        # Supersampled coordinates along each axis, in ASL voxel units
        axes_coords = []
        for ax in range(3):
            lo, hi = ranges[ax]
            offsets = (np.arange(subdiv[ax]) - (subdiv[ax] - 1) / 2) / subdiv[ax]
            axes_coords.append((np.arange(lo, hi)[:, None] + offsets[None, :]).ravel())

        grid = np.meshgrid(*axes_coords, indexing='ij')
        shape_ss = grid[0].shape
        points = np.vstack([g.ravel() for g in grid] + [np.ones(grid[0].size)])
        del grid
        anat_vox = (ref2anat_vox @ points)[:3]
        del points

        for i_tissue, tpm in enumerate(tpms):
            sampled = ndimage.map_coordinates(
                tpm, anat_vox, order=1, mode='constant', cval=0.0, prefilter=False
            ).reshape(shape_ss)
            np.clip(sampled, 0, 1, out=sampled)

            # Filter and decimate one axis at a time. The centre subvoxel of each ASL voxel
            # is at index (s - 1) / 2 within its block.
            for ax in range(3):
                sampled = ndimage.correlate1d(
                    sampled, kernels[ax], axis=ax, mode='constant', cval=0.0
                )
                centre = (subdiv[ax] - 1) // 2
                index = [slice(None)] * 3
                index[ax] = slice(centre, None, subdiv[ax])
                sampled = sampled[tuple(index)]

            # Drop the halos
            index = []
            for ax in range(3):
                index.append(slice(halo[ax], sampled.shape[ax] - halo[ax]))
            out_index = [slice(None)] * 3
            out_index[slice_axis] = slice(start, stop)
            pv_maps[i_tissue][tuple(out_index)] = sampled[tuple(index)]

        del anat_vox

    info = {
        'subdivision': [int(s) for s in subdiv],
        'subvoxel_mm': [float(s) for s in step],
        'sigma_mm': [float(s) for s in sigmas],
        'kernels': [k.tolist() for k in kernels],
        'halo_voxels': [int(h) for h in halo],
    }
    return pv_maps, info


def _box_sum(arr, size):
    """Sum over a centred box of ``size`` voxels, with zero padding beyond the image."""
    return ndimage.uniform_filter(arr, size=size, mode='constant', cval=0.0) * np.prod(size)


def _shifted(padded, offset, pad, shape):
    index = tuple(slice(p + o, p + o + n) for p, o, n in zip(pad, offset, shape, strict=True))
    return padded[index]


def local_regression_pvc(
    y,
    pv_gm,
    pv_wm,
    mask,
    valid,
    kernel_size=(5, 5, 1),
    tau=0.5,
    kappa_max=100.0,
    c_max=0.05,
    min_n_two=5,
    min_n_one=4,
):
    """Estimate gray and white matter signal amplitudes by local linear regression.

    For each voxel in ``mask``, the signal ``y`` of the valid voxels in a box around it is
    regressed on the tissue PV maps, without an intercept (Asllani et al., 2008).
    See Section 5 of the design spec for the status rules.

    Parameters
    ----------
    y : :obj:`numpy.ndarray`
        The signal to unmix (3D).
    pv_gm, pv_wm : :obj:`numpy.ndarray`
        Gray and white matter PV maps (3D).
    mask : :obj:`numpy.ndarray` of :obj:`bool`
        The brain mask.
    valid : :obj:`numpy.ndarray` of :obj:`bool`
        Voxels whose signal can be used, and which can be calibrated.
    kernel_size : :obj:`tuple` of :obj:`int`
        The box size along each array axis. All values must be odd.
    tau, kappa_max, c_max : :obj:`float`
        Minimum tissue support, maximum condition number of the normalized Gram matrix, and
        maximum contamination ratio for single-tissue fits.
    min_n_two, min_n_one : :obj:`int`
        Minimum valid voxels for two- and single-tissue fits.

    Returns
    -------
    :obj:`dict` of :obj:`numpy.ndarray`
        ``beta_gm``, ``beta_wm``, ``status``, ``condition``, ``noise_gm``, ``noise_wm``,
        ``rmse`` and ``n``.
    """
    kernel_size = tuple(int(k) for k in kernel_size)
    if any(k < 1 or k % 2 == 0 for k in kernel_size):
        raise ValueError(f'Kernel sizes must be odd and positive, got {kernel_size}.')

    y = np.asarray(y, dtype=np.float64)
    pv_gm = np.asarray(pv_gm, dtype=np.float64)
    pv_wm = np.asarray(pv_wm, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    valid = (
        np.asarray(valid, dtype=bool)
        & mask
        & np.isfinite(y)
        & np.isfinite(pv_gm)
        & np.isfinite(pv_wm)
        & ((np.nan_to_num(pv_gm) + np.nan_to_num(pv_wm)) > 0)
    )
    v = valid.astype(np.float64)
    yv = np.where(valid, y, 0.0)
    pg = np.where(valid, pv_gm, 0.0)
    pw = np.where(valid, pv_wm, 0.0)

    n = np.rint(_box_sum(v, kernel_size))
    s_g = _box_sum(pg, kernel_size)
    s_w = _box_sum(pw, kernel_size)
    s_gg = _box_sum(pg * pg, kernel_size)
    s_ww = _box_sum(pw * pw, kernel_size)
    s_gw = _box_sum(pg * pw, kernel_size)
    b_g = _box_sum(pg * yv, kernel_size)
    b_w = _box_sum(pw * yv, kernel_size)

    # The box filter can leave tiny negative values where the sums should be zero.
    tiny = 1e-12
    s_gg = np.where(s_gg > tiny, s_gg, 0.0)
    s_ww = np.where(s_ww > tiny, s_ww, 0.0)

    with np.errstate(divide='ignore', invalid='ignore'):
        both = (s_gg > 0) & (s_ww > 0)
        r = np.where(both, s_gw / np.sqrt(s_gg * s_ww), 1.0)
        r = np.clip(np.abs(r), 0, 1)
        kappa = np.where(r < 1, (1 + r) / (1 - r), np.inf)
        c_g = np.where(s_gg > 0, s_gw / s_gg, np.inf)
        c_w = np.where(s_ww > 0, s_gw / s_ww, np.inf)

    centre_ok = mask & valid
    status = np.zeros(y.shape, dtype=np.uint8)
    status[mask] = 4
    two = centre_ok & (n >= min_n_two) & (s_g >= tau) & (s_w >= tau) & both & (kappa <= kappa_max)
    gm_only = centre_ok & ~two & (n >= min_n_one) & (s_g >= tau) & (s_gg > 0) & (c_g <= c_max)
    wm_only = (
        centre_ok & ~two & ~gm_only & (n >= min_n_one) & (s_w >= tau) & (s_ww > 0) & (c_w <= c_max)
    )
    status[two] = 1
    status[gm_only] = 2
    status[wm_only] = 3

    beta_gm = np.zeros(y.shape)
    beta_wm = np.zeros(y.shape)
    noise_gm = np.zeros(y.shape)
    noise_wm = np.zeros(y.shape)

    det = s_gg[two] * s_ww[two] - s_gw[two] ** 2
    beta_gm[two] = (s_ww[two] * b_g[two] - s_gw[two] * b_w[two]) / det
    beta_wm[two] = (s_gg[two] * b_w[two] - s_gw[two] * b_g[two]) / det
    noise_gm[two] = np.sqrt(s_ww[two] / det)
    noise_wm[two] = np.sqrt(s_gg[two] / det)

    beta_gm[gm_only] = b_g[gm_only] / s_gg[gm_only]
    noise_gm[gm_only] = 1 / np.sqrt(s_gg[gm_only])
    beta_wm[wm_only] = b_w[wm_only] / s_ww[wm_only]
    noise_wm[wm_only] = 1 / np.sqrt(s_ww[wm_only])

    # Residual sum of squares, computed directly rather than from the expansion,
    # so near-exact fits are not lost to cancellation.
    fitted = status > 0
    fitted &= status < 4
    pad = [k // 2 for k in kernel_size]
    shape = y.shape
    y_pad = np.pad(yv, [(p, p) for p in pad])
    g_pad = np.pad(pg, [(p, p) for p in pad])
    w_pad = np.pad(pw, [(p, p) for p in pad])
    v_pad = np.pad(v, [(p, p) for p in pad])
    rss = np.zeros(shape)
    for offset in np.ndindex(*kernel_size):
        offset = tuple(o - p for o, p in zip(offset, pad, strict=True))
        resid = _shifted(y_pad, offset, pad, shape) - (
            beta_gm * _shifted(g_pad, offset, pad, shape)
            + beta_wm * _shifted(w_pad, offset, pad, shape)
        )
        rss += _shifted(v_pad, offset, pad, shape) * resid**2

    n_params = np.where(status == 1, 2, 1)
    rmse = np.zeros(shape)
    rmse[fitted] = np.sqrt(rss[fitted] / (n[fitted] - n_params[fitted]))

    condition = np.zeros(shape)
    defined = mask & both
    condition[defined] = np.minimum(kappa[defined], 1e6)

    return {
        'beta_gm': beta_gm,
        'beta_wm': beta_wm,
        'status': status,
        'condition': condition,
        'noise_gm': noise_gm,
        'noise_wm': noise_wm,
        'rmse': rmse,
        'n': n,
    }


def calibrate_pvc(fit, m0, pv_gm, pv_wm, min_pv=0.1):
    """Turn LR signal amplitudes into tissue CBF, calibrated by the centre voxel's M0.

    Parameters
    ----------
    fit : :obj:`dict`
        The output of :func:`local_regression_pvc`, fit to ``y = CBF * M0``.
    m0 : :obj:`numpy.ndarray`
        The M0 image used to quantify CBF.
    pv_gm, pv_wm : :obj:`numpy.ndarray`
        The PV maps used for the fit.
    min_pv : :obj:`float`
        A tissue's estimate is written only where that tissue's PV is at least this.

    Returns
    -------
    :obj:`dict` of :obj:`numpy.ndarray`
        ``cbf_gm``, ``cbf_wm``, ``se_gm``, ``se_wm`` (zero where not estimated), and the boolean
        masks ``valid_gm`` and ``valid_wm`` of where they are.
    """
    status = fit['status']
    m0 = np.asarray(m0, dtype=np.float64)
    calibratable = np.isfinite(m0) & (m0 > 0)
    valid_gm = np.isin(status, (1, 2)) & (pv_gm >= min_pv) & calibratable
    valid_wm = np.isin(status, (1, 3)) & (pv_wm >= min_pv) & calibratable

    out = {}
    for tissue, tissue_valid in (('gm', valid_gm), ('wm', valid_wm)):
        cbf = np.zeros(status.shape)
        se = np.zeros(status.shape)
        cbf[tissue_valid] = fit[f'beta_{tissue}'][tissue_valid] / m0[tissue_valid]
        se[tissue_valid] = (
            fit[f'noise_{tissue}'][tissue_valid] * fit['rmse'][tissue_valid] / m0[tissue_valid]
        )
        out[f'cbf_{tissue}'] = cbf
        out[f'se_{tissue}'] = se
        out[f'valid_{tissue}'] = tissue_valid
    return out


#: CSF fit status codes. Every fit includes CSF; the code says which other tissues it includes.
CSF_FIT_STATUS = {
    0: 'outside the brain mask',
    1: 'gray matter, white matter and CSF fit',
    2: 'gray matter and CSF fit (white matter contamination below the limit)',
    3: 'white matter and CSF fit (gray matter contamination below the limit)',
    4: 'not identifiable',
    5: 'CSF only fit (gray and white matter contamination below the limit)',
}

# Tissue indices (GM, WM, CSF) of each fit, in order of precedence, with its status code.
_CSF_SUBSETS = ((1, (0, 1, 2)), (2, (0, 2)), (3, (1, 2)), (5, (2,)))


def local_regression_csf(
    y,
    pv_gm,
    pv_wm,
    pv_csf,
    mask,
    valid,
    kernel_size=(5, 5, 1),
    tau=0.5,
    kappa_max=100.0,
    c_max=0.05,
):
    """Estimate the CSF signal amplitude by local regression on GM, WM and CSF, for QC.

    CSF has no perfusion, so its estimate should be near zero. A clearly nonzero estimate points
    at partial volume map errors (segmentation, registration, effective resolution),
    macrovascular or choroid plexus signal, or a poorly fitting kinetic model.
    The gray and white matter estimates of this fit are not used: the corrected GM and WM maps
    come from :func:`local_regression_pvc`, which assumes CSF has no perfusion.

    The fits are tried in the order of :data:`CSF_FIT_STATUS` (all three tissues, then CSF with
    one tissue, then CSF alone). A tissue is left out only when its contamination of the CSF
    estimate (the coefficient of its PV map on the CSF estimate, ``[G_S^-1 X_S^T p_o]_CSF``)
    is at most ``c_max``. Each fit needs ``tau`` of support for every included tissue, a
    normalized Gram matrix with condition number at most ``kappa_max``, and at least three
    residual degrees of freedom.

    Parameters
    ----------
    y, pv_gm, pv_wm, pv_csf, mask, valid, kernel_size, tau, kappa_max, c_max
        As for :func:`local_regression_pvc`, with the CSF PV map added.

    Returns
    -------
    :obj:`dict` of :obj:`numpy.ndarray`
        ``beta_csf``, ``status`` (see :data:`CSF_FIT_STATUS`), ``noise_csf``, ``rmse`` and
        ``n``.
    """
    kernel_size = tuple(int(k) for k in kernel_size)
    if any(k < 1 or k % 2 == 0 for k in kernel_size):
        raise ValueError(f'Kernel sizes must be odd and positive, got {kernel_size}.')

    y = np.asarray(y, dtype=np.float64)
    pvs = [np.asarray(p, dtype=np.float64) for p in (pv_gm, pv_wm, pv_csf)]
    mask = np.asarray(mask, dtype=bool)
    finite = np.isfinite(y)
    for p in pvs:
        finite &= np.isfinite(p)
    valid = (
        np.asarray(valid, dtype=bool) & mask & finite & (sum(np.nan_to_num(p) for p in pvs) > 0)
    )
    v = valid.astype(np.float64)
    yv = np.where(valid, y, 0.0)
    ps = [np.where(valid, p, 0.0) for p in pvs]

    n = np.rint(_box_sum(v, kernel_size))
    support = [_box_sum(p, kernel_size) for p in ps]
    gram = np.empty((*y.shape, 3, 3))
    for i in range(3):
        for j in range(i, 3):
            gram[..., i, j] = gram[..., j, i] = _box_sum(ps[i] * ps[j], kernel_size)
    rhs = np.stack([_box_sum(p * yv, kernel_size) for p in ps], axis=-1)

    status = np.zeros(y.shape, dtype=np.uint8)
    status[mask] = 4
    beta = np.zeros((*y.shape, 3))
    noise_csf = np.zeros(y.shape)
    undecided = mask & valid
    tiny = 1e-12

    for code, subset in _CSF_SUBSETS:
        subset = list(subset)
        omitted = [t for t in range(3) if t not in subset]
        k = len(subset)
        candidate = undecided & (n >= k + 3)
        for t in subset:
            candidate &= (support[t] >= tau) & (gram[..., t, t] > tiny)
        if not candidate.any():
            continue

        g = gram[candidate][:, subset][:, :, subset]
        d = np.sqrt(np.einsum('...ii->...i', g))
        normalized = g / (d[:, :, None] * d[:, None, :])
        eig = np.linalg.eigvalsh(normalized)
        with np.errstate(divide='ignore', invalid='ignore'):
            kappa = np.where(eig[:, 0] > tiny, eig[:, -1] / eig[:, 0], np.inf)
        ok = kappa <= kappa_max

        g_ok = g[ok]
        g_inv = np.linalg.inv(g_ok) if g_ok.size else g_ok
        csf_row = subset.index(2)
        # Contamination of the CSF coefficient by each omitted tissue
        for o in omitted:
            cross = gram[candidate][ok][:, subset, o]
            leak = np.einsum('vij,vj->vi', g_inv, cross)[:, csf_row]
            keep = np.abs(leak) <= c_max
            ok_index = np.flatnonzero(ok)
            ok[ok_index[~keep]] = False
            g_inv = g_inv[keep]

        voxels = tuple(idx[ok] for idx in np.nonzero(candidate))
        b = rhs[voxels][:, subset]
        coef = np.einsum('vij,vj->vi', g_inv, b)
        for row, t in enumerate(subset):
            beta[(*voxels, np.full(len(voxels[0]), t))] = coef[:, row]
        noise_csf[voxels] = np.sqrt(g_inv[:, csf_row, csf_row])
        status[voxels] = code
        undecided[voxels] = False

    # Residuals, computed directly over the kernel offsets
    fitted = np.isin(status, [code for code, _ in _CSF_SUBSETS])
    pad = [kk // 2 for kk in kernel_size]
    shape = y.shape
    y_pad = np.pad(yv, [(p, p) for p in pad])
    p_pad = [np.pad(p, [(q, q) for q in pad]) for p in ps]
    v_pad = np.pad(v, [(q, q) for q in pad])
    rss = np.zeros(shape)
    for offset in np.ndindex(*kernel_size):
        offset = tuple(o - p for o, p in zip(offset, pad, strict=True))
        fit_value = sum(beta[..., t] * _shifted(p_pad[t], offset, pad, shape) for t in range(3))
        resid = _shifted(y_pad, offset, pad, shape) - fit_value
        rss += _shifted(v_pad, offset, pad, shape) * resid**2

    n_params = np.zeros(shape)
    for code, subset in _CSF_SUBSETS:
        n_params[status == code] = len(subset)
    rmse = np.zeros(shape)
    rmse[fitted] = np.sqrt(rss[fitted] / (n[fitted] - n_params[fitted]))

    return {
        'beta_csf': beta[..., 2],
        'status': status,
        'noise_csf': noise_csf,
        'rmse': rmse,
        'n': n,
    }


def calibrate_csf(fit, m0, pv_csf, min_pv=0.1):
    """Calibrate the CSF amplitude of :func:`local_regression_csf` with the centre voxel's M0.

    Returns
    -------
    :obj:`dict` of :obj:`numpy.ndarray`
        ``cbf_csf``, ``se_csf`` (zero where not estimated) and the boolean mask ``valid_csf``.
    """
    m0 = np.asarray(m0, dtype=np.float64)
    calibratable = np.isfinite(m0) & (m0 > 0)
    fitted = np.isin(fit['status'], [code for code, _ in _CSF_SUBSETS])
    valid = fitted & (np.asarray(pv_csf) >= min_pv) & calibratable
    cbf = np.zeros(m0.shape)
    se = np.zeros(m0.shape)
    cbf[valid] = fit['beta_csf'][valid] / m0[valid]
    se[valid] = fit['noise_csf'][valid] * fit['rmse'][valid] / m0[valid]
    return {'cbf_csf': cbf, 'se_csf': se, 'valid_csf': valid}


def summarize_csf_qc(cbf_csf, valid_csf, pv_csf, min_pv=0.5):
    """QC summary of the CSF estimate: median and spread where CSF makes up most of the voxel.

    Returns
    -------
    :obj:`dict`
        ``median`` and ``iqr`` (mL/100 g/min; NaN without voxels) and ``n`` (the voxels used).
    """
    sel = np.asarray(valid_csf, dtype=bool) & (np.asarray(pv_csf) >= min_pv)
    values = np.asarray(cbf_csf)[sel]
    if not values.size:
        return {'median': float('nan'), 'iqr': float('nan'), 'n': 0}
    q25, q50, q75 = np.percentile(values, [25, 50, 75])
    return {'median': float(q50), 'iqr': float(q75 - q25), 'n': int(values.size)}


def normalized_resample_ratio(data, weight, threshold=0.5):
    """Finish a validity-normalized resampling: divide the warped data by the warped weight.

    ``data`` is the warp of (value * valid) and ``weight`` the warp of ``valid``, with the same
    linear interpolation. Where ``weight`` is below ``threshold`` the result is zero.
    """
    data = np.asarray(data, dtype=np.float64)
    weight = np.asarray(weight, dtype=np.float64)
    out = np.zeros(data.shape)
    ok = np.isfinite(data) & np.isfinite(weight) & (weight >= threshold)
    out[ok] = data[ok] / weight[ok]
    return out
