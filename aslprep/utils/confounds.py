"""Functions for calculating and collecting confounds."""

import itertools
import os

import nibabel as nb
import numpy as np
import pandas as pd
from nibabel.processing import smooth_image
from nipype.interfaces.base import isdefined

# Default cbf thresholds (mL/100 g/min) for computing voxel percentages
CBF_THRESH_DEFAULTS = (100, 150, 200)


def _less_breakable(a_string):
    """Harden the string to different environments, whatever that means."""
    return ''.join(a_string.split()).strip('#')


def _camel_to_snake(name):
    """Convert a camelCase or PascalCase string to snake_case.

    An underscore is inserted before an uppercase letter that follows a lowercase letter
    or digit, or that starts a new word after an acronym (e.g., "StdDVARS" -> "std_dvars").
    """
    chars = []
    for i, char in enumerate(name):
        if char.isupper() and i > 0:
            prev_char = name[i - 1]
            next_char = name[i + 1] if i + 1 < len(name) else ''
            if (
                prev_char.islower()
                or prev_char.isdigit()
                or (prev_char.isupper() and next_char.islower())
            ):
                chars.append('_')
        chars.append(char.lower())
    return ''.join(chars)


def _adjust_indices(left_df, right_df):
    """Force missing values to appear at the beginning of the DataFrame instead of the end."""
    index_diff = left_df.shape[0] - right_df.shape[0]
    if index_diff > 0:
        # right_df is shorter
        empty_df = pd.DataFrame(
            np.full((np.abs(index_diff), right_df.shape[1]), np.nan),
            columns=right_df.columns,
        )
        right_df = pd.concat((empty_df, right_df), axis=0).reset_index(drop=True)
    elif index_diff < 0:
        # left_df is shorter
        empty_df = pd.DataFrame(
            np.full((np.abs(index_diff), left_df.shape[1]), np.nan),
            columns=left_df.columns,
        )
        left_df = pd.concat((empty_df, left_df), axis=0).reset_index(drop=True)

    return left_df, right_df


def _gather_confounds(
    signals=None,
    dvars=None,
    std_dvars=None,
    fdisp=None,
    rmsd=None,
    motion=None,
    newpath=None,
    score=None,
):
    """Load confounds from the filenames, concatenate together horizontally, and save new file.

    For some confounds (e.g., FD), the number of rows in the file will be one less than the
    number of volumes. This will be adjusted automatically in this function.
    """
    all_files = []
    confounds_list = []
    for confound, name in (
        (signals, 'Global signals'),
        (std_dvars, 'Standardized DVARS'),
        (dvars, 'DVARS'),
        (fdisp, 'Framewise displacement'),
        (rmsd, 'RMSD'),
        (motion, 'Motion parameters'),
        (score, 'score_outlier_index'),
    ):
        if confound is not None and isdefined(confound):
            confounds_list.append(name)
            if os.path.exists(confound) and os.stat(confound).st_size > 0:
                all_files.append(confound)

    confounds_data = pd.DataFrame()
    for file_name in all_files:  # assumes they all have headings already
        new = pd.read_table(file_name)
        for column_name in new.columns:
            new.rename(
                columns={column_name: _camel_to_snake(_less_breakable(column_name))}, inplace=True
            )

        confounds_data, new = _adjust_indices(confounds_data, new)
        confounds_data = pd.concat((confounds_data, new), axis=1)

    if newpath is None:
        newpath = os.getcwd()

    combined_out = os.path.join(newpath, 'confounds.tsv')
    confounds_data.to_csv(combined_out, sep='\t', index=False, na_rep='n/a')

    return combined_out, confounds_list


def dice(input1, input2):
    r"""Calculate Dice coefficient between two arrays.

    Computes the Dice coefficient (also known as Sorensen index) between two binary images.

    The metric is defined as

    .. math::

        DC=\frac{2|A\cap B|}{|A|+|B|}

    , where :math:`A` is the first and :math:`B` the second set of samples (here: binary objects).
    This method was first proposed in :footcite:t:`dice1945measures` and
    :footcite:t:`sorensen1948method`.

    Parameters
    ----------
    input1/input2 : :obj:`numpy.ndarray`
        Numpy arrays to compare.
        Can be any type but will be converted into binary:
        False where 0, True everywhere else.

    Returns
    -------
    coef : :obj:`float`
        The Dice coefficient between ``input1`` and ``input2``.
        It ranges from 0 (no overlap) to 1 (perfect overlap).

    References
    ----------
    .. footbibliography::
    """
    input1 = np.atleast_1d(input1.astype(bool))
    input2 = np.atleast_1d(input2.astype(bool))

    intersection = np.count_nonzero(input1 & input2)

    size_i1 = np.count_nonzero(input1)
    size_i2 = np.count_nonzero(input2)

    if (size_i1 + size_i2) == 0:
        coef = 0
    else:
        coef = (2 * intersection) / (size_i1 + size_i2)

    return coef


def pearson(input1, input2):
    """Calculate Pearson product moment correlation between two mask images.

    Parameters
    ----------
    input1/input2 : :obj:`numpy.ndarray`
        Numpy arrays to compare.
        Can be any type but will be converted into binary:
        False where 0, True everywhere else.

    Returns
    -------
    coef : :obj:`float`
        Correlation between the two images.
    """
    input1 = np.atleast_1d(input1.astype(bool)).flatten()
    input2 = np.atleast_1d(input2.astype(bool)).flatten()

    return np.corrcoef(input1, input2)[0, 1]


def overlap(input1, input2):
    r"""Calculate overlap coefficient between two images.

    The metric is defined as

    .. math::

        DC=\frac{|A \cap B||}{min(|A|,|B|)}

    , where :math:`A` is the first and :math:`B` the second set of samples (here: binary objects).

    The overlap coefficient is also known as the Szymkiewicz-Simpson coefficient
    :footcite:p:`vijaymeena2016survey`.

    Parameters
    ----------
    input1/input2 : :obj:`numpy.ndarray`
        Numpy arrays to compare.
        Can be any type but will be converted into binary:
        False where 0, True everywhere else.

    Returns
    -------
    coef : :obj:`float`
        Coverage between two images.

    References
    ----------
    .. footbibliography::
    """
    input1 = np.atleast_1d(input1.astype(bool))
    input2 = np.atleast_1d(input2.astype(bool))

    intersection = np.count_nonzero(input1 & input2)
    smallv = np.minimum(np.sum(input1), np.sum(input2))

    return intersection / smallv


def average_cbf_by_tissue(cbf, gm, wm, csf, thresh):
    """Compute mean GM, WM, and CSF CBF values.

    Parameters
    ----------
    cbf : str
        Path to CBF file.
    gm, wm, csf : str
        Paths to GM, WM, and CSF tissue probability maps, in same space and resolution as cbf.
    thresh : float
        Threshold to apply to the TPMs. Default is 0.7.

    Returns
    -------
    mean_tissue_cbfs : list of float
        Mean CBF values from binarized versions of the tissue maps.
    """
    cbf = nb.load(cbf).get_fdata()

    mean_tissue_cbfs = []
    for tpm in [gm, wm, csf]:
        tpm_data = nb.load(tpm).get_fdata()
        mean_tissue_cbf = np.mean(cbf[tpm_data >= thresh])
        mean_tissue_cbfs.append(mean_tissue_cbf)

    return mean_tissue_cbfs


def compute_cbf_threshold_stats(cbf, mask, thresholds=CBF_THRESH_DEFAULTS):
    """Compute the percentage of in-mask voxels that exceed specified CBF thresholds.

    Args:
        cbf: path or array of cbf data
        mask: path or array of binary brain mask
        thresholds: cutoff to check against (default: 100, 150, 200)

    Returns:
        a dictionary of percentages like {'perc_voxels_cbf_gt_100': 45.2}
        returns NaN if mask is empty
    """
    if isinstance(cbf, str):
        cbf_data = nb.load(cbf).get_fdata()
    else:
        cbf_data = np.array(cbf, dtype=float)

    if cbf_data.ndim > 3:
        cbf_data = np.nanmean(cbf_data, axis=3)

    if isinstance(mask, str):
        mask_data = nb.load(mask).get_fdata()
    else:
        mask_data = np.asarray(mask)

    mask_bool = mask_data > 0
    masked_cbf = cbf_data[mask_bool]

    # Drop NaN voxels so they dont skew the stats
    masked_cbf = masked_cbf[~np.isnan(masked_cbf)]
    n_voxels = masked_cbf.size

    stats = {}
    for thresh in thresholds:
        if n_voxels == 0:
            stats[f'perc_voxels_cbf_gt_{int(thresh)}'] = np.nan
        else:
            n_above = np.sum(masked_cbf > thresh)
            stats[f'perc_voxels_cbf_gt_{int(thresh)}'] = 100.0 * n_above / n_voxels

    return stats


def compute_clipped_voxel_percentage(in_file, mask):
    """Compute the percentage of in-mask voxels that reach the image's maximum intensity.

    When the scanner's reconstruction clips the signal at the top of its dynamic range
    (for example, 4095 for 12-bit data), many voxels share exactly the same, maximal value.
    In data that are not clipped, at most a handful of voxels have the maximal value,
    so this percentage is close to zero.

    Parameters
    ----------
    in_file : :obj:`str`
        Path to a 3D or 4D image, as acquired (not resampled).
        Interpolation would smooth over the clipped values.
    mask : :obj:`str`
        Path to a binary brain mask.
        It is resampled to the grid of ``in_file`` if the two grids differ.

    Returns
    -------
    :obj:`float`
        Percentage of in-mask voxels with the maximal value in at least one volume.
        NaN if the mask is empty or the image has no finite values.
    """
    from nilearn import image

    img = nb.load(in_file)
    mask_img = nb.load(mask)
    if mask_img.shape[:3] != img.shape[:3] or not np.allclose(mask_img.affine, img.affine):
        mask_img = image.resample_to_img(
            mask_img,
            image.index_img(img, 0) if img.ndim > 3 else img,
            interpolation='nearest',
            force_resample=True,
            copy_header=True,
        )

    data = img.get_fdata()
    if data.ndim == 3:
        data = data[..., None]

    mask_bool = mask_img.get_fdata() > 0
    n_voxels = int(mask_bool.sum())
    finite = np.isfinite(data)
    if n_voxels == 0 or not finite.any():
        return np.nan

    max_value = data[finite].max()
    clipped = np.any(data == max_value, axis=3)
    return 100.0 * np.sum(clipped & mask_bool) / n_voxels


#: Upper edges of the CSF partial volume bins of the background-suppression profile.
BS_CSF_BINS = np.linspace(0, 1, 11)
#: Signal-to-floor ratio below which a voxel counts as near the noise floor.
#: At twice the floor, a magnitude image keeps about 87% of the label-control difference.
BS_NEAR_FLOOR_RATIO = 2.0


def _resample_to(img, target, interpolation):
    """Resample ``img`` to the grid of the 3D or 4D image ``target`` when the grids differ."""
    from nilearn import image

    if img.shape[:3] == target.shape[:3] and np.allclose(img.affine, target.affine):
        return img
    return image.resample_to_img(
        img,
        image.index_img(target, 0) if target.ndim > 3 else target,
        interpolation=interpolation,
        force_resample=True,
        copy_header=True,
    )


def estimate_noise_floor(data, mask, zooms, margin=20.0):
    """Estimate the magnitude noise floor from the background of an image series.

    In a magnitude image, a voxel with no signal does not read zero: it reads the
    magnitude of the noise (Rayleigh-distributed for one receive channel and noncentral
    chi for a sum-of-squares or adaptive coil combination).
    For any of these, the mean squared magnitude is the squared true signal plus the
    squared floor, :math:`E[M^2] = A^2 + F^2`, where :math:`F^2` is the mean squared
    magnitude of pure noise.

    Parameters
    ----------
    data : :obj:`numpy.ndarray` of shape (X, Y, Z, T)
        Magnitude images, as acquired.
    mask : :obj:`numpy.ndarray` of shape (X, Y, Z)
        Boolean brain mask.
    zooms : :obj:`tuple` of :obj:`float`
        Voxel size in mm.
    margin : :obj:`float`
        Background voxels are farther than this many mm from the brain mask,
        which leaves out the skull and scalp.

    Returns
    -------
    floor : :obj:`float`
        Root mean squared magnitude of the background, in the units of ``data``.
        NaN if fewer than 100 background voxels have nonzero values in every volume
        (for example, when the scanner sets the background to zero).
    background : :obj:`numpy.ndarray` of shape (X, Y, Z)
        The background voxels used.
    """
    from scipy import ndimage

    distance = ndimage.distance_transform_edt(~mask, sampling=zooms)
    background = (distance > margin) & np.all(data > 0, axis=3) & np.all(np.isfinite(data), axis=3)
    if background.sum() < 100:
        return np.nan, background

    # The median over voxels ignores bright structures outside the brain (face, neck, ghosts).
    # Each voxel's mean over volumes is close to symmetric, so its median is close to its mean.
    mean_square = np.mean(data[background] ** 2, axis=1)
    return float(np.sqrt(np.median(mean_square))), background


def compute_background_suppression_qc(
    asl_file,
    aslcontext,
    mask,
    gm=None,
    csf=None,
    tpm_threshold=0.7,
    min_bin_voxels=20,
):
    """Measure how close background-suppressed control and label images are to zero.

    Background suppression lowers the static tissue signal at readout.
    Scanners usually reconstruct magnitude images, which have two consequences
    when the suppressed signal is small:

    1.  Near the noise floor, the magnitude overestimates the signal, and the
        label-control difference (delta-M) is compressed.
        A voxel whose mean squared control magnitude is :math:`M^2` keeps about
        :math:`\\sqrt{1 - F^2 / M^2}` of its delta-M, where :math:`F` is the noise floor.
    2.  If the static magnetization is inverted at readout, the magnitude flips its sign.
        Labeling lowers the signed magnetization, so the label image is then *brighter*
        than the control image and delta-M is negative.
        Tissue and CSF have different T1s, so they cross zero at different readout times.
        When one is inverted and the other is not, the net magnetization of voxels that
        mix them passes through zero, and the magnitude dips to the noise floor.

    Both are summarized along the tissue-to-CSF partial volume axis:
    the median signal-to-floor ratio and the median delta-M in bins of CSF fraction.

    Parameters
    ----------
    asl_file : :obj:`str`
        ASL series, as acquired (not resampled or motion corrected).
    aslcontext : :obj:`str`
        The ``aslcontext.tsv`` file of ``asl_file``.
    mask : :obj:`str`
        Binary brain mask. Resampled to ``asl_file``'s grid if needed.
    gm, csf : :obj:`str` or None
        Gray matter and CSF probability maps. Resampled to ``asl_file``'s grid if needed.
        Without them, only the noise floor and the near-floor percentage are computed.
    tpm_threshold : :obj:`float`
        Probability threshold for the gray matter mask.
    min_bin_voxels : :obj:`int`
        Bins of the CSF profile with fewer voxels are left out.

    Returns
    -------
    :obj:`dict`
        ``metrics`` (scalar QC measures), ``profile`` (a DataFrame with one row per CSF
        bin), and the arrays ``signal_floor_ratio`` (root mean squared control magnitude
        over the floor), ``deltam`` (mean control minus mean label magnitude), ``mask``,
        and ``background``, plus the ``affine`` of ``asl_file``.
        Measures that cannot be computed are NaN: all of them for series without both
        control and label volumes, and all but the CSF delta-M profile when there is no
        background to estimate the floor from.
    """
    img = nb.load(asl_file)
    context = pd.read_table(aslcontext)['volume_type'].to_numpy()
    mask_arr = _resample_to(nb.load(mask), img, 'nearest').get_fdata() > 0
    shape = img.shape[:3]

    metrics = dict.fromkeys(
        (
            'bs_noise_floor',
            'bs_gm_signal_floor_ratio',
            'bs_gm_deltam_retained',
            'bs_perc_voxels_near_noise_floor',
            'bs_csf_profile_signal_floor_ratio_min',
            'bs_csf_profile_deltam_min',
        ),
        np.nan,
    )
    out = {
        'metrics': metrics,
        'profile': pd.DataFrame(),
        'signal_floor_ratio': np.full(shape, np.nan),
        'deltam': np.full(shape, np.nan),
        'mask': mask_arr,
        'background': np.zeros(shape, dtype=bool),
        'affine': img.affine,
    }
    is_control = context == 'control'
    is_label = context == 'label'
    if img.ndim != 4 or not is_control.any() or not is_label.any() or not mask_arr.any():
        return out

    data = img.get_fdata()
    control = data[..., is_control]
    label = data[..., is_label]
    floor, background = estimate_noise_floor(
        data[..., is_control | is_label],
        mask_arr,
        img.header.get_zooms()[:3],
    )
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.sqrt(np.mean(control**2, axis=-1)) / floor
        retained = np.sqrt(np.clip(1 - 1 / ratio**2, 0, 1))
    deltam = control.mean(axis=-1) - label.mean(axis=-1)
    out.update(signal_floor_ratio=ratio, deltam=deltam, background=background)
    metrics['bs_noise_floor'] = floor
    if np.isfinite(floor):
        metrics['bs_perc_voxels_near_noise_floor'] = 100 * float(
            np.mean(ratio[mask_arr] < BS_NEAR_FLOOR_RATIO),
        )

    if gm is None or csf is None:
        return out

    gm_arr, csf_arr = (
        _resample_to(nb.load(tpm), img, 'continuous').get_fdata() for tpm in (gm, csf)
    )
    gm_mask = mask_arr & (gm_arr >= tpm_threshold)
    if gm_mask.any() and np.isfinite(floor):
        metrics['bs_gm_signal_floor_ratio'] = float(np.median(ratio[gm_mask]))
        metrics['bs_gm_deltam_retained'] = float(np.median(retained[gm_mask]))

    rows = []
    for low, high in itertools.pairwise(BS_CSF_BINS):
        in_bin = mask_arr & (csf_arr >= low) & ((csf_arr < high) | (high == BS_CSF_BINS[-1]))
        n_voxels = int(in_bin.sum())
        row = {'csf_low': low, 'csf_high': high, 'n_voxels': n_voxels}
        enough = n_voxels >= min_bin_voxels
        for name, arr in (('signal_floor_ratio', ratio), ('deltam', deltam)):
            q25, q50, q75 = np.percentile(arr[in_bin], [25, 50, 75]) if enough else [np.nan] * 3
            row.update({f'{name}_q25': q25, f'{name}_median': q50, f'{name}_q75': q75})
        rows.append(row)
    profile = pd.DataFrame(rows)
    out['profile'] = profile

    if profile['signal_floor_ratio_median'].notna().any():
        metrics['bs_csf_profile_signal_floor_ratio_min'] = float(
            profile['signal_floor_ratio_median'].min(),
        )
    dm = profile['deltam_median']
    scale = dm.abs().max()
    if dm.notna().any() and scale > 0:
        metrics['bs_csf_profile_deltam_min'] = float(dm.min() / scale)

    return out


def _load_one_image(nii_file):
    if isinstance(nii_file, str):
        img = nb.load(nii_file)
        data = img.get_fdata()
    else:
        data = nii_file
    if data.ndim > 3:
        return data[:, :, :, 0]
    return data


def structural_pseudocbf_correlation(gm_probseg, wm_probseg, cbf_image):
    """Compute structural pseudocbf (rho_ss) from :footcite:t:`dolui2024automated`.

    Parameters
    ----------
    gm_probseg, wm_probseg : str
        Paths to GM, WM tissue probability maps, in same space and resolution as cbf.
    cbf_image : str or :obj:`nibabel.nifti1.Nifti1Image`
        Path to CBF file or Nifti1Image.
        This image should be smoothed prior use in this function.

    Returns
    -------
    structural_pseudocbf_correlation : :obj:`float`
        Structural pseudocbf correlation.

    References
    ----------
    .. footbibliography
    """
    gm_probseg_data = _load_one_image(gm_probseg)
    wm_probseg_data = _load_one_image(wm_probseg)
    cbf_data = _load_one_image(cbf_image)

    structural_pseudocbf = 2.5 * gm_probseg_data + wm_probseg_data
    msk = (cbf_data != 0) & (~np.isnan(cbf_data)) & (~np.isnan(structural_pseudocbf))
    return np.clip(np.corrcoef(cbf_data[msk], structural_pseudocbf[msk])[1, 0], 0, None)


def dispersion_index(wm_mask, gm_mask, csf_mask, cbf_image):
    """Calculate the dispersion index (DI) from :footcite:t:`dolui2024automated`.

    Parameters
    ----------
    wm_mask, gm_mask, csf_mask : :obj:`numpy.ndarray`
        Binary masks for WM, GM, and CSF.
    cbf_image : str or :obj:`nibabel.nifti1.Nifti1Image`
        Path to CBF file or Nifti1Image.
        This image should be smoothed prior use in this function.

    Returns
    -------
    di : :obj:`float`
        Dispersion index.

    References
    ----------
    .. footbibliography
    """
    cbf_data = _load_one_image(cbf_image)

    n_gm = np.sum(gm_mask)
    n_wm = np.sum(wm_mask)
    n_csf = np.sum(csf_mask)
    # TODO: possibly add a check to ensure that n_gm, n_wm, n_csf are all greater than 1

    V = (
        (n_gm - 1) * np.var(cbf_data[gm_mask])
        + (n_wm - 1) * np.var(cbf_data[wm_mask])
        + (n_csf - 1) * np.var(cbf_data[csf_mask])
    ) / (n_gm + n_wm + n_csf - 3)

    mean_gm_cbf = np.mean(cbf_data[gm_mask])
    return V / np.abs(mean_gm_cbf)


def compute_qei(
    gm,
    wm,
    csf,
    img,
    thresh,
    alpha=-3.0126,
    beta=2.4419,
    gamma=0.054,
    delta=0.9272,
    epsilon=2.8478,
    zeta=0.5196,
):
    r"""Compute quality evaluation index (QEI) using the individual components.

    The QEI is based on :footcite:t:`dolui2024automated`.

    .. math::

        QEI = \sqrt[3]{(1 - e^{\alpha\rho_{ss}^{\beta}})
              e^{-\left(\gamma\text{DI}^{\delta} + \epsilon\text{negGM}^{\zeta}\right)}}

    The constants used here differ slightlyfrom those in the paper, but match the actual values
    used in the original QEI implementation.

    Parameters
    ----------
    gm, wm, csf : str
        Paths to GM, WM, and CSF tissue probability maps, in same space and resolution as cbf.
    img : str or :obj:`nibabel.nifti1.Nifti1Image`
        Path to CBF file or Nifti1Image.
        This image should NOT be smoothed prior use in this function.
    thresh : float
        Threshold to apply to the TPMs. Default is 0.7.

    Returns
    -------
    qei : float
        Quality evaluation index.
    rho_ss : float
        Structural pseudocbf correlation.
    di : float
        Dispersion index (DI).
    neg_gm : float
        Percentage of negative voxels within grey matter mask.

    References
    ----------
    .. footbibliography
    """
    smoothed_cbf = _load_one_image(smooth_image(nb.load(img), fwhm=5).get_fdata())

    gm_probseg_data = _load_one_image(gm)
    wm_probseg_data = _load_one_image(wm)
    csf_probseg_data = _load_one_image(csf)
    gm_mask = gm_probseg_data > thresh
    wm_mask = wm_probseg_data > thresh
    csf_mask = csf_probseg_data > thresh

    # Calculate the components of the QEI equation
    di = dispersion_index(wm_mask, gm_mask, csf_mask, smoothed_cbf)
    neg_gm = np.sum(smoothed_cbf[gm_mask] < 0) / (np.sum(gm_mask))
    rho_ss = structural_pseudocbf_correlation(gm_probseg_data, wm_probseg_data, smoothed_cbf)

    # evaluate the QEI equation
    part1 = 1 - np.exp(alpha * np.power(rho_ss, beta))
    part2 = np.exp(-(gamma * np.power(di, delta) + epsilon * np.power(neg_gm, zeta)))
    qei = np.power(part1 * part2, 1 / 3)
    return qei, rho_ss, di, neg_gm * 100
