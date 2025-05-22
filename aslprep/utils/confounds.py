"""Functions for calculating and collecting confounds."""

import os
import re

import nibabel as nb
import numpy as np
import pandas as pd
from nibabel.processing import smooth_image
from nipype.interfaces.base import isdefined


def _less_breakable(a_string):
    """Harden the string to different environments, whatever that means."""
    return ''.join(a_string.split()).strip('#')


def _camel_to_snake(name):
    """Convert camelCase string to snake_case.

    Taken from https://stackoverflow.com/questions/1175208/.
    If we end up using it more than just here, probably worth pulling in a well-tested package.
    """
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


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
    """Compute structural pseudocbf (rho_ss) from :footcite:t:`dolui2017automated`.

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
    msk = (cbf_data != 0) & (cbf_data != np.nan) & (structural_pseudocbf != np.nan)
    return np.clip(np.corrcoef(cbf_data[msk], structural_pseudocbf[msk])[1, 0], 0, None)


def dispersion_index(wm_mask, gm_mask, csf_mask, cbf_image):
    """Calculate the dispersion index (DI) from :footcite:t:`dolui2017automated`.

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

    The QEI is based on :footcite:t:`dolui2017automated`.

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
    cbf_data = _load_one_image(img)
    smoothed_cbf = smooth_image(cbf_data, fwhm=5)

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
