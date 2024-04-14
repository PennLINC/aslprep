"""Functions for calculating and collecting confounds."""

import os
import re

import nibabel as nb
import numpy as np
import pandas as pd
from nibabel.processing import smooth_image
from nipype.interfaces.base import isdefined
from scipy.stats import gmean


def _less_breakable(a_string):
    """Harden the string to different environments, whatever that means."""
    return "".join(a_string.split()).strip("#")


def _camel_to_snake(name):
    """Convert camelCase string to snake_case.

    Taken from https://stackoverflow.com/questions/1175208/.
    If we end up using it more than just here, probably worth pulling in a well-tested package.
    """
    s1 = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub("([a-z0-9])([A-Z])", r"\1_\2", s1).lower()


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
        (signals, "Global signals"),
        (std_dvars, "Standardized DVARS"),
        (dvars, "DVARS"),
        (fdisp, "Framewise displacement"),
        (rmsd, "RMSD"),
        (motion, "Motion parameters"),
        (score, "score_outlier_index"),
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

    combined_out = os.path.join(newpath, "confounds.tsv")
    confounds_data.to_csv(combined_out, sep="\t", index=False, na_rep="n/a")

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
    """Calculate Pearson product moment correlation between two images.

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


def compute_qei(gm, wm, csf, img, thresh):
    """Compute quality evaluation index (QEI) of CBF.

    The QEI is based on :footcite:t:`dolui2017automated`.

    References
    ----------
    .. footbibliography
    """

    def fun1(x, xdata):
        d1 = np.exp(-(x[0]) * np.power(xdata, x[1]))
        return d1

    def fun2(x, xdata):
        d1 = 1 - np.exp(-(x[0]) * np.power(xdata, x[1]))
        return d1

    x1 = [0.054, 0.9272]
    x2 = [2.8478, 0.5196]
    x4 = [3.0126, 2.4419]
    scbf = smooth_image(nb.load(img), fwhm=5).get_fdata()

    # Only use first volume for time series data.
    if scbf.ndim > 3:
        scbf = scbf[:, :, :, 0]

    # load prob maps
    gmm = nb.load(gm).get_fdata()
    wmm = nb.load(wm).get_fdata()
    ccf = nb.load(csf).get_fdata()
    if gmm.ndim > 3:
        gmm = gmm[:, :, :, 0]
        wmm = wmm[:, :, :, 0]
        ccf = ccf[:, :, :, 0]

    pbcf = 2.5 * gmm + wmm  # gmm is 2.5 times wm
    msk = np.array((scbf != 0) & (scbf != np.nan) & (pbcf != np.nan)).astype(int)

    gm1 = np.array(gmm > thresh)
    wm1 = np.array(wmm > thresh)
    cc1 = np.array(ccf > thresh)
    r1 = np.array([0, np.corrcoef(scbf[msk == 1], pbcf[msk == 1])[1, 0]]).max()

    V = (
        (np.sum(gm1) - 1) * np.var(scbf[gm1 > 0])
        + (np.sum(wm1) - 1) * np.var(scbf[wm1 > 0])
        + (np.sum(cc1) - 1) * np.var(scbf[cc1 > 0])
    ) / (np.sum(gm1 > 0) + np.sum(wm1 > 0) + np.sum(cc1 > 0) - 3)

    negGM = np.sum(scbf[gm1] < 0) / (np.sum(gm1))
    GMCBF = np.mean(scbf[gm1])
    CV = V / np.abs(GMCBF)
    Q = [fun1(x1, CV), fun1(x2, negGM), fun2(x4, r1)]
    return gmean(Q)


def negativevoxel(cbf, gm, thresh):
    """Compute percentage of negative voxels within grey matter mask."""
    gm = nb.load(gm).get_fdata()
    cbf = nb.load(cbf).get_fdata()
    gm_bin = gm > thresh

    n_gm_voxels = np.sum(gm_bin)

    negative_cbf_bin = cbf < 0
    n_negative_cbf_voxels_in_gm = np.sum(negative_cbf_bin * gm_bin)

    return 100 * (n_negative_cbf_voxels_in_gm / n_gm_voxels)
