"""Functions for calculating CBF."""
from __future__ import annotations

import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import median_abs_deviation

from aslprep import config


def _get_wfun_weight(x, wfun="huber"):
    """Get wavelet function weight."""
    weights = {
        "andrews": (np.abs(x) < np.pi) * np.sin(x),
        "bisquare": (np.abs(x) < 1) * np.power((1 - np.power(x, 2)), 2),
        "cauchy": 1 / (1 + np.power(x, 2)),
        "logistic": np.tanh(x) / x,
        "ols": np.repeat(1, len(x)),
        "talwar": 1 * (np.abs(x) < 1),
        "welsch": np.exp(-(np.power(x, 2))),
    }
    return weights.get(wfun, 1 / np.abs(x))


def _get_wfun_tuner(wfun="huber"):
    """Get weight fun and tuner.

    But wait, you might say, the docstring makes no sense! Correct.
    """
    weights = {
        "andrews": 1.339,
        "bisquare": 4.685,
        "cauchy": 2.385,
        "logistic": 1.205,
        "ols": 1,
        "talwar": 2.795,
        "welsch": 2.985,
    }
    return weights.get(wfun, 1.345)


def _getchisquare(n):
    """Get chi-square value for a given sample size.

    TODO: Replace with scipy.stats call, assuming this function is correct in any way.
    """
    a = [
        0.000000,
        15.484663,
        8.886835,
        7.224733,
        5.901333,
        5.126189,
        4.683238,
        4.272937,
        4.079918,
        3.731612,
        3.515615,
        3.459711,
        3.280471,
        3.078046,
        3.037280,
        2.990761,
        2.837119,
        2.795526,
        2.785189,
        2.649955,
        2.637642,
        2.532700,
        2.505253,
        2.469810,
        2.496135,
        2.342210,
        2.384975,
        2.275019,
        2.244482,
        2.249109,
        2.271968,
        2.210340,
        2.179537,
        2.133762,
        2.174928,
        2.150072,
        2.142526,
        2.071512,
        2.091061,
        2.039329,
        2.053183,
        2.066396,
        1.998564,
        1.993568,
        1.991905,
        1.981837,
        1.950225,
        1.938580,
        1.937753,
        1.882911,
        1.892665,
        1.960767,
        1.915530,
        1.847124,
        1.947374,
        1.872383,
        1.852023,
        1.861169,
        1.843109,
        1.823870,
        1.809643,
        1.815038,
        1.848064,
        1.791687,
        1.768343,
        1.778231,
        1.779046,
        1.759597,
        1.774383,
        1.774876,
        1.751232,
        1.755293,
        1.757028,
        1.751388,
        1.739384,
        1.716395,
        1.730631,
        1.718389,
        1.693839,
        1.696862,
        1.691245,
        1.682541,
        1.702515,
        1.700991,
        1.674607,
        1.669986,
        1.688864,
        1.653713,
        1.641309,
        1.648462,
        1.630380,
        1.634156,
        1.660821,
        1.625298,
        1.643779,
        1.631554,
        1.643987,
        1.624604,
        1.606314,
        1.609462,
    ]
    b = [
        0,
        2.177715,
        1.446966,
        1.272340,
        1.190646,
        1.151953,
        1.122953,
        1.103451,
        1.089395,
        1.079783,
        1.071751,
        1.063096,
        1.058524,
        1.054137,
        1.049783,
        1.046265,
        1.043192,
        1.039536,
        1.038500,
        1.037296,
        1.033765,
        1.032317,
        1.031334,
        1.029551,
        1.028829,
        1.027734,
        1.024896,
        1.024860,
        1.025207,
        1.024154,
        1.022032,
        1.021962,
        1.021514,
        1.020388,
        1.019238,
        1.020381,
        1.019068,
        1.018729,
        1.018395,
        1.017134,
        1.016539,
        1.015676,
        1.015641,
        1.015398,
        1.015481,
        1.015566,
        1.014620,
        1.014342,
        1.013901,
        1.013867,
        1.013838,
        1.013602,
        1.013322,
        1.012083,
        1.013168,
        1.012667,
        1.011087,
        1.011959,
        1.011670,
        1.011494,
        1.010463,
        1.010269,
        1.010393,
        1.010004,
        1.010775,
        1.009399,
        1.011000,
        1.010364,
        1.009831,
        1.009563,
        1.010085,
        1.009149,
        1.008444,
        1.009455,
        1.009705,
        1.008597,
        1.008644,
        1.008051,
        1.008085,
        1.008550,
        1.008265,
        1.009141,
        1.008235,
        1.008002,
        1.008007,
        1.007660,
        1.007993,
        1.007184,
        1.008093,
        1.007816,
        1.007770,
        1.007932,
        1.007819,
        1.007063,
        1.006712,
        1.006752,
        1.006703,
        1.006650,
        1.006743,
        1.007087,
    ]
    return a[n - 1], b[n - 1]


def _score_cbf(cbf_ts, wm, gm, csf, mask, thresh=0.7):
    """Apply SCORE algorithm to remove noisy CBF volumes.

    Parameters
    ----------
    cbf_ts : 4D numpy.ndarray
        CBF time series.
    gm, wm, csf : 3D numpy.ndarray
        Gray matter, white matter, and CSF tissue probability maps.
    mask : 3D numpy.ndarray
        Binary brain mask.
    thresh : float
        Tissue probability threshold.
        Voxels in gm, wm, and csf with values greater than or equal to the threshold are
        considered to be of that tissue.

    Returns
    -------
    cbf_ts_good_vols : 4D numpy.ndarray
        CBF time series containing only non-outlier volumes. This time series is also brain-masked.
    outlier_idx : 1D numpy.ndarray
        An index of outlier volumes in the CBF time series. Values of 1 indicate outliers.
    """
    n_volumes = cbf_ts.shape[3]

    # Create binary tissue-type masks
    gm_bin = gm >= thresh
    wm_bin = wm >= thresh
    csf_bin = csf >= thresh
    tissue_type_masks = [gm_bin, wm_bin, csf_bin]

    # Get the number of voxels within each tissue type, then subtract 1 to get Nk-1.
    n_gm_voxels = np.sum(gm_bin) - 1
    n_wm_voxels = np.sum(wm_bin) - 1
    n_csf_voxels = np.sum(csf_bin) - 1
    n_voxels = [n_gm_voxels, n_wm_voxels, n_csf_voxels]

    # Create a mask only including GM, WM, and CSF
    tissue_type_mask = (gm_bin + wm_bin + csf_bin).astype(bool)

    # Perform initial outlier exclusion.
    # Remove volumes with mean GM-CBF outside 2.5 SDs of their means.
    mean_gb_cbf = np.squeeze(np.mean(cbf_ts[gm_bin, :], axis=0))
    median_gm_cbf = np.median(mean_gb_cbf)
    mad_gm_cbf = median_abs_deviation(mean_gb_cbf) / 0.675
    outlier_idx = np.abs(mean_gb_cbf - median_gm_cbf) > (2.5 * mad_gm_cbf)

    iter_mean_cbf = np.mean(cbf_ts[:, :, :, ~outlier_idx], axis=3)
    n_voxels_total = np.sum(n_voxels)
    iter_pooled_variance = (
        np.sum([n_voxels[k] * np.var(iter_mean_cbf[tissue_type_masks[k]]) for k in range(3)])
        / n_voxels_total
    )
    last_pooled_variance = iter_pooled_variance + 1  # add 1 to ensure it is >
    while iter_pooled_variance < last_pooled_variance:
        # If variance increases from the previous iteration, we interpret this as removing
        # random noise instead of artifact.
        last_pooled_variance = iter_pooled_variance
        volume_wise_correlations = np.full(n_volumes, np.nan)
        for i_vol in range(n_volumes):
            if outlier_idx[i_vol]:
                continue

            cbf_vol = cbf_ts[:, :, :, i_vol]
            # Correlate mean CBF's masked values with volume CBF's masked values.
            volume_wise_correlations[i_vol] = np.corrcoef(
                iter_mean_cbf[tissue_type_mask],
                cbf_vol[tissue_type_mask],
            )[0, 1]

        # The volume with the maximum correlation value is this iteration's outlier.
        max_correlation_idx = np.nanargmax(volume_wise_correlations)
        outlier_idx[max_correlation_idx] = True

        # Recalculate mean CBF map and pooled variance.
        iter_mean_cbf = np.mean(cbf_ts[:, :, :, ~outlier_idx], axis=3)
        iter_pooled_variance = (
            np.sum([n_voxels[k] * np.var(iter_mean_cbf[tissue_type_masks[k]]) for k in range(3)])
            / n_voxels_total
        )

    # Create CBF time series from just the non-outlier volumes and brain-mask it.
    cbf_ts_good_vols = cbf_ts[:, :, :, ~outlier_idx]
    cbf_ts_good_vols *= mask[..., None]
    cbf_ts_good_vols = np.nan_to_num(cbf_ts_good_vols)

    return cbf_ts_good_vols, outlier_idx.astype(int)


def _robust_fit(
    Y,
    mu,
    global_prior,
    modrobprior,
    lmd=0,
    localprior=0,
    wfun="huber",
    tune=1.345,
    flagstd=1,
    flagmodrobust=1,
):
    """Perform robust fit, whatever that means."""
    dimcbf = Y.shape
    priow = np.ones([dimcbf[0], dimcbf[1]])
    sw = 1
    X = priow
    b = (np.sum(X * Y, axis=0) + mu * global_prior + lmd * localprior) / (
        np.sum(X * X, axis=0) + mu + lmd
    )
    b0 = np.repeat(0, len(b))
    h1 = X / np.power(np.tile(np.sqrt(np.sum(X * X, axis=0)), (dimcbf[0], 1)), 2)
    h0 = 0.9999 * np.ones([dimcbf[0], dimcbf[1]])
    h = np.minimum(h0, h1)
    adjfactor = 1 / (np.sqrt(1 - h / priow))
    tiny_s = (1e-6) * (np.std(h, axis=0))
    tiny_s[tiny_s == 0] = 1
    D = np.sqrt(np.finfo(float).eps)

    iter_num, max_iters = 0, 10
    while iter_num < max_iters:
        print(f"iteration {iter_num}")
        iter_num += 1
        check1 = np.subtract(np.abs(b - b0), (D * np.maximum(np.abs(b), np.abs(b0))))
        check1[check1 > 0] = 0
        if any(check1):
            print(f"converged after {iter_num} iterations")
            break

        r = Y - X * (np.tile(b, (dimcbf[0], 1)))
        radj = r * adjfactor / sw
        if flagstd == 1:
            s = np.sqrt(np.mean(np.power(radj, 2), axis=0))
        else:
            rs = np.sort(np.abs(radj), axis=0)
            s = np.median(rs, axis=0) / 0.6745

        rx1 = radj * (1 - flagmodrobust * np.exp(-np.tile(modrobprior, (dimcbf[0], 1))))
        rx2 = np.tile(np.maximum(s, tiny_s) * tune, (dimcbf[0], 1))
        r1 = rx1 / rx2
        w = _get_wfun_weight(r1, wfun)
        b0 = b
        z = np.sqrt(w)
        x = X * z
        yz = Y * z
        b = (np.sum(x * yz, axis=0) + mu * global_prior + lmd * localprior) / (
            np.sum(x * x, axis=0) + mu + lmd
        )
        b = np.nan_to_num(b)

    return b


def _scrub_cbf(cbf_ts, gm, wm, csf, mask, wavelet_function="huber", thresh=0.7):
    """Apply SCRUB algorithm to CBF data.

    Parameters
    ----------
    cbf_ts : 4D numpy.ndarray
        CBF time series. Generally comes after SCORE censoring.
    gm, wm, csf : 3D numpy.ndarray
        Gray matter, white matter, and CSF tissue probability maps.
    mask : 3D numpy.ndarray
        Binary brain mask.
    wavelet_function : str
        Wavelet function
    thresh : float
        Tissue probability threshold.
        Voxels in gm, wm, and csf with values greater than or equal to the threshold are
        considered to be of that tissue.

    Returns
    -------
    mean_cbf_scrub : 3D numpy.ndarray
        Mean CBF after denoising with SCRUB.
    """
    n_volumes = cbf_ts.shape[3]

    # Mask out non-brain voxels from tissue probability maps
    mask = mask.astype(bool)
    gm = gm * mask
    wm = wm * mask
    csf = csf * mask

    gm_bin = gm >= thresh
    gm_idx = gm_bin[mask]

    mean_cbf = np.mean(cbf_ts, axis=3)
    masked_cbf_ts = cbf_ts[mask, :].T  # TxS array
    masked_var_map = np.var(masked_cbf_ts, axis=0)  # S array
    thresh1, thresh3 = _getchisquare(n_volumes)
    mu1 = masked_var_map / (np.median(masked_var_map[gm_idx == 1]) * thresh3)
    mu = (
        ((mu1 > thresh1) & (mu1 < 10 * thresh1)) * (mu1 - thresh1)
        + (mu1 >= 10 * thresh1) * (1 / (2 * thresh1 * 10) * np.power(mu1, 2))
        + (thresh1 * 10 / 2 - thresh1)
    )
    # Create 3D array from mu
    mu_map = np.zeros_like(mean_cbf)
    mu_map[mask] = mu

    modrobprior = mu / 10
    gm_idx2 = (
        [gm.flatten() > thresh] and [mu_map.flatten() == 0] and [wm.flatten() > csf.flatten()]
    )[0].astype(int)
    wm_idx2 = (
        [wm.flatten() > thresh] and [mu_map.flatten() == 0] and [gm.flatten() > csf.flatten()]
    )[0].astype(int)
    if np.sum(gm_idx2) == 0 or np.sum(wm_idx2) == 0:
        gm_idx2 = (gm.flatten() > thresh).astype(int)
        wm_idx2 = (wm.flatten() > thresh).astype(int)

    idxx = gm_idx2 + wm_idx2
    idxx[idxx > 0] = 1
    X = np.zeros([len(idxx), 2])
    X[:, 0] = gm.flatten()[gm.flatten() >= (0)] * idxx
    X[:, 1] = wm.flatten()[wm.flatten() >= (0)] * idxx
    A = (mean_cbf.flatten()[idxx >= 0]) * idxx
    c = np.linalg.lstsq(X, A)[0]
    global_prior_full = c[0] * gm.flatten() + c[1] * wm.flatten()
    global_prior = global_prior_full[mask.flatten() == 1]

    tune = _get_wfun_tuner(wfun=wavelet_function)
    masked_mean_cbf_scrub = _robust_fit(
        Y=masked_cbf_ts,
        mu=mu,
        global_prior=global_prior,
        modrobprior=modrobprior,
        lmd=0,
        localprior=0,
        wfun=wavelet_function,
        tune=tune,
        flagstd=1,
        flagmodrobust=1,
    )
    mean_cbf_scrub = np.zeros_like(mean_cbf)
    mean_cbf_scrub[mask] = masked_mean_cbf_scrub
    mean_cbf_scrub = np.nan_to_num(mean_cbf_scrub)
    return mean_cbf_scrub


def estimate_att_pcasl(deltam_arr, plds, lds, t1blood, t1tissue):
    """Estimate arterial transit time using the weighted average method.

    The weighted average method comes from :footcite:t:`dai2012reduced`.

    Parameters
    ----------
    deltam_arr : :obj:`numpy.ndarray` of shape (S, D)
        Delta-M array, averaged by PLD.
    plds : :obj:`numpy.ndarray` of shape (S, D)
        Post-labeling delays. w in Dai 2012.
        In case of a 2D acquisition, PLDs may vary by slice, and thus the plds array will vary
        in the spatial dimension.
        For 3D acquisitions, or 2D acquisitions without slice timing info, plds will only vary
        along the second dimension.
    lds : :obj:`numpy.ndarray`
        Labeling durations. tau in Dai 2012.
    t1blood : :obj:`float`
        T1 relaxation rate for blood.
    t1tissue : :obj:`float`
        T1 relaxation rate for tissue.

    Returns
    -------
    att_arr : :obj:`numpy.ndarray`
        Arterial transit time array.

    Notes
    -----
    This function was originally written in MATLAB by Jianxun Qu and William Tackett.
    It was translated to Python by Taylor Salo.
    Taylor Salo modified the code to loop over voxels, in order to account for slice timing-shifted
    post-labeling delays.

    Please see https://shorturl.at/wCO56 and https://shorturl.at/aKQU3 for the original MATLAB
    code.

    The code could probably be improved by operating on arrays, rather than looping over voxels.
    It is also overkill for 3D acquisitions, where PLD doesn't vary by voxel.

    License
    -------
    MIT License

    Copyright (c) 2023 willtack

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    References
    ----------
    .. footbibliography::
    """
    n_voxels, n_plds = plds.shape
    assert deltam_arr.shape == plds.shape, f"{deltam_arr.shape} != {plds.shape}"

    # Beginning of auxil_asl_gen_wsum
    assert lds.size == n_plds, f"{lds.size} != {n_plds}"

    att_arr = np.empty(n_voxels)
    for i_voxel in range(n_voxels):
        deltam_by_voxel = deltam_arr[i_voxel, :]
        plds_by_voxel = plds[i_voxel, :]

        # Define the possible transit times to evaluate for this voxel
        transit_times = np.arange(
            np.round(np.min(plds_by_voxel), 3),
            np.round(np.max(plds_by_voxel), 3) + 0.001,
            0.001,
        )

        sig_sum = np.zeros((transit_times.size))
        sig_pld_sum = np.zeros((transit_times.size))

        for j_pld in range(n_plds):
            pld = plds_by_voxel[j_pld]
            ld = lds[j_pld]

            # e ^ (-delta / T1a)
            exp_tt_T1b = np.exp(-transit_times / t1blood)
            # e ^ (-max(w - delta, 0) / T1t)
            exp_pld_tt_T1t = np.exp(-np.maximum(0, pld - transit_times) / t1tissue)
            # e ^ (-max(tau + w - delta, 0) / T1t)
            exp_pld_ld_tt_T1t = np.exp(-np.maximum(0, pld + ld - transit_times) / t1tissue)

            # The combination of exponent terms in Equation 1
            sig = exp_tt_T1b * (exp_pld_tt_T1t - exp_pld_ld_tt_T1t)

            # The elements in Equation 1 that aren't touched here are:
            # 2M0t (2 * equilibrium magnetization of brain tissue),
            # Beta (a term to compensate for static tissue signal loss caused by vessel supp.
            # pulses),
            # alpha (labeling efficiency),
            # T1t (t1tissue),
            # f (CBF; perfusion rate)
            # lambda (partition coefficient)
            # It seems like they are cancelled out though, so it's probably fine.

            # Numerator in Equation 4, for a range of transit times
            sig_pld_sum += sig * pld
            # Denominator in Equation 4, for a range of transit times
            sig_sum += sig

        # Predicted weighted delay values for a range of transit times
        weighted_delay_predicted = sig_pld_sum / sig_sum  # TT
        # End of auxil_asl_gen_wsum

        # Calculate the observed weighted delay for each voxel
        weighted_delay_denom = np.sum(deltam_by_voxel)
        weighted_delay_num = np.sum(deltam_by_voxel * plds_by_voxel)
        weighted_delay_observed = weighted_delay_num / (
            np.abs(weighted_delay_denom) + np.finfo(float).eps
        )

        # Truncate extreme transit time value to the PLD limits
        weighted_delay_min = min(weighted_delay_predicted)
        weighted_delay_max = max(weighted_delay_predicted)
        weighted_delay_observed = np.maximum(weighted_delay_observed, weighted_delay_min)
        weighted_delay_observed = np.minimum(weighted_delay_observed, weighted_delay_max)

        # Use linear interpolation to get the ATT for each weighted delay value (i.e., each voxel),
        # using the predicted weighted delay and associated transit time arrays.
        interp_func = interp1d(weighted_delay_predicted, transit_times)

        att_arr[i_voxel] = interp_func(weighted_delay_observed)

    return att_arr


def estimate_cbf_pcasl_multipld(
    deltam_arr,
    scaled_m0data,
    plds,
    tau,
    labeleff,
    t1blood,
    t1tissue,
    unit_conversion,
    partition_coefficient,
):
    """Estimate CBF and ATT for multi-delay PCASL data.

    Parameters
    ----------
    deltam_arr : :obj:`numpy.ndarray` of shape (S, P)
        Control-label values for each voxel and PLDs.
        S = sample (i.e., voxel).
        P = Post-labeling delay (i.e., volume).
    scaled_m0data : :obj:`numpy.ndarray` of shape (S,)
        The M0 volume, after scaling based on the M0-scale value.
    plds : :obj:`numpy.ndarray` of shape (S, P)
        Post-labeling delays. One value for each volume in ``deltam_arr``.
    tau : :obj:`numpy.ndarray` of shape (P,) or (0,)
        Label duration. May be a single value or may vary across volumes/PLDs.
    labeleff : :obj:`float`
        Estimated labeling efficiency.
    t1blood : :obj:`float`
        The longitudinal relaxation time of blood in seconds.
    t1tissue : :obj:`float`
        The longitudinal relaxation time of tissue in seconds.
    unit_conversion : :obj:`float`
        The factor to convert CBF units from mL/g/s to mL/ (100 g)/min. 6000.
    partition_coefficient : :obj:`float`
        The brain/blood partition coefficient in mL/g. Called lambda in the literature.

    Returns
    -------
    att_arr : :obj:`numpy.ndarray` of shape (S,)
        Arterial transit time map.
    cbf : :obj:`numpy.ndarray` of shape (S,)
        Estimated cerebrospinal fluid map, after estimating for each PLD and averaging across
        PLDs.

    Notes
    -----
    Delta-M values are first averaged over time for each unique post-labeling delay value.

    Arterial transit time is estimated on a voxel-wise basis according to
    :footcite:t:`dai2012reduced`.

    CBF is then calculated for each unique PLD value using the mean delta-M values and the
    estimated ATT.

    CBF is then averaged over PLDs according to :footcite:t:`juttukonda2021characterizing`,
    in which an unweighted average is calculated for each voxel across all PLDs in which
    PLD + tau > ATT.

    References
    ----------
    .. footbibliography::
    """
    n_voxels, n_volumes = deltam_arr.shape
    n_voxels_pld, n_plds = plds.shape
    assert n_voxels_pld == n_voxels
    first_voxel_plds = plds[0, :]

    if n_plds != n_volumes:
        raise ValueError(
            f"Number of PostLabelingDelays ({n_plds}) does not match "
            f"number of delta-M volumes ({n_volumes})."
        )

    # Formula from Fan 2017 (equation 2)
    # Determine unique original post-labeling delays (ignoring slice timing shifts)
    unique_first_voxel_plds, unique_pld_idx = np.unique(first_voxel_plds, return_index=True)
    unique_plds = plds[:, unique_pld_idx]  # S x unique PLDs
    n_unique_plds = unique_pld_idx.size

    # tau should be a 1D array, with one volume per unique PLD
    if tau.size > 1:
        if tau.size != n_plds:
            raise ValueError(
                f"Number of LabelingDurations ({tau.size}) != "
                f"number of PostLabelingDelays ({n_plds})"
            )

        tau = tau[unique_pld_idx]
    else:
        tau = np.full(n_unique_plds, tau)

    mean_deltam_by_pld = np.zeros((n_voxels, n_unique_plds))
    for i_pld, first_voxel_pld in enumerate(unique_first_voxel_plds):
        pld_idx = first_voxel_plds == first_voxel_pld
        mean_deltam_by_pld[:, i_pld] = np.mean(deltam_arr[:, pld_idx], axis=1)

    # Estimate ATT for each voxel
    att_arr = estimate_att_pcasl(
        deltam_arr=mean_deltam_by_pld,
        plds=unique_plds,
        lds=tau,
        t1blood=t1blood,
        t1tissue=t1tissue,
    )

    # Start calculating CBF
    num_factor = unit_conversion * partition_coefficient
    denom_factor = 2 * labeleff * scaled_m0data * t1blood

    # Loop over PLDs and calculate CBF for each, accounting for ATT.
    cbf_by_pld = np.zeros((n_voxels, n_unique_plds))
    for i_pld in range(n_unique_plds):
        pld_by_voxel = unique_plds[:, i_pld]
        tau_for_pld = tau[i_pld]

        pld_num_factor = num_factor * mean_deltam_by_pld[:, i_pld] * np.exp(att_arr / t1blood)
        pld_denom_factor = denom_factor * (
            np.exp(-np.maximum(pld_by_voxel - att_arr, 0))
            - np.exp(-np.maximum(tau_for_pld + pld_by_voxel - att_arr, 0))
        )
        cbf_by_pld[:, i_pld] = pld_num_factor / pld_denom_factor

    # Average CBF across PLDs, but only include PLDs where PLD + tau > ATT for that voxel,
    # per Juttukonda 2021 (section 2.6).
    cbf = np.zeros(n_voxels)  # mean CBF
    for i_voxel in range(n_voxels):
        plds_voxel = unique_plds[i_voxel, :]
        cbf_by_pld_voxel = cbf_by_pld[i_voxel, :]
        att_voxel = att_arr[i_voxel]
        cbf[i_voxel] = np.mean(cbf_by_pld_voxel[(plds_voxel + tau) > att_voxel])

    return att_arr, cbf


def estimate_t1(metadata):
    """Estimate the relaxation rates of blood and gray matter based on magnetic field strength.

    t1blood is set based on the scanner's field strength,
    according to :footcite:t:`zhang2013vivo,alsop_recommended_2015`.
    If recommended values from these publications cannot be used
    (i.e., if the field strength isn't 1.5T, 3T, 7T),
    then the formula from :footcite:t:`zhang2013vivo` will be applied.

    t1tissue is set based on the scanner's field strength as well,
    according to :footcite:t:`wright2008water`.
    At the moment, field strengths other than 1.5T, 3T, and 7T are not supported and will
    raise an exception.

    Parameters
    ----------
    metadata : :obj:`dict`
        Dictionary of metadata from the ASL file.

    Returns
    -------
    t1blood : :obj:`float`
        Estimated relaxation rate of blood based on magnetic field strength.
    t1tissue : :obj:`float`
        Estimated relaxation rate of gray matter based on magnetic field strength.

    References
    ----------
    .. footbibliography::
    """
    T1BLOOD_DICT = {
        1.5: 1.35,
        3: 1.65,
        7: 2.087,
    }
    t1blood = T1BLOOD_DICT.get(metadata["MagneticFieldStrength"])
    if not t1blood:
        config.loggers.interface.warning(
            f"T1blood cannot be inferred for {metadata['MagneticFieldStrength']}T data. "
            "Defaulting to formula from Zhang et al. (2013)."
        )
        t1blood = (110 * metadata["MagneticFieldStrength"] + 1316) / 1000

    # TODO: Supplement with formula for other field strengths
    T1TISSUE_DICT = {
        1.5: 1.197,
        3: 1.607,
        7: 1.939,
    }
    t1tissue = T1TISSUE_DICT.get(metadata["MagneticFieldStrength"])
    if not t1tissue:
        raise ValueError(
            f"T1tissue cannot be inferred for {metadata['MagneticFieldStrength']}T data."
        )

    return t1blood, t1tissue
