"""Functions for calculating CBF."""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import median_abs_deviation

from aslprep import config


def _weightfun(x, wfun='huber'):
    """Get weight fun and tuner."""
    if wfun == 'andrews':
        tuner = 1.339
        weight = (np.abs(x) < np.pi) * np.sin(x)
    elif wfun == 'bisquare':
        tuner = 4.685
        weight = (np.abs(x) < 1) * np.power((1 - np.power(x, 2)), 2)
    elif wfun == 'cauchy':
        tuner = 2.385
        weight = 1 / (1 + np.power(x, 2))
    elif wfun == 'logistic':
        tuner = 1.205
        weight = np.tanh(x) / x
    elif wfun == 'ols':
        tuner = 1
        weight = np.repeat(1, len(x))
    elif wfun == 'talwar':
        tuner = 2.795
        weight = 1 * (np.abs(x) < 1)
    elif wfun == 'welsch':
        tuner = 2.985
        weight = np.exp(-(np.power(x, 2)))
    else:
        tuner = 1.345
        weight = 1 / np.abs(x)
    return weight, tuner


def _tune(wfun='huber'):
    """Get weight fun and tuner.

    But wait, you might say, the docstring makes no sense! Correct.
    """
    if wfun == 'andrews':
        tuner = 1.339
    elif wfun == 'bisquare':
        tuner = 4.685
    elif wfun == 'cauchy':
        tuner = 2.385
    elif wfun == 'logistic':
        tuner = 1.205
    elif wfun == 'ols':
        tuner = 1
    elif wfun == 'talwar':
        tuner = 2.795
    elif wfun == 'welsch':
        tuner = 2.985
    else:
        tuner = 1.345
    return tuner


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


def _getcbfscore(cbfts, wm, gm, csf, mask, thresh=0.7):
    """Apply SCORE algorithm to remove noisy CBF volumes.

    Parameters
    ----------
    cbf_ts
       n-dimensional array of 3D or 4D computed cbf
    gm,wm,csf
       numpy array of grey matter, whitematter, and csf
    mask
       numpy array of mask
    """
    gm_bin = gm >= thresh
    wm_bin = wm >= thresh
    csf_bin = csf >= thresh

    mask1 = gm_bin + wm_bin + csf_bin
    # Tissue classes with fewer than two voxels have no within-tissue variance,
    # so they are left out of the pooled variance (as in OSIPI's SCORE_denoising.m)
    tissue_masks = [tissue_bin for tissue_bin in (gm_bin, wm_bin, csf_bin) if tissue_bin.sum() > 1]

    # mean  of times series cbf within greymatter
    gm_cbf_ts = np.squeeze(np.mean(cbfts[gm_bin, :], axis=0))

    # robust mean and median
    median_gm_cbf = np.median(gm_cbf_ts)
    mad_gm_cbf = median_abs_deviation(gm_cbf_ts) / 0.675
    index = 1 * (np.abs(gm_cbf_ts - median_gm_cbf) > (2.5 * mad_gm_cbf))

    def _pooled_variance(R):
        # Pooled within-tissue variance of the mean CBF map
        return sum((tissue_bin.sum() - 1) * np.var(R[tissue_bin]) for tissue_bin in tissue_masks)

    # Iteratively remove the volume most correlated with the mean CBF map,
    # for as long as doing so lowers the pooled within-tissue variance.
    R = np.mean(cbfts[:, :, :, index == 0], axis=3)
    V = _pooled_variance(R)
    while np.sum(index == 0) > 2:
        CC = np.full(cbfts.shape[3], np.nan)
        for s in range(cbfts.shape[3]):
            if index[s] != 0:
                continue

            tmp1 = cbfts[:, :, :, s]
            with np.errstate(divide='ignore', invalid='ignore'):
                CC[s] = np.corrcoef(R[mask1 > 0], tmp1[mask1 > 0])[0][1]

        # A volume that is constant within the tissue mask has an undefined correlation,
        # so it is never the most correlated one (MATLAB's max skips NaNs too)
        if np.all(np.isnan(CC)):
            break

        inx = np.nanargmax(CC)
        index[inx] = 2
        R = np.mean(cbfts[:, :, :, index == 0], axis=3)
        V_new = _pooled_variance(R)
        if not V_new < V:
            # The removal did not lower the variance, so undo it
            index[inx] = 0
            break

        V = V_new

    config.loggers.utils.warning(f'SCORE retains {np.sum(index == 0)}/{index.size} volumes')
    cbfts_recon = cbfts[:, :, :, index == 0]
    cbfts_recon1 = np.zeros_like(cbfts_recon)
    for i in range(cbfts_recon.shape[3]):
        cbfts_recon1[:, :, :, i] = cbfts_recon[:, :, :, i] * mask

    cbfts_recon1 = np.nan_to_num(cbfts_recon1)
    return cbfts_recon1, index


def _robust_fit(
    Y,
    mu,
    Globalprior,
    modrobprior,
    lmd=0,
    localprior=0,
    wfun='huber',
    tune=1.345,
    flagstd=1,
    flagmodrobust=1,
):
    """Perform robust fit, whatever that means."""
    dimcbf = Y.shape
    priow = np.ones([dimcbf[0], dimcbf[1]])
    sw = 1
    X = priow
    b = (np.sum(X * Y, axis=0) + mu * Globalprior + lmd * localprior) / (
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
    iter_num = 0
    interlim = 10
    while iter_num < interlim:
        print('iteration  ', iter_num, '\n')
        iter_num = iter_num + 1
        check1 = np.subtract(np.abs(b - b0), (D * np.maximum(np.abs(b), np.abs(b0))))
        check1[check1 > 0] = 0
        if any(check1):
            print(' \n converged after ', iter_num, 'iterations\n')
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
        w, _ = _weightfun(r1, wfun)
        b0 = b
        z = np.sqrt(w)
        x = X * z
        yz = Y * z
        b = (np.sum(x * yz, axis=0) + mu * Globalprior + lmd * localprior) / (
            np.sum(x * x, axis=0) + mu + lmd
        )
        b = np.nan_to_num(b)
    return b


def _scrubcbf(cbf_ts, gm, wm, csf, mask, wfun='huber', thresh=0.7):
    """Apply SCRUB algorithm to CBF data.

    Parameters
    ----------
    cbf_ts
       n-dimensional array of 3D or 4D computed cbf
    gm,wm,csf
       numpy array of grey matter, whitematter, and csf
    mask
       numpy array of mask
    wf
      wave function
    """
    gm = mask * gm
    gmidx = gm[mask == 1]
    gmidx[gmidx < thresh] = 0
    gmidx[gmidx > 0] = 1

    wm = mask * wm
    wmidx = wm[mask == 1]
    wmidx[wmidx < thresh] = 0
    wmidx[wmidx > 0] = 1

    csf = csf * mask
    csfidx = csf[mask == 1]
    csfidx[csfidx < thresh] = 0
    csfidx[csfidx > 0] = 1
    # midx = mask[mask==1]

    mean_cbf = np.mean(cbf_ts, axis=3)
    y = np.transpose(
        cbf_ts[
            mask == 1,
            :,
        ]
    )
    VV = np.var(y, axis=0)
    thresh1, thresh3 = _getchisquare(y.shape[0])
    mu1 = VV / (np.median(VV[gmidx == 1]) * thresh3)
    mu = (
        ((mu1 > thresh1) & (mu1 < 10 * thresh1)) * (mu1 - thresh1)
        + (mu1 >= 10 * thresh1) * (1 / (2 * thresh1 * 10) * np.power(mu1, 2))
        + (thresh1 * 10 / 2 - thresh1)
    )
    M = mean_cbf * mask
    M[mask == 1] = mu
    modrobprior = mu / 10
    gmidx2 = 1 * ((gm.flatten() > thresh) & (M.flatten() == 0) & (wm.flatten() > csf.flatten()))
    wmidx2 = 1 * ((wm.flatten() > thresh) & (M.flatten() == 0) & (gm.flatten() > csf.flatten()))
    if np.sum(gmidx2) == 0 or np.sum(wmidx2) == 0:
        gmidx2 = 1 * (gm.flatten() > thresh)
        wmidx2 = 1 * (wm.flatten() > thresh)
    idxx = gmidx2 + wmidx2
    idxx[idxx > 0] = 1
    X = np.zeros([len(idxx), 2])
    X[:, 0] = gm.flatten()[gm.flatten() >= (0)] * idxx
    X[:, 1] = wm.flatten()[wm.flatten() >= (0)] * idxx
    A = (mean_cbf.flatten()[idxx >= 0]) * idxx
    c = np.linalg.lstsq(X, A, rcond=None)[0]
    Globalpriorfull = c[0] * gm.flatten() + c[1] * wm.flatten()
    Globalprior = Globalpriorfull[mask.flatten() == 1]
    localprior = 0
    lmd = 0
    tune = _tune(wfun=wfun)
    bb = _robust_fit(
        Y=y,
        mu=mu,
        Globalprior=Globalprior,
        modrobprior=modrobprior,
        lmd=lmd,
        localprior=localprior,
        wfun=wfun,
        tune=tune,
        flagstd=1,
        flagmodrobust=1,
    )
    newcbf = mean_cbf * mask
    newcbf[mask == 1] = bb
    newcbf = np.nan_to_num(newcbf)
    return newcbf


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
    if deltam_arr.shape != plds.shape:
        raise ValueError(f'{deltam_arr.shape} != {plds.shape}')

    # Beginning of auxil_asl_gen_wsum
    if lds.size != n_plds:
        raise ValueError(f'{lds.size} != {n_plds}')

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

        sig_sum = np.zeros(transit_times.size)
        sig_pld_sum = np.zeros(transit_times.size)

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
        The blood-brain partition coefficient in mL/g. Called lambda in the literature.

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
    if n_voxels_pld != n_voxels:
        raise ValueError(f'{n_voxels_pld} != {n_voxels}')

    first_voxel_plds = plds[0, :]

    if n_plds != n_volumes:
        raise ValueError(
            f'Number of PostLabelingDelays ({n_plds}) does not match '
            f'number of delta-M volumes ({n_volumes}).'
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
                f'Number of LabelingDurations ({tau.size}) != '
                f'number of PostLabelingDelays ({n_plds})'
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
            np.exp(-(np.maximum(pld_by_voxel - att_arr, 0) / t1tissue))
            - np.exp(-(np.maximum(tau_for_pld + pld_by_voxel - att_arr, 0) / t1tissue))
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


def calculate_deltam_pcasl(X, cbf, att, abat, abv):
    """Calculate delta-M for one voxel of multi-PLD (P)CASL data with the general kinetic model.

    This scalar implementation is the reference for the vectorized model used by
    :func:`fit_deltam_multipld`.

    Parameters
    ----------
    cbf : :obj:`float`
        Cerebral blood flow.
        Used for deltam_tiss calculation.
        Estimated by the model.
    att : :obj:`float`
        Arterial transit time.
        Used for deltam_tiss calculation.
        Estimated by the model.
    abat : :obj:`float`
        The arterial bolus arrival time, describing the first arrival of the labeled blood
        water in the voxel within the arterial vessel. Used for deltam_art calculation.
        Estimated by the model.
    abv : :obj:`float`
        Arterial blood volume.
        The fraction of the voxel occupied by the arterial vessel.
        Used for deltam_art calculation.
        Estimated by the model.
    X : :obj:`numpy.ndarray` of shape (n_plds, 6)
        Dependent variables: labeleff, t1blood, m0_a, m0_b, plds, taus.

    Returns
    -------
    deltam : :obj:`numpy.ndarray` of shape (n_plds,)
        Delta-M values for the voxel.

    Notes
    -----
    X breaks down into six variables:

    labeleff : :obj:`float`
        Labeling efficiency.
        If background suppression was performed, then this is already adjusted and combines
        alpha_bs as well.
        Used for both deltam_tiss and deltam_art calculation.
    t1blood : :obj:`float`
        Longitudinal relaxation time of arterial blood in seconds.
        Used for both deltam_tiss and deltam_art calculation.
    m0_a : :obj:`float`
        Equilibrium magnetization of arterial blood (M0b in the OSIPI lexicon),
        calculated as M0a = SIpd / lambda,
        where SIpd is a proton density weighted image (the equilibrium magnetization of tissue)
        and lambda is the blood-brain partition coefficient.
        Used for deltam_tiss calculation.
    m0_b : :obj:`float`
        Equilibrium magnetization of arterial blood.
        ASLPrep uses the same value as ``m0_a``.
        Used for deltam_art calculation.
    plds : :obj:`numpy.ndarray`
        Post-labeling delays.
        Used for both deltam_tiss and deltam_art calculation.
    taus : :obj:`numpy.ndarray`
        Labeling durations. Also known as LD.
        Used for both deltam_tiss and deltam_art calculation.
    """
    # float parameters
    labeleff = X[0, 0]
    t1blood = X[0, 1]
    m0_a = X[0, 2]
    m0_b = X[0, 3]
    # array parameters
    plds = X[:, 4]
    taus = X[:, 5]

    deltam = np.zeros(plds.size)

    for i_pld, pld in enumerate(plds):
        tau = taus[i_pld]

        # Equation 2: the tissue component
        if 0 < (tau + pld) < att:
            # 0 < LD + PLD < ATT
            deltam_tiss = 0
        elif att <= (tau + pld) < (att + tau):
            # ATT < LD + PLD < ATT + LD
            deltam_tiss = (
                2
                * labeleff
                * t1blood
                * m0_a
                * cbf
                * np.exp(-att / t1blood)
                * (1 - np.exp(-(tau + pld - att) / t1blood))
            ) / 6000
        else:
            # 0 < ATT < PLD
            deltam_tiss = (
                2
                * labeleff
                * t1blood
                * m0_a
                * cbf
                * np.exp(-pld / t1blood)
                * (1 - np.exp(-tau / t1blood))
            ) / 6000

        # Equation 4: the intravascular component
        if 0 < (tau + pld) < abat:
            # 0 < LD + PLD < aBAT
            deltam_art = 0
        elif abat <= (tau + pld) < (abat + tau):
            # aBAT < LD + PLD < aBAT + LD
            deltam_art = 2 * labeleff * m0_b * abv * np.exp(-abat / t1blood)
        else:
            # aBAT < PLD
            deltam_art = 0

        deltam[i_pld] = deltam_tiss + deltam_art

    return deltam


def calculate_deltam_pasl(X, cbf, att, abat, abv):
    """Calculate delta-M for one voxel of multi-TI PASL data with the general kinetic model.

    This scalar implementation is the reference for the vectorized model used by
    :func:`fit_deltam_multipld`.

    Parameters
    ----------
    cbf : :obj:`float`
        Cerebral blood flow.
        Used for deltam_tiss calculation.
        Estimated by the model.
    att : :obj:`float`
        Arterial transit time.
        Used for deltam_tiss calculation.
        Estimated by the model.
    abat : :obj:`float`
        The arterial bolus arrival time, describing the first arrival of the labeled blood
        water in the voxel within the arterial vessel. Used for deltam_art calculation.
        Estimated by the model.
    abv : :obj:`float`
        Arterial blood volume.
        The fraction of the voxel occupied by the arterial vessel.
        Used for deltam_art calculation.
        Estimated by the model.
    X : :obj:`numpy.ndarray` of shape (n_plds/tis, 6)
        Dependent variables: labeleff, t1blood, m0_a, m0_b, tis, ti1s.

    Returns
    -------
    deltam : :obj:`numpy.ndarray` of shape (n_tis,)
        Delta-M values for the voxel.

    Notes
    -----
    X breaks down into six variables:

    labeleff : :obj:`float`
        Labeling efficiency.
        In the case of background suppression, this combines alpha and alpha_bs,
        so we don't need a separate term for alpha_bs.
        Used for both deltam_tiss and deltam_art calculation.
    t1blood : :obj:`float`
        Longitudinal relaxation time of arterial blood in seconds.
        Used for both deltam_tiss and deltam_art calculation.
    m0_a : :obj:`float`
        Equilibrium magnetization of arterial blood (M0b in the OSIPI lexicon),
        calculated as M0a = SIpd / lambda,
        where SIpd is a proton density weighted image (the equilibrium magnetization of tissue)
        and lambda is the blood-brain partition coefficient.
        SIpd is also already scaled by any user-provided scaling factor.
        Used for deltam_tiss calculation.
    m0_b : :obj:`float`
        Equilibrium magnetization of arterial blood.
        ASLPrep uses the same value as ``m0_a``.
        Used for deltam_art calculation.
    tis : :obj:`numpy.ndarray`
        Inversion times in seconds.
    ti1 : :obj:`float`
        Bolus cutoff delay time. For Q2TIPS, this is the first delay time.
        For QUIPSSII, it's the only one.
    """
    # float parameters
    labeleff = X[0, 0]
    t1blood = X[0, 1]
    m0_a = X[0, 2]
    m0_b = X[0, 3]
    tis = X[:, 4]  # array parameter
    ti1 = X[0, 5]

    deltam = np.zeros(tis.size)

    for i_ti, ti in enumerate(tis):
        # Equation 3: the tissue component
        if 0 < ti < att:
            # 0 < TI < ATT
            deltam_tiss = 0
        elif att <= ti < (att + ti1):
            # ATT < TI < ATT + TI1
            deltam_tiss = (2 * labeleff * m0_a * cbf * np.exp(-ti / t1blood) * (ti - att)) / 6000
        else:
            # ATT + TI1 < TI
            deltam_tiss = (2 * labeleff * m0_a * cbf * np.exp(-ti / t1blood) * ti1) / 6000

        # Equation 5: the intravascular component
        if 0 < ti < abat:
            # 0 < TI < aBAT
            deltam_art = 0
        elif abat < ti < (abat + ti1):
            # aBAT < TI < aBAT + TI1
            deltam_art = 2 * labeleff * m0_b * abv * np.exp(-ti / t1blood)
        else:
            # aBAT + TI1 < TI
            deltam_art = 0

        deltam[i_ti] = deltam_tiss + deltam_art

    return deltam


#: Lower and upper bounds for (CBF, ATT, aBAT, aBV) in the multi-delay general kinetic model.
#: Manuel Taso recommended 5, 5, and 0.2 for the ATT, aBAT, and aBV upper bounds.
MULTIPLD_BOUNDS = ((0.0, 300.0), (0.0, 5.0), (0.0, 5.0), (0.0, 0.1))

#: Number of evenly spaced ATT values evaluated in each interval on each search pass.
MULTIPLD_ATT_GRID_SIZE = 11

#: Number of ATT search passes. Each pass narrows the search to the neighbors of the previous
#: pass's best ATT, so the final ATT spacing is the interval width / (10 * 5 ** 5).
MULTIPLD_ATT_PASSES = 6

#: Number of (ATT interval, aBAT interval) combinations, ranked by the first pass,
#: that later passes refine in each voxel.
MULTIPLD_REFINED_INTERVALS = 8

#: Number of bisection steps used to locate each end of the range of equally good ATTs.
#: The search starts from an interval no wider than one first-pass grid step.
MULTIPLD_BOUND_STEPS = 30

#: Number of voxels fit at once with NumPy.
MULTIPLD_CHUNK_SIZE = 256

#: Number of elements in each of the largest working arrays when fitting on a GPU.
#: The fit holds a few dozen such arrays of doubles at once, so this uses 1-2 GB of GPU memory.
MULTIPLD_GPU_CHUNK_ELEMENTS = 15_000_000


def _gkm_tissue_basis(att, plds, bolus, is_casl, labeleff, t1blood, xp):
    """Calculate the tissue delta-M per unit CBF and per unit M0a.

    This is a vectorized version of the tissue terms in :func:`calculate_deltam_pcasl`
    (where ``bolus`` is the labeling duration) and :func:`calculate_deltam_pasl`
    (where ``plds`` are inversion times and ``bolus`` is TI1).
    ``att`` must broadcast against ``plds`` and ``bolus``.
    """
    if is_casl:
        end = bolus + plds
        partial = xp.exp(-att / t1blood) * (1 - xp.exp(-(end - att) / t1blood))
        full = xp.exp(-plds / t1blood) * (1 - xp.exp(-bolus / t1blood))
        basis = xp.where(
            (end > 0) & (end < att),
            0.0,
            xp.where((att <= end) & (end < att + bolus), partial, full),
        )
        return basis * (2 * labeleff * t1blood / 6000)

    basis = xp.where(
        (plds > 0) & (plds < att),
        0.0,
        xp.where((att <= plds) & (plds < att + bolus), plds - att, bolus),
    )
    return basis * xp.exp(-plds / t1blood) * (2 * labeleff / 6000)


def _gkm_arterial_basis(abat, plds, bolus, is_casl, labeleff, t1blood, xp):
    """Calculate the arterial delta-M per unit aBV and per unit M0b.

    This is a vectorized version of the intravascular terms in :func:`calculate_deltam_pcasl`
    and :func:`calculate_deltam_pasl`. ``abat`` must broadcast against ``plds`` and ``bolus``.
    """
    if is_casl:
        end = bolus + plds
        active = (abat <= end) & (end < abat + bolus)
        return xp.where(active, 2 * labeleff * xp.exp(-abat / t1blood), 0.0)

    active = (abat < plds) & (plds < abat + bolus)
    return xp.where(active, 2 * labeleff * xp.exp(-plds / t1blood), 0.0)


def _delay_edges(plds, bolus, is_casl, xp):
    """Find the parameter values at which the kinetic model's terms change form.

    The tissue term changes form when ATT, and the arterial term when aBAT, crosses one of these
    values: each PLD and PLD + labeling duration for (P)CASL, or each TI - TI1 and TI for PASL.
    Between consecutive edges, the model is smooth in ATT,
    and aBAT changes the arterial signal only by a common scale factor for (P)CASL,
    which aBV absorbs, and not at all for PASL.

    Returns the sorted edges, together with the parameter bounds, as an array of shape
    (n_voxels, n_edges). The ATT and aBAT bounds are the same.
    """
    lower, upper = MULTIPLD_BOUNDS[1]
    if is_casl:
        edges = [plds, plds + bolus]
    else:
        edges = [plds - bolus, plds]

    n_voxels = plds.shape[0]
    edges = xp.concatenate(
        [xp.broadcast_to(edge, plds.shape) for edge in edges]
        + [xp.full((n_voxels, 1), lower), xp.full((n_voxels, 1), upper)],
        axis=1,
    )
    return xp.sort(xp.clip(edges, lower, upper), axis=1)


def _solve_box_lsq(yt, ya, tt, aa, ta, tol, xp, v_hi=None):
    """Solve the bounded two-parameter linear least squares problem for CBF and aBV.

    With ``c`` as CBF, ``v`` as aBV, ``t`` as the tissue basis, ``a`` as the arterial basis,
    and ``y`` as the data, the inputs are the inner products of those vectors.
    The sum of squared errors, minus the constant ``y.y``, is minimized over the CBF and aBV
    bounds. The objective is convex, so the minimum is either the unconstrained solution,
    if it is in bounds, or the best bounded one-parameter solution on one of the box edges.
    A candidate only replaces an earlier one if its loss is lower by more than ``tol``,
    so that ties are not broken by rounding error.
    ``v_hi`` replaces the upper bound on aBV, and may be an array.
    """
    (c_lo, c_hi), _, _, (v_lo, default_v_hi) = MULTIPLD_BOUNDS
    v_hi = default_v_hi if v_hi is None else v_hi

    def loss(c, v):
        return c * (c * tt - 2 * yt) + v * (v * aa - 2 * ya) + 2 * c * v * ta

    det = tt * aa - ta * ta
    interior = det > 1e-10 * tt * aa
    det = xp.where(interior, det, 1.0)
    c = (aa * yt - ta * ya) / det
    v = (tt * ya - ta * yt) / det
    interior &= (c >= c_lo) & (c <= c_hi) & (v >= v_lo) & (v <= v_hi)
    best_c = xp.where(interior, c, c_lo)
    best_v = xp.where(interior, v, v_lo)
    best_loss = xp.where(interior, loss(c, v), xp.inf)

    def keep_better(c, v):
        nonlocal best_c, best_v, best_loss
        candidate_loss = loss(c, v)
        better = candidate_loss < best_loss - tol
        best_c = xp.where(better, c, best_c)
        best_v = xp.where(better, v, best_v)
        best_loss = xp.where(better, candidate_loss, best_loss)

    has_tissue = tt > 0
    has_arterial = aa > 0
    safe_tt = xp.where(has_tissue, tt, 1.0)
    safe_aa = xp.where(has_arterial, aa, 1.0)
    for v_edge in (v_lo, v_hi):
        c = xp.where(has_tissue, xp.clip((yt - v_edge * ta) / safe_tt, c_lo, c_hi), c_lo)
        keep_better(c, v_edge)

    for c_edge in (c_lo, c_hi):
        v = (ya - c_edge * ta) / safe_aa
        v = xp.where(has_arterial, xp.minimum(xp.maximum(v, v_lo), v_hi), v_lo)
        keep_better(c_edge, v)

    return best_c, best_v, best_loss


def _tie_rank(loss, tol, axis, xp):
    """Rank losses along an axis, treating losses within ``tol`` of each other as equal.

    Using these ranks with :func:`numpy.argmin` or a stable sort selects the first of several
    equally good candidates, rather than whichever one rounding error favors,
    so that the result does not depend on the array library or chunk size.
    """
    return xp.floor((loss - loss.min(axis=axis, keepdims=True)) / tol)


def _argmin_att(loss, tol, xp):
    """Find the best ATT along the last axis, at the middle of a run of equally good ATTs.

    Where a range of ATTs fits equally well, the data do not locate ATT within that range,
    so the middle of the first run of grid points tied for the minimum is returned.
    Without ties, this is the minimum.
    """
    tied = _tie_rank(loss, tol, -1, xp) == 0
    first = xp.argmax(tied, axis=-1)[..., None]
    # Grid points belong to the first run if no untied point separates them from its start.
    breaks = xp.cumsum(~tied, axis=-1)
    in_run = tied & (breaks == xp.take_along_axis(breaks, first, axis=-1))
    return first + (in_run.sum(axis=-1, keepdims=True) - 1) // 2


def _fit_gkm_chunk(y, plds, bolus, weights, is_casl, labeleff, t1blood, xp):
    """Fit the general kinetic model to a chunk of voxels by searching over ATT and aBAT.

    Parameters
    ----------
    y : array of shape (n_voxels, n_delays)
        Delta-M divided by M0a, then multiplied by ``weights``.
    plds : array of shape (n_voxels, n_delays)
        Post-labeling delays (PCASL) or inversion times (PASL).
    bolus : array of shape (1, n_delays) or float
        Labeling durations (PCASL) or TI1 (PASL).
    weights : array of shape (1, n_delays)
        Square root of the number of volumes averaged for each delay.

    Notes
    -----
    The search treats each combination of an ATT interval and an aBAT interval
    (see :func:`_delay_edges`) separately.
    aBAT is fixed at its interval's midpoint, because the data cannot localize it further.
    ATT is searched with a grid that spans its interval, including the edges.
    Where a range of ATTs fits equally well, the middle of that range is selected.
    The first pass evaluates every combination,
    and later passes refine the best :data:`MULTIPLD_REFINED_INTERVALS` combinations,
    each narrowing the grid to the neighbors of the previous pass's best ATT.

    aBV is solved as if aBAT were at its interval's midpoint.
    For (P)CASL, a shorter aBAT in the same interval gives more arterial signal per unit aBV,
    so the aBV bound is relaxed to allow the signal that the true bound allows at the interval's
    lower edge. aBAT is then reported as the midpoint of the part of the interval in which the
    fitted arterial signal respects the aBV bound, which is the whole interval unless aBV is
    close to its bound.

    Returns
    -------
    cbf, att, abat, abv : arrays of shape (n_voxels,)
        The estimates.
    att_lower, att_upper : arrays of shape (n_voxels,)
        The ends of the range of ATTs around ``att`` that fit as well as it does,
        with aBAT in the same interval, and with CBF and aBV free.
    abat_lower, abat_upper : arrays of shape (n_voxels,)
        The ends of the range of aBATs around ``abat`` that fit as well as it does,
        with ATT fixed at ``att``, and with CBF and aBV free.
    """
    n_voxels = y.shape[0]
    rows = xp.arange(n_voxels)
    basis_kwargs = {'is_casl': is_casl, 'labeleff': labeleff, 't1blood': t1blood, 'xp': xp}
    grid = xp.linspace(0, 1, MULTIPLD_ATT_GRID_SIZE)
    # Losses closer than this are treated as ties. It is well above the rounding error of the
    # losses (about 1e-15 * y.y), and the constant covers all-zero data.
    tol = 1e-12 * xp.sum(y * y, axis=1) + 1e-30

    edges = _delay_edges(plds, bolus, is_casl, xp)
    n_intervals = edges.shape[1] - 1

    # Arterial terms at each aBAT interval's midpoint: (voxels, aBATs, delays)
    abat_mid = (edges[:, :-1] + edges[:, 1:]) / 2
    all_arterial = _gkm_arterial_basis(
        abat_mid[:, :, None],
        plds[:, None, :],
        bolus[:, None, :] if is_casl else bolus,
        **basis_kwargs,
    )
    all_arterial *= weights
    # The largest midpoint-equivalent aBV in each aBAT interval: (voxels, aBATs)
    abv_max = MULTIPLD_BOUNDS[3][1]
    if is_casl:
        all_abv_cap = abv_max * xp.exp((abat_mid - edges[:, :-1]) / t1blood)
    else:
        all_abv_cap = xp.full(abat_mid.shape, abv_max)

    # First pass over every combination, with shape (voxels, ATTs, aBATs, ATT grid).
    # The tissue terms only depend on the ATT interval.
    att = edges[:, :-1, None] + (edges[:, 1:] - edges[:, :-1])[:, :, None] * grid
    tissue = _gkm_tissue_basis(
        att[..., None],
        plds[:, None, None, :],
        bolus[:, None, None, :] if is_casl else bolus,
        **basis_kwargs,
    )
    tissue *= weights
    _, _, loss = _solve_box_lsq(
        yt=xp.einsum('vd,vagd->vag', y, tissue)[:, :, None, :],
        ya=xp.einsum('vd,vbd->vb', y, all_arterial)[:, None, :, None],
        tt=xp.einsum('vagd,vagd->vag', tissue, tissue)[:, :, None, :],
        aa=xp.einsum('vbd,vbd->vb', all_arterial, all_arterial)[:, None, :, None],
        ta=xp.einsum('vagd,vbd->vabg', tissue, all_arterial),
        tol=tol[:, None, None, None],
        xp=xp,
        v_hi=all_abv_cap[:, None, :, None],
    )
    best = _argmin_att(loss, tol[:, None, None, None], xp)
    att = xp.broadcast_to(att[:, :, None, :], loss.shape)
    lower = xp.take_along_axis(att, xp.maximum(best - 1, 0), axis=3).reshape(n_voxels, -1)
    upper = xp.take_along_axis(att, xp.minimum(best + 1, grid.size - 1), axis=3)
    upper = upper.reshape(n_voxels, -1)

    # Keep the best combinations, with shape (voxels, combinations).
    loss = xp.take_along_axis(loss, best, axis=3).reshape(n_voxels, -1)
    keep = xp.argsort(_tie_rank(loss, tol[:, None], 1, xp), axis=1, kind='stable')
    keep = keep[:, :MULTIPLD_REFINED_INTERVALS]
    lower = xp.take_along_axis(lower, keep, axis=1)
    upper = xp.take_along_axis(upper, keep, axis=1)
    abat_idx = keep % n_intervals
    arterial = xp.take_along_axis(all_arterial, abat_idx[..., None], axis=1)
    abv_cap = xp.take_along_axis(all_abv_cap, abat_idx, axis=1)[..., None]
    ya = xp.einsum('vd,vcd->vc', y, arterial)[..., None]
    aa = xp.einsum('vcd,vcd->vc', arterial, arterial)[..., None]

    # Later passes, with shape (voxels, combinations, ATT grid)
    for _ in range(MULTIPLD_ATT_PASSES - 1):
        att = lower[..., None] + (upper - lower)[..., None] * grid
        tissue = _gkm_tissue_basis(
            att[..., None],
            plds[:, None, None, :],
            bolus[:, None, None, :] if is_casl else bolus,
            **basis_kwargs,
        )
        tissue *= weights
        cbf, abv, loss = _solve_box_lsq(
            yt=xp.einsum('vd,vcgd->vcg', y, tissue),
            ya=ya,
            tt=xp.einsum('vcgd,vcgd->vcg', tissue, tissue),
            aa=aa,
            ta=xp.einsum('vcgd,vcd->vcg', tissue, arterial),
            tol=tol[:, None, None],
            xp=xp,
            v_hi=abv_cap,
        )
        best = _argmin_att(loss, tol[:, None, None], xp)
        lower = xp.take_along_axis(att, xp.maximum(best - 1, 0), axis=2)[..., 0]
        upper = xp.take_along_axis(att, xp.minimum(best + 1, grid.size - 1), axis=2)[..., 0]

    att, cbf, abv, loss = (
        xp.take_along_axis(arr, best, axis=2)[..., 0] for arr in (att, cbf, abv, loss)
    )
    # Among tied combinations, prefer the one that comes first in the original order.
    tied = _tie_rank(loss, tol[:, None], 1, xp) == 0
    best = xp.argmin(xp.where(tied, keep, keep.max() + 1), axis=1)
    cbf, att, abv, loss, abat_idx = (arr[rows, best] for arr in (cbf, att, abv, loss, abat_idx))

    def solve_at(att, arterial, abv_cap):
        """Solve for CBF and aBV with ATTs of shape (voxels, n) and broadcastable arterial terms."""
        tissue = _gkm_tissue_basis(
            att[..., None],
            plds[:, None, :],
            bolus[:, None, :] if is_casl else bolus,
            **basis_kwargs,
        )
        tissue = tissue * weights
        return _solve_box_lsq(
            yt=(y[:, None, :] * tissue).sum(axis=-1),
            ya=(y[:, None, :] * arterial).sum(axis=-1),
            tt=(tissue * tissue).sum(axis=-1),
            aa=(arterial * arterial).sum(axis=-1),
            ta=(tissue * arterial).sum(axis=-1),
            tol=tol[:, None],
            xp=xp,
            v_hi=abv_cap,
        )

    # The range of ATTs that fit as well as the estimate, with aBAT in the same interval.
    # The first-pass grid brackets each end of the range, and bisection then locates it.
    tied_loss = (loss + tol)[:, None]
    best_arterial = all_arterial[rows, abat_idx][:, None, :]
    best_abv_cap = all_abv_cap[rows, abat_idx][:, None]
    grid_att = edges[:, :-1, None] + (edges[:, 1:] - edges[:, :-1])[:, :, None] * grid
    grid_att = grid_att.reshape(n_voxels, -1)
    n_points = grid_att.shape[1]
    points = xp.arange(n_points)
    tied = solve_at(grid_att, best_arterial, best_abv_cap)[2] <= tied_loss
    # The grid points on either side of the estimate
    below = xp.maximum((grid_att <= att[:, None]).sum(axis=1) - 1, 0)
    above = xp.minimum((grid_att < att[:, None]).sum(axis=1), n_points - 1)
    # The nearest grid points beyond the estimate that fit worse, if there are any
    worse_below = xp.where(~tied & (points <= below[:, None]), points, -1).max(axis=1)
    worse_above = xp.where(~tied & (points >= above[:, None]), points, n_points).min(axis=1)
    inside = xp.stack(
        [
            xp.where(tied[rows, below], grid_att[rows, worse_below + 1], att),
            xp.where(tied[rows, above], grid_att[rows, worse_above - 1], att),
        ],
        axis=1,
    )
    outside = xp.stack(
        [
            grid_att[rows, xp.maximum(worse_below, 0)],
            grid_att[rows, xp.minimum(worse_above, n_points - 1)],
        ],
        axis=1,
    )
    for _ in range(MULTIPLD_BOUND_STEPS):
        middle = (inside + outside) / 2
        is_tied = solve_at(middle, best_arterial, best_abv_cap)[2] <= tied_loss
        inside = xp.where(is_tied, middle, inside)
        outside = xp.where(is_tied, outside, middle)

    att_lower = xp.minimum(inside[:, 0], att)
    att_upper = xp.maximum(inside[:, 1], att)

    # The run of aBAT intervals that fit as well as the estimate, with ATT fixed at the estimate
    _, all_abv, all_loss = solve_at(att[:, None], all_arterial, all_abv_cap)
    intervals = xp.arange(n_intervals)
    is_best = intervals == abat_idx[:, None]
    tied = (all_loss <= tied_loss) | is_best
    all_abv = xp.where(is_best, abv[:, None], all_abv)
    all_upper = edges[:, 1:]
    if is_casl:
        # The longest aBAT in each interval at which the fitted arterial signal needs no more
        # than the largest allowed aBV
        longest = abat_mid + t1blood * xp.log(abv_max / xp.maximum(all_abv, 1e-300))
        all_upper = xp.maximum(xp.minimum(all_upper, longest), edges[:, :-1])

    # The range extends through neighboring intervals that fit as well over their whole width.
    # An interval whose longest aBATs are ruled out ends the range part of the way through it.
    partial = all_upper < edges[:, 1:]
    below = intervals < abat_idx[:, None]
    first = xp.where((~tied | partial) & below, intervals, -1).max(axis=1) + 1
    last = xp.minimum(
        xp.where(~tied & ~below & ~is_best, intervals, n_intervals).min(axis=1) - 1,
        xp.where(partial & ~below, intervals, n_intervals).min(axis=1),
    )

    abat = (edges[rows, abat_idx] + all_upper[rows, abat_idx]) / 2
    if is_casl:
        # An interval does not include its lower edge, where the arterial term changes form,
        # so an aBV at its bound leaves aBAT just above that edge. The margin of 1 ns is far
        # larger than the rounding error in the model's comparisons of aBAT with the delays.
        abat = xp.maximum(abat, xp.minimum(edges[rows, abat_idx] + 1e-9, abat_mid[rows, abat_idx]))
        # Convert the midpoint-equivalent aBV to the aBV at the reported aBAT
        abv = xp.minimum(abv * xp.exp((abat - abat_mid[rows, abat_idx]) / t1blood), abv_max)

    abat_upper = xp.maximum(all_upper[rows, last], abat)
    return cbf, att, abat, abv, att_lower, att_upper, edges[rows, first], abat_upper


def _to_numpy(arr):
    """Copy an array from any supported array library (e.g., CuPy) to NumPy."""
    return arr.get() if hasattr(arr, 'get') else np.asarray(arr)


def get_array_module(use_gpu=False):
    """Select the array library for the multi-delay model fit.

    Parameters
    ----------
    use_gpu : :obj:`bool`, optional
        Whether to run the fit on a GPU.

    Returns
    -------
    xp : module
        :mod:`cupy` if ``use_gpu`` is True, otherwise :mod:`numpy`.

    Raises
    ------
    RuntimeError
        If ``use_gpu`` is True, but CuPy cannot be imported or cannot find a CUDA device.
    """
    if not use_gpu:
        return np

    # CuPy caches its compiled kernels on disk, by default in the home directory,
    # which is read-only in some containers. Keep them in the working directory instead,
    # where later runs can reuse them.
    if config.execution.work_dir is not None:
        cache_dir = Path(config.execution.work_dir) / 'cupy_cache'
        os.environ.setdefault('CUPY_CACHE_DIR', str(cache_dir))

    try:
        import cupy
    except ImportError as exc:
        raise RuntimeError(
            'Running on a GPU requires CuPy, which could not be imported. '
            'It is included in the ASLPrep Docker image. Elsewhere, install the CuPy build '
            'that matches your CUDA version (e.g., cupy-cuda12x).'
        ) from exc

    try:
        n_devices = cupy.cuda.runtime.getDeviceCount()
    except Exception as exc:  # CuPy raises backend-specific errors here
        raise RuntimeError(f'CuPy could not query CUDA devices: {exc}') from exc

    if n_devices < 1:
        raise RuntimeError('CuPy did not find a CUDA device.')

    return cupy


def fit_deltam_multipld(
    deltam_arr,
    scaled_m0data,
    plds,
    labeleff,
    t1blood,
    partition_coefficient,
    is_casl,
    tau=None,
    ti1=None,
    xp=None,
    chunk_size=None,
    n_jobs=1,
):
    """Estimate CBF, ATT, aBAT, and aBV from multi-PLD PCASL data.

    This function uses the general kinetic model specified by
    :footcite:t:`woods2023recommendations`.

    Parameters
    ----------
    deltam_arr : :obj:`numpy.ndarray` of shape (n_voxels, n_volumes)
        Subtracted (control - label) signal.
    scaled_m0data : :obj:`numpy.ndarray` of shape (n_voxels,)
        M0 data, after any scaling requested by the user.
    plds : :obj:`numpy.ndarray` of shape (n_voxels, n_volumes)
        Delay times: post-labeling delays for (P)CASL data, inversion times for PASL data.
    labeleff : :obj:`float`
        Labeling efficiency, also referred to as alpha.
        In the case of background suppression, this combines alpha and alpha_bs,
        so we don't need a separate term for alpha_bs.
    t1blood : :obj:`float`
        Longitudinal relaxation time of arterial blood in seconds.
        Used for both deltam_tiss and deltam_art calculation.
    partition_coefficient : :obj:`float`
        Blood-brain partition coefficient, in mL/g (lambda in Alsop 2015). Generally 0.9.
    is_casl : :obj:`bool`
        True if is (P)CASL sequence.
    tau : :obj:`float` or :obj:`numpy.ndarray` or None, optional
        Label duration. May be a single value or may vary across volumes/PLDs.
        Only defined for (P)CASL data.
    ti1 : :obj:`float` or None, optional
        TI1. Only defined for PASL data.
    xp : module or None, optional
        NumPy-compatible array library used for the fit, such as :mod:`cupy` to run on a GPU.
        Defaults to NumPy.
    chunk_size : :obj:`int` or None, optional
        Number of voxels fit at once. Memory use scales linearly with this value.
        By default, NumPy uses :data:`MULTIPLD_CHUNK_SIZE`, which keeps the working arrays in
        the CPU cache, and other libraries use as many voxels as fit in
        :data:`MULTIPLD_GPU_CHUNK_ELEMENTS` elements per working array.
    n_jobs : :obj:`int`, optional
        Number of threads across which chunks are distributed.
        The array operations that dominate the fit release the GIL, so threads run in parallel.
        Memory use scales linearly with this value. The results do not depend on it.

    Returns
    -------
    cbf : :obj:`numpy.ndarray` of shape (n_voxels,)
        Cerebral blood flow.
    att : :obj:`numpy.ndarray` of shape (n_voxels,)
        Arterial transit time.
    abat : :obj:`numpy.ndarray` of shape (n_voxels,)
        Arterial bolus arrival time.
    abv : :obj:`numpy.ndarray` of shape (n_voxels,)
        Arterial blood volume.
    failures : :obj:`numpy.ndarray` of shape (n_voxels,)
        Boolean array indicating whether the model could not be fit for each voxel,
        because its M0 is not positive or its data are not finite.
        False indicates success, True indicates failure.
    bounds : :obj:`dict` of :obj:`numpy.ndarray` of shape (n_voxels,)
        The ranges of ATT and aBAT that fit the data as well as the estimates do,
        with keys ``att_lower``, ``att_upper``, ``abat_lower``, and ``abat_upper``.
        The ATT range is found with aBAT kept in its interval, and the aBAT range with ATT fixed
        at its estimate. CBF and aBV are free in both.

    Notes
    -----
    The model is linear in CBF and aBV once ATT and aBAT are fixed,
    so the least squares fit is found by searching over ATT and aBAT,
    solving for the bounded CBF and aBV exactly at each candidate.

    The tissue and arterial terms change form only where ATT or aBAT crosses a delay-dependent
    edge (see :func:`_delay_edges`), so each combination of an ATT interval and an aBAT interval
    between consecutive edges is searched separately.
    Within its interval, ATT is found with a grid search that is repeatedly narrowed around the
    best value (see :func:`_fit_gkm_chunk`), to a precision of well under 1 ms.
    aBAT is only identifiable up to its interval, so it is reported as the interval's midpoint.
    ATT is likewise not identifiable where a range of values fits equally well,
    most commonly when it is shorter than every delay, so that all of the labeled blood has
    arrived by the first one. It is then reported as the midpoint of that range.
    The ranges themselves are returned in ``bounds``. A range can extend past the interval in
    which the estimate was found, so the estimate is not always the midpoint of its bounds,
    and CBF and aBV can differ across it. Where the data determine ATT,
    its bounds differ from the estimate only by the precision of the search.
    Where aBV is zero, aBAT does not affect the model at all,
    and its bounds span every interval that fits equally well.
    For (P)CASL data, the fitted aBV depends on the reported aBAT,
    because only aBV * exp(-aBAT / T1b) is identifiable.
    If that product is too large for aBV to stay within its bound over the whole aBAT interval,
    aBAT is restricted to the lower part of the interval where it does,
    down to the interval's lower edge when aBV is at its bound.
    When solutions in different intervals fit equally well (to within rounding error),
    the one with the shortest ATT, then the shortest aBAT, then the smallest aBV is reported,
    so the results do not depend on the array library or the chunk size.

    Volumes that share the same delay (and labeling duration) in every voxel are averaged before
    fitting, with each average weighted by its number of volumes,
    which yields the same fit as using every volume.

    This model does not include a separate term for T1,tissue, but we could expand it to do so
    in the future.

    References
    ----------
    .. footbibliography::
    """
    n_voxels, n_volumes = deltam_arr.shape
    if deltam_arr.shape != plds.shape:
        raise ValueError(
            f'Number of PostLabelingDelays ({plds.shape}) does not match '
            f'number of delta-M volumes ({deltam_arr.shape}).'
        )

    if scaled_m0data.size != n_voxels:
        raise ValueError(
            f'Number of voxels in M0 data ({scaled_m0data.size}) does not match '
            f'number of delta-M voxels ({n_voxels}).'
        )

    if is_casl:
        if ti1 is not None:
            raise ValueError('ti1 should be None for PCASL data.')

        if tau is None:
            raise ValueError('tau should not be None for PCASL data.')

        if np.ndim(tau) == 0:
            tau = np.full((n_volumes,), float(tau))
        else:
            tau = np.array(tau, dtype=float)

            if tau.size != n_volumes:
                raise ValueError(
                    f'Number of values in tau ({tau.shape}) does not match number of '
                    f'delta-M volumes ({n_volumes}).'
                )
    else:
        if ti1 is None:
            raise ValueError('ti1 should not be None for PASL data.')

        if tau is not None:
            raise ValueError('tau should be None for PASL data.')

    xp = np if xp is None else xp

    # Average volumes that share the same delay (and labeling duration) in every voxel.
    delay_keys = np.vstack([plds, tau[None, :]]) if is_casl else plds
    _, first_volumes, delay_idx, counts = np.unique(
        delay_keys,
        axis=1,
        return_index=True,
        return_inverse=True,
        return_counts=True,
    )
    delay_idx = delay_idx.ravel()
    mean_deltam = np.stack(
        [deltam_arr[:, delay_idx == i_delay].mean(axis=1) for i_delay in range(counts.size)],
        axis=1,
    )
    unique_plds = plds[:, first_volumes]
    weights = np.sqrt(counts)[None, :]

    m0_a = scaled_m0data / partition_coefficient
    valid = (
        np.isfinite(m0_a)
        & (m0_a > 0)
        & np.all(np.isfinite(mean_deltam), axis=1)
        & np.all(np.isfinite(unique_plds), axis=1)
    )
    y = np.zeros_like(mean_deltam)
    y[valid] = mean_deltam[valid] / m0_a[valid, None] * weights

    bolus = xp.asarray(tau[first_volumes][None, :]) if is_casl else float(ti1)
    weights = xp.asarray(weights)
    if chunk_size is None and xp is np:
        chunk_size = MULTIPLD_CHUNK_SIZE
    elif chunk_size is None:
        # The largest working arrays hold one element per voxel, per combination of an ATT
        # interval and an aBAT interval, per ATT grid point.
        n_intervals = 2 * counts.size + 1
        elements_per_voxel = n_intervals**2 * MULTIPLD_ATT_GRID_SIZE
        chunk_size = max(1, MULTIPLD_GPU_CHUNK_ELEMENTS // elements_per_voxel)

    estimates = np.full((8, n_voxels), np.nan)
    valid_idx = np.flatnonzero(valid)
    chunks = [
        valid_idx[start : start + chunk_size] for start in range(0, valid_idx.size, chunk_size)
    ]

    def fit_chunk(chunk):
        chunk_estimates = _fit_gkm_chunk(
            y=xp.asarray(y[chunk]),
            plds=xp.asarray(unique_plds[chunk]),
            bolus=bolus,
            weights=weights,
            is_casl=is_casl,
            labeleff=labeleff,
            t1blood=t1blood,
            xp=xp,
        )
        # Chunks hold different voxels, so threads never write to the same elements.
        estimates[:, chunk] = [_to_numpy(estimate) for estimate in chunk_estimates]

    if n_jobs > 1 and len(chunks) > 1:
        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            # Consume the iterator so that exceptions raised in the threads propagate.
            list(executor.map(fit_chunk, chunks))
    else:
        for chunk in chunks:
            fit_chunk(chunk)

    failures = ~valid
    fail_count = int(failures.sum())
    if fail_count:
        fail_percent = 100 * fail_count / n_voxels
        config.loggers.workflow.warning(
            f'General kinetic model could not be fit on '
            f'{fail_count}/{n_voxels} ({fail_percent:.2f}%) voxel(s) '
            'with non-positive M0 or non-finite data.'
        )

    cbf, att, abat, abv, *bounds = estimates
    bounds = dict(zip(('att_lower', 'att_upper', 'abat_lower', 'abat_upper'), bounds, strict=True))
    return cbf, att, abat, abv, failures, bounds


#: M0 below this fraction of its 98th percentile over the field of view is treated as absent.
M0_THRESHOLD = 0.1

#: Blood-brain partition coefficient (lambda), in mL/g, from Alsop et al. (2015).
PARTITION_COEFFICIENT = 0.9


def zero_low_m0(m0data, threshold=M0_THRESHOLD):
    """Set M0 values with almost no signal to zero.

    CBF is proportional to delta-M divided by M0, so voxels with almost no M0 signal
    (air, mostly) would otherwise get arbitrarily large CBF values,
    which resampling then spreads into neighboring voxels.
    A voxel with an M0 of zero gets no CBF estimate.

    Parameters
    ----------
    m0data : :obj:`numpy.ndarray`
        M0 values.
    threshold : :obj:`float`, optional
        Fraction of the 98th percentile of the finite M0 values below which M0 is set to zero.
        The 98th percentile stands in for the brain's M0, without relying on a brain mask.

    Returns
    -------
    m0data : :obj:`numpy.ndarray`
        A copy of the M0 values, with the low ones set to zero.
        Non-finite values are left as they are.
    """
    m0data = np.array(m0data, dtype=float)
    finite = np.isfinite(m0data)
    if not finite.any():
        return m0data

    cutoff = threshold * np.percentile(m0data[finite], 98)
    m0data[finite & (m0data < cutoff)] = 0
    return m0data


# Oxygen saturation of arterial blood, as used by Hales et al. (2016) for arterial T1.
ARTERIAL_OXYGEN_SATURATION = 0.97
# Field strengths, in tesla, over which the Hales et al. (2016) model was fitted and validated.
HALES_FIELD_STRENGTH_RANGE = (1.5, 7.0)
# Haematocrit that the default T1blood values are taken to assume:
# the "typical human Hct" of Lu et al. (2004).
REFERENCE_HAEMATOCRIT = 0.42


def t1blood_from_haematocrit(haematocrit, field_strength, oxygen_saturation=None):
    """Estimate the T1 of blood from haematocrit, with the model of Hales et al. (2016).

    The model treats blood as erythrocytes and plasma in fast exchange,
    with relaxation rates that depend on the magnetic field strength,
    and on oxygen saturation for the erythrocytes.
    It was fitted to in vitro measurements between 1.5 and 7 T,
    and includes the authors' +108 ms offset for in vivo blood.

    Parameters
    ----------
    haematocrit : :obj:`float`
        Haematocrit, as a fraction of whole-blood volume (0-1).
    field_strength : :obj:`float`
        Magnetic field strength, in tesla.
    oxygen_saturation : :obj:`float`, optional
        Oxygen saturation of the blood (0-1).
        Defaults to 0.97, for arterial blood.

    Returns
    -------
    t1blood : :obj:`float`
        T1 of blood, in seconds.

    Notes
    -----
    From :footcite:t:`hales2016general`, equations 7 and 8:

    .. math::

        f_e = \\frac{0.70 \\cdot Hct}{0.70 \\cdot Hct + 0.95 \\cdot (1 - Hct)}

        R_{1,bl} = f_e \\cdot (1.099 - 0.057 B_0 + 0.033 [Hb] (1 - Y))
                   + (1 - f_e) \\cdot (0.496 - 0.023 B_0)

        T_{1,bl} = 1 / R_{1,bl} + 0.108

    where [Hb] is the mean corpuscular haemoglobin concentration (5.15 mmol/L).

    References
    ----------
    .. footbibliography::
    """
    if oxygen_saturation is None:
        oxygen_saturation = ARTERIAL_OXYGEN_SATURATION

    haemoglobin = 5.15  # mmol Hb tetramer/L plasma
    in_vivo_offset = 0.108  # seconds
    # Fraction of the water in whole blood that resides in erythrocytes
    fe = (0.70 * haematocrit) / (0.70 * haematocrit + 0.95 * (1 - haematocrit))
    r1_erythrocytes = (
        1.099 - 0.057 * field_strength + 0.033 * haemoglobin * (1 - oxygen_saturation)
    )
    r1_plasma = 0.496 - 0.023 * field_strength
    r1_blood = fe * r1_erythrocytes + (1 - fe) * r1_plasma
    return 1 / r1_blood + in_vivo_offset


def haematocrit_applies(metadata):
    """Determine whether T1blood should be calculated from the haematocrit in the metadata.

    The haematocrit is used when it is present and the field strength is within the range of
    the model of :footcite:t:`hales2016general` (1.5-7 T).

    References
    ----------
    .. footbibliography::
    """
    haematocrit = metadata.get('Haematocrit')
    if haematocrit is None:
        return False

    field_strength = metadata['MagneticFieldStrength']
    low, high = HALES_FIELD_STRENGTH_RANGE
    if not low <= field_strength <= high:
        config.loggers.interface.warning(
            f'A haematocrit of {haematocrit:g} is available, but the model that uses it is only '
            f'valid between {low} and {high} T, not at {field_strength} T. '
            'T1blood will be set from the field strength alone.'
        )
        return False

    return True


def estimate_t1(metadata):
    """Estimate the relaxation rates of blood and gray matter based on magnetic field strength.

    t1blood is set based on the scanner's field strength,
    according to :footcite:t:`zhang2013vivo,alsop_recommended_2015`.
    If recommended values from these publications cannot be used
    (i.e., if the field strength isn't 1.5T, 3T, 7T),
    then the formula from :footcite:t:`zhang2013vivo` will be applied.

    If the metadata include a ``Haematocrit`` (as a fraction),
    and the field strength is between 1.5 and 7 T,
    that value is then scaled for the participant's haematocrit,
    with the model of :footcite:t:`hales2016general` (see :func:`t1blood_from_haematocrit`):

    .. math::

        T_{1,blood} = T_{1,blood}^{default} \\cdot
        \\frac{T_{1,bl}^{Hales}(Hct)}{T_{1,bl}^{Hales}(0.42)}

    The defaults assume a typical haematocrit, taken as 0.42 (the "typical human Hct" of
    :footcite:t:`lu2004determining`, whose measurements underlie the 3 T default),
    so a haematocrit of 0.42 leaves them unchanged.
    The published model on its own gives 7-21% longer values than the defaults at a typical
    haematocrit, which would make CBF from participants with and without a measured haematocrit
    differ systematically.

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
    field_strength = metadata['MagneticFieldStrength']
    T1BLOOD_DICT = {
        1.5: 1.35,
        3: 1.65,
        7: 2.087,
    }
    t1blood = T1BLOOD_DICT.get(field_strength)
    if not t1blood:
        config.loggers.interface.warning(
            f'T1blood cannot be inferred for {metadata["MagneticFieldStrength"]}T data. '
            'Defaulting to formula from Zhang et al. (2013).'
        )
        t1blood = (110 * metadata['MagneticFieldStrength'] + 1316) / 1000

    if haematocrit_applies(metadata):
        t1blood *= t1blood_from_haematocrit(
            metadata['Haematocrit'], field_strength
        ) / t1blood_from_haematocrit(REFERENCE_HAEMATOCRIT, field_strength)

    # TODO: Supplement with formula for other field strengths
    T1TISSUE_DICT = {
        1.5: 1.197,
        3: 1.607,
        7: 1.939,
    }
    t1tissue = T1TISSUE_DICT.get(metadata['MagneticFieldStrength'])
    if not t1tissue:
        raise ValueError(
            f'T1tissue cannot be inferred for {metadata["MagneticFieldStrength"]}T data.'
        )

    return t1blood, t1tissue
