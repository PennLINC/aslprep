# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Miscellaneous utilities."""
from __future__ import annotations

import os
from typing import Any

import nibabel as nb
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import median_abs_deviation


def check_deps(workflow):
    """Make sure dependencies are present in this system."""
    from nipype.utils.filemanip import which

    return sorted(
        (node.interface.__class__.__name__, node.interface._cmd)
        for node in workflow._get_all_nodes()
        if (hasattr(node.interface, "_cmd") and which(node.interface._cmd.split()[0]) is None)
    )


def _get_series_len(asl_fname):
    """Determine the number of volumes in an image, after removing outlier volumes."""
    from aslprep.niworkflows.interfaces.registration import _get_vols_to_discard

    img = nb.load(asl_fname)
    if len(img.shape) < 4:
        return 1

    skip_vols = _get_vols_to_discard(img)

    return img.shape[3] - skip_vols


def _create_mem_gb(asl_fname):
    """Estimate the memory needed for different operations, based on the size of the data."""
    asl_size_gb = os.path.getsize(asl_fname) / (1024**3)
    asl_tlen = nb.load(asl_fname).shape[-1]
    mem_gb = {
        "filesize": asl_size_gb,
        "resampled": asl_size_gb * 4,
        "largemem": asl_size_gb * (max(asl_tlen / 100, 1.0) + 4),
    }

    return asl_tlen, mem_gb


def _get_wf_name(asl_fname):
    """Derive the workflow name for a supplied ASL file.

    >>> _get_wf_name('/completely/made/up/path/sub-01_task-nback_asl.nii.gz')
    'func_preproc_task_nback_wf'
    >>> _get_wf_name('/completely/made/up/path/sub-01_task-nback_run-01_echo-1_asl.nii.gz')
    'func_preproc_task_nback_run_01_echo_1_wf'

    """
    from nipype.utils.filemanip import split_filename

    fname = split_filename(asl_fname)[1]
    fname_nosub = "_".join(fname.split("_")[1:])
    name = "asl_preproc_" + fname_nosub.replace(".", "_").replace(" ", "").replace(
        "-", "_"
    ).replace("_asl", "_wf")

    return name


def get_n_volumes(fname):
    """Get the number of volumes in a niimg file."""
    img = nb.load(fname)
    if img.ndim == 3:
        n_volumes = 0
    elif img.ndim == 4:
        n_volumes = img.shape[3]
    else:
        raise ValueError(f"Image has {img.ndim} dimensions: {fname}")

    return n_volumes


def _split_spec(in_target):
    """Split space-resolution specification into space, template, and remaining info."""
    space, spec = in_target
    template = space.split(":")[0]
    return space, template, spec


def _select_template(template):
    """Select template file based on space/template specification."""
    from aslprep.niworkflows.utils.misc import get_template_specs

    template, specs = template
    template = template.split(":")[0]  # Drop any cohort modifier if present
    specs = specs.copy()
    specs["suffix"] = specs.get("suffix", "T1w")

    # Sanitize resolution
    res = specs.pop("res", None) or specs.pop("resolution", None) or "native"
    if res != "native":
        specs["resolution"] = res
        return get_template_specs(template, template_spec=specs)[0]

    # Map nonstandard resolutions to existing resolutions
    specs["resolution"] = 2
    try:
        out = get_template_specs(template, template_spec=specs)
    except RuntimeError:
        specs["resolution"] = 1
        out = get_template_specs(template, template_spec=specs)

    return out[0]


def _aslist(in_value):
    """Convert input to a list.

    TODO: Replace with listify from something like NiMARE.
    """
    if isinstance(in_value, list):
        return in_value
    return [in_value]


def _is_native(in_value):
    """Determine if input dictionary references native-space data."""
    return in_value.get("resolution") == "native" or in_value.get("res") == "native"


def _select_last_in_list(lst):
    """Select the last element in a list."""
    return lst[-1]


def _conditional_downsampling(in_file, in_mask, zoom_th=4.0):
    """Downsample the input dataset for sloppy mode."""
    from pathlib import Path

    import nibabel as nb
    import nitransforms as nt
    import numpy as np
    from scipy.ndimage.filters import gaussian_filter

    img = nb.load(in_file)

    zooms = np.array(img.header.get_zooms()[:3])
    if not np.any(zooms < zoom_th):
        return in_file, in_mask

    out_file = Path("desc-resampled_input.nii.gz").absolute()
    out_mask = Path("desc-resampled_mask.nii.gz").absolute()

    shape = np.array(img.shape[:3])
    scaling = zoom_th / zooms
    newrot = np.diag(scaling).dot(img.affine[:3, :3])
    newshape = np.ceil(shape / scaling).astype(int)
    old_center = img.affine.dot(np.hstack((0.5 * (shape - 1), 1.0)))[:3]
    offset = old_center - newrot.dot((newshape - 1) * 0.5)
    newaffine = nb.affines.from_matvec(newrot, offset)

    newref = nb.Nifti1Image(np.zeros(newshape, dtype=np.uint8), newaffine)
    nt.Affine(reference=newref).apply(img).to_filename(out_file)

    mask = nb.load(in_mask)
    mask.set_data_dtype(float)
    mdata = gaussian_filter(mask.get_fdata(dtype=float), scaling)
    floatmask = nb.Nifti1Image(mdata, mask.affine, mask.header)
    newmask = nt.Affine(reference=newref).apply(floatmask)
    hdr = newmask.header.copy()
    hdr.set_data_dtype(np.uint8)
    newmaskdata = (newmask.get_fdata(dtype=float) > 0.5).astype(np.uint8)
    nb.Nifti1Image(newmaskdata, newmask.affine, hdr).to_filename(out_mask)

    return str(out_file), str(out_mask)


def select_target(subject_id, space):
    """Get the target subject ID, given a source subject ID and a target space."""
    return subject_id if space == "fsnative" else space


def _select_first_in_list(lst):
    """Select the first element in a list."""
    return lst[0]


def _itk2lta(in_file, src_file, dst_file):
    """Convert ITK file to LTA file."""
    from pathlib import Path

    import nitransforms as nt

    out_file = Path("out.lta").absolute()
    lta_object = nt.linear.load(
        in_file,
        fmt="fs" if in_file.endswith(".lta") else "itk",
        reference=src_file,
    )
    lta_object.to_filename(out_file, moving=dst_file, fmt="fs")
    return str(out_file)


def _prefix(subid):
    """Add sub- prefix to subject ID, if necessary."""
    return subid if subid.startswith("sub-") else f"sub-{subid}"


def pcasl_or_pasl(metadata):
    """Determine if metadata indicates a PCASL or ASL scan."""
    aslt = metadata["ArterialSpinLabelingType"]

    if aslt in ["CASL", "PCASL"]:
        is_casl = True
    elif aslt == "PASL":
        is_casl = False
    else:
        raise ValueError(
            "Labeling type cannot be classified as (P)CASL or PASL based on "
            f"ArterialSpinLabelingType: '{aslt}'."
        )

    return is_casl


def readjson(jsonfile):
    """Read a JSON file into memory."""
    import json

    with open(jsonfile, "r") as f:
        data = json.load(f)
    return data


def get_tis(metadata: "dict[str, Any]") -> list:
    """Determine inversion times from metadata.

    Parameters
    ----------
    metadata : :obj:`dict`
        Dictionary of metadata associated with ASL file.

    Returns
    -------
    :obj:`list`
        List of PostLabelingDelay values.
    """
    if "CASL" in metadata["ArterialSpinLabelingType"]:
        return np.add(metadata["PostLabelingDelay"], metadata["LabelingDuration"]).tolist()
    else:
        return np.array(metadata["PostLabelingDelay"]).tolist()


def _weightfun(x, wfun="huber"):
    """Get weight fun and tuner."""
    if wfun == "andrews":
        tuner = 1.339
        weight = (np.abs(x) < np.pi) * np.sin(x)
    elif wfun == "bisquare":
        tuner = 4.685
        weight = (np.abs(x) < 1) * np.power((1 - np.power(x, 2)), 2)
    elif wfun == "cauchy":
        tuner = 2.385
        weight = 1 / (1 + np.power(x, 2))
    elif wfun == "logistic":
        tuner = 1.205
        weight = np.tanh(x) / x
    elif wfun == "ols":
        tuner = 1
        weight = np.repeat(1, len(x))
    elif wfun == "talwar":
        tuner = 2.795
        weight = 1 * (np.abs(x) < 1)
    elif wfun == "welsch":
        tuner = 2.985
        weight = np.exp(-(np.power(x, 2)))
    else:
        tuner = 1.345
        weight = 1 / np.abs(x)
    return weight, tuner


def _tune(wfun="huber"):
    """Get weight fun and tuner.

    But wait, you might say, the docstring makes no sense! Correct.
    """
    if wfun == "andrews":
        tuner = 1.339
    elif wfun == "bisquare":
        tuner = 4.685
    elif wfun == "cauchy":
        tuner = 2.385
    elif wfun == "logistic":
        tuner = 1.205
    elif wfun == "ols":
        tuner = 1
    elif wfun == "talwar":
        tuner = 2.795
    elif wfun == "welsch":
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
       nd array of 3D or 4D computed cbf
    gm,wm,csf
       numpy array of grey matter, whitematter, and csf
    mask
       numpy array of mask
    """
    gm_bin = gm >= thresh
    wm_bin = wm >= thresh
    csf_bin = csf >= thresh

    # get the total number of voxle within csf,gm and wm
    n_gm_voxels = np.sum(gm_bin) - 1
    n_wm_voxels = np.sum(wm_bin) - 1
    n_csf_voxels = np.sum(csf_bin) - 1
    mask1 = gm_bin + wm_bin + csf_bin

    # mean  of times series cbf within greymatter
    gm_cbf_ts = np.squeeze(np.mean(cbfts[gm_bin, :], axis=0))

    # robust mean and median
    median_gm_cbf = np.median(gm_cbf_ts)
    mad_gm_cbf = median_abs_deviation(gm_cbf_ts) / 0.675
    indx = 1 * (np.abs(gm_cbf_ts - median_gm_cbf) > (2.5 * mad_gm_cbf))
    R = np.mean(cbfts[:, :, :, indx == 0], axis=3)
    V = (
        n_gm_voxels * np.var(R[gm == 1])
        + n_wm_voxels * np.var(R[wm == 1])
        + n_csf_voxels * np.var(R[csf == 1])
    )
    V1 = V + 1
    while V < V1:
        V1 = V
        CC = np.zeros(cbfts.shape[3]) * (-2)
        for s in range(cbfts.shape[3]):
            if indx[s] != 0:
                break
            else:
                tmp1 = cbfts[:, :, :, s]
                CC[s] = np.corrcoef(R[mask1 > 0], tmp1[mask1 > 0])[0][1]

        inx = np.argmax(CC)
        indx[inx] = 2
        R = np.mean(cbfts[:, :, :, indx == 0], axis=3)
        V = (
            (n_gm_voxels * np.var(R[gm == 1]))
            + (n_wm_voxels * np.var(R[wm == 1]))
            + (n_csf_voxels * np.var(R[csf == 1]))
        )
    cbfts_recon = cbfts[:, :, :, indx == 0]
    cbfts_recon1 = np.zeros_like(cbfts_recon)
    for i in range(cbfts_recon.shape[3]):
        cbfts_recon1[:, :, :, i] = cbfts_recon[:, :, :, i] * mask

    cbfts_recon1 = np.nan_to_num(cbfts_recon1)
    return cbfts_recon1, indx


def _robust_fit(
    Y,
    mu,
    Globalprior,
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
        print("iteration  ", iter_num, "\n")
        iter_num = iter_num + 1
        check1 = np.subtract(np.abs(b - b0), (D * np.maximum(np.abs(b), np.abs(b0))))
        check1[check1 > 0] = 0
        if any(check1):
            print(" \n converged after ", iter_num, "iterations\n")
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


def _scrubcbf(cbf_ts, gm, wm, csf, mask, wfun="huber", thresh=0.7):
    """Apply SCRUB algorithm to CBF data.

    Parameters
    ----------
    cbf_ts
       nd array of 3D or 4D computed cbf
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
    gmidx2 = (
        1 * ([gm.flatten() > thresh] and [M.flatten() == 0] and [wm.flatten() > csf.flatten()])[0]
    )
    wmidx2 = (
        1 * ([wm.flatten() > thresh] and [M.flatten() == 0] and [gm.flatten() > csf.flatten()])[0]
    )
    if np.sum(gmidx2) == 0 or np.sum(wmidx2) == 0:
        gmidx2 = 1 * (gm.flatten() > thresh)
        wmidx2 = 1 * (wm.flatten() > thresh)
    idxx = gmidx2 + wmidx2
    idxx[idxx > 0] = 1
    X = np.zeros([len(idxx), 2])
    X[:, 0] = gm.flatten()[gm.flatten() >= (0)] * idxx
    X[:, 1] = wm.flatten()[wm.flatten() >= (0)] * idxx
    A = (mean_cbf.flatten()[idxx >= 0]) * idxx
    c = np.linalg.lstsq(X, A)[0]
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


def estimate_labeling_efficiency(metadata):
    """Estimate labeling efficiency based on the available metadata.

    Parameters
    ----------
    metadata : :obj:`dict`
        Dictionary of metadata from the ASL file.

    Returns
    -------
    labeleff : :obj:`float`
        Labeling efficiency.

    Notes
    -----
    If LabelingEfficiency is defined in the metadata, then this value will be used.
    Otherwise, efficiency will be estimated based on the ASL type and number of background
    suppression pulses (if any).
    PCASL and PASL values come from :footcite:t:`alsop_recommended_2015`.
    The CASL value comes from :footcite:t:`wang2005amplitude`.
    The adjustment based on number of background suppression pulses is not described in any papers,
    as far as we know, but is apparently widely used.

    References
    ----------
    .. footbibliography::
    """
    if "LabelingEfficiency" in metadata.keys():
        labeleff = metadata["LabelingEfficiency"]
    else:
        BASE_LABELEFF = {
            "CASL": 0.68,
            "PCASL": 0.85,
            "PASL": 0.98,
        }
        labeleff = BASE_LABELEFF[metadata["ArterialSpinLabelingType"]]

        if metadata.get("BackgroundSuppression", False):
            BS_PULSE_EFF = 0.95  # hardcoded BackgroundSuppressionPulse efficiency
            # We assume there was one pulse if suppression was applied,
            # but the number of pulses isn't defined.
            labeleff *= BS_PULSE_EFF ** metadata.get("BackgroundSuppressionNumberPulses", 1)

    return labeleff


def estimate_att_pcasl(deltam_arr, plds, lds, t1blood, t1tissue):
    """Estimate arterial transit time using the weighted average method.

    Parameters
    ----------
    deltam_arr : :obj:`numpy.ndarray`
        Delta-M array.
    plds : :obj:`numpy.ndarray`
        Post-labeling delays. w in Dai 2012.
    lds : :obj:`numpy.ndarray`
        Labeling durations. tau in Dai 2012.
    t1blood
        T1 relaxation rate for blood.
    t1tissue
        T1 relaxation rate for tissue.

    Returns
    -------
    att_arr : :obj:`numpy.ndarray`
        Arterial transit time array.

    Notes
    -----
    Written by Jianxun Qu, @jianxun_qu@163.com in MATLAB.
    Translated to Python by Taylor Salo.
    """
    n_plds = plds.size
    n_voxels, n_volumes = deltam_arr.shape
    assert n_volumes == n_plds, f"{n_plds} != {n_volumes}"

    # Beginning of auxil_asl_gen_wsum
    assert lds.size == n_plds, f"{lds.size} != {n_plds}"

    # Define the possible transit times to evaluate
    transit_times = np.arange(np.round(np.min(plds), 3), np.round(np.max(plds), 3) + 0.001, 0.001)

    sig_sum = np.zeros(transit_times.size)
    sig_pld_sum = np.zeros(transit_times.size)

    for i_pld in range(n_plds):
        pld = plds[i_pld]
        ld = lds[i_pld]

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
        # Beta (a term to compensate for static tissue signal loss caused by vessel supp. pulses),
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
    weighted_delay_predicted = sig_pld_sum / sig_sum
    # End of auxil_asl_gen_wsum

    # Calculate the observed weighted delay for each voxel
    weighted_delay_denom = np.sum(deltam_arr, axis=1)
    weighted_delay_num = np.zeros(n_voxels)
    for i_pld in range(n_plds):
        weighted_delay_num += deltam_arr[:, i_pld] * plds[i_pld]

    weighted_delay_observed = weighted_delay_num / (
        np.abs(weighted_delay_denom) + np.finfo(float).eps
    )

    # Truncate extreme transit time values to the PLD limits
    weighted_delay_min = min(weighted_delay_predicted)
    weighted_delay_max = max(weighted_delay_predicted)
    weighted_delay_observed[weighted_delay_observed < weighted_delay_min] = weighted_delay_min
    weighted_delay_observed[weighted_delay_observed > weighted_delay_max] = weighted_delay_max

    # Find the ATT for each voxel based on the predicted weighted delay value for each transit
    # time and the actual weighted delay value.
    interp_func = interp1d(weighted_delay_predicted, transit_times)
    return interp_func(weighted_delay_observed)


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
    """Estimate CBF and ATT for multi-PLD PCASL data."""
    n_voxels, n_volumes = deltam_arr.shape
    if plds.size != n_volumes:
        raise ValueError(
            f"Number of PostLabelingDelay values ({plds.size}) does not match number of delta-M "
            "volumes."
        )

    # Formula from Fan 2017 (equation 2)
    unique_plds, unique_pld_idx = np.unique(plds, return_index=True)
    if tau.size > 1:
        if tau.size != plds.size:
            raise ValueError(
                f"Number of LabelingDuration values {tau.size} != number of "
                f"PostLabelingDelay values {plds.size}"
            )

        tau = tau[unique_pld_idx]
    else:
        tau = np.full(unique_plds.size, tau)

    mean_deltam_by_pld = np.zeros((n_voxels, unique_plds.size))
    for i_pld, pld in enumerate(unique_plds):
        pld_idx = plds == pld
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
    print(f"num_factor: {num_factor}")
    denom_factor = 2 * labeleff * scaled_m0data * t1blood
    print(f"denom_factor: {denom_factor.shape}")

    # Loop over PLDs and calculate CBF for each, accounting for ATT.
    cbf_by_pld = np.zeros((n_voxels, unique_plds.size))
    for i_pld, pld in enumerate(unique_plds):
        tau_for_pld = tau[i_pld]

        pld_num_factor = (
            num_factor * mean_deltam_by_pld[:, i_pld] * np.exp(att_arr[:, None] / t1blood)
        )
        print(f"mean_deltam_by_pld[:, i_pld]: {mean_deltam_by_pld[:, i_pld].shape}")
        print(f"att_arr[:, None]: {att_arr[:, None].shape}")
        print(f"pld_num_factor: {pld_num_factor.shape}")

        pld_denom_factor = denom_factor * (
            np.exp(-np.maximum(pld - att_arr, 0))
            - np.exp(-np.maximum(tau_for_pld + pld - att_arr, 0))
        )
        print(
            "np.exp(-np.maximum(pld - att_arr, 0)): "
            f"{np.exp(-np.maximum(pld - att_arr, 0)).shape}"
        )
        print(
            "np.exp(-np.maximum(tau_for_pld + pld - att_arr, 0)): "
            f"{np.exp(-np.maximum(tau_for_pld + pld - att_arr, 0)).shape}"
        )
        print(f"pld_denom_factor: {pld_denom_factor.shape}")
        cbf_by_pld[:, i_pld] = pld_num_factor / pld_denom_factor

    # Average CBF across PLDs, but only include PLDs where PLD + tau > ATT for that voxel,
    # per Juttukonda 2021 (section 2.6).
    cbf = np.zeros(n_voxels)  # mean CBF
    for i_voxel in range(n_voxels):
        cbf_by_pld_voxel = cbf_by_pld[i_voxel, :]
        arr_voxel = att_arr[i_voxel]
        cbf[i_voxel] = np.mean(cbf_by_pld_voxel[(unique_plds + tau) > arr_voxel])

    return att_arr, cbf


def determine_multi_pld(metadata):
    """Determine if a run is multi-PLD or not."""
    plds = np.array(metadata["PostLabelingDelay"])
    return np.unique(plds).size > 1
