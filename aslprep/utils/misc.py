# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Miscellaneous utilities."""
from __future__ import annotations

import os

import nibabel as nb


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
    """Select the first element in a list.

    If the input isn't a list, then it will return the whole input.
    """
    return lst[0] if isinstance(lst, list) else lst


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


def readjson(jsonfile):
    """Read a JSON file into memory."""
    import json

    with open(jsonfile, "r") as f:
        data = json.load(f)
    return data


def get_template_str(template, kwargs):
    """Get template from templateflow, as a string."""
    from templateflow.api import get as get_template

    return str(get_template(template, **kwargs))
