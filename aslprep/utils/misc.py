# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Miscellaneous utilities."""
import os
from pathlib import Path

import nibabel as nb
import numpy as np
from nipype.utils.filemanip import fname_presuffix


def check_deps(workflow):
    """Make sure dependencies are present in this system."""
    from nipype.utils.filemanip import which

    return sorted(
        (node.interface.__class__.__name__, node.interface._cmd)
        for node in workflow._get_all_nodes()
        if (hasattr(node.interface, "_cmd") and which(node.interface._cmd.split()[0]) is None)
    )


def _get_series_len(asl_fname):
    from aslprep.niworkflows.interfaces.registration import _get_vols_to_discard

    img = nb.load(asl_fname)
    if len(img.shape) < 4:
        return 1

    skip_vols = _get_vols_to_discard(img)

    return img.shape[3] - skip_vols


def _create_mem_gb(asl_fname):
    asl_size_gb = os.path.getsize(asl_fname) / (1024**3)
    asl_tlen = nb.load(asl_fname).shape[-1]
    mem_gb = {
        "filesize": asl_size_gb,
        "resampled": asl_size_gb * 4,
        "largemem": asl_size_gb * (max(asl_tlen / 100, 1.0) + 4),
    }

    return asl_tlen, mem_gb


def _get_wf_name(asl_fname):
    """Derive the workflow name for supplied asl file.

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
    """Get number of volumes in niimg file."""
    img = nb.load(fname)
    if img.ndim == 3:
        n_volumes = 0
    elif img.ndim == 4:
        n_volumes = img.shape[3]
    else:
        raise ValueError(f"Image has {img.ndim} dimensions: {fname}")

    return n_volumes


def gen_reference(in_img, fwhm=5, newpath=None):
    """Generate reference for a GE scan with few volumes."""
    newpath = Path(newpath or ".")
    ss = get_n_volumes(in_img)
    if ss == 0:
        ref_data = nb.load(in_img).get_fdata()
    else:
        nii = nb.load(in_img).get_fdata()
        ref_data = np.mean(nii, axis=3)
    new_file = nb.Nifti1Image(
        dataobj=ref_data, header=nb.load(in_img).header, affine=nb.load(in_img).affine
    )

    new_file = nb.processing.smooth_image(new_file, fwhm=fwhm)
    out_file = fname_presuffix(
        "aslref",
        suffix="_reference.nii.gz",
        newpath=str(newpath.absolute()),
    )
    new_file.to_filename(out_file)
    return out_file


def _split_spec(in_target):
    space, spec = in_target
    template = space.split(":")[0]
    return space, template, spec


def _select_template(template):
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
    if isinstance(in_value, list):
        return in_value
    return [in_value]


def _is_native(in_value):
    return in_value.get("resolution") == "native" or in_value.get("res") == "native"


def _select_last_in_list(lst):
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
    return lst[0]


def _itk2lta(in_file, src_file, dst_file):
    from pathlib import Path

    import nitransforms as nt

    out_file = Path("out.lta").absolute()
    nt.linear.load(
        in_file, fmt="fs" if in_file.endswith(".lta") else "itk", reference=src_file
    ).to_filename(out_file, moving=dst_file, fmt="fs")
    return str(out_file)


def _prefix(subid):
    return subid if subid.startswith("sub-") else f"sub-{subid}"
