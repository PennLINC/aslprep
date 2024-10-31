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
        if (hasattr(node.interface, '_cmd') and which(node.interface._cmd.split()[0]) is None)
    )


def get_n_volumes(fname):
    """Get the number of volumes in a niimg file."""
    img = nb.load(fname)
    if img.ndim == 3:
        n_volumes = 0
    elif img.ndim == 4:
        n_volumes = img.shape[3]
    else:
        raise ValueError(f'Image has {img.ndim} dimensions: {fname}')

    return n_volumes


def _create_mem_gb(asl_fname):
    """Estimate the memory needed for different operations, based on the size of the data."""
    asl_size_gb = os.path.getsize(asl_fname) / (1024**3)
    asl_tlen = nb.load(asl_fname).shape[-1]
    mem_gb = {
        'filesize': asl_size_gb,
        'resampled': asl_size_gb * 4,
        'largemem': asl_size_gb * (max(asl_tlen / 100, 1.0) + 4),
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
    fname_nosub = '_'.join(fname.split('_')[1:])
    name = 'asl_preproc_' + fname_nosub.replace('.', '_').replace(' ', '').replace(
        '-', '_'
    ).replace('_asl', '_wf')

    return name


def _select_last_in_list(lst):
    """Select the last element in a list."""
    return lst[-1]


def _prefix(subid):
    """Add sub- prefix to subject ID, if necessary."""
    return subid if subid.startswith('sub-') else f'sub-{subid}'


def estimate_asl_mem_usage(asl_fname: str) -> tuple[int, dict]:
    """Estimate ASL memory usage."""
    import nibabel as nb
    import numpy as np

    img = nb.load(asl_fname)
    nvox = int(np.prod(img.shape, dtype='u8'))
    # Assume tools will coerce to 8-byte floats to be safe
    asl_size_gb = 8 * nvox / (1024**3)
    asl_tlen = img.shape[-1]
    mem_gb = {
        'filesize': asl_size_gb,
        'resampled': asl_size_gb * 4,
        'largemem': asl_size_gb * (max(asl_tlen / 100, 1.0) + 4),
    }

    return asl_tlen, mem_gb
