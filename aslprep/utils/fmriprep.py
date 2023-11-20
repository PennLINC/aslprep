"""Functions from fMRIPrep.

TODO: Replace with fMRIPrep imports once next is merged and released.
"""
import asyncio
import warnings
from pathlib import Path
from typing import Callable, TypeVar

import h5py
import nibabel as nb
import nitransforms as nt
import numpy as np
from nitransforms.io.itk import ITKCompositeH5

R = TypeVar("R")


async def worker(job: Callable[[], R], semaphore: asyncio.Semaphore) -> R:
    async with semaphore:
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, job)


def load_transforms(xfm_paths: list[Path], inverse: list[bool]) -> nt.base.TransformBase:
    """Load a series of transforms as a nitransforms TransformChain

    An empty list will return an identity transform
    """
    if len(inverse) == 1:
        inverse *= len(xfm_paths)
    elif len(inverse) != len(xfm_paths):
        raise ValueError("Mismatched number of transforms and inverses")

    chain = None
    for path, inv in zip(xfm_paths[::-1], inverse[::-1]):
        path = Path(path)
        if path.suffix == ".h5":
            xfm = load_ants_h5(path)
        else:
            xfm = nt.linear.load(path)
        if inv:
            xfm = ~xfm
        if chain is None:
            chain = xfm
        else:
            chain += xfm
    if chain is None:
        chain = nt.base.TransformBase()
    return chain


FIXED_PARAMS = np.array([
    193.0, 229.0, 193.0,  # Size
    96.0, 132.0, -78.0,   # Origin
    1.0, 1.0, 1.0,        # Spacing
    -1.0, 0.0, 0.0,       # Directions
    0.0, -1.0, 0.0,
    0.0, 0.0, 1.0,
])  # fmt:skip


def load_ants_h5(filename: Path) -> nt.base.TransformBase:
    """Load ANTs H5 files as a nitransforms TransformChain"""
    # Borrowed from https://github.com/feilong/process
    # process.resample.parse_combined_hdf5()
    #
    # Changes:
    #   * Tolerate a missing displacement field
    #   * Return the original affine without a round-trip
    #   * Always return a nitransforms TransformChain
    #
    # This should be upstreamed into nitransforms
    h = h5py.File(filename)
    xform = ITKCompositeH5.from_h5obj(h)

    # nt.Affine
    transforms = [nt.Affine(xform[0].to_ras())]

    if "2" not in h["TransformGroup"]:
        return transforms[0]

    transform2 = h["TransformGroup"]["2"]

    # Confirm these transformations are applicable
    if transform2["TransformType"][:][0] != b"DisplacementFieldTransform_float_3_3":
        msg = "Unknown transform type [2]\n"
        for i in h["TransformGroup"].keys():
            msg += f'[{i}]: {h["TransformGroup"][i]["TransformType"][:][0]}\n'
        raise ValueError(msg)

    fixed_params = transform2["TransformFixedParameters"][:]
    if not np.array_equal(fixed_params, FIXED_PARAMS):
        msg = "Unexpected fixed parameters\n"
        msg += f"Expected: {FIXED_PARAMS}\n"
        msg += f"Found: {fixed_params}"
        if not np.array_equal(fixed_params[6:], FIXED_PARAMS[6:]):
            raise ValueError(msg)
        warnings.warn(msg)

    shape = tuple(fixed_params[:3].astype(int))
    warp = h["TransformGroup"]["2"]["TransformParameters"][:]
    warp = warp.reshape((*shape, 3)).transpose(2, 1, 0, 3)
    warp *= np.array([-1, -1, 1])

    warp_affine = np.eye(4)
    warp_affine[:3, :3] = fixed_params[9:].reshape((3, 3))
    warp_affine[:3, 3] = fixed_params[3:6]
    lps_to_ras = np.eye(4) * np.array([-1, -1, 1, 1])
    warp_affine = lps_to_ras @ warp_affine
    if np.array_equal(fixed_params, FIXED_PARAMS):
        # Confirm that we construct the right affine when fixed parameters are known
        assert np.array_equal(
            warp_affine,
            np.array(
                [
                    [1.0, 0.0, 0.0, -96.0],
                    [0.0, 1.0, 0.0, -132.0],
                    [0.0, 0.0, 1.0, -78.0],
                    [0.0, 0.0, 0.0, 1.0],
                ]
            ),
        )
    transforms.insert(0, nt.DenseFieldTransform(nb.Nifti1Image(warp, warp_affine)))
    return nt.TransformChain(transforms)
