# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Interfaces for parcellating maps with atlases."""

import gc

import nibabel as nb
import numpy as np
import pandas as pd
from nilearn.image import concat_imgs
from nilearn.maskers import NiftiLabelsMasker
from nipype import logging
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    InputMultiObject,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.utils.filemanip import fname_presuffix

LOGGER = logging.getLogger('nipype.interface')


class _NiftiParcellateInputSpec(BaseInterfaceInputSpec):
    in_files = InputMultiObject(
        File(exists=True),
        mandatory=True,
        desc='3D maps to parcellate.',
    )
    masks = InputMultiObject(
        File(exists=True),
        mandatory=True,
        desc=(
            'The mask of each map in in_files. '
            'Only voxels in the mask enter the parcel means and coverage.'
        ),
    )
    names = traits.List(
        traits.Str(),
        mandatory=True,
        desc='The column name of each map in in_files.',
    )
    atlas = File(exists=True, mandatory=True, desc='atlas file')
    atlas_labels = File(exists=True, mandatory=True, desc='atlas labels file')
    min_coverage = traits.Float(
        0.5,
        usedefault=True,
        desc=(
            'Coverage threshold to apply to parcels. '
            'Any parcels with lower coverage than the threshold will be replaced with NaNs. '
            'Must be a value between zero and one. '
            'Default is 0.5.'
        ),
    )


class _NiftiParcellateOutputSpec(TraitedSpec):
    scalarstats = File(
        exists=True,
        desc='Mean of each map (columns) in each parcel (rows).',
    )
    coverage = File(
        exists=True,
        desc='Fraction of each parcel (rows) within the mask of each map (columns).',
    )


class NiftiParcellate(SimpleInterface):
    """Parcellate 3D maps with an atlas, using Nilearn's NiftiLabelsMasker.

    Maps that share a mask are concatenated and parcellated together.
    The outputs have a ``node`` column, with the name of the parcel in that row,
    followed by one column per map. This matches the tables of the BIDS Derivatives
    Toolbox (BDT).
    """

    input_spec = _NiftiParcellateInputSpec
    output_spec = _NiftiParcellateOutputSpec

    def _run_interface(self, runtime):
        in_files = self.inputs.in_files
        masks = self.inputs.masks
        names = self.inputs.names
        if not len(in_files) == len(masks) == len(names):
            raise ValueError(
                f'Got {len(in_files)} maps, {len(masks)} masks and {len(names)} names.'
            )
        if len(set(names)) != len(names):
            raise ValueError(f'Map names must be unique, got {names}.')

        node_labels_df = pd.read_table(self.inputs.atlas_labels)
        # The index tells us which row/column in the matrix the parcel is associated with.
        # The 'index' column tells us what that parcel's value in the atlas image is.
        # One requirement for later is that the index values are sorted in ascending order.
        node_labels_df = node_labels_df.sort_values(by='index').reset_index(drop=True)
        node_labels = node_labels_df['name'].tolist()

        # Parcellate the maps that share a mask together.
        maps_by_mask = {}
        for i_map, mask in enumerate(masks):
            maps_by_mask.setdefault(mask, []).append(i_map)

        values = np.full((len(node_labels), len(names)), np.nan, dtype=np.float32)
        coverage = np.zeros((len(node_labels), len(names)), dtype=np.float32)
        atlas_img = nb.load(self.inputs.atlas)
        for mask, map_idx in maps_by_mask.items():
            group_values, group_coverage = _parcellate(
                concat_imgs([in_files[i] for i in map_idx], ensure_ndim=4),
                mask=mask,
                atlas_img=atlas_img,
                node_labels_df=node_labels_df,
                min_coverage=self.inputs.min_coverage,
            )
            values[:, map_idx] = group_values.T
            coverage[:, map_idx] = group_coverage[:, None]

        for output, data in (('scalarstats', values), ('coverage', coverage)):
            self._results[output] = fname_presuffix(
                f'{output}.tsv',
                newpath=runtime.cwd,
                use_ext=True,
            )
            pd.DataFrame(data, index=node_labels, columns=names).to_csv(
                self._results[output], sep='\t', na_rep='n/a', index_label='node'
            )

        return runtime


def _parcellate(in_img, mask, atlas_img, node_labels_df, min_coverage):
    """Parcellate the volumes of a 4D image within a mask.

    Returns
    -------
    values : numpy.ndarray of shape (n_volumes, n_parcels)
        Mean of each volume in each parcel, or NaN for parcels with less coverage than
        ``min_coverage``.
    parcel_coverage : numpy.ndarray of shape (n_parcels,)
        Fraction of each parcel within the mask.
    """
    node_labels = node_labels_df['name'].tolist()
    n_volumes = in_img.shape[3]
    # Create a dictionary mapping df['index'] to df.index
    full_parcel_mapper = {v: k for k, v in enumerate(node_labels_df['index'].tolist())}
    masker_lut = node_labels_df[['index', 'name']].copy()
    atlas_values = np.unique(atlas_img.get_fdata())
    atlas_values = atlas_values[atlas_values != 0]
    atlas_values = atlas_values.astype(int)
    masker_lut = masker_lut.loc[masker_lut['index'].isin(atlas_values)].reset_index(drop=True)

    if not np.any(nb.load(mask).get_fdata()):
        # Nothing was estimated (possible for maps with their own validity masks, such as
        # partial-volume corrected CBF). Nilearn rejects empty masks, so return directly:
        # no parcel has coverage or a value.
        LOGGER.warning(f'The mask {mask} is empty, so every parcel has 0% coverage.')
        return np.full((n_volumes, len(node_labels)), np.nan), np.zeros(len(node_labels))

    # Before anything, we need to measure coverage
    atlas_img_bin = nb.Nifti1Image(
        (atlas_img.get_fdata() > 0).astype(np.uint8),
        atlas_img.affine,
        atlas_img.header,
    )

    sum_masker_masked = NiftiLabelsMasker(
        labels_img=atlas_img,
        lut=masker_lut,
        background_label=0,
        mask_img=mask,
        smoothing_fwhm=None,
        standardize=False,
        strategy='sum',
        resampling_target=None,  # they should be in the same space/resolution already
        keep_masked_labels=True,
    )
    sum_masker_unmasked = NiftiLabelsMasker(
        labels_img=atlas_img,
        lut=masker_lut,
        background_label=0,
        smoothing_fwhm=None,
        standardize=False,
        strategy='sum',
        resampling_target=None,  # they should be in the same space/resolution already
        keep_masked_labels=True,
    )
    n_voxels_in_masked_parcels = sum_masker_masked.fit_transform(atlas_img_bin)
    n_voxels_in_parcels = sum_masker_unmasked.fit_transform(atlas_img_bin)
    parcel_coverage = np.squeeze(n_voxels_in_masked_parcels / n_voxels_in_parcels)
    coverage_thresholded = parcel_coverage < min_coverage
    del sum_masker_masked, sum_masker_unmasked, n_voxels_in_masked_parcels, n_voxels_in_parcels
    gc.collect()

    n_nodes = node_labels_df.shape[0]
    n_found_nodes = coverage_thresholded.size
    n_bad_nodes = np.sum(parcel_coverage == 0)
    n_poor_parcels = np.sum(np.logical_and(parcel_coverage > 0, parcel_coverage < min_coverage))
    n_partial_parcels = np.sum(
        np.logical_and(parcel_coverage >= min_coverage, parcel_coverage < 1)
    )

    if n_found_nodes != n_nodes:
        LOGGER.warning(f'{n_nodes - n_found_nodes}/{n_nodes} of parcels not found in atlas file.')

    if n_bad_nodes:
        LOGGER.warning(f'{n_bad_nodes}/{n_nodes} of parcels have 0% coverage.')

    if n_poor_parcels:
        LOGGER.warning(
            f'{n_poor_parcels}/{n_nodes} of parcels have <{min_coverage * 100}% coverage. '
            "These parcels' values will be set to n/a."
        )

    if n_partial_parcels:
        LOGGER.warning(
            f'{n_partial_parcels}/{n_nodes} of parcels have at least one uncovered '
            'voxel, but have enough good voxels to be usable. '
            "The bad voxels will be ignored and the parcels' values will be "
            'calculated from the remaining voxels.'
        )

    masker = NiftiLabelsMasker(
        labels_img=atlas_img,
        lut=masker_lut,
        background_label=0,
        mask_img=mask,
        smoothing_fwhm=None,
        standardize=False,
        resampling_target=None,  # they should be in the same space/resolution already
        keep_masked_labels=True,
    )

    # Use nilearn to parcellate the image
    values = masker.fit_transform(in_img)
    if values.ndim == 1:
        # Add singleton first dimension representing the volumes.
        values = values[None, :]

    if values.shape[1] != n_found_nodes:
        raise ValueError(f'values.shape[1] ({values.shape[1]}) != n_found_nodes ({n_found_nodes})')
    # Map from atlas value to column index for parcels found in the atlas image
    # Keys are cols/rows in the matrix, values are atlas values
    masker_parcel_mapper = masker.region_ids_
    # Remove 'background' label
    masker_parcel_mapper = {k: v for k, v in masker_parcel_mapper.items() if k != 'background'}
    del masker
    gc.collect()

    # Apply the coverage mask
    values[:, coverage_thresholded] = np.nan

    if n_found_nodes != n_nodes:  # parcels lost by warping/downsampling atlas
        # Fill in any missing nodes in the values array with NaNs,
        # and in the coverage array with zero.
        new_values = np.full((n_volumes, n_nodes), fill_value=np.nan, dtype=values.dtype)
        new_parcel_coverage = np.zeros(n_nodes, dtype=parcel_coverage.dtype)
        for col_num, node_value in masker_parcel_mapper.items():
            full_col_num = full_parcel_mapper[node_value]
            new_values[:, full_col_num] = values[:, col_num]
            new_parcel_coverage[full_col_num] = parcel_coverage[col_num]

        values = new_values
        parcel_coverage = new_parcel_coverage

    return values, parcel_coverage
