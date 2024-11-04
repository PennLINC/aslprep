# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Handling functional connectvity."""

import nibabel as nb
import numpy as np
import pandas as pd
from nilearn.maskers import NiftiLabelsMasker
from nipype import logging
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.utils.filemanip import fname_presuffix

LOGGER = logging.getLogger('nipype.interface')


class _ParcellateCBFInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='File to be parcellated.')
    mask = File(exists=True, mandatory=True, desc='brain mask file')
    atlas = File(exists=True, mandatory=True, desc='atlas file')
    atlas_labels = File(exists=True, mandatory=True, desc='atlas labels file')
    min_coverage = traits.Float(
        default=0.5,
        usedefault=True,
        desc=(
            'Coverage threshold to apply to parcels. '
            'Any parcels with lower coverage than the threshold will be replaced with NaNs. '
            'Must be a value between zero and one. '
            'Default is 0.5.'
        ),
    )


class _ParcellateCBFOutputSpec(TraitedSpec):
    timeseries = File(exists=True, desc='Parcellated time series file.')
    coverage = File(exists=True, desc='Parcel-wise coverage file.')


class ParcellateCBF(SimpleInterface):
    """Extract timeseries and compute connectivity matrices.

    Write out time series using Nilearn's NiftiLabelMasker
    Then write out functional correlation matrix of
    timeseries using numpy.
    """

    input_spec = _ParcellateCBFInputSpec
    output_spec = _ParcellateCBFOutputSpec

    def _run_interface(self, runtime):
        mask = self.inputs.mask
        atlas = self.inputs.atlas
        min_coverage = self.inputs.min_coverage

        node_labels_df = pd.read_table(self.inputs.atlas_labels, index_col='index')

        # Fix any nonsequential values or mismatch between atlas and DataFrame.
        atlas_img, node_labels_df = _sanitize_nifti_atlas(atlas, node_labels_df)
        node_labels = node_labels_df['label'].tolist()
        # prepend "background" to node labels to satisfy NiftiLabelsMasker
        # The background "label" won't be present in the output timeseries.
        masker_labels = ['background'] + node_labels

        # Before anything, we need to measure coverage
        atlas_img_bin = nb.Nifti1Image(
            (atlas_img.get_fdata() > 0).astype(np.uint8),
            atlas_img.affine,
            atlas_img.header,
        )

        sum_masker_masked = NiftiLabelsMasker(
            labels_img=atlas_img,
            labels=masker_labels,
            background_label=0,
            mask_img=mask,
            smoothing_fwhm=None,
            standardize=False,
            strategy='sum',
            resampling_target=None,  # they should be in the same space/resolution already
        )
        sum_masker_unmasked = NiftiLabelsMasker(
            labels_img=atlas_img,
            labels=masker_labels,
            background_label=0,
            smoothing_fwhm=None,
            standardize=False,
            strategy='sum',
            resampling_target=None,  # they should be in the same space/resolution already
        )
        n_voxels_in_masked_parcels = sum_masker_masked.fit_transform(atlas_img_bin)
        n_voxels_in_parcels = sum_masker_unmasked.fit_transform(atlas_img_bin)
        parcel_coverage = np.squeeze(n_voxels_in_masked_parcels / n_voxels_in_parcels)
        coverage_thresholded = parcel_coverage < min_coverage
        del sum_masker_masked, sum_masker_unmasked, n_voxels_in_masked_parcels, n_voxels_in_parcels

        n_nodes = len(node_labels)
        n_found_nodes = coverage_thresholded.size
        n_bad_nodes = np.sum(parcel_coverage == 0)
        n_poor_parcels = np.sum(
            np.logical_and(parcel_coverage > 0, parcel_coverage < min_coverage)
        )
        n_partial_parcels = np.sum(
            np.logical_and(parcel_coverage >= min_coverage, parcel_coverage < 1)
        )

        if n_found_nodes != n_nodes:
            LOGGER.warning(
                f'{n_nodes - n_found_nodes}/{n_nodes} of parcels not found in atlas file.'
            )

        if n_bad_nodes:
            LOGGER.warning(f'{n_bad_nodes}/{n_nodes} of parcels have 0% coverage.')

        if n_poor_parcels:
            LOGGER.warning(
                f'{n_poor_parcels}/{n_nodes} of parcels have <50% coverage. '
                "These parcels' time series will be replaced with zeros."
            )

        if n_partial_parcels:
            LOGGER.warning(
                f'{n_partial_parcels}/{n_nodes} of parcels have at least one uncovered '
                'voxel, but have enough good voxels to be usable. '
                "The bad voxels will be ignored and the parcels' time series will be "
                'calculated from the remaining voxels.'
            )

        masker = NiftiLabelsMasker(
            labels_img=atlas_img,
            labels=masker_labels,
            background_label=0,
            mask_img=mask,
            smoothing_fwhm=None,
            standardize=False,
            resampling_target=None,  # they should be in the same space/resolution already
        )

        # Use nilearn to parcellate the file
        timeseries_arr = masker.fit_transform(self.inputs.in_file)
        if timeseries_arr.shape[1] != n_found_nodes:
            raise ValueError(
                f'Expected {n_found_nodes} parcels, but found {timeseries_arr.shape[1]}.'
            )
        masker_labels = masker.labels_[:]
        del masker

        # Apply the coverage mask
        timeseries_arr[:, coverage_thresholded] = np.nan

        # Region indices in the atlas may not be sequential, so we map them to sequential ints.
        seq_mapper = {idx: i for i, idx in enumerate(node_labels_df['sanitized_index'].tolist())}

        if n_found_nodes != n_nodes:  # parcels lost by warping/downsampling atlas
            # Fill in any missing nodes in the timeseries array with NaNs.
            new_timeseries_arr = np.full(
                (timeseries_arr.shape[0], n_nodes),
                fill_value=np.nan,
                dtype=timeseries_arr.dtype,
            )
            for col in range(timeseries_arr.shape[1]):
                label_col = seq_mapper[masker_labels[col]]
                new_timeseries_arr[:, label_col] = timeseries_arr[:, col]

            timeseries_arr = new_timeseries_arr
            del new_timeseries_arr

            # Fill in any missing nodes in the coverage array with zero.
            new_parcel_coverage = np.zeros(n_nodes, dtype=parcel_coverage.dtype)
            for row in range(parcel_coverage.shape[0]):
                label_row = seq_mapper[masker_labels[row]]
                new_parcel_coverage[label_row] = parcel_coverage[row]

            parcel_coverage = new_parcel_coverage
            del new_parcel_coverage

        # The time series file is tab-delimited, with node names included in the first row.
        self._results['timeseries'] = fname_presuffix(
            'timeseries.tsv',
            newpath=runtime.cwd,
            use_ext=True,
        )
        timeseries_df = pd.DataFrame(data=timeseries_arr, columns=node_labels)
        timeseries_df.to_csv(self._results['timeseries'], sep='\t', na_rep='n/a', index=False)

        # Save out the coverage tsv
        coverage_df = pd.DataFrame(
            data=parcel_coverage.astype(np.float32),
            index=node_labels,
            columns=['coverage'],
        )
        self._results['coverage'] = fname_presuffix(
            'coverage.tsv',
            newpath=runtime.cwd,
            use_ext=True,
        )
        coverage_df.to_csv(self._results['coverage'], sep='\t', na_rep='n/a', index_label='Node')

        return runtime


def _sanitize_nifti_atlas(atlas, df):
    atlas_img = nb.load(atlas)
    atlas_data = atlas_img.get_fdata()
    atlas_data = atlas_data.astype(np.int16)

    # Check that all labels in the DataFrame are present in the NIfTI file, and vice versa.
    if 0 in df.index:
        df = df.drop(index=[0])

    df.sort_index(inplace=True)  # ensure index is in order
    expected_values = df.index.values

    found_values = np.unique(atlas_data)
    found_values = found_values[found_values != 0]  # drop the background value
    if not np.all(np.isin(found_values, expected_values)):
        raise ValueError('Atlas file contains values that are not present in the DataFrame.')

    # Map the labels in the DataFrame to sequential values.
    label_mapper = {value: i + 1 for i, value in enumerate(expected_values)}
    df['sanitized_index'] = [label_mapper[i] for i in df.index.values]

    # Map the values in the atlas image to sequential values.
    new_atlas_data = np.zeros(atlas_data.shape, dtype=np.int16)
    for old_value, new_value in label_mapper.items():
        new_atlas_data[atlas_data == old_value] = new_value

    new_atlas_img = nb.Nifti1Image(new_atlas_data, atlas_img.affine, atlas_img.header)

    return new_atlas_img, df
