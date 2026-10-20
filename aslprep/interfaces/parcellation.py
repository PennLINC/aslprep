# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Handling functional connectvity."""

import gc

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


class _NiftiParcellateInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='input file')
    mask = File(exists=True, mandatory=True, desc='brain mask file')
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
    coverage = File(exists=True, desc='Parcel-wise coverage file.')
    timeseries = File(exists=True, desc='Parcellated time series file.')


class NiftiParcellate(SimpleInterface):
    """Extract timeseries and compute connectivity matrices.

    Write out time series using Nilearn's NiftiLabelMasker
    """

    input_spec = _NiftiParcellateInputSpec
    output_spec = _NiftiParcellateOutputSpec

    def _run_interface(self, runtime):
        mask = self.inputs.mask
        atlas_img = nb.load(self.inputs.atlas)
        min_coverage = self.inputs.min_coverage

        node_labels_df = pd.read_table(self.inputs.atlas_labels)
        # The index tells us which row/column in the matrix the parcel is associated with.
        # The 'index' column tells us what that parcel's value in the atlas image is.
        # One requirement for later is that the index values are sorted in ascending order.
        node_labels_df = node_labels_df.sort_values(by='index').reset_index(drop=True)
        node_labels = node_labels_df['name'].tolist()
        # Create a dictionary mapping df['index'] to df.index
        full_parcel_mapper = {v: k for k, v in enumerate(node_labels_df['index'].tolist())}
        masker_lut = node_labels_df.copy()
        masker_lut = masker_lut[['index', 'name']]
        atlas_values = np.unique(atlas_img.get_fdata())
        atlas_values = atlas_values[atlas_values != 0]
        atlas_values = atlas_values.astype(int)
        masker_lut = masker_lut.loc[masker_lut['index'].isin(atlas_values)].reset_index(drop=True)

        if not np.any(nb.load(mask).get_fdata()):
            # Nothing was estimated (possible for maps with their own validity masks, such as
            # partial-volume corrected CBF). Nilearn rejects empty masks, so write the outputs
            # directly: no parcel has coverage or a value.
            LOGGER.warning(f'The mask {mask} is empty, so every parcel has 0% coverage.')
            in_shape = nb.load(self.inputs.in_file).shape
            n_volumes = in_shape[3] if len(in_shape) == 4 else 1
            self._write_outputs(
                runtime,
                timeseries_arr=np.full((n_volumes, len(node_labels)), np.nan),
                parcel_coverage=np.zeros(len(node_labels)),
                node_labels=node_labels,
            )
            return runtime

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
            lut=masker_lut,
            background_label=0,
            mask_img=mask,
            smoothing_fwhm=None,
            standardize=False,
            resampling_target=None,  # they should be in the same space/resolution already
            keep_masked_labels=True,
        )

        # Use nilearn to parcellate the file
        timeseries_arr = masker.fit_transform(self.inputs.in_file)
        if timeseries_arr.ndim == 1:
            # Add singleton first dimension representing time.
            timeseries_arr = timeseries_arr[None, :]

        if timeseries_arr.shape[1] != n_found_nodes:
            raise ValueError(
                f'timeseries_arr.shape[1] ({timeseries_arr.shape[1]}) != '
                f'n_found_nodes ({n_found_nodes})'
            )
        # Map from atlas value to column index for parcels found in the atlas image
        # Keys are cols/rows in the matrix, values are atlas values
        masker_parcel_mapper = masker.region_ids_
        # Remove 'background' label
        masker_parcel_mapper = {k: v for k, v in masker_parcel_mapper.items() if k != 'background'}
        del masker
        gc.collect()

        # Apply the coverage mask
        timeseries_arr[:, coverage_thresholded] = np.nan

        if n_found_nodes != n_nodes:  # parcels lost by warping/downsampling atlas
            # Fill in any missing nodes in the timeseries array with NaNs.
            new_timeseries_arr = np.full(
                (timeseries_arr.shape[0], n_nodes),
                fill_value=np.nan,
                dtype=timeseries_arr.dtype,
            )
            for col_num, node_value in masker_parcel_mapper.items():
                full_col_num = full_parcel_mapper[node_value]
                new_timeseries_arr[:, full_col_num] = timeseries_arr[:, col_num]

            timeseries_arr = new_timeseries_arr
            del new_timeseries_arr
            gc.collect()

            # Fill in any missing nodes in the coverage array with zero.
            new_parcel_coverage = np.zeros(n_nodes, dtype=parcel_coverage.dtype)
            for col_num, node_value in masker_parcel_mapper.items():
                full_col_num = full_parcel_mapper[node_value]
                new_parcel_coverage[full_col_num] = parcel_coverage[col_num]

            parcel_coverage = new_parcel_coverage
            del new_parcel_coverage
            gc.collect()

        self._write_outputs(runtime, timeseries_arr, parcel_coverage, node_labels)
        return runtime

    def _write_outputs(self, runtime, timeseries_arr, parcel_coverage, node_labels):
        # The time series file is tab-delimited, with node names included in the first row.
        self._results['timeseries'] = fname_presuffix(
            'timeseries.tsv',
            newpath=runtime.cwd,
            use_ext=True,
        )
        timeseries_df = pd.DataFrame(data=timeseries_arr, columns=node_labels)
        timeseries_df.to_csv(self._results['timeseries'], sep='\t', na_rep='n/a', index=False)

        # The coverage file has the same layout: one row, with node names as columns.
        coverage_df = pd.DataFrame(
            data=parcel_coverage.astype(np.float32)[None, :],
            columns=node_labels,
        )
        self._results['coverage'] = fname_presuffix(
            'coverage.tsv',
            newpath=runtime.cwd,
            use_ext=True,
        )
        coverage_df.to_csv(self._results['coverage'], sep='\t', na_rep='n/a', index=False)


class _StackParcellatedTablesInputSpec(BaseInterfaceInputSpec):
    in_files = traits.List(
        traits.List(File(exists=True)),
        mandatory=True,
        desc=(
            'One-row parcellated tables, from NiftiParcellate. '
            'The outer list is over estimates, and the inner list is over atlases.'
        ),
    )
    estimates = traits.List(
        traits.Str(),
        mandatory=True,
        desc='Label of each estimate in in_files, used as its column name.',
    )


class _StackParcellatedTablesOutputSpec(TraitedSpec):
    out_files = traits.List(
        File(exists=True),
        desc='One table per atlas, with one row per parcel and one column per estimate.',
    )


class StackParcellatedTables(SimpleInterface):
    """Combine one-row parcellated tables into one table per atlas.

    Each output table has a ``node`` column, with the name of the parcel in that row,
    followed by one column per estimate. This matches the tables of the BIDS Derivatives
    Toolbox (BDT).
    """

    input_spec = _StackParcellatedTablesInputSpec
    output_spec = _StackParcellatedTablesOutputSpec

    def _run_interface(self, runtime):
        in_files = self.inputs.in_files
        estimates = self.inputs.estimates
        if len(in_files) != len(estimates):
            raise ValueError(
                f'Got {len(in_files)} lists of tables but {len(estimates)} estimate labels.'
            )
        if len(set(estimates)) != len(estimates):
            raise ValueError(f'Estimate labels must be unique, got {estimates}.')

        n_atlases = {len(atlas_files) for atlas_files in in_files}
        if len(n_atlases) != 1:
            raise ValueError(f'Every estimate needs one table per atlas, got {n_atlases}.')

        self._results['out_files'] = []
        for i_atlas, atlas_files in enumerate(zip(*in_files, strict=True)):
            columns = {}
            for estimate, in_file in zip(estimates, atlas_files, strict=True):
                df = pd.read_table(in_file)
                if df.shape[0] != 1:
                    raise ValueError(f'{in_file} has {df.shape[0]} rows, not 1.')
                if columns and list(df.columns) != list(next(iter(columns.values())).index):
                    raise ValueError(f'{in_file} has different parcels than {atlas_files[0]}.')
                columns[estimate] = df.iloc[0]

            out_file = fname_presuffix(
                f'atlas{i_atlas:02d}.tsv',
                newpath=runtime.cwd,
                use_ext=True,
            )
            out_df = pd.DataFrame(columns)
            out_df.to_csv(out_file, sep='\t', na_rep='n/a', index_label='node')
            self._results['out_files'].append(out_file)

        return runtime
