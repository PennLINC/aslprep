# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Handling functional connectvity."""
import nibabel as nb
import numpy as np
import pandas as pd
from nilearn.maskers import NiftiLabelsMasker
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.utils.filemanip import fname_presuffix

from aslprep import config


class _ParcellateCBFInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="File to be parcellated.")
    mask = File(exists=True, mandatory=True, desc="brain mask file")
    atlas = File(exists=True, mandatory=True, desc="atlas file")
    atlas_labels = File(exists=True, mandatory=True, desc="atlas labels file")
    min_coverage = traits.Float(
        default=0.5,
        usedefault=True,
        desc=(
            "Coverage threshold to apply to parcels. "
            "Any parcels with lower coverage than the threshold will be replaced with NaNs. "
            "Must be a value between zero and one. "
            "Default is 0.5."
        ),
    )


class _ParcellateCBFOutputSpec(TraitedSpec):
    timeseries = File(exists=True, mandatory=True, desc="Parcellated time series file.")
    coverage = File(exists=True, mandatory=True, desc="Parcel-wise coverage file.")


class ParcellateCBF(SimpleInterface):
    """Extract timeseries and compute connectivity matrices.

    Write out time series using Nilearn's NiftiLabelMasker
    Then write out functional correlation matrix of
    timeseries using numpy.
    """

    input_spec = _ParcellateCBFInputSpec
    output_spec = _ParcellateCBFOutputSpec

    def _run_interface(self, runtime):
        in_file = self.inputs.in_file
        mask = self.inputs.mask
        atlas = self.inputs.atlas
        atlas_labels = self.inputs.atlas_labels
        min_coverage = self.inputs.min_coverage

        node_labels_df = pd.read_table(atlas_labels, index_col="index")
        node_labels_df.sort_index(inplace=True)  # ensure index is in order

        # Explicitly remove label corresponding to background (index=0), if present.
        if 0 in node_labels_df.index:
            config.loggers.interface.warning(
                "Index value of 0 found in atlas labels file. "
                "Will assume this describes the background and ignore it."
            )
            node_labels_df = node_labels_df.drop(index=[0])

        node_labels = node_labels_df["label"].tolist()

        self._results["timeseries"] = fname_presuffix(
            "timeseries.tsv",
            newpath=runtime.cwd,
            use_ext=True,
        )
        self._results["coverage"] = fname_presuffix(
            "coverage.tsv",
            newpath=runtime.cwd,
            use_ext=True,
        )

        # Before anything, we need to measure coverage
        atlas_img = nb.load(atlas)
        atlas_data = atlas_img.get_fdata()
        atlas_data_bin = (atlas_data > 0).astype(np.float32)
        atlas_img_bin = nb.Nifti1Image(atlas_data_bin, atlas_img.affine, atlas_img.header)

        sum_masker_masked = NiftiLabelsMasker(
            labels_img=atlas,
            labels=node_labels,
            mask_img=mask,
            smoothing_fwhm=None,
            standardize=False,
            strategy="sum",
            resampling_target=None,  # they should be in the same space/resolution already
        )
        sum_masker_unmasked = NiftiLabelsMasker(
            labels_img=atlas,
            labels=node_labels,
            smoothing_fwhm=None,
            standardize=False,
            strategy="sum",
            resampling_target=None,  # they should be in the same space/resolution already
        )
        n_voxels_in_masked_parcels = sum_masker_masked.fit_transform(atlas_img_bin)
        n_voxels_in_parcels = sum_masker_unmasked.fit_transform(atlas_img_bin)
        parcel_coverage = np.squeeze(n_voxels_in_masked_parcels / n_voxels_in_parcels)
        coverage_thresholded = parcel_coverage < min_coverage

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
            config.loggers.interface.warning(
                f"{n_nodes - n_found_nodes}/{n_nodes} of parcels not found in atlas file."
            )

        if n_bad_nodes:
            config.loggers.interface.warning(
                f"{n_bad_nodes}/{n_nodes} of parcels have 0% coverage."
            )

        if n_poor_parcels:
            config.loggers.interface.warning(
                f"{n_poor_parcels}/{n_nodes} of parcels have <50% coverage. "
                "These parcels' time series will be replaced with zeros."
            )

        if n_partial_parcels:
            config.loggers.interface.warning(
                f"{n_partial_parcels}/{n_nodes} of parcels have at least one uncovered "
                "voxel, but have enough good voxels to be useable. "
                "The bad voxels will be ignored and the parcels' time series will be "
                "calculated from the remaining voxels."
            )

        masker = NiftiLabelsMasker(
            labels_img=atlas,
            labels=node_labels,
            mask_img=mask,
            smoothing_fwhm=None,
            standardize=False,
            resampling_target=None,  # they should be in the same space/resolution already
        )

        # Use nilearn for time_series
        timeseries_arr = masker.fit_transform(in_file)
        assert timeseries_arr.shape[1] == n_found_nodes

        # Apply the coverage mask
        timeseries_arr[:, coverage_thresholded] = np.nan

        # Region indices in the atlas may not be sequential, so we map them to sequential ints.
        seq_mapper = {idx: i for i, idx in enumerate(node_labels_df.index.tolist())}

        if n_found_nodes != n_nodes:  # parcels lost by warping/downsampling atlas
            # Fill in any missing nodes in the timeseries array with NaNs.
            new_timeseries_arr = np.full(
                (timeseries_arr.shape[0], n_nodes),
                fill_value=np.nan,
                dtype=timeseries_arr.dtype,
            )
            for col in range(timeseries_arr.shape[1]):
                label_col = seq_mapper[masker.labels_[col]]
                new_timeseries_arr[:, label_col] = timeseries_arr[:, col]

            timeseries_arr = new_timeseries_arr

            # Fill in any missing nodes in the coverage array with zero.
            new_parcel_coverage = np.zeros(n_nodes, dtype=parcel_coverage.dtype)
            for row in range(parcel_coverage.shape[0]):
                label_row = seq_mapper[masker.labels_[row]]
                new_parcel_coverage[label_row] = parcel_coverage[row]

            parcel_coverage = new_parcel_coverage

        # The time series file is tab-delimited, with node names included in the first row.
        timeseries_df = pd.DataFrame(data=timeseries_arr, columns=node_labels)
        coverage_df = pd.DataFrame(data=parcel_coverage, index=node_labels, columns=["coverage"])

        timeseries_df.to_csv(self._results["timeseries"], sep="\t", na_rep="n/a", index=False)
        coverage_df.to_csv(self._results["coverage"], sep="\t", index_label="Node")

        return runtime
