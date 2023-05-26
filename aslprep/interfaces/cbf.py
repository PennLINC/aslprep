"""Interfaces for calculating CBF."""
import os
from numbers import Number

import nibabel as nb
import numpy as np
import pandas as pd
from nibabel.processing import smooth_image
from nilearn import image, maskers
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    isdefined,
    traits,
)
from nipype.interfaces.fsl import MultiImageMaths
from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec
from nipype.utils.filemanip import fname_presuffix

from aslprep import config
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.utils.asl import (
    determine_multi_pld,
    estimate_labeling_efficiency,
    pcasl_or_pasl,
)
from aslprep.utils.cbf import (
    _score_cbf,
    _scrub_cbf,
    estimate_cbf_pcasl_multipld,
    estimate_t1,
)


class _RefineMaskInputSpec(BaseInterfaceInputSpec):
    t1w_mask = File(exists=True, mandatory=True, desc="t1 mask")
    asl_mask = File(exists=True, mandatory=True, desct="asl mask")
    transforms = File(exists=True, mandatory=True, desc="transfom")


class _RefineMaskOutputSpec(TraitedSpec):
    out_mask = File(exists=False, desc="output mask")
    out_tmp = File(exists=False, desc="tmp mask")


class RefineMask(SimpleInterface):
    """Reduce the ASL-derived brain mask using the associated T1w mask."""

    input_spec = _RefineMaskInputSpec
    output_spec = _RefineMaskOutputSpec

    def _run_interface(self, runtime):
        self._results["out_tmp"] = fname_presuffix(
            self.inputs.asl_mask,
            suffix="_tempmask",
            newpath=runtime.cwd,
        )
        self._results["out_mask"] = fname_presuffix(
            self.inputs.asl_mask,
            suffix="_refinemask",
            newpath=runtime.cwd,
        )

        refine_ref_mask(
            t1w_mask=self.inputs.t1w_mask,
            ref_asl_mask=self.inputs.asl_mask,
            t12ref_transform=self.inputs.transforms,
            tmp_mask=self._results["out_tmp"],
            refined_mask=self._results["out_mask"],
        )

        return runtime


class _ExtractCBFInputSpec(BaseInterfaceInputSpec):
    name_source = File(exists=True, mandatory=True, desc="raw asl file")
    asl_file = File(exists=True, mandatory=True, desc="preprocessed asl file")
    metadata = traits.Dict(mandatory=True, desc="metadata for ASL file")
    aslcontext = File(exists=True, mandatory=True, desc="aslcontext TSV file for run.")
    m0scan = traits.Either(
        File(exists=True),
        None,
        mandatory=True,
        desc="m0scan file associated with the ASL file. Only defined if M0Type is 'Separate'.",
    )
    m0scan_metadata = traits.Either(
        traits.Dict,
        None,
        mandatory=True,
        desc="metadata for M0 scan. Only defined if M0Type is 'Separate'.",
    )
    in_mask = File(exists=True, mandatory=True, desc="mask")
    dummy_vols = traits.Int(
        default_value=0,
        usedefault=True,
        mandatory=False,
        desc="remove first n volumes",
    )
    fwhm = traits.Float(default_value=5, usedefault=True, mandatory=False, desc="fwhm")


class _ExtractCBFOutputSpec(TraitedSpec):
    out_file = File(exists=False, desc="Either CBF or deltaM time series.")
    m0_file = File(exists=False, desc="Mean M0 image, after smoothing.")
    metadata = traits.Dict(
        desc=(
            "Metadata for the ASL run. "
            "The dictionary may be modified to only include metadata associated with the selected "
            "volumes."
        ),
    )
    m0tr = traits.Either(
        traits.Float,
        None,
        desc="RepetitionTimePreparation for M0 scans.",
    )


class ExtractCBF(SimpleInterface):
    """Extract CBF time series by subtracting label volumes from control volumes.

    TODO: Mock up test data and write tests to cover all of the branches in this interface.

    Notes
    -----
    The M0 information is extracted in the same way as GeReferenceFile,
    so there's duplication that could be reduced.
    """

    input_spec = _ExtractCBFInputSpec
    output_spec = _ExtractCBFOutputSpec

    def _run_interface(self, runtime):
        aslcontext = pd.read_table(self.inputs.aslcontext)
        metadata = self.inputs.metadata.copy()

        mask_data = nb.load(self.inputs.in_mask).get_fdata()

        # read the preprocessed ASL data
        asl_img = nb.load(self.inputs.asl_file)
        asl_data = asl_img.get_fdata()

        # get the control, tag, moscan or label
        vol_types = aslcontext["volume_type"].tolist()
        control_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "control"]
        label_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "label"]
        m0_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "m0scan"]
        deltam_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "deltam"]
        cbf_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "cbf"]

        # extract m0 file and register it to ASL if separate
        if metadata["M0Type"] == "Separate":
            m0file = self.inputs.m0scan

            m0_in_asl = fname_presuffix(self.inputs.asl_file, suffix="_m0file")
            m0_in_asl = regmotoasl(asl=self.inputs.asl_file, m0file=m0file, m02asl=m0_in_asl)
            m0data_smooth = smooth_image(nb.load(m0_in_asl), fwhm=self.inputs.fwhm).get_fdata()
            if len(m0data_smooth.shape) > 3:
                m0data = mask_data * np.mean(m0data_smooth, axis=3)
            else:
                m0data = mask_data * m0data_smooth

            m0tr = self.inputs.m0scan_metadata["RepetitionTimePreparation"]
            if np.array(m0tr).size > 1 and np.std(m0tr) > 0:
                raise ValueError("M0 scans have variable TR. ASLPrep does not support this.")

        elif metadata["M0Type"] == "Included":
            m0data = asl_data[:, :, :, m0_volume_idx]
            m0img = nb.Nifti1Image(m0data, asl_img.affine, asl_img.header)
            m0data_smooth = smooth_image(m0img, fwhm=self.inputs.fwhm).get_fdata()
            m0data = mask_data * np.mean(m0data_smooth, axis=3)

            if np.array(metadata["RepetitionTimePreparation"]).size > 1:
                m0tr = np.array(metadata["RepetitionTimePreparation"])[m0_volume_idx]
            else:
                m0tr = metadata["RepetitionTimePreparation"]

            if np.array(m0tr).size > 1 and np.std(m0tr) > 0:
                raise ValueError("M0 scans have variable TR. ASLPrep does not support this.")

        elif metadata["M0Type"] == "Estimate":
            m0data = metadata["M0Estimate"] * mask_data

            m0tr = None

        elif metadata["M0Type"] == "Absent":
            if control_volume_idx:
                # Estimate M0 using the smoothed mean control volumes.
                control_data = asl_data[:, :, :, control_volume_idx]
                control_img = nb.Nifti1Image(control_data, asl_img.affine, asl_img.header)
                control_img = smooth_image(control_img, fwhm=self.inputs.fwhm).get_fdata()
                m0data = mask_data * np.mean(control_img, axis=3)

            elif cbf_volume_idx:
                # If we have precalculated CBF data, we don't need M0, so we'll just use the mask.
                m0data = mask_data

            else:
                raise RuntimeError(
                    "m0scan is absent, "
                    "and there are no control volumes that can be used as a substitute"
                )

            m0tr = None

        else:
            raise RuntimeError("no pathway to m0scan")

        if deltam_volume_idx:
            config.loggers.interface.info("Extracting deltaM from ASL file.")
            metadata_idx = deltam_volume_idx
            out_data = asl_data[:, :, :, deltam_volume_idx]

        elif label_volume_idx:
            config.loggers.interface.info(
                "Calculating deltaM from label-control pairs in ASL file."
            )
            assert len(label_volume_idx) == len(control_volume_idx)
            metadata_idx = control_volume_idx
            control_data = asl_data[:, :, :, control_volume_idx]
            label_data = asl_data[:, :, :, label_volume_idx]
            out_data = control_data - label_data

        elif cbf_volume_idx:
            metadata_idx = cbf_volume_idx
            out_data = asl_data[:, :, :, cbf_volume_idx]

        else:
            raise RuntimeError("No valid ASL or CBF image.")

        # Update the metadata as necessary
        VOLUME_WISE_FIELDS = [
            "PostLabelingDelay",
            "VascularCrushingVENC",
            "LabelingDuration",
            "EchoTime",
            "FlipAngle",
            "RepetitionTimePreparation",
        ]

        for field in VOLUME_WISE_FIELDS:
            if field not in metadata:
                continue

            value = metadata[field]
            if isinstance(value, list) and len(value) != asl_data.shape[3]:
                raise ValueError(
                    f"{field} is an array, but the number of values ({len(value)}) "
                    f"does not match the number of volumes in the ASL data ({asl_data.shape[3]})."
                )
            elif isinstance(value, list):
                # Reduce to only the selected volumes
                value = [value[i] for i in metadata_idx]

                # Remove dummy volumes as well
                if self.inputs.dummy_vols != 0:
                    value = value[self.inputs.dummy_vols :]

                metadata[field] = value

        self._results["metadata"] = metadata
        self._results["m0tr"] = m0tr
        self._results["out_file"] = fname_presuffix(
            self.inputs.name_source,
            suffix="_DeltaMOrCBF",
            newpath=runtime.cwd,
        )
        self._results["m0_file"] = fname_presuffix(
            self.inputs.name_source,
            suffix="_m0file",
            newpath=runtime.cwd,
        )
        nb.Nifti1Image(out_data, asl_img.affine, asl_img.header).to_filename(
            self._results["out_file"]
        )
        nb.Nifti1Image(m0data, asl_img.affine, asl_img.header).to_filename(
            self._results["m0_file"]
        )

        return runtime


class _ExtractCBForDeltaMInputSpec(BaseInterfaceInputSpec):
    asl_file = File(exists=True, mandatory=True, desc="raw asl file")
    aslcontext = File(exists=True, mandatory=True, desc="aslcontext TSV file for run.")
    asl_mask = File(exists=True, mandatory=True, desct="asl mask")
    file_type = traits.Str(desc="file type, c for cbf, d for deltam", mandatory=True)


class _ExtractCBForDeltaMOutputSpec(TraitedSpec):
    out_file = File(exists=False, desc="cbf or deltam")


class ExtractCBForDeltaM(SimpleInterface):
    """Load an ASL file and grab the CBF or DeltaM volumes from it."""

    input_spec = _ExtractCBForDeltaMInputSpec
    output_spec = _ExtractCBForDeltaMOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = fname_presuffix(
            self.inputs.asl_mask,
            suffix="_cbfdeltam",
            newpath=runtime.cwd,
        )
        asl_img = nb.load(self.inputs.asl_file)
        asl_data = asl_img.get_fdata()

        aslcontext = pd.read_table(self.inputs.aslcontext)
        vol_types = aslcontext["volume_type"].tolist()
        control_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "control"]
        label_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "label"]
        deltam_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "deltam"]
        cbf_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == "cbf"]

        if self.inputs.file_type == "d":
            if len(control_volume_idx) > 0:
                # Grab control and label volumes from ASL file,
                # then calculate deltaM by subtracting label volumes from control volumes.
                deltam_data = (
                    asl_data[:, :, :, control_volume_idx] - asl_data[:, :, :, label_volume_idx]
                )
                out_img = nb.Nifti1Image(
                    dataobj=deltam_data,
                    affine=asl_img.affine,
                    header=asl_img.header,
                )
            else:
                # Grab deltaM volumes from ASL file.
                if len(asl_data.shape) < 4:
                    # 3D volume is written out without any changes.
                    # NOTE: Why not return the original file then?
                    out_img = nb.Nifti1Image(
                        dataobj=asl_data,
                        affine=asl_img.affine,
                        header=asl_img.header,
                    )
                else:
                    deltam_data = asl_data[:, :, :, deltam_volume_idx]
                    out_img = nb.Nifti1Image(
                        dataobj=deltam_data,
                        affine=asl_img.affine,
                        header=asl_img.header,
                    )

        elif self.inputs.file_type == "c":
            if len(asl_data.shape) < 4:
                # 3D volume is written out without any changes.
                # NOTE: Why not return the original file then?
                out_img = nb.Nifti1Image(
                    dataobj=asl_data,
                    affine=asl_img.affine,
                    header=asl_img.header,
                )
            else:
                # Grab CBF volumes from ASL file.
                cbf_data = asl_data[:, :, :, cbf_volume_idx]
                out_img = nb.Nifti1Image(
                    dataobj=cbf_data,
                    affine=asl_img.affine,
                    header=asl_img.header,
                )

        out_img.to_filename(self._results["out_file"])

        return runtime


class _ComputeCBFInputSpec(BaseInterfaceInputSpec):
    deltam = File(
        exists=True,
        mandatory=True,
        desc=(
            "NIfTI file containing raw CBF volume(s). "
            "These raw CBF values are the result of subtracting label volumes from "
            "control volumes, without any kind of additional scaling. "
            "This file may be 3D or 4D."
        ),
    )
    metadata = traits.Dict(
        mandatory=True,
        desc="Metadata for the raw CBF file, taken from the raw ASL data's sidecar JSON file.",
    )
    m0_scale = traits.Float(
        mandatory=True,
        desc="Relative scale between ASL and M0.",
    )
    m0_file = File(exists=True, mandatory=True, desc="M0 nifti file")
    mask = File(exists=True, mandatory=True, desc="Mask nifti file")
    cbf_only = traits.Bool(
        mandatory=True,
        desc="Whether data are deltam (False) or CBF (True).",
    )


class _ComputeCBFOutputSpec(TraitedSpec):
    cbf_ts = traits.Either(
        File(exists=True),
        None,
        desc="Quantitative CBF time series, in mL/100g/min. Only generated for single-delay data.",
    )
    mean_cbf = File(exists=True, desc="Quantified CBF, averaged over time.")
    att = traits.Either(
        File(exists=True),
        None,
        desc="Arterial transit time map, in seconds. Only generated for multi-delay data.",
    )


class ComputeCBF(SimpleInterface):
    """Calculate CBF time series and mean control.

    Notes
    -----
    This interface calculates CBF from deltam and M0 data.
    It can handle single-delay and multi-delay data, single CBF volumes and CBF time series,
    and PASL and (P)CASL data.

    Single-delay CBF, for both (P)CASL and QUIPSSII PASL
    is calculated according to :footcite:t:`alsop_recommended_2015`.
    Multi-delay CBF is handled using a weighted average,
    based on :footcite:t:`dai2012reduced,wang2013multi`.

    Multi-delay CBF is calculated according to :footcite:t:`fan2017long`,
    although CBF is averaged across PLDs according to the method in
    :footcite:t:`juttukonda2021characterizing`.
    Arterial transit time is estimated according to :footcite:t:`dai2012reduced`.

    If slice timing information is detected, then PLDs will be shifted by the slice times.

    See Also
    --------
    :func:`~aslprep.utils.asl.pcasl_or_pasl`
    :func:`~aslprep.utils.asl.determine_multi_pld`
    :func:`~aslprep.utils.cbf.estimate_t1`
    :func:`~aslprep.utils.asl.estimate_labeling_efficiency`
    :func:`~aslprep.utils.cbf.estimate_cbf_pcasl_multipld`

    References
    ----------
    .. footbibliography::
    """

    input_spec = _ComputeCBFInputSpec
    output_spec = _ComputeCBFOutputSpec

    def _run_interface(self, runtime):
        metadata = self.inputs.metadata
        m0_file = self.inputs.m0_file
        m0_scale = self.inputs.m0_scale
        mask_file = self.inputs.mask
        deltam_file = self.inputs.deltam  # control - label signal intensities

        if self.inputs.cbf_only:
            config.loggers.interface.debug("CBF data detected. Skipping CBF estimation.")
            self._results["cbf_ts"] = fname_presuffix(
                deltam_file,
                suffix="_cbf_ts",
                newpath=runtime.cwd,
            )
            cbf_img = nb.load(deltam_file)
            cbf_img.to_filename(self._results["cbf_ts"])
            self._results["mean_cbf"] = fname_presuffix(
                deltam_file,
                suffix="_meancbf",
                newpath=runtime.cwd,
            )
            mean_cbf_img = image.mean_img(cbf_img)
            mean_cbf_img.to_filename(self._results["mean_cbf"])

            # No ATT available for pre-calculated CBF
            self._results["att"] = None

            return runtime

        is_casl = pcasl_or_pasl(metadata=metadata)
        is_multi_pld = determine_multi_pld(metadata=metadata)
        t1blood, t1tissue = estimate_t1(metadata=metadata)

        # PostLabelingDelay is either a single number or an array of numbers.
        # If it is an array of numbers, then there should be one value for every volume in the
        # time series, with any M0 volumes having a value of 0.
        plds = np.atleast_1d(metadata["PostLabelingDelay"])

        # Get labeling efficiency (alpha in Alsop 2015).
        labeleff = estimate_labeling_efficiency(metadata=metadata)

        UNIT_CONV = 6000  # convert units from mL/g/s to mL/(100 g)/min
        PARTITION_COEF = 0.9  # brain partition coefficient (lambda in Alsop 2015)

        # NOTE: Nilearn will still add a singleton time dimension for 3D imgs with
        # NiftiMasker.transform, until 0.12.0, so the arrays will currently be 2D no matter what.
        masker = maskers.NiftiMasker(mask_img=mask_file)
        deltam_arr = masker.fit_transform(deltam_file).T  # Transpose to SxT
        assert deltam_arr.ndim == 2, f"deltam is {deltam_arr.ndim}"

        # Load the M0 map and average over time, in case there's more than one map in the file.
        m0data = masker.transform(m0_file)
        m0data = np.mean(m0data, axis=0)
        scaled_m0data = m0_scale * m0data

        if "SliceTiming" in metadata:
            # Offset PLD(s) by slice times
            # This step builds a voxel-wise array of post-labeling delay values,
            # where voxels from each slice have the appropriately-shifted PLD value.
            # If there are multiple PLDs, then the second dimension of the PLD array will
            # correspond to volumes in the time series.
            config.loggers.interface.info(
                "2D acquisition with slice timing information detected. "
                "Shifting post-labeling delay values across the brain by slice times."
            )
            slice_times = np.array(metadata["SliceTiming"])

            # Determine which axis slices come from.
            # ASL data typically acquires along z axis, from inferior to superior.
            slice_encoding_direction = metadata.get("SliceEncodingDirection", "k")
            slice_encoding_axis = "ijk".index(slice_encoding_direction[0])

            deltam_img = nb.load(deltam_file)
            shape = deltam_img.shape[:3]
            if slice_times.size != shape[slice_encoding_axis]:
                raise ValueError(
                    f"Number of slices ({shape[slice_encoding_axis]}) != "
                    f"slice times ({slice_times.size})"
                )

            # Reverse the slice times if slices go from maximum index to zero.
            # This probably won't occur with ASL data though, since I --> S makes more sense than
            # S --> I.
            if slice_encoding_direction.endswith("-"):
                slice_times = slice_times[::-1]

            # Determine which dimensions to add to the slice times array,
            # so that all 4 dims are 1, except for the slice encoding axis,
            # which will have the slice times.
            new_dims = [0, 1, 2, 3]
            new_dims.pop(slice_encoding_axis)
            slice_times = np.expand_dims(slice_times, new_dims)

            # Create a 4D array of PLDs, matching shape of ASL data (except only one volume).
            pld_brain = np.tile(plds, list(shape) + [1])

            # Shift the PLDs by the appropriate slice times.
            pld_brain = pld_brain + slice_times

            # Mask the PLD array to go from (X, Y, Z, delay) to (S, delay)
            pld_img = nb.Nifti1Image(pld_brain, deltam_img.affine, deltam_img.header)
            plds = masker.transform(pld_img).T

            # Write out the slice-shifted PLDs to the working directory, for debugging.
            pld_file = fname_presuffix(
                deltam_file,
                suffix="_plds",
                newpath=runtime.cwd,
            )
            pld_img.to_filename(pld_file)

        elif is_multi_pld:
            # Broadcast PLDs to voxels by PLDs
            plds = np.dot(plds[:, None], np.ones((1, deltam_arr.shape[0]))).T

        if is_casl:
            tau = np.array(metadata["LabelingDuration"])

        if is_multi_pld:
            if is_casl:
                att, mean_cbf = estimate_cbf_pcasl_multipld(
                    deltam_arr,
                    scaled_m0data,
                    plds,
                    tau,
                    labeleff,
                    t1blood=t1blood,
                    t1tissue=t1tissue,
                    unit_conversion=UNIT_CONV,
                    partition_coefficient=PARTITION_COEF,
                )

            else:
                # Dai's approach can't be used on PASL data, so we'll need another method.
                raise ValueError(
                    "Multi-delay data are not supported for PASL sequences at the moment."
                )

            mean_cbf_img = masker.inverse_transform(mean_cbf)
            att_img = masker.inverse_transform(att)

            # Multi-delay data won't produce a CBF time series
            self._results["cbf_ts"] = None
            self._results["att"] = fname_presuffix(
                self.inputs.deltam,
                suffix="_att",
                newpath=runtime.cwd,
            )
            att_img.to_filename(self._results["att"])

        else:  # Single-delay
            if is_casl:
                denom_factor = t1blood * (1 - np.exp(-(tau / t1blood)))

            elif not metadata["BolusCutOffFlag"]:
                raise ValueError(
                    "PASL without a bolus cut-off technique is not supported in ASLPrep."
                )

            elif metadata["BolusCutOffTechnique"] == "QUIPSS":
                # PASL + QUIPSS
                # Only one BolusCutOffDelayTime allowed.
                assert isinstance(metadata["BolusCutOffDelayTime"], Number)
                denom_factor = plds - metadata["BolusCutOffDelayTime"]  # delta_TI, per Wong 1998

            elif metadata["BolusCutOffTechnique"] == "QUIPSSII":
                # PASL + QUIPSSII
                # Per SD, use PLD as TI for PASL, so we will just use 'plds' in the numerator when
                # calculating the perfusion factor.
                # Only one BolusCutOffDelayTime allowed.
                assert isinstance(metadata["BolusCutOffDelayTime"], Number)
                denom_factor = metadata["BolusCutOffDelayTime"]  # called TI1 in Alsop 2015

            elif metadata["BolusCutOffTechnique"] == "Q2TIPS":
                # PASL + Q2TIPS
                # Q2TIPS should have two BolusCutOffDelayTimes.
                assert len(metadata["BolusCutOffDelayTime"]) == 2
                denom_factor = metadata["BolusCutOffDelayTime"][0]  # called TI1 in Noguchi 2015

            else:
                raise ValueError(
                    f"Unknown BolusCutOffTechnique {metadata['BolusCutOffTechnique']}"
                )

            # Q2TIPS uses TI2 instead of w (PLD), see Noguchi 2015 for this info.
            exp_numerator = (
                metadata["BolusCutOffDelayTime"][1]
                if metadata.get("BolusCutOffTechnique") == "Q2TIPS"
                else plds
            )

            # Scale difference signal to absolute CBF units by dividing by PD image (M0 * M0scale).
            deltam_scaled = deltam_arr / scaled_m0data[:, None]

            perfusion_factor = (UNIT_CONV * PARTITION_COEF * np.exp(exp_numerator / t1blood)) / (
                denom_factor * 2 * labeleff
            )

            cbf_ts = deltam_scaled * perfusion_factor
            cbf_ts = np.nan_to_num(cbf_ts)

            cbf_ts_img = masker.inverse_transform(cbf_ts.T)
            mean_cbf_img = image.mean_img(cbf_ts_img)
            self._results["cbf_ts"] = fname_presuffix(
                self.inputs.deltam,
                suffix="_cbf",
                newpath=runtime.cwd,
            )
            cbf_ts_img.to_filename(self._results["cbf_ts"])
            # Single-delay data won't produce an ATT image
            self._results["att"] = None

        # Mean CBF is returned no matter what
        self._results["mean_cbf"] = fname_presuffix(
            self.inputs.deltam,
            suffix="_meancbf",
            newpath=runtime.cwd,
        )
        mean_cbf_img.to_filename(self._results["mean_cbf"])

        return runtime


class _ScoreAndScrubCBFInputSpec(BaseInterfaceInputSpec):
    cbf_ts = File(exists=True, mandatory=True, desc="computed CBF from ComputeCBF")
    gm_tpm = File(exists=True, mandatory=True, desc="Gray matter tissue probability map.")
    wm_tpm = File(exists=True, mandatory=True, desc="White matter tissue probability map.")
    csf_tpm = File(exists=True, mandatory=True, desc="Cerebrospinal fluid tissue probability map.")
    brain_mask = File(exists=True, mandatory=True, desc="Brain mask.")
    tpm_threshold = traits.Float(
        default_value=0.7,
        mandatory=False,
        usedefault=True,
        desc="Threshold for tissue probability maps.",
    )
    cost_function = traits.Str(
        mandatory=False,
        default_value="huber",
        option=["bisquare", "andrews", "cauchy", "fair", "logistics", "ols", "talwar", "welsch"],
        desc="Wavelet function used in SCRUB.",
        usedefault=True,
    )


class _ScoreAndScrubCBFOutputSpec(TraitedSpec):
    cbf_ts_score = File(
        exists=True,
        desc="CBF time series after removing outlier volumes flagged by SCORE algorithm.",
    )
    score_outlier_index = File(exists=True, desc="Index of removed volumes, in a CSV file.")
    mean_cbf_score = File(
        exists=True,
        desc="Mean CBF image calculated from SCORE-censored CBF time series.",
    )
    mean_cbf_scrub = File(exists=True, desc="average scrub")


class ScoreAndScrubCBF(SimpleInterface):
    """Apply the SCORE and SCRUB algorithms.

    The Structural Correlation-based Outlier Rejection (SCORE) algorithm is applied to the CBF
    time series to discard CBF volumes with outlying values :footcite:p:`dolui2017structural`
    before computing the mean CBF.
    The Structural Correlation with RobUst Bayesian (SCRUB) algorithm is then applied to the CBF
    maps using structural tissue probability maps to reweight the mean CBF
    :footcite:p:`dolui2016scrub`.

    References
    ----------
    .. footbibliography::
    """

    input_spec = _ScoreAndScrubCBFInputSpec
    output_spec = _ScoreAndScrubCBFOutputSpec

    def _run_interface(self, runtime):
        cbf_img = nb.load(self.inputs.cbf_ts)
        cbf_arr = cbf_img.get_fdata()

        if cbf_arr.ndim > 3:
            mask_arr = nb.load(self.inputs.brain_mask).get_fdata()
            gm_tpm_arr = nb.load(self.inputs.gm_tpm).get_fdata()
            wm_tpm_arr = nb.load(self.inputs.wm_tpm).get_fdata()
            csf_tpm_arr = nb.load(self.inputs.csf_tpm).get_fdata()

            cbf_ts_score, score_outliers_idx = _score_cbf(
                cbf_ts=cbf_arr,
                gm=gm_tpm_arr,
                wm=wm_tpm_arr,
                csf=csf_tpm_arr,
                mask=mask_arr,
                thresh=self.inputs.tpm_threshold,
            )
            mean_cbf_score = np.mean(cbf_ts_score, axis=3)
            mean_cbf_scrub = _scrub_cbf(
                cbf_ts=cbf_ts_score,
                gm=gm_tpm_arr,
                wm=wm_tpm_arr,
                csf=csf_tpm_arr,
                mask=mask_arr,
                cost_function=self.inputs.cost_function,
                thresh=self.inputs.tpm_threshold,
            )
        else:
            config.loggers.interface.warning(
                f"CBF time series is only {cbf_arr.ndim}D. Skipping SCORE and SCRUB."
            )
            cbf_ts_score = cbf_arr
            mean_cbf_score = cbf_arr
            score_outliers_idx = np.array([0])
            mean_cbf_scrub = cbf_arr

        self._results["cbf_ts_score"] = fname_presuffix(
            self.inputs.cbf_ts,
            suffix="_cbfscorets",
            newpath=runtime.cwd,
        )
        self._results["mean_cbf_score"] = fname_presuffix(
            self.inputs.cbf_ts,
            suffix="_meancbfscore",
            newpath=runtime.cwd,
        )
        self._results["mean_cbf_scrub"] = fname_presuffix(
            self.inputs.cbf_ts,
            suffix="_cbfscrub",
            newpath=runtime.cwd,
        )
        self._results["score_outlier_index"] = fname_presuffix(
            self.inputs.cbf_ts,
            suffix="_scoreindex.csv",
            newpath=runtime.cwd,
            use_ext=False,
        )

        nb.Nifti1Image(
            dataobj=cbf_ts_score,
            affine=cbf_img.affine,
            header=cbf_img.header,
        ).to_filename(self._results["cbf_ts_score"])
        nb.Nifti1Image(
            dataobj=mean_cbf_score,
            affine=cbf_img.affine,
            header=cbf_img.header,
        ).to_filename(self._results["mean_cbf_score"])
        np.savetxt(self._results["score_outlier_index"], score_outliers_idx, delimiter=",")
        nb.Nifti1Image(
            dataobj=mean_cbf_scrub,
            affine=cbf_img.affine,
            header=cbf_img.header,
        ).to_filename(self._results["mean_cbf_scrub"])

        return runtime


class _BASILCBFInputSpec(FSLCommandInputSpec):
    # We use position args here as list indices - so a negative number
    # will put something on the end
    deltam = File(
        exists=True,
        desc=(
            "ASL data after subtracting tag-control or control-tag. "
            "This matches with ``--iaf diff``, which is the default."
        ),
        argstr="-i %s",
        position=0,
        mandatory=True,
    )
    mask = File(
        exists=True,
        argstr="-m %s",
        desc="mask in the same space as deltam",
        mandatory=True,
    )
    mzero = File(exists=True, argstr="-c %s", desc="m0 scan", mandatory=False)
    m0_scale = traits.Float(desc="calibration of asl", argstr="--cgain %.2f", mandatory=True)
    m0tr = traits.Float(
        desc="The repetition time for the calibration image (the M0 scan).",
        argstr="--tr %.2f",
        mandatory=True,
    )
    tis = traits.Either(
        traits.Float(),
        traits.List(traits.Float()),
        desc=(
            "The list of inflow times (TIs), a comma separated list of values should be provided "
            "(that matches the order in the data).\n\n"
            "Note, the inflow time is the PLD plus bolus duration for pcASL (and cASL), "
            "it equals the inversion time for pASL. "
            "If the data contains multiple repeats of the same set of TIs then it is only "
            "necessary to list the unique TIs.\n\n"
            "When using the ``--tis=`` you can specify a full list of all TIs/PLDs in the data "
            "(i.e., as many entries as there are label-control pairs). "
            "Or, if you have a number of TIs/PLDs repeated multiple times you can just list the "
            "unique TIs in order and ``oxford_asl`` will automatically replicate that list to "
            "match the number of repeated measurements in the data. "
            "If you have a variable number of repeats at each TI/PLD then either list all TIs "
            "or use the ``--rpts=<csv>`` option (see below)."
        ),
        argstr="--tis %s",
        mandatory=True,
        sep=",",
    )
    pcasl = traits.Bool(
        desc=(
            "Data were acquired using cASL or pcASL labelling "
            "(pASL labeling is assumed by default)."
        ),
        argstr="--casl",
        mandatory=False,
        default_value=False,
    )
    bolus = traits.Either(
        traits.Float(),
        traits.List(traits.Float()),
        desc="bolus or tau: label duration",
        argstr="--bolus %s",
        mandatory=True,
        sep=",",
    )
    slice_spacing = traits.Float(
        desc="Slice times",
        argstr="--slicedt %s",
        mandatory=False,
    )
    pvc = traits.Bool(
        desc="Do partial volume correction.",
        mandatory=False,
        argstr="--pvcorr",
        default_value=True,
        usedefault=True,
    )
    gm_tpm = File(
        exists=True,
        mandatory=False,
        desc="Partial volume estimates for GM. This is just a GM tissue probability map.",
        argstr="--pvgm %s",
    )
    wm_tpm = File(
        exists=True,
        mandatory=False,
        desc="Partial volume estimates for WM. This is just a WM tissue probability map.",
        argstr="--pvwm %s",
    )
    alpha = traits.Float(
        desc=(
            "Inversion efficiency - [default: 0.98 (pASL); 0.85 (cASL)]. "
            "This is equivalent to the BIDS metadata field 'LabelingEfficiency'."
        ),
        argstr="--alpha %.2f",
    )
    out_basename = File(desc="base name of output files", argstr="-o %s", mandatory=True)


class _BASILCBFOutputSpec(TraitedSpec):
    mean_cbf_basil = File(exists=True, desc="CBF map in absolute units (ml/100g/min).")
    mean_cbf_gm_basil = File(
        exists=True,
        desc=(
            "CBF map with gray matter partial volume correction. "
            "This means that the map contains 'pure' GM perfusion estimates, in absolute units."
        ),
    )
    mean_cbf_wm_basil = File(
        exists=True,
        desc=(
            "CBF map with white matter partial volume correction. "
            "This means that the map contains 'pure' WM perfusion estimates, in absolute units."
        ),
    )
    att_basil = File(exists=True, desc="Arterial transit time map.")


class BASILCBF(FSLCommand):
    """Apply Bayesian Inference for Arterial Spin Labeling (BASIL).

    This interface calculates:
    (1) arterial transit time,
    (2) CBF with spatial correction,
    (3) CBF with spatial partial volume white matter correction, and
    (4) CBF with spatial partial volume correction.

    See https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BASIL and https://asl-docs.readthedocs.io.
    """

    _cmd = "oxford_asl"
    input_spec = _BASILCBFInputSpec
    output_spec = _BASILCBFOutputSpec

    def _run_interface(self, runtime):
        runtime = super(BASILCBF, self)._run_interface(runtime)
        return runtime

    def _gen_outfilename(self, suffix):
        if isdefined(self.inputs.deltam):
            out_file = self._gen_fname(self.inputs.deltam, suffix=suffix)

        return os.path.abspath(out_file)

    def _list_outputs(self):
        basename = self.inputs.out_basename

        outputs = self.output_spec().get()

        outputs["mean_cbf_basil"] = os.path.join(basename, "native_space/perfusion_calib.nii.gz")
        outputs["att_basil"] = os.path.join(basename, "native_space/arrival.nii.gz")
        outputs["mean_cbf_gm_basil"] = os.path.join(
            basename,
            "native_space/pvcorr/perfusion_calib.nii.gz",
        )
        outputs["mean_cbf_wm_basil"] = os.path.join(
            basename,
            "native_space/pvcorr/perfusion_wm_calib.nii.gz",
        )

        return outputs


def regmotoasl(asl, m0file, m02asl):
    """Calculate mean M0 image and mean ASL image, then FLIRT M0 image to ASL space.

    TODO: This should not be a function. It uses interfaces, so it should be a workflow.
    """
    from nipype.interfaces import fsl

    meanasl = fsl.MeanImage()
    meanasl.inputs.in_file = asl
    meanasl.inputs.out_file = fname_presuffix(asl, suffix="_meanasl")
    meanasl.run()
    meanm0 = fsl.MeanImage()
    meanm0.inputs.in_file = m0file
    meanm0.inputs.out_file = fname_presuffix(asl, suffix="_meanm0")
    meanm0.run()
    flt = fsl.FLIRT(bins=640, cost_func="mutualinfo")
    flt.inputs.in_file = meanm0.inputs.out_file
    flt.inputs.reference = meanasl.inputs.out_file
    flt.inputs.out_file = m02asl
    flt.run()
    return m02asl


def refine_ref_mask(t1w_mask, ref_asl_mask, t12ref_transform, tmp_mask, refined_mask):
    """Warp T1w mask to ASL space, then use it to mask the ASL mask.

    TODO: This should not be a function. It uses interfaces, so it should be a workflow.
    """
    warp_t1w_mask_to_asl = ApplyTransforms(
        dimension=3,
        float=True,
        input_image=t1w_mask,
        interpolation="NearestNeighbor",
        reference_image=ref_asl_mask,
        transforms=[t12ref_transform],
        input_image_type=3,
        output_image=tmp_mask,
    )
    results = warp_t1w_mask_to_asl.run()

    modify_asl_mask = MultiImageMaths(
        in_file=results.outputs.output_image,
        op_string="-mul %s -bin",
        operand_files=ref_asl_mask,
        out_file=refined_mask,
    )
    results = modify_asl_mask.run()

    return results.outputs.out_file
