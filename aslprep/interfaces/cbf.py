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
from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec
from nipype.utils.filemanip import fname_presuffix

from aslprep import config
from aslprep.utils.asl import (
    check_m0tr,
    determine_multi_pld,
    pcasl_or_pasl,
)
from aslprep.utils.calibration import (
    CBF_UNIT_CONVERSION,
    calibration_value,
    m0_partition_coefficient,
    resolve_calibration,
)
from aslprep.utils.cbf import (
    _getcbfscore,
    _scrubcbf,
    att_status,
    fit_deltam_multipld,
    get_array_module,
    multipld_bolus,
    slice_shifted_delays,
    voxel_delays,
    zero_low_m0,
)


class _RefineMaskInputSpec(BaseInterfaceInputSpec):
    t1w_mask = File(exists=True, mandatory=True, desc='t1 mask')
    asl_mask = File(exists=True, mandatory=True, desc='asl mask')
    m0_mask = File(exists=True, mandatory=False, desc='M0 mask (if available)')


class _RefineMaskOutputSpec(TraitedSpec):
    out_mask = File(exists=False, desc='output mask')


class RefineMask(SimpleInterface):
    """Reduce the ASL-derived brain mask using the associated T1w mask and possibly an M0 mask."""

    input_spec = _RefineMaskInputSpec
    output_spec = _RefineMaskOutputSpec

    def _run_interface(self, runtime):
        from nilearn import image

        self._results['out_mask'] = fname_presuffix(
            self.inputs.asl_mask,
            suffix='_refinemask',
            newpath=runtime.cwd,
        )

        img1 = nb.load(self.inputs.asl_mask)
        img1 = nb.funcs.squeeze_image(img1)
        img2 = nb.load(self.inputs.t1w_mask)
        img2 = nb.funcs.squeeze_image(img2)
        if isdefined(self.inputs.m0_mask):
            img3 = nb.load(self.inputs.m0_mask)
            img3 = nb.funcs.squeeze_image(img3)
            out_mask = image.math_img(
                'img1 * img2 * img3',
                img1=img1,
                img2=img2,
                img3=img3,
            )
        else:
            out_mask = image.math_img(
                'img1 * img2',
                img1=img1,
                img2=img2,
            )

        out_mask.to_filename(self._results['out_mask'])

        return runtime


class _ExtractCBFInputSpec(BaseInterfaceInputSpec):
    name_source = File(exists=True, mandatory=True, desc='raw asl file')
    asl_file = File(exists=True, mandatory=True, desc='preprocessed asl file')
    metadata = traits.Dict(mandatory=True, desc='metadata for ASL file')
    aslcontext = File(exists=True, mandatory=True, desc='aslcontext TSV file for run.')
    m0scan = traits.Either(
        File(exists=True),
        None,
        mandatory=False,
        desc="m0scan file associated with the ASL file. Only defined if M0Type is 'Separate'.",
    )
    m0scan_metadata = traits.Either(
        traits.Dict,
        None,
        mandatory=True,
        desc="metadata for M0 scan. Only defined if M0Type is 'Separate'.",
    )
    fwhm = traits.Float(default_value=5, usedefault=True, mandatory=False, desc='fwhm')


class _ExtractCBFOutputSpec(TraitedSpec):
    out_file = File(exists=False, desc='Either CBF or deltaM time series.')
    m0_file = File(exists=False, desc='Mean M0 image, after smoothing.')
    metadata = traits.Dict(
        desc=(
            'Metadata for the ASL run. '
            'The dictionary may be modified to only include metadata associated with the selected '
            'volumes.'
        ),
    )
    m0tr = traits.Either(
        traits.Float,
        None,
        desc='RepetitionTimePreparation for M0 scans.',
    )


class ExtractCBF(SimpleInterface):
    """Extract CBF time series by subtracting label volumes from control volumes.

    TODO: Mock up test data and write tests to cover all of the branches in this interface.
    """

    input_spec = _ExtractCBFInputSpec
    output_spec = _ExtractCBFOutputSpec

    def _run_interface(self, runtime):
        aslcontext = pd.read_table(self.inputs.aslcontext)
        metadata = self.inputs.metadata.copy()

        # read the preprocessed ASL data
        asl_img = nb.load(self.inputs.asl_file)
        asl_data = asl_img.get_fdata()

        dims = asl_data.shape[:3]

        if aslcontext.shape[0] != asl_img.shape[3]:
            raise ValueError(
                f'Number of rows in aslcontext ({aslcontext.shape[0]}) != '
                f'number of volumes in ASL file ({asl_img.shape[3]})'
            )

        # get the control, label, m0scan, deltam, or cbf volumes
        vol_types = aslcontext['volume_type'].tolist()
        control_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == 'control']
        label_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == 'label']
        m0_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == 'm0scan']
        deltam_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == 'deltam']
        cbf_volume_idx = [i for i, vol_type in enumerate(vol_types) if vol_type == 'cbf']

        # extract m0 file and register it to ASL if separate
        if metadata['M0Type'] == 'Separate':
            m0file = self.inputs.m0scan
            m0data = smooth_image(nb.load(m0file), fwhm=self.inputs.fwhm).get_fdata()
            if len(m0data.shape) > 3:
                m0data = np.mean(m0data, axis=3)

            m0tr = self.inputs.m0scan_metadata['RepetitionTimePreparation']
            if np.array(m0tr).size > 1 and np.std(m0tr) > 0:
                raise ValueError('M0 scans have variable TR. ASLPrep does not support this.')

        elif metadata['M0Type'] == 'Included':
            m0data = asl_data[:, :, :, m0_volume_idx]
            m0img = nb.Nifti1Image(m0data, asl_img.affine, asl_img.header)
            m0data = smooth_image(m0img, fwhm=self.inputs.fwhm).get_fdata()
            m0data = np.mean(m0data, axis=3)

            if np.array(metadata['RepetitionTimePreparation']).size > 1:
                m0tr = np.array(metadata['RepetitionTimePreparation'])[m0_volume_idx]
            else:
                m0tr = metadata['RepetitionTimePreparation']

            if np.array(m0tr).size > 1 and np.std(m0tr) > 0:
                raise ValueError('M0 scans have variable TR. ASLPrep does not support this.')

        elif metadata['M0Type'] == 'Estimate':
            m0data = np.full(dims, metadata['M0Estimate'])

            m0tr = None

        elif metadata['M0Type'] == 'Absent':
            if control_volume_idx and not cbf_volume_idx:
                # BackgroundSuppression is required, so no need to use get().
                if metadata['BackgroundSuppression']:
                    raise ValueError(
                        'Background-suppressed control volumes cannot be used for calibration.'
                    )

            if control_volume_idx:
                # Estimate M0 using the smoothed mean control volumes.
                control_data = asl_data[:, :, :, control_volume_idx]
                control_img = nb.Nifti1Image(control_data, asl_img.affine, asl_img.header)
                control_img = smooth_image(control_img, fwhm=self.inputs.fwhm).get_fdata()
                m0data = np.mean(control_img, axis=3)

                # Use the control volumes' TR as the M0 TR.
                if np.array(metadata['RepetitionTimePreparation']).size > 1:
                    m0tr = np.array(metadata['RepetitionTimePreparation'])[control_volume_idx[0]]
                else:
                    m0tr = metadata['RepetitionTimePreparation']

            elif cbf_volume_idx:
                # If we have precalculated CBF data, we don't need M0, so we'll just use the mask.
                m0data = np.ones(dims)

                m0tr = None

            else:
                raise RuntimeError(
                    'm0scan is absent, '
                    'and there are no control volumes that can be used as a substitute'
                )

        else:
            raise RuntimeError('no pathway to m0scan')

        check_m0tr(m0tr)

        if deltam_volume_idx:
            config.loggers.interface.info('Extracting deltaM from ASL file.')
            metadata_idx = deltam_volume_idx
            out_data = asl_data[:, :, :, deltam_volume_idx]

        elif label_volume_idx:
            config.loggers.interface.info(
                'Calculating deltaM from label-control pairs in ASL file.'
            )
            if len(label_volume_idx) != len(control_volume_idx):
                raise ValueError(
                    f'Number of label volumes ({len(label_volume_idx)}) != '
                    f'number of control volumes ({len(control_volume_idx)})'
                )
            metadata_idx = control_volume_idx
            control_data = asl_data[:, :, :, control_volume_idx]
            label_data = asl_data[:, :, :, label_volume_idx]
            out_data = control_data - label_data

        elif cbf_volume_idx:
            metadata_idx = cbf_volume_idx
            out_data = asl_data[:, :, :, cbf_volume_idx]

        else:
            raise RuntimeError('No valid ASL or CBF image.')

        # Remove volume-wise metadata for M0 scans as necessary
        VOLUME_WISE_FIELDS = [
            'PostLabelingDelay',
            'VascularCrushingVENC',
            'LabelingDuration',
            'EchoTime',
            'FlipAngle',
            'RepetitionTimePreparation',
        ]

        for field in VOLUME_WISE_FIELDS:
            if field not in metadata:
                continue

            value = metadata[field]
            if isinstance(value, list) and len(value) != asl_data.shape[3]:
                raise ValueError(
                    f'{field} is an array, but the number of values ({len(value)}) '
                    f'does not match the number of volumes in the ASL data ({asl_data.shape[3]}).'
                )
            elif isinstance(value, list):
                # Reduce to only the selected volumes
                value = [value[i] for i in metadata_idx]

                metadata[field] = value

        self._results['metadata'] = metadata
        self._results['m0tr'] = np.mean(m0tr) if m0tr is not None else None
        self._results['out_file'] = fname_presuffix(
            self.inputs.name_source,
            suffix='_DeltaMOrCBF',
            newpath=runtime.cwd,
        )
        self._results['m0_file'] = fname_presuffix(
            self.inputs.name_source,
            suffix='_m0file',
            newpath=runtime.cwd,
        )
        nb.Nifti1Image(out_data, asl_img.affine, asl_img.header).to_filename(
            self._results['out_file']
        )
        nb.Nifti1Image(m0data, asl_img.affine, asl_img.header).to_filename(
            self._results['m0_file']
        )

        return runtime


class _ComputeCBFInputSpec(BaseInterfaceInputSpec):
    deltam = File(
        exists=True,
        mandatory=True,
        desc=(
            'NIfTI file containing raw CBF volume(s). '
            'These raw CBF values are the result of subtracting label volumes from '
            'control volumes, without any kind of additional scaling. '
            'This file may be 3D or 4D.'
        ),
    )
    metadata = traits.Dict(
        exists=True,
        mandatory=True,
        desc="Metadata for the raw CBF file, taken from the raw ASL data's sidecar JSON file.",
    )
    calibration = traits.Dict(
        desc=(
            'Calibration specification, from aslprep.utils.calibration.resolve_calibration. '
            'If not defined, it is resolved from metadata, m0_scale and m0tr.'
        ),
        xor=['m0_scale', 'm0tr'],
    )
    m0_scale = traits.Float(
        desc='Relative scale between ASL and M0. Only used without calibration (default 1).',
        xor=['calibration'],
    )
    m0_file = File(exists=True, mandatory=True, desc='M0 nifti file')
    m0tr = traits.Float(
        desc='M0 TR, in seconds. Only used without calibration.',
        xor=['calibration'],
    )
    cbf_only = traits.Bool(
        mandatory=True,
        desc='Whether data are deltam (False) or CBF (True).',
    )
    num_threads = traits.Int(
        1,
        usedefault=True,
        nohash=True,
        desc='Number of threads used to fit the multi-delay model. Ignored for single-delay data.',
    )
    use_gpu = traits.Bool(
        False,
        usedefault=True,
        desc='Fit the multi-delay model on a GPU, with CuPy. Ignored for single-delay data.',
    )
    write_calibration = traits.Bool(
        False,
        usedefault=True,
        desc=(
            'Also write the calibration M0 and the quantification-valid mask '
            '(m0_calib, quant_valid), which partial volume correction needs.'
        ),
    )


class _ComputeCBFOutputSpec(TraitedSpec):
    cbf_ts = traits.Either(
        File(exists=True),
        None,
        desc='Quantitative CBF time series, in mL/100g/min. Only generated for single-delay data.',
    )
    mean_cbf = File(exists=True, desc='Quantified CBF, averaged over time.')
    att = traits.Either(
        File(exists=True),
        None,
        desc='Arterial transit time map, in seconds. Only generated for multi-delay data.',
    )
    abat = traits.Either(
        File(exists=True),
        None,
        desc='Arterial bolus arrival time map, in seconds. Only generated for multi-delay data.',
    )
    abv = traits.Either(
        File(exists=True),
        None,
        desc='Arterial blood volume map. Only generated for multi-delay data.',
    )
    att_lower = traits.Either(
        File(exists=True),
        None,
        desc=(
            'Shortest arterial transit time that fits as well as the estimate, in seconds. '
            'Only generated for multi-delay data.'
        ),
    )
    att_upper = traits.Either(
        File(exists=True),
        None,
        desc=(
            'Longest arterial transit time that fits as well as the estimate, in seconds. '
            'Only generated for multi-delay data.'
        ),
    )
    abat_lower = traits.Either(
        File(exists=True),
        None,
        desc=(
            'Shortest arterial bolus arrival time that fits as well as the estimate, in seconds. '
            'Only generated for multi-delay data.'
        ),
    )
    abat_upper = traits.Either(
        File(exists=True),
        None,
        desc=(
            'Longest arterial bolus arrival time that fits as well as the estimate, in seconds. '
            'Only generated for multi-delay data.'
        ),
    )
    att_lower_ci = traits.Either(
        File(exists=True),
        None,
        desc=(
            'Lower end of the noise-aware arterial transit time interval, in seconds. '
            'Only generated for multi-delay data.'
        ),
    )
    att_upper_ci = traits.Either(
        File(exists=True),
        None,
        desc=(
            'Upper end of the noise-aware arterial transit time interval, in seconds. '
            'Only generated for multi-delay data.'
        ),
    )
    att_status = traits.Either(
        File(exists=True),
        None,
        desc=(
            'How well the fit located the arterial transit time, from the noise-aware interval '
            '(see aslprep.utils.cbf.ATT_STATUS). Only generated for multi-delay data.'
        ),
    )
    earliest_delay = traits.Either(
        File(exists=True),
        None,
        desc=(
            'Earliest delay of each voxel, including its slice-timing offset: the PLD for '
            '(P)CASL, or TI - TI1 for PASL. Only generated for multi-delay data.'
        ),
    )
    abv_decayed = traits.Either(
        File(exists=True),
        None,
        desc=(
            'Arterial blood volume multiplied by exp(-aBAT / T1blood), '
            'which does not depend on where aBAT falls within its interval between delays. '
            'Only generated for multi-delay (P)CASL data.'
        ),
    )
    plds = traits.Either(
        File(exists=True),
        None,
        desc=(
            'Delay times (post-labeling delays for (P)CASL, inversion times for PASL). '
            'Only defined if slice-timing correction is applied.'
        ),
    )
    failures = traits.Either(
        File(exists=True),
        None,
        desc='Boolean map of failed model fits. Only generated for multi-delay data.',
    )
    m0_calib = traits.Either(
        File(exists=True),
        None,
        desc=(
            'The M0 image that delta-M was divided by to quantify CBF: scaled, corrected for '
            'a short TR, and with low values set to zero. '
            'Only generated with write_calibration, and not for scanner-provided CBF.'
        ),
    )
    quant_valid = traits.Either(
        File(exists=True),
        None,
        desc=(
            'Mask of voxels where CBF could be quantified: M0 is finite and positive, and the '
            'delta-M values and quantification factor are finite. '
            'Only generated with write_calibration, and not for scanner-provided CBF.'
        ),
    )


def _image_finite(in_file, masker):
    """Whether every value of each voxel in an image is finite, in the masker's voxel order.

    NiftiMasker replaces non-finite values with zeros, so this is taken from the image itself.
    """
    img = nb.load(in_file)
    finite = np.isfinite(img.get_fdata())
    if finite.ndim == 4:
        finite = finite.all(axis=3)
    return masker.transform(nb.Nifti1Image(finite.astype(np.uint8), img.affine)).ravel() > 0


class ComputeCBF(SimpleInterface):
    """Calculate CBF time series and mean control.

    Notes
    -----
    This interface calculates CBF from deltam and M0 data.
    It can handle single-delay and multi-delay data, single CBF volumes and CBF time series,
    and PASL and (P)CASL data.

    Single-delay CBF, for both (P)CASL and QUIPSSII PASL
    is calculated according to :footcite:t:`alsop_recommended_2015`.
    Multi-delay CBF is calculated according to :footcite:t:`woods2023recommendations`.

    If slice timing information is detected, then PLDs will be shifted by the slice times.

    See Also
    --------
    :func:`~aslprep.utils.asl.pcasl_or_pasl`
    :func:`~aslprep.utils.asl.determine_multi_pld`
    :func:`~aslprep.utils.calibration.resolve_calibration`
    :func:`~aslprep.utils.cbf.calculate_deltam_pasl`
    :func:`~aslprep.utils.cbf.calculate_deltam_pcasl`
    :func:`~aslprep.utils.cbf.fit_deltam_multipld`

    References
    ----------
    .. footbibliography::
    """

    input_spec = _ComputeCBFInputSpec
    output_spec = _ComputeCBFOutputSpec

    def _run_interface(self, runtime):
        metadata = self.inputs.metadata
        m0_file = self.inputs.m0_file
        deltam_file = self.inputs.deltam  # control - label signal intensities

        deltam_img = nb.load(deltam_file)
        dims = deltam_img.shape[:3]
        mask_img = nb.Nifti1Image(np.ones(dims), deltam_img.affine, deltam_img.header)
        self._results['m0_calib'] = None
        self._results['quant_valid'] = None

        if self.inputs.cbf_only:
            config.loggers.interface.debug('CBF data detected. Skipping CBF estimation.')
            self._results['cbf_ts'] = fname_presuffix(
                deltam_file,
                suffix='_cbf_ts',
                newpath=runtime.cwd,
            )
            cbf_img = nb.load(deltam_file)
            cbf_img.to_filename(self._results['cbf_ts'])
            self._results['mean_cbf'] = fname_presuffix(
                deltam_file,
                suffix='_meancbf',
                newpath=runtime.cwd,
            )
            mean_cbf_img = image.mean_img(cbf_img)
            mean_cbf_img.to_filename(self._results['mean_cbf'])

            # No ATT available for pre-calculated CBF
            self._results['att'] = None

            return runtime

        is_casl = pcasl_or_pasl(metadata=metadata)
        is_multi_pld = determine_multi_pld(metadata=metadata)
        if isdefined(self.inputs.calibration):
            calibration = self.inputs.calibration
        else:
            calibration = resolve_calibration(
                metadata,
                m0tr=self.inputs.m0tr if isdefined(self.inputs.m0tr) else None,
                m0_scale=self.inputs.m0_scale if isdefined(self.inputs.m0_scale) else 1.0,
            )
            for warning in calibration['Warnings']:
                config.loggers.interface.warning(warning)

        t1blood = calibration_value(calibration, 'BloodT1')
        # 1 for an M0Estimate, which is already the M0 of blood
        partition_coefficient = m0_partition_coefficient(calibration)
        # Labeling efficiency, including the loss from background suppression (alpha in Alsop 2015)
        labeleff = calibration_value(calibration, 'EffectiveLabelingEfficiency')
        # mL/g/s to mL/100 g/min. A fixed constant, which the multi-delay model and BASIL share,
        # so it is recorded in the specification but not read from it.
        UNIT_CONV = CBF_UNIT_CONVERSION

        # PostLabelingDelay is either a single number or an array of numbers.
        # If it is an array of numbers, then there should be one value for every volume in the
        # time series, with any M0 volumes having a value of 0.
        plds = np.atleast_1d(metadata['PostLabelingDelay'])

        # NOTE: Nilearn will still add a singleton time dimension for 3D imgs with
        # NiftiMasker.transform, until 0.12.0, so the arrays will currently be 2D no matter what.
        masker = maskers.NiftiMasker(mask_img=mask_img)
        deltam_arr = masker.fit_transform(deltam_file).T  # Transpose to SxT
        if deltam_arr.ndim != 2:
            raise ValueError(f'deltam is {deltam_arr.ndim}')

        # Load the M0 map and average over time, in case there's more than one map in the file.
        m0data = masker.transform(m0_file)
        # TODO: Scale each M0 volume by its TR separately instead of averaging.
        if m0data.ndim != 1:
            m0data = np.mean(m0data, axis=0)

        recovery_factor = calibration_value(calibration, 'M0RecoveryFactor')
        if recovery_factor is not None:
            config.loggers.interface.warning(
                'M0 TR is less than 5 seconds '
                f'({calibration_value(calibration, "M0RepetitionTime")}), '
                'so ASLPrep will scale the M0 data according to Alsop 2015.'
            )
            # Alsop 2015, page 113
            m0data = m0data * recovery_factor

        # Voxels with almost no M0 signal get no CBF estimate (zero for single-delay data,
        # a failed fit for multi-delay data), rather than an arbitrarily large one.
        scaled_m0data = calibration_value(calibration, 'M0Scale') * m0data
        m0_threshold = calibration_value(calibration, 'M0Threshold')
        if m0_threshold is not None:
            scaled_m0data = zero_low_m0(scaled_m0data, threshold=m0_threshold)

        self._results['plds'] = None
        plds, pld_img = voxel_delays(
            metadata, nb.load(deltam_file), masker, deltam_arr.shape[0], is_multi_pld
        )
        if pld_img is not None:
            # Write out the slice-shifted PLDs to the working directory, for debugging.
            pld_file = fname_presuffix(
                deltam_file,
                suffix='_plds',
                newpath=runtime.cwd,
            )
            pld_img.to_filename(pld_file)
            self._results['plds'] = pld_file

        # tau is defined for (P)CASL
        tau = metadata['LabelingDuration'] if is_casl else None

        # Now estimate CBF and any other metrics
        if is_multi_pld:
            tau, ti1 = multipld_bolus(metadata, is_casl)

            cbf, att, abat, abv, failures, bounds = fit_deltam_multipld(
                deltam_arr=deltam_arr,
                scaled_m0data=scaled_m0data,
                plds=plds,
                labeleff=labeleff,
                t1blood=t1blood,
                partition_coefficient=partition_coefficient,
                is_casl=is_casl,
                tau=tau,  # defined for (P)CASL
                ti1=ti1,  # defined for PASL
                xp=get_array_module(self.inputs.use_gpu),
                n_jobs=self.inputs.num_threads,
                att_interval=True,
            )

            mean_cbf_img = masker.inverse_transform(cbf)
            att_img = masker.inverse_transform(att)
            abat_img = masker.inverse_transform(abat)
            abv_img = masker.inverse_transform(abv)
            failures_img = masker.inverse_transform(failures.astype(int))

            # Multi-delay data won't produce a CBF time series
            self._results['cbf_ts'] = None
            self._results['att'] = fname_presuffix(
                self.inputs.deltam,
                suffix='_att',
                newpath=runtime.cwd,
            )
            att_img.to_filename(self._results['att'])
            self._results['abat'] = fname_presuffix(
                self.inputs.deltam,
                suffix='_abat',
                newpath=runtime.cwd,
            )
            abat_img.to_filename(self._results['abat'])
            self._results['abv'] = fname_presuffix(
                self.inputs.deltam,
                suffix='_abv',
                newpath=runtime.cwd,
            )
            abv_img.to_filename(self._results['abv'])
            for name, bound in bounds.items():
                self._results[name] = fname_presuffix(
                    self.inputs.deltam,
                    suffix=f'_{name}',
                    newpath=runtime.cwd,
                )
                masker.inverse_transform(bound).to_filename(self._results[name])

            # The earliest ATT at which the tissue signal no longer depends on ATT
            earliest = np.min(plds, axis=1) if is_casl else np.min(plds, axis=1) - ti1
            self._results['earliest_delay'] = fname_presuffix(
                self.inputs.deltam,
                suffix='_earliest_delay',
                newpath=runtime.cwd,
            )
            masker.inverse_transform(np.maximum(earliest, 0)).to_filename(
                self._results['earliest_delay']
            )

            status = att_status(cbf, bounds['att_lower_ci'], bounds['att_upper_ci'], failures)
            # The masker returns floats, but the status is categorical
            status_img = masker.inverse_transform(status)
            status_img = nb.Nifti1Image(
                np.asarray(status_img.dataobj).astype(np.uint8), status_img.affine
            )
            self._results['att_status'] = fname_presuffix(
                self.inputs.deltam,
                suffix='_attstatus',
                newpath=runtime.cwd,
            )
            status_img.to_filename(self._results['att_status'])

            # For (P)CASL, the data only determine the product of aBV and the decay of the
            # label over aBAT. For PASL, the arterial signal does not depend on aBAT.
            self._results['abv_decayed'] = None
            if is_casl:
                self._results['abv_decayed'] = fname_presuffix(
                    self.inputs.deltam,
                    suffix='_abv_decayed',
                    newpath=runtime.cwd,
                )
                abv_decayed_img = masker.inverse_transform(abv * np.exp(-abat / t1blood))
                abv_decayed_img.to_filename(self._results['abv_decayed'])
            self._results['failures'] = fname_presuffix(
                self.inputs.deltam,
                suffix='_failures',
                newpath=runtime.cwd,
            )
            failures_img.to_filename(self._results['failures'])

        else:  # Single-delay
            if is_casl:
                denom_factor = t1blood * (1 - np.exp(-(tau / t1blood)))

            elif not metadata['BolusCutOffFlag']:
                raise ValueError(
                    'PASL without a bolus cut-off technique is not supported in ASLPrep.'
                )

            elif metadata['BolusCutOffTechnique'] == 'QUIPSS':
                # PASL + QUIPSS
                # Only one BolusCutOffDelayTime allowed.
                if not isinstance(metadata['BolusCutOffDelayTime'], Number):
                    raise ValueError(
                        f'Expected a single BolusCutOffDelayTime, but got '
                        f'{metadata["BolusCutOffDelayTime"]}'
                    )
                denom_factor = plds - metadata['BolusCutOffDelayTime']  # delta_TI, per Wong 1998

            elif metadata['BolusCutOffTechnique'] == 'QUIPSSII':
                # PASL + QUIPSSII
                # Per SD, use PLD as TI for PASL, so we will just use 'plds' in the numerator when
                # calculating the perfusion factor.
                # Only one BolusCutOffDelayTime allowed.
                if not isinstance(metadata['BolusCutOffDelayTime'], Number):
                    raise ValueError(
                        f'Expected a single BolusCutOffDelayTime, but got '
                        f'{metadata["BolusCutOffDelayTime"]}'
                    )
                denom_factor = metadata['BolusCutOffDelayTime']  # called TI1 in Alsop 2015

            elif metadata['BolusCutOffTechnique'] == 'Q2TIPS':
                # PASL + Q2TIPS
                # Q2TIPS should have two BolusCutOffDelayTimes.
                if len(metadata['BolusCutOffDelayTime']) != 2:
                    raise ValueError(
                        f'Expected two BolusCutOffDelayTimes, but got '
                        f'{metadata["BolusCutOffDelayTime"]}'
                    )
                # The first BolusCutOffDelayTime is the bolus duration (TI1 in Alsop 2015).
                # The second is the end of the saturation pulse train, not the inversion time
                # (TI in Alsop 2015, TI2 in Noguchi 2015), which is the PLD, as for QUIPSS II.
                denom_factor = metadata['BolusCutOffDelayTime'][0]

            else:
                raise ValueError(
                    f'Unknown BolusCutOffTechnique {metadata["BolusCutOffTechnique"]}'
                )

            # Scale difference signal to absolute CBF units by dividing by PD image (M0 * M0scale).
            with np.errstate(divide='ignore', invalid='ignore'):
                deltam_scaled = deltam_arr / scaled_m0data[:, None]

            perfusion_factor = (UNIT_CONV * partition_coefficient * np.exp(plds / t1blood)) / (
                denom_factor * 2 * labeleff
            )

            cbf_ts = deltam_scaled * perfusion_factor
            cbf_ts = np.nan_to_num(cbf_ts, nan=0, posinf=0, neginf=0)
            cbf_ts_img = masker.inverse_transform(cbf_ts.T)
            mean_cbf_img = image.mean_img(cbf_ts_img)
            self._results['cbf_ts'] = fname_presuffix(
                self.inputs.deltam,
                suffix='_cbf',
                newpath=runtime.cwd,
            )
            cbf_ts_img.to_filename(self._results['cbf_ts'])
            # Single-delay data won't produce an ATT image
            self._results['att'] = None

        if self.inputs.write_calibration:
            # The calibration image and where quantification can succeed, for partial volume
            # correction, which unmixes the M0-free signal and calibrates with this M0 again.
            # The masker replaces non-finite values with zeros, so finiteness is taken from the
            # images themselves, mapped through the same masker to keep the voxel order.
            quant_valid = (
                np.isfinite(scaled_m0data)
                & (scaled_m0data > 0)
                & _image_finite(deltam_file, masker)
                & _image_finite(m0_file, masker)
                & np.all(np.isfinite(np.broadcast_to(plds, deltam_arr.shape)), axis=1)
            )
            if not is_multi_pld:
                quant_valid &= np.all(
                    np.isfinite(np.broadcast_to(perfusion_factor, deltam_arr.shape)), axis=1
                )
            self._results['m0_calib'] = fname_presuffix(
                self.inputs.deltam,
                suffix='_m0calib',
                newpath=runtime.cwd,
            )
            masker.inverse_transform(np.nan_to_num(scaled_m0data)).to_filename(
                self._results['m0_calib']
            )
            self._results['quant_valid'] = fname_presuffix(
                self.inputs.deltam,
                suffix='_quantvalid',
                newpath=runtime.cwd,
            )
            masker.inverse_transform(quant_valid.astype(np.uint8)).to_filename(
                self._results['quant_valid']
            )

        # Mean CBF is returned no matter what
        self._results['mean_cbf'] = fname_presuffix(
            self.inputs.deltam,
            suffix='_meancbf',
            newpath=runtime.cwd,
        )
        mean_cbf_img.to_filename(self._results['mean_cbf'])

        return runtime


class _WriteCalibrationInputSpec(BaseInterfaceInputSpec):
    calibration = traits.Dict(
        mandatory=True,
        desc='Calibration specification, from aslprep.utils.calibration.resolve_calibration.',
    )


class _WriteCalibrationOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='The calibration specification, as a JSON file.')


class WriteCalibration(SimpleInterface):
    """Write a calibration specification to a JSON file."""

    input_spec = _WriteCalibrationInputSpec
    output_spec = _WriteCalibrationOutputSpec

    def _run_interface(self, runtime):
        import json

        out_file = os.path.join(runtime.cwd, 'calibration_parameters.json')
        with open(out_file, 'w') as fobj:
            json.dump(self.inputs.calibration, fobj, indent=2, sort_keys=False)
        self._results['out_file'] = out_file
        return runtime


class _ScoreAndScrubCBFInputSpec(BaseInterfaceInputSpec):
    cbf_ts = File(exists=True, mandatory=True, desc='Computed CBF from ComputeCBF.')
    mask = File(exists=True, mandatory=True, desc='mask')
    gm_tpm = File(exists=True, mandatory=True, desc='Gray matter tissue probability map.')
    wm_tpm = File(exists=True, mandatory=True, desc='White matter tissue probability map.')
    csf_tpm = File(exists=True, mandatory=True, desc='CSF tissue probability map.')
    tpm_threshold = traits.Float(
        default_value=0.7,
        usedefault=True,
        mandatory=False,
        desc='Tissue probability threshold for binarizing GM, WM, and CSF masks.',
    )
    wavelet_function = traits.Enum(
        'huber',
        'andrews',
        'bisquare',
        'cauchy',
        'fair',
        'logistic',
        'ols',
        'talwar',
        'welsch',
        usedefault=True,
        mandatory=False,
        desc="Robust weight function for SCRUB's data fidelity term.",
    )


class _ScoreAndScrubCBFOutputSpec(TraitedSpec):
    cbf_ts_score = File(exists=False, mandatory=False, desc='score timeseries data')
    mean_cbf_score = File(exists=False, mandatory=False, desc='average score')
    mean_cbf_scrub = File(exists=False, mandatory=False, desc='average scrub')
    score_outlier_index = File(exists=False, mandatory=False, desc='index of volume remove ')


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
        cbf_ts = nb.load(self.inputs.cbf_ts).get_fdata()
        mask = nb.load(self.inputs.mask).get_fdata()
        greym = nb.load(self.inputs.gm_tpm).get_fdata()
        whitem = nb.load(self.inputs.wm_tpm).get_fdata()
        csf = nb.load(self.inputs.csf_tpm).get_fdata()
        if cbf_ts.ndim > 3:
            cbf_scorets, index_score = _getcbfscore(
                cbfts=cbf_ts,
                wm=whitem,
                gm=greym,
                csf=csf,
                mask=mask,
                thresh=self.inputs.tpm_threshold,
            )
            cbfscrub = _scrubcbf(
                cbf_ts=cbf_scorets,
                gm=greym,
                wm=whitem,
                csf=csf,
                mask=mask,
                wfun=self.inputs.wavelet_function,
                thresh=self.inputs.tpm_threshold,
            )
            mean_cbf_score = np.mean(cbf_scorets, axis=3)
        else:
            config.loggers.interface.warning(
                f'CBF time series is only {cbf_ts.ndim}D. Skipping SCORE and SCRUB.'
            )
            cbf_scorets = cbf_ts
            index_score = np.array([0])
            cbfscrub = cbf_ts
            mean_cbf_score = cbf_ts

        self._results['cbf_ts_score'] = fname_presuffix(
            self.inputs.cbf_ts,
            suffix='_cbfscorets',
            newpath=runtime.cwd,
        )
        self._results['mean_cbf_score'] = fname_presuffix(
            self.inputs.cbf_ts,
            suffix='_meancbfscore',
            newpath=runtime.cwd,
        )
        self._results['mean_cbf_scrub'] = fname_presuffix(
            self.inputs.cbf_ts,
            suffix='_cbfscrub',
            newpath=runtime.cwd,
        )
        self._results['score_outlier_index'] = fname_presuffix(
            self.inputs.cbf_ts,
            suffix='_scoreindex.tsv',
            newpath=runtime.cwd,
            use_ext=False,
        )
        # nibabel writes an image in its header's data type, and the mask is uint8,
        # so the CBF maps get a float32 copy of the mask's header.
        mask_img = nb.load(self.inputs.mask)
        header = mask_img.header.copy()
        header.set_data_dtype(np.float32)
        for out_name, data in (
            ('cbf_ts_score', cbf_scorets),
            ('mean_cbf_score', mean_cbf_score),
            ('mean_cbf_scrub', cbfscrub),
        ):
            nb.Nifti1Image(data.astype(np.float32), mask_img.affine, header).to_filename(
                self._results[out_name]
            )

        score_outlier_df = pd.DataFrame(columns=['score_outlier_index'], data=index_score)
        score_outlier_df.to_csv(self._results['score_outlier_index'], sep='\t', index=False)

        return runtime


class _BASILInflowTimeImageInputSpec(BaseInterfaceInputSpec):
    deltam = File(
        exists=True,
        mandatory=True,
        desc='The deltaM file passed to BASIL, on the acquired slice grid.',
    )
    tis = traits.Either(
        traits.Float(),
        traits.List(traits.Float()),
        mandatory=True,
        desc='The inflow times BASIL would otherwise get with ``--tis``.',
    )
    metadata = traits.Dict(
        mandatory=True,
        desc='ASL metadata, with SliceTiming.',
    )


class _BASILInflowTimeImageOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='4D image of the inflow time of each voxel.')


class BASILInflowTimeImage(SimpleInterface):
    """Write the voxelwise inflow time image for BASIL's ``--tiimg``.

    Each volume holds one of the inflow times, plus the acquisition time of each voxel's slice.
    BASIL reads these values as inflow times (TIs), and does not add the labeling duration,
    so for (P)CASL data they must already include it.
    """

    input_spec = _BASILInflowTimeImageInputSpec
    output_spec = _BASILInflowTimeImageOutputSpec

    def _run_interface(self, runtime):
        tiimg = slice_shifted_delays(
            self.inputs.tis,
            self.inputs.metadata,
            nb.load(self.inputs.deltam),
        )
        tiimg = nb.Nifti1Image(tiimg.get_fdata(dtype=np.float32), tiimg.affine, tiimg.header)
        tiimg.set_data_dtype(np.float32)
        self._results['out_file'] = fname_presuffix(
            self.inputs.deltam,
            suffix='_tiimg',
            newpath=runtime.cwd,
        )
        tiimg.to_filename(self._results['out_file'])
        return runtime


class _BASILCBFInputSpec(FSLCommandInputSpec):
    # We use position args here as list indices - so a negative number
    # will put something on the end
    deltam = File(
        exists=True,
        desc=(
            'ASL difference (deltaM) data, after subtracting the label images from the control '
            'images. '
            'This matches with ``--iaf diff``, which is the default.'
        ),
        argstr='-i %s',
        position=0,
        mandatory=True,
    )
    mask = File(
        exists=True,
        argstr='-m %s',
        desc='mask in the same space as deltam',
        mandatory=True,
    )
    mzero = File(
        exists=True,
        argstr='-c %s',
        desc='Calibration (M0) image, which oxford_asl divides by the partition coefficient.',
        xor=['m0_value'],
    )
    m0_value = traits.Float(
        desc=(
            'A single M0 value of arterial blood (an M0Estimate, multiplied by M0Scale), '
            'used in place of a calibration image. oxford_asl does not divide it by the '
            'partition coefficient, or apply --cgain or --tr to it.'
        ),
        argstr='--M0 %s',
        xor=['mzero', 'm0_scale'],
    )
    # oxford_asl only passes --cgain to fslmaths, which reads exponents, so small scales keep
    # their precision. The other values also go through bc, which does not read exponents.
    m0_scale = traits.Float(
        desc='Relative scale between ASL and M0 (the M0Scale calibration parameter).',
        argstr='--cgain %.10g',
        xor=['m0_value'],
    )
    m0tr = traits.Float(
        desc=(
            'The repetition time for the calibration image (the M0 scan). '
            'oxford_asl corrects the calibration image for T1 recovery if it is below 5 s, '
            'and assumes 3.2 s if it is not given.'
        ),
        argstr='--tr %.6f',
        mandatory=False,
    )
    t1tissue_calib = traits.Float(
        desc=(
            'T1 of tissue, in seconds, used to correct the calibration image for T1 recovery '
            '(the TissueT1 calibration parameter) [default: 1.3].'
        ),
        argstr='--t1t %.6f',
    )
    t1tissue = traits.Float(
        desc="T1 of tissue, in seconds, in BASIL's kinetic model [default: 1.3].",
        argstr='--t1 %.6f',
    )
    tis = traits.Either(
        traits.Float(),
        traits.List(traits.Float()),
        desc=(
            'The list of inflow times (TIs), a comma separated list of values should be provided '
            '(that matches the order in the data).\n\n'
            'Note, the inflow time is the PLD plus bolus duration for pcASL (and cASL), '
            'it equals the inversion time for pASL. '
            'If the data contains multiple repeats of the same set of TIs then it is only '
            'necessary to list the unique TIs.\n\n'
            'When using the ``--tis=`` you can specify a full list of all TIs/PLDs in the data '
            '(i.e., as many entries as there are label-control pairs). '
            'Or, if you have a number of TIs/PLDs repeated multiple times you can just list the '
            'unique TIs in order and ``oxford_asl`` will automatically replicate that list to '
            'match the number of repeated measurements in the data. '
            'If you have a variable number of repeats at each TI/PLD then either list all TIs '
            'or use the ``--rpts=<csv>`` option (see below).'
        ),
        argstr='--tis %s',
        mandatory=True,
        xor=['tis', 'tiimg'],
        sep=',',
    )
    tiimg = File(
        exists=True,
        desc=(
            'A 4D image of the inflow time (TI) of each voxel, with one volume per entry of '
            '``tis``. Used in place of ``tis`` when slice timing cannot be described by '
            '``--slicedt``, such as for interleaved slice orders.'
        ),
        argstr='--tiimg %s',
        mandatory=True,
        xor=['tis', 'tiimg', 'slice_spacing', 'sliceband'],
    )
    pcasl = traits.Bool(
        desc=(
            'Data were acquired using CASL or PCASL labeling '
            '(PASL labeling is assumed by default).'
        ),
        argstr='--casl',
        mandatory=False,
        default_value=False,
    )
    bolus = traits.Either(
        traits.Float(),
        traits.List(traits.Float()),
        desc='bolus or tau: label duration',
        argstr='--bolus %s',
        mandatory=True,
        sep=',',
    )
    fixbolus = traits.Bool(
        desc=(
            'The bolus duration is fixed, e.g., by QUIPSS II, Q2TIPS or (P)CASL labeling. '
            'Without this flag, oxford_asl estimates the bolus duration for multi-delay PASL data, '
            'with ``--bolus`` only setting its prior.'
        ),
        argstr='--fixbolus',
        mandatory=False,
    )
    slice_spacing = traits.Float(
        desc='Increase in inflow time per slice, in seconds.',
        argstr='--slicedt %s',
        mandatory=False,
        xor=['tiimg'],
    )
    sliceband = traits.Int(
        desc='Number of slices per band in a multi-band acquisition.',
        argstr='--sliceband %d',
        mandatory=False,
        xor=['tiimg'],
    )
    pvc = traits.Bool(
        desc='Do partial volume correction.',
        mandatory=False,
        argstr='--pvcorr',
        default_value=True,
    )
    gm_tpm = File(
        exists=True,
        mandatory=False,
        desc='Partial volume estimates for GM. This is just a GM tissue probability map.',
        argstr='--pvgm %s',
    )
    wm_tpm = File(
        exists=True,
        mandatory=False,
        desc='Partial volume estimates for WM. This is just a WM tissue probability map.',
        argstr='--pvwm %s',
    )
    alpha = traits.Float(
        desc=(
            'Inversion efficiency - [default: 0.98 (pASL); 0.85 (cASL)]. '
            'ASLPrep passes the EffectiveLabelingEfficiency calibration parameter: '
            "the BIDS metadata field 'LabelingEfficiency', reduced by background suppression."
        ),
        argstr='--alpha %.6f',
    )
    t1blood = traits.Float(
        desc=(
            'T1 of blood, in seconds [default: 1.65]. '
            'ASLPrep passes the BloodT1 calibration parameter, which its own models use.'
        ),
        argstr='--t1b %.6f',
    )


class _BASILCBFOutputSpec(TraitedSpec):
    mean_cbf_basil = File(exists=True, desc='cbf with spatial correction')
    mean_cbf_gm_basil = File(exists=True, desc='cbf with spatial correction')
    mean_cbf_wm_basil = File(
        exists=True,
        desc='cbf with spatial partial volume white matter correction',
    )
    att_basil = File(exists=True, desc='arterial transit time')
    abv_basil = File(exists=True, desc='arterial blood volume')
    cbf_sd_basil = File(exists=True, desc='posterior standard deviation of mean_cbf_basil')
    cbf_gm_sd_basil = File(
        exists=True,
        desc='posterior standard deviation of mean_cbf_gm_basil',
    )
    cbf_wm_sd_basil = File(
        exists=True,
        desc='posterior standard deviation of mean_cbf_wm_basil',
    )
    att_sd_basil = File(exists=True, desc='posterior standard deviation of att_basil')


# Posterior variance maps written by oxford_asl, and the standard-deviation outputs made from
# them. The CBF variances are calibrated by oxford_asl, so they are in (mL/100 g/min)^2.
_BASIL_VARIANCE_MAPS = {
    'cbf_sd_basil': 'native_space/perfusion_var_calib.nii.gz',
    'cbf_gm_sd_basil': 'native_space/pvcorr/perfusion_var_calib.nii.gz',
    'cbf_wm_sd_basil': 'native_space/pvcorr/perfusion_wm_var_calib.nii.gz',
    'att_sd_basil': 'native_space/arrival_var.nii.gz',
}


def _variance_to_sd(variance_file, mask_file, out_file):
    """Write the square root of a BASIL variance map, with zeros outside the analysis mask.

    oxford_asl fills voxels it did not fit with a variance of 1e12,
    which would otherwise show up as a standard deviation of 1e6.
    """
    img = nb.load(variance_file)
    variance = np.asanyarray(img.dataobj, dtype=np.float64)
    mask = np.asanyarray(nb.load(mask_file).dataobj) > 0
    sd = np.zeros_like(variance, dtype=np.float32)
    sd[mask] = np.sqrt(np.clip(variance[mask], 0, None))
    out_img = nb.Nifti1Image(sd, img.affine, img.header)
    out_img.set_data_dtype(np.float32)
    out_img.to_filename(out_file)


class BASILCBF(FSLCommand):
    """Apply Bayesian Inference for Arterial Spin Labeling (BASIL).

    This interface calculates:
    (1) arterial transit time,
    (2) CBF with spatial correction,
    (3) CBF with spatial partial volume white matter correction, and
    (4) CBF with spatial partial volume correction.

    See https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BASIL and https://asl-docs.readthedocs.io.
    """

    _cmd = 'oxford_asl'
    input_spec = _BASILCBFInputSpec
    output_spec = _BASILCBFOutputSpec

    def _format_arg(self, name, spec, value):
        if name == 'm0_value':
            # oxford_asl multiplies the value with bc, which does not read exponents
            value = np.format_float_positional(value, precision=10, fractional=False, trim='-')
            return spec.argstr % value
        return super()._format_arg(name, spec, value)

    def _run_interface(self, runtime):
        # oxford_asl keeps its intermediate files in a temporary directory named by FSL's
        # tmpnam, which creates it in TMPDIR, or in /tmp if TMPDIR is not set.
        # Inside a container, /tmp may be too small for these files,
        # so point TMPDIR to the node's working directory.
        runtime.environ['TMPDIR'] = runtime.cwd
        runtime = super()._run_interface(runtime)

        # Convert BASIL's posterior variance maps to standard deviations,
        # which share the units of the estimates.
        mask_file = os.path.join(runtime.cwd, 'native_space', 'mask.nii.gz')
        for output, variance_file in _BASIL_VARIANCE_MAPS.items():
            variance_file = os.path.join(runtime.cwd, variance_file)
            if os.path.isfile(variance_file) and os.path.isfile(mask_file):
                _variance_to_sd(
                    variance_file,
                    mask_file,
                    os.path.join(runtime.cwd, f'{output}.nii.gz'),
                )

        return runtime

    def _gen_outfilename(self, suffix):
        if isdefined(self.inputs.deltam):
            out_file = self._gen_fname(self.inputs.deltam, suffix=suffix)
        return os.path.abspath(out_file)

    def _list_outputs(self):
        outputs = self.output_spec().get()

        outputs['mean_cbf_basil'] = os.path.abspath('native_space/perfusion_calib.nii.gz')
        # oxford_asl only writes the partial-volume corrected maps with --pvcorr
        if self.inputs.pvc is True:
            outputs['mean_cbf_gm_basil'] = os.path.abspath(
                'native_space/pvcorr/perfusion_calib.nii.gz',
            )
            outputs['mean_cbf_wm_basil'] = os.path.abspath(
                'native_space/pvcorr/perfusion_wm_calib.nii.gz',
            )

        # These are only meaningful for multi-PLD data, and oxford_asl only writes them when
        # it infers arterial transit time / the macrovascular (arterial) component. For
        # single-PLD data these files are not produced, so only report them when they exist
        # to avoid failing the mandatory existence check during output aggregation. Whether
        # they are surfaced as derivatives is gated on multi-PLD status in the CBF workflow
        # (see init_cbf_wf).
        att_basil = os.path.abspath('native_space/arrival.nii.gz')
        if os.path.isfile(att_basil):
            outputs['att_basil'] = att_basil

        abv_basil = os.path.abspath('native_space/aCBV_calib.nii.gz')
        if os.path.isfile(abv_basil):
            outputs['abv_basil'] = abv_basil

        # Standard deviations exist only for the variance maps oxford_asl wrote.
        # oxford_asl writes no variance for aCBV.
        for output in _BASIL_VARIANCE_MAPS:
            if output in ('cbf_gm_sd_basil', 'cbf_wm_sd_basil') and self.inputs.pvc is not True:
                continue
            sd_file = os.path.abspath(f'{output}.nii.gz')
            if os.path.isfile(sd_file):
                outputs[output] = sd_file

        return outputs
