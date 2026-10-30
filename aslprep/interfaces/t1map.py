"""Interfaces that use a quantitative T1 map in the tissue-T1 kinetic model."""

import os

import nibabel as nb
import numpy as np
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    isdefined,
    traits,
)

from aslprep import config


class _PrepareT1MapInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='The T1 map.')
    metadata = traits.Dict(desc="The T1 map's sidecar, whose Units decide its scale if present.")


class _PrepareT1MapOutputSpec(TraitedSpec):
    t1 = File(exists=True, desc='T1, in seconds. Non-finite and non-positive values are zero.')
    r1 = File(
        exists=True,
        desc='R1 (1 / T1), in 1/s, zero where T1 is zero, for registering the map itself.',
    )
    units = traits.Str(desc='How the units of the map were decided.')


class PrepareT1Map(SimpleInterface):
    """Convert a T1 map to seconds, and compute its R1 map.

    The R1 map has a T1-weighted contrast, so the T1 map can be registered to the anatomical image
    when there is no T1-weighted image on its grid.
    """

    input_spec = _PrepareT1MapInputSpec
    output_spec = _PrepareT1MapOutputSpec

    def _run_interface(self, runtime):
        from aslprep.utils.t1map import t1_scale

        img = nb.load(self.inputs.in_file)
        values = np.asanyarray(img.dataobj, dtype=np.float64)
        if values.ndim == 4 and values.shape[3] == 1:
            values = values[..., 0]

        if values.ndim != 3:
            raise ValueError(f'The T1 map must be a 3D image, but it has shape {values.shape}.')

        metadata = self.inputs.metadata if isdefined(self.inputs.metadata) else {}
        scale, how = t1_scale(values, metadata)
        t1 = np.where(np.isfinite(values) & (values > 0), values * scale, 0.0)
        with np.errstate(divide='ignore'):
            r1 = np.where(t1 > 0, 1 / t1, 0.0)

        header = img.header.copy()
        header.set_data_dtype(np.float32)
        for name, data in (('t1', t1), ('r1', r1)):
            out_file = os.path.join(runtime.cwd, f'{name}.nii.gz')
            nb.Nifti1Image(data.astype(np.float32), img.affine, header).to_filename(out_file)
            self._results[name] = out_file

        self._results['units'] = how
        return runtime


class _MeasureTissueT1InputSpec(BaseInterfaceInputSpec):
    t1 = File(exists=True, mandatory=True, desc='T1, in seconds, on the grid of the tissue maps.')
    tpms = traits.List(
        File(exists=True),
        minlen=3,
        mandatory=True,
        desc='Gray matter, white matter and CSF probability maps, in that order.',
    )
    t1map_name = traits.Str(mandatory=True, desc='The T1 map, as named in the provenance.')
    units = traits.Str(mandatory=True, desc='How the units of the map were decided.')
    registration = traits.Str(
        mandatory=True,
        desc='How the T1 map was aligned with the anatomical image, as named in the provenance.',
    )
    # None first, so that it is the default. A default_value argument would be stored as trait
    # metadata, which hides the trait's default_value method from nipype's documentation.
    field_strength = traits.Either(
        None,
        traits.Float,
        usedefault=True,
        desc='Field strength of the T1 map, in tesla, to check the measured values against.',
    )


class _MeasureTissueT1OutputSpec(TraitedSpec):
    updates = traits.Dict(
        desc=(
            'GrayMatterT1 and WhiteMatterT1 entries for the calibration specification '
            '(see aslprep.utils.calibration.update_calibration).'
        ),
    )
    warnings = traits.List(traits.Str, desc='Warnings about the measured values.')
    medians = traits.Dict(desc='The medians (see aslprep.utils.t1map.tissue_t1_medians).')


class MeasureTissueT1(SimpleInterface):
    """Measure the T1 of gray and white matter from a T1 map, as medians over pure voxels.

    See :func:`~aslprep.utils.t1map.tissue_t1_medians`.
    """

    input_spec = _MeasureTissueT1InputSpec
    output_spec = _MeasureTissueT1OutputSpec

    def _run_interface(self, runtime):
        from aslprep.utils.calibration import measured_tissue_t1_updates
        from aslprep.utils.t1map import (
            MAX_CSF_PROBABILITY,
            PURE_TISSUE_PROBABILITY,
            implausible_t1_warnings,
            tissue_t1_medians,
        )

        t1_img = nb.load(self.inputs.t1)
        t1 = np.asanyarray(t1_img.dataobj, dtype=np.float64)
        tpms = []
        for tpm_file in self.inputs.tpms[:3]:
            tpm_img = nb.load(tpm_file)
            if tpm_img.shape[:3] != t1.shape or not np.allclose(
                tpm_img.affine,
                t1_img.affine,
                atol=1e-3,
            ):
                raise ValueError(f'Tissue probability map {tpm_file} is not on the T1 grid.')

            tpms.append(np.asanyarray(tpm_img.dataobj, dtype=np.float64))

        medians = tissue_t1_medians(t1, *tpms)
        origin = {
            tissue: (
                f'Median ({medians[tissue]["t1"]:.4g} s, interquartile range '
                f'{medians[tissue]["iqr"]:.3g} s) over {medians[tissue]["n"]} {tissue} matter '
                f'voxels (probability >= {PURE_TISSUE_PROBABILITY:g}, CSF probability < '
                f'{MAX_CSF_PROBABILITY:g}, eroded by one voxel) of the T1 map '
                f'{self.inputs.t1map_name} ({self.inputs.units}), '
                f'{self.inputs.registration}'
            )
            for tissue in ('gray', 'white')
        }
        self._results['updates'] = measured_tissue_t1_updates(
            medians['gray']['t1'],
            medians['white']['t1'],
            origin,
        )
        warnings = implausible_t1_warnings(medians, self.inputs.field_strength)
        for warning in warnings:
            config.loggers.interface.warning(warning)

        self._results['warnings'] = warnings
        self._results['medians'] = medians
        return runtime


class _UnmixTissueR1InputSpec(BaseInterfaceInputSpec):
    t1 = File(exists=True, mandatory=True, desc='T1, in seconds, on the grid of the tissue maps.')
    tpms = traits.List(
        File(exists=True),
        minlen=3,
        mandatory=True,
        desc='Gray matter, white matter and CSF probability maps, in that order.',
    )
    mask = File(exists=True, mandatory=True, desc='Brain mask on the grid of the tissue maps.')
    kernel_mm = traits.Float(mandatory=True, desc="Edge of the regression's cube, in mm.")
    medians = traits.Dict(
        mandatory=True,
        desc=(
            'The median T1s (see aslprep.utils.t1map.tissue_t1_medians), used where there is no '
            'local estimate.'
        ),
    )
    updates = traits.Dict(
        mandatory=True,
        desc='The median T1s, as calibration updates, whose provenance is completed.',
    )


class _UnmixTissueR1OutputSpec(TraitedSpec):
    gm_r1 = File(exists=True, desc='Gray matter R1 (1/s) in each voxel.')
    wm_r1 = File(exists=True, desc='White matter R1 (1/s) in each voxel.')
    status = File(
        exists=True,
        desc='Fit status of each voxel (see aslprep.utils.t1map.R1_FIT_STATUS).',
    )
    updates = traits.Dict(desc='The calibration updates, with the method in their provenance.')
    summary = traits.Dict(desc='How the maps were made, for their sidecars.')


class UnmixTissueR1(SimpleInterface):
    """Estimate the R1 of gray and white matter in each voxel of a T1 map, by local regression.

    See :func:`~aslprep.utils.t1map.local_regression_r1`.
    """

    input_spec = _UnmixTissueR1InputSpec
    output_spec = _UnmixTissueR1OutputSpec

    def _run_interface(self, runtime):
        import copy

        from aslprep.utils.t1map import (
            MAX_CONDITION,
            MAX_CONTAMINATION,
            MAX_R1_RATIO,
            MIN_SUPPORT_FRACTION,
            R1_FIT_STATUS,
            kernel_voxels,
            local_regression_r1,
        )

        t1_img = nb.load(self.inputs.t1)
        t1 = np.asanyarray(t1_img.dataobj, dtype=np.float64)
        images = {}
        for name, in_file in (
            ('gm', self.inputs.tpms[0]),
            ('wm', self.inputs.tpms[1]),
            ('csf', self.inputs.tpms[2]),
            ('mask', self.inputs.mask),
        ):
            img = nb.load(in_file)
            if img.shape[:3] != t1.shape or not np.allclose(img.affine, t1_img.affine, atol=1e-3):
                raise ValueError(f'{in_file} is not on the grid of the T1 map.')

            images[name] = np.asanyarray(img.dataobj, dtype=np.float64)

        kernel = kernel_voxels(self.inputs.kernel_mm, t1_img.affine)
        voxel_size = np.sqrt((t1_img.affine[:3, :3] ** 2).sum(axis=0))
        span = np.asarray(kernel) * voxel_size
        if np.any(span > self.inputs.kernel_mm + voxel_size):
            config.loggers.interface.warning(
                f'The anatomical voxels ({" x ".join(f"{v:.3g}" for v in voxel_size)} mm) are '
                f'too large for a {self.inputs.kernel_mm:g} mm cube, so the T1 map was unmixed in '
                f'{" x ".join(str(k) for k in kernel)} voxels '
                f'({" x ".join(f"{s:.3g}" for s in span)} mm) instead.'
            )

        result = local_regression_r1(
            t1,
            images['gm'],
            images['wm'],
            images['csf'],
            images['mask'] > 0,
            kernel,
            self.inputs.medians,
        )

        header = t1_img.header.copy()
        header.set_data_dtype(np.float32)
        for name in ('gm_r1', 'wm_r1'):
            out_file = os.path.join(runtime.cwd, f'{name}.nii.gz')
            nb.Nifti1Image(result[name].astype(np.float32), t1_img.affine, header).to_filename(
                out_file
            )
            self._results[name] = out_file

        status_file = os.path.join(runtime.cwd, 'r1_fit_status.nii.gz')
        status_img = nb.Nifti1Image(result['status'], t1_img.affine)
        status_img.set_data_dtype(np.uint8)
        status_img.to_filename(status_file)
        self._results['status'] = status_file

        fractions = result['estimated_fraction']
        rule = (
            f'Estimated in each voxel by local linear regression of R1 (1/T1) on the gray matter, '
            f'white matter and CSF probability maps, in a {self.inputs.kernel_mm:g} mm cube '
            f'({" x ".join(str(k) for k in kernel)} voxels); a tissue was left out of a fit only '
            f'when its contamination of the gray and white matter estimates was at most '
            f'{MAX_CONTAMINATION:g}. Where the tissue was not in the fit, or its estimate was more '
            f'than a factor of {MAX_R1_RATIO:g} from its median, the median was used.'
        )
        updates = copy.deepcopy(self.inputs.updates)
        for name, tissue in (('GrayMatterT1', 'gray'), ('WhiteMatterT1', 'white')):
            updates[name]['Origin'] = (
                f'{updates[name]["Origin"]}. The tissue-T1 model used this value only where the '
                f'per-voxel {tissue} matter T1 map did not estimate one: '
                f'{rule} The map estimated {100 * fractions[tissue]:.3g}% of the {tissue} matter '
                '(by probability, within the brain mask).'
            )

        self._results['updates'] = updates
        self._results['summary'] = {
            'Method': 'local-regression',
            'Description': rule,
            'KernelSize': float(self.inputs.kernel_mm),
            'KernelVoxels': list(kernel),
            'MinimumSupportFraction': MIN_SUPPORT_FRACTION,
            'MaximumConditionNumber': MAX_CONDITION,
            'MaximumContamination': MAX_CONTAMINATION,
            'MaximumRatioToMedian': MAX_R1_RATIO,
            'EstimatedFraction': {
                'GrayMatter': fractions['gray'],
                'WhiteMatter': fractions['white'],
            },
            'FitStatus': {str(code): label for code, label in R1_FIT_STATUS.items()},
        }
        return runtime
