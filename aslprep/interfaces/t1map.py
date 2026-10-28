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
    field_strength = traits.Either(
        traits.Float,
        None,
        default_value=None,
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
        return runtime
