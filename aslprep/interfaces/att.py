"""Interfaces that bring an external ATT map onto a single-delay run's grid."""

import os

import nibabel as nb
import numpy as np
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)

from aslprep import config


class _ExternalATTInputSpec(BaseInterfaceInputSpec):
    att = File(exists=True, mandatory=True, desc="The source run's ATT, on its ASL grid.")
    att_status = File(exists=True, mandatory=True, desc="The source run's ATT status.")
    att_ci_lower = File(exists=True, mandatory=True, desc="The source run's ATT interval, low.")
    att_ci_upper = File(exists=True, mandatory=True, desc="The source run's ATT interval, high.")
    source_aslref2anat = File(
        exists=True,
        mandatory=True,
        desc="The source run's ASL-reference-to-anatomical transform (ITK).",
    )
    source_t1w = File(
        exists=True,
        mandatory=True,
        desc='The anatomical image the source run was registered to.',
    )
    t1w_preproc = File(exists=True, mandatory=True, desc="This run's anatomical image.")
    aslref2anat = File(
        exists=True,
        mandatory=True,
        desc="This run's ASL-reference-to-anatomical transform (ITK).",
    )
    reference = File(
        exists=True,
        mandatory=True,
        desc="This run's ASL brain mask, which sets the grid.",
    )
    gm_pv = File(exists=True, mandatory=True, desc="This run's gray matter partial volume.")
    wm_pv = File(exists=True, mandatory=True, desc="This run's white matter partial volume.")
    nominal_cbf = traits.List(
        traits.Float,
        minlen=2,
        maxlen=2,
        mandatory=True,
        desc='Gray and white matter nominal CBFs of the tissue-T1 model, in mL/100 g/min.',
    )
    fill = traits.Enum('median', 'none', usedefault=True, desc='--att-fill.')
    source_name = traits.Str(mandatory=True, desc='The ATT map, as named in the provenance.')


class _ExternalATTOutputSpec(TraitedSpec):
    att = File(exists=True, desc="ATT (s) on this run's grid, zero where there is none.")
    att_status = File(
        exists=True,
        desc='Where the ATT came from (see aslprep.utils.att.EXTERNAL_ATT_STATUS).',
    )
    summary = traits.Dict(desc='The tissue medians and the fraction of voxels with each status.')


class ExternalATT(SimpleInterface):
    """Resample a multi-delay run's ATT onto a single-delay run's grid, and fill it.

    See :func:`~aslprep.utils.att.resample_external_att`.
    """

    input_spec = _ExternalATTInputSpec
    output_spec = _ExternalATTOutputSpec

    def _run_interface(self, runtime):
        from aslprep.utils.att import (
            EXTERNAL_ATT_STATUS,
            check_same_anatomy,
            resample_external_att,
            target_to_source_voxels,
        )

        check_same_anatomy(self.inputs.source_t1w, self.inputs.t1w_preproc)

        source_img = nb.load(self.inputs.att)
        source = {}
        for name in ('att', 'att_status', 'att_ci_lower', 'att_ci_upper'):
            in_file = getattr(self.inputs, name)
            img = nb.load(in_file)
            same_grid = img.shape[:3] == source_img.shape[:3] and np.allclose(
                img.affine,
                source_img.affine,
                atol=1e-3,
            )
            if not same_grid:
                raise ValueError(f'{in_file} is not on the grid of the ATT map.')

            source[name] = np.asanyarray(img.dataobj, dtype=np.float64).reshape(img.shape[:3])

        ref_img = nb.load(self.inputs.reference)
        shape = ref_img.shape[:3]
        target = {}
        for name in ('reference', 'gm_pv', 'wm_pv'):
            img = nb.load(getattr(self.inputs, name))
            if img.shape[:3] != shape or not np.allclose(img.affine, ref_img.affine, atol=1e-3):
                raise ValueError(f'{getattr(self.inputs, name)} is not on the ASL grid.')

            target[name] = np.asanyarray(img.dataobj, dtype=np.float64).reshape(-1)

        coords = target_to_source_voxels(
            ref_img.affine,
            shape,
            source_img.affine,
            self.inputs.aslref2anat,
            self.inputs.source_aslref2anat,
        )
        att, status, summary = resample_external_att(
            source['att'],
            source['att_status'],
            source['att_ci_lower'],
            source['att_ci_upper'],
            coords,
            target['gm_pv'],
            target['wm_pv'],
            target['reference'] > 0,
            tuple(self.inputs.nominal_cbf),
            fill=self.inputs.fill,
        )

        att_file = os.path.join(runtime.cwd, 'external_att.nii.gz')
        att_data = np.nan_to_num(att).reshape(shape).astype(np.float32)
        att_img = nb.Nifti1Image(att_data, ref_img.affine)
        att_img.set_data_dtype(np.float32)
        att_img.to_filename(att_file)
        status_file = os.path.join(runtime.cwd, 'external_attstatus.nii.gz')
        status_img = nb.Nifti1Image(status.reshape(shape), ref_img.affine)
        status_img.set_data_dtype(np.uint8)
        status_img.to_filename(status_file)

        fractions = summary['ResamplingStatusFractions']
        config.loggers.interface.info(
            'External ATT: %.1f%% of the brain mask measured, %.1f%% filled, %.1f%% without ATT.',
            100 * fractions['1'],
            100 * (fractions['2'] + fractions['3']),
            100 * (fractions['0'] + fractions['4'] + fractions['5'] + fractions['7']),
        )
        self._results['att'] = att_file
        self._results['att_status'] = status_file
        self._results['summary'] = {
            'Source': self.inputs.source_name,
            **summary,
            'Fill': self.inputs.fill,
            'StatusLevels': {str(code): label for code, label in EXTERNAL_ATT_STATUS.items()},
        }
        return runtime
