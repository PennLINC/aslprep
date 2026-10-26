"""Interfaces for building reference images."""

import nibabel as nb
import numpy as np
import pandas as pd
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    isdefined,
    traits,
)
from nipype.utils.filemanip import fname_presuffix

from aslprep import config


class _SelectHighestContrastVolumesInputSpec(BaseInterfaceInputSpec):
    asl_file = File(exists=True, mandatory=True, desc='ASL file.')
    aslcontext = File(exists=True, mandatory=True, desc='ASL context file.')
    m0scan = File(
        exists=True,
        mandatory=False,
        desc="M0 scan file. Only provided if M0Type is 'separate'.",
    )
    prioritize_m0 = traits.Bool(
        mandatory=True,
        desc='Whether to prioritize the M0 scan (useful for GE data) or not.',
    )


class _SelectHighestContrastVolumesOutputSpec(TraitedSpec):
    selected_volumes_file = File(desc='File containing the highest-contrast available volumes.')


class SelectHighestContrastVolumes(SimpleInterface):
    """Select the highest-contrast volume type for reference volume generation.

    This interface selects the highest-contrast volume type from the ASL context file and
    returns a file containing the highest-contrast volumes.
    The choice of volume type is based on a hard-coded priority list:
    - Included M0 scan (if prioritized)
    - Separate M0 scan (if prioritized)
    - CBF
    - DeltaM
    - Included M0 scan (if not prioritized)
    - Separate M0 scan (if not prioritized)
    - Control
    """

    input_spec = _SelectHighestContrastVolumesInputSpec
    output_spec = _SelectHighestContrastVolumesOutputSpec

    def _run_interface(self, runtime):
        asl_img = nb.load(self.inputs.asl_file)
        aslcontext_df = pd.read_table(self.inputs.aslcontext)
        if aslcontext_df.shape[0] != asl_img.shape[3]:
            raise ValueError(
                f'Number of volumes in {self.inputs.asl_file} ({asl_img.shape[3]}) '
                f'does not match the number of rows in {self.inputs.aslcontext} '
                f'({aslcontext_df.shape[0]})'
            )

        if 'm0scan' in aslcontext_df['volume_type'].tolist() and self.inputs.prioritize_m0:
            target_type = 'm0scan'
        elif isdefined(self.inputs.m0scan) and self.inputs.prioritize_m0:
            target_type = 'separate_m0scan'
        elif 'cbf' in aslcontext_df['volume_type'].tolist():
            target_type = 'cbf'
        elif 'deltam' in aslcontext_df['volume_type'].tolist():
            target_type = 'deltam'
        elif 'm0scan' in aslcontext_df['volume_type'].tolist():
            target_type = 'm0scan'
        elif isdefined(self.inputs.m0scan):
            target_type = 'separate_m0scan'
        else:
            target_type = 'control'

        config.loggers.interface.info(
            f'Selecting {target_type} as highest-contrast volume type for reference volume '
            'generation.'
        )

        if target_type == 'separate_m0scan':
            self._results['selected_volumes_file'] = self.inputs.m0scan
            return runtime

        # Otherwise, split up the ASL file based on the volume type with the highest contrast.
        target_idx = aslcontext_df.loc[aslcontext_df['volume_type'] == target_type].index.values
        if target_idx.size == 0:
            raise ValueError(f"Volume type '{target_type}' missing from {self.inputs.aslcontext}")

        asl_data = asl_img.get_fdata()
        highest_contrast_data = asl_data[:, :, :, target_idx]

        self._results['selected_volumes_file'] = fname_presuffix(
            self.inputs.asl_file,
            suffix='_contrast.nii.gz',
            newpath=runtime.cwd,
            use_ext=False,
        )

        nb.Nifti1Image(
            dataobj=highest_contrast_data,
            affine=asl_img.affine,
            header=asl_img.header,
        ).to_filename(self._results['selected_volumes_file'])

        return runtime


#: Sidecar metadata of the anatomical brain mask's quality-control measures.
ANAT_MASK_QC_METADATA = {
    'anat_mask_synthstrip_dice': {
        'LongName': 'Dice Between the Anatomical and SynthStrip ASL Brain Masks',
        'Description': (
            'Sørensen-Dice coefficient between the anatomical brain mask, resampled to the ASL '
            'reference through the coregistration, and the SynthStrip mask of the ASL reference. '
            'Low values mean that one of the two masks is wrong: a failed coregistration '
            '(the anatomical mask is displaced) or a failed SynthStrip mask.'
        ),
        'Term URL': 'https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient',
    },
    'anat_mask_synthstrip_volume_ratio': {
        'LongName': 'SynthStrip-to-Anatomical ASL Brain Mask Volume Ratio',
        'Description': (
            'Volume of the SynthStrip mask of the ASL reference divided by the volume of the '
            'anatomical brain mask resampled to the ASL reference.'
        ),
    },
    'anat_mask_used': {
        'LongName': 'Anatomical ASL Brain Mask Used',
        'Description': (
            '1 if the anatomical brain mask was used as the ASL brain mask, 0 if the SynthStrip '
            'mask was used instead because the two masks disagreed in a way that indicates a '
            'failed coregistration.'
        ),
    },
}


class _GuardAnatomicalMaskInputSpec(BaseInterfaceInputSpec):
    anat_mask = File(
        exists=True,
        mandatory=True,
        desc='Anatomical brain mask, resampled to the ASL reference.',
    )
    asl_mask = File(
        exists=True,
        mandatory=True,
        desc='Brain mask estimated from the ASL reference (SynthStrip), on the same grid.',
    )
    min_dice = traits.Float(
        0.9,
        usedefault=True,
        desc='Below this Dice between the two masks, they disagree.',
    )
    min_volume_ratio = traits.Float(
        0.95,
        usedefault=True,
        desc=(
            'When the masks disagree, the ASL mask is used if its volume is at least this '
            'fraction of the anatomical mask volume.'
        ),
    )


class _GuardAnatomicalMaskOutputSpec(TraitedSpec):
    out_mask = File(exists=True, desc='The brain mask to use.')
    qc = traits.Dict(desc='Quality-control measures (see ANAT_MASK_QC_METADATA).')
    metadata = traits.Dict(desc='Sidecar metadata of the brain mask.')


class GuardAnatomicalMask(SimpleInterface):
    """Use the anatomical brain mask, unless it disagrees with the ASL-derived mask.

    The anatomical mask is only as good as the coregistration. When the two masks disagree
    (Dice below ``min_dice``), one of them is wrong, and their volumes tell which:
    a failed SynthStrip mask misses brain, so it is the smaller one,
    while a failed coregistration displaces the anatomical mask (often partly out of the
    field of view) without shrinking the ASL mask.
    So when they disagree, the ASL mask is used if its volume is at least ``min_volume_ratio``
    of the anatomical mask's, and the anatomical mask otherwise.
    """

    input_spec = _GuardAnatomicalMaskInputSpec
    output_spec = _GuardAnatomicalMaskOutputSpec

    def _run_interface(self, runtime):
        from aslprep.utils.confounds import dice

        anat_img = nb.load(self.inputs.anat_mask)
        anat = np.asanyarray(anat_img.dataobj) > 0
        asl = np.asanyarray(nb.load(self.inputs.asl_mask).dataobj) > 0
        if anat.shape != asl.shape:
            raise ValueError(
                f'The anatomical {anat.shape} and ASL {asl.shape} masks are not on one grid.'
            )

        overlap = float(dice(anat, asl))
        ratio = float(asl.sum() / anat.sum()) if anat.any() else np.inf
        use_anat = overlap >= self.inputs.min_dice or ratio < self.inputs.min_volume_ratio
        if not use_anat:
            config.loggers.interface.warning(
                'The anatomical brain mask disagrees with the SynthStrip mask of the ASL '
                'reference (Dice %.3f; SynthStrip/anatomical volume %.2f), which indicates a '
                'failed coregistration. Using the SynthStrip mask. Check the coregistration.',
                overlap,
                ratio,
            )

        mask = anat if use_anat else asl
        self._results['out_mask'] = fname_presuffix(
            self.inputs.anat_mask, suffix='_guarded', newpath=runtime.cwd
        )
        out_img = nb.Nifti1Image(mask.astype(np.uint8), anat_img.affine, anat_img.header)
        out_img.set_data_dtype(np.uint8)
        out_img.to_filename(self._results['out_mask'])

        self._results['qc'] = {
            'anat_mask_synthstrip_dice': overlap,
            'anat_mask_synthstrip_volume_ratio': ratio,
            'anat_mask_used': int(use_anat),
        }
        self._results['metadata'] = {
            'SkullStripASL': 'anatomical',
            'BrainMaskSource': 'anatomical' if use_anat else 'synthstrip',
            'AnatomicalSynthStripDice': overlap,
        }
        return runtime
