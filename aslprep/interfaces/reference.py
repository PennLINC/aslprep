"""Interfaces for building reference images."""

import nibabel as nb
import pandas as pd
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    CommandLine,
    CommandLineInputSpec,
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


class _BrainChopInputSpec(CommandLineInputSpec):
    in_file = File(
        exists=True,
        mandatory=True,
        desc='Input file.',
        argstr='%s',
        position=0,
    )
    skullstripped_file = File(
        'skullstripped.nii.gz',
        usedefault=True,
        desc='Output file.',
        argstr='-o %s',
        position=1,
    )
    mask_file = File(
        'mask.nii.gz',
        usedefault=True,
        desc='Output file.',
        argstr='--mask %s',
        position=2,
    )


class _BrainChopOutputSpec(TraitedSpec):
    skullstripped_file = File(desc='Output file.', exists=True)
    mask_file = File(desc='Output file.', exists=True)


class BrainChop(CommandLine):
    """Brain chop interface."""

    input_spec = _BrainChopInputSpec
    output_spec = _BrainChopOutputSpec
    _cmd = 'brainchop'
