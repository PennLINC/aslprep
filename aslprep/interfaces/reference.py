"""Interfaces for building reference images."""

import nibabel as nb
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
