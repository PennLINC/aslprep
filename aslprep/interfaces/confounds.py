# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Interfaces for calculating and collecting confounds."""

import json
import os

import nibabel as nb
import numpy as np
import pandas as pd
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.utils.filemanip import fname_presuffix
from nipype.utils.misc import normalize_mc_params

from aslprep.utils.confounds import (
    _gather_confounds,
    average_cbf_by_tissue,
    compute_qei,
    dice,
    overlap,
    pearson,
)


class _GatherConfoundsInputSpec(BaseInterfaceInputSpec):
    signals = File(exists=True, desc='input signals')
    dvars = File(exists=True, desc='file containing DVARS')
    rmsd = File(exists=True, desc='input RMS framewise displacement')
    std_dvars = File(exists=True, desc='file containing standardized DVARS')
    fd = File(exists=True, desc='input framewise displacement')
    motion = File(exists=True, desc='input motion parameters')


class _GatherConfoundsOutputSpec(TraitedSpec):
    confounds_file = File(exists=True, desc='output confounds file')
    confounds_list = traits.List(traits.Str, desc='list of headers')


class GatherConfounds(SimpleInterface):
    """Combine various sources of confounds in one TSV file."""

    input_spec = _GatherConfoundsInputSpec
    output_spec = _GatherConfoundsOutputSpec

    def _run_interface(self, runtime):
        combined_out, confounds_list = _gather_confounds(
            signals=self.inputs.signals,
            dvars=self.inputs.dvars,
            rmsd=self.inputs.rmsd,
            std_dvars=self.inputs.std_dvars,
            fdisp=self.inputs.fd,
            motion=self.inputs.motion,
            newpath=runtime.cwd,
        )
        self._results['confounds_file'] = combined_out
        self._results['confounds_list'] = confounds_list
        return runtime


class _GatherCBFConfoundsInputSpec(BaseInterfaceInputSpec):
    signals = File(exists=True, desc='input signals')
    score = File(exists=True, desc='SCORE outlier index')


class _GatherCBFConfoundsOutputSpec(TraitedSpec):
    confounds_file = File(exists=True, desc='output confounds file')
    confounds_list = traits.List(traits.Str, desc='list of headers')


class GatherCBFConfounds(SimpleInterface):
    """Combine various sources of confounds in one TSV file."""

    input_spec = _GatherCBFConfoundsInputSpec
    output_spec = _GatherCBFConfoundsOutputSpec

    def _run_interface(self, runtime):
        combined_out, confounds_list = _gather_confounds(
            signals=self.inputs.signals,
            dvars=None,
            rmsd=None,
            std_dvars=None,
            fdisp=None,
            motion=None,
            score=self.inputs.score,
            newpath=runtime.cwd,
        )
        self._results['confounds_file'] = combined_out
        self._results['confounds_list'] = confounds_list
        return runtime


class _NormalizeMotionParamsInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='the input parameters file')
    format = traits.Enum('FSL', 'AFNI', 'FSFAST', 'NIPY', usedefault=True, desc='output format')


class _NormalizeMotionParamsOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='written file path')


class NormalizeMotionParams(SimpleInterface):
    """Convert input motion parameters into the designated convention.

    From niworkflows.interfaces.confounds.NormalizeMotionParams.

    NOTE: Modified by Taylor Salo to support 1D arrays.
    """

    input_spec = _NormalizeMotionParamsInputSpec
    output_spec = _NormalizeMotionParamsOutputSpec

    def _run_interface(self, runtime):
        mpars = np.loadtxt(self.inputs.in_file)  # mpars is N_t x 6

        # Support single-volume motion parameters
        if mpars.ndim == 1:
            mpars = mpars[None, :]

        mpars = np.apply_along_axis(
            func1d=normalize_mc_params, axis=1, arr=mpars, source=self.inputs.format
        )
        self._results['out_file'] = os.path.join(runtime.cwd, 'motion_params.txt')
        np.savetxt(self._results['out_file'], mpars)
        return runtime


class _ComputeCBFQCInputSpec(BaseInterfaceInputSpec):
    name_source = File(
        exists=True,
        mandatory=True,
        desc='Original asl_file. Used to extract entity information.',
    )
    mean_cbf = File(exists=True, mandatory=True, desc='Mean CBF from standard CBF calculation.')
    # SCORE/SCRUB inputs
    mean_cbf_score = File(exists=True, mandatory=False, desc='Mean CBF after SCORE censoring.')
    mean_cbf_scrub = File(exists=True, mandatory=False, desc='Mean CBF after SCRUB denoising.')
    # BASIL inputs
    mean_cbf_basil = File(exists=True, mandatory=False, desc='Mean CBF produced by BASIL.')
    mean_cbf_gm_basil = File(
        exists=True,
        mandatory=False,
        desc='GM partial volume corrected CBF with BASIL.',
    )
    # Tissue probability maps and masks
    gm_tpm = File(exists=True, mandatory=True, desc='Gray matter tissue probability map')
    wm_tpm = File(exists=True, mandatory=True, desc='White matter tissue probability map')
    csf_tpm = File(exists=True, mandatory=True, desc='CSF tissue probability map')
    asl_mask = File(exists=True, mandatory=True, desc='ASL mask in native ASL reference space')
    t1w_mask = File(exists=True, mandatory=True, desc='T1w mask in native space')
    asl_mask_std = File(exists=True, mandatory=False, desc='ASL mask in standard space')
    template_mask = File(exists=True, mandatory=False, desc='template mask or image')
    tpm_threshold = traits.Float(
        default_value=0.7,
        usedefault=True,
        mandatory=False,
        desc='Tissue probability threshold for binarizing GM, WM, and CSF masks.',
    )
    # Non-GE-only inputs
    confounds_file = File(
        exists=True,
        mandatory=False,
        desc='Confounds file. Will not be defined for GE data.',
    )
    rmsd_file = File(
        exists=True,
        mandatory=False,
        desc='RMSD file. Will not be defined for GE data.',
    )


class _ComputeCBFQCOutputSpec(TraitedSpec):
    qc_file = File(exists=True, desc='qc file')
    qc_metadata = File(exists=True, desc='qc metadata')


class ComputeCBFQC(SimpleInterface):
    """Calculate a series of CBF quality control metrics for GE data.

    compute qc from confound regressors
    and cbf maps,
    coregistration and registration indexes
    """

    input_spec = _ComputeCBFQCInputSpec
    output_spec = _ComputeCBFQCOutputSpec

    def _run_interface(self, runtime):
        thresh = self.inputs.tpm_threshold

        confounds_df = pd.read_table(self.inputs.confounds_file)
        confounds_df.fillna(0, inplace=True)
        if 'framewise_displacement' in confounds_df.columns:
            # FD and RMSD only available for multi-volume datasets
            mean_fd = np.mean(confounds_df['framewise_displacement'])
            mean_rms = np.mean(confounds_df['rmsd'])
        else:
            mean_fd = np.nan
            mean_rms = np.nan

        asl_mask_arr = nb.load(self.inputs.asl_mask).get_fdata()
        t1w_mask_arr = nb.load(self.inputs.t1w_mask).get_fdata()
        coreg_dice = dice(asl_mask_arr, t1w_mask_arr)
        coreg_correlation = pearson(asl_mask_arr, t1w_mask_arr)
        coreg_overlap = overlap(asl_mask_arr, t1w_mask_arr)

        if self.inputs.asl_mask_std and self.inputs.template_mask:
            asl_mask_std_arr = nb.load(self.inputs.asl_mask_std).get_fdata()
            template_mask_arr = nb.load(self.inputs.template_mask).get_fdata()
            norm_dice = dice(asl_mask_std_arr, template_mask_arr)
            norm_correlation = pearson(asl_mask_std_arr, template_mask_arr)
            norm_overlap = overlap(asl_mask_std_arr, template_mask_arr)

        qei_cbf, qei_cbf_rho_ss, qei_cbf_di, percentage_negative_cbf = compute_qei(
            gm=self.inputs.gm_tpm,
            wm=self.inputs.wm_tpm,
            csf=self.inputs.csf_tpm,
            img=self.inputs.mean_cbf,
            thresh=thresh,
        )
        mean_cbf_mean = average_cbf_by_tissue(
            gm=self.inputs.gm_tpm,
            wm=self.inputs.wm_tpm,
            csf=self.inputs.csf_tpm,
            cbf=self.inputs.mean_cbf,
            thresh=thresh,
        )

        if self.inputs.mean_cbf_score:
            (
                qei_cbf_score,
                qei_cbf_score_rho_ss,
                qei_cbf_score_di,
                percentage_negative_cbf_score,
            ) = compute_qei(
                gm=self.inputs.gm_tpm,
                wm=self.inputs.wm_tpm,
                csf=self.inputs.csf_tpm,
                img=self.inputs.mean_cbf_score,
                thresh=thresh,
            )
            (
                qei_cbf_scrub,
                qei_cbf_scrub_rho_ss,
                qei_cbf_scrub_di,
                percentage_negative_cbf_scrub,
            ) = compute_qei(
                gm=self.inputs.gm_tpm,
                wm=self.inputs.wm_tpm,
                csf=self.inputs.csf_tpm,
                img=self.inputs.mean_cbf_scrub,
                thresh=thresh,
            )
        else:
            print('no score inputs, setting to np.nan')
            qei_cbf_score = np.nan
            qei_cbf_score_rho_ss = np.nan
            qei_cbf_score_di = np.nan
            percentage_negative_cbf_score = np.nan
            qei_cbf_scrub = np.nan
            qei_cbf_scrub_rho_ss = np.nan
            qei_cbf_scrub_di = np.nan
            percentage_negative_cbf_scrub = np.nan

        if self.inputs.mean_cbf_basil:
            (
                qei_cbf_basil,
                qei_cbf_basil_rho_ss,
                qei_cbf_basil_di,
                percentage_negative_cbf_basil,
            ) = compute_qei(
                gm=self.inputs.gm_tpm,
                wm=self.inputs.wm_tpm,
                csf=self.inputs.csf_tpm,
                img=self.inputs.mean_cbf_basil,
                thresh=thresh,
            )
            (
                qei_cbf_basil_gm,
                qei_cbf_basil_gm_rho_ss,
                qei_cbf_basil_gm_di,
                percentage_negative_cbf_basil_gm,
            ) = compute_qei(
                gm=self.inputs.gm_tpm,
                wm=self.inputs.wm_tpm,
                csf=self.inputs.csf_tpm,
                img=self.inputs.mean_cbf_gm_basil,
                thresh=thresh,
            )
        else:
            print('no basil inputs, setting to np.nan')
            qei_cbf_basil = np.nan
            qei_cbf_basil_rho_ss = np.nan
            qei_cbf_basil_di = np.nan
            qei_cbf_basil_gm = np.nan
            qei_cbf_basil_gm_rho_ss = np.nan
            qei_cbf_basil_gm_di = np.nan
            percentage_negative_cbf_basil = np.nan
            percentage_negative_cbf_basil_gm = np.nan

        ratio_gm_wm_cbf = np.divide(mean_cbf_mean[0], mean_cbf_mean[1])

        metrics_dict = {
            # Motion metrics
            'mean_fd': [mean_fd],
            'rmsd': [mean_rms],
            # Coregistration metrics
            'coreg_dice': [coreg_dice],
            'coreg_correlation': [coreg_correlation],
            'coreg_overlap': [coreg_overlap],
            # Standard CBF QEI metrics
            'qei_cbf': [qei_cbf],
            'qei_cbf_rho_ss': [qei_cbf_rho_ss],
            'qei_cbf_di': [qei_cbf_di],
            'percentage_negative_cbf': [percentage_negative_cbf],
            # SCORE QEI metrics
            'qei_cbf_score': [qei_cbf_score],
            'qei_cbf_score_rho_ss': [qei_cbf_score_rho_ss],
            'qei_cbf_score_di': [qei_cbf_score_di],
            'percentage_negative_cbf_score': [percentage_negative_cbf_score],
            # SCRUB QEI metrics
            'qei_cbf_scrub': [qei_cbf_scrub],
            'qei_cbf_scrub_rho_ss': [qei_cbf_scrub_rho_ss],
            'qei_cbf_scrub_di': [qei_cbf_scrub_di],
            'percentage_negative_cbf_scrub': [percentage_negative_cbf_scrub],
            # BASIL QEI metrics
            'qei_cbf_basil': [qei_cbf_basil],
            'qei_cbf_basil_rho_ss': [qei_cbf_basil_rho_ss],
            'qei_cbf_basil_di': [qei_cbf_basil_di],
            'percentage_negative_cbf_basil': [percentage_negative_cbf_basil],
            # BASIL Partial Volume Corrected QEI metrics
            'qei_cbf_basil_gm': [qei_cbf_basil_gm],
            'qei_cbf_basil_gm_rho_ss': [qei_cbf_basil_gm_rho_ss],
            'qei_cbf_basil_gm_di': [qei_cbf_basil_gm_di],
            'percentage_negative_cbf_basil_gm': [percentage_negative_cbf_basil_gm],
            # Mean CBF values and ratios
            'mean_gm_cbf': [mean_cbf_mean[0]],
            'mean_wm_cbf': [mean_cbf_mean[1]],
            'ratio_gm_wm_cbf': [ratio_gm_wm_cbf],
        }

        qc_metadata = {
            # Motion metrics
            'mean_fd': {
                'LongName': 'Mean Framewise Displacement',
                'Description': (
                    'Average framewise displacement without any motion parameter filtering. '
                    'This value includes high-motion outliers. '
                    'FD is calculated according to the Power definition.'
                ),
                'Units': 'mm',
                'Term URL': 'https://doi.org/10.1016/j.neuroimage.2011.10.018',
            },
            'rmsd': {
                'LongName': 'Mean Relative Root Mean Squared',
                'Description': (
                    'Average relative root mean squared calculated from motion parameters. '
                    "Relative in this case means 'relative to the previous scan'."
                ),
                'Units': 'arbitrary',
            },
            # Coregistration metrics
            'coreg_dice': {
                'LongName': 'Coregistration Sørensen-Dice Coefficient',
                'Description': (
                    'The Sørensen-Dice coefficient calculated between the binary brain masks from '
                    'the coregistered anatomical and ASL reference images. '
                    'Values are bounded between 0 and 1, '
                    'with higher values indicating better coregistration.'
                ),
                'Term URL': 'https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient',
            },
            'coreg_correlation': {
                'LongName': 'Coregistration Pearson Correlation',
                'Description': (
                    'The Pearson correlation coefficient calculated between the binary brain '
                    'masks from the coregistered anatomical and ASL reference images. '
                    'Values are bounded between -1 and 1, '
                    'with higher values indicating better coregistration.'
                ),
                'Term URL': 'https://en.wikipedia.org/wiki/Pearson_correlation_coefficient',
            },
            'coreg_overlap': {
                'LongName': 'Coregistration Overlap Coefficient',
                'Description': (
                    'The Szymkiewicz-Simpson overlap coefficient calculated between the binary '
                    'brain masks from the coregistered anatomical and ASL reference images. '
                    'Higher values indicate better normalization.'
                ),
                'Term URL': 'https://en.wikipedia.org/wiki/Overlap_coefficient',
            },
            # Standard CBF QEI metrics
            'qei_cbf': {
                'LongName': 'Cerebral Blood Flow Quality Evaluation Index',
                'Description': 'QEI calculated on mean CBF image.',
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'qei_cbf_rho_ss': {
                'LongName': 'Cerebral Blood Flow Quality Evaluation Index Rho SS',
                'Description': 'QEI rho SS calculated on mean CBF image.',
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'qei_cbf_di': {
                'LongName': 'Cerebral Blood Flow Quality Evaluation Index DI',
                'Description': 'QEI DI calculated on mean CBF image.',
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'percentage_negative_cbf': {
                'LongName': 'Percentage of Negative Cerebral Blood Flow Values',
                'Description': (
                    'Percentage of negative CBF values, calculated on the mean CBF image.'
                ),
                'Units': 'percent',
            },
            # SCORE QEI metrics
            'qei_cbf_score': {
                'LongName': 'SCORE-Denoised Cerebral Blood Flow Quality Evaluation Index',
                'Description': 'QEI calculated on mean SCORE-denoised CBF image.',
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'qei_cbf_score_rho_ss': {
                'LongName': 'SCORE-Denoised Cerebral Blood Flow Quality Evaluation Index Rho SS',
                'Description': 'QEI rho SS calculated on mean SCORE-denoised CBF image.',
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'qei_cbf_score_di': {
                'LongName': 'SCORE-Denoised Cerebral Blood Flow Quality Evaluation Index DI',
                'Description': 'QEI DI calculated on mean SCORE-denoised CBF image.',
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'percentage_negative_cbf_score': {
                'LongName': 'Percentage of Negative SCORE-Denoised Cerebral Blood Flow Values',
                'Description': (
                    'Percentage of negative CBF values, calculated on the SCORE-denoised '
                    'CBF image.'
                ),
                'Units': 'percent',
            },
            # SCRUB QEI metrics
            'qei_cbf_scrub': {
                'LongName': 'SCRUB-Denoised Cerebral Blood Flow Quality Evaluation Index',
                'Description': 'QEI calculated on mean SCRUB-denoised CBF image.',
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'qei_cbf_scrub_rho_ss': {
                'LongName': 'SCRUB-Denoised Cerebral Blood Flow Quality Evaluation Index Rho SS',
                'Description': 'QEI rho SS calculated on mean SCRUB-denoised CBF image.',
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'qei_cbf_scrub_di': {
                'LongName': 'SCRUB-Denoised Cerebral Blood Flow Quality Evaluation Index DI',
                'Description': 'QEI DI calculated on mean SCRUB-denoised CBF image.',
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'percentage_negative_cbf_scrub': {
                'LongName': 'Percentage of Negative SCRUB-Denoised Cerebral Blood Flow Values',
                'Description': (
                    'Percentage of negative CBF values, calculated on the SCRUB-denoised '
                    'CBF image.'
                ),
                'Units': 'percent',
            },
            # BASIL QEI metrics
            'qei_cbf_basil': {
                'LongName': 'BASIL Cerebral Blood Flow Quality Evaluation Index',
                'Description': 'QEI calculated on CBF image produced by BASIL.',
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'qei_cbf_basil_rho_ss': {
                'LongName': 'BASIL Cerebral Blood Flow Quality Evaluation Index Rho SS',
                'Description': 'QEI rho SS calculated on CBF image produced by BASIL.',
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'qei_cbf_basil_di': {
                'LongName': 'BASIL Cerebral Blood Flow Quality Evaluation Index DI',
                'Description': 'QEI DI calculated on CBF image produced by BASIL.',
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'percentage_negative_cbf_basil': {
                'LongName': 'Percentage of Negative BASIL Cerebral Blood Flow Values',
                'Description': (
                    'Percentage of negative CBF values, calculated on CBF image produced by BASIL.'
                ),
                'Units': 'percent',
            },
            # BASIL Partial Volume Corrected QEI metrics
            'qei_cbf_basil_gm': {
                'LongName': (
                    'BASIL Partial Volume Corrected Cerebral Blood Flow Quality Evaluation Index'
                ),
                'Description': (
                    'QEI calculated on partial volume-corrected CBF image produced by BASIL.'
                ),
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'qei_cbf_basil_gm_rho_ss': {
                'LongName': (
                    'BASIL Partial Volume Corrected Cerebral Blood Flow Quality Evaluation '
                    'Index Rho SS'
                ),
                'Description': (
                    'QEI rho SS calculated on partial volume-corrected CBF image produced by '
                    'BASIL.'
                ),
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'qei_cbf_basil_gm_di': {
                'LongName': (
                    'BASIL Partial Volume Corrected Cerebral Blood Flow Quality Evaluation '
                    'Index DI'
                ),
                'Description': (
                    'QEI DI calculated on partial volume-corrected CBF image produced by BASIL.'
                ),
                'Term URL': 'http://indexsmart.mirasmart.com/ISMRM2017/PDFfiles/0682.html',
            },
            'percentage_negative_cbf_basil_gm': {
                'LongName': (
                    'Percentage of Negative BASIL Partial Volume Corrected Cerebral Blood Flow '
                    'Values'
                ),
                'Description': (
                    'Percentage of negative CBF values, calculated on partial volume-corrected '
                    'CBF image produced by BASIL.'
                ),
                'Units': 'percent',
            },
            # Mean CBF values and ratios
            'mean_gm_cbf': {
                'LongName': 'Mean Cerebral Blood Flow of Gray Matter',
                'Description': 'Mean CBF value of gray matter.',
                'Units': 'mL/100 g/min',
            },
            'mean_wm_cbf': {
                'LongName': 'Mean Cerebral Blood Flow of White Matter',
                'Description': 'Mean CBF value of white matter.',
                'Units': 'mL/100 g/min',
            },
            'ratio_gm_wm_cbf': {
                'LongName': 'Mean Gray Matter-White Matter Cerebral Blood Flow Ratio',
                'Description': (
                    'The ratio between the mean gray matter and mean white matter CBF values.'
                ),
            },
        }

        if self.inputs.asl_mask_std and self.inputs.template_mask:
            metrics_dict.update(
                {
                    'norm_dice': [norm_dice],
                    'norm_correlation': [norm_correlation],
                    'norm_overlap': [norm_overlap],
                }
            )

            qc_metadata.update(
                {
                    'norm_dice': {
                        'LongName': 'Normalization Sørensen-Dice Coefficient',
                        'Description': (
                            'The Sørensen-Dice coefficient calculated between the binary brain '
                            'masks from the normalized ASL reference image and the associated '
                            'template. '
                            'Values are bounded between 0 and 1, '
                            'with higher values indicating better normalization.'
                        ),
                        'Term URL': (
                            'https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient'
                        ),
                    },
                    'norm_correlation': {
                        'LongName': 'Normalization Pearson Correlation',
                        'Description': (
                            'The Pearson correlation coefficient calculated between the binary '
                            'brain masks from the normalized ASL reference image and the '
                            'associated template. '
                            'Values are bounded between -1 and 1, '
                            'with higher values indicating better coregistration.'
                        ),
                        'Term URL': (
                            'https://en.wikipedia.org/wiki/Pearson_correlation_coefficient'
                        ),
                    },
                    'norm_overlap': {
                        'LongName': 'Normalization Overlap Coefficient',
                        'Description': (
                            'The Szymkiewicz-Simpson overlap coefficient calculated between the '
                            'binary brain masks from the normalized ASL reference image and the '
                            'associated template. '
                            'Higher values indicate better normalization.'
                        ),
                        'Term URL': 'https://en.wikipedia.org/wiki/Overlap_coefficient',
                    },
                }
            )

        # Extract entities from the input file.
        # Useful for identifying ASL files after concatenating the QC files across runs.
        base_file = os.path.basename(self.inputs.name_source)
        entities = base_file.split('_')[:-1]
        entities_dict = {ent.split('-')[0]: ent.split('-')[1] for ent in entities}

        # Combine the dictionaries and convert to a DataFrame.
        qc_dict = {**entities_dict, **metrics_dict}
        qc_df = pd.DataFrame(qc_dict)

        self._results['qc_file'] = fname_presuffix(
            self.inputs.mean_cbf,
            suffix='qc_cbf.tsv',
            newpath=runtime.cwd,
            use_ext=False,
        )
        qc_df.to_csv(self._results['qc_file'], index=False, header=True, sep='\t', na_rep='n/a')

        self._results['qc_metadata'] = fname_presuffix(
            self.inputs.mean_cbf,
            suffix='qc_cbf.json',
            newpath=runtime.cwd,
            use_ext=False,
        )
        with open(self._results['qc_metadata'], 'w') as fobj:
            json.dump(qc_metadata, fobj, indent=4, sort_keys=True)

        return runtime
