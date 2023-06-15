"""Interfaces for calculating CBF."""
import json
import os

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

from aslprep.utils.qc import (
    average_cbf_by_tissue,
    compute_qei,
    coverage,
    crosscorr,
    dice,
    jaccard,
    negativevoxel,
)


class _ComputeCBFQCInputSpec(BaseInterfaceInputSpec):
    name_source = File(
        exists=True,
        mandatory=True,
        desc="Original asl_file. Used to extract entity information.",
    )
    mean_cbf = File(exists=True, mandatory=True, desc="Mean CBF from standard CBF calculation.")
    # SCORE/SCRUB inputs
    mean_cbf_score = File(exists=True, mandatory=False, desc="Mean CBF after SCORE censoring.")
    mean_cbf_scrub = File(exists=True, mandatory=False, desc="Mean CBF after SCRUB denoising.")
    # BASIL inputs
    mean_cbf_basil = File(exists=True, mandatory=False, desc="Mean CBF produced by BASIL.")
    mean_cbf_gm_basil = File(
        exists=True,
        mandatory=False,
        desc="GM partial volume corrected CBF with BASIL.",
    )
    # Tissue probability maps and masks
    gm_tpm = File(exists=True, mandatory=True, desc="Gray matter tissue probability map")
    wm_tpm = File(exists=True, mandatory=True, desc="White matter tissue probability map")
    csf_tpm = File(exists=True, mandatory=True, desc="CSF tissue probability map")
    asl_mask = File(exists=True, mandatory=True, desc="ASL mask in native ASL reference space")
    t1w_mask = File(exists=True, mandatory=True, desc="T1w mask in native space")
    asl_mask_std = File(exists=True, mandatory=False, desc="ASL mask in standard space")
    template_mask = File(exists=True, mandatory=False, desc="template mask or image")
    tpm_threshold = traits.Float(
        default_value=0.7,
        usedefault=True,
        mandatory=False,
        desc="Tissue probability threshold for binarizing GM, WM, and CSF masks.",
    )
    # Non-GE-only inputs
    confounds_file = File(
        exists=True,
        mandatory=False,
        desc="Confounds file. Will not be defined for GE data.",
    )
    rmsd_file = File(
        exists=True,
        mandatory=False,
        desc="RMSD file. Will not be defined for GE data.",
    )


class _ComputeCBFQCOutputSpec(TraitedSpec):
    qc_file = File(exists=False, desc="qc file ")


class ComputeCBFQC(SimpleInterface):
    """Calculate a series of CBF quality control metrics for GE data.

    compute qc from confound regressors
    and cbf maps,
    coregistration and regsitration indexes
    """

    input_spec = _ComputeCBFQCInputSpec
    output_spec = _ComputeCBFQCOutputSpec

    def _run_interface(self, runtime):
        thresh = self.inputs.tpm_threshold

        if isdefined(self.inputs.confounds_file):
            confounds_df = pd.read_table(self.inputs.confounds_file)
            confounds_df.fillna(0, inplace=True)
            mean_fd = np.mean(confounds_df["framewise_displacement"])
            mean_rms = pd.read_csv(self.inputs.rmsd_file, header=None).mean().values[0]
        else:
            mean_fd = np.nan
            mean_rms = np.nan

        coreg_dice = dice(self.inputs.asl_mask, self.inputs.t1w_mask)
        coreg_jaccard = jaccard(self.inputs.asl_mask, self.inputs.t1w_mask)
        coreg_crosscorr = crosscorr(self.inputs.asl_mask, self.inputs.t1w_mask)
        coreg_coverage = coverage(self.inputs.asl_mask, self.inputs.t1w_mask)

        if self.inputs.asl_mask_std and self.inputs.template_mask:
            norm_dice = dice(self.inputs.asl_mask_std, self.inputs.template_mask)
            norm_jaccard = jaccard(self.inputs.asl_mask_std, self.inputs.template_mask)
            norm_crosscorr = crosscorr(self.inputs.asl_mask_std, self.inputs.template_mask)
            norm_coverage = coverage(self.inputs.asl_mask_std, self.inputs.template_mask)

        mean_cbf_qei = compute_qei(
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
            mean_cbf_score_qei = compute_qei(
                gm=self.inputs.gm_tpm,
                wm=self.inputs.wm_tpm,
                csf=self.inputs.csf_tpm,
                img=self.inputs.mean_cbf_score,
                thresh=thresh,
            )
            mean_cbf_scrub_qei = compute_qei(
                gm=self.inputs.gm_tpm,
                wm=self.inputs.wm_tpm,
                csf=self.inputs.csf_tpm,
                img=self.inputs.mean_cbf_scrub,
                thresh=thresh,
            )
            mean_cbf_score_negvox = negativevoxel(
                cbf=self.inputs.mean_cbf_score,
                gm=self.inputs.gm_tpm,
                thresh=thresh,
            )
            mean_cbf_scrub_negvox = negativevoxel(
                cbf=self.inputs.mean_cbf_scrub,
                gm=self.inputs.gm_tpm,
                thresh=thresh,
            )
        else:
            print("no score inputs, setting to np.nan")
            mean_cbf_score_qei = np.nan
            mean_cbf_scrub_qei = np.nan
            mean_cbf_score_negvox = np.nan
            mean_cbf_scrub_negvox = np.nan

        if self.inputs.mean_cbf_basil:
            mean_cbf_basil_qei = compute_qei(
                gm=self.inputs.gm_tpm,
                wm=self.inputs.wm_tpm,
                csf=self.inputs.csf_tpm,
                img=self.inputs.mean_cbf_basil,
                thresh=thresh,
            )
            mean_cbf_gm_basil_qei = compute_qei(
                gm=self.inputs.gm_tpm,
                wm=self.inputs.wm_tpm,
                csf=self.inputs.csf_tpm,
                img=self.inputs.mean_cbf_gm_basil,
                thresh=thresh,
            )
            mean_cbf_basil_negvox = negativevoxel(
                cbf=self.inputs.mean_cbf_basil,
                gm=self.inputs.gm_tpm,
                thresh=thresh,
            )
            mean_cbf_gm_basil_negvox = negativevoxel(
                cbf=self.inputs.mean_cbf_gm_basil,
                gm=self.inputs.gm_tpm,
                thresh=thresh,
            )
        else:
            print("no basil inputs, setting to np.nan")
            mean_cbf_basil_qei = np.nan
            mean_cbf_gm_basil_qei = np.nan
            mean_cbf_basil_negvox = np.nan
            mean_cbf_gm_basil_negvox = np.nan

        gm_wm_ratio = np.divide(mean_cbf_mean[0], mean_cbf_mean[1])
        mean_cbf_negvox = negativevoxel(
            cbf=self.inputs.mean_cbf,
            gm=self.inputs.gm_tpm,
            thresh=thresh,
        )

        metrics_dict = {
            "FD": [mean_fd],
            "rmsd": [mean_rms],
            "coregDC": [coreg_dice],
            "coregJC": [coreg_jaccard],
            "coregCC": [coreg_crosscorr],
            "coregCOV": [coreg_coverage],
            "cbfQEI": [mean_cbf_qei],
            "scoreQEI": [mean_cbf_score_qei],
            "scrubQEI": [mean_cbf_scrub_qei],
            "basilQEI": [mean_cbf_basil_qei],
            "pvcQEI": [mean_cbf_gm_basil_qei],
            "GMmeanCBF": [mean_cbf_mean[0]],
            "WMmeanCBF": [mean_cbf_mean[1]],
            "Gm_Wm_CBF_ratio": [gm_wm_ratio],
            "NEG_CBF_PERC": [mean_cbf_negvox],
            "NEG_SCORE_PERC": [mean_cbf_score_negvox],
            "NEG_SCRUB_PERC": [mean_cbf_scrub_negvox],
            "NEG_BASIL_PERC": [mean_cbf_basil_negvox],
            "NEG_PVC_PERC": [mean_cbf_gm_basil_negvox],
        }

        qc_metadata = {
            "FD": {
                "LongName": "Mean Framewise Displacement",
                "Description": (
                    "Average framewise displacement without any motion parameter filtering. "
                    "This value includes high-motion outliers, but not dummy volumes. "
                    "FD is calculated according to the Power definition."
                ),
                "Units": "mm",
                "Term URL": "https://doi.org/10.1016/j.neuroimage.2011.10.018",
            },
            "rmsd": {
                "LongName": "Mean Relative Root Mean Squared",
                "Description": (
                    "Average relative root mean squared calculated from motion parameters, "
                    "after removal of dummy volumes and high-motion outliers. "
                    "Relative in this case means 'relative to the previous scan'."
                ),
                "Units": "arbitrary",
            },
            "coregDC": {
                "LongName": "Coregistration Sørensen-Dice Coefficient",
                "Description": "",
                "Units": "",
                "Term URL": "https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient",
            },
            "coregJC": {
                "LongName": "Coregistration Jaccard Index",
                "Description": "",
                "Units": "",
                "Term URL": "https://en.wikipedia.org/wiki/Jaccard_index",
            },
            "coregCC": {
                "LongName": "",
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "coregCOV": {
                "LongName": "",
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "cbfQEI": {
                "LongName": "Cerebral Blood Flow Quality Evaluation Index",
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "scoreQEI": {
                "LongName": "SCORE-Denoised Cerebral Blood Flow Quality Evaluation Index",
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "scrubQEI": {
                "LongName": "SCRUB-Denoised Cerebral Blood Flow Quality Evaluation Index",
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "basilQEI": {
                "LongName": "BASIL Cerebral Blood Flow Quality Evaluation Index",
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "pvcQEI": {
                "LongName": (
                    "BASIL Partial Volume Corrected Cerebral Blood Flow Quality Evaluation Index"
                ),
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "GMmeanCBF": {
                "LongName": "Mean Cerebral Blood Flow of Gray Matter",
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "WMmeanCBF": {
                "LongName": "Mean Cerebral Blood Flow of White Matter",
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "Gm_Wm_CBF_ratio": {
                "LongName": "Mean Gray Matter-White Matter Cerebral Blood Flow Ratio",
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "NEG_CBF_PERC": {
                "LongName": "Percentage of Negative Cerebral Blood Flow Values",
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "NEG_SCORE_PERC": {
                "LongName": "Percentage of Negative SCORE-Denoised Cerebral Blood Flow Values",
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "NEG_SCRUB_PERC": {
                "LongName": "Percentage of Negative SCRUB-Denoised Cerebral Blood Flow Values",
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "NEG_BASIL_PERC": {
                "LongName": "Percentage of Negative BASIL Cerebral Blood Flow Values",
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
            "NEG_PVC_PERC": {
                "LongName": (
                    "Percentage of Negative BASIL Partial Volume Corrected Cerebral Blood Flow "
                    "Values"
                ),
                "Description": "",
                "Units": "",
                "Term URL": "",
            },
        }

        if self.inputs.asl_mask_std and self.inputs.template_mask:
            metrics_dict.update(
                {
                    "normDC": [norm_dice],
                    "normJC": [norm_jaccard],
                    "normCC": [norm_crosscorr],
                    "normCOV": [norm_coverage],
                }
            )

            qc_metadata.update(
                {
                    "normDC": {
                        "LongName": "",
                        "Description": "",
                        "Units": "",
                        "Term URL": (
                            "https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient"
                        ),
                    },
                    "normJC": {
                        "LongName": "",
                        "Description": "",
                        "Units": "",
                        "Term URL": "https://en.wikipedia.org/wiki/Jaccard_index",
                    },
                    "normCC": {
                        "LongName": "",
                        "Description": "",
                        "Units": "",
                        "Term URL": "",
                    },
                    "normCOV": {
                        "LongName": "",
                        "Description": "",
                        "Units": "",
                        "Term URL": "",
                    },
                }
            )

        # Extract entities from the input file.
        # Useful for identifying ASL files after concatenating the QC files across runs.
        base_file = os.path.basename(self.inputs.name_source)
        entities = base_file.split("_")[:-1]
        entities_dict = {ent.split("-")[0]: ent.split("-")[1] for ent in entities}

        # Combine the dictionaries and convert to a DataFrame.
        qc_dict = {**entities_dict, **metrics_dict}
        qc_df = pd.DataFrame(qc_dict)

        self._results["qc_file"] = fname_presuffix(
            self.inputs.mean_cbf,
            suffix="qc_cbf.csv",
            newpath=runtime.cwd,
            use_ext=False,
        )
        qc_df.to_csv(self._results["qc_file"], index=False, header=True, na_rep="n/a")

        self._results["qc_metadata"] = fname_presuffix(
            self.inputs.mean_cbf,
            suffix="qc_cbf.json",
            newpath=runtime.cwd,
            use_ext=False,
        )
        with open(self._results["qc_metadata"], "w") as fo:
            json.dump(qc_metadata, fo, indent=4, sort_keys=True)

        return runtime
