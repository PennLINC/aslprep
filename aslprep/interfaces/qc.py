"""Interfaces for calculating CBF."""
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
    tpm_threshold = traits.Float(desc="Typically 0.7 for non-GE data and 0.8 for GE data.")
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

        normalization_metrics_dict = {}
        if self.inputs.asl_mask_std and self.inputs.template_mask:
            normalization_metrics_dict = {
                "normDC": [norm_dice],
                "normJC": [norm_jaccard],
                "normCC": [norm_crosscorr],
                "normCOV": [norm_coverage],
            }

        # Extract entities from the input file.
        # Useful for identifying ASL files after concatenating the QC files across runs.
        base_file = os.path.basename(self.inputs.name_source)
        entities = base_file.split("_")[:-1]
        entities_dict = {ent.split("-")[0]: ent.split("-")[1] for ent in entities}

        # Combine the dictionaries and convert to a DataFrame.
        qc_dict = {**entities_dict, **metrics_dict, **normalization_metrics_dict}
        qc_df = pd.DataFrame(qc_dict)

        self._results["qc_file"] = fname_presuffix(
            self.inputs.mean_cbf,
            suffix="qc_cbf.csv",
            newpath=runtime.cwd,
            use_ext=False,
        )
        qc_df.to_csv(self._results["qc_file"], index=False, header=True, na_rep="n/a")

        return runtime
