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
    cbf_qei,
    coverage,
    crosscorr,
    dice,
    globalcbf,
    jaccard,
    negativevoxel,
)


class _ComputeCBFQCInputSpec(BaseInterfaceInputSpec):
    name_source = File(
        exists=True,
        mandatory=True,
        desc="Original asl_file. Used to extract entity information.",
    )
    in_meancbf = File(exists=True, mandatory=True, desc="cbf img")
    in_avgscore = File(exists=True, mandatory=False, desc="cbf img")
    in_scrub = File(exists=True, mandatory=False, desc="cbf img")
    in_basil = File(exists=True, mandatory=False, desc="cbf img")
    in_pvc = File(exists=True, mandatory=False, desc="cbf img")
    in_greyM = File(exists=True, mandatory=True, desc="grey matter tissue probability map")
    in_whiteM = File(exists=True, mandatory=True, desc="white matter tissue probability map")
    in_csf = File(exists=True, mandatory=True, desc="csf tissue probability map")
    in_aslmask = File(exists=True, mandatory=True, desc="asl mask in native space")
    in_t1mask = File(exists=True, mandatory=True, desc="t1wmask in native space ")
    in_aslmaskstd = File(exists=True, mandatory=False, desc="asl mask in native space")
    in_templatemask = File(exists=True, mandatory=False, desc="template mask or image")
    tpm_threshold = traits.Float(desc="Typically 0.7 for non-GE data and 0.8 for GE data.")
    # non-GE-only inputs
    in_confmat = File(
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

        if isdefined(self.inputs.in_confmat):
            time1 = pd.read_table(self.inputs.in_confmat)
            time1.fillna(0, inplace=True)
            fd = np.mean(time1["framewise_displacement"])
            rms = pd.read_csv(self.inputs.rmsd_file, header=None).mean().values[0]
        else:
            fd = np.nan
            rms = np.nan

        regDC = dice(self.inputs.in_aslmask, self.inputs.in_t1mask)
        regJC = jaccard(self.inputs.in_aslmask, self.inputs.in_t1mask)
        regCC = crosscorr(self.inputs.in_aslmask, self.inputs.in_t1mask)
        regCov = coverage(self.inputs.in_aslmask, self.inputs.in_t1mask)

        if self.inputs.in_aslmaskstd and self.inputs.in_templatemask:
            normDC = dice(self.inputs.in_aslmaskstd, self.inputs.in_templatemask)
            normJC = jaccard(self.inputs.in_aslmaskstd, self.inputs.in_templatemask)
            normCC = crosscorr(self.inputs.in_aslmaskstd, self.inputs.in_templatemask)
            normCov = coverage(self.inputs.in_aslmaskstd, self.inputs.in_templatemask)

        meancbf_qei = cbf_qei(
            gm=self.inputs.in_greyM,
            wm=self.inputs.in_whiteM,
            csf=self.inputs.in_csf,
            img=self.inputs.in_meancbf,
            thresh=thresh,
        )
        meancbf = globalcbf(
            gm=self.inputs.in_greyM,
            wm=self.inputs.in_whiteM,
            csf=self.inputs.in_csf,
            cbf=self.inputs.in_meancbf,
            thresh=thresh,
        )

        if self.inputs.in_avgscore:
            scorecbf_qei = cbf_qei(
                gm=self.inputs.in_greyM,
                wm=self.inputs.in_whiteM,
                csf=self.inputs.in_csf,
                img=self.inputs.in_avgscore,
                thresh=thresh,
            )
            scrub_qei = cbf_qei(
                gm=self.inputs.in_greyM,
                wm=self.inputs.in_whiteM,
                csf=self.inputs.in_csf,
                img=self.inputs.in_scrub,
                thresh=thresh,
            )
            negscore = negativevoxel(
                cbf=self.inputs.in_avgscore,
                gm=self.inputs.in_greyM,
                thresh=thresh,
            )
            negscrub = negativevoxel(
                cbf=self.inputs.in_scrub,
                gm=self.inputs.in_greyM,
                thresh=thresh,
            )
        else:
            print("no score inputs, setting to np.nan")
            scorecbf_qei = np.nan
            scrub_qei = np.nan
            negscore = np.nan
            negscrub = np.nan

        if self.inputs.in_basil:
            basilcbf_qei = cbf_qei(
                gm=self.inputs.in_greyM,
                wm=self.inputs.in_whiteM,
                csf=self.inputs.in_csf,
                img=self.inputs.in_basil,
                thresh=thresh,
            )
            pvcbf_qei = cbf_qei(
                gm=self.inputs.in_greyM,
                wm=self.inputs.in_whiteM,
                csf=self.inputs.in_csf,
                img=self.inputs.in_pvc,
                thresh=thresh,
            )
            negbasil = negativevoxel(
                cbf=self.inputs.in_basil,
                gm=self.inputs.in_greyM,
                thresh=thresh,
            )
            negpvc = negativevoxel(
                cbf=self.inputs.in_pvc,
                gm=self.inputs.in_greyM,
                thresh=thresh,
            )
        else:
            print("no basil inputs, setting to np.nan")
            basilcbf_qei = np.nan
            pvcbf_qei = np.nan
            negbasil = np.nan
            negpvc = np.nan

        gwratio = np.divide(meancbf[0], meancbf[1])
        negcbf = negativevoxel(
            cbf=self.inputs.in_meancbf,
            gm=self.inputs.in_greyM,
            thresh=thresh,
        )

        metrics_dict = {
            "FD": [fd],
            "rmsd": [rms],
            "coregDC": [regDC],
            "coregJC": [regJC],
            "coregCC": [regCC],
            "coregCOV": [regCov],
            "cbfQEI": [meancbf_qei],
            "scoreQEI": [scorecbf_qei],
            "scrubQEI": [scrub_qei],
            "basilQEI": [basilcbf_qei],
            "pvcQEI": [pvcbf_qei],
            "GMmeanCBF": [meancbf[0]],
            "WMmeanCBF": [meancbf[1]],
            "Gm_Wm_CBF_ratio": [gwratio],
            "NEG_CBF_PERC": [negcbf],
            "NEG_SCORE_PERC": [negscore],
            "NEG_SCRUB_PERC": [negscrub],
            "NEG_BASIL_PERC": [negbasil],
            "NEG_PVC_PERC": [negpvc],
        }

        normalization_metrics_dict = {}
        if self.inputs.in_aslmaskstd and self.inputs.in_templatemask:
            normalization_metrics_dict = {
                "normDC": [normDC],
                "normJC": [normJC],
                "normCC": [normCC],
                "normCOV": [normCov],
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
            self.inputs.in_meancbf,
            suffix="qc_cbf.csv",
            newpath=runtime.cwd,
            use_ext=False,
        )
        qc_df.to_csv(self._results["qc_file"], index=False, header=True, na_rep="n/a")

        return runtime
