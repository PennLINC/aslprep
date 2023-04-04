"""Interfaces for calculating QC measures."""
import os

import numpy as np
import pandas as pd
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
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
    in_file = File(exists=True, mandatory=True, desc="original asl_file")
    in_meancbf = File(exists=True, mandatory=True, desc="cbf img")
    in_avgscore = File(exists=True, mandatory=False, desc="cbf img")
    in_scrub = File(exists=True, mandatory=False, desc="cbf img")
    in_basil = File(exists=True, mandatory=False, desc="cbf img")
    in_pvc = File(exists=True, mandatory=False, desc="cbf img")
    in_greyM = File(exists=True, mandatory=True, desc="grey  matter")
    in_whiteM = File(exists=True, mandatory=True, desc="white  matter")
    in_csf = File(exists=True, mandatory=True, desc="csf")
    in_confmat = File(exists=True, mandatory=False, desc=" cofnound matrix")
    in_aslmask = File(exists=True, mandatory=True, desc="asl mask in native space")
    in_t1mask = File(exists=True, mandatory=True, desc="t1wmask in native space ")
    in_aslmaskstd = File(exists=True, mandatory=False, desc="asl mask in native space")
    in_templatemask = File(exists=True, mandatory=False, desc="template mask or image")
    rmsd_file = File(exists=True, mandatory=True, desc="rmsd file")


class _ComputeCBFQCOutputSpec(TraitedSpec):
    qc_file = File(exists=False, desc="qc file ")


class ComputeCBFQC(SimpleInterface):
    """Calculate a series of CBF quality control metrics for non-GE data.

    compute qc from confound regressors
    and cbf maps,
    coregistration and regsitration indexes
    """

    input_spec = _ComputeCBFQCInputSpec
    output_spec = _ComputeCBFQCOutputSpec

    def _run_interface(self, runtime):
        time1 = pd.read_table(self.inputs.in_confmat)
        time1.fillna(0, inplace=True)
        fd = np.mean(time1["framewise_displacement"])
        # rms = time1[['rot_x', 'rot_y', 'rot_z']]
        # rms1 = rms.pow(2)
        # rms = np.mean(np.sqrt(rms1.sum(axis=1)/3))
        print(self.inputs.rmsd_file)
        rms = pd.read_csv(self.inputs.rmsd_file, header=None).mean().values[0]
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
            thresh=0.7,
        )
        meancbf = globalcbf(
            gm=self.inputs.in_greyM,
            wm=self.inputs.in_whiteM,
            csf=self.inputs.in_csf,
            cbf=self.inputs.in_meancbf,
            thresh=0.7,
        )

        if self.inputs.in_avgscore:
            scorecbf_qei = cbf_qei(
                gm=self.inputs.in_greyM,
                wm=self.inputs.in_whiteM,
                csf=self.inputs.in_csf,
                img=self.inputs.in_avgscore,
                thresh=0.7,
            )
            scrub_qei = cbf_qei(
                gm=self.inputs.in_greyM,
                wm=self.inputs.in_whiteM,
                csf=self.inputs.in_csf,
                img=self.inputs.in_scrub,
                thresh=0.7,
            )
            negscore = negativevoxel(
                cbf=self.inputs.in_avgscore, gm=self.inputs.in_greyM, thresh=0.7
            )
            negscrub = negativevoxel(cbf=self.inputs.in_scrub, gm=self.inputs.in_greyM, thresh=0.7)
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
                thresh=0.7,
            )
            pvcbf_qei = cbf_qei(
                gm=self.inputs.in_greyM,
                wm=self.inputs.in_whiteM,
                csf=self.inputs.in_csf,
                img=self.inputs.in_pvc,
                thresh=0.7,
            )
            negbasil = negativevoxel(cbf=self.inputs.in_basil, gm=self.inputs.in_greyM, thresh=0.7)
            negpvc = negativevoxel(cbf=self.inputs.in_pvc, gm=self.inputs.in_greyM, thresh=0.7)
        else:
            print("no basil inputs, setting to np.nan")
            basilcbf_qei = np.nan
            pvcbf_qei = np.nan
            negbasil = np.nan
            negpvc = np.nan

        gwratio = np.divide(meancbf[0], meancbf[1])
        negcbf = negativevoxel(cbf=self.inputs.in_meancbf, gm=self.inputs.in_greyM, thresh=0.7)

        if self.inputs.in_aslmaskstd and self.inputs.in_templatemask:
            dict1 = {
                "FD": [fd],
                "rmsd": [rms],
                "coregDC": [regDC],
                "coregJC": [regJC],
                "coregCC": [regCC],
                "coregCOV": [regCov],
                "normDC": [normDC],
                "normJC": [normJC],
                "normCC": [normCC],
                "normCOV": [normCov],
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
        else:
            dict1 = {
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
        _, file1 = os.path.split(self.inputs.in_file)
        bb = file1.split("_")
        dict2 = {}
        for i in range(len(bb) - 1):
            dict2.update({bb[i].split("-")[0]: bb[i].split("-")[1]})
        dict2.update(dict1)

        df = pd.DataFrame(dict2)

        self._results["qc_file"] = fname_presuffix(
            self.inputs.in_meancbf,
            suffix="qc_cbf.csv",
            newpath=runtime.cwd,
            use_ext=False,
        )
        df.to_csv(self._results["qc_file"], index=False, header=True)

        return runtime


class _ComputeCBFQCforGEInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="original asl_file")
    in_meancbf = File(exists=True, mandatory=True, desc="cbf img")
    in_avgscore = File(exists=True, mandatory=False, desc="cbf img")
    in_scrub = File(exists=True, mandatory=False, desc="cbf img")
    in_basil = File(exists=True, mandatory=False, desc="cbf img")
    in_pvc = File(exists=True, mandatory=False, desc="cbf img")
    in_greyM = File(exists=True, mandatory=True, desc="grey  matter")
    in_whiteM = File(exists=True, mandatory=True, desc="white  matter")
    in_csf = File(exists=True, mandatory=True, desc="csf")
    in_aslmask = File(exists=True, mandatory=True, desc="asl mask in native space")
    in_t1mask = File(exists=True, mandatory=True, desc="t1wmask in native space ")
    in_aslmaskstd = File(exists=True, mandatory=False, desc="asl mask in native space")
    in_templatemask = File(exists=True, mandatory=False, desc="template mask or image")


class _ComputeCBFQCforGEOutputSpec(TraitedSpec):
    qc_file = File(exists=False, desc="qc file ")


class ComputeCBFQCforGE(SimpleInterface):
    """Calculate a series of CBF quality control metrics for GE data.

    compute qc from confound regressors
    and cbf maps,
    coregistration and regsitration indexes
    """

    input_spec = _ComputeCBFQCforGEInputSpec
    output_spec = _ComputeCBFQCforGEOutputSpec

    def _run_interface(self, runtime):
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
            thresh=0.8,
        )
        meancbf = globalcbf(
            gm=self.inputs.in_greyM,
            wm=self.inputs.in_whiteM,
            csf=self.inputs.in_csf,
            cbf=self.inputs.in_meancbf,
            thresh=0.8,
        )

        if self.inputs.in_avgscore:
            scorecbf_qei = cbf_qei(
                gm=self.inputs.in_greyM,
                wm=self.inputs.in_whiteM,
                csf=self.inputs.in_csf,
                img=self.inputs.in_avgscore,
                thresh=0.8,
            )
            scrub_qei = cbf_qei(
                gm=self.inputs.in_greyM,
                wm=self.inputs.in_whiteM,
                csf=self.inputs.in_csf,
                img=self.inputs.in_scrub,
                thresh=0.8,
            )
            negscore = negativevoxel(
                cbf=self.inputs.in_avgscore, gm=self.inputs.in_greyM, thresh=0.8
            )
            negscrub = negativevoxel(cbf=self.inputs.in_scrub, gm=self.inputs.in_greyM, thresh=0.8)
        else:
            scorecbf_qei = 0
            scrub_qei = 0
            negscore = 0
            negscrub = 0

        if self.inputs.in_basil:
            basilcbf_qei = cbf_qei(
                gm=self.inputs.in_greyM,
                wm=self.inputs.in_whiteM,
                csf=self.inputs.in_csf,
                img=self.inputs.in_basil,
                thresh=0.8,
            )
            pvcbf_qei = cbf_qei(
                gm=self.inputs.in_greyM,
                wm=self.inputs.in_whiteM,
                csf=self.inputs.in_csf,
                img=self.inputs.in_pvc,
                thresh=0.8,
            )
            negbasil = negativevoxel(cbf=self.inputs.in_basil, gm=self.inputs.in_greyM, thresh=0.8)
            negpvc = negativevoxel(cbf=self.inputs.in_pvc, gm=self.inputs.in_greyM, thresh=0.8)
        else:
            basilcbf_qei = 0
            pvcbf_qei = 0
            negbasil = 0
            negpvc = 0
        gwratio = np.divide(meancbf[0], meancbf[1])
        negcbf = negativevoxel(cbf=self.inputs.in_meancbf, gm=self.inputs.in_greyM, thresh=0.8)

        if self.inputs.in_aslmaskstd and self.inputs.in_templatemask:
            dict1 = {
                "FD": 0,
                "relRMS": 0,
                "rmsd": np.nan,
                "coregDC": [regDC],
                "coregJC": [regJC],
                "coregCC": [regCC],
                "coregCOV": [regCov],
                "normDC": [normDC],
                "normJC": [normJC],
                "normCC": [normCC],
                "normCOV": [normCov],
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
        else:
            dict1 = {
                "FD": 0,
                "relRMS": 0,
                "rmsd": np.nan,
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
        _, file1 = os.path.split(self.inputs.in_file)
        bb = file1.split("_")
        dict2 = {}
        for i in range(len(bb) - 1):
            dict2.update({bb[i].split("-")[0]: bb[i].split("-")[1]})
        dict2.update(dict1)

        df = pd.DataFrame(dict2)

        self._results["qc_file"] = fname_presuffix(
            self.inputs.in_meancbf,
            suffix="qc_cbf.csv",
            newpath=runtime.cwd,
            use_ext=False,
        )
        df.to_csv(self._results["qc_file"], index=False, header=True)

        return runtime
