"""Interfaces for calculating CBF."""
import os

import nibabel as nb
import numpy as np
import pandas as pd
from nibabel.processing import smooth_image
from nilearn import image, maskers
from nipype import logging
from nipype.interfaces.ants import ApplyTransforms
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    isdefined,
    traits,
)
from nipype.interfaces.fsl import MultiImageMaths
from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec
from nipype.utils.filemanip import fname_presuffix

from aslprep.utils.misc import (
    _getcbfscore,
    _scrubcbf,
    parcellate_cbf,
    pcasl_or_pasl,
    readjson,
)
from aslprep.utils.qc import (
    cbf_qei,
    coverage,
    crosscorr,
    dice,
    globalcbf,
    jaccard,
    negativevoxel,
)

LOGGER = logging.getLogger("nipype.interface")


class _RefineMaskInputSpec(BaseInterfaceInputSpec):
    in_t1mask = File(exists=True, mandatory=True, desc="t1 mask")
    in_aslmask = File(exists=True, mandatory=True, desct="asl mask")
    transforms = File(exists=True, mandatory=True, desc="transfom")


class _RefineMaskOutputSpec(TraitedSpec):
    out_mask = File(exists=False, desc="output mask")
    out_tmp = File(exists=False, desc="tmp mask")


class RefineMask(SimpleInterface):
    """Reduce the ASL-derived brain mask using the associated T1w mask."""

    input_spec = _RefineMaskInputSpec
    output_spec = _RefineMaskOutputSpec

    def _run_interface(self, runtime):
        self._results["out_tmp"] = fname_presuffix(
            self.inputs.in_aslmask,
            suffix="_tempmask",
            newpath=runtime.cwd,
        )
        self._results["out_mask"] = fname_presuffix(
            self.inputs.in_aslmask,
            suffix="_refinemask",
            newpath=runtime.cwd,
        )

        refine_ref_mask(
            t1w_mask=self.inputs.in_t1mask,
            ref_asl_mask=self.inputs.in_aslmask,
            t12ref_transform=self.inputs.transforms,
            tmp_mask=self._results["out_tmp"],
            refined_mask=self._results["out_mask"],
        )

        return runtime


class _ExtractCBFInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="raw asl file")
    asl_file = File(exists=True, mandatory=True, desc="preprocessed asl file")
    in_mask = File(exists=True, mandatory=True, desc="mask")
    dummy_vols = traits.Int(
        default_value=0, exit=False, mandatory=False, desc="remove first n volumes"
    )
    in_metadata = traits.Dict(exists=True, mandatory=True, desc="metadata for asl or deltam ")
    bids_dir = traits.Str(exits=True, mandatory=True, desc=" bids directory")
    fwhm = traits.Float(default_value=5, exists=True, mandatory=False, desc="fwhm")


class _ExtractCBFOutputSpec(TraitedSpec):
    out_file = File(exists=False, desc="cbf time series data")
    out_avg = File(exists=False, desc="average control")


class ExtractCBF(SimpleInterface):
    """Extract CBF time series by subtracting label volumes from control volumes.

    TODO: Mock up test data and write tests to cover all of the branches in this interface.
    """

    input_spec = _ExtractCBFInputSpec
    output_spec = _ExtractCBFOutputSpec

    def _run_interface(self, runtime):
        file1 = os.path.abspath(self.inputs.in_file)
        # check if there is m0 file
        # m0num = 0
        m0file = []
        aslfile_linkedM0 = []
        mask = nb.load(self.inputs.in_mask).get_fdata()
        aslcontext1 = file1.replace("_asl.nii.gz", "_aslcontext.tsv")
        idasl = pd.read_csv(aslcontext1)["volume_type"].tolist()

        # read the data
        allasl = nb.load(self.inputs.asl_file)
        dataasl = allasl.get_fdata()

        # get the control,tag,moscan or label
        controllist = [i for i in range(0, len(idasl)) if idasl[i] == "control"]
        labellist = [i for i in range(0, len(idasl)) if idasl[i] == "label"]
        m0list = [i for i in range(0, len(idasl)) if idasl[i] == "m0scan"]
        deltamlist = [i for i in range(0, len(idasl)) if idasl[i] == "deltam"]
        cbflist = [i for i in range(0, len(idasl)) if idasl[i] == "CBF"]

        # extcract m0 file and register it to ASL if separate
        if self.inputs.in_metadata["M0Type"] == "Separate":
            m0file = self.inputs.in_file.replace("asl.nii.gz", "m0scan.nii.gz")
            m0file_metadata = readjson(m0file.replace("nii.gz", "json"))
            aslfile_linkedM0 = os.path.abspath(
                self.inputs.bids_dir + "/" + m0file_metadata["IntendedFor"]
            )
            if self.inputs.in_file not in aslfile_linkedM0:
                raise RuntimeError("there is no separate m0scan for the asl data")

            newm0 = fname_presuffix(self.inputs.asl_file, suffix="_m0file")
            newm0 = regmotoasl(asl=self.inputs.asl_file, m0file=m0file, m02asl=newm0)
            m0data_smooth = smooth_image(nb.load(newm0), fwhm=self.inputs.fwhm).get_fdata()
            if len(m0data_smooth.shape) > 3:
                m0dataf = mask * np.mean(m0data_smooth, axis=3)
            else:
                m0dataf = mask * m0data_smooth

        elif self.inputs.in_metadata["M0Type"] == "Included":
            modata2 = dataasl[:, :, :, m0list]
            con2 = nb.Nifti1Image(modata2, allasl.affine, allasl.header)
            m0data_smooth = smooth_image(con2, fwhm=self.inputs.fwhm).get_fdata()
            if len(m0data_smooth.shape) > 3:
                m0dataf = mask * np.mean(m0data_smooth, axis=3)
            else:
                m0dataf = mask * m0data_smooth

        elif self.inputs.in_metadata["M0Type"] == "Estimate":
            moestimate = self.inputs.in_metadata["M0Estimate"]
            m0dataf = moestimate * mask

        elif self.inputs.in_metadata["M0Type"] == "Absent":
            if len(controllist) > 0:
                control_img = dataasl[:, :, :, controllist]
                con = nb.Nifti1Image(control_img, allasl.affine, allasl.header)
                control_img1 = smooth_image(con, fwhm=self.inputs.fwhm).get_fdata()
                m0dataf = mask * np.mean(control_img1, axis=3)
            elif len(cbflist) > 0:
                m0dataf = mask
            else:
                raise RuntimeError("m0scan is absent")

        else:
            raise RuntimeError("no pathway to m0scan")

        if len(dataasl.shape) == 5:
            raise RuntimeError("Input image (%s) is 5D.", self.inputs.asl_file)

        if len(deltamlist) > 0:
            cbf_data = dataasl[:, :, :, deltamlist]

        if len(cbflist) > 0:
            cbf_data = dataasl[:, :, :, cbflist]
        elif len(labellist) > 0:
            control_img = dataasl[:, :, :, controllist]
            label_img = dataasl[:, :, :, labellist]
            cbf_data = np.subtract(control_img, label_img)
        else:
            raise RuntimeError("no valid asl or cbf image.")

        if self.inputs.dummy_vols != 0:
            cbf_data = np.delete(cbf_data, range(0, self.inputs.dummy_vols), axis=3)
            # control_img = np.delete(control_img, range(0, self.inputs.dummy_vols), axis=3)

        self._results["out_file"] = fname_presuffix(
            self.inputs.in_file,
            suffix="_cbftimeseries",
            newpath=runtime.cwd,
        )
        self._results["out_avg"] = fname_presuffix(
            self.inputs.in_file,
            suffix="_m0file",
            newpath=runtime.cwd,
        )
        nb.Nifti1Image(cbf_data, allasl.affine, allasl.header).to_filename(
            self._results["out_file"]
        )
        nb.Nifti1Image(m0dataf, allasl.affine, allasl.header).to_filename(self._results["out_avg"])

        return runtime


class _ComputeCBFInputSpec(BaseInterfaceInputSpec):
    in_cbf = File(exists=True, mandatory=True, desc="cbf nifti")
    in_metadata = traits.Dict(exists=True, mandatory=True, desc="metadata for CBF ")
    in_m0scale = traits.Float(
        exists=True,
        mandatory=True,
        desc="relative scale between asl and m0",
    )
    in_m0file = File(exists=True, mandatory=False, desc="M0 nifti file")
    in_mask = File(exists=True, mandatory=False, desc="mask")


class _ComputeCBFOutputSpec(TraitedSpec):
    out_cbf = File(exists=False, desc="cbf timeries data")
    out_mean = File(exists=False, desc="average control")


class ComputeCBF(SimpleInterface):
    """Calculate CBF time series, mean control, and arterial transit time."""

    input_spec = _ComputeCBFInputSpec
    output_spec = _ComputeCBFOutputSpec

    def _run_interface(self, runtime):
        metadata = self.inputs.in_metadata
        m0file = self.inputs.in_m0file
        m0scale = self.inputs.in_m0scale
        mask = self.inputs.in_mask
        cl_diff_file = self.inputs.in_cbf  # control - label signal intensities

        is_casl = pcasl_or_pasl(metadata=metadata)
        plds = np.array(metadata["PostLabelingDelay"])

        # https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.24550
        t1blood = (110 * int(metadata["MagneticFieldStrength"]) + 1316) / 1000

        # Get labeling efficiency (alpha in Alsop 2015)
        if "LabelingEfficiency" in metadata.keys():
            labeleff = metadata["LabelingEfficiency"]
        elif is_casl:
            labeleff = 0.72
        else:
            labeleff = 0.8

        PARTITION_COEF = 0.9  # brain partition coefficient (lambda in Alsop 2015)

        if is_casl:
            tau = metadata["LabelingDuration"]
            perfusion_factor = np.exp(plds / t1blood) / (t1blood * (1 - np.exp(-(tau / t1blood))))
        else:
            inversiontime = plds  # As per BIDS: inversiontime for PASL == PostLabelingDelay
            inversiontime1 = metadata["BolusCutOffDelayTime"]  # called TI1 in Alsop 2015
            perfusion_factor = np.exp(inversiontime / t1blood) / inversiontime1

        perfusion_factor *= (6000 * PARTITION_COEF) / (2 * labeleff)

        # NOTE: Nilearn will still add a singleton time dimension for 3D imgs with
        # NiftiMasker.transform, until 0.12.0, so the arrays will currently be 2D no matter what.
        masker = maskers.NiftiMasker(mask_img=mask)
        cl_diff_arr = masker.fit_transform(cl_diff_file).T
        n_voxels, n_volumes = cl_diff_arr.shape
        m0data = masker.transform(m0file).T
        scaled_m0data = m0scale * m0data

        # Scale difference signal to absolute CBF units by dividing by PD image (M0)
        cl_diff_scaled = cl_diff_arr / scaled_m0data

        if (perfusion_factor.size > 1) and (n_volumes > 1):
            # CBF is a time series and there are multiple PLDs.
            permfactor = np.tile(
                perfusion_factor,
                int(n_volumes / len(perfusion_factor)),
            )

            cbf_data_ts = cl_diff_scaled * permfactor

            cbf = np.zeros([n_voxels, int(n_volumes / len(perfusion_factor))])
            cbf_xx = np.split(cbf_data_ts, int(n_volumes / len(perfusion_factor)), axis=1)

            # calculate weighted cbf with multiple PostLabelingDelays
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3791289/
            # https://pubmed.ncbi.nlm.nih.gov/22084006/
            for k in range(len(cbf_xx)):
                cbf_plds = cbf_xx[k]
                pldx = np.zeros(cbf_plds.shape)
                for j in range(cbf_plds.shape[0]):
                    pldx[:, j] = np.array(np.multiply(cbf_plds[:, j], plds[j]))

                cbf[:, k] = np.sum(pldx, axis=1) / np.sum(plds)

        elif (perfusion_factor.size > 1) and (n_volumes == 1):
            # CBF is not a time series, but there are multiple PLDs.
            cbf_ts = np.zeros(cl_diff_arr.shape, len(perfusion_factor))
            for i in len(perfusion_factor):
                cbf_ts[:, i] = cl_diff_scaled * perfusion_factor[i]

            cbf = np.sum(cbf_ts, axis=2) / np.sum(perfusion_factor)

        else:
            # There is only a single PLD.
            cbf = cl_diff_scaled * perfusion_factor

        # return cbf to nifti shape
        cbf = np.nan_to_num(cbf)
        tcbf = masker.inverse_transform(cbf.T)
        meancbf = image.mean_img(tcbf)

        self._results["out_cbf"] = fname_presuffix(
            self.inputs.in_cbf,
            suffix="_cbf",
            newpath=runtime.cwd,
        )
        self._results["out_mean"] = fname_presuffix(
            self.inputs.in_cbf,
            suffix="_meancbf",
            newpath=runtime.cwd,
        )
        tcbf.to_filename(self._results["out_cbf"])
        meancbf.to_filename(self._results["out_mean"])

        return runtime


class _ScoreAndScrubCBFInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="computed CBF from ComputeCBF")
    in_greyM = File(exists=True, mandatory=True, desc="grey  matter")
    in_whiteM = File(exists=True, mandatory=True, desc="white  matter")
    in_mask = File(exists=True, mandatory=True, desc="mask")
    in_csf = File(exists=True, mandatory=True, desc="csf")
    in_thresh = traits.Float(
        default_value=0.7, exists=True, mandatory=False, desc="threshold of propbaility matter"
    )
    in_wfun = traits.Str(
        exists=True,
        mandatory=False,
        default_value="huber",
        option=["bisquare", "andrews", "cauchy", "fair", "logistics", "ols", "talwar", "welsch"],
        desc="wavelet fun ",
    )
    out_score = File(exists=False, desc="score timeseries data")
    out_avgscore = File(exists=False, desc="average score")
    out_scrub = File(exists=False, desc="average scrub")
    out_scoreindex = File(exists=False, desc="index of volume remove or leave by score")


class _ScoreAndScrubCBFOutputSpec(TraitedSpec):
    out_score = File(exists=False, mandatory=False, desc="score timeseries data")
    out_avgscore = File(exists=False, mandatory=False, desc="average score")
    out_scrub = File(exists=False, mandatory=False, desc="average scrub")
    out_scoreindex = File(exists=False, mandatory=False, desc="index of volume remove ")


class ScoreAndScrubCBF(SimpleInterface):
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

    input_spec = _ScoreAndScrubCBFInputSpec
    output_spec = _ScoreAndScrubCBFOutputSpec

    def _run_interface(self, runtime):
        cbf_ts = nb.load(self.inputs.in_file).get_fdata()
        mask = nb.load(self.inputs.in_mask).get_fdata()
        greym = nb.load(self.inputs.in_greyM).get_fdata()
        whitem = nb.load(self.inputs.in_whiteM).get_fdata()
        csf = nb.load(self.inputs.in_csf).get_fdata()
        if cbf_ts.ndim > 3:
            cbf_scorets, index_score = _getcbfscore(
                cbfts=cbf_ts,
                wm=whitem,
                gm=greym,
                csf=csf,
                mask=mask,
                thresh=self.inputs.in_thresh,
            )
            cbfscrub = _scrubcbf(
                cbf_ts=cbf_scorets,
                gm=greym,
                wm=whitem,
                csf=csf,
                mask=mask,
                wfun=self.inputs.in_wfun,
                thresh=self.inputs.in_thresh,
            )
            avgscore = np.mean(cbf_scorets, axis=3)
        else:
            LOGGER.warning(f"CBF time series is only {cbf_ts.ndim}D. Skipping SCORE and SCRUB.")
            cbf_scorets = cbf_ts
            index_score = np.array([0])
            cbfscrub = cbf_ts
            avgscore = cbf_ts

        self._results["out_score"] = fname_presuffix(
            self.inputs.in_file,
            suffix="_cbfscorets",
            newpath=runtime.cwd,
        )
        self._results["out_avgscore"] = fname_presuffix(
            self.inputs.in_file,
            suffix="_meancbfscore",
            newpath=runtime.cwd,
        )
        self._results["out_scrub"] = fname_presuffix(
            self.inputs.in_file,
            suffix="_cbfscrub",
            newpath=runtime.cwd,
        )
        self._results["out_scoreindex"] = fname_presuffix(
            self.inputs.in_file,
            suffix="_scoreindex.txt",
            newpath=runtime.cwd,
            use_ext=False,
        )
        samplecbf = nb.load(self.inputs.in_mask)

        nb.Nifti1Image(
            dataobj=cbf_scorets,
            affine=samplecbf.affine,
            header=samplecbf.header,
        ).to_filename(self._results["out_score"])
        nb.Nifti1Image(
            dataobj=avgscore,
            affine=samplecbf.affine,
            header=samplecbf.header,
        ).to_filename(self._results["out_avgscore"])
        nb.Nifti1Image(
            dataobj=cbfscrub,
            affine=samplecbf.affine,
            header=samplecbf.header,
        ).to_filename(self._results["out_scrub"])

        np.savetxt(self._results["out_scoreindex"], index_score, delimiter=",")

        self.inputs.out_score = os.path.abspath(self._results["out_score"])
        self.inputs.out_avgscore = os.path.abspath(self._results["out_avgscore"])
        self.inputs.out_scrub = os.path.abspath(self._results["out_scrub"])
        self.inputs.out_scoreindex = os.path.abspath(self._results["out_scoreindex"])
        return runtime


class _BASILCBFInputSpec(FSLCommandInputSpec):
    # We use position args here as list indices - so a negative number
    # will put something on the end
    in_file = File(
        exists=True,
        desc=(
            "ASL data after subtracting tag-control or control-tag. "
            "This matches with ``--iaf diff``, which is the default."
        ),
        argstr="-i %s",
        position=0,
        mandatory=True,
    )
    mask = File(
        exists=True,
        argstr="-m %s",
        desc="mask in the same space as in_infile",
        mandatory=True,
    )
    mzero = File(exists=True, argstr="-c %s", desc="m0 scan", mandatory=False)
    m0scale = traits.Float(desc="calibration of asl", argstr="--cgain %.2f", mandatory=True)
    m0tr = traits.Float(
        desc="The repetition time for the calibration image (the M0 scan).",
        argstr="--tr %.2f",
        mandatory=True,
    )
    tis = traits.Either(
        traits.Float(),
        traits.List(traits.Float()),
        desc=(
            "The list of inflow times (TIs), a comma separated list of values should be provided "
            "(that matches the order in the data).\n\n"
            "Note, the inflow time is the PLD plus bolus duration for pcASL (and cASL), "
            "it equals the inversion time for pASL. "
            "If the data contains multiple repeats of the same set of TIs then it is only "
            "necessary to list the unique TIs.\n\n"
            "When using the ``--tis=`` you can specify a full list of all TIs/PLDs in the data "
            "(i.e., as many entries as there are label-control pairs). "
            "Or, if you have a number of TIs/PLDs repeated multiple times you can just list the "
            "unique TIs in order and ``oxford_asl`` will automatically replicate that list to "
            "match the number of repeated measurements in the data. "
            "If you have a variable number of repeats at each TI/PLD then either list all TIs "
            "or use the ``--rpts=<csv>`` option (see below)."
        ),
        argstr="--tis %s",
        mandatory=True,
        sep=",",
    )
    pcasl = traits.Bool(
        desc=(
            "Data were acquired using cASL or pcASL labelling "
            "(pASL labeling is assumed by default)."
        ),
        argstr="--casl",
        mandatory=False,
        default_value=False,
    )
    bolus = traits.Either(
        traits.Float(),
        traits.List(traits.Float()),
        desc="bolus or tau: label duration",
        argstr="--bolus %s",
        mandatory=True,
        sep=",",
    )
    pvc = traits.Bool(
        desc="Do partial volume correction.",
        mandatory=False,
        argstr="--pvcorr",
        default_value=True,
    )
    pvgm = File(
        exists=True,
        mandatory=False,
        desc="Partial volume estimates for GM. This is just a GM tissue probability map.",
        argstr="--pvgm %s",
    )
    pvwm = File(
        exists=True,
        mandatory=False,
        desc="Partial volume estimates for WM. This is just a WM tissue probability map.",
        argstr="--pvwm %s",
    )
    alpha = traits.Float(
        desc=(
            "Inversion efficiency - [default: 0.98 (pASL); 0.85 (cASL)]. "
            "This is equivalent the BIDS metadata field 'LabelingEfficiency'."
        ),
        argstr="--alpha %.2f",
    )
    out_basename = File(desc="base name of output files", argstr="-o %s", mandatory=True)


class _BASILCBFOutputSpec(TraitedSpec):
    out_cbfb = File(exists=False, desc="cbf with spatial correction")
    out_cbfpv = File(exists=False, desc="cbf with spatial correction")
    out_cbfpvwm = File(
        exists=False, desc="cbf with spatial partial volume white matter correction"
    )
    out_att = File(exists=False, desc="aretrial transist time")


class BASILCBF(FSLCommand):
    """Apply Bayesian Inference for Arterial Spin Labeling (BASIL).

    This interface calculates:
    (1) arterial transit time,
    (2) CBF with spatial correction,
    (3) CBF with spatial partial volume white matter correction, and
    (4) CBF with spatial partial volume correction.

    See https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BASIL and https://asl-docs.readthedocs.io.
    """

    _cmd = "oxford_asl"
    input_spec = _BASILCBFInputSpec
    output_spec = _BASILCBFOutputSpec

    def _run_interface(self, runtime):
        runtime = super(BASILCBF, self)._run_interface(runtime)
        return runtime

    def _gen_outfilename(self, suffix):
        if isdefined(self.inputs.in_file):
            out_file = self._gen_fname(self.inputs.in_file, suffix=suffix)
        return os.path.abspath(out_file)

    def _list_outputs(self):
        basename = self.inputs.out_basename

        outputs = self.output_spec().get()

        outputs["out_cbfb"] = os.path.join(basename, "native_space/perfusion_calib.nii.gz")
        outputs["out_att"] = os.path.join(basename, "native_space/arrival.nii.gz")
        outputs["out_cbfpv"] = os.path.join(basename, "native_space/pvcorr/perfusion_calib.nii.gz")
        outputs["out_cbfpvwm"] = os.path.join(
            basename,
            "native_space/pvcorr/perfusion_wm_calib.nii.gz",
        )

        return outputs


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
        time1 = pd.read_csv(self.inputs.in_confmat, sep="\t")
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


class _ParcellateCBFInputSpec(BaseInterfaceInputSpec):
    in_cbf = File(exists=True, mandatory=True, desc="cbf img")
    atlasfile = File(exists=True, mandatory=True, desc="data")
    atlasdata = File(exists=True, mandatory=True, desc="data")
    atlaslabel = File(exists=True, mandatory=True, desc="data")


class _ParcellateCBFOutputSpec(TraitedSpec):
    atlascsv = File(exists=False, desc="harvard output csv")


class ParcellateCBF(SimpleInterface):
    """Parcellate CBF time series according to a given atlas."""

    input_spec = _ParcellateCBFInputSpec
    output_spec = _ParcellateCBFOutputSpec

    def _run_interface(self, runtime):
        self._results["atlascsv"] = fname_presuffix(
            self.inputs.in_cbf,
            suffix="atlas.csv",
            newpath=runtime.cwd,
            use_ext=False,
        )
        roiquant = parcellate_cbf(
            roi_label=self.inputs.atlaslabel,
            roi_file=self.inputs.atlasfile,
            cbfmap=self.inputs.in_cbf,
        )
        data1 = pd.read_table(self.inputs.atlasdata, header=None, index_col=None, sep="\t")
        bb = list(data1.values.tolist())
        flattened = [val for sublist in bb for val in sublist]
        datat = pd.DataFrame([flattened, roiquant])
        datat.to_csv(self._results["atlascsv"], header=None, index=None)
        return runtime


class _ExtractCBForDeltaMInputSpec(BaseInterfaceInputSpec):
    in_asl = File(exists=True, mandatory=True, desc="raw asl file")
    in_aslmask = File(exists=True, mandatory=True, desct="asl mask")
    file_type = traits.Str(desc="file type, c for cbf, d for deltam", mandatory=True)


class _ExtractCBForDeltaMOutputSpec(TraitedSpec):
    out_file = File(exists=False, desc="cbf or deltam")


class ExtractCBForDeltaM(SimpleInterface):
    """Load an ASL file and grab the CBF or DeltaM volumes from it."""

    input_spec = _ExtractCBForDeltaMInputSpec
    output_spec = _ExtractCBForDeltaMOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = fname_presuffix(
            self.inputs.in_aslmask,
            suffix="_cbfdeltam",
            newpath=runtime.cwd,
        )
        filex = self.inputs.in_asl
        # NOTE: Not a good way to find the aslcontext file.
        aslcontext = pd.read_csv(filex.replace("_asl.nii.gz", "_aslcontext.tsv"))
        idasl = aslcontext["volume_type"].tolist()
        fdata = nb.load(filex).get_fdata()
        img = nb.load(filex)

        controllist = [i for i in range(0, len(idasl)) if idasl[i] == "control"]
        labelist = [i for i in range(0, len(idasl)) if idasl[i] == "label"]

        if self.inputs.file_type == "d":
            if len(controllist) > 0:
                # Grab control and label volumes from ASL file,
                # then calculate deltaM by subtracting label volumes from control volumes.
                ffdata = fdata[:, :, :, controllist] - fdata[:, :, :, labelist]
                newdata = nb.Nifti1Image(dataobj=ffdata, affine=img.affine, header=img.header)
            else:
                # Grab deltaM volumes from ASL file.
                dlist = [i for i in range(0, len(idasl)) if idasl[i] == "deltam"]
                if len(fdata.shape) < 4:
                    newdata = nb.Nifti1Image(dataobj=fdata, affine=img.affine, header=img.header)
                else:
                    ffdata = fdata[:, :, :, dlist]
                    newdata = nb.Nifti1Image(dataobj=ffdata, affine=img.affine, header=img.header)

        elif self.inputs.file_type == "c":
            dlist = [i for i in range(0, len(idasl)) if idasl[i] == "CBF"]
            if len(fdata.shape) < 4:
                # 3D volume is written out without any changes.
                # NOTE: Why not return the original file then?
                newdata = nb.Nifti1Image(dataobj=fdata, affine=img.affine, header=img.header)
            else:
                # Grab CBF volumes from ASL file.
                ffdata = fdata[:, :, :, dlist]
                newdata = nb.Nifti1Image(dataobj=ffdata, affine=img.affine, header=img.header)

        newdata.to_filename(self._results["out_file"])

        return runtime


def regmotoasl(asl, m0file, m02asl):
    """Calculate mean M0 image and mean ASL image, then FLIRT M0 image to ASL space.

    TODO: This should not be a function. It uses interfaces, so it should be a workflow.
    """
    from nipype.interfaces import fsl

    meanasl = fsl.MeanImage()
    meanasl.inputs.in_file = asl
    meanasl.inputs.out_file = fname_presuffix(asl, suffix="_meanasl")
    meanasl.run()
    meanm0 = fsl.MeanImage()
    meanm0.inputs.in_file = m0file
    meanm0.inputs.out_file = fname_presuffix(asl, suffix="_meanm0")
    meanm0.run()
    flt = fsl.FLIRT(bins=640, cost_func="mutualinfo")
    flt.inputs.in_file = meanm0.inputs.out_file
    flt.inputs.reference = meanasl.inputs.out_file
    flt.inputs.out_file = m02asl
    flt.run()
    return m02asl


def refine_ref_mask(t1w_mask, ref_asl_mask, t12ref_transform, tmp_mask, refined_mask):
    """Warp T1w mask to ASL space, then use it to mask the ASL mask.

    TODO: This should not be a function. It uses interfaces, so it should be a workflow.
    """
    b1 = ApplyTransforms()
    b1.inputs.dimension = 3
    b1.inputs.float = True
    b1.inputs.input_image = t1w_mask
    b1.inputs.interpolation = "NearestNeighbor"
    b1.inputs.reference_image = ref_asl_mask
    b1.inputs.transforms = t12ref_transform
    b1.inputs.input_image_type = 3
    b1.inputs.output_image = tmp_mask
    b1.run()

    mat1 = MultiImageMaths()
    mat1.inputs.in_file = tmp_mask
    mat1.inputs.op_string = " -mul  %s -bin"
    mat1.inputs.operand_files = ref_asl_mask
    mat1.inputs.out_file = refined_mask
    mat1.run()

    return refined_mask
