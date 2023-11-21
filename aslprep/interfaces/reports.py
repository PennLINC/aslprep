# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Interfaces to generate reportlets."""

import os
import re
import time

import pandas as pd
from fmriprep.interfaces.reports import get_world_pedir
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    Directory,
    File,
    InputMultiObject,
    SimpleInterface,
    Str,
    TraitedSpec,
    isdefined,
    traits,
)
from smriprep.interfaces.freesurfer import ReconAll

SUBJECT_TEMPLATE = """\
\t<ul class="elem-desc">
\t\t<li>Subject ID: {subject_id}</li>
\t\t<li>Structural images: {n_t1s:d} T1-weighted {t2w}</li>
\t\t<li>ASL series: {n_asl:d}</li>
\t\t<li>Standard output spaces: {std_spaces}</li>
\t\t<li>Non-standard output spaces: {nstd_spaces}</li>
\t\t<li>FreeSurfer reconstruction: {freesurfer_status}</li>
\t</ul>
"""

FUNCTIONAL_TEMPLATE = """\t\t<h3 class="elem-title">Summary</h3>
\t\t<ul class="elem-desc">
\t\t\t<li>Repetition time (TR): {tr:.03g}s</li>
\t\t\t<li>Phase-encoding (PE) direction: {pedir}</li>
\t\t\t<li>Susceptibility distortion correction: {sdc}</li>
\t\t\t<li>Registration: {registration}</li>

\t\t</ul>
"""

QC_TEMPLATE = """\t\t<h3 class="elem-title">QC Summary</h3>
\t\t<ul class="elem-desc">
\t\t\t<li>Confounds collected: {confounds}</li>
\t\t\t<li>Motion summary measures: {motionparam}</li>
\t\t\t<li>Coregistration quality: {coregindex}</li>
\t\t\t<li>Normalization quality: {normindex}</li>
\t\t\t<li>Quality evaluation index : {qei}</li>
\t\t\t<li>Mean CBF (mL 100/g/min) : {mean_cbf}</li>
\t\t\t<li>Percentage of negative voxel : {negvoxel}</li>

\t\t</ul>
"""

ABOUT_TEMPLATE = """\t<ul>
\t\t<li>ASLPrep version: {version}</li>
\t\t<li>ASLPrep command: <code>{command}</code></li>
\t\t<li>Date preprocessed: {date}</li>
\t</ul>
</div>
"""


class _SummaryOutputSpec(TraitedSpec):
    out_report = File(exists=True, desc="HTML segment containing summary")


class SummaryInterface(SimpleInterface):
    """A basic summary interface."""

    output_spec = _SummaryOutputSpec

    def _run_interface(self, runtime):
        segment = self._generate_segment()
        fname = os.path.join(runtime.cwd, "report.html")
        with open(fname, "w") as fobj:
            fobj.write(segment)
        self._results["out_report"] = fname
        return runtime

    def _generate_segment(self):
        raise NotImplementedError


class _SubjectSummaryInputSpec(BaseInterfaceInputSpec):
    t1w = InputMultiObject(File(exists=True), desc="T1w structural images")
    t2w = InputMultiObject(File(exists=True), desc="T2w structural images")
    subjects_dir = Directory(desc="FreeSurfer subjects directory")
    subject_id = Str(desc="Subject ID")
    asl = InputMultiObject(
        traits.Either(File(exists=True), traits.List(File(exists=True))),
        desc="ASL functional series",
    )
    std_spaces = traits.List(Str, desc="list of standard spaces")
    nstd_spaces = traits.List(Str, desc="list of non-standard spaces")


class _SubjectSummaryOutputSpec(_SummaryOutputSpec):
    # This exists to ensure that the summary is run prior to the first ReconAll
    # call, allowing a determination whether there is a pre-existing directory
    subject_id = Str(desc="FreeSurfer subject ID")


class SubjectSummary(SummaryInterface):
    """A summary describing the subject's data as a whole."""

    input_spec = _SubjectSummaryInputSpec
    output_spec = _SubjectSummaryOutputSpec

    def _run_interface(self, runtime):
        if isdefined(self.inputs.subject_id):
            self._results["subject_id"] = self.inputs.subject_id
        return super(SubjectSummary, self)._run_interface(runtime)

    def _generate_segment(self):
        if not isdefined(self.inputs.subjects_dir):
            freesurfer_status = "Not run"
        else:
            recon = ReconAll(
                subjects_dir=self.inputs.subjects_dir,
                subject_id="sub-" + self.inputs.subject_id,
                T1_files=self.inputs.t1w,
                flags="-noskullstrip",
            )
            if recon.cmdline.startswith("echo"):
                freesurfer_status = "Pre-existing directory"
            else:
                freesurfer_status = "Run by ASLPrep"

        t2w_seg = ""
        if self.inputs.t2w:
            t2w_seg = f"(+ {len(self.inputs.t2w):d} T2-weighted)"

        # Add list of tasks with number of runs
        asl_series = self.inputs.asl if isdefined(self.inputs.asl) else []
        asl_series = [s[0] if isinstance(s, list) else s for s in asl_series]

        return SUBJECT_TEMPLATE.format(
            subject_id=self.inputs.subject_id,
            n_t1s=len(self.inputs.t1w),
            t2w=t2w_seg,
            n_asl=len(asl_series),
            std_spaces=", ".join(self.inputs.std_spaces),
            nstd_spaces=", ".join(self.inputs.nstd_spaces),
            freesurfer_status=freesurfer_status,
        )


class _FunctionalSummaryInputSpec(BaseInterfaceInputSpec):
    distortion_correction = traits.Str(
        desc="Susceptibility distortion correction method",
        mandatory=True,
    )
    pe_direction = traits.Enum(
        None,
        "i",
        "i-",
        "j",
        "j-",
        mandatory=True,
        desc="Phase-encoding direction detected",
    )
    registration = traits.Enum(
        "FSL",
        "FreeSurfer",
        mandatory=True,
        desc="Functional/anatomical registration method",
    )
    fallback = traits.Bool(desc="Boundary-based registration rejected")
    registration_dof = traits.Enum(
        6,
        9,
        12,
        desc="Registration degrees of freedom",
        mandatory=True,
    )
    registration_init = traits.Enum(
        "register",
        "header",
        mandatory=True,
        desc='Whether to initialize registration with the "header"'
        ' or by centering the volumes ("register")',
    )
    confounds_file = File(exists=True, mandatory=False, desc="Confounds file")
    qc_file = File(exists=True, desc="qc file")
    tr = traits.Float(desc="Repetition time", mandatory=True)
    orientation = traits.Str(mandatory=True, desc="Orientation of the voxel axes")


class FunctionalSummary(SummaryInterface):
    """A summary of a functional run, with QC measures included."""

    input_spec = _FunctionalSummaryInputSpec

    def _generate_segment(self):
        dof = self.inputs.registration_dof
        reg = {
            "FSL": [
                "FSL <code>flirt</code> with boundary-based registration"
                f" (BBR) metric - {dof} dof",
                "FSL <code>flirt</code> rigid registration - 6 dof",
            ],
            "FreeSurfer": [
                "FreeSurfer <code>bbregister</code> "
                f"(boundary-based registration, BBR) - {dof} dof",
                f"FreeSurfer <code>mri_coreg</code> - {dof} dof",
            ],
        }[self.inputs.registration][self.inputs.fallback]

        pedir = get_world_pedir(self.inputs.orientation, self.inputs.pe_direction)

        if self.inputs.pe_direction is None:
            pedir = "MISSING - Assuming Anterior-Posterior"
        else:
            pedir = {"i": "Left-Right", "j": "Anterior-Posterior"}[self.inputs.pe_direction[0]]

        # the number of dummy scans was specified by the user and
        # it is not equal to the number detected by the algorithm

        # the number of dummy scans was not specified by the user

        return FUNCTIONAL_TEMPLATE.format(
            pedir=pedir,
            sdc=self.inputs.distortion_correction,
            registration=reg,
            tr=self.inputs.tr,
        )


class _CBFSummaryInputSpec(BaseInterfaceInputSpec):
    confounds_file = File(exists=True, mandatory=False, desc="Confounds file")
    qc_file = File(exists=True, desc="qc file")


class CBFSummary(SummaryInterface):
    """A summary of a functional run, with QC measures included."""

    input_spec = _CBFSummaryInputSpec

    def _generate_segment(self):
        qcfile = pd.read_csv(self.inputs.qc_file)
        motionparam = f"FD : {round(qcfile['FD'][0], 4)}, rmsd: {round(qcfile['rmsd'][0], 4)} "
        coregindex = (
            f" Dice Index: {round(qcfile['coregDC'][0], 4)}, "
            f"Jaccard Index: {round(qcfile['coregJC'][0], 4)}, "
            f"Cross Cor.: {round(qcfile['coregCC'][0], 4)}, "
            f"Coverage: {round(qcfile['coregCOV'][0], 4)} "
        )
        normindex = (
            f" Dice Index: {round(qcfile['normDC'][0], 4)}, "
            f"Jaccard Index: {round(qcfile['normJC'][0], 4)}, "
            f"Cross Cor.: {round(qcfile['normCC'][0], 4)}, "
            f"Coverage: {round(qcfile['normCOV'][0], 4)} "
        )
        qei = (
            f"cbf: {round(qcfile['cbfQEI'][0], 4)}, "
            f"score: {round(qcfile['scoreQEI'][0], 4)}, "
            f"scrub: {round(qcfile['scrubQEI'][0], 4)}, "
            f"basil: {round(qcfile['basilQEI'][0], 4)}, "
            f"pvc: {round(qcfile['pvcQEI'][0], 4)} "
        )
        mean_cbf = (
            f"GM CBF: {round(qcfile['GMmeanCBF'][0], 2)}, "
            f"WM CBF: {round(qcfile['WMmeanCBF'][0], 2)}, "
            f"GM/WM CBF ratio: {round(qcfile['Gm_Wm_CBF_ratio'][0], 2)} "
        )
        negvoxel = (
            f"cbf: {round(qcfile['NEG_CBF_PERC'][0], 2)}, "
            f"score: {round(qcfile['NEG_SCORE_PERC'][0], 2)}, "
            f"scrub: {round(qcfile['NEG_SCRUB_PERC'][0], 2)}, "
            f"basil: {round(qcfile['NEG_BASIL_PERC'][0], 2)}, "
            f"pvc: {round(qcfile['NEG_PVC_PERC'][0], 2)} "
        )

        if isdefined(self.inputs.confounds_file):
            with open(self.inputs.confounds_file) as cfh:
                conflist = cfh.readline().strip("\n").strip()
        else:
            conflist = "None"
        # the number of dummy scans was specified by the user and
        # it is not equal to the number detected by the algorithm

        # the number of dummy scans was not specified by the user

        return QC_TEMPLATE.format(
            confounds=re.sub(r"[\t ]+", ", ", conflist),
            motionparam=motionparam,
            qei=qei,
            coregindex=coregindex,
            normindex=normindex,
            mean_cbf=mean_cbf,
            negvoxel=negvoxel,
        )


class _AboutSummaryInputSpec(BaseInterfaceInputSpec):
    version = Str(desc="ASLPREP version")
    command = Str(desc="ASLPREP command")
    # Date not included - update timestamp only if version or command changes


class AboutSummary(SummaryInterface):
    """A basic summary of the ASLPrep run."""

    input_spec = _AboutSummaryInputSpec

    def _generate_segment(self):
        return ABOUT_TEMPLATE.format(
            version=self.inputs.version,
            command=self.inputs.command,
            date=time.strftime("%Y-%m-%d %H:%M:%S %z"),
        )
