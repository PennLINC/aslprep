# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Interfaces to generate reportlets."""

import os
import time
import re

from collections import Counter
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec,
    File, Directory, InputMultiObject, Str, isdefined,
    SimpleInterface)
from ..smriprep.interfaces.freesurfer import ReconAll


SUBJECT_TEMPLATE = """\
\t<ul class="elem-desc">
\t\t<li>Subject ID: {subject_id}</li>
\t\t<li>Structural images: {n_t1s:d} T1-weighted {t2w}</li>
\t\t<li>ASL series: {n_asl:d}</li>
{tasks}
\t\t<li>Standard output spaces: {std_spaces}</li>
\t\t<li>Non-standard output spaces: {nstd_spaces}</li>
\t\t<li>FreeSurfer reconstruction: {freesurfer_status}</li>
\t</ul>
"""

FUNCTIONAL_TEMPLATE = """\t\t<h3 class="elem-title">Summary</h3>
\t\t<ul class="elem-desc">
\t\t\t<li>Repetition time (TR): {tr:.03g}s</li>
\t\t\t<li>Phase-encoding (PE) direction: {pedir}</li>
\t\t\t<li>Slice timing correction: {stc}</li>
\t\t\t<li>Susceptibility distortion correction: {sdc}</li>
\t\t\t<li>Registration: {registration}</li>
\t\t\t<li>Confounds collected: {confounds}</li>
\t\t\t<li>Motion summary measures: {motionparam}</li>
\t\t\t<li>Coregistration quality: {coregindex}</li>
\t\t\t<li>Normalization quality: {normindex}</li>
\t\t\t<li>Quality evaluation index : {qei}</li>
\t\t\t<li>Mean CBF (mL 100/g/min) : {meancbf}</li>
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


class SummaryOutputSpec(TraitedSpec):
    out_report = File(exists=True, desc='HTML segment containing summary')


class SummaryInterface(SimpleInterface):
    output_spec = SummaryOutputSpec

    def _run_interface(self, runtime):
        segment = self._generate_segment()
        fname = os.path.join(runtime.cwd, 'report.html')
        with open(fname, 'w') as fobj:
            fobj.write(segment)
        self._results['out_report'] = fname
        return runtime

    def _generate_segment(self):
        raise NotImplementedError


class SubjectSummaryInputSpec(BaseInterfaceInputSpec):
    t1w = InputMultiObject(File(exists=True), desc='T1w structural images')
    t2w = InputMultiObject(File(exists=True), desc='T2w structural images')
    subjects_dir = Directory(desc='FreeSurfer subjects directory')
    subject_id = Str(desc='Subject ID')
    asl = InputMultiObject(traits.Either(
        File(exists=True), traits.List(File(exists=True))),
        desc='ASL functional series')
    std_spaces = traits.List(Str, desc='list of standard spaces')
    nstd_spaces = traits.List(Str, desc='list of non-standard spaces')


class SubjectSummaryOutputSpec(SummaryOutputSpec):
    # This exists to ensure that the summary is run prior to the first ReconAll
    # call, allowing a determination whether there is a pre-existing directory
    subject_id = Str(desc='FreeSurfer subject ID')


class SubjectSummary(SummaryInterface):
    input_spec = SubjectSummaryInputSpec
    output_spec = SubjectSummaryOutputSpec

    def _run_interface(self, runtime):
        if isdefined(self.inputs.subject_id):
            self._results['subject_id'] = self.inputs.subject_id
        return super(SubjectSummary, self)._run_interface(runtime)

    def _generate_segment(self):
        BIDS_NAME = re.compile(
            r'^(.*\/)?'
            '(?P<subject_id>sub-[a-zA-Z0-9]+)'
            '(_(?P<session_id>ses-[a-zA-Z0-9]+))?'
            '(_(?P<task_id>task-[a-zA-Z0-9]+))?'
            '(_(?P<acq_id>acq-[a-zA-Z0-9]+))?'
            '(_(?P<rec_id>rec-[a-zA-Z0-9]+))?'
            '(_(?P<run_id>run-[a-zA-Z0-9]+))?')

        if not isdefined(self.inputs.subjects_dir):
            freesurfer_status = 'Not run'
        else:
            recon = ReconAll(subjects_dir=self.inputs.subjects_dir,
                             subject_id='sub-' + self.inputs.subject_id,
                             T1_files=self.inputs.t1w,
                             flags='-noskullstrip')
            if recon.cmdline.startswith('echo'):
                freesurfer_status = 'Pre-existing directory'
            else:
                freesurfer_status = 'Run by ASLPrep'

        t2w_seg = ''
        if self.inputs.t2w:
            t2w_seg = '(+ {:d} T2-weighted)'.format(len(self.inputs.t2w))

        # Add list of tasks with number of runs
        asl_series = self.inputs.asl if isdefined(self.inputs.asl) else []
        asl_series = [s[0] if isinstance(s, list) else s for s in asl_series]

        counts = Counter(BIDS_NAME.search(series).groupdict()['task_id'][5:]
                         for series in asl_series)

        tasks = ''
        if counts:
            header = '\t\t<ul class="elem-desc">'
            footer = '\t\t</ul>'
            lines = ['\t\t\t<li>Task: {task_id} ({n_runs:d} run{s})</li>'.format(
                     task_id=task_id, n_runs=n_runs, s='' if n_runs == 1 else 's')
                     for task_id, n_runs in sorted(counts.items())]
            tasks = '\n'.join([header] + lines + [footer])

        return SUBJECT_TEMPLATE.format(
            subject_id=self.inputs.subject_id,
            n_t1s=len(self.inputs.t1w),
            t2w=t2w_seg,
            n_asl=len(asl_series),
            tasks=tasks,
            std_spaces=', '.join(self.inputs.std_spaces),
            nstd_spaces=', '.join(self.inputs.nstd_spaces),
            freesurfer_status=freesurfer_status)


class FunctionalSummaryInputSpec(BaseInterfaceInputSpec):
    slice_timing = traits.Enum(False, True, 'TooShort', usedefault=True,
                               desc='Slice timing correction used')
    distortion_correction = traits.Str(desc='Susceptibility distortion correction method',
                                       mandatory=True)
    pe_direction = traits.Enum(None, 'i', 'i-', 'j', 'j-', mandatory=True,
                               desc='Phase-encoding direction detected')
    registration = traits.Enum('FSL', 'FreeSurfer', mandatory=True,
                               desc='Functional/anatomical registration method')
    fallback = traits.Bool(desc='Boundary-based registration rejected')
    registration_dof = traits.Enum(6, 9, 12, desc='Registration degrees of freedom',
                                   mandatory=True)
    registration_init = traits.Enum('register', 'header', mandatory=True,
                                    desc='Whether to initialize registration with the "header"'
                                         ' or by centering the volumes ("register")')
    confounds_file = File(exists=True, mandatory=False, desc='Confounds file')
    qc_file = File(exists=True, desc='qc file')
    tr = traits.Float(desc='Repetition time', mandatory=True)
    #dummy_scans = traits.Either(traits.Int(), None, desc='number of dummy scans specified by user')
    #algo_dummy_scans = traits.Int(desc='number of dummy scans determined by algorithm')


class FunctionalSummary(SummaryInterface):
    input_spec = FunctionalSummaryInputSpec

    def _generate_segment(self):
        dof = self.inputs.registration_dof
        stc = {True: 'Applied',
               False: 'Not applied',
               'TooShort': 'Skipped (too few volumes)'}[self.inputs.slice_timing]
        reg = {
            'FSL': [
                'FSL <code>flirt</code> with boundary-based registration'
                ' (BBR) metric - %d dof' % dof,
                'FSL <code>flirt</code> rigid registration - 6 dof'],
            'FreeSurfer': [
                'FreeSurfer <code>bbregister</code> '
                '(boundary-based registration, BBR) - %d dof' % dof,
                'FreeSurfer <code>mri_coreg</code> - %d dof' % dof],
        }[self.inputs.registration][self.inputs.fallback]
        import pandas as pd
        qcfile = pd.read_csv(self.inputs.qc_file)
        motionparam = "FD : {}, relRMS: {} ".format(round(qcfile['FD'][0], 4),
                                                    round(qcfile['relRMS'][0], 4))
        coregindex = " Dice Index: {}, Jaccard Index: {}, Cross Cor.: {}, Coverage: {} ".format(
                        round(qcfile['coregDC'][0], 4), round(qcfile['coregJC'][0], 4),
                        round(qcfile['coregCC'][0], 4), round(qcfile['coregCOV'][0], 4))
        normindex = " Dice Index: {}, Jaccard Index: {}, Cross Cor.: {}, Coverage: {} ".format(
                        round(qcfile['normDC'][0], 4), round(qcfile['normJC'][0], 4),
                        round(qcfile['normCC'][0], 4), round(qcfile['normCOV'][0], 4))
        qei = "cbf: {},score: {},scrub: {}, basil: {}, pvc: {} ".format(
                round(qcfile['cbfQEI'][0], 4), round(qcfile['scoreQEI'][0], 4),
                round(qcfile['scrubQEI'][0], 4), round(qcfile['basilQEI'][0], 4),
                round(qcfile['pvcQEI'][0], 4))
        meancbf = "GM CBF: {}, WM CBF: {}, GM/WM CBF ratio: {} ".format(
                    round(qcfile['GMmeanCBF'][0], 2), round(qcfile['WMmeanCBF'][0], 2),
                    round(qcfile['Gm_Wm_CBF_ratio'][0], 2))
        negvoxel = "cbf: {}, score: {}, scrub: {}, basil: {}, pvc: {} ".format(
                round(qcfile['NEG_CBF_PERC'][0], 2), round(qcfile['NEG_SCORE_PERC'][0], 2),
                round(qcfile['NEG_SCRUB_PERC'][0], 2), round(qcfile['NEG_BASIL_PERC'][0], 2),
                round(qcfile['NEG_PVC_PERC'][0], 2))
        if self.inputs.pe_direction is None:
            pedir = 'MISSING - Assuming Anterior-Posterior'
        else:
            pedir = {'i': 'Left-Right', 'j': 'Anterior-Posterior'}[self.inputs.pe_direction[0]]

        if isdefined(self.inputs.confounds_file):
            with open(self.inputs.confounds_file) as cfh:
                conflist = cfh.readline().strip('\n').strip()
        else:
            conflist = 'None'
        # the number of dummy scans was specified by the user and
        # it is not equal to the number detected by the algorithm
        
        # the number of dummy scans was not specified by the user
       
        return FUNCTIONAL_TEMPLATE.format(
            pedir=pedir, stc=stc, sdc=self.inputs.distortion_correction, registration=reg,
            confounds=re.sub(r'[\t ]+', ', ', conflist), tr=self.inputs.tr,
            motionparam=motionparam, qei=qei,
            coregindex=coregindex, normindex=normindex, meancbf=meancbf,
            negvoxel=negvoxel)


class AboutSummaryInputSpec(BaseInterfaceInputSpec):
    version = Str(desc='ASLPREP version')
    command = Str(desc='ASLPREP command')
    # Date not included - update timestamp only if version or command changes


class AboutSummary(SummaryInterface):
    input_spec = AboutSummaryInputSpec

    def _generate_segment(self):
        return ABOUT_TEMPLATE.format(version=self.inputs.version,
                                     command=self.inputs.command,
                                     date=time.strftime("%Y-%m-%d %H:%M:%S %z"))