# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
ASLprep base processing workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_aslprep_wf
.. autofunction:: init_single_subject_wf

"""

import sys
import os
from copy import deepcopy

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from .. import config
from ..interfaces import SubjectSummary, AboutSummary, DerivativesDataSink
from .asl import init_asl_preproc_wf
from .asl import init_asl_gepreproc_wf


def init_aslprep_wf():
    """
    Build *ASLPrep*'s pipeline.

    This workflow organizes the execution of aslprep, with a sub-workflow for
    each subject.

    If FreeSurfer's ``recon-all`` is to be run, a corresponding folder is created
    and populated with any needed template subjects under the derivatives folder.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.tests import mock_config
            from aslprep.workflows.base import init_aslprep_wf
            with mock_config():
                wf = init_aslprep_wf()

    """
    from ..niworkflows.engine.workflows import LiterateWorkflow as Workflow
    #from ..niworkflows.interfaces.bids import BIDSFreeSurferDir

    aslprep_wf = Workflow(name='aslprep_wf')
    aslprep_wf.base_dir = config.execution.work_dir

    for subject_id in config.execution.participant_label:
        single_subject_wf = init_single_subject_wf(subject_id)

        single_subject_wf.config['execution']['crashdump_dir'] = str(
            config.execution.output_dir / "aslprep" / "-".join(("sub", subject_id))
            / "log" / config.execution.run_uuid
        )
        for node in single_subject_wf._get_all_nodes():
            node.config = deepcopy(single_subject_wf.config)
        
        aslprep_wf.add_nodes([single_subject_wf])

        # Dump a copy of the config file into the log directory
        log_dir = config.execution.output_dir / 'aslprep' / 'sub-{}'.format(subject_id) \
            / 'log' / config.execution.run_uuid
        log_dir.mkdir(exist_ok=True, parents=True)
        config.to_filename(log_dir / 'aslprep.toml')

    return aslprep_wf


def init_single_subject_wf(subject_id):
    """
    Organize the preprocessing pipeline for a single subject.

    It collects and reports information about the subject, and prepares
    sub-workflows to perform anatomical and functional preprocessing.
    Anatomical preprocessing is performed in a single workflow, regardless of
    the number of sessions.
    Functional preprocessing is performed using a separate workflow for each
    individual ASL series.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.tests import mock_config
            from aslprep.workflows.base import init_single_subject_wf
            with mock_config():
                wf = init_single_subject_wf('01')

    Parameters
    ----------
    subject_id : :obj:`str`
        Subject label for this single-subject workflow.

    Inputs
    ------
    subjects_dir : :obj:`str`
        FreeSurfer's ``$SUBJECTS_DIR``.

    """
    from ..niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from ..niworkflows.interfaces.bids import BIDSInfo, BIDSDataGrabber
    from ..niworkflows.interfaces.nilearn import NILEARN_VERSION
    from ..niworkflows.utils.bids import collect_data
    from ..niworkflows.utils.misc import fix_multi_T1w_source_name
    from ..niworkflows.utils.spaces import Reference
    from ..smriprep.workflows.anatomical import init_anat_preproc_wf

    name = "single_subject_%s_wf" % subject_id
    subject_data = collect_data(
        config.execution.bids_dir,
        subject_id,
        config.execution.task_id,
        config.execution.echo_idx,
        bids_filters=config.execution.bids_filters)[0]

    if 'flair' in config.workflow.ignore:
        subject_data['flair'] = []
    if 't2w' in config.workflow.ignore:
        subject_data['t2w'] = []

    anat_only = config.workflow.anat_only
    # Make sure we always go through these two checks
    if not anat_only and not subject_data['asl']:
        task_id = config.execution.task_id
        raise RuntimeError(
            "No ASL images found for participant {} and task {}. "
            "All workflows require ASL images.".format(
                subject_id, task_id if task_id else '<all>')
        )

    if not subject_data['t1w']:
        raise Exception("No T1w images found for participant {}. "
                        "All workflows require T1w images.".format(subject_id))

    workflow = Workflow(name=name)
    spaces = config.workflow.spaces
    output_dir = str(config.execution.output_dir)

    inputnode = pe.Node(niu.IdentityInterface(fields=['subjects_dir']),
                        name='inputnode')

    bidssrc = pe.Node(BIDSDataGrabber(subject_data=subject_data,
                                      anat_only=anat_only,
                                      subject_id=subject_id),
                      name='bidssrc')

    bids_info = pe.Node(BIDSInfo(
        bids_dir=config.execution.bids_dir, bids_validate=False), name='bids_info')

    summary = pe.Node(SubjectSummary(std_spaces=spaces.get_spaces(nonstandard=False),
                                     nstd_spaces=spaces.get_spaces(standard=False)),
                      name='summary', run_without_submitting=True)

    about = pe.Node(AboutSummary(version=config.environment.version,
                                 command=' '.join(sys.argv)),
                    name='about', run_without_submitting=True)

    ds_report_summary = pe.Node(
        DerivativesDataSink(base_directory=output_dir, desc='summary', datatype="figures",
                            dismiss_entities=("echo",)),
        name='ds_report_summary', run_without_submitting=True)

    ds_report_about = pe.Node(
        DerivativesDataSink(base_directory=output_dir, desc='about', datatype="figures",
                            dismiss_entities=("echo",)),
        name='ds_report_about', run_without_submitting=True)
    anat_derivatives = config.execution.anat_derivatives
    if anat_derivatives:
        from ..smriprep.utils.bids import collect_derivatives
        std_spaces = spaces.get_spaces(nonstandard=False, dim=(3,))
        anat_derivatives = collect_derivatives(
            derivative_dir=anat_derivatives.absolute(),
            subject_id=subject_id,
            std_spaces=std_spaces,
            freesurfer=None
        )
        if anat_derivatives is None:
            config.loggers.workflow.warning(f"""\
Attempted to access pre-existing anatomical derivatives at \
<{config.execution.anat_derivatives}>, however not all expectations of aslprep \
were met (for participant <{subject_id}>, spaces <{', '.join(std_spaces)}>, \
""")

    workflow.__desc__ = """
Results included in this manuscript come from preprocessing
performed using *aslprep* {aslprep_ver},
which is based on *Nipype* {nipype_ver}
(@nipype1; @nipype2; RRID:SCR_002502).

""".format(aslprep_ver=config.environment.version,
           nipype_ver=config.environment.nipype_version)  
    workflow.__postdesc__ = """


Many internal operations of *aslprep* use
*Nilearn* {nilearn_ver} [@nilearn, RRID:SCR_001362],
mostly within the functional processing workflow.
For more details of the pipeline, see [the section corresponding
to workflows in *aslprep*'s documentation]\
(https://aslprep.readthedocs.io/en/latest/workflows.html \
"aslprep's documentation").


### Copyright Waiver

The above boilerplate text was automatically generated by aslprep
with the express intention that users should copy and paste this
text into their manuscripts *unchanged*.
It is released under the [CC0]\
(https://creativecommons.org/publicdomain/zero/1.0/) license.

### References

""".format(nilearn_ver=NILEARN_VERSION)

    # Preprocessing of T1w (includes registration to MNI)
    anat_preproc_wf = init_anat_preproc_wf(
        bids_root=str(config.execution.bids_dir),
        freesurfer=None,
        debug=config.execution.debug is True,
        existing_derivatives=anat_derivatives,
        hires=config.workflow.hires,
        longitudinal=config.workflow.longitudinal,
        omp_nthreads=config.nipype.omp_nthreads,
        output_dir=output_dir,
        skull_strip_fixed_seed=config.workflow.skull_strip_fixed_seed,
        skull_strip_mode=config.workflow.skull_strip_t1w,
        skull_strip_template=Reference.from_string(
            config.workflow.skull_strip_template)[0],
        spaces=spaces,
        t1w=subject_data['t1w'],
    )

    workflow.connect([
        (inputnode, anat_preproc_wf, [('subjects_dir', 'inputnode.subjects_dir')]),
        (bidssrc, bids_info, [(('t1w', fix_multi_T1w_source_name), 'in_file')]),
        (inputnode, summary, [('subjects_dir', 'subjects_dir')]),
        (bidssrc, summary, [('t1w', 't1w'),
                            ('t2w', 't2w'),
                            ('asl', 'asl')]),
        (bids_info, summary, [('subject', 'subject_id')]),
        (bids_info, anat_preproc_wf, [(('subject', _prefix), 'inputnode.subject_id')]),
        (bidssrc, anat_preproc_wf, [('t1w', 'inputnode.t1w'),
                                    ('t2w', 'inputnode.t2w'),
                                    ('roi', 'inputnode.roi'),
                                    ('flair', 'inputnode.flair')]),
        (bidssrc, ds_report_summary, [(('t1w', fix_multi_T1w_source_name), 'source_file')]),
        (summary, ds_report_summary, [('out_report', 'in_file')]),
        (bidssrc, ds_report_about, [(('t1w', fix_multi_T1w_source_name), 'source_file')]),
        (about, ds_report_about, [('out_report', 'in_file')]),
    ])

    # Overwrite ``out_path_base`` of smriprep's DataSinks
    for node in workflow.list_node_names():
        if node.split('.')[-1].startswith('ds_'):
            workflow.get_node(node).interface.out_path_base = 'aslprep'

    if anat_only:
        return workflow

    # Append the functional section to the existing anatomical exerpt
    # That way we do not need to stream down the number of asl datasets
    anat_preproc_wf.__postdesc__ = (anat_preproc_wf.__postdesc__ or '') + """

Functional data preprocessing

: For each of the {num_asl} ASL runs found per subject (across all
tasks and sessions), the following preprocessing was performed.
""".format(num_asl=len(subject_data['asl']))

    for asl_file in subject_data['asl']:
        if check_img(img=asl_file) > 5: # if number of volume of ASL is less than 5. motion slicetiming etc will be skipped
            asl_preproc_wf = init_asl_preproc_wf(asl_file)
            workflow.connect([
             (anat_preproc_wf, asl_preproc_wf,
             [('outputnode.t1w_preproc', 'inputnode.t1w_preproc'),
              ('outputnode.t1w_mask', 'inputnode.t1w_mask'),
              ('outputnode.t1w_dseg', 'inputnode.t1w_dseg'),
              ('outputnode.t1w_tpms', 'inputnode.t1w_tpms'),
              ('outputnode.template', 'inputnode.template'),
              ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
              ('outputnode.std2anat_xfm', 'inputnode.std2anat_xfm'),
              ]),
            ])
        else:
            asl_preproc_wf = init_asl_gepreproc_wf(asl_file)
            workflow.connect([
             (anat_preproc_wf, asl_preproc_wf,
             [('outputnode.t1w_preproc', 'inputnode.t1w_preproc'),
              ('outputnode.t1w_mask', 'inputnode.t1w_mask'),
              ('outputnode.t1w_dseg', 'inputnode.t1w_dseg'),
              ('outputnode.t1w_tpms', 'inputnode.t1w_tpms'),
              ('outputnode.template', 'inputnode.template'),
              ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
              ('outputnode.std2anat_xfm', 'inputnode.std2anat_xfm')]),
            ])

    return workflow


def _prefix(subid):
    return subid if subid.startswith('sub-') else f'sub-{subid}'


def _pop(inlist):
    if isinstance(inlist, (list, tuple)):
        return inlist[0]
    return inlist

def check_img(img):
    # get the 4th dimension
    import numpy as np 
    import nibabel as nb 
    ss=nb.load(img).get_fdata().shape
    if len(ss) == 3:
        ss=np.hstack([ss,0])
    return ss[3]