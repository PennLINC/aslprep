# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
ASLPrep base processing workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_aslprep_wf
.. autofunction:: init_single_subject_wf

"""

import sys
import os
from copy import deepcopy

from nipype import __version__ as nipype_ver
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from ..niworkflows.interfaces.nilearn import NILEARN_VERSION

from ..niworkflows.engine.workflows import LiterateWorkflow as Workflow
from ..niworkflows.interfaces.bids import (
    BIDSInfo, BIDSDataGrabber, BIDSFreeSurferDir
)
from ..niworkflows.utils.bids import collect_data
from ..niworkflows.utils.misc import fix_multi_T1w_source_name
from ..smriprep.workflows.anatomical import init_anat_preproc_wf

from ..interfaces import SubjectSummary, AboutSummary, DerivativesDataSink
from ..__about__ import __version__
from .bold import init_func_preproc_wf


def init_aslprep_wf(
    anat_only,
    pcasl,
    #aroma_melodic_dim,
    bold2t1w_dof,
    cifti_output,
    debug,
    dummy_scans,
    dummy_vols,
    echo_idx,
    #err_on_aroma_warn,
    fmap_bspline,
    fmap_demean,
    force_syn,
    freesurfer,
    fs_subjects_dir,
    hires,
    ignore,
    layout,
    longitudinal,
    low_mem,
    medial_surface_nan,
    omp_nthreads,
    output_dir,
    #regressors_all_comps,
    #regressors_dvars_th,
    #regressors_fd_th,
    run_uuid,
    skull_strip_fixed_seed,
    skull_strip_template,
    spaces,
    subject_list,
    t2s_coreg,
    task_id,
    #use_aroma,
    use_bbr,
    use_syn,
    work_dir,
    smooth_kernel,
    bids_filters,
):
    """
    Build *ASLPrep*'s pipeline.

    This workflow organizes the execution of ASLPREP, with a sub-workflow for
    each subject.

    If FreeSurfer's ``recon-all`` is to be run, a corresponding folder is created
    and populated with any needed template subjects under the derivatives folder.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            import os
            from collections import namedtuple, OrderedDict
            BIDSLayout = namedtuple('BIDSLayout', ['root'])
            from aslprep.workflows.base import init_aslprep_wf
            from ..niworkflows.utils.spaces import Reference, SpatialReferences

            os.environ['FREESURFER_HOME'] = os.getcwd()
            wf = init_aslprep_wf(
                anat_only=False,
                aroma_melodic_dim=-200,
                bold2t1w_dof=9,
                cifti_output=False,
                debug=False,
                dummy_scans=None,
                echo_idx=None,
                err_on_aroma_warn=False,
                fmap_bspline=False,
                fmap_demean=True,
                force_syn=True,
                freesurfer=True,
                fs_subjects_dir=None,
                hires=True,
                ignore=[],
                layout=BIDSLayout('.'),
                longitudinal=False,
                low_mem=False,
                medial_surface_nan=False,
                omp_nthreads=1,
                output_dir='.',
                regressors_all_comps=False,
                regressors_dvars_th=1.5,
                regressors_fd_th=0.5,
                run_uuid='X',
                skull_strip_fixed_seed=False,
                skull_strip_template=Reference('OASIS30ANTs'),
                spaces=SpatialReferences(
                    spaces=['MNI152Lin',
                            ('fsaverage', {'density': '10k'}),
                            'T1w',
                            'fsnative'],
                    checkpoint=True),
                subject_list=['aslpreptest'],
                t2s_coreg=False,
                task_id='',
                use_aroma=False,
                use_bbr=True,
                use_syn=True,
                work_dir='.',
                bids_filters=None,
            )


    Parameters
    ----------
    anat_only : :obj:`bool`
        Disable functional workflows
    bold2t1w_dof : 6, 9 or 12
        Degrees-of-freedom for BOLD-T1w registration
    cifti_output : :obj:`bool`
        Generate bold CIFTI file in output spaces
    debug : :obj:`bool`
        Enable debugging outputs
    dummy_scans : :obj:`int` or None
        Number of volumes to consider as non steady state
    echo_idx : :obj:`int` or None
        Index of echo to preprocess in multiecho BOLD series,
        or ``None`` to preprocess all
    err_on_aroma_warn : :obj:`bool`
        Do not fail on ICA-AROMA errors
    fmap_bspline : :obj:`bool`
        **Experimental**: Fit B-Spline field using least-squares
    fmap_demean : :obj:`bool`
        Demean voxel-shift map during unwarp
    force_syn : :obj:`bool`
        **Temporary**: Always run SyN-based SDC
    freesurfer : :obj:`bool`
        Enable FreeSurfer surface reconstruction (may increase runtime)
    hires : :obj:`bool`
        Enable sub-millimeter preprocessing in FreeSurfer
    ignore : :obj:`list`
        Preprocessing steps to skip (may include "slicetiming", "fieldmaps")
    layout : :py:class:`~bids.layout.BIDSLayout`
        BIDS dataset layout
    longitudinal : :obj:`bool`
        Treat multiple sessions as longitudinal (may increase runtime)
        See sub-workflows for specific differences
    low_mem : :obj:`bool`
        Write uncompressed .nii files in some cases to reduce memory usage
    medial_surface_nan : :obj:`bool`
        Replace medial wall values with NaNs on functional GIFTI files
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    output_dir : :obj:`str`
        Directory in which to save derivatives
    regressors_all_comps : :obj:`bool`
        Return all CompCor component time series instead of the top fraction
    regressors_dvars_th : :obj:`float`
        Criterion for flagging DVARS outliers
    regressors_fd_th : :obj:`float`
        Criterion for flagging framewise displacement outliers
    run_uuid : :obj:`str`
        Unique identifier for execution instance
    skull_strip_template : tuple
        Name of target template for brain extraction with ANTs' ``antsBrainExtraction``,
        and corresponding dictionary of output-space modifiers.
    skull_strip_fixed_seed : :obj:`bool`
        Do not use a random seed for skull-stripping - will ensure
        run-to-run replicability when used with --omp-nthreads 1
    spaces : :py:class:`~..niworkflows.utils.spaces.SpatialReferences`
        A container for storing, organizing, and parsing spatial normalizations. Composed of
        :py:class:`~..niworkflows.utils.spaces.Reference` objects representing spatial references.
        Each ``Reference`` contains a space, which is a string of either TemplateFlow template IDs
        (e.g., ``MNI152Lin``, ``MNI152NLin6Asym``, ``MNIPediatricAsym``), nonstandard references
        (e.g., ``T1w`` or ``anat``, ``sbref``, ``run``, etc.), or a custom template located in
        the TemplateFlow root directory. Each ``Reference`` may also contain a spec, which is a
        dictionary with template specifications (e.g., a specification of ``{'resolution': 2}``
        would lead to resampling on a 2mm resolution of the space).
    subject_list : :obj:`list`
        List of subject labels
    t2s_coreg : :obj:`bool`
        For multi-echo EPI, use the calculated T2*-map for T2*-driven coregistration
    task_id : :obj:`str` or None
        Task ID of BOLD series to preprocess, or ``None`` to preprocess all
    use_aroma : :obj:`bool`
        Perform ICA-AROMA on MNI-resampled functional series
    use_bbr : :obj:`bool` or None
        Enable/disable boundary-based registration refinement.
        If ``None``, test BBR result for distortion before accepting.
    use_syn : :obj:`bool`
        **Experimental**: Enable ANTs SyN-based susceptibility distortion correction (SDC).
        If fieldmaps are present and enabled, this is not run, by default.
    work_dir : :obj:`str`
        Directory in which to store workflow execution state and temporary files
    bids_filters : :obj:`dict`
        Provides finer specification of the pipeline input files using pybids entities filters.
        A dict with the following structure {<suffix>:{<entity>:<filter>,...},...}

    """
    aslprep_wf = Workflow(name='aslprep_wf')
    aslprep_wf.base_dir = work_dir

    if freesurfer:
        fsdir = pe.Node(
            BIDSFreeSurferDir(
                derivatives=output_dir,
                freesurfer_home=os.getenv('FREESURFER_HOME'),
                spaces=spaces.get_fs_spaces()),
            name='fsdir_run_' + run_uuid.replace('-', '_'), run_without_submitting=True)
        if fs_subjects_dir is not None:
            fsdir.inputs.subjects_dir = str(fs_subjects_dir.absolute())

    reportlets_dir = os.path.join(work_dir, 'reportlets')
    for subject_id in subject_list:
        single_subject_wf = init_single_subject_wf(
            anat_only=anat_only,
            #aroma_melodic_dim=aroma_melodic_dim,
            bold2t1w_dof=bold2t1w_dof,
            cifti_output=cifti_output,
            debug=debug,
            dummy_scans=dummy_scans,
            echo_idx=echo_idx,
            #err_on_aroma_warn=err_on_aroma_warn,
            fmap_bspline=fmap_bspline,
            fmap_demean=fmap_demean,
            force_syn=force_syn,
            freesurfer=freesurfer,
            hires=hires,
            ignore=ignore,
            dummy_vols=dummy_vols,
            layout=layout,
            longitudinal=longitudinal,
            low_mem=low_mem,
            medial_surface_nan=medial_surface_nan,
            name="single_subject_" + subject_id + "_wf",
            omp_nthreads=omp_nthreads,
            output_dir=output_dir,
            #regressors_all_comps=regressors_all_comps,
            #regressors_dvars_th=regressors_dvars_th,
            #regressors_fd_th=regressors_fd_th,
            reportlets_dir=reportlets_dir,
            skull_strip_fixed_seed=skull_strip_fixed_seed,
            skull_strip_template=skull_strip_template,
            spaces=spaces,
            subject_id=subject_id,
            t2s_coreg=t2s_coreg,
            task_id=task_id,
            pcasl=pcasl,
            smooth_kernel=smooth_kernel,
            #use_aroma=use_aroma,
            use_bbr=use_bbr,
            use_syn=use_syn,
            bids_filters=bids_filters,
        )

        single_subject_wf.config['execution']['crashdump_dir'] = (
            os.path.join(output_dir, "aslprep", "sub-" + subject_id, 'log', run_uuid)
        )
        for node in single_subject_wf._get_all_nodes():
            node.config = deepcopy(single_subject_wf.config)
        if freesurfer:
            aslprep_wf.connect(fsdir, 'subjects_dir',
                                single_subject_wf, 'inputnode.subjects_dir')
        else:
            aslprep_wf.add_nodes([single_subject_wf])

    return aslprep_wf


def init_single_subject_wf(
    anat_only,
    #aroma_melodic_dim,
    bold2t1w_dof,
    cifti_output,
    debug,
    dummy_scans,
    echo_idx,
    #err_on_aroma_warn,
    fmap_bspline,
    fmap_demean,
    force_syn,
    freesurfer,
    hires,
    ignore,
    layout,
    dummy_vols,
    longitudinal,
    low_mem,
    medial_surface_nan,name,
    omp_nthreads,
    output_dir,
    reportlets_dir,
    smooth_kernel,
    #regressors_all_comps,
    #regressors_dvars_th,
    #regressors_fd_th,
    skull_strip_fixed_seed,
    skull_strip_template,
    spaces,
    subject_id,
    t2s_coreg,
    task_id,
    #use_aroma,
    use_bbr,
    use_syn,
    pcasl,
    bids_filters,
):
    """
    This workflow organizes the preprocessing pipeline for a single subject.

    It collects and reports information about the subject, and prepares
    sub-workflows to perform anatomical and functional preprocessing.
    Anatomical preprocessing is performed in a single workflow, regardless of
    the number of sessions.
    Functional preprocessing is performed using a separate workflow for each
    individual BOLD series.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from collections import namedtuple
            from ..niworkflows.utils.spaces import Reference, SpatialReferences
            from aslprep.workflows.base import init_single_subject_wf

            BIDSLayout = namedtuple('BIDSLayout', ['root'])
            wf = init_single_subject_wf(
                anat_only=False,
                aroma_melodic_dim=-200,
                bold2t1w_dof=9,
                cifti_output=False,
                debug=False,
                dummy_scans=None,
                echo_idx=None,
                err_on_aroma_warn=False,
                fmap_bspline=False,
                fmap_demean=True,
                force_syn=True,
                freesurfer=True,
                hires=True,
                ignore=[],
                layout=BIDSLayout('.'),
                longitudinal=False,
                low_mem=False,
                medial_surface_nan=False,
                name='single_subject_wf',
                omp_nthreads=1,
                output_dir='.',
                reportlets_dir='.',
                regressors_all_comps=False,
                regressors_dvars_th=1.5,
                regressors_fd_th=0.5,
                skull_strip_fixed_seed=False,
                skull_strip_template=Reference('OASIS30ANTs'),
                spaces=SpatialReferences(
                    spaces=['MNI152Lin',
                            ('fsaverage', {'density': '10k'}),
                            'T1w',
                            'fsnative'],
                    checkpoint=True),
                subject_id='test',
                t2s_coreg=False,
                task_id='',
                use_aroma=False,
                use_bbr=True,
                use_syn=True,
                bids_filters=None,
            )

    Parameters
    ----------
    anat_only : :obj:`bool`
        Disable functional workflows
    aroma_melodic_dim : :obj:`int`
        Maximum number of components identified by MELODIC within ICA-AROMA
        (default is -200, i.e., no limitation).
    bold2t1w_dof : 6, 9 or 12
        Degrees-of-freedom for BOLD-T1w registration
    cifti_output : :obj:`bool`
        Generate bold CIFTI file in output spaces
    debug : :obj:`bool`
        Enable debugging outputs
    dummy_scans : :obj:`int` or None
        Number of volumes to consider as non steady state
    echo_idx : :obj:`int` or None
        Index of echo to preprocess in multiecho BOLD series,
        or ``None`` to preprocess all
    err_on_aroma_warn : :obj:`bool`
        Do not fail on ICA-AROMA errors
    fmap_bspline : :obj:`bool`
        **Experimental**: Fit B-Spline field using least-squares
    fmap_demean : :obj:`bool`
        Demean voxel-shift map during unwarp
    force_syn : :obj:`bool`
        **Temporary**: Always run SyN-based SDC
    freesurfer : :obj:`bool`
        Enable FreeSurfer surface reconstruction (may increase runtime)
    hires : :obj:`bool`
        Enable sub-millimeter preprocessing in FreeSurfer
    ignore : :obj:`list`
        Preprocessing steps to skip (may include "slicetiming", "fieldmaps")
    layout : :py:class:`~bids.layout.BIDSLayout`
        BIDS dataset layout
    longitudinal : :obj:`bool`
        Treat multiple sessions as longitudinal (may increase runtime)
        See sub-workflows for specific differences
    low_mem : :obj:`bool`
        Write uncompressed .nii files in some cases to reduce memory usage
    medial_surface_nan : :obj:`bool`
        Replace medial wall values with NaNs on functional GIFTI files
    name : :obj:`str`
        Name of workflow
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    output_dir : :obj:`str`
        Directory in which to save derivatives
    reportlets_dir : :obj:`str`
        Directory in which to save reportlets
    regressors_all_comps : :obj:`bool`
        Return all CompCor component time series instead of the top fraction
    regressors_dvars_th : :obj:`float`
        Criterion for flagging DVARS outliers
    regressors_fd_th : :obj:`float`
        Criterion for flagging framewise displacement outliers
    skull_strip_fixed_seed : :obj:`bool`
        Do not use a random seed for skull-stripping - will ensure
        run-to-run replicability when used with --omp-nthreads 1
    skull_strip_template : tuple
        Name of target template for brain extraction with ANTs' ``antsBrainExtraction``,
        and corresponding dictionary of output-space modifiers.
    subject_id : :obj:`str`
        List of subject labels
    t2s_coreg : :obj:`bool`
        For multi-echo EPI, use the calculated T2*-map for T2*-driven coregistration
    spaces : :py:class:`~..niworkflows.utils.spaces.SpatialReferences`
        A container for storing, organizing, and parsing spatial normalizations. Composed of
        :py:class:`~..niworkflows.utils.spaces.Reference` objects representing spatial references.
        Each ``Reference`` contains a space, which is a string of either TemplateFlow template IDs
        (e.g., ``MNI152Lin``, ``MNI152NLin6Asym``, ``MNIPediatricAsym``), nonstandard references
        (e.g., ``T1w`` or ``anat``, ``sbref``, ``run``, etc.), or a custom template located in
        the TemplateFlow root directory. Each ``Reference`` may also contain a spec, which is a
        dictionary with template specifications (e.g., a specification of ``{'resolution': 2}``
        would lead to resampling on a 2mm resolution of the space).
    task_id : :obj:`str` or None
        Task ID of BOLD series to preprocess, or ``None`` to preprocess all
    use_aroma : :obj:`bool`
        Perform ICA-AROMA on MNI-resampled functional series
    use_bbr : :obj:`bool` or None
        Enable/disable boundary-based registration refinement.
        If ``None``, test BBR result for distortion before accepting.
    use_syn : :obj:`bool`
        **Experimental**: Enable ANTs SyN-based susceptibility distortion correction (SDC).
        If fieldmaps are present and enabled, this is not run, by default.
    bids_filters : :obj:`dict`
        Provides finer specification of the pipeline input files using pybids entities filters.
        A dict with the following structure {<suffix>:{<entity>:<filter>,...},...}

    Inputs
    ------
    subjects_dir : :obj:`str`
        FreeSurfer's ``$SUBJECTS_DIR``.

    """
    if name in ('single_subject_wf', 'single_subject_aslpreptest_wf'):
        # for documentation purposes
        subject_data = {
            't1w': ['/completely/made/up/path/sub-01_T1w.nii.gz'],
            'asl': ['/completely/made/up/path/sub-01_task-nback_bold.nii.gz']
        }
    else:
        subject_data = collect_data(layout, subject_id, task_id, echo_idx,
                                    bids_filters=None)[0]

    # Make sure we always go through these two checks
    if not anat_only and subject_data['asl'] == []:
        raise Exception("No ASL images found for participant {} and task {}. "
                        "All workflows require ASL images.".format(
                            subject_id, task_id if task_id else '<all>'))

    if not subject_data['t1w']:
        raise Exception("No T1w images found for participant {}. "
                        "All workflows require T1w images.".format(subject_id))

    workflow = Workflow(name=name)
    
    inputnode = pe.Node(niu.IdentityInterface(fields=['subjects_dir']),
                        name='inputnode')

    bidssrc = pe.Node(BIDSDataGrabber(subject_data=subject_data,
                                      anat_only=anat_only,
                                      subject_id=subject_id),
                      name='bidssrc')

    bids_info = pe.Node(BIDSInfo(
        bids_dir=layout.root, bids_validate=False), name='bids_info')

    summary = pe.Node(SubjectSummary(std_spaces=spaces.get_spaces(nonstandard=False),
                                     nstd_spaces=spaces.get_spaces(standard=False)),
                      name='summary', run_without_submitting=True)

    about = pe.Node(AboutSummary(version=__version__,
                                 command=' '.join(sys.argv)),
                    name='about', run_without_submitting=True)

    ds_report_summary = pe.Node(
        DerivativesDataSink(base_directory=reportlets_dir,
                            desc='summary', keep_dtype=True),
        name='ds_report_summary', run_without_submitting=True)

    ds_report_about = pe.Node(
        DerivativesDataSink(base_directory=reportlets_dir,
                            desc='about', keep_dtype=True),
        name='ds_report_about', run_without_submitting=True)

    workflow.__desc__ = """
Results included in this manuscript come from preprocessing
performed using *ASLPrep* {aslprep_ver}
(RRID:SCR_016216),
which is based on *Nipype* {nipype_ver}
(@nipype1; @nipype2; RRID:SCR_002502).

""".format(aslprep_ver=__version__, nipype_ver=nipype_ver)
    workflow.__postdesc__ = """

Many internal operations of *ASLPrep* use
*Nilearn* {nilearn_ver} [@nilearn, RRID:SCR_001362],
mostly within the functional processing workflow.
For more details of the pipeline, see [the section corresponding
to workflows in *ASLPrep*'s documentation]\
(https://aslprep.readthedocs.io/en/latest/workflows.html \
"ASLPrep's documentation").


### Copyright Waiver

The above boilerplate text was automatically generated by ASLPrep
with the express intention that users should copy and paste this
text into their manuscripts *unchanged*.
It is released under the [CC0]\
(https://creativecommons.org/publicdomain/zero/1.0/) license.

### References

""".format(nilearn_ver=NILEARN_VERSION)

    # Preprocessing of T1w (includes registration to MNI)
    anat_preproc_wf = init_anat_preproc_wf(
        bids_root=layout.root,
        debug=debug,
        freesurfer=freesurfer,
        hires=hires,
        longitudinal=longitudinal,
        name="anat_preproc_wf",
        num_t1w=len(subject_data['t1w']),
        omp_nthreads=omp_nthreads,
        output_dir=output_dir,
        reportlets_dir=reportlets_dir,
        spaces=spaces,
        skull_strip_fixed_seed=skull_strip_fixed_seed,
        skull_strip_template=skull_strip_template,
    )

    workflow.connect([
        (inputnode, anat_preproc_wf, [('subjects_dir', 'inputnode.subjects_dir')]),
        (bidssrc, bids_info, [(('t1w', fix_multi_T1w_source_name), 'in_file')]),
        (inputnode, summary, [('subjects_dir', 'subjects_dir')]),
        (bidssrc, summary, [('t1w', 't1w'),
                            ('t2w', 't2w'),
                            ('asl', 'bold')]),
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

    for bold_file in subject_data['asl']:
        #import os 
        #file1=os.path.abspath(bold_file)
        #aslcontext=file1.replace('.nii.gz','_ASLContext.tsv')
        func_preproc_wf = init_func_preproc_wf(
            #aroma_melodic_dim=aroma_melodic_dim,
            bold2t1w_dof=bold2t1w_dof,
            bold_file=bold_file,
            cifti_output=cifti_output,
            debug=debug,
            dummy_scans=dummy_scans,
            #err_on_aroma_warn=err_on_aroma_warn,
            fmap_bspline=fmap_bspline,
            fmap_demean=fmap_demean,
            force_syn=force_syn,
            freesurfer=freesurfer,
            ignore=ignore,
            layout=layout,
            low_mem=low_mem,
            dummy_vols=dummy_vols,
            medial_surface_nan=medial_surface_nan,
            num_bold=len(subject_data['asl']),
            omp_nthreads=omp_nthreads,
            output_dir=output_dir,
            reportlets_dir=reportlets_dir,
            pcasl=False,
            smooth_kernel=smooth_kernel,
            #aslcontext=aslcontext,
            #regressors_all_comps=regressors_all_comps,
            #regressors_fd_th=regressors_fd_th,
            #regressors_dvars_th=regressors_dvars_th,
            spaces=spaces,
            t2s_coreg=t2s_coreg,
            #use_aroma=use_aroma,
            use_bbr=use_bbr,
            use_syn=use_syn,
        )

        workflow.connect([
            (anat_preproc_wf, func_preproc_wf,
             [(('outputnode.t1w_preproc', _pop), 'inputnode.t1w_preproc'),
              ('outputnode.t1w_brain', 'inputnode.t1w_brain'),
              ('outputnode.t1w_mask', 'inputnode.t1w_mask'),
              ('outputnode.t1w_dseg', 'inputnode.t1w_dseg'),
              ('outputnode.t1w_aseg', 'inputnode.t1w_aseg'),
              ('outputnode.t1w_aparc', 'inputnode.t1w_aparc'),
              ('outputnode.t1w_tpms', 'inputnode.t1w_tpms'),
              ('outputnode.template', 'inputnode.template'),
              ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
              ('outputnode.std2anat_xfm', 'inputnode.std2anat_xfm'),
              ('outputnode.joint_template', 'inputnode.joint_template'),
              ('outputnode.joint_anat2std_xfm', 'inputnode.joint_anat2std_xfm'),
              ('outputnode.joint_std2anat_xfm', 'inputnode.joint_std2anat_xfm'),
              # Undefined if --fs-no-reconall, but this is safe
              ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
              ('outputnode.subject_id', 'inputnode.subject_id'),
              ('outputnode.t1w2fsnative_xfm', 'inputnode.t1w2fsnative_xfm'),
              ('outputnode.fsnative2t1w_xfm', 'inputnode.fsnative2t1w_xfm')]),
        ])

    return workflow


def _prefix(subid):
    if subid.startswith('sub-'):
        return subid
    return '-'.join(('sub', subid))


def _pop(inlist):
    if isinstance(inlist, (list, tuple)):
        return inlist[0]
    return inlist
