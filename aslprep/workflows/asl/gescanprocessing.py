# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

from ... import config

import os

import nibabel as nb
from nipype.interfaces.fsl import Split as FSLSplit
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from ...utils.meepi import combine_meepi_source

from ...interfaces import DerivativesDataSink
from ...interfaces.reports import FunctionalSummary

# asl workflows
from .confounds import init_asl_confs_wf, init_carpetplot_wf
from .hmc import init_asl_hmc_wf
from .stc import init_asl_stc_wf
from .t2s import init_asl_t2s_wf
from .registration import init_asl_t1_trans_wf, init_asl_reg_wf
from .resampling import (
    init_asl_surf_wf,
    init_asl_std_trans_wf,
    init_asl_preproc_trans_wf,
)
from .cbf import (
    init_cbf_compt_wf,
    init_cbfqc_compt_wf,
    init_cbfplot_wf,
    init_cbfroiquant_wf,
    init_gecbf_compt_wf)
from .outputs import init_asl_derivatives_wf

from ...interfaces.cbf_computation import refinemask
from .ge_utils import ( init_asl_geref_wf, init_asl_gereg_wf,
             init_asl_t1_getrans_wf,init_asl_gestd_trans_wf) 


def init_asl_gepreproc_wf(asl_file):
    """
    This workflow controls the functional preprocessing stages of *aslprep*.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.tests import mock_config
            from aslprep import config
            from aslprep.workflows.asl.base import init_asl_preproc_wf
            with mock_config():
                asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_task-restEyesOpen_asl.nii.gz'
                wf = init_asl_preproc_wf(str(asl_file))

    Parameters
    ----------
    asl_file
        asl series NIfTI file

    Inputs
    ------
    asl_file
        asl series NIfTI file
    t1w_preproc
        Bias-corrected structural template image
    t1w_mask
        Mask of the skull-stripped template image
    t1w_dseg
        Segmentation of preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF)
    t1w_asec
        Segmentation of structural image, done with FreeSurfer.
    t1w_aparc
        Parcellation of structural image, done with FreeSurfer.
    t1w_tpms
        List of tissue probability maps in T1w space
    template
        List of templates to target
    anat2std_xfm
        List of transform files, collated with templates
    std2anat_xfm
        List of inverse transform files, collated with templates
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    t1w2fsnative_xfm
        LTA-style affine matrix translating from T1w to FreeSurfer-conformed subject space
    fsnative2t1w_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed subject space to T1w

    Outputs
    -------
    asl_t1
        asl series, resampled to T1w space
    asl_mask_t1
        asl series mask in T1w space
    asl_std
        asl series, resampled to template space
    asl_mask_std
        asl series mask in template space
    confounds
        TSV of confounds
    surfaces
        asl series, resampled to FreeSurfer surfaces
    asl_cifti
        asl CIFTI image
    cifti_variant
        combination of target spaces for `asl_cifti`
    cbf_t1
        cbf times series in T1w space
    meancbf_t1
        mean cbf   in T1w space
    scorecbf_t1
        scorecbf times series in T1w space
    avgscorecbf_t1
        mean score cbf  in T1w space
    scrub_t1, pv_t1, basil_t1
        scrub, parital volume corrected and basil cbf   in T1w space
    cbf_std
        cbf times series in template space
    meancbf_std
        mean cbf   in template space
    scorecbf_std
        scorecbf times series in template space
    avgscorecbf_std
        mean score cbf  in template space
    scrub_std, pv_std, basil_std
        scrub, parital volume corrected and basil cbf   in template space
    qc_file
        quality control meausres 

    See Also
    --------

    * :py:func:`~aslprep.niworkflows.func.util.init_asl_reference_wf`
    * :py:func:`~aslprep.workflows.asl.stc.init_asl_stc_wf`
    * :py:func:`~aslprep.workflows.asl.hmc.init_asl_hmc_wf`
    * :py:func:`~aslprep.workflows.asl.t2s.init_asl_t2s_wf`
    * :py:func:`~aslprep.workflows.asl.registration.init_asl_t1_trans_wf`
    * :py:func:`~aslprep.workflows.asl.registration.init_asl_reg_wf`
    * :py:func:`~aslprep.workflows.asl.confounds.init_asl_confounds_wf`
    * :py:func:`~aslprep.workflows.asl.confounds.init_ica_aroma_wf`
    * :py:func:`~aslprep.workflows.asl.resampling.init_asl_std_trans_wf`
    * :py:func:`~aslprep.workflows.asl.resampling.init_asl_preproc_trans_wf`
    * :py:func:`~aslprep.workflows.asl.resampling.init_asl_surf_wf`
    * :py:func:`~sdcflows.workflows.fmap.init_fmap_wf`
    * :py:func:`~sdcflows.workflows.pepolar.init_pepolar_unwarp_wf`
    * :py:func:`~sdcflows.workflows.phdiff.init_phdiff_wf`
    * :py:func:`~sdcflows.workflows.syn.init_syn_sdc_wf`
    * :py:func:`~sdcflows.workflows.unwarp.init_sdc_unwarp_wf`

    """
    from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from ...niworkflows.func.util import init_asl_reference_wf
    from ...niworkflows.interfaces.nibabel import ApplyMask
    from ...niworkflows.interfaces.utility import KeySelect
    from sdcflows.workflows.base import init_sdc_estimate_wf, fieldmap_wrangler

    ref_file = asl_file
    mem_gb = {'filesize': 1, 'resampled': 1, 'largemem': 1}
    asl_tlen = 10

    # Have some options handy
    layout = config.execution.layout
    omp_nthreads = config.nipype.omp_nthreads
    spaces = config.workflow.spaces
    output_dir = str(config.execution.output_dir)
    mscale = config.workflow.m0_scale


    if os.path.isfile(ref_file):
        asl_tlen, mem_gb = _create_mem_gb(ref_file)

    wf_name = _get_wf_name(ref_file)
    config.loggers.workflow.debug(
        'Creating asl processing workflow for "%s" (%.2f GB / %d TRs). '
        'Memory resampled/largemem=%.2f/%.2f GB.',
        ref_file, mem_gb['filesize'], asl_tlen, mem_gb['resampled'], mem_gb['largemem'])

    
    # Find associated sbref, if possible
    refbase = os.path.basename(ref_file)
    config.loggers.workflow.info("No single-band-reference found for %s.",
                                     refbase)
    metadata = layout.get_metadata(ref_file)
    # Build workflow
    workflow = Workflow(name=wf_name)
    workflow.__postdesc__ = """\
All resamplings can be performed with *a single interpolation
step* by composing all the pertinent transformations (i.e. head-motion
transform matrices, susceptibility distortion correction when available,
and co-registrations to anatomical and output spaces).
Gridded (volumetric) resamplings were performed using `antsApplyTransforms` (ANTs),
configured with Lanczos interpolation to minimize the smoothing
effects of other kernels [@lanczos].
Non-gridded (surface) resamplings were performed using `mri_vol2surf`
(FreeSurfer).
"""

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['asl_file','t1w_preproc', 't1w_mask', 't1w_dseg', 't1w_tpms',
                'anat2std_xfm', 'std2anat_xfm', 'template']),
        name='inputnode')
    inputnode.inputs.asl_file = asl_file
    subj_dir=str(config.execution.bids_dir) + '/sub-' + str(config.execution.participant_label[0])
    

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['asl_t1', 'asl_t1_ref', 'asl_mask_t1','asl_std', 'asl_std_ref', 'asl_mask_std',
                'asl_native','cbf_t1', 'cbf_std', 'meancbf_t1', 'meancbf_std', 'score_t1', 'score_std',
                'avgscore_t1', 'avgscore_std', 'avgscore_cifti', ' scrub_t1', 'scrub_std',
                'basil_t1', 'basil_std', 'pv_t1', 'pv_std', 'pv_native','att','att_t1','att_std',
                'qc_file']),
        name='outputnode')
      
    # Generate a brain-masked conversion of the t1w
    t1w_brain = pe.Node(ApplyMask(), name='t1w_brain')

    # asl buffer: an identity used as a pointer to either the original asl
    # or the STC'ed one for further use.
    aslbuffer = pe.Node(niu.IdentityInterface(fields=['asl_file']), name='aslbuffer')
    aslbuffer.inputs.asl_file=ref_file
    summary = pe.Node(
        FunctionalSummary(
            registration=('FSL'),
            registration_dof=config.workflow.asl2t1w_dof,
            registration_init=config.workflow.asl2t1w_init,
            pe_direction=metadata.get("PhaseEncodingDirection"),
            tr=metadata.get("RepetitionTime")),
        name='summary', mem_gb=config.DEFAULT_MEMORY_MIN_GB, run_without_submitting=True)
    summary.inputs.dummy_scans = 0

    asl_derivatives_wf = init_asl_derivatives_wf(
        bids_root=layout.root,
        cifti_output=None,
        freesurfer=None,
        metadata=metadata,
        output_dir=output_dir,
        spaces=spaces,
    )

    workflow.connect([
        (outputnode, asl_derivatives_wf, [
            ('asl_t1', 'inputnode.asl_t1'),
            ('asl_t1_ref', 'inputnode.asl_t1_ref'),
            ('asl_mask_t1', 'inputnode.asl_mask_t1'),
            ('asl_native', 'inputnode.asl_native'),
        ]),
    ])
     
     # begin workflow 
    gen_ref_wf = init_asl_geref_wf(omp_nthreads=omp_nthreads,mem_gb=mem_gb['filesize'],
                              metadata=metadata,bids_dir=subj_dir,brainmask_thresh=0.5,
                              pre_mask=False, name="asl_gereference_wf",gen_report=False)
    
    reg_ge_wf = init_asl_gereg_wf(use_bbr=config,
                asl2t1w_dof=config.workflow.asl2t1w_dof,
                asl2t1w_init=config.workflow.asl2t1w_init,
                mem_gb=2, omp_nthreads=omp_nthreads, name='asl_reg_wf',
                sloppy=False, use_compression=True, write_report=True)

    t1w_gereg_wf = init_asl_t1_getrans_wf(mem_gb=3, omp_nthreads=omp_nthreads, cbft1space=True,
                          use_compression=True, name='asl_t1_trans_wf')

    std_gereg_wf = init_asl_gestd_trans_wf(mem_gb=4,omp_nthreads=omp_nthreads,spaces=spaces,
                            name='asl_gestd_trans_wf', use_compression=True)

    cbf_compt_wf = init_gecbf_compt_wf(mem_gb=mem_gb, metadata=3,bids_dir=subj_dir,omp_nthreads=omp_nthreads,
                                M0Scale=1,smooth_kernel=5,name='cbf_compt_wf')
    
    
    workflow.connect([
         (inputnode, gen_ref_wf,[('asl_file','inputnode.asl_file')]),
         (gen_ref_wf, reg_ge_wf,[('ref_image_brain','ref_asl_brain')]),
         (inputnode, t1w_brain,[('t1w_preproc', 'in_file'),
                                ('t1w_mask', 'in_mask')]),
         (t1w_brain,reg_ge_wf,[('out_file','inputnode.t1w_brain')]),
         (inputnode,reg_ge_wf,[('t1w_dseg','inputnode.t1w_dseg')]),
    
     ])
    
    

    

def _get_series_len(asl_fname):
    from ...niworkflows.interfaces.registration import _get_vols_to_discard
    img = nb.load(asl_fname)
    if len(img.shape) < 4:
        return 1

    skip_vols = _get_vols_to_discard(img)

    return img.shape[3] - skip_vols


def _create_mem_gb(asl_fname):
    asl_size_gb = os.path.getsize(asl_fname) / (1024**3)
    asl_tlen = nb.load(asl_fname).shape[-1]
    mem_gb = {
        'filesize': asl_size_gb,
        'resampled': asl_size_gb * 4,
        'largemem': asl_size_gb * (max(asl_tlen / 100, 1.0) + 4),
    }

    return asl_tlen, mem_gb


def _get_wf_name(asl_fname):
    """
    Derive the workflow name for supplied asl file.

    >>> _get_wf_name('/completely/made/up/path/sub-01_task-nback_asl.nii.gz')
    'func_preproc_task_nback_wf'
    >>> _get_wf_name('/completely/made/up/path/sub-01_task-nback_run-01_echo-1_asl.nii.gz')
    'func_preproc_task_nback_run_01_echo_1_wf'

    """
    from nipype.utils.filemanip import split_filename
    fname = split_filename(asl_fname)[1]
    fname_nosub = '_'.join(fname.split("_")[1:])
    # if 'echo' in fname_nosub:
    #     fname_nosub = '_'.join(fname_nosub.split("_echo-")[:1]) + "_asl"
    name = "asl_preproc_" + fname_nosub.replace(
        ".", "_").replace(" ", "").replace("-", "_").replace("_asl", "_wf")

    return name


def _to_join(in_file, join_file):
    """Join two tsv files if the join_file is not ``None``."""
    from ...niworkflows.interfaces.utils import JoinTSVColumns
    if join_file is None:
        return in_file
    res = JoinTSVColumns(in_file=in_file, join_file=join_file).run()
    return res.outputs.out_file


def check_img(img):
    # get the 4th dimension
    import numpy as np 
    ss=nb.load(img).get_fdata().shape
    if len(ss) == 3:
        ss=np.hstack([ss,0])
    return ss[3]

def _getTR(img):
    import nibabel as nib
    import numpy as np
    asl_img = nib.load(img)
    tr = asl_img.header.get_zooms()[-1]
    return np.float(tr)