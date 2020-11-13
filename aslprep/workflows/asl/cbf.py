# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
from ...niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from ...interfaces.cbf_computation import (extractCBF, computeCBF, scorescrubCBF, BASILCBF,
                                           refinemask, qccbf, cbfqroiquant,qccbfge,extractCB)
from ...niworkflows.interfaces.plotting import (CBFSummary, CBFtsSummary)
from ...interfaces import DerivativesDataSink
import numpy as np
import pandas as pd
import os
import tempfile
from ...config import DEFAULT_MEMORY_MIN_GB


def init_cbf_compt_wf(mem_gb, metadata,bids_dir,dummy_vols, omp_nthreads,
                     scorescrub=False,basil=False,M0Scale=1,smooth_kernel=5,
                      name='cbf_compt_wf'):
    """
    Create a workflow for :abbr:`CCBF ( compute cbf)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.cbf import init_cbf_compt_wf
            wf = init_cbf_compt_wf(mem_gb=0.1,smooth_kernel=5,dummy_vols=0)

    Parameters
    ----------
    metadata : :obj:`dict`
        BIDS metadata for asl file
    name : :obj:`str`
        Name of workflow (default: ``cbf_compt_wf``)

    Inputs
    ------
    in_file
        raw asl file 
    bids_dir
        bids dir of the subject
    asl_file
        asl series NIfTI file
    asl_mask
        asl mask NIFTI file 
    t1w_tpms
        t1w probability maps 
    t1w_mask
        t1w mask Nifti
    t1_asl_xform
        t1w to asl transfromation file 
    itk_asl_to_t1
        asl to t1q transfromation file
    scorescrub
        compute score and scrub 
    basil
        compute basil and PVC

    Outputs
    -------
    *cbf
       all cbf outputs
       cbf,score, scrub, pv, and basil

    """

                    
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The CBF was quantified from  *preproccessed* ASL data using a relatively basic
model [@detre_perfusion] [@alsop_recommended].
"""
    inputnode = pe.Node(niu.IdentityInterface(fields=['asl_file', 'in_file', 'asl_mask',
                                                      't1w_tpms', 't1w_mask', 't1_asl_xform',
                                                      'itk_asl_to_t1']),
                        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_cbf', 'out_mean', 'out_score',
                                                       'out_avgscore', 'out_scrub', 'out_cbfb',
                                                       'out_scoreindex', 'out_cbfpv','out_att']),
                         name='outputnode')
    # convert tmps to asl_space
    csf_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                      name='csf_tfm', mem_gb=0.1)
    wm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='wm_tfm', mem_gb=0.1)
    gm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='gm_tfm', mem_gb=0.1)
    
    

    tiscbf=np.add(metadata["PostLabelingDelay"],metadata["LabelingDuration"])
    if hasattr(tiscbf, '__len__'):
        tisasl = ",".join([str(i) for i in tiscbf])
    else:
        tisasl = str(tiscbf)

    def pcaslorasl(metadata):
        if 'CASL' in metadata["LabelingType"]:
            pcasl1 = True
        elif 'PASL' in metadata["LabelingType"]:
            pcasl1 = False
        return pcasl1
  
    extractcbf = pe.Node(extractCBF(dummy_vols=dummy_vols,bids_dir=bids_dir,
                  fwhm=smooth_kernel,in_metadata=metadata), mem_gb=0.2,
                         run_without_submitting=True, name="extractcbf")
    computecbf = pe.Node(computeCBF(in_metadata=metadata,in_m0scale=M0Scale), mem_gb=0.2,
                         run_without_submitting=True, name="computecbf")
    scorescrub = pe.Node(scorescrubCBF(in_thresh=0.7, in_wfun='huber'), mem_gb=0.2,
                         name='scorescrub', run_without_submitting=True)
    basilcbf = pe.Node(BASILCBF(m0scale=M0Scale, bolus=metadata["LabelingDuration"],
                                m0tr=metadata['RepetitionTime'], pvc=True,tis = tisasl,
                                pcasl = pcaslorasl(metadata = metadata)), name='basilcbf',
                       run_without_submitting=True, mem_gb=0.2)

    refinemaskj = pe.Node(refinemask(), mem_gb=0.2, run_without_submitting=True, name="refinemask")

    def _pick_csf(files):
        return files[-1]

    def _pick_gm(files):
        return files[0]

    def _pick_wm(files):
        return files[1]

    def _getfiledir(file):
        import os
        return os.path.dirname(file)

    workflow.connect([
        # extract CBF data and compute cbf
        (inputnode,  extractcbf, [('in_file', 'in_file'), ('asl_file', 'asl_file')]),
        (extractcbf, computecbf, [('out_file', 'in_cbf'), ('out_avg', 'in_m0file')]),
        # (inputnode,computecbf,[('asl_mask','in_mask')]),
        (inputnode, refinemaskj, [('t1w_mask', 'in_t1mask'), ('asl_mask', 'in_aslmask'),
                                  ('t1_asl_xform', 'transforms')]),
        (refinemaskj, extractcbf, [('out_mask', 'in_mask')]),
        (refinemaskj, computecbf, [('out_mask', 'in_mask')]),
        (computecbf, outputnode, [('out_cbf', 'out_cbf'),
                                  ('out_mean', 'out_mean')]),
        (inputnode, csf_tfm, [('asl_mask', 'reference_image'),
                              ('t1_asl_xform', 'transforms')]),
        (inputnode, csf_tfm, [(('t1w_tpms', _pick_csf), 'input_image')]),
        (inputnode, wm_tfm, [('asl_mask', 'reference_image'),
                             ('t1_asl_xform', 'transforms')]),
        (inputnode, wm_tfm, [(('t1w_tpms', _pick_wm), 'input_image')]),
        (inputnode, gm_tfm, [('asl_mask', 'reference_image'),
                             ('t1_asl_xform', 'transforms')]),
        (inputnode, gm_tfm, [(('t1w_tpms', _pick_gm), 'input_image')]),
    ])

    if scorescrub:
        workflow.__desc__ = workflow.__desc__ +  """\
CBF are susceptible to artifactsdue to low signal to noise ratio  and  sensitivity to  motion, Structural Correlation
based Outlier Rejection (SCORE) algothim was applied to the CBF to discard few extreme outliers [@score_dolui]. 
Furthermore,Structural Correlation with RobUst Bayesian (SCRUB) algorithms was applied to the CBF by iteratively 
reweighted  CBF  with structural tissues probalility maps [@scrub_dolui].
"""      
        workflow.connect([
        (refinemaskj, scorescrub, [('out_mask', 'in_mask')]),
        (computecbf, scorescrub, [('out_cbf', 'in_file')]),
        (gm_tfm, scorescrub, [('output_image', 'in_greyM')]),
        (wm_tfm, scorescrub, [('output_image', 'in_whiteM')]),
        (csf_tfm, scorescrub, [('output_image', 'in_csf')]),
        (scorescrub, outputnode, [('out_score', 'out_score'), 
                                  ('out_scoreindex', 'out_scoreindex'),
                                  ('out_avgscore', 'out_avgscore'), 
                                  ('out_scrub', 'out_scrub')]),
        ])
    if basil:
        workflow.__desc__ = workflow.__desc__ +  """\
Alternate method of CBF computation is Bayesian Inference
for Arterial Spin Labeling (BASIL) as implmented in FSL which is  based on Bayeisan inference
principles [@chappell_basil]. BASIL computed the CBF from ASL incoporating natural varaibility
of other model parameters and spatial regularization of the estimated perfusion image. BASIL
also included correction for partial volume effects [@chappell_pvc].
"""
        workflow.connect([
         (refinemaskj, basilcbf, [('out_mask', 'mask')]),
         (extractcbf, basilcbf, [(('out_avg', _getfiledir), 'out_basename')]),
         (extractcbf, basilcbf, [('out_file', 'in_file')]),
         (gm_tfm, basilcbf, [('output_image', 'pvgm')]),
         (wm_tfm, basilcbf, [('output_image', 'pvwm')]),
         (extractcbf, basilcbf, [('out_avg', 'mzero')]),
         (basilcbf, outputnode, [('out_cbfb', 'out_cbfb'),
                                ('out_cbfpv', 'out_cbfpv'),
                                ('out_att', 'out_att')]),
        ])
            
    return workflow




def init_cbfqc_compt_wf(mem_gb,asl_file, metadata,omp_nthreads,scorescrub=False,basil=False,name='cbfqc_compt_wf'):
    """
    Create a workflow for :abbr:`cbfqc( compute cbf)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.cbf import init_cbfqc_compt_wf
            wf = init_cbfqc_compt_wf(mem_gb=0.1)

    Parameters
    ----------
    metadata : :obj:`dict`
        BIDS metadata for asl file
    name : :obj:`str`
        Name of workflow (default: ``cbfqc_compt_wf'``)

    Inputs
    ------
    *cbf
        all cbf 
    asl_mask
        asl mask NIFTI file 
    t1w_tpms
        t1w probability maps 
    t1_asl_xform
        t1w to asl transfromation file 

    Outputs
    -------
    qc_file
       qc measures in tsv

    """

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The following quality control (qc) measures was estimated: framewise displacement and relative
root mean square dice index. Other qc meaure include dice and jaccard indices, cross-correlation
and coverage that estimate the coregistration quality of  ASL and T1W images and  normalization
quality of ASL to template. Quality evaluation index (QEI) was also computed for CBF [@cbfqc].
The  QEI is  automated for objective quality evaluation of CBF maps and measured the CBF quality
based on structural similarity,spatial variability and the percentatge  of voxels with  negtaive
CBF within Grey matter
"""
    inputnode = pe.Node(niu.IdentityInterface(fields=['meancbf', 'avgscore', 'scrub', 'basil',
                                                      'asl_mask', 't1w_tpms',  'confmat',
                                                      'asl_mask_std', 't1_asl_xform', 'pv',
                                                      't1w_mask']),
                        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['qc_file']), name='outputnode')

    def _pick_csf(files):
        return files[-1]

    def _pick_gm(files):
        return files[0]

    def _pick_wm(files):
        return files[1]

    csf_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                      name='csf_tfm', mem_gb=0.1)
    wm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='wm_tfm', mem_gb=0.1)
    gm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='gm_tfm', mem_gb=0.1)

    mask_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                       name='masktonative', mem_gb=0.1)

    from templateflow.api import get as get_template
    brain_mask = str(get_template('MNI152NLin2009cAsym', resolution=2, desc='brain',
                                  suffix='mask'))

    from nipype.interfaces.afni import Resample
    resample = pe.Node(Resample(in_file=brain_mask, outputtype='NIFTI_GZ'),
                       name='resample', mem_gb=0.1)

    qccompute = pe.Node(qccbf(in_file=asl_file), name='qccompute',
                        run_without_submitting=True, mem_gb=0.2)

    workflow.connect([(inputnode, csf_tfm, [('asl_mask', 'reference_image'),
                                            ('t1_asl_xform', 'transforms')]),
                      (inputnode, csf_tfm, [(('t1w_tpms', _pick_csf), 'input_image')]),
                      (inputnode, wm_tfm, [('asl_mask', 'reference_image'),
                                           ('t1_asl_xform', 'transforms')]),
                      (inputnode, wm_tfm, [(('t1w_tpms', _pick_wm), 'input_image')]),
                      (inputnode, gm_tfm, [('asl_mask', 'reference_image'),
                                           ('t1_asl_xform', 'transforms')]),
                      (inputnode, gm_tfm, [(('t1w_tpms', _pick_gm), 'input_image')]),
                      (inputnode, mask_tfm, [('asl_mask', 'reference_image'),
                                             ('t1_asl_xform', 'transforms'),
                                             ('t1w_mask', 'input_image')]),
                      (mask_tfm, qccompute, [('output_image', 'in_t1mask')]),
                      (inputnode, qccompute, [('asl_mask', 'in_aslmask'),
                                              ('confmat', 'in_confmat')]),
                      (inputnode, qccompute, [(('asl_mask_std', _pick_csf), 'in_aslmaskstd')]),
                      (inputnode, resample, [(('asl_mask_std', _pick_csf), 'master')]),
                      (resample, qccompute, [('out_file', 'in_templatemask')]),
                      (gm_tfm, qccompute, [('output_image', 'in_greyM')]),
                      (wm_tfm, qccompute, [('output_image', 'in_whiteM')]),
                      (csf_tfm, qccompute, [('output_image', 'in_csf')]),
                      (inputnode, qccompute, [('meancbf', 'in_meancbf')],)
       ])


    if scorescrub:
        workflow.connect([ 
                      (inputnode, qccompute, [('scrub', 'in_scrub'),
                                              ('avgscore', 'in_avgscore'),]),
        ]) 
     
    if basil: 
        workflow.connect([ 
                      (inputnode, qccompute, [('basil', 'in_basil'), ('pv', 'in_pvc')]),
        ]),
                                              
    workflow.connect([(qccompute, outputnode, [('qc_file', 'qc_file')]),
                      ])

    return workflow

def init_cbfgeqc_compt_wf(mem_gb, asl_file, metadata, omp_nthreads,scorescrub=False,basil=False,name='cbfqc_compt_wf'):
    """
    Create a workflow for :abbr:`cbfqc( compute cbf)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.cbf import init_cbfqc_compt_wf
            wf = init_cbfqc_compt_wf(mem_gb=0.1)

    Parameters
    ----------
    metadata : :obj:`dict`
        BIDS metadata for asl file
    name : :obj:`str`
        Name of workflow (default: ``cbfqc_compt_wf'``)

    Inputs
    ------
    *cbf
        all cbf 
    asl_mask
        asl mask NIFTI file 
    t1w_tpms
        t1w probability maps 
    t1_asl_xform
        t1w to asl transfromation file 

    Outputs
    -------
    qc_file
       qc measures in tsv

    """

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The following quality control (qc) measures was estimated: dice and jaccard indices, cross-correlation
and coverage that estimate the coregistration quality of  ASL and T1W images and  normalization
quality of ASL to template. Quality evaluation index (QEI) was also computed for CBF [@cbfqc].
The  QEI is  automated for objective quality evaluation of CBF maps and measured the CBF quality
based on structural similarity,spatial variability and the percentatge  of voxels with  negtaive
CBF within Grey matter
"""
    inputnode = pe.Node(niu.IdentityInterface(fields=['meancbf', 'avgscore', 'scrub', 'basil',
                                                      'asl_mask', 't1w_tpms',
                                                      'asl_mask_std', 't1_asl_xform', 'pv',
                                                      't1w_mask']),
                        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['qc_file']), name='outputnode')

    def _pick_csf(files):
        return files[-1]

    def _pick_gm(files):
        return files[0]

    def _pick_wm(files):
        return files[1]

    csf_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                      name='csf_tfm', mem_gb=0.1)
    wm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='wm_tfm', mem_gb=0.1)
    gm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='gm_tfm', mem_gb=0.1)

    mask_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                       name='masktonative', mem_gb=0.1)

    from templateflow.api import get as get_template
    brain_mask = str(get_template('MNI152NLin2009cAsym', resolution=2, desc='brain',
                                  suffix='mask'))

    from nipype.interfaces.afni import Resample
    resample = pe.Node(Resample(in_file=brain_mask, outputtype='NIFTI_GZ'),
                       name='resample', mem_gb=0.1)

    qccompute = pe.Node(qccbfge(in_file=asl_file), name='qccompute',
                        run_without_submitting=True, mem_gb=0.2)

    workflow.connect([
                    (inputnode, csf_tfm, [('asl_mask', 'reference_image'),
                                            ('t1_asl_xform', 'transforms')]),
                      (inputnode, csf_tfm, [(('t1w_tpms', _pick_csf), 'input_image')]),
                      (inputnode, wm_tfm, [('asl_mask', 'reference_image'),
                                           ('t1_asl_xform', 'transforms')]),
                      (inputnode, wm_tfm, [(('t1w_tpms', _pick_wm), 'input_image')]),
                      (inputnode, gm_tfm, [('asl_mask', 'reference_image'),
                                           ('t1_asl_xform', 'transforms')]),
                      (inputnode, gm_tfm, [(('t1w_tpms', _pick_gm), 'input_image')]),
                      (inputnode, mask_tfm, [('asl_mask', 'reference_image'),
                                             ('t1_asl_xform', 'transforms'),
                                             ('t1w_mask', 'input_image')]),
                      (mask_tfm, qccompute, [('output_image', 'in_t1mask')]),
                      (inputnode, qccompute, [('asl_mask', 'in_aslmask'),]),
                      (inputnode, qccompute, [(('asl_mask_std', _pick_csf), 'in_aslmaskstd')]),
                      (inputnode, resample, [(('asl_mask_std', _pick_csf), 'master')]),
                      (resample, qccompute, [('out_file', 'in_templatemask')]),
                      (gm_tfm, qccompute, [('output_image', 'in_greyM')]),
                      (wm_tfm, qccompute, [('output_image', 'in_whiteM')]),
                      (csf_tfm, qccompute, [('output_image', 'in_csf')]),
                      (inputnode, qccompute, [('meancbf', 'in_meancbf')]),
                     
                      ])
    if scorescrub:
        workflow.connect([ 
                      (inputnode, qccompute, [('scrub', 'in_scrub'),
                                              ('avgscore', 'in_avgscore'),]),
        ]) 
     
    if basil: 
        workflow.connect([ 
                      (inputnode, qccompute, [('basil', 'in_basil'), ('pv', 'in_pvc')]),
        ]),
                                              
    workflow.connect([
          (qccompute, outputnode, [('qc_file', 'qc_file')]),
                      ])
    return workflow



def init_cbfplot_wf(mem_gb, metadata, omp_nthreads, scorescrub=False,basil=False,name='cbf_plot'):
    workflow = Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=['cbf', 'cbf_ts', 'score_ts', 'score',
                                                      'scrub', 'asl_ref', 'basil', 'pvc',
                                                      'asl_mask', 't1_asl_xform', 'std2anat_xfm',
                                                      'confounds_file', 'scoreindex']),
                        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['cbf_carpetplot', 'score_carpetplot',
                                                       'cbf_summary_plot', 'cbf_summary_plot',
                                                       'score_summary_plot', 'scrub_summary_plot',
                                                       'basil_summary_plot', 'pvc_summary_plot']),
                         name='outputnode')
    mrg_xfms = pe.Node(niu.Merge(2), name='mrg_xfms')
    from templateflow.api import get as get_template
    seg = get_template(
            'MNI152NLin2009cAsym', resolution=1, desc='carpet',
            suffix='dseg')
    print(seg)
    resample_parc = pe.Node(ApplyTransforms(
        float=True,
        input_image=str(seg),
        dimension=3, default_value=0, interpolation='MultiLabel'),
        name='resample_parc')

    cbftssummary = pe.Node(CBFtsSummary(tr=metadata['RepetitionTime']),
                           name='cbf_ts_summary', mem_gb=2)
    cbfsummary = pe.Node(CBFSummary(label='cbf'), name='cbf_summary', mem_gb=1)
    ds_report_cbftsplot = pe.Node(
        DerivativesDataSink(desc='cbftsplot', datatype="figures",  keep_dtype=True),
        name='ds_report_cbftsplot', run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)
    ds_report_cbfplot = pe.Node(
        DerivativesDataSink(desc='cbfplot', datatype="figures", keep_dtype=True),
        name='ds_report_cbfplot', run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)
    workflow.connect(
        [(inputnode, mrg_xfms, [('t1_asl_xform', 'in1'), ('std2anat_xfm', 'in2')]),
                      (inputnode, resample_parc, [('asl_mask', 'reference_image')]),
                      (mrg_xfms, resample_parc, [('out', 'transforms')]),
                      (resample_parc, cbftssummary, [('output_image', 'seg_file')]),
                      (inputnode, cbftssummary, [('cbf_ts', 'cbf_ts'),
                                                 ('confounds_file', 'conf_file'),
                                                 ('scoreindex', 'score_file')]),
                      (cbftssummary, ds_report_cbftsplot, [('out_file', 'in_file')]),
                      (cbftssummary, outputnode, [('out_file', 'cbf_carpetplot')]),
                      (inputnode, cbfsummary, [('cbf', 'cbf'), ('asl_ref', 'ref_vol')]),
                      (cbfsummary, ds_report_cbfplot, [('out_file', 'in_file')]),
                      (cbfsummary, outputnode, [('out_file', 'cbf_summary_plot')]),
       ])

    if scorescrub:
        scoresummary = pe.Node(CBFSummary(label='score'), name='score_summary', mem_gb=1)
        scrubsummary = pe.Node(CBFSummary(label='scrub'), name='scrub_summary', mem_gb=1)
        ds_report_scoreplot = pe.Node(
           DerivativesDataSink(desc='scoreplot', datatype="figures",  keep_dtype=True),
          name='ds_report_scoreplot', run_without_submitting=True,
           mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_report_scrubplot = pe.Node(
           DerivativesDataSink(desc='scrubplot', datatype="figures",  keep_dtype=True),
           name='ds_report_scrubplot', run_without_submitting=True,
           mem_gb=DEFAULT_MEMORY_MIN_GB)
        workflow.connect([
            (inputnode, scoresummary, [('score', 'cbf'), ('asl_ref', 'ref_vol')]),
            (scoresummary, ds_report_scoreplot, [('out_file', 'in_file')]),
            (scoresummary, outputnode, [('out_file', 'score_summary_plot')]),
            (inputnode, scrubsummary, [('scrub', 'cbf'), ('asl_ref', 'ref_vol')]),
            (scrubsummary, ds_report_scrubplot, [('out_file', 'in_file')]),
            (scrubsummary, outputnode, [('out_file', 'scrub_summary_plot')]),

        ])

    if basil:
        basilsummary = pe.Node(CBFSummary(label='basil'), name='basil_summary', mem_gb=1)
        pvcsummary = pe.Node(CBFSummary(label='pvc'), name='pvc_summary', mem_gb=1)
        ds_report_basilplot = pe.Node(
           DerivativesDataSink(desc='basilplot', datatype="figures",  keep_dtype=True),
            name='ds_report_basilplot', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_report_pvcplot = pe.Node(
            DerivativesDataSink(desc='pvcplot', datatype="figures", keep_dtype=True),
            name='ds_report_pvcplot', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
    
        workflow.connect([
                      (inputnode, basilsummary, [('basil', 'cbf'), ('asl_ref', 'ref_vol')]),
                      (basilsummary, ds_report_basilplot, [('out_file', 'in_file')]),
                      (basilsummary, outputnode, [('out_file', 'basil_summary_plot')]),
                      (inputnode, pvcsummary, [('pvc', 'cbf'), ('asl_ref', 'ref_vol')]),
                      (pvcsummary, ds_report_pvcplot, [('out_file', 'in_file')]),
                      (pvcsummary, outputnode, [('out_file', 'pvc_summary_plot')]),
                      ])
    return workflow

def init_gecbfplot_wf(mem_gb, metadata, omp_nthreads,scorescrub=False,basil=False, name='cbf_plot'):
    workflow = Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=['cbf', 'score',
                                                      'scrub', 'asl_ref', 'basil', 'pvc'
                                                      ]),
                        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=[
                                                       'cbf_summary_plot', 
                                                       'score_summary_plot', 'scrub_summary_plot',
                                                       'basil_summary_plot', 'pvc_summary_plot']),
                         name='outputnode')
    mrg_xfms = pe.Node(niu.Merge(2), name='mrg_xfms')
    from templateflow.api import get as get_template
    seg = get_template(
            'MNI152NLin2009cAsym', resolution=1, desc='carpet',
            suffix='dseg')
    print(seg)
    resample_parc = pe.Node(ApplyTransforms(
        float=True,
        input_image=str(seg),
        dimension=3, default_value=0, interpolation='MultiLabel'),
        name='resample_parc')



    cbfsummary = pe.Node(CBFSummary(label='cbf'), name='cbf_summary', mem_gb=1)
    ds_report_cbfplot = pe.Node(
        DerivativesDataSink(desc='cbfplot', datatype="figures", keep_dtype=True),
        name='ds_report_cbfplot', run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)
    workflow.connect([
                      (inputnode, cbfsummary, [('cbf', 'cbf'), ('asl_ref', 'ref_vol')]),
                      (cbfsummary, ds_report_cbfplot, [('out_file', 'in_file')]),
                      (cbfsummary, outputnode, [('out_file', 'cbf_summary_plot')]),
       ])

    if scorescrub:
        scoresummary = pe.Node(CBFSummary(label='score'), name='score_summary', mem_gb=1)
        scrubsummary = pe.Node(CBFSummary(label='scrub'), name='scrub_summary', mem_gb=1)
        ds_report_scoreplot = pe.Node(
           DerivativesDataSink(desc='scoreplot', datatype="figures",  keep_dtype=True),
          name='ds_report_scoreplot', run_without_submitting=True,
           mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_report_scrubplot = pe.Node(
           DerivativesDataSink(desc='scrubplot', datatype="figures",  keep_dtype=True),
           name='ds_report_scrubplot', run_without_submitting=True,
           mem_gb=DEFAULT_MEMORY_MIN_GB)
        workflow.connect([
            (inputnode, scoresummary, [('score', 'cbf'), ('asl_ref', 'ref_vol')]),
            (scoresummary, ds_report_scoreplot, [('out_file', 'in_file')]),
            (scoresummary, outputnode, [('out_file', 'score_summary_plot')]),
            (inputnode, scrubsummary, [('scrub', 'cbf'), ('asl_ref', 'ref_vol')]),
            (scrubsummary, ds_report_scrubplot, [('out_file', 'in_file')]),
            (scrubsummary, outputnode, [('out_file', 'scrub_summary_plot')]),

        ])

    if basil:
        basilsummary = pe.Node(CBFSummary(label='basil'), name='basil_summary', mem_gb=1)
        pvcsummary = pe.Node(CBFSummary(label='pvc'), name='pvc_summary', mem_gb=1)
        ds_report_basilplot = pe.Node(
           DerivativesDataSink(desc='basilplot', datatype="figures",  keep_dtype=True),
            name='ds_report_basilplot', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
        ds_report_pvcplot = pe.Node(
            DerivativesDataSink(desc='pvcplot', datatype="figures", keep_dtype=True),
            name='ds_report_pvcplot', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)
    
        workflow.connect([
                      (inputnode, basilsummary, [('basil', 'cbf'), ('asl_ref', 'ref_vol')]),
                      (basilsummary, ds_report_basilplot, [('out_file', 'in_file')]),
                      (basilsummary, outputnode, [('out_file', 'basil_summary_plot')]),
                      (inputnode, pvcsummary, [('pvc', 'cbf'), ('asl_ref', 'ref_vol')]),
                      (pvcsummary, ds_report_pvcplot, [('out_file', 'in_file')]),
                      (pvcsummary, outputnode, [('out_file', 'pvc_summary_plot')]),
                      ])
    return workflow


def init_cbfroiquant_wf(mem_gb, omp_nthreads,scorescrub=False,basil=False,name='cbf_roiquant'):
    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['cbf', 'score', 'scrub', 'basil', 'pvc',
                                                      'aslmask', 't1_asl_xform',
                                                      'std2anat_xfm']),
                        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['cbf_hvoxf', 'score_hvoxf', 'scrub_hvoxf', 'basil_hvoxf', 'pvc_hvoxf',
                'cbf_sc207', 'score_sc207', 'scrub_sc207', 'basil_sc207', 'pvc_sc207',
                'cbf_sc217', 'score_sc217', 'scrub_sc217', 'basil_sc217', 'pvc_sc217',
                'cbf_sc407', 'score_sc407', 'scrub_sc407', 'basil_sc407', 'pvc_sc407',
                'cbf_sc417', 'score_sc417', 'scrub_sc417', 'basil_sc417', 'pvc_sc417']),
        name='outputnode')

    from ...interfaces.cbf_computation import get_atlas
    hvoxfile, hvoxdata, hvoxlabel = get_atlas(atlasname='HarvardOxford')
    sc207file, sc207data, sc207label = get_atlas(atlasname='schaefer200x7')
    sc217file, sc217data, sc217label = get_atlas(atlasname='schaefer200x17')
    sc407file, sc407data, sc407label = get_atlas(atlasname='schaefer400x7')
    sc417file, sc417data, sc417label = get_atlas(atlasname='schaefer400x17')

    mrg_xfms = pe.Node(niu.Merge(2), name='mrg_xfms')
    hvoftrans = pe.Node(ApplyTransforms(float=True, input_image=hvoxfile,
                        dimension=3, default_value=0, interpolation='NearestNeighbor'),
                        name='hvoftrans')
    sc207trans = pe.Node(ApplyTransforms(float=True, input_image=sc207file,
                         dimension=3, default_value=0, interpolation='NearestNeighbor'),
                         name='sc207trans')
    sc217trans = pe.Node(ApplyTransforms(float=True, input_image=sc217file,
                         dimension=3, default_value=0, interpolation='NearestNeighbor'),
                         name='sc217trans')
    sc407trans = pe.Node(ApplyTransforms(float=True, input_image=sc407file,
                         dimension=3, default_value=0, interpolation='NearestNeighbor'),
                         name='sc407trans')
    sc417trans = pe.Node(ApplyTransforms(float=True, input_image=sc417file,
                         dimension=3, default_value=0, interpolation='NearestNeighbor'),
                         name='sc417trans')

    cbfroihv = pe.Node(cbfqroiquant(atlaslabel=hvoxlabel, atlasdata=hvoxdata), name='cbfroihv')
    cbfroi207 = pe.Node(cbfqroiquant(atlaslabel=sc207label, atlasdata=sc207data), name='cbf207')
    cbfroi217 = pe.Node(cbfqroiquant(atlaslabel=sc217label, atlasdata=sc217data), name='cbf217')
    cbfroi407 = pe.Node(cbfqroiquant(atlaslabel=sc407label, atlasdata=sc407data), name='cbf407')
    cbfroi417 = pe.Node(cbfqroiquant(atlaslabel=sc417label, atlasdata=sc417data), name='cbf417')
    if scorescrub:
        scorehv = pe.Node(cbfqroiquant(atlaslabel=hvoxlabel, atlasdata=hvoxdata), name='scorehv')
        score207 = pe.Node(cbfqroiquant(atlaslabel=sc207label, atlasdata=sc207data), name='score207')
        score217 = pe.Node(cbfqroiquant(atlaslabel=sc217label, atlasdata=sc217data), name='score217')
        score407 = pe.Node(cbfqroiquant(atlaslabel=sc407label, atlasdata=sc407data), name='score407')
        score417 = pe.Node(cbfqroiquant(atlaslabel=sc417label, atlasdata=sc417data), name='score417')
        scrubhv = pe.Node(cbfqroiquant(atlaslabel=hvoxlabel, atlasdata=hvoxdata), name='scrubhv')
        scrub207 = pe.Node(cbfqroiquant(atlaslabel=sc207label, atlasdata=sc207data), name='scrub207')
        scrub217 = pe.Node(cbfqroiquant(atlaslabel=sc217label, atlasdata=sc217data), name='scrub217')
        scrub407 = pe.Node(cbfqroiquant(atlaslabel=sc407label, atlasdata=sc407data), name='scrub407')
        scrub417 = pe.Node(cbfqroiquant(atlaslabel=sc417label, atlasdata=sc417data), name='scrub417')
    
    if basil:
        basilhv = pe.Node(cbfqroiquant(atlaslabel=hvoxlabel, atlasdata=hvoxdata), name='basilhv')
        basil207 = pe.Node(cbfqroiquant(atlaslabel=sc207label, atlasdata=sc207data), name='basil207')
        basil217 = pe.Node(cbfqroiquant(atlaslabel=sc217label, atlasdata=sc217data), name='basil217')
        basil407 = pe.Node(cbfqroiquant(atlaslabel=sc407label, atlasdata=sc407data), name='basil407')
        basil417 = pe.Node(cbfqroiquant(atlaslabel=sc417label, atlasdata=sc417data), name='basil417')
        pvchv = pe.Node(cbfqroiquant(atlaslabel=hvoxlabel, atlasdata=hvoxdata), name='pvchv')
        pvc207 = pe.Node(cbfqroiquant(atlaslabel=sc207label, atlasdata=sc207data), name='pvc207')
        pvc217 = pe.Node(cbfqroiquant(atlaslabel=sc217label, atlasdata=sc217data), name='pvc217')
        pvc407 = pe.Node(cbfqroiquant(atlaslabel=sc407label, atlasdata=sc407data), name='pvc407')
        pvc417 = pe.Node(cbfqroiquant(atlaslabel=sc417label, atlasdata=sc417data), name='pvc417')
    

    
    
    workflow.connect([(inputnode, mrg_xfms, [('t1_asl_xform', 'in1'),
                                             ('std2anat_xfm', 'in2')]),
                      (inputnode, hvoftrans, [('aslmask', 'reference_image')]),
                      (mrg_xfms, hvoftrans, [('out', 'transforms')]),
                      (inputnode, sc207trans, [('aslmask', 'reference_image')]),
                      (mrg_xfms, sc207trans, [('out', 'transforms')]),
                      (inputnode, sc217trans, [('aslmask', 'reference_image')]),
                      (mrg_xfms, sc217trans, [('out', 'transforms')]),
                      (inputnode, sc407trans, [('aslmask', 'reference_image')]),
                      (mrg_xfms, sc407trans, [('out', 'transforms')]),
                      (inputnode, sc417trans, [('aslmask', 'reference_image')]),
                      (mrg_xfms, sc417trans, [('out', 'transforms')]),
                      (hvoftrans, cbfroihv, [('output_image', 'atlasfile')]),
                      (sc207trans, cbfroi207, [('output_image', 'atlasfile')]),
                      (sc217trans, cbfroi217, [('output_image', 'atlasfile')]),
                      (sc407trans, cbfroi407, [('output_image', 'atlasfile')]),
                      (sc417trans, cbfroi417, [('output_image', 'atlasfile')]),
                      (inputnode, cbfroihv, [('cbf', 'in_cbf')]),
                      (inputnode, cbfroi207, [('cbf', 'in_cbf')]),
                      (inputnode, cbfroi217, [('cbf', 'in_cbf')]),
                      (inputnode, cbfroi407, [(('cbf', 'in_cbf'))]),
                      (inputnode, cbfroi417, [(('cbf', 'in_cbf'))]),
                      (cbfroihv, outputnode, [('atlascsv', 'cbf_hvoxf')]),
                      (cbfroi207, outputnode, [('atlascsv', 'cbf_sc207')]),
                      (cbfroi217, outputnode, [('atlascsv', 'cbf_sc217')]),
                      (cbfroi407, outputnode, [('atlascsv', 'cbf_sc407')]),
                      (cbfroi417, outputnode, [('atlascsv', 'cbf_sc417')]),
                      ])

    if scorescrub: 
        workflow.connect([
            (hvoftrans, scorehv, [('output_image', 'atlasfile')]),
            (hvoftrans, scrubhv, [('output_image', 'atlasfile')]),
            (sc207trans, score207, [('output_image', 'atlasfile')]),
            (sc207trans, scrub207, [('output_image', 'atlasfile')]),
            (sc217trans, score217, [('output_image', 'atlasfile')]),
            (sc217trans, scrub217, [('output_image', 'atlasfile')]),
            (sc407trans, score407, [('output_image', 'atlasfile')]),
            (sc407trans, scrub407, [('output_image', 'atlasfile')]),
            (sc417trans, score417, [('output_image', 'atlasfile')]),
            (sc417trans, scrub417, [('output_image', 'atlasfile')]),
            (inputnode, scorehv, [('score', 'in_cbf')]),
            (inputnode, scrubhv, [('scrub', 'in_cbf')]),
            (inputnode, score207, [('score', 'in_cbf')]),
            (inputnode, scrub207, [('scrub', 'in_cbf')]),
            (inputnode, score217, [('score', 'in_cbf')]),
            (inputnode, scrub217, [(('scrub', 'in_cbf'))]),
            (inputnode, score407, [(('score', 'in_cbf'))]),
            (inputnode, scrub407, [(('scrub', 'in_cbf'))]),
            (inputnode, score417, [(('score', 'in_cbf'))]),
            (inputnode, scrub417, [(('scrub', 'in_cbf'))]),
            (scorehv, outputnode, [('atlascsv', 'score_hvoxf')]),
            (score207, outputnode, [('atlascsv', 'score_sc207')]),
            (score217, outputnode, [('atlascsv', 'score_sc217')]),
            (score407, outputnode, [('atlascsv', 'score_sc407')]),
            (score417, outputnode, [('atlascsv', 'score_sc417')]),
            (scrubhv, outputnode, [('atlascsv', 'scrub_hvoxf')]),
            (scrub207, outputnode, [('atlascsv', 'scrub_sc207')]),
            (scrub217, outputnode, [('atlascsv', 'scrub_sc217')]),
            (scrub407, outputnode, [('atlascsv', 'scrub_sc407')]),
            (scrub417, outputnode, [('atlascsv', 'scrub_sc417')]),
        ])

    if basil:
        workflow.connect([
           (hvoftrans, basilhv, [('output_image', 'atlasfile')]),
           (hvoftrans, pvchv, [('output_image', 'atlasfile')]),
           (sc207trans, pvc207, [('output_image', 'atlasfile')]),
            (sc207trans, basil207, [('output_image', 'atlasfile')]),
            (sc217trans, pvc217, [('output_image', 'atlasfile')]),
            (sc217trans, basil217, [('output_image', 'atlasfile')]),
            (sc407trans, pvc407, [('output_image', 'atlasfile')]),
            (sc407trans, basil407, [('output_image', 'atlasfile')]),
            (sc417trans, pvc417, [('output_image', 'atlasfile')]),
            (sc417trans, basil417, [('output_image', 'atlasfile')]),
           (inputnode, basilhv, [('basil', 'in_cbf')]),
           (inputnode, pvchv, [('pvc', 'in_cbf')]),
           (inputnode, basil207, [('basil', 'in_cbf')]),
           (inputnode, pvc207, [('pvc', 'in_cbf')]),
           (inputnode, basil217, [(('basil', 'in_cbf'))]),
           (inputnode, pvc217, [(('pvc', 'in_cbf'))]),
           (inputnode, basil407, [(('basil', 'in_cbf'))]),
           (inputnode, pvc407, [(('pvc', 'in_cbf'))]),
           (inputnode, basil417, [('basil', 'in_cbf')]),
           (inputnode, pvc417, [('pvc', 'in_cbf')]),
           (basilhv, outputnode, [('atlascsv', 'basil_hvoxf')]),
            (basil207, outputnode, [('atlascsv', 'basil_sc207')]),
            (basil217, outputnode, [('atlascsv', 'basil_sc217')]),
            (basil407, outputnode, [('atlascsv', 'basil_sc407')]),
            (basil417, outputnode, [('atlascsv', 'basil_sc417')]),
            (pvchv, outputnode, [('atlascsv', 'pvc_hvoxf')]),
            (pvc207, outputnode, [('atlascsv', 'pvc_sc207')]),
            (pvc217, outputnode, [('atlascsv', 'pvc_sc217')]),
            (pvc407, outputnode, [('atlascsv', 'pvc_sc407')]),
            (pvc417, outputnode, [('atlascsv', 'pvc_sc417')]),
        ])
    return workflow


def init_gecbf_compt_wf(metadata, asl_file,mem_gb, bids_dir,omp_nthreads,M0Scale=1,smooth_kernel=5,
                       scorescrub=False,basil=False,name='cbf_compt_wf'):
    """
    be back
    """

                    
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The CBF was quantified from  *preproccessed* ASL data  using a relatively basic
model [@detre_perfusion] [@alsop_recommended]. 
"""
    inputnode = pe.Node(niu.IdentityInterface(fields=['asl_file', 'in_file', 'asl_mask',
                                                      't1w_tpms', 't1w_mask', 't1_asl_xform',
                                                      'itk_asl_to_t1','m0_file']),
                        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_cbf', 'out_mean', 'out_score',
                                                       'out_avgscore', 'out_scrub', 'out_cbfb',
                                                       'out_scoreindex', 'out_cbfpv','out_att']),
                         name='outputnode')
    # convert tmps to asl_space
    csf_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                      name='csf_tfm', mem_gb=0.1)
    wm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='wm_tfm', mem_gb=0.1)
    gm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='gm_tfm', mem_gb=0.1)
    


    filex = os.path.abspath(asl_file)
    aslcontext1 = filex.replace('_asl.nii.gz', '_aslcontext.tsv')
    aslcontext = pd.read_csv(aslcontext1)
    idasl = aslcontext['volume_type'].tolist()
    deltamlist = [i for i in range(0, len(idasl)) if idasl[i] == 'deltam']
    cbflist = [i for i in range(0, len(idasl)) if idasl[i] == 'CBF']
    
    tiscbf = np.add(metadata['PostLabelingDelay'],metadata['LabelingDuration'])

    if hasattr(tiscbf, '__len__'):
        tisasl = ",".join([str(i) for i in tiscbf])
    else:
        tisasl = str(tiscbf)
    def pcaslorasl(metadata):
        if 'CASL' in metadata["LabelingType"]:
            pcasl1 = True
        elif 'PASL' in metadata["LabelingType"]:
            pcasl1 = False
        return pcasl1

    computecbf = pe.Node(computeCBF(in_metadata=metadata,in_m0scale=M0Scale), mem_gb=mem_gb,
                         run_without_submitting=True, name="computecbf")
    
    
    
    

    refinemaskj = pe.Node(refinemask(), mem_gb=1, run_without_submitting=True, name="refinemask")

    def _pick_csf(files):
        return files[-1]

    def _pick_gm(files):
        return files[0]

    def _pick_wm(files):
        return files[1]

    def _getfiledir(file):
        import os
        return os.path.dirname(file)

    if  len(deltamlist) > 0 :
        extractcbf1 = pe.Node(extractCB(file_type='d'),mem_gb=1,
                         run_without_submitting=True, name="extract_d")
        
        workflow.connect([
        # extract CBF data and compute cbf
        (inputnode, extractcbf1,[('asl_file','in_asl'),('asl_mask','in_aslmask')]),
        (extractcbf1, computecbf,[('out_file','in_cbf')]),
        (inputnode, computecbf, [('m0_file', 'in_m0file')]),
        # (inputnode,computecbf,[('asl_mask','in_mask')]),
        (inputnode, refinemaskj, [('t1w_mask', 'in_t1mask'), ('asl_mask', 'in_aslmask'),
                                  ('t1_asl_xform', 'transforms')]),
        (refinemaskj, computecbf, [('out_mask', 'in_mask')]),
        
        # extract probability maps
        (inputnode, csf_tfm, [('asl_mask', 'reference_image'),
                              ('t1_asl_xform', 'transforms')]),
        (inputnode, csf_tfm, [(('t1w_tpms', _pick_csf), 'input_image')]),
        (inputnode, wm_tfm, [('asl_mask', 'reference_image'),
                             ('t1_asl_xform', 'transforms')]),
        (inputnode, wm_tfm, [(('t1w_tpms', _pick_wm), 'input_image')]),
        (inputnode, gm_tfm, [('asl_mask', 'reference_image'),
                             ('t1_asl_xform', 'transforms')]),
        (inputnode, gm_tfm, [(('t1w_tpms', _pick_gm), 'input_image')]),
        (computecbf, outputnode, [('out_cbf', 'out_cbf'),
                                  ('out_mean', 'out_mean')]),
        ])
        if scorescrub:
            workflow.__desc__ = workflow.__desc__ +  """\
CBF are susceptible to artifacts due to low signal to noise ratio  and  sensitivity to  motion, Structural Correlation
based Outlier Rejection (SCORE) algothim was applied to the CBF to discard few extreme
outliers [@score_dolui]. Furthermore,Structural Correlation with RobUst Bayesian (SCRUB)
algorithms was applied to the CBF by iteratively reweighted  CBF  with structural tissues
probalility maps [@scrub_dolui].i 
"""
            scorescrub1 = pe.Node(scorescrubCBF(in_thresh=0.7, in_wfun='huber'), mem_gb=mem_gb,
                         name='scorescrub', run_without_submitting=True)
            workflow.connect([
            (computecbf, scorescrub1, [('out_cbf', 'in_file')]),
            (gm_tfm, scorescrub1, [('output_image', 'in_greyM')]),
            (refinemaskj, scorescrub1, [('out_mask', 'in_mask')]),
            (wm_tfm, scorescrub1, [('output_image', 'in_whiteM')]),
            (csf_tfm, scorescrub1, [('output_image', 'in_csf')]),
            (scorescrub1, outputnode, [('out_score', 'out_score'), ('out_scoreindex', 'out_scoreindex'),
                                  ('out_avgscore', 'out_avgscore'), ('out_scrub', 'out_scrub')]),
            ])
        if basil:
            workflow.__desc__ = workflow.__desc__ + """\
Alternate method of CBF computation is Bayesian Inference
for Arterial Spin Labeling (BASIL) as implmented in FSL which is  based on Bayeisan inference
principles [@chappell_basil]. BASIL computed the CBF from ASL incoporating natural varaibility
of other model parameters and spatial regularization of the estimated perfusion image. BASIL
also included correction for partial volume effects [@chappell_pvc].
"""         
            basilcbf = pe.Node(BASILCBF(m0scale=M0Scale, bolus=metadata["LabelingDuration"],
                                m0tr=metadata['RepetitionTime'], pvc=True,tis = tisasl,
                                pcasl = pcaslorasl(metadata = metadata)), name='basilcbf',
                       run_without_submitting=True, mem_gb=mem_gb)
            workflow.connect([
             (inputnode, basilcbf, [('asl_file', 'in_file')]),
             (gm_tfm, basilcbf, [('output_image', 'pvgm')]),
             (wm_tfm, basilcbf, [('output_image', 'pvwm')]),
             (inputnode, basilcbf, [(('asl_mask', _getfiledir), 'out_basename')]),
             (inputnode, basilcbf, [('m0_file', 'mzero')]),
             (refinemaskj, basilcbf, [('out_mask', 'mask')]),
             (basilcbf, outputnode, [('out_cbfb', 'out_cbfb'),
                                ('out_cbfpv', 'out_cbfpv'),
                                ('out_att', 'out_att')]),
             ])
    elif len(cbflist) > 0: 
        basil = False
        from nipype.interfaces.fsl import MultiImageMaths
        mask_cbf = pe.Node(MultiImageMaths(op_string=" -mul  %s "),mem_gb=1, run_without_submitting=True, name="mask_cbf")
         
        extractcbf2 = pe.Node(extractCB(file_type='c'),mem_gb=1,
                         run_without_submitting=True, name="extract_c")
        workflow.connect([
        (inputnode, refinemaskj, [('t1w_mask', 'in_t1mask'), ('asl_mask', 'in_aslmask'),
                                  ('t1_asl_xform', 'transforms')]),
        (refinemaskj, mask_cbf,[('out_mask','in_file')]),
        (inputnode, extractcbf2,[('asl_file','in_asl'),('asl_mask','in_aslmask')]),
        (extractcbf2, mask_cbf,[('out_file','operand_files')]),
        (mask_cbf, outputnode, [('out_file', 'out_cbf'),
                                  ('out_file', 'out_mean')]),
        ])

        if scorescrub:
            workflow.__desc__ = workflow.__desc__ + """\
CBF are susceptible to artifacts due to low signal to noise ratio  and  sensitivity to  motion, Structural Correlation
based Outlier Rejection (SCORE) algothim was applied to the CBF to discard few extreme
outliers [@score_dolui]. Furthermore,Structural Correlation with RobUst Bayesian (SCRUB)
algorithms was applied to the CBF by iteratively reweighted  CBF  with structural tissues
probalility maps [@scrub_dolui]. 
"""         
            scorescrub1 = pe.Node(scorescrubCBF(in_thresh=0.7, in_wfun='huber'), mem_gb=mem_gb,
                         name='scorescrub', run_without_submitting=True)
            workflow.connect([
            
             (inputnode, csf_tfm, [('asl_mask', 'reference_image'),
                              ('t1_asl_xform', 'transforms')]),
             (inputnode, csf_tfm, [(('t1w_tpms', _pick_csf), 'input_image')]),
             (inputnode, wm_tfm, [('asl_mask', 'reference_image'),
                             ('t1_asl_xform', 'transforms')]),
             (inputnode, wm_tfm, [(('t1w_tpms', _pick_wm), 'input_image')]),
             (inputnode, gm_tfm, [('asl_mask', 'reference_image'),
                             ('t1_asl_xform', 'transforms')]),
             (inputnode, gm_tfm, [(('t1w_tpms', _pick_gm), 'input_image')]),
             (inputnode, scorescrub1, [('asl_file', 'in_file')]),
             (gm_tfm, scorescrub1, [('output_image', 'in_greyM')]),
             (wm_tfm, scorescrub1, [('output_image', 'in_whiteM')]),
             (csf_tfm, scorescrub1, [('output_image', 'in_csf')]),
             (scorescrub1, outputnode, [('out_score', 'out_score'), ('out_scoreindex', 'out_scoreindex'),
                                  ('out_avgscore', 'out_avgscore'), ('out_scrub', 'out_scrub')]),
             ])
    return workflow


