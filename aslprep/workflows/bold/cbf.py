from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl
from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
from ...niworkflows.interfaces import NormalizeMotionParams
from ...niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from ...niworkflows.interfaces.itk import MCFLIRT2ITK
from ...niworkflows.interfaces.cbf_computation import (extractCBF,computeCBF
       ,scorescrubCBF,BASILCBF,refinemask,qccbf)
from ...niworkflows.interfaces.utility import KeySelect
from ...niworkflows.interfaces.plotting import (CBFSummary,CBFtsSummary)
from ...interfaces import  DerivativesDataSink
import nibabel as nb 
import numpy as np
import os,sys
from ...config import DEFAULT_MEMORY_MIN_GB


def init_cbf_compt_wf(mem_gb,metadata,aslcontext,pcasl,omp_nthreads, name='cbf_compt_wf'):
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The CBF was quantified from  *preproccessed* ASL data using a relatively basic model 
[@detre_perfusion] [@alsop_recommended]. CBF are susceptible to artifacts due to low signal to noise ratio  and  sensitivity 
to  motion, Structural Correlation based Outlier Rejection (SCORE) algothim was applied to the CBF to 
discard few extreme outliers [@score_dolui]. Furthermore,Structural Correlation with RobUst Bayesian (SCRUB)
algorithms was applied to the CBF by iteratively reweighted  CBF  with structural tissues probalility maps 
[@scrub_dolui].  Alternate method of CBF computation is Bayesian Inference for Arterial Spin Labeling (BASIL) 
as implmented in FSL which is  based on Bayeisan inference principles [@chappell_basil]. 
BASIL computed the CBF from ASL incoporating natural varaibility of other model parameters and spatial regularization 
of the estimated perfusion image. BASIL also included correction for partial volume effects [@chappell_pvc].  
"""
    


    inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold', 'bold_mask','t1w_tpms','t1w_mask','t1_bold_xform','itk_bold_to_t1']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_cbf', 'out_mean','out_score','out_avgscore','out_scrub',
             'out_scoreindex','out_cbfb','out_cbfpv']),
        name='outputnode')


    
    # convert tmps to bold_space
    csf_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                      name='csf_tfm', mem_gb=0.1)
    wm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='wm_tfm', mem_gb=0.1)
    gm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='gm_tfm', mem_gb=0.1)
    
     
    labeltype=metadata['LabelingType']
    if 'CASL' in labeltype: 
        pcasl=True
    elif 'PASL' in labeltype:
        pcasl=False
    else:
        print('unknown label type')

        
    extractcbf = pe.Node(extractCBF(in_ASLcontext=aslcontext),mem_gb=0.2,run_without_submitting=True,name="extractcbf") 
    computecbf = pe.Node(computeCBF(in_metadata=metadata),mem_gb=0.2,
              run_without_submitting=True,name="computecbf")
    scorescrub= pe.Node(scorescrubCBF(in_thresh=0.7,in_wfun='huber'),
              name='scorescrub',run_without_submitting=True,mem_gb=0.2)
    basilcbf= pe.Node(BASILCBF(m0scale=metadata["M0"],
               bolus=metadata["InitialPostLabelDelay"],m0tr=metadata['RepetitionTime'],pvc=True,
               tis=np.add(metadata["InitialPostLabelDelay"],metadata["LabelingDuration"]),
               pcasl=pcasl,out_basename=os.getcwd()),
              name='basilcbf',run_without_submitting=True,mem_gb=0.2) 
    
   
    refinemaskj=pe.Node(refinemask(),mem_gb=0.2,run_without_submitting=True,name="refinemask")
    
    #def _getTR(file):
        #import nibabel as nb
        #motr=nb.load(file).header.get_zooms()[3]
        #return motr
    
    def _pick_csf(files):
        return files[0]
    
    def _pick_gm(files):
        return files[1]

    def _pick_wm(files):
        return files[-1]

    
    workflow.connect([
        # extract CBF data and compute cbf
        (inputnode,  extractcbf, [('bold','in_file')]),
        (extractcbf, computecbf, [('out_file','in_cbf'),('out_avg','in_m0file')]),
        #(inputnode,computecbf,[('bold_mask','in_mask')]),
        (inputnode,refinemaskj,[('t1w_mask','in_t1mask'),('bold_mask','in_boldmask'),
                                ('t1_bold_xform','transforms')]),
        (inputnode,computecbf,[('bold_mask','in_mask')]),
        (inputnode,scorescrub,[('bold_mask','in_mask')]),
        (inputnode,basilcbf,[('bold_mask','mask')]),

        # extract probability maps
        (inputnode, csf_tfm, [('bold_mask', 'reference_image'),
                              ('t1_bold_xform', 'transforms')]),
        (inputnode, csf_tfm, [(('t1w_tpms', _pick_csf), 'input_image')]),
        (inputnode, wm_tfm, [('bold_mask', 'reference_image'),
                              ('t1_bold_xform', 'transforms')]),
        (inputnode, wm_tfm, [(('t1w_tpms', _pick_wm), 'input_image')]),
        (inputnode, gm_tfm, [('bold_mask', 'reference_image'),
                              ('t1_bold_xform', 'transforms')]),
        (inputnode, gm_tfm, [(('t1w_tpms', _pick_gm), 'input_image')]),
        (computecbf,scorescrub,[('out_cbf','in_file')]),
        (gm_tfm,scorescrub,[('output_image','in_greyM')]),
        (wm_tfm,scorescrub,[('output_image','in_whiteM')]),
        (csf_tfm,scorescrub,[('output_image','in_csf')]),
        #(inputnode,scorescrub,[('bold_mask','in_mask')]),
        (extractcbf,basilcbf,[('out_file','in_file')]),
        (gm_tfm,basilcbf,[('output_image','pvgm')]),
        (wm_tfm,basilcbf,[('output_image','pvwm')]),
        #(inputnode,basilcbf,[('bold_mask','mask')]),
        (extractcbf,basilcbf,[('out_avg','mzero')]),
        (basilcbf,outputnode,[('out_cbfb','out_cbfb'),
                ('out_cbfpv','out_cbfpv')]),
        (computecbf,outputnode,[('out_cbf','out_cbf'),
                     ('out_mean','out_mean')]),
        (scorescrub,outputnode,[('out_score','out_score'),('out_scoreindex','out_scoreindex'),
                    ('out_avgscore','out_avgscore'),('out_scrub','out_scrub')]),
        
        ])
    return workflow


def init_cbfqc_compt_wf(mem_gb,bold_file,metadata,omp_nthreads, name='cbfqc_compt_wf'):
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The following quality control (qc) measures was estimated: framewise displacement and relative root mean square dice index. 
Other qc meaure include dice and jaccard indices, cross-correlation and coverage that estimate the coregistration  
quality of  ASL and T1W images and  normalization quality of ASL to template. Quality evaluation index (QEI) 
was also computed for CBF [@cbfqc]. The  QEI is  automated for objective quality evaluation of CBF maps and measured 
the CBF quality based on structural similarity,spatial variability and the percentatge  of voxels with  negtaive  CBF within Grey matter 

"""
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['meancbf','avgscore','scrub','basil','pv','bold_mask','t1w_tpms','t1w_mask','t1_bold_xform','bold_mask_std',
        'confmat']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['qc_file']),
        name='outputnode')
    
    def _pick_csf(files):
        return files[0]
    
    def _pick_gm(files):
        return files[1]

    def _pick_wm(files):
        return files[-1]
    
    csf_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                      name='csf_tfm', mem_gb=0.1)
    wm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='wm_tfm', mem_gb=0.1)
    gm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='gm_tfm', mem_gb=0.1)

    mask_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='masktonative', mem_gb=0.1)

    from templateflow.api import get as get_template  
    brain_mask = str(get_template(
            'MNI152NLin2009cAsym', resolution=2, desc='brain', suffix='mask'))

    from nipype.interfaces.afni  import Resample
    resample = pe.Node(Resample(in_file=brain_mask,outputtype='NIFTI_GZ'),name='resample', mem_gb=0.1)

    #template_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True,input_image=brain_mask),
                     #name='template_tfm', mem_gb=0.1)
    qccompute=pe.Node(qccbf(in_file=bold_file),name='qccompute',run_without_submitting=True,mem_gb=0.2)
    
    workflow.connect([(inputnode, csf_tfm, [('bold_mask', 'reference_image'),
                              ('t1_bold_xform', 'transforms')]),
                      (inputnode, csf_tfm, [(('t1w_tpms', _pick_csf), 'input_image')]),
                      (inputnode, wm_tfm, [('bold_mask', 'reference_image'),
                              ('t1_bold_xform', 'transforms')]),
                      (inputnode, wm_tfm, [(('t1w_tpms', _pick_wm), 'input_image')]),
                      (inputnode, gm_tfm, [('bold_mask', 'reference_image'),
                              ('t1_bold_xform', 'transforms')]),
                      (inputnode, gm_tfm, [(('t1w_tpms', _pick_gm), 'input_image')]),
                      (inputnode, mask_tfm, [('bold_mask', 'reference_image'),
                              ('t1_bold_xform', 'transforms'),('t1w_mask', 'input_image')]),
                      (mask_tfm,qccompute,[('output_image','in_t1mask')]),
                      (inputnode,qccompute,[('bold_mask','in_boldmask'),
                            ('confmat','in_confmat')]),
                      (inputnode,qccompute,[(('bold_mask_std',_pick_csf),'in_boldmaskstd')]),
                     (inputnode,resample,[(('bold_mask_std',_pick_csf),'master')]),
                     (resample,qccompute,[('out_file','in_templatemask')]),
                     (gm_tfm,qccompute,[('output_image','in_greyM')]),
                     (wm_tfm,qccompute,[('output_image','in_whiteM')]),
                     (csf_tfm,qccompute,[('output_image','in_csf')]),    
                     (inputnode,qccompute,[('scrub','in_scrub'),
                                ('meancbf','in_meancbf'),('avgscore','in_avgscore'),
                                ('basil','in_basil'),('pv','in_pvc')]),
                     (qccompute,outputnode,[('qc_file','qc_file')]), 
                    ])
    return workflow


def init_cbfplot_wf(mem_gb,metadata,omp_nthreads, name='cbf_plot'):
    workflow = Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['cbf', 'cbf_ts','score_ts','score','scrub','bold_ref',
            'basil','pvc','bold_mask','t1_bold_xform','std2anat_xfm']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['cbf_carpetplot','cbf_summary_plot']),
        name='outputnode')
    mrg_xfms = pe.Node(niu.Merge(2), name='mrg_xfms')
     
    from templateflow.api import get as get_template
    resample_parc = pe.Node(ApplyTransforms(
        float=True,
        input_image=str(get_template(
            'MNI152NLin2009cAsym', resolution=1, desc='carpet',
            suffix='dseg', extension=['.nii', '.nii.gz'])),
        dimension=3, default_value=0, interpolation='MultiLabel'),
        name='resample_parc')
    


    cbfsummary=pe.Node(CBFSummary(),name='cbf_summary',mem_gb=0.2)
    cbftssummary=pe.Node(CBFtsSummary(tr=metadata['RepetitionTime']),name='cbf_ts_summary',mem_gb=0.2)

    ds_report_cbfplot = pe.Node(
        DerivativesDataSink(desc='cbfplot', keep_dtype=True),
        name='ds_report_cbfplot', run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)
    ds_report_cbftsplot = pe.Node(
        DerivativesDataSink(desc='cbftsplot', keep_dtype=True),
        name='ds_report_cbftsplot', run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)

    
    workflow.connect([(inputnode, mrg_xfms, [('t1_bold_xform', 'in1'),
                               ('std2anat_xfm', 'in2')]),
                      (inputnode, resample_parc, [('bold_mask', 'reference_image'),
                              ('t1_bold_xform', 'transforms')]),
                      (resample_parc,cbftssummary,[('output_image','seg_file')]),
                      (inputnode,cbftssummary,[('cbf_ts','cbf_ts'),('score_ts','score_ts')]),
                      (cbftssummary,ds_report_cbftsplot,[('out_file','in_file')]),
                      (cbftssummary,outputnode,[('out_file','cbf_carpetplot')]),
                      (inputnode ,cbfsummary,[('cbf','cbf'),('score','score'),
                                  ('scrub','scrub'),('basil','basil'),('pvc','pvc'),
                                ('bold_ref','ref_vol')]),
                      (cbfsummary,ds_report_cbfplot,[('out_file','in_file')]),
                      (cbfsummary,outputnode,[('out_file','cbf_summary_plot')]),

    ])
    return workflow
        
 
        


       
        
              
        


    


    


