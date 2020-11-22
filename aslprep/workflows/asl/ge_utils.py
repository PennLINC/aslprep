# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
utils code to process GE scan with fewer volumes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

"""

from pathlib import Path
import nibabel as nb
import numpy as np
import os 
import pandas as pd
import os
import os.path as op
import pkg_resources as pkgr
from nipype.utils.filemanip import fname_presuffix
from ...niworkflows.func.util import init_enhance_and_skullstrip_asl_wf
from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
from ...niworkflows.interfaces.masks import SimpleShowMaskRPT 
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl, c3
from nipype.interfaces import fsl
from ... import config

DEFAULT_MEMORY_MIN_GB = config.DEFAULT_MEMORY_MIN_GB
LOGGER = config.loggers.workflow


def init_asl_geref_wf(omp_nthreads,mem_gb,metadata,bids_dir,brainmask_thresh=0.5,pre_mask=False, name="asl_gereference_wf",
    gen_report=False):

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
First, a reference volume and its skull-stripped version were generated.
        """
    
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file"
            ]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "raw_ref_image",
                "ref_image_brain",
                "asl_mask",
                "m0_file",
                "mask_report",
            ]
        ),
        name="outputnode",
    )

    gen_ref = pe.Node(GeReferenceFile(bids_dir=bids_dir, in_metadata=metadata),
               omp_nthreads=1,mem_gb=1,name='gen_ge_ref')
    gen_ref.base_dir=os.getcwd()
    skull_strip_wf =  pe.Node(
        fsl.BET(frac=0.5, mask=True), name="fslbet")
    apply_mask = pe.Node(fsl.ApplyMask(), name="apply_mask")
    mask_reportlet = pe.Node(SimpleShowMaskRPT(), name="mask_reportlet")

    workflow.connect([
         (inputnode,gen_ref,[('asl_file','in_file')]),
         (gen_ref,skull_strip_wf,[('ref_file','in_file')]),
         (gen_ref, outputnode, [
            ("ref_file", "raw_ref_image")]),
        (gen_ref, apply_mask, [
                ("ref_file", "in_file")]),
        (skull_strip_wf, outputnode, [
            ("mask_file", "asl_mask")]),
        (skull_strip_wf, apply_mask, [
            ("mask_file", "mask_file")]),
        (apply_mask,outputnode, [
                ("out_file", "ref_image_brain")]),
         (gen_ref,mask_reportlet,[("ref_file", "background_file")]),
         (skull_strip_wf,mask_reportlet,[('mask_file','mask_file')]),
         (gen_ref,outputnode,[("m0_file", "m0_file")]),
         ])
    return workflow



def init_asl_gereg_wf(use_bbr,asl2t1w_dof,asl2t1w_init,
        mem_gb, omp_nthreads, name='asl_reg_wf',
        sloppy=False, use_compression=True, write_report=True):
    
    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['ref_asl_brain', 't1w_brain', 't1w_dseg']),name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['itk_asl_to_t1','itk_t1_to_asl',
             'fallback']),name='outputnode')
    
    from .registration import init_fsl_bbr_wf

    bbr_wf = init_fsl_bbr_wf(use_bbr=use_bbr, asl2t1w_dof=asl2t1w_dof,
                                 asl2t1w_init=asl2t1w_init, sloppy=sloppy)
    #bbr_wf.base_dir=os.getcwd()
    from ...interfaces import DerivativesDataSink

    workflow.connect([
        (inputnode, bbr_wf, [
            ('ref_asl_brain', 'inputnode.in_file'),
            ('t1w_dseg', 'inputnode.t1w_dseg'),
            ('t1w_brain', 'inputnode.t1w_brain')]),
        (bbr_wf, outputnode, [('outputnode.itk_asl_to_t1', 'itk_asl_to_t1'),
                              ('outputnode.itk_t1_to_asl', 'itk_t1_to_asl'),
                              ('outputnode.fallback', 'fallback')]),
    ])

    if write_report:
        ds_report_reg = pe.Node(
            DerivativesDataSink(datatype="figures", dismiss_entities=("echo",)),
            name='ds_report_reg', run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB)

        def _asl_reg_suffix(fallback):
            return 'flirtbbr' 

        workflow.connect([
            (bbr_wf, ds_report_reg, [
                ('outputnode.out_report', 'in_file'),
                (('outputnode.fallback', _asl_reg_suffix), 'desc')]),
        ])

    return workflow

def init_asl_t1_getrans_wf(mem_gb, omp_nthreads, cbft1space=False,scorescrub=False, basil=False,
                          use_compression=True, name='asl_t1_trans_wf'):
    """
    Co-register the reference ASL image to T1w-space.

    The workflow uses :abbr:`BBR (boundary-based registration)`.

    

    """
    from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from ...niworkflows.func.util import init_asl_reference_wf
    from ...niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
    from ...niworkflows.interfaces.itk import MultiApplyTransforms
    from ...niworkflows.interfaces.nilearn import Merge
    from ...niworkflows.interfaces.utils import GenerateSamplingReference

    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['name_source', 'ref_asl_brain', 'ref_asl_mask','asl_file',
                    't1w_brain', 't1w_mask','cbf', 'meancbf','att',
                    'score', 'avgscore', 'scrub', 'basil', 'pv', 'itk_asl_to_t1']),
        name='inputnode'
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=[
            'asl_t1', 'asl_t1_ref', 'asl_mask_t1','att_t1','cbf_t1', 'meancbf_t1', 
            'score_t1', 'avgscore_t1', 'scrub_t1', 'basil_t1', 'pv_t1']),
        name='outputnode'
    )

    #gen_ref = pe.Node(GenerateSamplingReference(), name='gen_ref',
                      #mem_gb=0.3)  # 256x256x256 * 64 / 8 ~ 150MB

    mask_t1w_tfm = pe.Node(ApplyTransforms(interpolation='MultiLabel'),
                           name='mask_t1w_tfm', mem_gb=0.1)

    workflow.connect([
        (inputnode, mask_t1w_tfm, [('ref_asl_mask', 'input_image')]),
        (inputnode, mask_t1w_tfm, [('t1w_brain', 'reference_image')]),
        (inputnode, mask_t1w_tfm, [('itk_asl_to_t1', 'transforms')]),
        (mask_t1w_tfm, outputnode, [('output_image', 'asl_mask_t1')]),
    ])

    asl_to_t1w_transform = pe.Node(
                  ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3,
                        dimension=3),
                  name='asl_to_t1w_transform', mem_gb=mem_gb)
   
    # Generate a reference on the target T1w space
   
    
    workflow.connect([
            (inputnode, asl_to_t1w_transform, [('ref_asl_brain', 'input_image')]),
            (inputnode, asl_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
            (inputnode, asl_to_t1w_transform, [('t1w_brain', 'reference_image')]),
            (asl_to_t1w_transform, outputnode, [('output_image', 'asl_t1')]),
        ])

    if cbft1space:
        cbf_to_t1w_transform = pe.Node(
               ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3,
                        dimension=3),
              name='cbf_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        meancbf_to_t1w_transform = pe.Node(
                         ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
                         name='meancbf_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        workflow.connect([
         
        (asl_to_t1w_transform, outputnode, [('output_image', 'asl_t1_ref')]),
        (inputnode, cbf_to_t1w_transform, [('cbf', 'input_image')]),
        (cbf_to_t1w_transform, outputnode, [('output_image', 'cbf_t1')]),
        (inputnode, cbf_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
        (inputnode, cbf_to_t1w_transform, [('t1w_brain', 'reference_image')]),
        (inputnode, meancbf_to_t1w_transform, [('meancbf', 'input_image')]),
        (meancbf_to_t1w_transform, outputnode, [('output_image', 'meancbf_t1')]),
        (inputnode, meancbf_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
        (inputnode, meancbf_to_t1w_transform, [('t1w_brain', 'reference_image')]),
        ])

    if cbft1space and scorescrub:

        score_to_t1w_transform = pe.Node(
             ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3,
                        dimension=3),
             name='score_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        avgscore_to_t1w_transform = pe.Node(
            ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
           name='avgscore_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        scrub_to_t1w_transform = pe.Node(
              ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
              name='scrub_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
    if cbft1space and basil:
        basil_to_t1w_transform = pe.Node(
               ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
            name='basil_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        pv_to_t1w_transform = pe.Node(
               ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
               name='pv_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        att_to_t1w_transform = pe.Node(
               ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3),
               name='att_to_t1w_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        
        workflow.connect([
         (inputnode, score_to_t1w_transform, [('score', 'input_image')]),
         (score_to_t1w_transform, outputnode, [('output_image', 'score_t1')]),
         (inputnode, score_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (inputnode, score_to_t1w_transform, [('t1w_brain', 'reference_image')]),

         (inputnode, avgscore_to_t1w_transform, [('avgscore', 'input_image')]),
         (avgscore_to_t1w_transform, outputnode, [('output_image', 'avgscore_t1')]),
         (inputnode, avgscore_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (inputnode, avgscore_to_t1w_transform, [('t1w_brain', 'reference_image')]),

         (inputnode, scrub_to_t1w_transform, [('scrub', 'input_image')]),
         (scrub_to_t1w_transform, outputnode, [('output_image', 'scrub_t1')]),
         (inputnode, scrub_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (inputnode, scrub_to_t1w_transform, [('t1w_brain', 'reference_image')]),
          ])
        
        workflow.connect([
         (inputnode, basil_to_t1w_transform, [('basil', 'input_image')]),
         (basil_to_t1w_transform, outputnode, [('output_image', 'basil_t1')]),
         (inputnode, basil_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (inputnode, basil_to_t1w_transform, [('t1w_brain', 'reference_image')]),

         (inputnode, pv_to_t1w_transform, [('pv', 'input_image')]),
         (pv_to_t1w_transform, outputnode, [('output_image', 'pv_t1')]),
         (inputnode, pv_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (inputnode, pv_to_t1w_transform, [('t1w_brain', 'reference_image')]),

         (inputnode, att_to_t1w_transform, [('att', 'input_image')]),
         (att_to_t1w_transform, outputnode, [('output_image', 'att_t1')]),
         (inputnode, att_to_t1w_transform, [('itk_asl_to_t1', 'transforms')]),
         (inputnode, att_to_t1w_transform, [('t1w_brain', 'reference_image')]),
         ])

    return workflow


def init_asl_gestd_trans_wf(
    mem_gb,
    omp_nthreads,
    spaces,
    scorescrub=False,
    basil=False,
    name='asl_gestd_trans_wf',
    use_compression=True,
):
    """

    """
    from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from ...niworkflows.func.util import init_asl_reference_wf
    from ...niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
    from ...niworkflows.interfaces.itk import MultiApplyTransforms
    from ...niworkflows.interfaces.utility import KeySelect
    from ...niworkflows.interfaces.utils import GenerateSamplingReference
    from ...niworkflows.interfaces.nilearn import Merge
    from ...niworkflows.utils.spaces import format_reference

    workflow = Workflow(name=name)
    output_references = spaces.cached.get_spaces(nonstandard=False, dim=(3,))
    std_vol_references = [
        (s.fullname, s.spec) for s in spaces.references if s.standard and s.dim == 3
    ]

    if len(output_references) == 1:
        workflow.__desc__ = """\
The ASL and CBF dreivatives  were resampled into standard space,
generating a *preprocessed ASL and computed CBF in {tpl} space*.
""".format(tpl=output_references[0])
    elif len(output_references) > 1:
        workflow.__desc__ = """\
The ASL and CBF dreivatives were resampled into several standard spaces,
correspondingly generating the following *spatially-normalized,
preprocessed ASL runs*: {tpl}.
""".format(tpl=', '.join(output_references))

    inputnode = pe.Node(
        niu.IdentityInterface(fields=[
            'anat2std_xfm',
            'cbf','meancbf','att','asl_file',
            'score','avgscore','scrub',
            'basil','pv','asl_mask',
            'itk_asl_to_t1',
            'name_source','templates',
        ]),
        name='inputnode'
    )

    iterablesource = pe.Node(
        niu.IdentityInterface(fields=['std_target']), name='iterablesource'
    )
    # Generate conversions for every template+spec at the input
    iterablesource.iterables = [('std_target', std_vol_references)]

    split_target = pe.Node(niu.Function(
        function=_split_spec, input_names=['in_target'],
        output_names=['space', 'template', 'spec']),
        run_without_submitting=True, name='split_target')

    select_std = pe.Node(KeySelect(fields=['anat2std_xfm']),
                         name='select_std', run_without_submitting=True)

    select_tpl = pe.Node(niu.Function(function=_select_template),
                         name='select_tpl', run_without_submitting=True)

    mask_std_tfm = pe.Node(ApplyTransforms(interpolation='MultiLabel'),
                           name='mask_std_tfm', mem_gb=1)

    # Write corrected file in the designated output dir
    mask_merge_tfms = pe.Node(niu.Merge(2), name='mask_merge_tfms', run_without_submitting=True,
                              mem_gb=DEFAULT_MEMORY_MIN_GB)
    nxforms = 3 
    merge_xforms = pe.Node(niu.Merge(nxforms), name='merge_xforms',
                           run_without_submitting=True, mem_gb=DEFAULT_MEMORY_MIN_GB)

    asl_to_std_transform = pe.Node(
        ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3,
                        dimension=3),
        name='asl_to_std_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
    cbf_to_std_transform = pe.Node(
        ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3,
                        dimension=3),
        name='cbf_to_std_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
    meancbf_to_std_transform = pe.Node(
        ApplyTransforms(interpolation="LanczosWindowedSinc", float=True,input_image_type=3),
        name='meancbf_to_std_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
    
    if scorescrub:
        score_to_std_transform = pe.Node(
            ApplyTransforms(interpolation="LanczosWindowedSinc", float=True, input_image_type=3,
                        dimension=3),
            name='score_to_std_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)

        avgscore_to_std_transform = pe.Node(
               ApplyTransforms(interpolation="LanczosWindowedSinc", float=True,input_image_type=3),
            name='avgscore_to_std_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)

        scrub_to_std_transform = pe.Node(
               ApplyTransforms(interpolation="LanczosWindowedSinc", float=True,input_image_type=3),
            name='scrub_to_std_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)

    if basil:
        basil_to_std_transform = pe.Node(
          ApplyTransforms(interpolation="LanczosWindowedSinc", float=True,input_image_type=3),
          name='basil_to_std_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)

        pv_to_std_transform = pe.Node(
          ApplyTransforms(interpolation="LanczosWindowedSinc", float=True,input_image_type=3),
          name='pv_to_std_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
        att_to_std_transform = pe.Node(
           ApplyTransforms(interpolation="LanczosWindowedSinc", float=True,input_image_type=3),
           name='att_to_std_transform', mem_gb=mem_gb * 3 * omp_nthreads, n_procs=omp_nthreads)
    #merge = pe.Node(Merge(compress=use_compression), name='merge',
                    #mem_gb=mem_gb * 3)
    mask_merge_tfms = pe.Node(niu.Merge(2), name='mask_merge_tfms', run_without_submitting=True,
                              mem_gb=DEFAULT_MEMORY_MIN_GB)
    # Generate a reference on the target standard space
    gen_ref = pe.Node(GenerateSamplingReference(), name='gen_ref',
                      mem_gb=0.3) 
    #gen_final_ref = init_asl_reference_wf(omp_nthreads=omp_nthreads, pre_mask=True)

    workflow.connect([
        (iterablesource, split_target, [('std_target', 'in_target')]),
        (iterablesource, select_tpl, [('std_target', 'template')]),
        (inputnode, select_std, [('anat2std_xfm', 'anat2std_xfm'),
                                 ('templates', 'keys')]),
        (inputnode, mask_std_tfm, [('asl_mask', 'input_image')]),
        (inputnode, gen_ref, [('asl_file', 'moving_image')]),
        (inputnode, merge_xforms, [
            (('itk_asl_to_t1', _aslist), 'in2')]),
        #(inputnode, merge, [('name_source', 'header_source')]),
        (inputnode, mask_merge_tfms, [(('itk_asl_to_t1', _aslist), 'in2')]),
        (inputnode, asl_to_std_transform, [('asl_file', 'input_image')]),
        (split_target, select_std, [('space', 'key')]),
        (select_std, merge_xforms, [('anat2std_xfm', 'in1')]),
        (select_std, mask_merge_tfms, [('anat2std_xfm', 'in1')]),
        (split_target, gen_ref, [(('spec', _is_native), 'keep_native')]),
        (select_tpl, gen_ref, [('out', 'fixed_image')]),
        (merge_xforms, asl_to_std_transform, [('out', 'transforms')]),
        (gen_ref, asl_to_std_transform, [('out_file', 'reference_image')]),
        (gen_ref, mask_std_tfm, [('out_file', 'reference_image')]),
        (mask_merge_tfms, mask_std_tfm, [('out', 'transforms')]),
        #(mask_std_tfm, gen_final_ref, [('output_image', 'inputnode.asl_mask')]),
        #(asl_to_std_transform, merge, [('output_image', 'in_files')]),
        #(inputnode, gen_final_ref, [('asl_file', 'inputnode.asl_file')]),
    ])


    output_names = [
        'asl_mask_std',
        'asl_std',
        'asl_std_ref',
        'spatial_reference',
        'template',
        'cbf_std',
        'meancbf_std',
    ] 

    if scorescrub:
        output_names = output_names +['score_std','avgscore_std','scrub_std']
    if basil:
        output_names = output_names + ['basil_std', 'pv_std','att_std']

    poutputnode = pe.Node(niu.IdentityInterface(fields=output_names),
                          name='poutputnode')


    workflow.connect([
        # Connecting outputnode
        (iterablesource, poutputnode, [
            (('std_target', format_reference), 'spatial_reference')]),
        (asl_to_std_transform, poutputnode, [('output_image', 'asl_std')]),
        (asl_to_std_transform, poutputnode, [('output_image', 'asl_std_ref')]),
        (mask_std_tfm, poutputnode, [('output_image', 'asl_mask_std')]),
        (select_std, poutputnode, [('key', 'template')]),

        (mask_merge_tfms, cbf_to_std_transform, [('out', 'transforms')]),
        (gen_ref, cbf_to_std_transform, [('out_file', 'reference_image')]),
        (inputnode, cbf_to_std_transform, [('cbf', 'input_image')]),
        (cbf_to_std_transform, poutputnode, [('output_image', 'cbf_std')]),

        (mask_merge_tfms, meancbf_to_std_transform, [('out', 'transforms')]),
        (gen_ref, meancbf_to_std_transform, [('out_file', 'reference_image')]),
        (inputnode, meancbf_to_std_transform, [('cbf', 'input_image')]),
        (meancbf_to_std_transform, poutputnode, [('output_image', 'meancbf_std')]),
    
        
      ])

    if scorescrub:
        workflow.connect([
         (mask_merge_tfms, avgscore_to_std_transform, [('out', 'transforms')]),
         (gen_ref, avgscore_to_std_transform, [('out_file', 'reference_image')]),
         (inputnode, avgscore_to_std_transform, [('avgscore', 'input_image')]),
         (avgscore_to_std_transform, poutputnode, [('output_image', 'avgscore_std')]),

        (mask_merge_tfms, score_to_std_transform, [('out', 'transforms')]),
        (gen_ref, score_to_std_transform, [('out_file', 'reference_image')]),
        (inputnode, score_to_std_transform, [('score', 'input_image')]),
        (score_to_std_transform, poutputnode, [('output_image', 'score_std')]),

         (mask_merge_tfms, scrub_to_std_transform, [('out', 'transforms')]),
         (gen_ref, scrub_to_std_transform, [('out_file', 'reference_image')]),
         (inputnode, scrub_to_std_transform, [('scrub', 'input_image')]),
         (scrub_to_std_transform, poutputnode, [('output_image', 'scrub_std')]),

         ])
    if basil:
        workflow.connect([
        (mask_merge_tfms, basil_to_std_transform, [('out', 'transforms')]),
        (gen_ref, basil_to_std_transform, [('out_file', 'reference_image')]),
        (inputnode, basil_to_std_transform, [('basil', 'input_image')]),
        (basil_to_std_transform, poutputnode, [('output_image', 'basil_std')]),

        (mask_merge_tfms, pv_to_std_transform, [('out', 'transforms')]),
        (gen_ref, pv_to_std_transform, [('out_file', 'reference_image')]),
        (inputnode, pv_to_std_transform, [('pv', 'input_image')]),
        (pv_to_std_transform, poutputnode, [('output_image', 'pv_std')]),

        (mask_merge_tfms, att_to_std_transform, [('out', 'transforms')]),
        (gen_ref, att_to_std_transform, [('out_file', 'reference_image')]),
        (inputnode, att_to_std_transform, [('att', 'input_image')]),
        (att_to_std_transform, poutputnode, [('output_image', 'att_std')]),
         ])
    
    # Connect parametric outputs to a Join outputnode
    outputnode = pe.JoinNode(niu.IdentityInterface(fields=output_names),
                             name='outputnode', joinsource='iterablesource')
    workflow.connect([
        (poutputnode, outputnode, [(f, f) for f in output_names]),
    ])
    return workflow

from nipype.interfaces.base import (
    traits,
    isdefined,
    File,
    InputMultiPath,
    TraitedSpec,
    BaseInterfaceInputSpec,
    SimpleInterface,
    DynamicTraitedSpec,
)

class _GenerateReferenceInputSpec(BaseInterfaceInputSpec):
    input_image = File(
        exists=True, mandatory=True, desc="input images"
    )
    

class _GenerateReferenceOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="one file with all inputs flattened")


class GenerateReference(SimpleInterface):
    """
    Generates a reference grid for resampling one image keeping original resolution,
    but moving data to a different space (e.g. MNI).

    """

    input_spec = _GenerateReferenceInputSpec
    output_spec = _GenerateReferenceOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = gen_reference(
            in_img=self.inputs.input_image)
        return runtime




def gen_reference(in_img, newpath=None):

    """generate reference for a GE scan with few volumes."""
    import nibabel as nb 
    import numpy as np
    import os 
    newpath = Path(newpath or ".")
    ss=check_img(in_img)
    if ss == 0: 
        ref_data=nb.load(in_img).get_fdata()
    else: 
        nii = nb.load(in_img).get_fdata()
        ref_data=np.mean(nii,axis=3)
    
    new_file = nb.Nifti1Image(dataobj=ref_data,header=nb.load(in_img).header,
             affine=nb.load(in_img).affine)
    out_file = fname_presuffix('aslref', suffix="_reference.nii.gz", newpath=str(newpath.absolute()))
    new_file.to_filename(out_file)
    return out_file

def check_img(img):
    # get the 4th dimension 
    ss=nb.load(img).get_fdata().shape
    if len(ss) == 3:
        ss=np.hstack([ss,0])
    return ss[3]

def _split_spec(in_target):
    space, spec = in_target
    template = space.split(':')[0]
    return space, template, spec


def _select_template(template):
    from niworkflows.utils.misc import get_template_specs
    template, specs = template
    template = template.split(':')[0]  # Drop any cohort modifier if present
    specs = specs.copy()
    specs['suffix'] = specs.get('suffix', 'T1w')

    # Sanitize resolution
    res = specs.pop('res', None) or specs.pop('resolution', None) or 'native'
    if res != 'native':
        specs['resolution'] = res
        return get_template_specs(template, template_spec=specs)[0]

    # Map nonstandard resolutions to existing resolutions
    specs['resolution'] = 2
    try:
        out = get_template_specs(template, template_spec=specs)
    except RuntimeError:
        specs['resolution'] = 1
        out = get_template_specs(template, template_spec=specs)

    return out[0]

def _first(inlist):
    return inlist[0]


def _aslist(in_value):
    if isinstance(in_value, list):
        return in_value
    return [in_value]


def _is_native(in_value):
    return (
        in_value.get('resolution') == 'native'
        or in_value.get('res') == 'native'
    )


class _GeReferenceFileInputSpec(BaseInterfaceInputSpec):
    in_file = File(
        exists=True, mandatory=True, desc="asl_file"
    )
    in_metadata = traits.Dict(exists=True, mandatory=True,
                              desc='metadata for asl or deltam ')
    bids_dir=traits.Str(exits=True,mandatory=True,desc=' bids directory')
    ref_file = File(exists=False,mandatory=False, desc="ref file")
    m0_file = File(exists=False,mandatory=False, desc="m0 file")
    

class _GeReferenceFileOutputSpec(TraitedSpec):
    ref_file = File(exists=True,mandatory=True,desc="ref file")
    m0_file = File(exists=True,mandatory=True, desc="m0 file")


class GeReferenceFile(SimpleInterface):
    """
    Generates a reference grid for resampling one image keeping original resolution,
    but moving data to a different space (e.g. MNI).

    """

    input_spec = _GeReferenceFileInputSpec
    output_spec = _GeReferenceFileOutputSpec

    def _run_interface(self, runtime):
        import os 
        filex = os.path.abspath(self.inputs.in_file)
        aslcontext1 = filex.replace('_asl.nii.gz', '_aslcontext.tsv')
        aslcontext = pd.read_csv(aslcontext1)
        idasl = aslcontext['volume_type'].tolist()
        m0list = [i for i in range(0, len(idasl)) if idasl[i] == 'm0scan']
        deltamlist = [i for i in range(0, len(idasl)) if idasl[i] == 'deltam']
        cbflist = [i for i in range(0, len(idasl)) if idasl[i] == 'CBF']
        allasl = nb.load(self.inputs.in_file)
        dataasl = allasl.get_fdata()

        if self.inputs.in_metadata['M0'] != "True" and type(self.inputs.in_metadata['M0']) != int :
            m0file=os.path.abspath(self.inputs.bids_dir+'/'+self.inputs.in_metadata['M0'])
            reffile = gen_reference(m0file,newpath=runtime.cwd)
            m0file = reffile
        if self.inputs.in_metadata['M0'] == "True":
            modata2 = dataasl[:, :, :, m0list]
            m0filename=fname_presuffix(self.inputs.in_file,
                                                    suffix='_mofile', newpath=os.getcwd())
            m0obj = nb.Nifti1Image(modata2, allasl.affine, allasl.header)
            m0obj.to_filename(m0filename)
            reffile = gen_reference(m0filename,newpath=runtime.cwd)
            m0file = reffile

        elif type(self.inputs.in_metadata['M0']) == int or  type(self.inputs.in_metadata['M0']) == float :
            m0num=np.float(self.inputs.in_metadata['M0'])
            modata2 = dataasl[:, :, :, deltamlist]
            m0filename=fname_presuffix(self.inputs.in_file,
                                                    suffix='_mofile', newpath=os.getcwd())
            m0obj = nb.Nifti1Image(modata2, allasl.affine, allasl.header)
            m0obj.to_filename(m0filename)
            reffile = gen_reference(m0filename,newpath=runtime.cwd)
            m0file_data=m0num * np.ones_like(nb.load(reffile).get_fdata())

            m0filename1=fname_presuffix(self.inputs.in_file,
                                                    suffix='_mofile', newpath=os.getcwd())
            m0obj1 = nb.Nifti1Image(m0file_data, allasl.affine, allasl.header)
            m0obj1.to_filename(m0filename1)
            m0file = gen_reference(m0filename1,newpath=runtime.cwd)
            reffile = m0file

        elif len(cbflist) > 0 :
            reffile=gen_reference(self.inputs.in_file,newpath=runtime.cwd)
            m0file = reffile 
        
        self._results['ref_file']=reffile
        self._results['m0_file']=m0file
        self.inputs.ref_file = os.path.abspath(self._results['ref_file'])
        self.inputs.m0_file = os.path.abspath(self._results['m0_file'])
        return runtime


def readjson(jsonfile):
    import json
    with open(jsonfile) as f:
        data = json.load(f)
    return data


def init_fsl_gebbr_wf(use_bbr, asl2t1w_dof, asl2t1w_init, sloppy=False, name='fsl_bbr_wf'):
    """
    

    """
    from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from ...niworkflows.utils.images import dseg_label as _dseg_label
    from ...niworkflows.interfaces.freesurfer import PatchedLTAConvert as LTAConvert
    from ...niworkflows.interfaces.registration import FLIRTRPT
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The ASL reference was then co-registered to the T1w reference using
`flirt` [FSL {fsl_ver}, @flirt] with the boundary-based registration [@bbr]
cost-function.
Co-registration was configured with nine degrees of freedom to account
for distortions remaining in the ASL reference.
""".format(fsl_ver=FLIRTRPT().version or '<ver>')

    inputnode = pe.Node(
        niu.IdentityInterface([
            'in_file',
            't1w_dseg', 't1w_brain']),  # FLIRT BBR
        name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(['itk_asl_to_t1', 'itk_t1_to_asl', 'out_report', 'fallback']),
        name='outputnode')

    wm_mask = pe.Node(niu.Function(function=_dseg_label), name='wm_mask')
    wm_mask.inputs.label = 2  # BIDS default is WM=2
    flt_bbr_init = pe.Node(FLIRTRPT(dof=6, generate_report=not use_bbr,
                                    uses_qform=True), name='flt_bbr_init')

    if asl2t1w_init not in ("register", "header"):
        raise ValueError(f"Unknown ASL-T1w initialization option: {asl2t1w_init}")

    if asl2t1w_init == "header":
        raise NotImplementedError("Header-based registration initialization not supported for FSL")

    invt_bbr = pe.Node(fsl.ConvertXFM(invert_xfm=True), name='invt_bbr',
                       mem_gb=DEFAULT_MEMORY_MIN_GB)

    # ASL to T1 transform matrix is from fsl, using c3 tools to convert to
    # something ANTs will like.
    fsl2itk_fwd = pe.Node(c3.C3dAffineTool(fsl2ras=True, itk_transform=True),
                          name='fsl2itk_fwd', mem_gb=DEFAULT_MEMORY_MIN_GB)
    fsl2itk_inv = pe.Node(c3.C3dAffineTool(fsl2ras=True, itk_transform=True),
                          name='fsl2itk_inv', mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow.connect([
        (inputnode, flt_bbr_init, [('in_file', 'in_file'),
                                   ('t1w_brain', 'reference')]),
        (inputnode, fsl2itk_fwd, [('t1w_brain', 'reference_file'),
                                  ('in_file', 'source_file')]),
        (inputnode, fsl2itk_inv, [('in_file', 'reference_file'),
                                  ('t1w_brain', 'source_file')]),
        (invt_bbr, fsl2itk_inv, [('out_file', 'transform_file')]),
        (fsl2itk_fwd, outputnode, [('itk_transform', 'itk_asl_to_t1')]),
        (fsl2itk_inv, outputnode, [('itk_transform', 'itk_t1_to_asl')]),
    ])

    outputnode.inputs.fallback = True

    # Short-circuit workflow building, use rigid registration
    if use_bbr is False:
        workflow.connect([
            (flt_bbr_init, invt_bbr, [('out_matrix_file', 'in_file')]),
            (flt_bbr_init, fsl2itk_fwd, [('out_matrix_file', 'transform_file')]),
            (flt_bbr_init, outputnode, [('out_report', 'out_report')]),
        ])
        outputnode.inputs.fallback = True

        return workflow

    flt_bbr = pe.Node(
        FLIRTRPT(cost_func='bbr', dof=asl2t1w_dof, generate_report=True),
        name='flt_bbr')

    FSLDIR = os.getenv('FSLDIR')
    if FSLDIR:
        flt_bbr.inputs.schedule = op.join(FSLDIR, 'etc/flirtsch/bbr.sch')
    else:
        # Should mostly be hit while building docs
        LOGGER.warning("FSLDIR unset - using packaged BBR schedule")
        flt_bbr.inputs.schedule = pkgr.resource_filename('aslprep', 'data/flirtsch/bbr.sch')

    workflow.connect([
        (inputnode, wm_mask, [('t1w_dseg', 'in_seg')]),
        (inputnode, flt_bbr, [('in_file', 'in_file')]),
        (flt_bbr_init, flt_bbr, [('out_matrix_file', 'in_matrix_file')]),
    ])

    if sloppy is True:
        downsample = pe.Node(niu.Function(
            function=_conditional_downsampling, output_names=["out_file", "out_mask"]),
            name='downsample')
        workflow.connect([
            (inputnode, downsample, [("t1w_brain", "in_file")]),
            (wm_mask, downsample, [("out", "in_mask")]),
            (downsample, flt_bbr, [('out_file', 'reference'),
                                   ('out_mask', 'wm_seg')]),
        ])
    else:
        workflow.connect([
            (inputnode, flt_bbr, [('t1w_brain', 'reference')]),
            (wm_mask, flt_bbr, [('out', 'wm_seg')]),
        ])

    # Short-circuit workflow building, use boundary-based registration
    if use_bbr is True:
        workflow.connect([
            (flt_bbr, invt_bbr, [('out_matrix_file', 'in_file')]),
            (flt_bbr, fsl2itk_fwd, [('out_matrix_file', 'transform_file')]),
            (flt_bbr, outputnode, [('out_report', 'out_report')]),
        ])
        outputnode.inputs.fallback = False

        return workflow

def _conditional_downsampling(in_file, in_mask, zoom_th=4.0):
    """Downsamples the input dataset for sloppy mode."""
    from pathlib import Path
    import numpy as np
    import nibabel as nb
    import nitransforms as nt
    from scipy.ndimage.filters import gaussian_filter

    img = nb.load(in_file)

    zooms = np.array(img.header.get_zooms()[:3])
    if not np.any(zooms < zoom_th):
        return in_file, in_mask

    out_file = Path('desc-resampled_input.nii.gz').absolute()
    out_mask = Path('desc-resampled_mask.nii.gz').absolute()

    shape = np.array(img.shape[:3])
    scaling = zoom_th / zooms
    newrot = np.diag(scaling).dot(img.affine[:3, :3])
    newshape = np.ceil(shape / scaling).astype(int)
    old_center = img.affine.dot(np.hstack((0.5 * (shape - 1), 1.0)))[:3]
    offset = old_center - newrot.dot((newshape - 1) * 0.5)
    newaffine = nb.affines.from_matvec(newrot, offset)

    newref = nb.Nifti1Image(np.zeros(newshape, dtype=np.uint8), newaffine)
    nt.Affine(reference=newref).apply(img).to_filename(out_file)

    mask = nb.load(in_mask)
    mask.set_data_dtype(float)
    mdata = gaussian_filter(mask.get_fdata(dtype=float), scaling)
    floatmask = nb.Nifti1Image(mdata, mask.affine, mask.header)
    newmask = nt.Affine(reference=newref).apply(floatmask)
    hdr = newmask.header.copy()
    hdr.set_data_dtype(np.uint8)
    newmaskdata = (newmask.get_fdata(dtype=float) > 1.5).astype(np.uint8)
    nb.Nifti1Image(newmaskdata, newmask.affine, hdr).to_filename(out_mask)

    return str(out_file), str(out_mask)