# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Calculate asl confounds
^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_asl_confs_wf


"""
from os import getenv
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl
from nipype.algorithms import confounds as nac

from templateflow.api import get as get_template
from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
from ...niworkflows.interfaces.confounds import ExpandModel, SpikeRegressors
from ...niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from ...niworkflows.interfaces.images import SignalExtraction
from ...niworkflows.interfaces.masks import ROIsPlot
from ...niworkflows.interfaces.utility import KeySelect
from ...niworkflows.interfaces.patches import (
    RobustACompCor as ACompCor,
    RobustTCompCor as TCompCor,
)
from ...niworkflows.interfaces.plotting import (
    CompCorVariancePlot, ConfoundsCorrelationPlot
)
from ...niworkflows.interfaces.segmentation import ICA_AROMARPT
from ...niworkflows.interfaces.utils import (
    TPM2ROI, AddTPMs, AddTSVHeader, TSV2JSON, DictMerge
)

from ...config import DEFAULT_MEMORY_MIN_GB
from ...interfaces import (
    GatherConfounds, 
    ASLSummary, DerivativesDataSink
)


def init_asl_confs_wf(
    mem_gb,
    metadata,
    name="asl_confs_wf",
):
    """
    Build a workflow to generate and write out confounding signals.

    This workflow calculates confounds for a asl series, and aggregates them
    into a :abbr:`TSV (tab-separated value)` file, for use as nuisance
    regressors in a :abbr:`GLM (general linear model)`.
    The following confounds are calculated, with column headings in parentheses:

    #. DVARS - original and standardized variants (``dvars``, ``std_dvars``)
    #. Framewise displacement, based on head-motion parameters
       (``framewise_displacement``)
    #. Estimated head-motion parameters, in mm and rad
       (``trans_x``, ``trans_y``, ``trans_z``, ``rot_x``, ``rot_y``, ``rot_z``)

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.confounds import init_asl_confs_wf
            wf = init_asl_confs_wf(
                mem_gb=1,
                metadata={},
            )

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of asl file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    metadata : :obj:`dict`
        BIDS metadata for asl file
    name : :obj:`str`
        Name of workflow (default: ``asl_confs_wf``)


    Inputs
    ------
    asl
        asl image, after the prescribed corrections (STC, HMC and SDC)
        when available.
    asl_mask
        asl series mask
    movpar_file
        SPM-formatted motion parameters file
    skip_vols
        number of non steady state volumes
    t1w_mask
        Mask of the skull-stripped template image
    t1w_tpms
        List of tissue probability maps in T1w space
    t1_asl_xform
        Affine matrix that maps the T1w space into alignment with
        the native asl space

    Outputs
    -------
    confounds_file
        TSV of all aggregated confounds
    confounds_metadata
        Confounds metadata dictionary.

    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
Several confounding time-series were calculated based on the
*preprocessed ASL*: framewise displacement (FD) and DVARS. 
FD and DVARS are calculated for each ASL run, both using their
implementations in *Nipype* [following the definitions by @power_fd_dvars].
The head-motion estimates calculated in the correction step were also
placed within the corresponding confounds file.

"""
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['asl', 'asl_mask', 'movpar_file', 'skip_vols',
                't1w_mask', 't1w_tpms', 't1_asl_xform']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['confounds_file', 'confounds_metadata']),
        name='outputnode')

    # DVARS
    dvars = pe.Node(nac.ComputeDVARS(save_nstd=True, save_std=True, remove_zerovariance=True),
                    name="dvars", mem_gb=mem_gb)

    # Frame displacement
    fdisp = pe.Node(nac.FramewiseDisplacement(parameter_source="SPM"),
                    name="fdisp", mem_gb=mem_gb)


    # Global and segment regressors
    #signals_class_labels = ["csf", "white_matter", "global_signal"]

    # Arrange confounds
    add_dvars_header = pe.Node(
        AddTSVHeader(columns=["dvars"]),
        name="add_dvars_header", mem_gb=0.01, run_without_submitting=True)
    add_std_dvars_header = pe.Node(
        AddTSVHeader(columns=["std_dvars"]),
        name="add_std_dvars_header", mem_gb=0.01, run_without_submitting=True)
    add_motion_headers = pe.Node(
        AddTSVHeader(columns=["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"]),
        name="add_motion_headers", mem_gb=0.01, run_without_submitting=True)
    concat = pe.Node(GatherConfounds(), name="concat", mem_gb=0.01, run_without_submitting=True)



    # Expand model to include derivatives and quadratics

    workflow.connect([
        # connect inputnode to each non-anatomical confound node
        (inputnode, dvars, [('asl', 'in_file'),
                            ('asl_mask', 'in_mask')]),
        (inputnode, fdisp, [('movpar_file', 'in_file')]),
        # Collate computed confounds together
        (inputnode, add_motion_headers, [('movpar_file', 'in_file')]),
        (dvars, add_dvars_header, [('out_nstd', 'in_file')]),
        (dvars, add_std_dvars_header, [('out_std', 'in_file')]),
        (fdisp, concat, [('out_file', 'fd')]),
        (add_motion_headers, concat, [('out_file', 'motion')]),
        (add_dvars_header, concat, [('out_file', 'dvars')]),
        (add_std_dvars_header, concat, [('out_file', 'std_dvars')]),

        # Expand the model with derivatives, quadratics, and spikes

        # Set outputs
        (concat, outputnode, [('confounds_file', 'confounds_file')]),
    
    ])

    return workflow


def init_carpetplot_wf(mem_gb, metadata, name="asl_carpet_wf"):
    """
    Build a workflow to generate *carpet* plots.

    Resamples the MNI parcellation (ad-hoc parcellation derived from the
    Harvard-Oxford template and others).

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of ASL file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    metadata : :obj:`dict`
        BIDS metadata for ASL file
    name : :obj:`str`
        Name of workflow (default: ``asl_carpet_wf``)

    Inputs
    ------
    asl
        asl image, after the prescribed corrections (STC, HMC and SDC)
        when available.
    asl_mask
        ASL series mask
    confounds_file
        TSV of all aggregated confounds
    t1_asl_xform
        Affine matrix that maps the T1w space into alignment with
        the native ASL space
    std2anat_xfm
        ANTs-compatible affine-and-warp transform file

    Outputs
    -------
    out_carpetplot
        Path of the generated SVG file

    """
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['asl', 'asl_mask', 'confounds_file',
                't1_asl_xform', 'std2anat_xfm']),
        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_carpetplot']), name='outputnode')

    # List transforms
    mrg_xfms = pe.Node(niu.Merge(2), name='mrg_xfms')

    # Warp segmentation into EPI space
    resample_parc = pe.Node(ApplyTransforms(
        float=True,
        input_image=str(get_template(
            'MNI152NLin2009cAsym', resolution=1, desc='carpet',
            suffix='dseg', extension=['.nii', '.nii.gz'])),
        dimension=3, default_value=0, interpolation='MultiLabel'),
        name='resample_parc')

    # Carpetplot and confounds plot
    conf_plot = pe.Node(ASLSummary(
        tr=metadata['RepetitionTime'],
        confounds_list=[
            ('std_dvars', None, 'DVARS'),
            ('framewise_displacement', 'mm', 'FD')]),
        name='conf_plot', mem_gb=mem_gb)
    ds_report_asl_conf = pe.Node(
        DerivativesDataSink(desc='carpetplot',datatype="figures", 
        keep_dtype=True),
         name='ds_report_asl_conf', run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow = Workflow(name=name)
    workflow.connect([
        (inputnode, mrg_xfms, [('t1_asl_xform', 'in1'),
                               ('std2anat_xfm', 'in2')]),
        (inputnode, resample_parc, [('asl_mask', 'reference_image')]),
        (mrg_xfms, resample_parc, [('out', 'transforms')]),
        # Carpetplot
        (inputnode, conf_plot, [
            ('asl', 'in_func'),
            ('asl_mask', 'in_mask'),
            ('confounds_file', 'confounds_file')]),
        (resample_parc, conf_plot, [('output_image', 'in_segm')]),
        (conf_plot, ds_report_asl_conf, [('out_file', 'in_file')]),
        (conf_plot, outputnode, [('out_file', 'out_carpetplot')]),
    ])
    return workflow


def _remove_volumes(asl_file, skip_vols):
    """Remove skip_vols from asl_file."""
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix

    if skip_vols == 0:
        return asl_file

    out = fname_presuffix(asl_file, suffix='_cut')
    asl_img = nb.load(asl_file)
    asl_img.__class__(asl_img.dataobj[..., skip_vols:],
                       asl_img.affine, asl_img.header).to_filename(out)

    return out


def _add_volumes(asl_file, asl_cut_file, skip_vols):
    """Prepend skip_vols from asl_file onto asl_cut_file."""
    import nibabel as nb
    import numpy as np
    from nipype.utils.filemanip import fname_presuffix

    if skip_vols == 0:
        return asl_cut_file

    asl_img = nb.load(asl_file)
    asl_cut_img = nb.load(asl_cut_file)

    asl_data = np.concatenate((asl_img.dataobj[..., :skip_vols],
                                asl_cut_img.dataobj), axis=3)

    out = fname_presuffix(asl_cut_file, suffix='_addnonsteady')
    asl_img.__class__(asl_data, asl_img.affine, asl_img.header).to_filename(out)

    return out


def _maskroi(in_mask, roi_file):
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix

    roi = nb.load(roi_file)
    roidata = roi.get_data().astype(np.uint8)
    msk = nb.load(in_mask).get_data().astype(bool)
    roidata[~msk] = 0
    roi.set_data_dtype(np.uint8)

    out = fname_presuffix(roi_file, suffix='_aslmsk')
    roi.__class__(roidata, roi.affine, roi.header).to_filename(out)
    return out