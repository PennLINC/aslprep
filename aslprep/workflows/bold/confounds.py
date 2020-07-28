# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Calculate BOLD confounds
^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_bold_confs_wf

"""
from os import getenv
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl
from nipype.algorithms import confounds as nac

from templateflow.api import get as get_template
from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
from ...niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from ...niworkflows.interfaces.utility import KeySelect

from ...niworkflows.interfaces.segmentation import ICA_AROMARPT
from ...niworkflows.interfaces.utils import (AddTSVHeader, TSV2JSON)

from ...config import DEFAULT_MEMORY_MIN_GB
from ...interfaces import (
    GatherConfounds, ICAConfounds,
    FMRISummary, DerivativesDataSink
)


def init_bold_confs_wf(
    mem_gb,
    metadata,
    name="bold_confs_wf",
):
    """
    Build a workflow to generate and write out confounding signals.

    This workflow calculates confounds for a BOLD series, and aggregates them
    into a :abbr:`TSV (tab-separated value)` file, for use as nuisance
    regressors in a :abbr:`GLM (general linear model)`.
    The following confounds are calculated, with column headings in parentheses:

    #. DVARS - original and standardized variants (``dvars``, ``std_dvars``)
    #. Framewise displacement, based on head-motion parameters
       (``framewise_displacement``)
    #. Non-steady-state volumes (``non_steady_state_XX``)
    #. Estimated head-motion parameters, in mm and rad
       (``trans_x``, ``trans_y``, ``trans_z``, ``rot_x``, ``rot_y``, ``rot_z``)

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.bold.confounds import init_bold_confs_wf
            wf = init_bold_confs_wf(
                mem_gb=1,
                metadata={},
                regressors_all_comps=False,
                regressors_dvars_th=1.5,
                regressors_fd_th=0.5,
            )

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    metadata : :obj:`dict`
        BIDS metadata for BOLD file
    name : :obj:`str`
        Name of workflow (default: ``bold_confs_wf``)
    regressors_all_comps : :obj:`bool`
        Indicates whether CompCor decompositions should return all
        components instead of the minimal number of components necessary
        to explain 50 percent of the variance in the decomposition mask.
    regressors_dvars_th : :obj:`float`
        Criterion for flagging DVARS outliers
    regressors_fd_th : :obj:`float`
        Criterion for flagging framewise displacement outliers

    Inputs
    ------
    bold
        BOLD image, after the prescribed corrections (STC, HMC and SDC)
        when available.
    bold_mask
        BOLD series mask
    movpar_file
        SPM-formatted motion parameters file
    skip_vols
        number of non steady state volumes
    t1w_mask
        Mask of the skull-stripped template image
    t1w_tpms
        List of tissue probability maps in T1w space
    t1_bold_xform
        Affine matrix that maps the T1w space into alignment with
        the native BOLD space

    Outputs
    -------
    confounds_file
        TSV of all aggregated confounds
    rois_report
        Reportlet visualizing white-matter/CSF mask used for aCompCor,
        the ROI for tCompCor and the BOLD brain mask.
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
        fields=['bold', 'bold_mask', 'movpar_file', 'skip_vols',
                't1w_mask', 't1w_tpms', 't1_bold_xform']),
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

    workflow.connect([
        # connect inputnode to each non-anatomical confound node
        (inputnode, dvars, [('bold', 'in_file'),
                            ('bold_mask', 'in_mask')]),
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


def init_carpetplot_wf(mem_gb, metadata, name="bold_carpet_wf"):
    """
    Build a workflow to generate *carpet* plots.

    Resamples the MNI parcellation (ad-hoc parcellation derived from the
    Harvard-Oxford template and others).

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    metadata : :obj:`dict`
        BIDS metadata for BOLD file
    name : :obj:`str`
        Name of workflow (default: ``bold_carpet_wf``)

    Inputs
    ------
    bold
        BOLD image, after the prescribed corrections (STC, HMC and SDC)
        when available.
    bold_mask
        BOLD series mask
    confounds_file
        TSV of all aggregated confounds
    t1_bold_xform
        Affine matrix that maps the T1w space into alignment with
        the native BOLD space
    std2anat_xfm
        ANTs-compatible affine-and-warp transform file

    Outputs
    -------
    out_carpetplot
        Path of the generated SVG file

    """
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold', 'bold_mask', 'confounds_file',
                't1_bold_xform', 'std2anat_xfm']),
        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_carpetplot']), name='outputnode')

    # List transforms
    mrg_xfms = pe.Node(niu.Merge(2), name='mrg_xfms')

    # Warp segmentation into EPI space
    seg= get_template(
            'MNI152NLin2009cAsym', resolution=1, desc='carpet',
            suffix='dseg')
    resample_parc = pe.Node(ApplyTransforms(
        float=True,
        input_image=str(seg),
        dimension=3, default_value=0, interpolation='MultiLabel'),
        name='resample_parc')

    # Carpetplot and confounds plot
    conf_plot = pe.Node(FMRISummary(
        tr=metadata['RepetitionTime'],
        confounds_list=[
            ('std_dvars', None, 'DVARS'),
            ('framewise_displacement', 'mm', 'FD')]),
        name='conf_plot', mem_gb=mem_gb)
    ds_report_bold_conf = pe.Node(
        DerivativesDataSink(desc='carpetplot', datatype="figures", keep_dtype=True),
        name='ds_report_bold_conf', run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow = Workflow(name=name)
    workflow.connect([
        (inputnode, mrg_xfms, [('t1_bold_xform', 'in1'),
                               ('std2anat_xfm', 'in2')]),
        (inputnode, resample_parc, [('bold_mask', 'reference_image')]),
        (mrg_xfms, resample_parc, [('out', 'transforms')]),
        # Carpetplot
        (inputnode, conf_plot, [
            ('bold', 'in_func'),
            ('bold_mask', 'in_mask'),
            ('confounds_file', 'confounds_file')]),
        (resample_parc, conf_plot, [('output_image', 'in_segm')]),
        (conf_plot, ds_report_bold_conf, [('out_file', 'in_file')]),
        (conf_plot, outputnode, [('out_file', 'out_carpetplot')]),
    ])
    return workflow

def _remove_volumes(bold_file, skip_vols):
    """Remove skip_vols from bold_file."""
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix

    if skip_vols == 0:
        return bold_file

    out = fname_presuffix(bold_file, suffix='_cut')
    bold_img = nb.load(bold_file)
    bold_img.__class__(bold_img.dataobj[..., skip_vols:],
                       bold_img.affine, bold_img.header).to_filename(out)

    return out


def _add_volumes(bold_file, bold_cut_file, skip_vols):
    """Prepend skip_vols from bold_file onto bold_cut_file."""
    import nibabel as nb
    import numpy as np
    from nipype.utils.filemanip import fname_presuffix

    if skip_vols == 0:
        return bold_cut_file

    bold_img = nb.load(bold_file)
    bold_cut_img = nb.load(bold_cut_file)

    bold_data = np.concatenate((bold_img.dataobj[..., :skip_vols],
                                bold_cut_img.dataobj), axis=3)

    out = fname_presuffix(bold_cut_file, suffix='_addnonsteady')
    bold_img.__class__(bold_data, bold_img.affine, bold_img.header).to_filename(out)

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

    out = fname_presuffix(roi_file, suffix='_boldmsk')
    roi.__class__(roidata, roi.affine, roi.header).to_filename(out)
    return out
