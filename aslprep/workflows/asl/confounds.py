# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for calculating confounds for ASL data."""

from fmriprep.workflows.bold.confounds import _carpet_parcellation
from nipype.algorithms import confounds as nac
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.utility import AddTSVHeader
from templateflow.api import get as get_template

from aslprep import config
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.bids import DerivativesDataSink
from aslprep.interfaces.confounds import ComputeCBFQC, GatherConfounds
from aslprep.interfaces.plotting import ASLCarpetPlot


def init_asl_confounds_wf(
    n_volumes: int,
    mem_gb: float,
    freesurfer: bool = False,
    name: str = 'asl_confounds_wf',
):
    """Build a workflow to generate and write out confounding signals.

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

            from aslprep.tests.tests import mock_config
            from aslprep.workflows.asl.confounds import init_asl_confounds_wf

            with mock_config():
                wf = init_asl_confounds_wf(
                    n_volumes=50,
                    mem_gb=1,
                    freesurfer=True,
                )

    Parameters
    ----------
    n_volumes : :obj:`int`
        Number of volumes in the ASL file.
        Some relative measures (e.g., FD, DVARS, RMSD) will produce empty arrays for single-volume
        datasets, so we must skip those.
    mem_gb : :obj:`float`
        Size of asl file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    freesurfer : :obj:`bool`
        True if FreeSurfer derivatives were used.
    name : :obj:`str`
        Name of workflow (default: ``asl_confounds_wf``)

    Inputs
    ------
    asl
        asl image, after the prescribed corrections (HMC and SDC)
        when available.
    asl_mask
        asl series mask
    movpar_file
        SPM-formatted motion parameters file
    rmsd_file
        Framewise displacement as measured by ``fsl_motion_outliers``.
    skip_vols
        number of non steady state volumes
    t1w_mask
        Mask of the skull-stripped template image
    t1w_tpms
        List of tissue probability maps in T1w space
    aslref2anat_xfm
        Affine matrix that maps the T1w space into alignment with
        the native asl space

    Outputs
    -------
    confounds_file
        TSV of all aggregated confounds
    confounds_metadata
        Confounds metadata dictionary.
    crown_mask
        Mask of brain edge voxels
    acompcor_masks
    """
    from fmriprep.interfaces.confounds import aCompCorMasks
    from niworkflows.interfaces.images import SignalExtraction
    from niworkflows.interfaces.morphology import BinaryDilation, BinarySubtraction
    from niworkflows.interfaces.nibabel import ApplyMask, Binarize

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
Several confounding timeseries were calculated, including both framewise displacement
(FD) and DVARS. FD and DVARS are calculated using the implementations in in *Nipype*
(following the definition by [@power_fd_dvars]) for each ASL run.  ASLPrep summarizes
in-scanner motion as the mean framewise displacement and relative root-mean square displacement.
"""
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'asl',
                'asl_mask',
                'movpar_file',
                'rmsd_file',
                'skip_vols',
                't1w_mask',
                't1w_tpms',
                'aslref2anat_xfm',
            ],
        ),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'confounds_file',
                'confounds_metadata',
                'acompcor_masks',
                'crown_mask',
            ],
        ),
        name='outputnode',
    )

    add_motion_headers = pe.Node(
        AddTSVHeader(columns=['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']),
        name='add_motion_headers',
        mem_gb=0.01,
        run_without_submitting=True,
    )
    workflow.connect([(inputnode, add_motion_headers, [('movpar_file', 'in_file')])])

    if n_volumes > 2:  # set to 2 bc relative arrays will be 1D instead of 2D for 1-volume data
        # DVARS
        dvars = pe.Node(
            nac.ComputeDVARS(save_nstd=True, save_std=True, remove_zerovariance=True),
            name='dvars',
            mem_gb=mem_gb,
        )
        workflow.connect([
            (inputnode, dvars, [
                ('asl', 'in_file'),
                ('asl_mask', 'in_mask'),
            ]),
        ])  # fmt:skip

        # Frame displacement
        fdisp = pe.Node(
            nac.FramewiseDisplacement(parameter_source='SPM'),
            name='fdisp',
            mem_gb=mem_gb,
        )
        workflow.connect([(inputnode, fdisp, [('movpar_file', 'in_file')])])

        # Arrange confounds
        add_dvars_header = pe.Node(
            AddTSVHeader(columns=['dvars']),
            name='add_dvars_header',
            mem_gb=0.01,
            run_without_submitting=True,
        )
        add_std_dvars_header = pe.Node(
            AddTSVHeader(columns=['std_dvars']),
            name='add_std_dvars_header',
            mem_gb=0.01,
            run_without_submitting=True,
        )
        add_rmsd_header = pe.Node(
            AddTSVHeader(columns=['rmsd']),
            name='add_rmsd_header',
            mem_gb=0.01,
            run_without_submitting=True,
        )

        workflow.connect([
            # Collate computed confounds together
            (inputnode, add_rmsd_header, [('rmsd_file', 'in_file')]),
            (dvars, add_dvars_header, [('out_nstd', 'in_file')]),
            (dvars, add_std_dvars_header, [('out_std', 'in_file')]),
        ])  # fmt:skip

    # Project T1w mask into BOLD space and merge with BOLD brainmask
    t1w_mask_tfm = pe.Node(
        ApplyTransforms(interpolation='GenericLabel', invert_transform_flags=[True], args='-v'),
        name='t1w_mask_tfm',
    )
    union_mask = pe.Node(niu.Function(function=_binary_union), name='union_mask')

    # Create the crown mask
    dilated_mask = pe.Node(BinaryDilation(), name='dilated_mask')
    subtract_mask = pe.Node(BinarySubtraction(), name='subtract_mask')

    workflow.connect([
        # Brain mask
        (inputnode, t1w_mask_tfm, [
            ('t1w_mask', 'input_image'),
            ('asl_mask', 'reference_image'),
            ('aslref2anat_xfm', 'transforms'),
        ]),
        (inputnode, union_mask, [('asl_mask', 'mask1')]),
        (t1w_mask_tfm, union_mask, [('output_image', 'mask2')]),
        (union_mask, dilated_mask, [('out', 'in_mask')]),
        (union_mask, subtract_mask, [('out', 'in_subtract')]),
        (dilated_mask, subtract_mask, [('out_mask', 'in_base')]),
        (subtract_mask, outputnode, [('out_mask', 'crown_mask')]),
    ])  # fmt:skip

    # Generate aCompCor probseg maps
    acc_masks = pe.Node(aCompCorMasks(is_aseg=freesurfer), name='acc_masks')
    workflow.connect([
        (inputnode, acc_masks, [
            ('t1w_tpms', 'in_vfs'),
            (('asl', _get_zooms), 'bold_zooms'),
        ]),
    ])  # fmt:skip

    # Resample probseg maps in BOLD space via T1w-to-BOLD transform
    acc_msk_tfm = pe.MapNode(
        ApplyTransforms(interpolation='Gaussian', invert_transform_flags=[True], args='-v'),
        iterfield=['input_image'],
        name='acc_msk_tfm',
        mem_gb=0.1,
    )
    workflow.connect([
        (inputnode, acc_msk_tfm, [
            ('aslref2anat_xfm', 'transforms'),
            ('asl_mask', 'reference_image'),
        ]),
        (acc_masks, acc_msk_tfm, [('out_masks', 'input_image')]),
    ])  # fmt:skip

    acc_msk_brain = pe.MapNode(ApplyMask(), name='acc_msk_brain', iterfield=['in_file'])
    workflow.connect([
        (inputnode, acc_msk_brain, [('asl_mask', 'in_mask')]),
        (acc_msk_tfm, acc_msk_brain, [('output_image', 'in_file')]),
    ])  # fmt:skip

    acc_msk_bin = pe.MapNode(Binarize(thresh_low=0.99), name='acc_msk_bin', iterfield=['in_file'])
    workflow.connect([
        (acc_msk_brain, acc_msk_bin, [('out_file', 'in_file')]),
        (acc_msk_bin, outputnode, [('out_file', 'acompcor_masks')]),
    ])  # fmt:skip

    # Global and segment regressors
    signals_class_labels = [
        'global_signal',
        'csf',
        'white_matter',
        'csf_wm',
    ]
    merge_rois = pe.Node(
        niu.Merge(2, ravel_inputs=True),
        name='merge_rois',
        run_without_submitting=True,
    )
    signals = pe.Node(
        SignalExtraction(class_labels=signals_class_labels),
        name='signals',
        mem_gb=mem_gb,
    )
    workflow.connect([
        (inputnode, merge_rois, [('asl_mask', 'in1')]),
        (acc_msk_bin, merge_rois, [('out_file', 'in2')]),
        (inputnode, signals, [('asl', 'in_file')]),
        (merge_rois, signals, [('out', 'label_files')]),
    ])  # fmt:skip

    concat = pe.Node(
        GatherConfounds(),
        name='concat',
        mem_gb=0.01,
        run_without_submitting=True,
    )
    workflow.connect([
        (add_motion_headers, concat, [('out_file', 'motion')]),
        (signals, concat, [('out_file', 'signals')]),
        (concat, outputnode, [('confounds_file', 'confounds_file')]),
    ])  # fmt:skip
    if n_volumes > 2:
        workflow.connect([
            (fdisp, concat, [('out_file', 'fd')]),
            (add_rmsd_header, concat, [('out_file', 'rmsd')]),
            (add_dvars_header, concat, [('out_file', 'dvars')]),
            (add_std_dvars_header, concat, [('out_file', 'std_dvars')]),
        ])  # fmt:skip

    return workflow


def init_carpetplot_wf(
    mem_gb: float,
    confounds_list: list,
    metadata: dict,
    cifti_output: bool,
    suffix: str = 'asl',
    name: str = 'asl_carpet_wf',
):
    """Build a workflow to generate *carpet* plots.

    Resamples the MNI parcellation (ad-hoc parcellation derived from the
    Harvard-Oxford template and others).

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB - please note that this size
        should be calculated after resamplings that may extend
        the FoV
    confounds_list : :obj:`list` of length-3 :obj:`tuple` of :obj:`str`
    metadata : :obj:`dict`
        BIDS metadata for BOLD file
    name : :obj:`str`
        Name of workflow (default: ``asl_carpet_wf``)

    Inputs
    ------
    asl
        BOLD image, after the prescribed corrections (STC, HMC and SDC)
        when available.
    asl_mask
        BOLD series mask
    confounds_file
        TSV of all aggregated confounds
    aslref2anat_xfm
        Affine matrix that maps the T1w space into alignment with
        the native BOLD space
    std2anat_xfm
        ANTs-compatible affine-and-warp transform file
    cifti_asl
        BOLD image in CIFTI format, to be used in place of volumetric BOLD
    crown_mask
        Mask of brain edge voxels
    acompcor_mask
        Mask of deep WM+CSF
    dummy_scans
        Number of nonsteady states to be dropped at the beginning of the timeseries.

    Outputs
    -------
    out_carpetplot
        Path of the generated SVG file

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    from aslprep.interfaces.ants import ApplyTransforms

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'asl',
                'asl_mask',
                'confounds_file',
                'aslref2anat_xfm',
                'std2anat_xfm',
                'cifti_asl',
                'crown_mask',
                'acompcor_mask',
                'dummy_scans',
            ]
        ),
        name='inputnode',
    )

    outputnode = pe.Node(niu.IdentityInterface(fields=['out_carpetplot']), name='outputnode')

    # Carpetplot and confounds plot
    conf_plot = pe.Node(
        ASLCarpetPlot(
            tr=metadata['RepetitionTime'],
            confounds_list=confounds_list,
        ),
        name='conf_plot',
        mem_gb=mem_gb,
    )
    ds_report_asl_conf = pe.Node(
        DerivativesDataSink(
            desc='carpetplot',
            datatype='figures',
            suffix=suffix,
            extension='svg',
            dismiss_entities=('echo',),
        ),
        name='ds_report_asl_conf',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    parcels = pe.Node(niu.Function(function=_carpet_parcellation), name='parcels')
    parcels.inputs.nifti = not cifti_output
    # List transforms
    mrg_xfms = pe.Node(niu.Merge(2), name='mrg_xfms')

    # Warp segmentation into EPI space
    resample_parc = pe.Node(
        ApplyTransforms(
            dimension=3,
            input_image=str(
                get_template(
                    'MNI152NLin2009cAsym',
                    resolution=1,
                    desc='carpet',
                    suffix='dseg',
                    extension=['.nii', '.nii.gz'],
                ),
            ),
            interpolation='GenericLabel',
            invert_transform_flags=[True, False],
            args='-u int -v',
        ),
        name='resample_parc',
    )

    workflow = Workflow(name=name)
    if cifti_output:
        workflow.connect(inputnode, 'cifti_asl', conf_plot, 'in_cifti')

    workflow.connect([
        (inputnode, mrg_xfms, [
            ('aslref2anat_xfm', 'in1'),
            ('std2anat_xfm', 'in2'),
        ]),
        (inputnode, resample_parc, [('asl_mask', 'reference_image')]),
        (inputnode, parcels, [
            ('crown_mask', 'crown_mask'),
            ('acompcor_mask', 'acompcor_mask'),
        ]),
        (inputnode, conf_plot, [
            ('asl', 'in_nifti'),
            ('confounds_file', 'confounds_file'),
            ('dummy_scans', 'drop_trs'),
        ]),
        (mrg_xfms, resample_parc, [('out', 'transforms')]),
        (resample_parc, parcels, [('output_image', 'segmentation')]),
        (parcels, conf_plot, [('out', 'in_segm')]),
        (conf_plot, ds_report_asl_conf, [('out_file', 'in_file')]),
        (conf_plot, outputnode, [('out_file', 'out_carpetplot')]),
    ])  # fmt:skip
    return workflow


def init_cbf_confounds_wf(
    scorescrub=False,
    basil=False,
    name='cbf_confounds_wf',
):
    """Create a workflow for :abbr:`dolui2017automated (compute cbf)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep.workflows.asl.confounds import init_cbf_confounds_wf

            with mock_config():
                wf = init_cbf_confounds_wf(
                    scorescrub=True,
                    basil=True,
                    name="cbf_confounds_wf",
                )

    Parameters
    ----------
    scorescrub : bool
    basil : bool
    name : :obj:`str`
        Name of workflow (default: "cbf_qc_wf")

    Inputs
    ------
    *cbf
        all cbf
    asl_mask
        asl mask NIFTI file
    t1w_tpms
        t1w probability maps
    aslref2anat_xfm
        aslref to t1w transformation file
    asl_mask_std : list
        Since ASLPrep always includes MNI152NLin2009cAsym as a standard space,
        this should always be provided.

    Outputs
    -------
    qc_file
        qc measures in tsv
    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """
The quality evaluation index (QEI) was computed for each CBF map [@dolui2017automated].
QEI is based on the similarity between the CBF and the structural images, the spatial
variability of the CBF image, and the percentage of grey matter voxels containing
negative CBF values.
"""
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'name_source',
                'asl_mask',
                't1w_mask',
                't1w_tpms',
                'aslref2anat_xfm',
                'anat2mni2009c_xfm',
                # CBF inputs
                'mean_cbf',
                # SCORE/SCRUB inputs
                'mean_cbf_score',
                'mean_cbf_scrub',
                # BASIL inputs
                'mean_cbf_basil',
                'mean_cbf_gm_basil',
                'mean_cbf_wm_basil',
                # non-GE inputs
                'confounds_file',
                'rmsd_file',
            ],
        ),
        name='inputnode',
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=['qc_file']), name='outputnode')

    def _pick_gm(files):
        return files[0]

    def _pick_wm(files):
        return files[1]

    def _pick_csf(files):
        return files[2]

    gm_tfm = pe.Node(
        ApplyTransforms(
            interpolation='NearestNeighbor',
            float=True,
            invert_transform_flags=[True],
            args='-v',
        ),
        name='gm_tfm',
        mem_gb=0.1,
    )

    workflow.connect([
        (inputnode, gm_tfm, [
            ('asl_mask', 'reference_image'),
            ('aslref2anat_xfm', 'transforms'),
            (('t1w_tpms', _pick_gm), 'input_image'),
        ]),
    ])  # fmt:skip

    wm_tfm = pe.Node(
        ApplyTransforms(
            interpolation='NearestNeighbor',
            float=True,
            invert_transform_flags=[True],
            args='-v',
        ),
        name='wm_tfm',
        mem_gb=0.1,
    )

    workflow.connect([
        (inputnode, wm_tfm, [
            ('asl_mask', 'reference_image'),
            ('aslref2anat_xfm', 'transforms'),
            (('t1w_tpms', _pick_wm), 'input_image'),
        ]),
    ])  # fmt:skip

    csf_tfm = pe.Node(
        ApplyTransforms(
            interpolation='NearestNeighbor',
            float=True,
            invert_transform_flags=[True],
            args='-v',
        ),
        name='csf_tfm',
        mem_gb=0.1,
    )

    workflow.connect([
        (inputnode, csf_tfm, [
            ('asl_mask', 'reference_image'),
            ('aslref2anat_xfm', 'transforms'),
            (('t1w_tpms', _pick_csf), 'input_image'),
        ]),
    ])  # fmt:skip

    warp_t1w_mask_to_aslref = pe.Node(
        ApplyTransforms(
            interpolation='NearestNeighbor',
            float=True,
            invert_transform_flags=[True],
            args='-v',
        ),
        name='warp_t1w_mask_to_aslref',
        mem_gb=0.1,
    )

    workflow.connect([
        (inputnode, warp_t1w_mask_to_aslref, [
            ('asl_mask', 'reference_image'),
            ('aslref2anat_xfm', 'transforms'),
            ('t1w_mask', 'input_image'),
        ]),
    ])  # fmt:skip

    template_brain_mask = str(
        get_template('MNI152NLin2009cAsym', resolution=2, desc='brain', suffix='mask')
    )

    aslref2mni152nlin2009casym = pe.Node(niu.Merge(2), name='aslref2mni152nlin2009casym')
    workflow.connect([
        (inputnode, aslref2mni152nlin2009casym, [
            ('aslref2anat_xfm', 'in1'),
            ('anat2mni2009c_xfm', 'in2'),
        ]),
    ])  # fmt:skip

    warp_asl_mask_to_mni152nlin2009casym = pe.Node(
        ApplyTransforms(
            interpolation='NearestNeighbor',
            float=True,
            reference_image=template_brain_mask,
            args='-v',
        ),
        name='warp_asl_mask_to_mni152nlin2009casym',
        mem_gb=0.1,
    )

    workflow.connect([
        (inputnode, warp_asl_mask_to_mni152nlin2009casym, [('asl_mask', 'input_image')]),
        (aslref2mni152nlin2009casym, warp_asl_mask_to_mni152nlin2009casym, [
            ('out', 'transforms'),
        ])
    ])  # fmt:skip

    compute_qc_metrics = pe.Node(
        ComputeCBFQC(
            tpm_threshold=0.7,
            template_mask=template_brain_mask,
        ),
        name='compute_qc_metrics',
        run_without_submitting=True,
        mem_gb=0.2,
    )

    workflow.connect([
        (warp_t1w_mask_to_aslref, compute_qc_metrics, [('output_image', 't1w_mask')]),
        (inputnode, compute_qc_metrics, [
            ('name_source', 'name_source'),
            ('asl_mask', 'asl_mask'),
            ('mean_cbf', 'mean_cbf'),
            ('confounds_file', 'confounds_file'),
            ('rmsd_file', 'rmsd_file'),
        ]),
        (warp_asl_mask_to_mni152nlin2009casym, compute_qc_metrics, [
            ('output_image', 'asl_mask_std'),
        ]),
        (gm_tfm, compute_qc_metrics, [('output_image', 'gm_tpm')]),
        (wm_tfm, compute_qc_metrics, [('output_image', 'wm_tpm')]),
        (csf_tfm, compute_qc_metrics, [('output_image', 'csf_tpm')]),
        (compute_qc_metrics, outputnode, [('qc_file', 'qc_file')]),
    ])  # fmt:skip

    ds_qc = pe.Node(
        DerivativesDataSink(
            base_directory=config.execution.aslprep_dir,
            desc='qualitycontrol',
            suffix='cbf',
            compress=False,
        ),
        name='ds_qc',
        run_without_submitting=True,
    )

    workflow.connect([
        (inputnode, ds_qc, [('name_source', 'source_file')]),
        (compute_qc_metrics, ds_qc, [('qc_file', 'in_file')]),
    ])  # fmt:skip

    ds_qc_metadata = pe.Node(
        DerivativesDataSink(
            base_directory=config.execution.aslprep_dir,
            dismiss_entities=list(DerivativesDataSink._allowed_entities),
            allowed_entities=[],
            desc='qualitycontrol',
            suffix='cbf',
            extension='.json',
        ),
        name='ds_qc_metadata',
        run_without_submitting=True,
    )

    workflow.connect([
        (inputnode, ds_qc_metadata, [('name_source', 'source_file')]),
        (compute_qc_metrics, ds_qc_metadata, [('qc_metadata', 'in_file')]),
    ])  # fmt:skip

    if scorescrub:
        workflow.connect([
            (inputnode, compute_qc_metrics, [
                ('mean_cbf_scrub', 'mean_cbf_scrub'),
                ('mean_cbf_score', 'mean_cbf_score'),
            ]),
        ])  # fmt:skip

    if basil:
        workflow.connect([
            (inputnode, compute_qc_metrics, [
                ('mean_cbf_basil', 'mean_cbf_basil'),
                ('mean_cbf_gm_basil', 'mean_cbf_gm_basil'),
            ]),
        ])  # fmt:skip

    return workflow


def _binary_union(mask1, mask2):
    """Generate the union of two masks."""
    from pathlib import Path

    import nibabel as nb
    import numpy as np

    img = nb.load(mask1)
    mskarr1 = np.asanyarray(img.dataobj, dtype=int) > 0
    mskarr2 = np.asanyarray(nb.load(mask2).dataobj, dtype=int) > 0
    out = img.__class__(mskarr1 | mskarr2, img.affine, img.header)
    out.set_data_dtype('uint8')
    out_name = Path('mask_union.nii.gz').absolute()
    out.to_filename(out_name)
    return str(out_name)


def _get_zooms(in_file):
    import nibabel as nb

    return tuple(nb.load(in_file).header.get_zooms()[:3])
