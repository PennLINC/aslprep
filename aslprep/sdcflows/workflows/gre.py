# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Processing phase-difference (aka :abbr:`GRE (gradient-recalled echo)`) fieldmaps.

.. _gre-fieldmaps:

Workflows for processing :abbr:`GRE (gradient recalled echo)` fieldmaps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Workflows for preparing the magnitude part of :abbr:`GRE (gradient-recalled echo)` fieldmap
images and cleaning up the fieldmaps created from the phases or phasediff.

"""

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl, ants
from niflow.nipype1.workflows.dmri.fsl.utils import cleanup_edge_pipeline
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.images import IntraModalMerge
from niworkflows.interfaces.masks import BETRPT


def init_magnitude_wf(omp_nthreads, name='magnitude_wf'):
    """
    Prepare the magnitude part of :abbr:`GRE (gradient-recalled echo)` fieldmaps.

    Average (if not done already) the magnitude part of the
    :abbr:`GRE (gradient recalled echo)` images, run N4 to
    correct for B1 field nonuniformity, and skull-strip the
    preprocessed magnitude.

    Workflow Graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from sdcflows.workflows.fmap import init_magnitude_wf
            wf = init_magnitude_wf(omp_nthreads=6)

    Parameters
    ----------
    omp_nthreads : int
        Maximum number of threads an individual process may use
    name : str
        Name of workflow (default: ``prepare_magnitude_w``)

    Inputs
    ------
    magnitude : pathlike
        Path to the corresponding magnitude path(s).

    Outputs
    -------
    fmap_ref : pathlike
        Path to the fieldmap reference calculated in this workflow.
    fmap_mask : pathlike
        Path to a binary brain mask corresponding to the reference above.

    """
    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=['magnitude']), name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['fmap_ref', 'fmap_mask', 'mask_report']),
        name='outputnode')

    # Merge input magnitude images
    # Do not reorient to RAS to preserve the validity of PhaseEncodingDirection
    magmrg = pe.Node(IntraModalMerge(hmc=False, to_ras=False), name='magmrg')

    # de-gradient the fields ("bias/illumination artifact")
    n4_correct = pe.Node(ants.N4BiasFieldCorrection(dimension=3, copy_header=True),
                         name='n4_correct', n_procs=omp_nthreads)
    bet = pe.Node(BETRPT(generate_report=True, frac=0.6, mask=True),
                  name='bet')

    workflow.connect([
        (inputnode, magmrg, [('magnitude', 'in_files')]),
        (magmrg, n4_correct, [('out_avg', 'input_image')]),
        (n4_correct, bet, [('output_image', 'in_file')]),
        (bet, outputnode, [('mask_file', 'fmap_mask'),
                           ('out_file', 'fmap_ref'),
                           ('out_report', 'mask_report')]),
    ])
    return workflow


def init_fmap_postproc_wf(omp_nthreads, fmap_bspline, median_kernel_size=5,
                          name='fmap_postproc_wf'):
    """
    Postprocess a B0 map estimated elsewhere.

    This workflow denoises (mostly via smoothing) a B0 fieldmap.

    Workflow Graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from sdcflows.workflows.fmap import init_fmap_postproc_wf
            wf = init_fmap_postproc_wf(omp_nthreads=6, fmap_bspline=False)

    Parameters
    ----------
    omp_nthreads : int
        Maximum number of threads an individual process may use
    fmap_bspline : bool
        Whether the fieldmap should be smoothed and extrapolated to off-brain regions
        using B-Spline basis.
    median_kernel_size : int
        Size of the kernel when smoothing is done with a median filter.
    name : str
        Name of workflow (default: ``fmap_postproc_wf``)

    Inputs
    ------
    fmap_mask : pathlike
        A brain binary mask corresponding to this fieldmap.
    fmap_ref : pathlike
        A preprocessed magnitude/reference image for the fieldmap.
    fmap : pathlike
        A B0-field nonuniformity map (aka fieldmap) estimated elsewhere.

    Outputs
    -------
    out_fmap : pathlike
        Postprocessed fieldmap.

    """
    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['fmap_mask', 'fmap_ref', 'fmap', 'metadata']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_fmap', 'metadata']),
                         name='outputnode')
    if fmap_bspline:
        from ..interfaces.fmap import FieldEnhance
        # despike_threshold=1.0, mask_erode=1),
        fmapenh = pe.Node(
            FieldEnhance(unwrap=False, despike=False),
            name='fmapenh', mem_gb=4, n_procs=omp_nthreads)

        workflow.connect([
            (inputnode, fmapenh, [('fmap_mask', 'in_mask'),
                                  ('fmap_ref', 'in_magnitude'),
                                  ('fmap_hz', 'in_file')]),
            (fmapenh, outputnode, [('out_file', 'out_fmap')]),
        ])

    else:
        recenter = pe.Node(niu.Function(function=_recenter),
                           name='recenter', run_without_submitting=True)
        denoise = pe.Node(fsl.SpatialFilter(
                          operation='median', kernel_shape='sphere',
                          kernel_size=median_kernel_size), name='denoise')
        demean = pe.Node(niu.Function(function=_demean), name='demean')
        cleanup_wf = cleanup_edge_pipeline(name="cleanup_wf")

        workflow.connect([
            (inputnode, cleanup_wf, [('fmap_mask', 'inputnode.in_mask')]),
            (inputnode, recenter, [(('fmap', _pop), 'in_file')]),
            (recenter, denoise, [('out', 'in_file')]),
            (denoise, demean, [('out_file', 'in_file')]),
            (demean, cleanup_wf, [('out', 'inputnode.in_file')]),
            (cleanup_wf, outputnode, [('outputnode.out_file', 'out_fmap')]),
            (inputnode, outputnode, [(('metadata', _pop), 'metadata')]),
        ])

    return workflow


def _recenter(in_file):
    """Recenter the phase-map distribution to the -pi..pi range."""
    from os import getcwd
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix

    nii = nb.load(in_file)
    data = nii.get_fdata(dtype='float32')
    msk = data != 0
    msk[data == 0] = False
    data[msk] -= np.median(data[msk])

    out_file = fname_presuffix(in_file, suffix='_recentered',
                               newpath=getcwd())
    nb.Nifti1Image(data, nii.affine, nii.header).to_filename(out_file)
    return out_file


def _demean(in_file, in_mask=None, usemode=True):
    """
    Subtract the median (since it is robuster than the mean) from a map.

    Parameters
    ----------
    usemode : bool
        Use the mode instead of the median (should be even more robust
        against outliers).

    """
    from os import getcwd
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix

    nii = nb.load(in_file)
    data = nii.get_fdata(dtype='float32')

    msk = np.ones_like(data, dtype=bool)
    if in_mask is not None:
        msk[nb.load(in_mask).get_fdata(dtype='float32') < 1e-4] = False

    if usemode:
        from scipy.stats import mode
        data[msk] -= mode(data[msk], axis=None)[0][0]
    else:
        data[msk] -= np.median(data[msk], axis=None)

    out_file = fname_presuffix(in_file, suffix='_demean',
                               newpath=getcwd())
    nb.Nifti1Image(data, nii.affine, nii.header).to_filename(out_file)
    return out_file


def _pop(inlist):
    if isinstance(inlist, (tuple, list)):
        return inlist[0]
    return inlist
