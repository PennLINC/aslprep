# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Fieldmap-based estimation of susceptibility distortions.

.. _sdc_direct_b0 :

Direct B0 mapping sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~
When the fieldmap is directly measured with a prescribed sequence (such as
:abbr:`SE (spiral echo)`), we only need to calculate the corresponding
displacements field that accounts for the distortions.
This procedure is described with more detail `here
<https://cni.stanford.edu/wiki/GE_Processing#Fieldmaps>`__.

This corresponds to `this section of the BIDS specification
<https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#a-real-fieldmap-image>`__.

"""
import pkg_resources as pkgr

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.images import IntraModalMerge
from niworkflows.interfaces.registration import ANTSApplyTransformsRPT, ANTSRegistrationRPT

from ..interfaces.fmap import get_ees as _get_ees, FieldToRadS, FUGUEvsm2ANTSwarp
from .gre import init_fmap_postproc_wf, init_magnitude_wf


def init_fmap_wf(omp_nthreads, fmap_bspline, name='fmap_wf'):
    """
    Estimate the fieldmap based on a field-mapping MRI acquisition.

    When we have a sequence that directly measures the fieldmap,
    we just need to mask it (using the corresponding magnitude image)
    to remove the noise in the surrounding air region, and ensure that
    units are Hz.

    Workflow Graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from sdcflows.workflows.fmap import init_fmap_wf
            wf = init_fmap_wf(omp_nthreads=6, fmap_bspline=False)

    Parameters
    ----------
    omp_nthreads : int
        Maximum number of threads an individual process may use.
    fmap_bspline : bool
        Whether the fieldmap estimate will be smoothed using BSpline basis.
    name : str
        Unique name of this workflow.

    Inputs
    ------
    magnitude : str
        Path to the corresponding magnitude image for anatomical reference.
    fieldmap : str
        Path to the fieldmap acquisition (``*_fieldmap.nii[.gz]`` of BIDS).

    Outputs
    -------
    fmap : str
        Path to the estimated fieldmap.
    fmap_ref : str
        Path to a preprocessed magnitude image reference.
    fmap_mask : str
        Path to a binary brain mask corresponding to the ``fmap`` and ``fmap_ref``
        pair.

    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
A B0-nonuniformity map (or *fieldmap*) was directly measured with an MRI scheme
designed with that purpose (typically, a spiral pulse sequence).
"""
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['magnitude', 'fieldmap']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['fmap', 'fmap_ref', 'fmap_mask']),
                         name='outputnode')

    magnitude_wf = init_magnitude_wf(omp_nthreads=omp_nthreads)
    workflow.connect([
        (inputnode, magnitude_wf, [('magnitude', 'inputnode.magnitude')]),
        (magnitude_wf, outputnode, [('outputnode.fmap_mask', 'fmap_mask'),
                                    ('outputnode.fmap_ref', 'fmap_ref')]),
    ])

    # Merge input fieldmap images
    fmapmrg = pe.Node(IntraModalMerge(zero_based_avg=False, hmc=False),
                      name='fmapmrg')
    applymsk = pe.Node(fsl.ApplyMask(), name='applymsk')
    fmap_postproc_wf = init_fmap_postproc_wf(omp_nthreads=omp_nthreads,
                                             fmap_bspline=fmap_bspline)

    workflow.connect([
        (inputnode, fmapmrg, [('fieldmap', 'in_files')]),
        (fmapmrg, applymsk, [('out_avg', 'in_file')]),
        (magnitude_wf, applymsk, [('outputnode.fmap_mask', 'mask_file')]),
        (applymsk, fmap_postproc_wf, [('out_file', 'inputnode.fmap')]),
        (magnitude_wf, fmap_postproc_wf, [
            ('outputnode.fmap_mask', 'inputnode.fmap_mask'),
            ('outputnode.fmap_ref', 'inputnode.fmap_ref')]),
        (fmap_postproc_wf, outputnode, [('outputnode.out_fmap', 'fmap')]),
    ])
    return workflow


def init_fmap2field_wf(omp_nthreads, debug, name='fmap2field_wf',
                       generate_report=True):
    """
    Convert the estimated fieldmap in Hz into a displacements field.

    This workflow takes in a fieldmap and calculates the corresponding
    displacements field (in other words, an ANTs-compatible warp file).

    Workflow Graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from sdcflows.workflows.fmap import init_fmap2field_wf
            wf = init_fmap2field_wf(omp_nthreads=8,
                                    debug=False)

    Parameters
    ----------
    omp_nthreads : int
        Maximum number of threads an individual process may use.
    debug : bool
        Run fast configurations of registrations.
    name : str
        Unique name of this workflow.

    Inputs
    ------
    in_reference
        the reference image
    in_reference_brain
        the reference image (skull-stripped)
    metadata
        metadata associated to the ``in_reference`` EPI input
    fmap
        the fieldmap in Hz
    fmap_ref
        the reference (anatomical) image corresponding to ``fmap``
    fmap_mask
        a brain mask corresponding to ``fmap``


    Outputs
    -------
    out_reference
        the ``in_reference`` after unwarping
    out_reference_brain
        the ``in_reference`` after unwarping and skullstripping
    out_warp
        the corresponding :abbr:`DFM (displacements field map)` compatible with
        ANTs
    out_jacobian
        the jacobian of the field (for drop-out alleviation)
    out_mask
        mask of the unwarped input file

    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The *fieldmap* was then co-registered to the target EPI (echo-planar imaging)
reference run and converted to a displacements field map (amenable to registration
tools such as ANTs) with FSL's `fugue` and other *SDCflows* tools.
"""
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_reference', 'in_reference_brain', 'metadata',
                'fmap_ref', 'fmap_mask', 'fmap']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_warp']), name='outputnode')

    # Register the reference of the fieldmap to the reference
    # of the target image (the one that shall be corrected)
    ants_settings = pkgr.resource_filename('sdcflows', 'data/fmap-any_registration.json')
    if debug:
        ants_settings = pkgr.resource_filename(
            'sdcflows', 'data/fmap-any_registration_testing.json')

    fmap2ref_reg = pe.Node(
        ANTSRegistrationRPT(generate_report=False, from_file=ants_settings,
                            output_inverse_warped_image=True, output_warped_image=True),
        name='fmap2ref_reg', n_procs=omp_nthreads)

    # Map the VSM into the EPI space
    fmap2ref_apply = pe.Node(ANTSApplyTransformsRPT(
        generate_report=False, dimension=3, interpolation='BSpline', float=True),
        name='fmap2ref_apply')

    fmap_mask2ref_apply = pe.Node(ANTSApplyTransformsRPT(
        generate_report=False, dimension=3, interpolation='MultiLabel',
        float=True),
        name='fmap_mask2ref_apply')

    # Fieldmap to rads and then to voxels (VSM - voxel shift map)
    torads = pe.Node(FieldToRadS(fmap_range=0.5), name='torads')

    get_ees = pe.Node(niu.Function(function=_get_ees, output_names=['ees']), name='get_ees')

    gen_vsm = pe.Node(fsl.FUGUE(save_unmasked_shift=True), name='gen_vsm')
    # Convert the VSM into a DFM (displacements field map)
    # or: FUGUE shift to ANTS warping.
    vsm2dfm = pe.Node(FUGUEvsm2ANTSwarp(), name='vsm2dfm')

    workflow.connect([
        (inputnode, fmap2ref_reg, [('fmap_ref', 'moving_image'),
                                   ('in_reference_brain', 'fixed_image')]),
        (inputnode, fmap2ref_apply, [('fmap', 'input_image'),
                                     ('in_reference', 'reference_image')]),
        (inputnode, fmap_mask2ref_apply, [('in_reference', 'reference_image'),
                                          ('fmap_mask', 'input_image')]),
        (inputnode, get_ees, [('in_reference', 'in_file'),
                              ('metadata', 'in_meta')]),
        (inputnode, gen_vsm, [(('metadata', _get_pedir_fugue), 'unwarp_direction')]),
        (inputnode, vsm2dfm, [(('metadata', _get_pedir_bids), 'pe_dir')]),
        (fmap2ref_reg, fmap2ref_apply, [
            ('composite_transform', 'transforms')]),
        (fmap2ref_reg, fmap_mask2ref_apply, [
            ('composite_transform', 'transforms')]),
        (fmap2ref_apply, torads, [('output_image', 'in_file')]),
        (fmap_mask2ref_apply, gen_vsm, [('output_image', 'mask_file')]),
        (gen_vsm, vsm2dfm, [('shift_out_file', 'in_file')]),
        (get_ees, gen_vsm, [('ees', 'dwell_time')]),
        (torads, gen_vsm, [('out_file', 'fmap_in_file')]),
        (vsm2dfm, outputnode, [('out_file', 'out_warp')]),
    ])

    if generate_report:
        from niworkflows.interfaces.bids import DerivativesDataSink
        from ..interfaces.reportlets import FieldmapReportlet

        fmap_rpt = pe.Node(FieldmapReportlet(
            reference_label='EPI Reference',
            moving_label='Magnitude', show='both'), name='fmap_rpt')
        ds_report_sdc = pe.Node(
            DerivativesDataSink(desc='fieldmap', suffix='bold', datatype='figures'),
            name='ds_report_fmap', mem_gb=0.01, run_without_submitting=True
        )

        workflow.connect([
            (inputnode, fmap_rpt, [('in_reference', 'reference')]),
            (fmap2ref_reg, fmap_rpt, [('warped_image', 'moving')]),
            (fmap_mask2ref_apply, fmap_rpt, [('output_image', 'mask')]),
            (vsm2dfm, fmap_rpt, [('fieldmap', 'fieldmap')]),
            (fmap_rpt, ds_report_sdc, [('out_report', 'in_file')]),
        ])

    return workflow


# Helper functions
# ------------------------------------------------------------
def _get_pedir_bids(in_dict):
    return in_dict['PhaseEncodingDirection']


def _get_pedir_fugue(in_dict):
    return in_dict['PhaseEncodingDirection'].replace('i', 'x').replace('j', 'y').replace('k', 'z')
