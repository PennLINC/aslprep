# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Datasets with multiple phase encoded directions.

.. _sdc_pepolar :

:abbr:`PEPOLAR (Phase Encoding POLARity)` techniques
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This corresponds to `this section of the BIDS specification
<https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#case-4-multiple-phase-encoded-directions-pepolar>`__.

"""

import pkg_resources as pkgr

from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces import CopyHeader
from niworkflows.interfaces.freesurfer import StructuralReference
from niworkflows.interfaces.registration import ANTSApplyTransformsRPT
from niworkflows.func.util import init_enhance_and_skullstrip_bold_wf

from nipype.pipeline import engine as pe
from nipype.interfaces import afni, ants, utility as niu


def init_pepolar_unwarp_wf(omp_nthreads=1, matched_pe=False,
                           name="pepolar_unwarp_wf"):
    """
    Create the PE-Polar field estimation workflow.

    This workflow takes in a set of EPI files with opposite phase encoding
    direction than the target file and calculates a displacements field
    (in other words, an ANTs-compatible warp file).

    This procedure works if there is only one ``_epi`` file is present
    (as long as it has the opposite phase encoding direction to the target
    file). The target file will be used to estimate the field distortion.
    However, if there is another ``_epi`` file present with a matching
    phase encoding direction to the target it will be used instead.

    Currently, different phase encoding directions in the target file and the
    ``_epi`` file(s) (for example, ``i`` and ``j``) is not supported.

    The warp field correcting for the distortions is estimated using AFNI's
    ``3dQwarp``, with displacement estimation limited to the target file phase
    encoding direction.

    It also calculates a new mask for the input dataset that takes into
    account the distortions.

    Workflow Graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from sdcflows.workflows.pepolar import init_pepolar_unwarp_wf
            wf = init_pepolar_unwarp_wf()

    Parameters
    ----------
    matched_pe : bool
        Whether the input ``fmaps_epi`` will contain images with matched
        PE blips or not. Please use :func:`sdcflows.workflows.pepolar.check_pes`
        to determine whether they exist or not.
    name : str
        Name for this workflow
    omp_nthreads : int
        Parallelize internal tasks across the number of CPUs given by this option.

    Inputs
    ------
    fmaps_epi : list of tuple(pathlike, str)
        The list of EPI images that will be used in PE-Polar correction, and
        their corresponding ``PhaseEncodingDirection`` metadata.
        The workflow will use the ``epi_pe_dir`` input to separate out those
        EPI acquisitions with opposed PE blips and those with matched PE blips
        (the latter could be none, and ``in_reference_brain`` would then be
        used). The workflow raises a ``ValueError`` when no images with
        opposed PE blips are found.
    epi_pe_dir : str
        The baseline PE direction.
    in_reference : pathlike
        The baseline reference image (must correspond to ``epi_pe_dir``).
    in_reference_brain : pathlike
        The reference image above, but skullstripped.
    in_mask : pathlike
        Not used, present only for consistency across fieldmap estimation
        workflows.

    Outputs
    -------
    out_reference : pathlike
        The ``in_reference`` after unwarping
    out_reference_brain : pathlike
        The ``in_reference`` after unwarping and skullstripping
    out_warp : pathlike
        The corresponding :abbr:`DFM (displacements field map)` compatible with
        ANTs.
    out_mask : pathlike
        Mask of the unwarped input file

    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
A B0-nonuniformity map (or *fieldmap*) was estimated based on two (or more)
echo-planar imaging (EPI) references with opposing phase-encoding
directions, with `3dQwarp` @afni (AFNI {afni_ver}).
""".format(afni_ver=''.join(['%02d' % v for v in afni.Info().version() or []]))

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['fmaps_epi', 'in_reference', 'in_reference_brain',
                'in_mask', 'epi_pe_dir']), name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_reference', 'out_reference_brain', 'out_warp', 'out_mask']),
        name='outputnode')

    prepare_epi_wf = init_prepare_epi_wf(omp_nthreads=omp_nthreads,
                                         matched_pe=matched_pe,
                                         name="prepare_epi_wf")

    qwarp = pe.Node(afni.QwarpPlusMinus(
        pblur=[0.05, 0.05], blur=[-1, -1], noweight=True, minpatch=9, nopadWARP=True,
        environ={'OMP_NUM_THREADS': '%d' % omp_nthreads}),
        name='qwarp', n_procs=omp_nthreads)

    to_ants = pe.Node(niu.Function(function=_fix_hdr), name='to_ants',
                      mem_gb=0.01)

    cphdr_warp = pe.Node(CopyHeader(), name='cphdr_warp', mem_gb=0.01)

    unwarp_reference = pe.Node(ANTSApplyTransformsRPT(dimension=3,
                                                      generate_report=False,
                                                      float=True,
                                                      interpolation='LanczosWindowedSinc'),
                               name='unwarp_reference')

    enhance_and_skullstrip_bold_wf = init_enhance_and_skullstrip_bold_wf(
        omp_nthreads=omp_nthreads)

    workflow.connect([
        (inputnode, qwarp, [(('epi_pe_dir', _qwarp_args), 'args')]),
        (inputnode, cphdr_warp, [('in_reference', 'hdr_file')]),
        (inputnode, prepare_epi_wf, [
            ('fmaps_epi', 'inputnode.maps_pe'),
            ('epi_pe_dir', 'inputnode.epi_pe'),
            ('in_reference_brain', 'inputnode.ref_brain')]),
        (prepare_epi_wf, qwarp, [('outputnode.opposed_pe', 'base_file'),
                                 ('outputnode.matched_pe', 'in_file')]),
        (qwarp, cphdr_warp, [('source_warp', 'in_file')]),
        (cphdr_warp, to_ants, [('out_file', 'in_file')]),
        (to_ants, unwarp_reference, [('out', 'transforms')]),
        (inputnode, unwarp_reference, [('in_reference', 'reference_image'),
                                       ('in_reference', 'input_image')]),
        (unwarp_reference, enhance_and_skullstrip_bold_wf, [
            ('output_image', 'inputnode.in_file')]),
        (unwarp_reference, outputnode, [('output_image', 'out_reference')]),
        (enhance_and_skullstrip_bold_wf, outputnode, [
            ('outputnode.mask_file', 'out_mask'),
            ('outputnode.skull_stripped_file', 'out_reference_brain')]),
        (to_ants, outputnode, [('out', 'out_warp')]),
    ])

    return workflow


def init_prepare_epi_wf(omp_nthreads, matched_pe=False,
                        name="prepare_epi_wf"):
    """
    Prepare opposed-PE EPI images for PE-POLAR SDC.

    This workflow takes in a set of EPI files and returns two 3D volumes with
    matching and opposed PE directions, ready to be used in field distortion
    estimation.

    The procedure involves: estimating a robust template using FreeSurfer's
    ``mri_robust_template``, bias field correction using ANTs ``N4BiasFieldCorrection``
    and AFNI ``3dUnifize``, skullstripping using FSL BET and AFNI ``3dAutomask``,
    and rigid coregistration to the reference using ANTs.

    Workflow Graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from sdcflows.workflows.pepolar import init_prepare_epi_wf
            wf = init_prepare_epi_wf(omp_nthreads=8)

    Parameters
    ----------
    matched_pe : bool
        Whether the input ``fmaps_epi`` will contain images with matched
        PE blips or not. Please use :func:`sdcflows.workflows.pepolar.check_pes`
        to determine whether they exist or not.
    name : str
        Name for this workflow
    omp_nthreads : int
        Parallelize internal tasks across the number of CPUs given by this option.

    Inputs
    ------
    epi_pe : str
        Phase-encoding direction of the EPI image to be corrected.
    maps_pe : list of tuple(pathlike, str)
        list of 3D or 4D NIfTI images
    ref_brain
        coregistration reference (skullstripped and bias field corrected)

    Outputs
    -------
    opposed_pe : pathlike
        single 3D NIfTI file
    matched_pe : pathlike
        single 3D NIfTI file

    """
    inputnode = pe.Node(niu.IdentityInterface(fields=['epi_pe', 'maps_pe', 'ref_brain']),
                        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(fields=['opposed_pe', 'matched_pe']),
                         name='outputnode')

    ants_settings = pkgr.resource_filename('sdcflows',
                                           'data/translation_rigid.json')

    split = pe.Node(niu.Function(function=_split_epi_lists), name='split')

    merge_op = pe.Node(
        StructuralReference(auto_detect_sensitivity=True,
                            initial_timepoint=1,
                            fixed_timepoint=True,  # Align to first image
                            intensity_scaling=True,
                            # 7-DOF (rigid + intensity)
                            no_iteration=True,
                            subsample_threshold=200,
                            out_file='template.nii.gz'),
        name='merge_op')

    ref_op_wf = init_enhance_and_skullstrip_bold_wf(
        omp_nthreads=omp_nthreads, name='ref_op_wf')

    op2ref_reg = pe.Node(ants.Registration(
        from_file=ants_settings, output_warped_image=True),
        name='op2ref_reg', n_procs=omp_nthreads)

    workflow = Workflow(name=name)
    workflow.connect([
        (inputnode, split, [('maps_pe', 'in_files'),
                            ('epi_pe', 'pe_dir')]),
        (split, merge_op, [(('out', _front), 'in_files')]),
        (merge_op, ref_op_wf, [('out_file', 'inputnode.in_file')]),
        (ref_op_wf, op2ref_reg, [
            ('outputnode.skull_stripped_file', 'moving_image')]),
        (inputnode, op2ref_reg, [('ref_brain', 'fixed_image')]),
        (op2ref_reg, outputnode, [('warped_image', 'opposed_pe')]),
    ])

    if not matched_pe:
        workflow.connect([
            (inputnode, outputnode, [('ref_brain', 'matched_pe')]),
        ])
        return workflow

    merge_ma = pe.Node(
        StructuralReference(auto_detect_sensitivity=True,
                            initial_timepoint=1,
                            fixed_timepoint=True,  # Align to first image
                            intensity_scaling=True,
                            # 7-DOF (rigid + intensity)
                            no_iteration=True,
                            subsample_threshold=200,
                            out_file='template.nii.gz'),
        name='merge_ma')

    ref_ma_wf = init_enhance_and_skullstrip_bold_wf(
        omp_nthreads=omp_nthreads, name='ref_ma_wf')

    ma2ref_reg = pe.Node(ants.Registration(
        from_file=ants_settings, output_warped_image=True),
        name='ma2ref_reg', n_procs=omp_nthreads)

    workflow.connect([
        (split, merge_ma, [(('out', _last), 'in_files')]),
        (merge_ma, ref_ma_wf, [('out_file', 'inputnode.in_file')]),
        (ref_ma_wf, ma2ref_reg, [
            ('outputnode.skull_stripped_file', 'moving_image')]),
        (inputnode, ma2ref_reg, [('ref_brain', 'fixed_image')]),
        (ma2ref_reg, outputnode, [('warped_image', 'matched_pe')]),
    ])
    return workflow


def _front(inlist):
    if isinstance(inlist, (list, tuple)):
        return inlist[0]
    return inlist


def _last(inlist):
    if isinstance(inlist, (list, tuple)):
        return inlist[-1]
    return inlist


def _fix_hdr(in_file, newpath=None):
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix

    nii = nb.load(in_file)
    hdr = nii.header.copy()
    hdr.set_data_dtype('<f4')
    hdr.set_intent('vector', (), '')
    out_file = fname_presuffix(in_file, "_warpfield", newpath=newpath)
    nb.Nifti1Image(nii.get_fdata(dtype='float32'), nii.affine, hdr).to_filename(
        out_file)
    return out_file


def check_pes(epi_fmaps, pe_dir):
    """Check whether there are images with matched PE."""
    opposed_pe = False
    matched_pe = False

    for _, fmap_pe in epi_fmaps:
        if fmap_pe == pe_dir:
            matched_pe = True
        elif fmap_pe[0] == pe_dir[0]:
            opposed_pe = True

    if not opposed_pe:
        raise ValueError("""\
None of the discovered fieldmaps has the right phase encoding direction. \
This is possibly a problem with metadata. If not, rerun with \
``--ignore fieldmaps`` to skip the distortion correction step.""")

    return matched_pe


def _qwarp_args(pe_dir):
    return {'i': '-noYdis -noZdis', 'j': '-noXdis -noZdis',
            'k': '-noXdis -noYdis'}[pe_dir[0]]


def _split_epi_lists(in_files, pe_dir, max_trs=50):
    """
    Split input EPIs and generate an output list of PEs.

    Inputs
    ------
    in_files : list of ``BIDSFile``s
        The EPI images that will be pooled into field estimation.
    pe_dir : str
        The phase-encoding direction of the IntendedFor EPI scan.
    max_trs : int
        Index of frame after which all volumes will be discarded
        from the input EPI images.

    """
    from os import path as op
    import nibabel as nb

    matched_pe = []
    opposed_pe = []

    for i, (epi_path, epi_pe) in enumerate(in_files):
        if epi_pe[0] == pe_dir[0]:
            img = nb.load(epi_path)
            if len(img.shape) == 3:
                splitnii = [img]
            else:
                splitnii = nb.four_to_three(img.slicer[:, :, :, :max_trs])

            for j, nii in enumerate(splitnii):
                out_name = op.abspath(
                    'dir-%s_tstep-%03d_pe-%03d.nii.gz' % (epi_pe, i, j))
                nii.to_filename(out_name)

                if epi_pe == pe_dir:
                    matched_pe.append(out_name)
                else:
                    opposed_pe.append(out_name)

    if matched_pe:
        return [opposed_pe, matched_pe]

    return [opposed_pe]
