# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# Copyright 2021 The NiPreps Developers <nipreps@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# We support and encourage derived works from this project, please read
# about our expectations at
#
#     https://www.nipreps.org/community/licensing/
#
"""Workflows for generating reference images."""

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.header import ValidateImage

from aslprep import config
from aslprep.interfaces.reference import SelectHighestContrastVolumes
from aslprep.interfaces.utility import Smooth


def init_raw_aslref_wf(
    *,
    asl_file=None,
    reference_volume_type=None,
    m0scan=False,
    use_ge=False,
    name='raw_aslref_wf',
):
    """Build a workflow that generates reference ASL images for a series.

    CBF images have the best GM-WM contrast, but raw ASL files rarely have precomputed CBF volumes.
    As such, ASLPrep will select volumes in the following order: CBF, M0, deltam, control.

    The raw reference image is the target of :abbr:`HMC (head motion correction)`, and a
    contrast-enhanced reference is the subject of distortion correction, as well as
    boundary-based registration to T1w and template spaces.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.reference import init_raw_aslref_wf

            wf = init_raw_aslref_wf()

    Parameters
    ----------
    asl_file : :obj:`str`
        ASL series NIfTI file
    reference_volume_type : :obj:`str`
        Type of reference volume to use.
    m0scan : :obj:`bool`
        True if a separate M0 file is available. False if not.
    use_ge : :obj:`bool`
        If True, the M0 scan (when available) will be prioritized as the volume type for the
        reference image, as GE deltam volumes exhibit extreme background noise.
    name : :obj:`str`
        Name of workflow (default: ``asl_reference_wf``)

    Inputs
    ------
    asl_file : str
        ASL series NIfTI file

    Outputs
    -------
    asl_file : str
        Validated ASL series NIfTI file
    aslref : str
        Reference image to which ASL series is motion corrected
    algo_dummy_scans : int
        Number of non-steady-state volumes algorithmically detected at
        beginning of ``asl_file``
    """
    from niworkflows.interfaces.images import RobustAverage

    from aslprep.interfaces.ants import RegridImage
    from aslprep.interfaces.utility import GetImageType, IndexImage

    workflow = Workflow(name=name)
    workflow.__desc__ = f"""\
First, a reference volume was generated from the {reference_volume_type} volumes using a custom
methodology of *ASLPrep*, for use in head motion correction.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'asl_file',
                'aslcontext',
                'm0scan',
            ],
        ),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'asl_file',
                'aslref',
                'validation_report',
            ],
        ),
        name='outputnode',
    )

    # Simplify manually setting input image
    if asl_file is not None:
        inputnode.inputs.asl_file = asl_file

    val_asl = pe.Node(
        ValidateImage(),
        name='val_asl',
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([
        (inputnode, val_asl, [('asl_file', 'in_file')]),
        (val_asl, outputnode, [
            ('out_file', 'asl_file'),
            ('out_report', 'validation_report'),
        ]),
    ])  # fmt:skip

    if reference_volume_type == 'separate_m0scan':
        val_m0scan = pe.Node(
            ValidateImage(),
            name='val_m0scan',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(inputnode, val_m0scan, [('m0scan', 'in_file')])])

        # Grab first volume of ASL file
        extract_3d_asl = pe.Node(
            IndexImage(index=0),
            name='extract_3d_asl',
        )
        workflow.connect([(val_asl, extract_3d_asl, [('out_file', 'in_file')])])

        # In some cases, the separate M0 scan may have a different resolution than the ASL scan.
        # We resample the M0 scan to match the ASL scan at this point.
        get_image_type = pe.Node(GetImageType(), name='get_image_type')
        workflow.connect([(val_m0scan, get_image_type, [('out_file', 'image')])])

        # An M0 scan already on the ASL grid is passed through, as Gaussian interpolation would
        # only smooth it.
        resample_m0scan_to_asl = pe.Node(
            RegridImage(
                interpolation='Gaussian',
                transforms=['identity'],
                args='--verbose',
            ),
            name='resample_m0scan_to_asl',
        )
        workflow.connect([
            (extract_3d_asl, resample_m0scan_to_asl, [('out_file', 'reference_image')]),
            (val_m0scan, resample_m0scan_to_asl, [('out_file', 'input_image')]),
            (get_image_type, resample_m0scan_to_asl, [('image_type', 'input_image_type')]),
        ])  # fmt:skip

    select_highest_contrast_volumes = pe.Node(
        SelectHighestContrastVolumes(prioritize_m0=use_ge),
        name='select_highest_contrast_volumes',
        mem_gb=1,
    )
    workflow.connect([
        (inputnode, select_highest_contrast_volumes, [('aslcontext', 'aslcontext')]),
        (val_asl, select_highest_contrast_volumes, [('out_file', 'asl_file')]),
    ])  # fmt:skip

    if m0scan and reference_volume_type == 'separate_m0scan':
        workflow.connect([
            (resample_m0scan_to_asl, select_highest_contrast_volumes, [
                ('output_image', 'm0scan'),
            ]),
        ])  # fmt:skip

    gen_avg = pe.Node(RobustAverage(), name='gen_avg', mem_gb=1)
    workflow.connect([
        (select_highest_contrast_volumes, gen_avg, [('selected_volumes_file', 'in_file')]),
    ])  # fmt:skip

    if use_ge and (config.workflow.smooth_kernel > 0):
        workflow.__desc__ += (
            'The reference image was then smoothed with a Gaussian kernel '
            f'(FWHM = {config.workflow.smooth_kernel} mm).'
        )
        smooth_reference = pe.Node(
            Smooth(fwhm=config.workflow.smooth_kernel),
            name='smooth_reference',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([
            (gen_avg, smooth_reference, [('out_file', 'in_file')]),
            (smooth_reference, outputnode, [('out_file', 'aslref')]),
        ])  # fmt:skip
    else:
        workflow.connect([(gen_avg, outputnode, [('out_file', 'aslref')])])

    return workflow


def init_synthstrip_aslref_wf(
    disable_n4=False,
    name='synthstrip_aslref_wf',
):
    """Enhance and run brain extraction on an ASL image.

    This workflow uses SynthStrip to extract the brain from an ASLRef image.
    It is used when ``--skull-strip-asl synthstrip`` (the default) is selected.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep.workflows.asl.reference import init_synthstrip_aslref_wf

            with mock_config():
                wf = init_synthstrip_aslref_wf()

    Parameters
    ----------
    disable_n4 : :obj:`bool`
        If True, N4 bias field correction will be disabled.
    name : str
        Name of workflow (default: ``synthstrip_aslref_wf``)

    Inputs
    ------
    in_file : str
        ASLRef image

    Outputs
    -------
    bias_corrected_file : str
        the ``in_file`` after N4 bias field correction.
        If ``disable_n4`` is True, this will be the same as ``in_file``.
    mask_file : str
        SynthStrip brain mask of ``in_file``
    """
    from nipype.interfaces import utility as niu
    from nipype.pipeline import engine as pe
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.fixes import FixN4BiasFieldCorrection as N4BiasFieldCorrection

    from aslprep.workflows.segmentation import init_synthstrip_wf

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
A brain mask was computed from the reference image with SynthStrip [@synthstrip].
"""
    if not disable_n4:
        workflow.__desc__ += """\
The reference image was corrected for intensity non-uniformity with `N4BiasFieldCorrection`
[@n4].
"""

    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file']), name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['mask_file', 'bias_corrected_file']),
        name='outputnode',
    )

    n4_buffer = pe.Node(niu.IdentityInterface(fields=['bias_corrected_file']), name='n4_buffer')

    if not disable_n4:
        # Run N4 normally, force num_threads=1 for stability (images are small, no need for >1)
        n4_correct = pe.Node(
            N4BiasFieldCorrection(dimension=3, copy_header=True, bspline_fitting_distance=200),
            shrink_factor=2,
            name='n4_correct',
            n_procs=1,
        )
        n4_correct.inputs.rescale_intensities = True

        workflow.connect([
            (inputnode, n4_correct, [('in_file', 'input_image')]),
            (n4_correct, n4_buffer, [('output_image', 'bias_corrected_file')]),
        ])  # fmt:skip
    else:
        workflow.connect([(inputnode, n4_buffer, [('in_file', 'bias_corrected_file')])])

    synthstrip_wf = init_synthstrip_wf(name='synthstrip_wf')

    workflow.connect([
        (inputnode, synthstrip_wf, [('in_file', 'inputnode.original_image')]),
        (n4_buffer, outputnode, [('bias_corrected_file', 'bias_corrected_file')]),
        (synthstrip_wf, outputnode, [('outputnode.brain_mask', 'mask_file')]),
    ])  # fmt: skip

    return workflow


def init_skullstrip_aslref_wf(name='skullstrip_aslref_wf'):
    """Compute a brain mask for an ASL reference image that needs no further enhancement.

    This is used for references that are already final, such as a distortion-corrected
    or precomputed reference.
    The method follows ``--skull-strip-asl``:
    SynthStrip (:func:`init_synthstrip_aslref_wf` without N4) for ``synthstrip``
    and ``anatomical`` (which checks the anatomical mask against it,
    see :func:`init_anat_aslmask_wf`),
    and niworkflows' :func:`~niworkflows.func.util.init_skullstrip_bold_wf`
    (BET and ``3dAutomask``) for ``legacy``, as *ASLPrep* did before SynthStrip.

    Parameters
    ----------
    name : str
        Name of workflow (default: ``skullstrip_aslref_wf``)

    Inputs
    ------
    in_file : str
        ASLRef image

    Outputs
    -------
    mask_file : str
        brain mask of ``in_file``
    """
    if config.workflow.skull_strip_asl == 'legacy':
        from niworkflows.func.util import init_skullstrip_bold_wf

        return init_skullstrip_bold_wf(name=name)

    return init_synthstrip_aslref_wf(disable_n4=True, name=name)


def init_anat_aslmask_wf(name='anat_aslmask_wf'):
    """Take the ASL brain mask from the anatomical brain mask, checked against SynthStrip's.

    This is used when ``--skull-strip-asl anatomical`` is selected.
    The anatomical brain mask is resampled to the ASL reference through the inverse of the
    coregistration transform.
    The coregistration does not use the ASL brain mask, so there is no cycle.

    The anatomical mask is only as good as the coregistration, so it is checked against the
    SynthStrip mask of the ASL reference (:class:`~aslprep.interfaces.reference.GuardAnatomicalMask`):
    when the two disagree (Dice below 0.9) and SynthStrip's mask has at least 95% of the
    anatomical mask's volume, the coregistration has most likely failed, and the SynthStrip mask
    is used instead.
    The Dice coefficient, the volume ratio and which mask was used are added to the
    quality-control measures.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep.workflows.asl.reference import init_anat_aslmask_wf

            with mock_config():
                wf = init_anat_aslmask_wf()

    Parameters
    ----------
    name : str
        Name of workflow (default: ``anat_aslmask_wf``)

    Inputs
    ------
    aslref
        ASL reference image (the coregistration reference)
    asl_mask
        SynthStrip brain mask of ``aslref``
    t1w_mask
        Anatomical brain mask
    aslref2anat_xfm
        Transform from ``aslref`` to the anatomical space

    Outputs
    -------
    mask_file
        The brain mask to use for ``aslref``
    mask_metadata
        Sidecar metadata of ``mask_file``, recording which mask was used
    qc
        Quality-control measures of the check
    """
    from aslprep.interfaces.ants import ApplyTransforms
    from aslprep.interfaces.reference import GuardAnatomicalMask

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
The brain mask of the ASL reference was then taken from the anatomical brain mask, resampled
into the space of the ASL reference with the inverse of the ASL-to-anatomical transform.
Because this mask depends on the coregistration, it was checked against the SynthStrip mask:
where the two masks disagreed (Dice coefficient below 0.9) and the SynthStrip mask was not the
smaller one (at least 95% of the anatomical mask's volume), which indicates a failed
coregistration, the SynthStrip mask was used instead.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['aslref', 'asl_mask', 't1w_mask', 'aslref2anat_xfm']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['mask_file', 'mask_metadata', 'qc']),
        name='outputnode',
    )

    anat_mask_to_aslref = pe.Node(
        ApplyTransforms(interpolation='GenericLabel', invert_transform_flags=[True], args='-v'),
        name='anat_mask_to_aslref',
        mem_gb=0.1,
    )
    guard = pe.Node(GuardAnatomicalMask(), name='guard')

    workflow.connect([
        (inputnode, anat_mask_to_aslref, [
            ('t1w_mask', 'input_image'),
            ('aslref', 'reference_image'),
            ('aslref2anat_xfm', 'transforms'),
        ]),
        (anat_mask_to_aslref, guard, [('output_image', 'anat_mask')]),
        (inputnode, guard, [('asl_mask', 'asl_mask')]),
        (guard, outputnode, [
            ('out_mask', 'mask_file'),
            ('metadata', 'mask_metadata'),
            ('qc', 'qc'),
        ]),
    ])  # fmt:skip

    return workflow


def init_enhance_and_skullstrip_asl_wf(
    disable_n4=False,
    omp_nthreads=1,
    name='enhance_and_skullstrip_asl_wf',
):
    """Enhance and run brain extraction on an ASL image with BET and 3dAutomask.

    This is a modified version of the ``enhance_and_skullstrip_bold_wf`` workflow from
    niworkflows, with an added ``disable_n4`` parameter to disable N4 bias field correction,
    which tends to be unstable on data with prescan normalization enabled.
    It was *ASLPrep*'s brain extraction method before SynthStrip,
    and is used when ``--skull-strip-asl legacy`` is selected.

    The steps are:

    1. If N4 is enabled, register the reference image to fMRIPrep's MNI152NLin2009cAsym
       BOLD reference template and map the template's brain probability map back,
       to weight N4 bias field correction.
    2. Compute a liberal brain mask with FSL's BET (``-f 0.2``) on the bias-corrected image.
    3. Dilate that mask, apply it, and normalize the intensities of the masked image with
       AFNI's ``3dUnifize``.
    4. Compute a second mask with AFNI's ``3dAutomask`` on the normalized image.
    5. Take the intersection of the BET and ``3dAutomask`` masks.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep.workflows.asl.reference import init_enhance_and_skullstrip_asl_wf

            with mock_config():
                wf = init_enhance_and_skullstrip_asl_wf()

    Parameters
    ----------
    disable_n4 : :obj:`bool`
        If True, N4 bias field correction will be disabled,
        along with the template registration that only serves to weight it.
    omp_nthreads : :obj:`int`
        Number of threads available to the template registration.
    name : str
        Name of workflow (default: ``enhance_and_skullstrip_asl_wf``)

    Inputs
    ------
    in_file : str
        ASLRef image

    Outputs
    -------
    bias_corrected_file : str
        the ``in_file`` after N4 bias field correction.
        If ``disable_n4`` is True, this will be the same as ``in_file``.
    mask_file : str
        Intersection of the BET and ``3dAutomask`` brain masks of ``in_file``
    """
    from nipype.interfaces import afni, fsl
    from niworkflows.interfaces.header import CopyXForm
    from niworkflows.interfaces.nibabel import ApplyMask, BinaryDilation

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
A brain mask was computed from the reference image as the intersection of two masks:
a liberal mask from FSL's `bet` [@bet], and a mask from AFNI's `3dAutomask` [@afni],
computed on the BET-masked image after its intensities were normalized with AFNI's `3dUnifize`.
"""
    if not disable_n4:
        workflow.__desc__ += """\
Before brain extraction, the reference image was corrected for intensity non-uniformity with
`N4BiasFieldCorrection` [@n4], weighted by a brain mask mapped from the MNI152NLin2009cAsym
template [@mni152nlin2009casym] after registering the reference image to the template with
`antsRegistration` [@ants].
"""

    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file']), name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['mask_file', 'bias_corrected_file']),
        name='outputnode',
    )

    n4_buffer = pe.Node(niu.IdentityInterface(fields=['bias_corrected_file']), name='n4_buffer')

    if not disable_n4:
        from nipype.interfaces.ants.utils import AI
        from niworkflows import data
        from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
        from niworkflows.interfaces.fixes import FixHeaderRegistration as Registration
        from niworkflows.interfaces.fixes import (
            FixN4BiasFieldCorrection as N4BiasFieldCorrection,
        )
        from niworkflows.interfaces.header import CopyHeader
        from packaging.version import Version
        from packaging.version import parse as parseversion
        from templateflow.api import get as get_template

        asl_template = get_template(
            'MNI152NLin2009cAsym',
            resolution=2,
            desc='fMRIPrep',
            suffix='boldref',
        )
        brain_mask = get_template('MNI152NLin2009cAsym', resolution=2, desc='brain', suffix='mask')

        # Initialize transforms with antsAI
        init_aff = pe.Node(
            AI(
                fixed_image=str(asl_template),
                fixed_image_mask=str(brain_mask),
                metric=('Mattes', 32, 'Regular', 0.2),
                transform=('Affine', 0.1),
                search_factor=(20, 0.12),
                principal_axes=False,
                convergence=(10, 1e-6, 10),
                verbose=True,
            ),
            name='init_aff',
            n_procs=omp_nthreads,
        )

        # Registration().version may be None
        if parseversion(Registration().version or '0.0.0') > Version('2.2.0'):
            init_aff.inputs.search_grid = (40, (0, 40, 40))

        # Set up spatial normalization
        norm = pe.Node(
            Registration(from_file=data.load('epi_atlasbased_brainmask.json')),
            name='norm',
            n_procs=omp_nthreads,
        )
        norm.inputs.fixed_image = str(asl_template)
        map_brainmask = pe.Node(
            ApplyTransforms(
                interpolation='Linear',
                # Use the higher resolution and probseg for numerical stability in rounding
                input_image=str(
                    get_template(
                        'MNI152NLin2009cAsym',
                        resolution=1,
                        label='brain',
                        suffix='probseg',
                    )
                ),
            ),
            name='map_brainmask',
        )
        # Ensure mask's header matches reference's
        fix_header = pe.Node(CopyHeader(), name='fix_header', run_without_submitting=True)

        # Run N4 normally, force num_threads=1 for stability (images are small, no need for >1)
        n4_correct = pe.Node(
            N4BiasFieldCorrection(dimension=3, copy_header=True, bspline_fitting_distance=200),
            shrink_factor=2,
            name='n4_correct',
            n_procs=1,
        )
        n4_correct.inputs.rescale_intensities = True

        workflow.connect([
            (inputnode, fix_header, [('in_file', 'hdr_file')]),
            (inputnode, init_aff, [('in_file', 'moving_image')]),
            (inputnode, map_brainmask, [('in_file', 'reference_image')]),
            (inputnode, norm, [('in_file', 'moving_image')]),
            (init_aff, norm, [('output_transform', 'initial_moving_transform')]),
            (norm, map_brainmask, [
                ('reverse_invert_flags', 'invert_transform_flags'),
                ('reverse_transforms', 'transforms'),
            ]),
            (map_brainmask, fix_header, [('output_image', 'in_file')]),
            (inputnode, n4_correct, [('in_file', 'input_image')]),
            (fix_header, n4_correct, [('out_file', 'weight_image')]),
            (n4_correct, n4_buffer, [('output_image', 'bias_corrected_file')]),
        ])  # fmt:skip
    else:
        workflow.connect([(inputnode, n4_buffer, [('in_file', 'bias_corrected_file')])])

    # Create a generous BET mask out of the bias-corrected EPI
    skullstrip_first_pass = pe.Node(fsl.BET(frac=0.2, mask=True), name='skullstrip_first_pass')
    first_dilate = pe.Node(BinaryDilation(radius=6), name='first_dilate')
    first_mask = pe.Node(ApplyMask(), name='first_mask')

    # Use AFNI's unifize for T2 contrast & fix header
    unifize = pe.Node(
        afni.Unifize(
            t2=True,
            outputtype='NIFTI_GZ',
            # Default -clfrac is 0.1, 0.4 was too conservative
            # -rbt because I'm a Jedi AFNI Master (see 3dUnifize's documentation)
            args='-clfrac 0.2 -rbt 18.3 65.0 90.0',
            out_file='uni.nii.gz',
        ),
        name='unifize',
    )
    fixhdr_unifize = pe.Node(CopyXForm(), name='fixhdr_unifize', mem_gb=0.1)

    # Run AFNI's 3dAutomask to extract a refined brain mask
    skullstrip_second_pass = pe.Node(
        afni.Automask(dilate=1, outputtype='NIFTI_GZ'),
        name='skullstrip_second_pass',
    )
    fixhdr_skullstrip2 = pe.Node(CopyXForm(), name='fixhdr_skullstrip2', mem_gb=0.1)

    # Take intersection of both masks
    combine_masks = pe.Node(fsl.BinaryMaths(operation='mul'), name='combine_masks')

    workflow.connect([
        (inputnode, fixhdr_unifize, [('in_file', 'hdr_file')]),
        (inputnode, fixhdr_skullstrip2, [('in_file', 'hdr_file')]),
        (n4_buffer, skullstrip_first_pass, [('bias_corrected_file', 'in_file')]),
        (skullstrip_first_pass, first_dilate, [('mask_file', 'in_file')]),
        (first_dilate, first_mask, [('out_file', 'in_mask')]),
        (skullstrip_first_pass, first_mask, [('out_file', 'in_file')]),
        (first_mask, unifize, [('out_file', 'in_file')]),
        (unifize, fixhdr_unifize, [('out_file', 'in_file')]),
        (fixhdr_unifize, skullstrip_second_pass, [('out_file', 'in_file')]),
        (skullstrip_first_pass, combine_masks, [('mask_file', 'in_file')]),
        (skullstrip_second_pass, fixhdr_skullstrip2, [('out_file', 'in_file')]),
        (fixhdr_skullstrip2, combine_masks, [('out_file', 'operand_file')]),
        (combine_masks, outputnode, [('out_file', 'mask_file')]),
        (n4_buffer, outputnode, [('bias_corrected_file', 'bias_corrected_file')]),
    ])  # fmt:skip

    return workflow
