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
from aslprep.interfaces.reference import BrainChop, SelectHighestContrastVolumes
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

    from aslprep.interfaces.ants import ApplyTransforms
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

        resample_m0scan_to_asl = pe.Node(
            ApplyTransforms(
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
    if m0scan:
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


def init_enhance_and_skullstrip_asl_wf(
    name='enhance_and_skullstrip_bold_wf',
    bias_correct=False,
):
    """
    Enhance and run brain extraction on a BOLD EPI image.

    NOTE: This is a modified version of the bold workflow from niworkflows.
    I have added an option to skip the N4 bias field correction step,
    and have swapped out BET for SynthStrip.

    This workflow takes in a :abbr:`BOLD (blood-oxygen level-dependant)`
    :abbr:`fMRI (functional MRI)` average/summary (e.g., a reference image
    averaging non-steady-state timepoints), and sharpens the histogram
    with the application of the N4 algorithm for removing the
    :abbr:`INU (intensity non-uniformity)` bias field and calculates a signal
    mask.

    Steps of this workflow are:

      1. Calculate a tentative mask by registering (9-parameters) to *fMRIPrep*'s
         :abbr:`EPI (echo-planar imaging)` -*boldref* template, which
         is in MNI space.
         The tentative mask is obtained by resampling the MNI template's
         brainmask into *boldref*-space.
      2. Binary dilation of the tentative mask with a sphere of 3mm diameter.
      3. Run ANTs' ``N4BiasFieldCorrection`` on the input
         :abbr:`BOLD (blood-oxygen level-dependant)` average, using the
         mask generated in 1) instead of the internal Otsu thresholding.
      4. Calculate a loose mask using FSL's ``bet``, with one mathematical morphology
         dilation of one iteration and a sphere of 6mm as structuring element.
      5. Mask the :abbr:`INU (intensity non-uniformity)`-corrected image
         with the latest mask calculated in 3), then use AFNI's ``3dUnifize``
         to *standardize* the T2* contrast distribution.
      6. Calculate a mask using AFNI's ``3dAutomask`` after the contrast
         enhancement of 4).
      7. Calculate a final mask as the intersection of 4) and 6).
      8. Apply final mask on the enhanced reference.

    Step 1 can be skipped if the ``pre_mask`` argument is set to ``True`` and
    a tentative mask is passed in to the workflow through the ``pre_mask``
    Nipype input.


    Workflow graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from niworkflows.func.util import init_enhance_and_skullstrip_bold_wf
            wf = init_enhance_and_skullstrip_bold_wf(omp_nthreads=1)

    .. _N4BiasFieldCorrection: https://hdl.handle.net/10380/3053

    Parameters
    ----------
    brainmask_thresh: :obj:`float`
        Lower threshold for the probabilistic brainmask to obtain
        the final binary mask (default: 0.5).
    name : str
        Name of workflow (default: ``enhance_and_skullstrip_bold_wf``)
    omp_nthreads : int
        number of threads available to parallel nodes
    pre_mask : bool
        Indicates whether the ``pre_mask`` input will be set (and thus, step 1
        should be skipped).

    Inputs
    ------
    in_file : str
        BOLD image (single volume)
    pre_mask : bool
        A tentative brain mask to initialize the workflow (requires ``pre_mask``
        parameter set ``True``).


    Outputs
    -------
    bias_corrected_file : str
        the ``in_file`` after `N4BiasFieldCorrection`_
    skull_stripped_file : str
        the ``bias_corrected_file`` after skull-stripping
    mask_file : str
        mask of the skull-stripped input file
    out_report : str
        reportlet for the skull-stripping

    """
    from niworkflows.interfaces.fixes import FixN4BiasFieldCorrection as N4BiasFieldCorrection

    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file', 'pre_mask']), name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['mask_file', 'skull_stripped_file', 'bias_corrected_file']),
        name='outputnode',
    )

    n4_buffer = pe.Node(niu.IdentityInterface(fields=['bias_corrected_file']), name='n4_buffer')

    if bias_correct:
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
        ])  # fmt: skip
    else:
        workflow.connect([(inputnode, n4_buffer, [('in_file', 'bias_corrected_file')])])

    workflow.connect([(n4_buffer, outputnode, [('bias_corrected_file', 'bias_corrected_file')])])

    skullstrip = pe.Node(BrainChop(), name='skullstrip')
    workflow.connect([
        (n4_buffer, skullstrip, [('bias_corrected_file', 'in_file')]),
        (skullstrip, outputnode, [
            ('skullstripped_file', 'skull_stripped_file'),
            ('mask_file', 'mask_file'),
        ]),
    ])  # fmt: skip

    return workflow


def init_skullstrip_asl_wf(
    name='skullstrip_asl_wf',
):
    """
    Run brain extraction on a ASL reference image.

    This workflow takes in a ASL reference image and runs brain extraction
    using SynthStrip.

    Steps of this workflow are:

      1. Calculate a tentative mask by registering (9-parameters) to *fMRIPrep*'s
         :abbr:`EPI (echo-planar imaging)` -*boldref* template, which
         is in MNI space.
         The tentative mask is obtained by resampling the MNI template's
         brainmask into *boldref*-space.
      2. Binary dilation of the tentative mask with a sphere of 3mm diameter.
      3. Run ANTs' ``N4BiasFieldCorrection`` on the input
         :abbr:`BOLD (blood-oxygen level-dependant)` average, using the
         mask generated in 1) instead of the internal Otsu thresholding.
      4. Calculate a loose mask using FSL's ``bet``, with one mathematical morphology
         dilation of one iteration and a sphere of 6mm as structuring element.
      5. Mask the :abbr:`INU (intensity non-uniformity)`-corrected image
         with the latest mask calculated in 3), then use AFNI's ``3dUnifize``
         to *standardize* the T2* contrast distribution.
      6. Calculate a mask using AFNI's ``3dAutomask`` after the contrast
         enhancement of 4).
      7. Calculate a final mask as the intersection of 4) and 6).
      8. Apply final mask on the enhanced reference.

    Step 1 can be skipped if the ``pre_mask`` argument is set to ``True`` and
    a tentative mask is passed in to the workflow through the ``pre_mask``
    Nipype input.


    Workflow graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from niworkflows.func.util import init_enhance_and_skullstrip_bold_wf
            wf = init_enhance_and_skullstrip_bold_wf(omp_nthreads=1)

    .. _N4BiasFieldCorrection: https://hdl.handle.net/10380/3053

    Parameters
    ----------
    brainmask_thresh: :obj:`float`
        Lower threshold for the probabilistic brainmask to obtain
        the final binary mask (default: 0.5).
    name : str
        Name of workflow (default: ``enhance_and_skullstrip_bold_wf``)
    omp_nthreads : int
        number of threads available to parallel nodes
    pre_mask : bool
        Indicates whether the ``pre_mask`` input will be set (and thus, step 1
        should be skipped).

    Inputs
    ------
    in_file : str
        BOLD image (single volume)
    pre_mask : bool
        A tentative brain mask to initialize the workflow (requires ``pre_mask``
        parameter set ``True``).


    Outputs
    -------
    bias_corrected_file : str
        the ``in_file`` after `N4BiasFieldCorrection`_
    skull_stripped_file : str
        the ``bias_corrected_file`` after skull-stripping
    mask_file : str
        mask of the skull-stripped input file
    out_report : str
        reportlet for the skull-stripping

    """
    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file']), name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['mask_file', 'skull_stripped_file']),
        name='outputnode',
    )

    skullstrip = pe.Node(BrainChop(), name='skullstrip')
    workflow.connect([
        (inputnode, skullstrip, [('in_file', 'in_file')]),
        (skullstrip, outputnode, [
            ('skullstripped_file', 'skull_stripped_file'),
            ('mask_file', 'mask_file'),
        ]),
    ])  # fmt: skip

    return workflow
