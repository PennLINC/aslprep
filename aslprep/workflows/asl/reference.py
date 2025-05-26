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

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
First, a reference volume was generated using a custom methodology of *ASLPrep*,
for use in head motion correction.
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

    if m0scan:
        val_m0scan = pe.Node(
            ValidateImage(),
            name='val_m0scan',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(inputnode, val_m0scan, [('m0scan', 'in_file')])])

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
        workflow.connect([(val_m0scan, select_highest_contrast_volumes, [('out_file', 'm0scan')])])

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
