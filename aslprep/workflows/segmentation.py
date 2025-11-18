# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# Many portions of this were copied from pre-apache fmriprep, but some are
# copied post-apache.
#
# Copyright The NiPreps Developers <nipreps@gmail.com>
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
"""
Anatomical reference preprocessing workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_anat_preproc_wf
.. autofunction:: init_skullstrip_ants_wf

"""

from nipype.interfaces import afni, ants
from nipype.interfaces import utility as niu
from nipype.interfaces.ants import BrainExtraction
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow

from .. import config
from ..interfaces.freesurfer import (
    FixHeaderSynthStrip,
    MockSynthStrip,
    PrepareSynthStripGrid,
)

ANTS_VERSION = BrainExtraction().version or '<ver>'
FS_VERSION = '7.3.1'


def init_dl_prep_wf(name='dl_prep_wf') -> Workflow:
    """Prepare images for use in the FreeSurfer deep learning functions"""
    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['image']), name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['padded_image']),
        name='outputnode',
    )
    skulled_1mm_resample = pe.Node(
        afni.Resample(outputtype='NIFTI_GZ', voxel_size=(1.0, 1.0, 1.0)),
        name='skulled_1mm_resample',
    )
    skulled_autobox = pe.Node(
        afni.Autobox(outputtype='NIFTI_GZ', padding=3),
        name='skulled_autobox',
    )
    prepare_synthstrip_reference = pe.Node(
        PrepareSynthStripGrid(),
        name='prepare_synthstrip_reference',
    )
    resample_skulled_to_reference = pe.Node(
        ants.ApplyTransforms(
            dimension=3,
            interpolation='BSpline',
            transforms=['identity'],
        ),
        name='resample_skulled_to_reference',
    )

    workflow.connect([
        (inputnode, skulled_1mm_resample, [('image', 'in_file')]),
        (skulled_1mm_resample, skulled_autobox, [('out_file', 'in_file')]),
        (skulled_autobox, prepare_synthstrip_reference, [('out_file', 'input_image')]),
        (prepare_synthstrip_reference, resample_skulled_to_reference, [
            ('prepared_image', 'reference_image'),
        ]),
        (inputnode, resample_skulled_to_reference, [('image', 'input_image')]),
        (resample_skulled_to_reference, outputnode, [('output_image', 'padded_image')]),
    ])  # fmt:skip
    return workflow


def init_synthstrip_wf(do_padding=True, name='synthstrip_wf') -> Workflow:
    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=['padded_image', 'original_image']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['brain_image', 'brain_mask']),
        name='outputnode',
    )

    if not config.execution.sloppy:
        synthstrip = pe.Node(
            FixHeaderSynthStrip(),  # Threads are always fixed to 1 in the run
            name='synthstrip',
            n_procs=config.nipype.omp_nthreads,
        )
    else:
        synthstrip = pe.Node(
            MockSynthStrip(),
            name='mocksynthstrip',
        )

    mask_to_original_grid = pe.Node(
        ants.ApplyTransforms(
            dimension=3, transforms=['identity'], interpolation='NearestNeighbor'
        ),
        name='mask_to_original_grid',
    )
    mask_brain = pe.Node(
        ants.MultiplyImages(dimension=3, output_product_image='masked_brain.nii'),
        name='mask_brain',
    )

    # If the input image isn't already padded, do it here
    if do_padding:
        padding_wf = init_dl_prep_wf(name='pad_before_' + name)
        workflow.connect([
            (inputnode, padding_wf, [('original_image', 'inputnode.image')]),
            (padding_wf, synthstrip, [('outputnode.padded_image', 'input_image')]),
        ])  # fmt:skip
    else:
        workflow.connect([(inputnode, synthstrip, [('padded_image', 'input_image')])])

    workflow.connect([
        (synthstrip, mask_to_original_grid, [('out_brain_mask', 'input_image')]),
        (inputnode, mask_to_original_grid, [('original_image', 'reference_image')]),
        (mask_to_original_grid, outputnode, [('output_image', 'brain_mask')]),
        (inputnode, mask_brain, [('original_image', 'first_input')]),
        (mask_to_original_grid, mask_brain, [('output_image', 'second_input')]),
        (mask_brain, outputnode, [('output_product_image', 'brain_image')]),
    ])  # fmt:skip

    return workflow
