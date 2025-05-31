# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# Copyright 2023 The NiPreps Developers <nipreps@gmail.com>
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
"""Resampling workflows for ASLPrep.

TODO: Remove once fMRIPrep releases 23.2.0.
"""

from __future__ import annotations

from nipype.interfaces import freesurfer as fs
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.interfaces.freesurfer import MedialNaNs

from aslprep import config
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.bids import DerivativesDataSink


def init_asl_surf_wf(
    *,
    mem_gb: float,
    surface_spaces: list[str],
    medial_surface_nan: bool,
    metadata: dict,  # noqa: U100
    cbf_3d: list[str],
    cbf_4d: list[str],
    att: list[str],
    output_dir: str,
    name: str = 'asl_surf_wf',
):
    """Sample functional images to FreeSurfer surfaces.

    For each vertex, the cortical ribbon is sampled at six points (spaced 20% of thickness apart)
    and averaged.

    Outputs are in GIFTI format.

    The two main changes for ASLPrep are: prepare_timing_parameters is dropped and
    DerivativesDataSink is imported outside the function.
    The former is because prepare_timing_parameters relies on *fMRIPrep's* config,
    which will be uninitialized when called by ASLPrep.
    ASLPrep can work around this with a context manager though.
    The latter is because, when DerivativesDataSink is imported within the function,
    ASLPrep can't use a context manager to override it with its own version.
    TODO: Replace with fMRIPrep workflow once DerivativesDataSink import is moved out.

    I've made a bunch of further changes to write out CBF maps instead.

    Workflow Graph
        .. workflow::
            :graph2use: colored
            :simple_form: yes

            from aslprep.workflows.asl.resampling import init_asl_surf_wf

            wf = init_asl_surf_wf(
                mem_gb=0.1,
                surface_spaces=["fsnative", "fsaverage5"],
                medial_surface_nan=False,
                metadata={},
                output_dir=".",
            )

    Parameters
    ----------
    surface_spaces : :obj:`list`
        List of FreeSurfer surface-spaces (either ``fsaverage{3,4,5,6,}`` or ``fsnative``)
        the functional images are to be resampled to.
        For ``fsnative``, images will be resampled to the individual subject's
        native surface.
    medial_surface_nan : :obj:`bool`
        Replace medial wall values with NaNs on functional GIFTI files

    Inputs
    ------
    source_file
        Original ASL series
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    fsnative2t1w_xfm
        ITK-style affine matrix translating from FreeSurfer-conformed subject space to T1w

    Outputs
    -------
    surfaces
        ASL series, resampled to FreeSurfer surfaces

    """
    from nipype.interfaces.io import FreeSurferSource
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.nitransforms import ConcatenateXFMs
    from niworkflows.interfaces.surf import GiftiSetAnatomicalStructure

    from aslprep.workflows.asl.outputs import (
        BASE_INPUT_FIELDS,
        prepare_timing_parameters,
    )

    timing_parameters = prepare_timing_parameters(metadata)

    workflow = Workflow(name=name)
    out_spaces_str = ', '.join([f'*{s}*' for s in surface_spaces])
    workflow.__desc__ = f"""\
The CBF maps were resampled onto the following surfaces (FreeSurfer reconstruction nomenclature):
{out_spaces_str}.
"""
    inputnode_fields = [
        'source_file',
        'anat',
        'aslref2anat_xfm',
        'subject_id',
        'subjects_dir',
        'fsnative2t1w_xfm',
    ]
    inputnode_fields += cbf_3d
    inputnode_fields += cbf_4d
    inputnode_fields += att
    inputnode = pe.Node(
        niu.IdentityInterface(fields=inputnode_fields),
        name='inputnode',
    )

    itersource = pe.Node(niu.IdentityInterface(fields=['target']), name='itersource')
    itersource.iterables = [('target', surface_spaces)]

    get_fsnative = pe.Node(FreeSurferSource(), name='get_fsnative', run_without_submitting=True)
    workflow.connect([
        (inputnode, get_fsnative, [
            ('subject_id', 'subject_id'),
            ('subjects_dir', 'subjects_dir')
        ]),
    ])  # fmt:skip

    def select_target(subject_id, space):
        """Get the target subject ID, given a source subject ID and a target space."""
        return subject_id if space == 'fsnative' else space

    targets = pe.Node(
        niu.Function(function=select_target),
        name='targets',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([
        (inputnode, targets, [('subject_id', 'subject_id')]),
        (itersource, targets, [('target', 'space')]),
    ])  # fmt:skip

    for cbf_deriv in cbf_4d + cbf_3d + att:
        fields = BASE_INPUT_FIELDS[cbf_deriv]

        kwargs = {}
        if cbf_deriv in cbf_4d:
            kwargs['dimension'] = 3

        warp_cbf_to_anat = pe.Node(
            ApplyTransforms(
                interpolation='LanczosWindowedSinc',
                float=True,
                input_image_type=3,
                args='-v',
                **kwargs,
            ),
            name=f'warp_{cbf_deriv}_to_anat',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([
            (inputnode, warp_cbf_to_anat, [
                (cbf_deriv, 'input_image'),
                ('anat', 'reference_image'),
                ('aslref2anat_xfm', 'transforms'),
            ]),
        ])  # fmt:skip

        itk2lta = pe.Node(
            ConcatenateXFMs(out_fmt='fs', inverse=True),
            name=f'itk2lta_{cbf_deriv}',
            run_without_submitting=True,
        )
        workflow.connect([
            (inputnode, itk2lta, [('fsnative2t1w_xfm', 'in_xfms')]),
            (warp_cbf_to_anat, itk2lta, [('output_image', 'moving')]),
            (get_fsnative, itk2lta, [('T1', 'reference')]),
        ])  # fmt:skip

        sampler = pe.MapNode(
            fs.SampleToSurface(
                interp_method='trilinear',
                out_type='gii',
                override_reg_subj=True,
                sampling_method='average',
                sampling_range=(0, 1, 0.2),
                sampling_units='frac',
            ),
            iterfield=['hemi'],
            name=f'sampler_{cbf_deriv}',
            mem_gb=mem_gb * 3,
        )
        sampler.inputs.hemi = ['lh', 'rh']
        workflow.connect([
            (inputnode, sampler, [
                ('subjects_dir', 'subjects_dir'),
                ('subject_id', 'subject_id'),
            ]),
            (warp_cbf_to_anat, sampler, [('output_image', 'source_file')]),
            (itk2lta, sampler, [('out_inv', 'reg_file')]),
            (targets, sampler, [('out', 'target_subject')]),
        ])  # fmt:skip

        update_metadata = pe.MapNode(
            GiftiSetAnatomicalStructure(),
            iterfield=['in_file'],
            name=f'update_{cbf_deriv}_metadata',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        ds_surfs = pe.MapNode(
            DerivativesDataSink(
                base_directory=output_dir,
                extension='.func.gii',
                **timing_parameters,
                **fields,
            ),
            iterfield=['in_file', 'hemi'],
            name=f'ds_{cbf_deriv}_surfs',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        ds_surfs.inputs.hemi = ['L', 'R']

        workflow.connect([
            (inputnode, ds_surfs, [('source_file', 'source_file')]),
            (itersource, ds_surfs, [('target', 'space')]),
            (update_metadata, ds_surfs, [('out_file', 'in_file')]),
        ])  # fmt:skip

        # Refine if medial vertices should be NaNs
        medial_nans = pe.MapNode(
            MedialNaNs(),
            iterfield=['in_file'],
            name=f'medial_nans_{cbf_deriv}',
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )

        if medial_surface_nan:
            # fmt: off
            workflow.connect([
                (inputnode, medial_nans, [('subjects_dir', 'subjects_dir')]),
                (sampler, medial_nans, [('out_file', 'in_file')]),
                (medial_nans, update_metadata, [('out_file', 'in_file')]),
            ])
            # fmt: on
        else:
            workflow.connect([(sampler, update_metadata, [('out_file', 'in_file')])])

    return workflow
