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

import typing as ty

from fmriprep.interfaces.workbench import MetricDilate, MetricMask, MetricResample
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
        Original BOLD series
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    fsnative2t1w_xfm
        ITK-style affine matrix translating from FreeSurfer-conformed subject space to T1w

    Outputs
    -------
    surfaces
        BOLD series, resampled to FreeSurfer surfaces

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

        if cbf_deriv in att:
            meta = {'Units': 's'}
        else:
            meta = {'Units': 'mL/100 g/min'}

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
                **meta,
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


def init_bold_fsLR_resampling_wf(  # noqa: N802
    grayord_density: ty.Literal['91k', '170k'],
    omp_nthreads: int,
    mem_gb: float,
    name: str = 'bold_fsLR_resampling_wf',
):
    """Resample BOLD time series to fsLR surface.

    This is copied from fMRIPrep and can be replaced with an fMRIPrep import once 23.2.0 is
    released.

    This workflow is derived heavily from three scripts within the DCAN-HCP pipelines scripts

    Line numbers correspond to the locations of the code in the original scripts, found at:
    https://github.com/DCAN-Labs/DCAN-HCP/tree/9291324/

    Workflow Graph
        .. workflow::
            :graph2use: colored
            :simple_form: yes

            from fmriprep.workflows.bold.resampling import init_bold_fsLR_resampling_wf
            wf = init_bold_fsLR_resampling_wf(
                grayord_density="92k",
                omp_nthreads=1,
                mem_gb=1,
            )

    Parameters
    ----------
    grayord_density : :class:`str`
        Either ``"91k"`` or ``"170k"``, representing the total *grayordinates*.
    omp_nthreads : :class:`int`
        Maximum number of threads an individual process may use
    mem_gb : :class:`float`
        Size of BOLD file in GB
    name : :class:`str`
        Name of workflow (default: ``bold_fsLR_resampling_wf``)

    Inputs
    ------
    bold_file : :class:`str`
        Path to BOLD file resampled into T1 space
    white : :class:`list` of :class:`str`
        Path to left and right hemisphere white matter GIFTI surfaces.
    pial : :class:`list` of :class:`str`
        Path to left and right hemisphere pial GIFTI surfaces.
    midthickness : :class:`list` of :class:`str`
        Path to left and right hemisphere midthickness GIFTI surfaces.
    midthickness_fsLR : :class:`list` of :class:`str`
        Path to left and right hemisphere midthickness GIFTI surfaces in fsLR space.
    sphere_reg_fsLR : :class:`list` of :class:`str`
        Path to left and right hemisphere sphere.reg GIFTI surfaces, mapping from subject to fsLR
    cortex_mask : :class:`list` of :class:`str`
        Path to left and right hemisphere cortical masks.
    volume_roi : :class:`str` or Undefined
        Pre-calculated goodvoxels mask. Not required.

    Outputs
    -------
    bold_fsLR : :class:`list` of :class:`str`
        Path to BOLD series resampled as functional GIFTI files in fsLR space

    """
    import templateflow.api as tf
    from fmriprep.interfaces.workbench import VolumeToSurfaceMapping
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.utility import KeySelect
    from smriprep import data as smriprep_data

    fslr_density = '32k' if grayord_density == '91k' else '59k'

    workflow = Workflow(name=name)

    workflow.__desc__ = """\
The BOLD time-series were resampled onto the left/right-symmetric template
"fsLR" using the Connectome Workbench [@hcppipelines].
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'bold_file',
                'white',
                'pial',
                'midthickness',
                'midthickness_fsLR',
                'sphere_reg_fsLR',
                'cortex_mask',
                'volume_roi',
            ]
        ),
        name='inputnode',
    )

    hemisource = pe.Node(
        niu.IdentityInterface(fields=['hemi']),
        name='hemisource',
        iterables=[('hemi', ['L', 'R'])],
    )

    joinnode = pe.JoinNode(
        niu.IdentityInterface(fields=['bold_fsLR']),
        name='joinnode',
        joinsource='hemisource',
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['bold_fsLR']),
        name='outputnode',
    )

    # select white, midthickness and pial surfaces based on hemi
    select_surfaces = pe.Node(
        KeySelect(
            fields=[
                'white',
                'pial',
                'midthickness',
                'midthickness_fsLR',
                'sphere_reg_fsLR',
                'template_sphere',
                'cortex_mask',
                'template_roi',
            ],
            keys=['L', 'R'],
        ),
        name='select_surfaces',
        run_without_submitting=True,
    )
    select_surfaces.inputs.template_sphere = [
        str(sphere)
        for sphere in tf.get(
            template='fsLR',
            density=fslr_density,
            suffix='sphere',
            space=None,
            extension='.surf.gii',
        )
    ]
    atlases = smriprep_data.load_resource('atlases')
    select_surfaces.inputs.template_roi = [
        str(atlases / f'L.atlasroi.{fslr_density}_fs_LR.shape.gii'),
        str(atlases / f'R.atlasroi.{fslr_density}_fs_LR.shape.gii'),
    ]

    # RibbonVolumeToSurfaceMapping.sh
    # Line 85 thru ...
    volume_to_surface = pe.Node(
        VolumeToSurfaceMapping(method='ribbon-constrained'),
        name='volume_to_surface',
        mem_gb=mem_gb * 3,
        n_procs=omp_nthreads,
    )
    metric_dilate = pe.Node(
        MetricDilate(distance=10, nearest=True),
        name='metric_dilate',
        mem_gb=1,
        n_procs=omp_nthreads,
    )
    mask_native = pe.Node(MetricMask(), name='mask_native')
    resample_to_fsLR = pe.Node(
        MetricResample(method='ADAP_BARY_AREA', area_surfs=True),
        name='resample_to_fsLR',
        mem_gb=1,
        n_procs=omp_nthreads,
    )
    # ... line 89
    mask_fsLR = pe.Node(MetricMask(), name='mask_fsLR')

    workflow.connect([
        (inputnode, select_surfaces, [
            ('white', 'white'),
            ('pial', 'pial'),
            ('midthickness', 'midthickness'),
            ('midthickness_fsLR', 'midthickness_fsLR'),
            ('sphere_reg_fsLR', 'sphere_reg_fsLR'),
            ('cortex_mask', 'cortex_mask'),
        ]),
        (hemisource, select_surfaces, [('hemi', 'key')]),
        # Resample BOLD to native surface, dilate and mask
        (inputnode, volume_to_surface, [
            ('bold_file', 'volume_file'),
            ('volume_roi', 'volume_roi'),
        ]),
        (select_surfaces, volume_to_surface, [
            ('midthickness', 'surface_file'),
            ('white', 'inner_surface'),
            ('pial', 'outer_surface'),
        ]),
        (select_surfaces, metric_dilate, [('midthickness', 'surf_file')]),
        (select_surfaces, mask_native, [('cortex_mask', 'mask')]),
        (volume_to_surface, metric_dilate, [('out_file', 'in_file')]),
        (metric_dilate, mask_native, [('out_file', 'in_file')]),
        # Resample BOLD to fsLR and mask
        (select_surfaces, resample_to_fsLR, [
            ('sphere_reg_fsLR', 'current_sphere'),
            ('template_sphere', 'new_sphere'),
            ('midthickness', 'current_area'),
            ('midthickness_fsLR', 'new_area'),
            ('cortex_mask', 'roi_metric'),
        ]),
        (mask_native, resample_to_fsLR, [('out_file', 'in_file')]),
        (select_surfaces, mask_fsLR, [('template_roi', 'mask')]),
        (resample_to_fsLR, mask_fsLR, [('out_file', 'in_file')]),
        # Output
        (mask_fsLR, joinnode, [('out_file', 'bold_fsLR')]),
        (joinnode, outputnode, [('bold_fsLR', 'bold_fsLR')]),
    ])  # fmt:skip

    return workflow
