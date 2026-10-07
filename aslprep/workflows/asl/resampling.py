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
    metadata: dict,
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
    source_files
        Files to list as Sources in the output metadata
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

    from aslprep.interfaces.bids import BIDSURI
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
        'source_files',
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

    sources = pe.Node(
        BIDSURI(
            numinputs=3,
            dataset_links=config.execution.dataset_links,
            out_dir=str(output_dir),
        ),
        name='sources',
    )
    workflow.connect([
        (inputnode, sources, [
            ('source_files', 'in1'),
            ('aslref2anat_xfm', 'in2'),
            ('fsnative2t1w_xfm', 'in3'),
        ]),
    ])  # fmt:skip

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
            (sources, ds_surfs, [('out', 'Sources')]),
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


def init_asl_wb_surf_wf(
    *,
    source_file: str,
    surface_spaces: list,
    metadata: dict,
    output_dir: str,
    cbf_3d: list[str],
    cbf_4d: list[str],
    att: list[str],
    omp_nthreads: int,
    mem_gb: float,
    name: str = 'asl_wb_surf_wf',
):
    """Resample CBF derivatives to surface templates using the Connectome Workbench.

    This is the ASL counterpart of the Workbench-based surface resampling in fMRIPrep's
    ``init_bold_wf`` (nipreps/fmriprep#3461).
    Each derivative is warped to the anatomical reference, sampled onto the subject's native
    surface with the "ribbon-constrained" method, dilated, and then resampled to each surface
    template whose spheres are registered to fsLR.

    Workflow Graph
        .. workflow::
            :graph2use: colored
            :simple_form: yes

            from niworkflows.utils.spaces import Reference

            from aslprep.tests.tests import mock_config
            from aslprep.workflows.asl.resampling import init_asl_wb_surf_wf

            with mock_config():
                wf = init_asl_wb_surf_wf(
                    source_file='sub-01_asl.nii.gz',
                    surface_spaces=[Reference('fsLR', {'den': '32k'})],
                    metadata={'RepetitionTime': 4.0},
                    output_dir='.',
                    cbf_3d=['mean_cbf'],
                    cbf_4d=[],
                    att=[],
                    omp_nthreads=1,
                    mem_gb=1,
                )

    Parameters
    ----------
    source_file
        Original ASL series, used to name the outputs.
    surface_spaces
        :class:`~niworkflows.utils.spaces.Reference` objects for the target surface templates.
        Each one must specify a density.
    metadata
        BIDS metadata for the ASL series.
    output_dir
        Directory in which to save derivatives.
    cbf_3d, cbf_4d, att
        Names of the CBF derivatives to resample.
    omp_nthreads
        Maximum number of threads an individual process may use.
    mem_gb
        Size of the resampled ASL file in GB.
    name
        Name of workflow (default: ``asl_wb_surf_wf``).

    Inputs
    ------
    source_files
        Files to list as Sources in the output metadata.
    anat_ref_file
        ASL-resolution reference image in anatomical space.
    aslref2anat_xfm
        Affine transform from the ASL reference to the anatomical reference.
    white, pial, midthickness
        Left and right hemisphere GIFTI surfaces.
    sphere_reg_fsLR
        Left and right hemisphere registration spheres to fsLR.
    goodvoxels_mask
        Pre-computed goodvoxels mask in anatomical space. Only used if
        ``--project-goodvoxels`` is enabled.
    """
    import templateflow.api as tf
    from fmriprep.workflows.bold.resampling import init_wb_surf_surf_wf, init_wb_vol_surf_wf
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from smriprep.workflows.surfaces import init_resample_surfaces_wf

    from aslprep.interfaces.bids import BIDSURI
    from aslprep.workflows.asl.outputs import BASE_INPUT_FIELDS, prepare_timing_parameters

    targets = []
    for ref in surface_spaces:
        density = ref.spec.get('density') or ref.spec.get('den')
        if density is None:
            config.loggers.workflow.warning(
                f'Cannot resample to surface space {ref} without a density. Skipping.'
            )
            continue
        targets.append((ref.space, density))

    workflow = Workflow(name=name)

    template_strs = []
    for template, density in targets:
        template_meta = tf.get_metadata(template)
        template_refs = ['@onavg'] if template == 'onavg' else []
        if template_meta.get('RRID'):
            template_refs.append(f'RRID:{template_meta["RRID"]}')
        template_refs.append(f'TemplateFlow ID: {template}')
        template_strs.append(
            f'*{template_meta.get("Name", template)}* [{"; ".join(template_refs)}] '
            f'({density} density)'
        )
    goodvoxels_str = (
        ', excluding voxels marked by the "goodvoxels" mask,'
        if config.workflow.project_goodvoxels
        else ''
    )
    workflow.__desc__ = f"""\
The CBF maps were resampled onto the native surface of the subject{goodvoxels_str}
using the "ribbon-constrained" method of the Connectome Workbench [@hcppipelines],
dilated by 10 mm, and then resampled to the following surface templates:
{', '.join(template_strs)}.
"""

    inputnode_fields = [
        'source_files',
        'anat_ref_file',
        'aslref2anat_xfm',
        'white',
        'pial',
        'midthickness',
        'sphere_reg_fsLR',
        'goodvoxels_mask',
    ]
    inputnode_fields += cbf_3d
    inputnode_fields += cbf_4d
    inputnode_fields += att
    inputnode = pe.Node(
        niu.IdentityInterface(fields=inputnode_fields),
        name='inputnode',
    )

    sources = pe.Node(
        BIDSURI(
            numinputs=2,
            dataset_links=config.execution.dataset_links,
            out_dir=str(output_dir),
        ),
        name='sources',
    )
    workflow.connect([
        (inputnode, sources, [
            ('source_files', 'in1'),
            ('aslref2anat_xfm', 'in2'),
        ]),
    ])  # fmt:skip

    # Subject midthickness surfaces resampled to each template, shared by all derivatives
    resample_surfaces_wfs = {}
    for template, density in targets:
        resample_surfaces_wf = init_resample_surfaces_wf(
            surfaces=['midthickness'],
            template=template,
            density=density,
            name=f'resample_surfaces_{template}_{density}_wf',
        )
        workflow.connect([
            (inputnode, resample_surfaces_wf, [
                ('midthickness', 'inputnode.midthickness'),
                ('sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR'),
            ]),
        ])  # fmt:skip
        resample_surfaces_wfs[(template, density)] = resample_surfaces_wf

    timing_parameters = prepare_timing_parameters(metadata)
    for cbf_deriv in cbf_4d + cbf_3d + att:
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

        wb_vol_surf_wf = init_wb_vol_surf_wf(
            omp_nthreads=omp_nthreads,
            mem_gb=mem_gb,
            dilate=True,
            name=f'{cbf_deriv}_wb_vol_surf_wf',
        )
        # The parent workflow describes the resampling
        wb_vol_surf_wf.__desc__ = None
        workflow.connect([
            (inputnode, warp_cbf_to_anat, [
                (cbf_deriv, 'input_image'),
                ('anat_ref_file', 'reference_image'),
                ('aslref2anat_xfm', 'transforms'),
            ]),
            (inputnode, wb_vol_surf_wf, [
                ('white', 'inputnode.white'),
                ('pial', 'inputnode.pial'),
                ('midthickness', 'inputnode.midthickness'),
            ]),
            (warp_cbf_to_anat, wb_vol_surf_wf, [('output_image', 'inputnode.bold_file')]),
        ])  # fmt:skip
        if config.workflow.project_goodvoxels:
            workflow.connect([
                (inputnode, wb_vol_surf_wf, [('goodvoxels_mask', 'inputnode.volume_roi')]),
            ])  # fmt:skip

        for template, density in targets:
            wb_surf_surf_wf = init_wb_surf_surf_wf(
                template=template,
                density=density,
                omp_nthreads=omp_nthreads,
                mem_gb=mem_gb,
                name=f'{cbf_deriv}_wb_surf_{template}_{density}_wf',
            )
            wb_surf_surf_wf.__desc__ = None

            ds_cbf_surf = pe.MapNode(
                DerivativesDataSink(
                    source_file=source_file,
                    base_directory=output_dir,
                    space=template,
                    density=density,
                    extension='.func.gii',
                    **timing_parameters,
                    **BASE_INPUT_FIELDS[cbf_deriv],
                ),
                iterfield=['in_file', 'hemi'],
                name=f'ds_{cbf_deriv}_{template}_{density}',
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            ds_cbf_surf.inputs.hemi = ['L', 'R']

            workflow.connect([
                (inputnode, wb_surf_surf_wf, [
                    ('midthickness', 'inputnode.midthickness'),
                    ('sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR'),
                ]),
                (wb_vol_surf_wf, wb_surf_surf_wf, [
                    ('outputnode.bold_fsnative', 'inputnode.bold_fsnative'),
                ]),
                (resample_surfaces_wfs[(template, density)], wb_surf_surf_wf, [
                    (f'outputnode.midthickness_{template}', 'inputnode.midthickness_resampled'),
                ]),
                (wb_surf_surf_wf, ds_cbf_surf, [('outputnode.bold_resampled', 'in_file')]),
                (sources, ds_cbf_surf, [('out', 'Sources')]),
            ])  # fmt:skip

    return workflow
