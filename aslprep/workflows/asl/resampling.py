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
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe


def init_bold_fsLR_resampling_wf(
    grayord_density: ty.Literal["91k", "170k"],
    omp_nthreads: int,
    mem_gb: float,
    name: str = "bold_fsLR_resampling_wf",
):
    """Resample BOLD time series to fsLR surface.

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
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.utility import KeySelect
    from smriprep import data as smriprep_data

    from fmriprep.interfaces.workbench import VolumeToSurfaceMapping

    fslr_density = "32k" if grayord_density == "91k" else "59k"

    workflow = Workflow(name=name)

    workflow.__desc__ = """\
The BOLD time-series were resampled onto the left/right-symmetric template
"fsLR" using the Connectome Workbench [@hcppipelines].
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_file",
                "white",
                "pial",
                "midthickness",
                "midthickness_fsLR",
                "sphere_reg_fsLR",
                "cortex_mask",
                "volume_roi",
            ]
        ),
        name="inputnode",
    )

    hemisource = pe.Node(
        niu.IdentityInterface(fields=["hemi"]),
        name="hemisource",
        iterables=[("hemi", ["L", "R"])],
    )

    joinnode = pe.JoinNode(
        niu.IdentityInterface(fields=["bold_fsLR"]),
        name="joinnode",
        joinsource="hemisource",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=["bold_fsLR"]),
        name="outputnode",
    )

    # select white, midthickness and pial surfaces based on hemi
    select_surfaces = pe.Node(
        KeySelect(
            fields=[
                "white",
                "pial",
                "midthickness",
                "midthickness_fsLR",
                "sphere_reg_fsLR",
                "template_sphere",
                "cortex_mask",
                "template_roi",
            ],
            keys=["L", "R"],
        ),
        name="select_surfaces",
        run_without_submitting=True,
    )
    select_surfaces.inputs.template_sphere = [
        str(sphere)
        for sphere in tf.get(
            template="fsLR",
            density=fslr_density,
            suffix="sphere",
            space=None,
            extension=".surf.gii",
        )
    ]
    atlases = smriprep_data.load_resource("atlases")
    select_surfaces.inputs.template_roi = [
        str(atlases / "L.atlasroi.32k_fs_LR.shape.gii"),
        str(atlases / "R.atlasroi.32k_fs_LR.shape.gii"),
    ]

    # RibbonVolumeToSurfaceMapping.sh
    # Line 85 thru ...
    volume_to_surface = pe.Node(
        VolumeToSurfaceMapping(method="ribbon-constrained"),
        name="volume_to_surface",
        mem_gb=mem_gb * 3,
        n_procs=omp_nthreads,
    )
    metric_dilate = pe.Node(
        MetricDilate(distance=10, nearest=True),
        name="metric_dilate",
        mem_gb=1,
        n_procs=omp_nthreads,
    )
    mask_native = pe.Node(MetricMask(), name="mask_native")
    resample_to_fsLR = pe.Node(
        MetricResample(method="ADAP_BARY_AREA", area_surfs=True),
        name="resample_to_fsLR",
        mem_gb=1,
        n_procs=omp_nthreads,
    )
    # ... line 89
    mask_fsLR = pe.Node(MetricMask(), name="mask_fsLR")

    workflow.connect([
        (inputnode, select_surfaces, [
            ("white", "white"),
            ("pial", "pial"),
            ("midthickness", "midthickness"),
            ("midthickness_fsLR", "midthickness_fsLR"),
            ("sphere_reg_fsLR", "sphere_reg_fsLR"),
            ("cortex_mask", "cortex_mask"),
        ]),
        (hemisource, select_surfaces, [("hemi", "key")]),
        # Resample BOLD to native surface, dilate and mask
        (inputnode, volume_to_surface, [
            ("bold_file", "volume_file"),
            ("volume_roi", "volume_roi"),
        ]),
        (select_surfaces, volume_to_surface, [
            ("midthickness", "surface_file"),
            ("white", "inner_surface"),
            ("pial", "outer_surface"),
        ]),
        (select_surfaces, metric_dilate, [("midthickness", "surf_file")]),
        (select_surfaces, mask_native, [("cortex_mask", "mask")]),
        (volume_to_surface, metric_dilate, [("out_file", "in_file")]),
        (metric_dilate, mask_native, [("out_file", "in_file")]),
        # Resample BOLD to fsLR and mask
        (select_surfaces, resample_to_fsLR, [
            ("sphere_reg_fsLR", "current_sphere"),
            ("template_sphere", "new_sphere"),
            ("midthickness", "current_area"),
            ("midthickness_fsLR", "new_area"),
            ("cortex_mask", "roi_metric"),
        ]),
        (mask_native, resample_to_fsLR, [("out_file", "in_file")]),
        (select_surfaces, mask_fsLR, [("template_roi", "mask")]),
        (resample_to_fsLR, mask_fsLR, [("out_file", "in_file")]),
        # Output
        (mask_fsLR, joinnode, [("out_file", "bold_fsLR")]),
        (joinnode, outputnode, [("bold_fsLR", "bold_fsLR")]),
    ])  # fmt:skip

    return workflow
