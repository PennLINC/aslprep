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
"""Fit workflows for ASLPrep."""
import typing as ty

import bids
import nibabel as nb
import numpy as np
from fmriprep.workflows.bold.registration import init_bold_reg_wf
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.interfaces.header import ValidateImage
from niworkflows.interfaces.nitransforms import ConcatenateXFMs
from niworkflows.interfaces.utility import KeySelect
from sdcflows.workflows.apply.correction import init_unwarp_wf
from sdcflows.workflows.apply.registration import init_coeff2epi_wf

# ASL workflows
# TODO: Replace with fMRIPrep imports once new release is made
from aslprep import config
from aslprep.interfaces.reports import FunctionalSummary
from aslprep.interfaces.resampling import (
    DistortionParameters,
    ReconstructFieldmap,
    ResampleSeries,
)
from aslprep.interfaces.utility import ReduceASLFiles
from aslprep.utils.asl import select_processing_target
from aslprep.utils.bids import extract_entities  # TODO: Replace with fMRIPrep import
from aslprep.workflows.asl.hmc import init_asl_hmc_wf
from aslprep.workflows.asl.outputs import (
    init_asl_fit_reports_wf,
    init_ds_aslref_wf,
    init_ds_hmc_wf,
    init_ds_registration_wf,
)
from aslprep.workflows.asl.reference import init_raw_aslref_wf
from aslprep.workflows.asl.util import init_enhance_and_skullstrip_asl_wf


def get_sbrefs(
    asl_file: str,
    entity_overrides: ty.Dict[str, ty.Any],
    layout: bids.BIDSLayout,
) -> ty.List[str]:
    """Find single-band reference(s) associated with ASL file.

    Parameters
    ----------
    asl_file
        List of absolute paths to ASL files
    entity_overrides
        Query parameters to override defaults
    layout
        :class:`~bids.layout.BIDSLayout` to query

    Returns
    -------
    sbref_file
        List of absolute paths to sbref files associated with input ASL files,
        sorted by EchoTime
    """
    entities = extract_entities(asl_file)
    entities.update(suffix="sbref", extension=[".nii", ".nii.gz"], **entity_overrides)

    return layout.get(return_type="file", **entities)


def init_asl_fit_wf(
    *,
    asl_file: str,
    m0scan: ty.Union[str, None],
    fieldmap_id: ty.Optional[str] = None,
    omp_nthreads: int = 1,
    name: str = "asl_fit_wf",
) -> pe.Workflow:
    """Control the minimal estimation steps for functional preprocessing.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep import config
            from aslprep.workflows.asl.fit import init_asl_fit_wf

            with mock_config():
                asl_file = (
                    config.execution.bids_dir / "sub-01" / "perf" /
                    "sub-01_asl.nii.gz"
                )
                wf = init_asl_fit_wf(asl_file=str(asl_file))

    Parameters
    ----------
    asl_series
        List of paths to NIfTI files.
    fieldmap_id
        ID of the fieldmap to use to correct this ASL series. If :obj:`None`,
        no correction will be applied.

    Inputs
    ------
    asl_file
        ASL series NIfTI file
    t1w_preproc
        Bias-corrected structural template image
    t1w_mask
        Mask of the skull-stripped template image
    t1w_dseg
        Segmentation of preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF)
    anat2std_xfm
        List of transform files, collated with templates
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    fsnative2t1w_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed subject space to T1w
    fmap_id
        Unique identifiers to select fieldmap files
    fmap
        List of estimated fieldmaps (collated with fmap_id)
    fmap_ref
        List of fieldmap reference files (collated with fmap_id)
    fmap_coeff
        List of lists of spline coefficient files (collated with fmap_id)
    fmap_mask
        List of fieldmap masks (collated with fmap_id)
    sdc_method
        List of fieldmap correction method names (collated with fmap_id)

    Outputs
    -------
    hmc_aslref
        ASL reference image used for head motion correction.
        Minimally processed to ensure consistent contrast with ASL series.
    coreg_aslref
        ASL reference image used for coregistration. Contrast-enhanced
        and fieldmap-corrected for greater anatomical fidelity, and aligned
        with ``hmc_aslref``.
    asl_mask
        Mask of ``coreg_aslref``.
    motion_xfm
        Affine transforms from each ASL volume to ``hmc_aslref``, written
        as concatenated ITK affine transforms.
    aslref2anat_xfm
        Affine transform mapping from ASL reference space to the anatomical
        space.
    aslref2fmap_xfm
        Affine transform mapping from ASL reference space to the fieldmap
        space, if applicable.
    movpar_file
        MCFLIRT motion parameters, normalized to SPM format (X, Y, Z, Rx, Ry, Rz)
    rmsd_file
        Root mean squared deviation as measured by ``fsl_motion_outliers`` [Jenkinson2002]_.
    dummy_scans
        The number of dummy scans declared or detected at the beginning of the series.
    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    from aslprep.utils.misc import estimate_asl_mem_usage

    layout = config.execution.layout

    # Collect asl and sbref files, sorted by EchoTime
    sbref_file = get_sbrefs(
        asl_file,
        entity_overrides=config.execution.get().get("bids_filters", {}).get("sbref", {}),
        layout=layout,
    )

    # Get metadata from ASL file(s)
    metadata = layout.get_metadata(asl_file)
    # Patch RepetitionTimePreparation into RepetitionTime,
    # for the sake of BOLD-based interfaces and workflows.
    # This value shouldn't be used for anything except figures and reportlets.
    metadata["RepetitionTime"] = metadata.get(
        "RepetitionTime",
        np.mean(metadata["RepetitionTimePreparation"]),
    )

    orientation = "".join(nb.aff2axcodes(nb.load(asl_file).affine))

    _, mem_gb = estimate_asl_mem_usage(asl_file)

    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "aslcontext",
                # Fieldmap registration
                "fmap",
                "fmap_ref",
                "fmap_coeff",
                "fmap_mask",
                "fmap_id",
                "sdc_method",
                # Anatomical coregistration
                "t1w_preproc",
                "t1w_mask",
                "t1w_dseg",
                "subjects_dir",
                "subject_id",
                "fsnative2t1w_xfm",
                # Other things
                "dummy_scans",
            ],
        ),
        name="inputnode",
    )
    inputnode.inputs.asl_file = asl_file
    inputnode.inputs.dummy_scans = config.workflow.dummy_scans

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "dummy_scans",
                "hmc_aslref",
                "coreg_aslref",
                "asl_mask",
                "motion_xfm",
                "aslref2anat_xfm",
                "aslref2fmap_xfm",
                "movpar_file",  # motion parameters file, for confounds and plots
                "rmsd_file",
            ],
        ),
        name="outputnode",
    )

    # If all derivatives exist, inputnode could go unconnected, so add explicitly
    workflow.add_nodes([inputnode])

    hmcref_buffer = pe.Node(
        niu.IdentityInterface(fields=["aslref", "asl_file", "dummy_scans"]),
        name="hmcref_buffer",
    )
    fmapref_buffer = pe.Node(niu.Function(function=_select_ref), name="fmapref_buffer")
    hmc_buffer = pe.Node(niu.IdentityInterface(fields=["hmc_xforms"]), name="hmc_buffer")
    fmapreg_buffer = pe.Node(
        niu.IdentityInterface(fields=["aslref2fmap_xfm"]),
        name="fmapreg_buffer",
    )
    regref_buffer = pe.Node(
        niu.IdentityInterface(fields=["aslref", "aslmask"]),
        name="regref_buffer",
    )

    summary = pe.Node(
        FunctionalSummary(
            distortion_correction="None",  # Can override with connection
            registration=("FSL", "FreeSurfer")[config.workflow.run_reconall],
            registration_dof=config.workflow.asl2t1w_dof,
            registration_init=config.workflow.asl2t1w_init,
            pe_direction=metadata.get("PhaseEncodingDirection"),
            tr=metadata["RepetitionTime"],
            orientation=orientation,
        ),
        name="summary",
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )
    # workflow.connect([(inputnode, summary, [("dummy_scans", "dummy_scans")])])

    asl_fit_reports_wf = init_asl_fit_reports_wf(
        sdc_correction=not (fieldmap_id is None),
        freesurfer=config.workflow.run_reconall,
        output_dir=config.execution.aslprep_dir,
    )

    # fmt:off
    workflow.connect([
        # XXX: Was from hmc_aslref_wf
        (inputnode, hmcref_buffer, [("dummy_scans", "dummy_scans")]),
        (hmcref_buffer, outputnode, [
            ("aslref", "hmc_aslref"),
            ("dummy_scans", "dummy_scans"),
        ]),
        (regref_buffer, outputnode, [
            ("aslref", "coreg_aslref"),
            ("aslmask", "asl_mask"),
        ]),
        (fmapreg_buffer, outputnode, [("aslref2fmap_xfm", "aslref2fmap_xfm")]),
        (hmc_buffer, outputnode, [
            ("hmc_xforms", "motion_xfm"),
            ("movpar_file", "movpar_file"),
            ("rmsd_file", "rmsd_file"),
        ]),
        (inputnode, asl_fit_reports_wf, [
            ("asl_file", "inputnode.source_file"),
            ("t1w_preproc", "inputnode.t1w_preproc"),
            # May not need all of these
            ("t1w_mask", "inputnode.t1w_mask"),
            ("t1w_dseg", "inputnode.t1w_dseg"),
            ("subjects_dir", "inputnode.subjects_dir"),
            ("subject_id", "inputnode.subject_id"),
        ]),
        (outputnode, asl_fit_reports_wf, [
            ("coreg_aslref", "inputnode.coreg_aslref"),
            ("aslref2anat_xfm", "inputnode.aslref2anat_xfm"),
        ]),
        (summary, asl_fit_reports_wf, [("out_report", "inputnode.summary_report")]),
    ])
    # fmt:on

    # Stage 1: Generate motion correction aslref
    config.loggers.workflow.info("Stage 1: Adding HMC aslref workflow")
    hmc_aslref_wf = init_raw_aslref_wf(
        name="hmc_aslref_wf",
        asl_file=asl_file,
        m0scan=(metadata["M0Type"] == "Separate"),
    )
    hmc_aslref_wf.inputs.inputnode.m0scan = m0scan
    hmc_aslref_wf.inputs.inputnode.dummy_scans = config.workflow.dummy_scans

    workflow.connect([(inputnode, hmc_aslref_wf, [("aslcontext", "inputnode.aslcontext")])])

    ds_hmc_aslref_wf = init_ds_aslref_wf(
        bids_root=layout.root,
        output_dir=config.execution.aslprep_dir,
        desc="hmc",
        name="ds_hmc_aslref_wf",
    )
    ds_hmc_aslref_wf.inputs.inputnode.source_files = [asl_file]

    # fmt:off
    workflow.connect([
        (hmc_aslref_wf, hmcref_buffer, [
            ("outputnode.asl_file", "asl_file"),
            ("outputnode.aslref", "aslref"),
        ]),
        (hmcref_buffer, ds_hmc_aslref_wf, [("aslref", "inputnode.aslref")]),
        (hmc_aslref_wf, asl_fit_reports_wf, [
            ("outputnode.validation_report", "inputnode.validation_report"),
        ]),
    ])
    # fmt:on

    # Reduce the ASL series to only include volumes that need to be processed.
    processing_target = pe.Node(
        niu.Function(
            function=select_processing_target,
            input_names=["aslcontext"],
            output_names=["processing_target"],
        ),
        name="processing_target",
    )

    reduce_asl_file = pe.Node(
        ReduceASLFiles(metadata=metadata),
        name="reduce_asl_file",
    )
    # fmt:off
    workflow.connect([
        (inputnode, processing_target, [("aslcontext", "aslcontext")]),
        (inputnode, reduce_asl_file, [("aslcontext", "aslcontext")]),
        (processing_target, reduce_asl_file, [("processing_target", "processing_target")]),
        (hmcref_buffer, reduce_asl_file, [("asl_file", "asl_file")]),
    ])
    # fmt:on

    # Stage 2: Estimate head motion
    config.loggers.workflow.info("Stage 2: Adding motion correction workflow")
    asl_hmc_wf = init_asl_hmc_wf(
        name="asl_hmc_wf",
        mem_gb=mem_gb["filesize"],
        omp_nthreads=omp_nthreads,
    )

    ds_hmc_wf = init_ds_hmc_wf(
        bids_root=layout.root,
        output_dir=config.execution.aslprep_dir,
    )
    ds_hmc_wf.inputs.inputnode.source_files = [asl_file]

    # fmt:off
    workflow.connect([
        (hmcref_buffer, asl_hmc_wf, [("aslref", "inputnode.raw_ref_image")]),
        (reduce_asl_file, asl_hmc_wf, [
            ("asl_file", "inputnode.asl_file"),
            ("aslcontext", "inputnode.aslcontext"),
        ]),
        (asl_hmc_wf, ds_hmc_wf, [("outputnode.xforms", "inputnode.xforms")]),
        (asl_hmc_wf, hmc_buffer, [
            ("outputnode.xforms", "hmc_xforms"),
            ("outputnode.movpar_file", "movpar_file"),
            ("outputnode.rmsd_file", "rmsd_file"),
        ]),
    ])
    # fmt:on

    # Stage 3: Create coregistration reference
    # Fieldmap correction only happens during fit if this stage is needed
    config.loggers.workflow.info("Stage 3: Adding coregistration aslref workflow")

    # Select initial aslref, enhance contrast, and generate mask
    fmapref_buffer.inputs.sbref_file = sbref_file
    enhance_aslref_wf = init_enhance_and_skullstrip_asl_wf(pre_mask=False)

    ds_coreg_aslref_wf = init_ds_aslref_wf(
        bids_root=layout.root,
        output_dir=config.execution.aslprep_dir,
        desc="coreg",
        name="ds_coreg_aslref_wf",
    )

    # fmt:off
    workflow.connect([
        (hmcref_buffer, fmapref_buffer, [("aslref", "aslref_files")]),
        (fmapref_buffer, enhance_aslref_wf, [("out", "inputnode.in_file")]),
        (fmapref_buffer, ds_coreg_aslref_wf, [("out", "inputnode.source_files")]),
        (ds_coreg_aslref_wf, regref_buffer, [("outputnode.aslref", "aslref")]),
        (fmapref_buffer, asl_fit_reports_wf, [("out", "inputnode.sdc_aslref")]),
    ])
    # fmt:on

    if fieldmap_id:
        fmap_select = pe.Node(
            KeySelect(
                fields=["fmap_ref", "fmap_coeff", "fmap_mask", "sdc_method"],
                key=fieldmap_id,
            ),
            name="fmap_select",
            run_without_submitting=True,
        )

        fmapreg_wf = init_coeff2epi_wf(
            debug="fieldmaps" in config.execution.debug,
            omp_nthreads=config.nipype.omp_nthreads,
            sloppy=config.execution.sloppy,
            name="fmapreg_wf",
        )

        itk_mat2txt = pe.Node(ConcatenateXFMs(out_fmt="itk"), name="itk_mat2txt")

        ds_fmapreg_wf = init_ds_registration_wf(
            bids_root=layout.root,
            output_dir=config.execution.aslprep_dir,
            source="aslref",
            dest=fieldmap_id.replace("_", ""),
            name="ds_fmapreg_wf",
        )

        # fmt:off
        workflow.connect([
            (enhance_aslref_wf, fmapreg_wf, [
                ("outputnode.bias_corrected_file", "inputnode.target_ref"),
                ("outputnode.mask_file", "inputnode.target_mask"),
            ]),
            (fmap_select, fmapreg_wf, [
                ("fmap_ref", "inputnode.fmap_ref"),
                ("fmap_mask", "inputnode.fmap_mask"),
            ]),
            (fmapreg_wf, itk_mat2txt, [("outputnode.target2fmap_xfm", "in_xfms")]),
            (itk_mat2txt, ds_fmapreg_wf, [("out_xfm", "inputnode.xform")]),
            (fmapref_buffer, ds_fmapreg_wf, [("out", "inputnode.source_files")]),
            (ds_fmapreg_wf, fmapreg_buffer, [("outputnode.xform", "aslref2fmap_xfm")]),
        ])
        # fmt:on

        unwarp_wf = init_unwarp_wf(
            free_mem=config.environment.free_mem,
            debug="fieldmaps" in config.execution.debug,
            omp_nthreads=config.nipype.omp_nthreads,
        )
        unwarp_wf.inputs.inputnode.metadata = layout.get_metadata(asl_file)

        # fmt:off
        workflow.connect([
            (inputnode, fmap_select, [
                ("fmap_ref", "fmap_ref"),
                ("fmap_coeff", "fmap_coeff"),
                ("fmap_mask", "fmap_mask"),
                ("sdc_method", "sdc_method"),
                ("fmap_id", "keys"),
            ]),
            (fmap_select, unwarp_wf, [("fmap_coeff", "inputnode.fmap_coeff")]),
            (fmapreg_buffer, unwarp_wf, [
                # This looks backwards, but unwarp_wf describes transforms in
                # terms of points while we (and init_coeff2epi_wf) describe them
                # in terms of images. Mapping fieldmap coordinates into aslref
                # coordinates maps the aslref image onto the fieldmap image.
                ("aslref2fmap_xfm", "inputnode.fmap2data_xfm"),
            ]),
            (enhance_aslref_wf, unwarp_wf, [
                ("outputnode.bias_corrected_file", "inputnode.distorted"),
            ]),
            (unwarp_wf, ds_coreg_aslref_wf, [("outputnode.corrected", "inputnode.aslref")]),
            (unwarp_wf, regref_buffer, [("outputnode.corrected_mask", "aslmask")]),
            (fmap_select, asl_fit_reports_wf, [("fmap_ref", "inputnode.fmap_ref")]),
            (fmap_select, summary, [("sdc_method", "distortion_correction")]),
            (fmapreg_buffer, asl_fit_reports_wf, [
                ("aslref2fmap_xfm", "inputnode.aslref2fmap_xfm"),
            ]),
            (unwarp_wf, asl_fit_reports_wf, [("outputnode.fieldmap", "inputnode.fieldmap")]),
        ])
        # fmt:on
    else:
        # fmt:off
        workflow.connect([
            (enhance_aslref_wf, ds_coreg_aslref_wf, [
                ("outputnode.bias_corrected_file", "inputnode.aslref"),
            ]),
            (enhance_aslref_wf, regref_buffer, [("outputnode.mask_file", "aslmask")]),
        ])
        # fmt:on

    # calculate ASL registration to T1w
    asl_reg_wf = init_bold_reg_wf(
        bold2t1w_dof=config.workflow.asl2t1w_dof,
        bold2t1w_init=config.workflow.asl2t1w_init,
        freesurfer=config.workflow.run_reconall,
        mem_gb=mem_gb["resampled"],
        name="asl_reg_wf",
        omp_nthreads=omp_nthreads,
        sloppy=config.execution.sloppy,
        use_bbr=config.workflow.use_bbr,
        use_compression=False,
        write_report=False,
    )

    ds_aslreg_wf = init_ds_registration_wf(
        bids_root=layout.root,
        output_dir=config.execution.aslprep_dir,
        source="aslref",
        dest="T1w",
        name="ds_aslreg_wf",
    )

    # fmt:off
    workflow.connect([
        (inputnode, asl_reg_wf, [
            ("t1w_preproc", "inputnode.t1w_brain"),  # TODO: rename to t1w_preproc for next
            # ("t1w_mask", "inputnode.t1w_mask"),  # TODO: Add back for next
            ("t1w_dseg", "inputnode.t1w_dseg"),
            # Undefined if --fs-no-reconall, but this is safe
            ("subjects_dir", "inputnode.subjects_dir"),
            ("subject_id", "inputnode.subject_id"),
            ("fsnative2t1w_xfm", "inputnode.fsnative2t1w_xfm"),
        ]),
        (regref_buffer, asl_reg_wf, [("aslref", "inputnode.ref_bold_brain")]),
        # Incomplete sources
        (regref_buffer, ds_aslreg_wf, [("aslref", "inputnode.source_files")]),
        (asl_reg_wf, ds_aslreg_wf, [("outputnode.itk_bold_to_t1", "inputnode.xform")]),
        (ds_aslreg_wf, outputnode, [("outputnode.xform", "aslref2anat_xfm")]),
        (asl_reg_wf, summary, [("outputnode.fallback", "fallback")]),
    ])
    # fmt:on

    return workflow


def init_asl_native_wf(
    *,
    asl_file: str,
    fieldmap_id: ty.Optional[str] = None,
    omp_nthreads: int = 1,
    name: str = "asl_native_wf",
) -> pe.Workflow:
    r"""Apply minimal resampling workflow.

    This workflow resamples to aslref space with head motion and susceptibility distortion
    correction. It also selects the transforms needed to perform further resampling.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep import config
            from aslprep.workflows.asl.fit import init_asl_native_wf

            with mock_config():
                asl_file = (
                    config.execution.bids_dir / "sub-01" / "perf" /
                    "sub-01_asl.nii.gz"
                )
                wf = init_asl_native_wf(asl_file=str(asl_file))

    Parameters
    ----------
    asl_file
        List of paths to NIfTI files.
    fieldmap_id
        ID of the fieldmap to use to correct this ASL series. If :obj:`None`,
        no correction will be applied.

    Inputs
    ------
    aslref
        ASL reference file
    asl_mask
        Mask of ASL reference file
    m0scan
        If M0Type is 'Separate', then the M0 file will resampled into ASL reference space.
        Otherwise, this field will be undefined.
    motion_xfm
        Affine transforms from each ASL volume to ``hmc_aslref``, written
        as concatenated ITK affine transforms.
    aslref2fmap_xfm
        Affine transform mapping from ASL reference space to the fieldmap
        space, if applicable.
    fmap_id
        Unique identifiers to select fieldmap files
    fmap_ref
        List of fieldmap reference files (collated with fmap_id)
    fmap_coeff
        List of lists of spline coefficient files (collated with fmap_id)

    Outputs
    -------
    asl_minimal
        ASL series ready for further resampling.
        This ASL file will only contain the volumes needed for processing.
    asl_native
        ASL series resampled into ASL reference space.
        Head motion and susceptibility distortion correction will be applied to each file.
    m0scan_native
        If M0Type is 'Separate', then the M0 file will resampled into ASL reference space.
        Otherwise, this field will be undefined.
    metadata
        Metadata dictionary of ASL series with the shortest echo
    motion_xfm
        Motion correction transforms for further correcting asl_minimal.
        For multi-echo data, motion correction has already been applied, so
        this will be undefined.
    """
    layout = config.execution.layout
    metadata = layout.get_metadata(asl_file)

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "aslcontext",
                # ASL fit
                "aslref",
                "asl_mask",
                "m0scan",
                "motion_xfm",
                "aslref2fmap_xfm",
                "dummy_scans",
                # Fieldmap fit
                "fmap_ref",
                "fmap_coeff",
                "fmap_id",
            ],
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_minimal",
                "asl_native",
                "m0scan_native",
                "aslcontext",
                "metadata",
                # Transforms
                "motion_xfm",
            ],
        ),
        name="outputnode",
    )

    aslbuffer = pe.Node(
        niu.IdentityInterface(fields=["asl_file", "ro_time", "pe_dir"]),
        name="aslbuffer",
    )

    # Validate the ASL file
    validate_asl = pe.Node(ValidateImage(in_file=asl_file), name="validate_asl")

    # Drop volumes in the ASL file that won't be used
    # (e.g., precalculated CBF volumes if control-label pairs are available).
    processing_target = pe.Node(
        niu.Function(
            function=select_processing_target,
            input_names=["aslcontext"],
            output_names=["processing_target"],
        ),
        name="processing_target",
    )
    reduce_asl_file = pe.Node(
        ReduceASLFiles(metadata=metadata),
        name="reduce_asl_file",
    )
    # fmt:off
    workflow.connect([
        (inputnode, processing_target, [("aslcontext", "aslcontext")]),
        (processing_target, reduce_asl_file, [("processing_target", "processing_target")]),
        (inputnode, reduce_asl_file, [("aslcontext", "aslcontext")]),
        (validate_asl, reduce_asl_file, [("out_file", "asl_file")]),
        (reduce_asl_file, aslbuffer, [("asl_file", "asl_file")]),
        (reduce_asl_file, outputnode, [
            ("aslcontext", "aslcontext"),
            ("metadata", "metadata"),
        ]),
    ])
    # fmt:on

    # Prepare fieldmap metadata
    if fieldmap_id:
        fmap_select = pe.Node(
            KeySelect(fields=["fmap_ref", "fmap_coeff"], key=fieldmap_id),
            name="fmap_select",
            run_without_submitting=True,
        )

        distortion_params = pe.Node(
            DistortionParameters(),
            name="distortion_params",
            run_without_submitting=True,
        )
        workflow.connect([
            (inputnode, fmap_select, [
                ("fmap_ref", "fmap_ref"),
                ("fmap_coeff", "fmap_coeff"),
                ("fmap_id", "keys"),
            ]),
            (reduce_asl_file, distortion_params, [
                ("metadata", "metadata"),
                ("asl_file", "in_file"),
            ]),
            (distortion_params, aslbuffer, [
                ("readout_time", "ro_time"),
                ("pe_direction", "pe_dir"),
            ]),
        ])  # fmt:skip

    # Resample to aslref
    aslref_asl = pe.Node(ResampleSeries(), name="aslref_asl", n_procs=omp_nthreads)

    workflow.connect([
        (inputnode, aslref_asl, [
            ("aslref", "ref_file"),
            ("motion_xfm", "transforms"),
        ]),
        (aslbuffer, aslref_asl, [
            ("asl_file", "in_file"),
            ("ro_time", "ro_time"),
            ("pe_dir", "pe_dir"),
        ]),
    ])  # fmt:skip

    if fieldmap_id:
        aslref_fmap = pe.Node(ReconstructFieldmap(inverse=[True]), name="aslref_fmap")
        workflow.connect([
            (inputnode, aslref_fmap, [
                ("aslref", "target_ref_file"),
                ("aslref2fmap_xfm", "transforms"),
            ]),
            (fmap_select, aslref_fmap, [
                ("fmap_coeff", "in_coeffs"),
                ("fmap_ref", "fmap_ref_file"),
            ]),
            (aslref_fmap, aslref_asl, [("out_file", "fieldmap")]),
        ])  # fmt:skip

    workflow.connect([
        (inputnode, outputnode, [("motion_xfm", "motion_xfm")]),
        (aslbuffer, outputnode, [("asl_file", "asl_minimal")]),
        (aslref_asl, outputnode, [("out_file", "asl_native")]),
    ])  # fmt:skip

    return workflow


def _select_ref(sbref_file, aslref_files):
    """Select first sbref or aslref file, preferring sbref if available."""
    from niworkflows.utils.connections import listify

    refs = sbref_file or aslref_files
    return listify(refs)[0]
