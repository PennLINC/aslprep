"""Workflows to apply transforms to data."""
from __future__ import annotations

import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe
from fmriprep.interfaces.resampling import (
    DistortionParameters,
    ReconstructFieldmap,
    ResampleSeries,
)
from niworkflows.interfaces.nibabel import GenerateSamplingReference
from niworkflows.interfaces.utility import KeySelect


def init_asl_volumetric_resample_wf(
    *,
    metadata: dict,
    fieldmap_id: str | None = None,
    omp_nthreads: int = 1,
    name: str = "asl_volumetric_resample_wf",
) -> pe.Workflow:
    """Resample a ASL series to a volumetric target space.

    This workflow collates a sequence of transforms to resample a ASL series
    in a single shot, including motion correction and fieldmap correction, if
    requested.

    .. workflow::

        from fmriprep.workflows.asl.resampling import init_asl_volumetric_resample_wf
        wf = init_asl_volumetric_resample_wf(
            metadata={
                "RepetitionTime": 2.0,
                "PhaseEncodingDirection": "j-",
                "TotalReadoutTime": 0.03
            },
            fieldmap_id="my_fieldmap",
        )

    Parameters
    ----------
    metadata
        BIDS metadata for ASL file.
    fieldmap_id
        Fieldmap identifier, if fieldmap correction is to be applied.
    omp_nthreads
        Maximum number of threads an individual process may use.
    name
        Name of workflow (default: ``asl_volumetric_resample_wf``)

    Inputs
    ------
    asl_file
        ASL series to resample.
    asl_ref_file
        Reference image to which ASL series is aligned.
    target_ref_file
        Reference image defining the target space.
    target_mask
        Brain mask corresponding to ``target_ref_file``.
        This is used to define the field of view for the resampled ASL series.
    motion_xfm
        List of affine transforms aligning each volume to ``asl_ref_file``.
        If undefined, no motion correction is performed.
    aslref2fmap_xfm
        Affine transform from ``asl_ref_file`` to the fieldmap reference image.
    fmap_ref
        Fieldmap reference image defining the valid field of view for the fieldmap.
    fmap_coeff
        B-Spline coefficients for the fieldmap.
    fmap_id
        Fieldmap identifier, to select correct fieldmap in case there are multiple.
    aslref2anat_xfm
        Affine transform from ``asl_ref_file`` to the anatomical reference image.
    anat2std_xfm
        Affine transform from the anatomical reference image to standard space.
        Leave undefined to resample to anatomical reference space.

    Outputs
    -------
    asl_file
        The ``asl_file`` input, resampled to ``target_ref_file`` space.
    resampling_reference
        An empty reference image with the correct affine and header for resampling
        further images into the ASL series' space.

    """
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "asl_file",
                "asl_ref_file",
                "target_ref_file",
                "target_mask",
                # HMC
                "motion_xfm",
                # SDC
                "aslref2fmap_xfm",
                "fmap_ref",
                "fmap_coeff",
                "fmap_id",
                # Anatomical
                "aslref2anat_xfm",
                # Template
                "anat2std_xfm",
                # Entity for selecting target resolution
                "resolution",
            ],
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=["asl_file", "resampling_reference"]),
        name="outputnode",
    )

    gen_ref = pe.Node(GenerateSamplingReference(), name="gen_ref", mem_gb=0.3)

    aslref2target = pe.Node(niu.Merge(2), name="aslref2target")
    asl2target = pe.Node(niu.Merge(2), name="asl2target")
    resample = pe.Node(ResampleSeries(), name="resample", n_procs=omp_nthreads)

    workflow.connect([
        (inputnode, gen_ref, [
            ("asl_ref_file", "moving_image"),
            ("target_ref_file", "fixed_image"),
            ("target_mask", "fov_mask"),
            (("resolution", _is_native), "keep_native"),
        ]),
        (inputnode, aslref2target, [
            ("aslref2anat_xfm", "in1"),
            ("anat2std_xfm", "in2"),
        ]),
        (inputnode, asl2target, [("motion_xfm", "in1")]),
        (inputnode, resample, [("asl_file", "in_file")]),
        (gen_ref, resample, [("out_file", "ref_file")]),
        (aslref2target, asl2target, [("out", "in2")]),
        (asl2target, resample, [("out", "transforms")]),
        (gen_ref, outputnode, [("out_file", "resampling_reference")]),
        (resample, outputnode, [("out_file", "asl_file")]),
    ])  # fmt:skip

    if not fieldmap_id:
        return workflow

    fmap_select = pe.Node(
        KeySelect(fields=["fmap_ref", "fmap_coeff"], key=fieldmap_id),
        name="fmap_select",
        run_without_submitting=True,
    )
    distortion_params = pe.Node(
        DistortionParameters(metadata=metadata),
        name="distortion_params",
        run_without_submitting=True,
    )
    fmap2target = pe.Node(niu.Merge(2), name="fmap2target")
    inverses = pe.Node(niu.Function(function=_gen_inverses), name="inverses")

    fmap_recon = pe.Node(ReconstructFieldmap(), name="fmap_recon")

    workflow.connect([
        (inputnode, fmap_select, [
            ("fmap_ref", "fmap_ref"),
            ("fmap_coeff", "fmap_coeff"),
            ("fmap_id", "keys"),
        ]),
        (inputnode, distortion_params, [("asl_file", "in_file")]),
        (inputnode, fmap2target, [("aslref2fmap_xfm", "in1")]),
        (gen_ref, fmap_recon, [("out_file", "target_ref_file")]),
        (aslref2target, fmap2target, [("out", "in2")]),
        (aslref2target, inverses, [("out", "inlist")]),
        (fmap_select, fmap_recon, [
            ("fmap_coeff", "in_coeffs"),
            ("fmap_ref", "fmap_ref_file"),
        ]),
        (fmap2target, fmap_recon, [("out", "transforms")]),
        (inverses, fmap_recon, [("out", "inverse")]),
        # Inject fieldmap correction into resample node
        (distortion_params, resample, [
            ("readout_time", "ro_time"),
            ("pe_direction", "pe_dir"),
        ]),
        (fmap_recon, resample, [("out_file", "fieldmap")]),
    ])  # fmt:skip

    return workflow


def _gen_inverses(inlist: list) -> list[bool]:
    """Create a list indicating the first transform should be inverted.

    The input list is the collection of transforms that follow the
    inverted one.
    """
    from niworkflows.utils.connections import listify

    if not inlist:
        return [True]
    return [True] + [False] * len(listify(inlist))


def _is_native(value):
    return value == "native"
