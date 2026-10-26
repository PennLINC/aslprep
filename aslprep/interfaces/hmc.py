# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Interfaces for head-motion correction."""

import json
import os
import subprocess
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    OutputMultiObject,
    SimpleInterface,
    TraitedSpec,
    isdefined,
    traits,
)

#: Delays (seconds; inversion times for PASL) below which labels may still hold arterial signal.
#: At those delays the label mean differs from the control mean by bright vessels as well as
#: perfusion, and registering one to the other reads the vessels as an offset of the head, so
#: their labels are not corrected. Generous on purpose: transit times reach a second or more,
#: and a delay wrongly left uncorrected only keeps the type-wise transforms.
ARTERIAL_DELAY = {'PASL': 0.8, 'default': 1.0}


def plan_label_offset_groups(volume_types, plds, asl_type):
    """Group label volumes by delay, and say which delays are corrected.

    Each delay at or above :data:`ARTERIAL_DELAY` is corrected with the offset estimated from
    its own label and control volumes. Delays below it are not corrected: their labels keep the
    type-wise transforms. Each delay's own estimate is needed because, in a series that acquires
    the delays one after another, it also carries the residual misregistration of that delay's
    pairs; an estimate borrowed from another delay moves the labels by the wrong amount.

    Parameters
    ----------
    volume_types : list of str
        Volume type of each volume of the (reduced) series.
    plds : list of float
        Post-labeling delay (inversion time for PASL) of each volume.
    asl_type : str
        ``ArterialSpinLabelingType``.

    Returns
    -------
    estimated : dict
        Delay -> (label indices, control indices), for the delays whose offset is estimated from
        their own label and control volumes and applied to their labels.
    skipped : list of float
        Delays with labels that are not corrected (below the cutoff, or without controls).
    """
    volume_types = list(volume_types)
    plds = [float(p) for p in plds]
    cut = ARTERIAL_DELAY.get(str(asl_type).upper(), ARTERIAL_DELAY['default'])
    groups = {}
    for d in sorted({p for p, t in zip(plds, volume_types, strict=True) if t == 'label'}):
        labels = [
            i
            for i, (t, p) in enumerate(zip(volume_types, plds, strict=True))
            if t == 'label' and p == d
        ]
        controls = [
            i
            for i, (t, p) in enumerate(zip(volume_types, plds, strict=True))
            if t == 'control' and p == d
        ]
        groups[d] = (labels, controls)
    estimated = {d: g for d, g in groups.items() if d >= cut and g[0] and g[1]}
    skipped = sorted(d for d, g in groups.items() if g[0] and d not in estimated)
    return estimated, skipped


class _CorrectLabelOffsetInputSpec(BaseInterfaceInputSpec):
    asl_file = File(
        exists=True, mandatory=True, desc='ASL series the transforms were estimated on.'
    )
    aslcontext = File(exists=True, mandatory=True, desc='Volume types of asl_file.')
    metadata = traits.Dict(mandatory=True, desc='Metadata of asl_file (reduced to its volumes).')
    ref_file = File(exists=True, mandatory=True, desc='Reference the transforms map into.')
    xforms = File(exists=True, mandatory=True, desc='ITK transforms, one per volume.')


class _CorrectLabelOffsetOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='Corrected ITK transforms, one per volume.')
    out_report = File(exists=True, desc='JSON record of the estimated offsets.')


class CorrectLabelOffset(SimpleInterface):
    """Remove the systematic offset between motion-corrected label and control volumes.

    Volume-type-wise motion correction registers label and control volumes to the same
    reference separately, and leaves a small rigid offset between the corrected labels and
    controls, which shows up at the edges of the brain after subtraction. For each delay at or
    above :data:`ARTERIAL_DELAY`, the mean of the corrected label volumes is registered rigidly
    to the mean of the corrected control volumes (AFNI ``3dvolreg``; at those delays the two
    images should differ mainly by the perfusion signal), and the correction is composed into the
    transforms of that delay's label volumes. Labels at shorter delays, which may still hold
    arterial signal, keep their transforms (:func:`plan_label_offset_groups`), as do control,
    M0 and derived volumes. Series without both labels and controls pass through unchanged.
    """

    input_spec = _CorrectLabelOffsetInputSpec
    output_spec = _CorrectLabelOffsetOutputSpec

    def _run_interface(self, runtime):
        import nitransforms as nt
        from nitransforms.io.afni import AFNILinearTransformArray

        out_file = os.path.join(runtime.cwd, 'hmc_xfms_label_offset.txt')
        out_report = os.path.join(runtime.cwd, 'label_offset.json')
        volume_types = pd.read_table(self.inputs.aslcontext)['volume_type'].tolist()
        img = nb.load(self.inputs.asl_file)
        n_vols = img.shape[3] if img.ndim == 4 else 1
        if len(volume_types) != n_vols:
            raise ValueError(f'{len(volume_types)} volume types for {n_vols} volumes')
        plds = np.atleast_1d(np.asarray(self.inputs.metadata.get('PostLabelingDelay', 0.0), float))
        plds = np.repeat(plds, n_vols) if plds.size == 1 else plds
        if plds.size != n_vols:
            raise ValueError(f'{plds.size} post-labeling delays for {n_vols} volumes')
        asl_type = self.inputs.metadata.get('ArterialSpinLabelingType', 'PCASL')

        ref_img = nb.load(self.inputs.ref_file)
        xfms = nt.linear.load(self.inputs.xforms, fmt='itk', reference=ref_img)
        pulls = np.asarray(xfms.matrix).reshape(-1, 4, 4)  # reference point -> volume point
        if len(pulls) != n_vols:
            raise ValueError(f'{len(pulls)} transforms for {n_vols} volumes')
        estimated, skipped = plan_label_offset_groups(volume_types, plds, asl_type)
        new = pulls.copy()
        record = {'asl_type': str(asl_type), 'delays': {}}
        corrections = {}
        data = np.asanyarray(img.dataobj, np.float32)
        data = data[..., None] if data.ndim == 3 else data
        for d, (labels, controls) in estimated.items():
            means = {}
            for name, sel in (('label', labels), ('control', controls)):
                moved = [_resample(data[..., i], img.affine, pulls[i], ref_img) for i in sel]
                means[name] = os.path.join(runtime.cwd, f'{name}_mean_{d:g}.nii.gz')
                mean_img = nb.Nifti1Image(np.mean(moved, axis=0), ref_img.affine, ref_img.header)
                mean_img.to_filename(means[name])

            aff = os.path.join(runtime.cwd, f'offset_{d:g}.aff12.1D')
            prefix = os.path.join(runtime.cwd, f'offset_{d:g}.nii.gz')
            cmd = ['3dvolreg', '-base', means['control'], '-prefix', prefix]
            cmd += ['-1Dmatrix_save', aff, means['label']]
            subprocess.run(cmd, check=True, capture_output=True)
            # Pull from a control-mean point to the label-mean point (AFNI deobliques implicitly;
            # nitransforms reads its matrices in the plumb frames of the two headers).
            pull = AFNILinearTransformArray.from_filename(aff).to_ras(
                moving=nb.load(means['label']),
                reference=nb.load(means['control']),
            )
            corrections[d] = np.asarray(pull).reshape(-1, 4, 4)[0]
            record['delays'][f'{d:g}'] = {
                'n_labels': len(labels),
                'n_controls': len(controls),
                **_magnitude(corrections[d], ref_img),
            }

        for d in skipped:
            record['delays'][f'{d:g}'] = {'corrected': False}

        for d, (labels, _) in estimated.items():
            for i in labels:
                # The corrected label shows at y what the label mean shows at P y.
                new[i] = pulls[i] @ corrections[d]
        nt.linear.LinearTransformsMapping(new, reference=ref_img).to_filename(out_file, fmt='itk')
        Path(out_report).write_text(json.dumps(record, indent=2) + '\n')
        self._results['out_file'] = out_file
        self._results['out_report'] = out_report
        return runtime


class _IdentityTransformsInputSpec(BaseInterfaceInputSpec):
    aslcontext = File(exists=True, mandatory=True, desc='Volume types of the (reduced) series.')
    ref_file = File(exists=True, mandatory=True, desc='Reference image of the transforms.')


class _IdentityTransformsOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='Identity ITK transforms, one per volume.')


class IdentityTransforms(SimpleInterface):
    """Write an identity transform for every volume, in place of head-motion correction.

    Used for short series without label-control pairs (the GE-specific path), whose few volumes,
    an M0 with delta-M or CBF images, share too little contrast for registration to be more
    accurate than none.
    """

    input_spec = _IdentityTransformsInputSpec
    output_spec = _IdentityTransformsOutputSpec

    def _run_interface(self, runtime):
        import nitransforms as nt

        n_vols = len(pd.read_table(self.inputs.aslcontext)['volume_type'])
        out_file = os.path.join(runtime.cwd, 'hmc_xfms_identity.txt')
        nt.linear.LinearTransformsMapping(
            np.tile(np.eye(4), (n_vols, 1, 1)), reference=nb.load(self.inputs.ref_file)
        ).to_filename(out_file, fmt='itk')
        self._results['out_file'] = out_file
        return runtime


class _VolregSeriesInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='Series (3D or 4D) to register.')
    ref_file = File(
        exists=True,
        desc=(
            'Reference image, on the same grid as in_file. If not given, the middle volume '
            'of in_file (MCFLIRT default).'
        ),
    )


class _VolregSeriesOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='The series resampled onto the reference.')
    mat_file = OutputMultiObject(
        File(exists=True),
        desc='FSL matrices (volume to reference), one per volume, as MCFLIRT.',
    )
    par_file = File(
        exists=True,
        desc=(
            'Rotations (radians) and translations (mm) per volume, in the column order of '
            "MCFLIRT's .par file."
        ),
    )


class VolregSeries(SimpleInterface):
    """Register every volume of a series to a reference with AFNI ``3dvolreg``.

    ``3dvolreg`` fits each volume by least squares, the most precise cost between images of
    one contrast. Its outputs are written as MCFLIRT writes them, so it can take MCFLIRT's place:
    one FSL matrix per volume (``MAT_0000``, ...) and a ``.par`` file. The ``.par`` columns are
    in MCFLIRT's order, but decomposed from the matrices about the reference's geometric centre,
    not MCFLIRT's centre of gravity. They reach only the HMC workflow's ``movpar_file``; the
    motion confounds are computed from the final transforms. The series and the reference must
    share one grid, as MCFLIRT2ITK, downstream, assumes.
    """

    input_spec = _VolregSeriesInputSpec
    output_spec = _VolregSeriesOutputSpec

    def _run_interface(self, runtime):
        from nitransforms.io.afni import AFNILinearTransformArray

        img = nb.load(self.inputs.in_file)
        data = np.asanyarray(img.dataobj, np.float32)
        data = data[..., None] if data.ndim == 3 else data
        series = os.path.join(runtime.cwd, 'series.nii.gz')
        nb.Nifti1Image(data, img.affine, img.header).to_filename(series)
        ref = self.inputs.ref_file if isdefined(self.inputs.ref_file) else None
        if ref is None:
            ref = os.path.join(runtime.cwd, 'middle.nii.gz')
            nb.Nifti1Image(data[..., data.shape[-1] // 2], img.affine).to_filename(ref)

        series_img, ref_img = nb.load(series), nb.load(ref)
        same_shape = series_img.shape[:3] == ref_img.shape[:3]
        same_grid = same_shape and np.allclose(series_img.affine, ref_img.affine, atol=1e-3)
        if not same_grid:
            raise ValueError(
                'The series and the reference must share one grid (shape and affine): the FSL '
                'matrices and motion parameters are written in that one frame'
            )

        out_file = os.path.join(runtime.cwd, 'volreg.nii.gz')
        aff = os.path.join(runtime.cwd, 'volreg.aff12.1D')
        cmd = ['3dvolreg', '-base', ref, '-prefix', out_file, '-1Dmatrix_save', aff, series]
        subprocess.run(cmd, check=True, capture_output=True)
        # Per volume, the RAS pull from a reference point to the volume point (AFNI deobliques
        # implicitly; nitransforms reads its matrices in the plumb frames of the two headers).
        pulls = AFNILinearTransformArray.from_filename(aff).to_ras(
            moving=series_img,
            reference=ref_img,
        )
        pulls = np.asarray(pulls).reshape(-1, 4, 4)

        mat_dir = os.path.join(runtime.cwd, 'volreg.mat')
        os.makedirs(mat_dir, exist_ok=True)
        src_to_world, ref_to_world = _fsl_to_world(series_img), _fsl_to_world(ref_img)
        centre = np.r_[(np.array(ref_img.shape[:3]) - 1) / 2 * ref_img.header.get_zooms()[:3], 1]
        mat_files, params = [], []
        for k, pull in enumerate(pulls):
            mat = pull_to_fsl(pull, src_to_world, ref_to_world)
            mat_files.append(os.path.join(mat_dir, f'MAT_{k:04d}'))
            np.savetxt(mat_files[-1], mat, fmt='%.10f')
            params.append(_fsl_params(mat, centre))
        par_file = os.path.join(runtime.cwd, 'volreg.par')
        np.savetxt(par_file, np.array(params), fmt='%.8f')
        self._results['out_file'] = out_file
        self._results['mat_file'] = mat_files
        self._results['par_file'] = par_file
        return runtime


def _fsl_to_world(img):
    """FSL scaled-voxel coordinates -> the image's world (RAS mm) coordinates."""
    zooms = np.array(img.header.get_zooms()[:3], float)
    scale = np.diag(np.r_[zooms, 1.0])
    if np.linalg.det(img.affine[:3, :3]) > 0:  # FSL flips x for neurological storage
        scale[0, 0] = -zooms[0]
        scale[0, 3] = (img.shape[0] - 1) * zooms[0]
    return img.affine @ np.linalg.inv(scale)


def pull_to_fsl(pull, src_to_world, ref_to_world):
    """A RAS pull (reference point -> volume point) as an FSL matrix (volume -> reference).

    Parameters
    ----------
    pull : (4, 4) array
        RAS matrix mapping a reference point to the corresponding volume point.
    src_to_world, ref_to_world : (4, 4) array
        FSL scaled-voxel to world matrices of the volume and the reference.
    """
    return np.linalg.inv(ref_to_world) @ np.linalg.inv(pull) @ src_to_world


def _fsl_params(mat, centre):
    """Rotations (radians, x y z) and translation at ``centre`` (mm) of an FSL matrix."""
    r = mat[:3, :3]
    rotations = [
        np.arctan2(r[1, 2], r[2, 2]),
        -np.arcsin(np.clip(r[0, 2], -1, 1)),
        np.arctan2(r[0, 1], r[0, 0]),
    ]
    shift = (mat @ centre - centre)[:3]
    return [*rotations, *shift]


def _resample(volume, affine, pull, ref_img):
    """One volume resampled onto the reference grid through a RAS pull (cubic B-spline)."""
    from scipy.ndimage import map_coordinates

    shape = ref_img.shape[:3]
    ijk = np.indices(shape).reshape(3, -1)
    world = ref_img.affine @ np.vstack([ijk, np.ones((1, ijk.shape[1]))])
    src = np.linalg.inv(affine) @ pull @ world
    # Zero outside the field of view, as in the probe's validated implementation
    # (nitransforms' default), so the means match those C9-auto was judged on.
    return map_coordinates(volume, src[:3], order=3, mode='constant', cval=0.0).reshape(shape)


def _magnitude(pull, ref_img):
    """Translation at the reference's centre (mm) and rotation angle (degrees) of a transform."""
    centre = ref_img.affine @ np.r_[(np.array(ref_img.shape[:3]) - 1) / 2, 1]
    shift = (pull @ centre - centre)[:3]
    # From both the antisymmetric part and the trace: 3dvolreg writes its matrices with few
    # decimals, so for small rotations the diagonal rounds to 1 and the trace alone gives 0.
    r = pull[:3, :3]
    sin = np.linalg.norm([r[2, 1] - r[1, 2], r[0, 2] - r[2, 0], r[1, 0] - r[0, 1]]) / 2
    angle = np.degrees(np.arctan2(sin, (np.trace(r) - 1) / 2))
    return {'translation_mm': float(np.linalg.norm(shift)), 'rotation_deg': float(angle)}
