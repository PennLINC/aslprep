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
    SimpleInterface,
    TraitedSpec,
    traits,
)

#: Labels at shorter delays than these (seconds) still hold arterial signal, so their mean
#: differs from the control mean in structure and biases an offset estimated from them.
ARTERIAL_DELAY = {'PASL': 0.8, 'default': 0.5}


def plan_label_offset_groups(volume_types, plds, asl_type):
    """Group label volumes by delay, and say which delays borrow another delay's estimate.

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
        their own label and control volumes.
    borrowed : dict
        Delay -> (label indices, source delay): delays whose labels still hold arterial signal
        take the estimate of the estimated delay nearest in acquisition order.
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
    early = {d for d in groups if d < cut}
    if early >= set(groups):
        early = set()  # every delay is early: each is estimated from its own volumes
    estimated = {d: g for d, g in groups.items() if d not in early and g[0] and g[1]}
    borrowed = {}
    # Only the early delays borrow; a later delay without controls keeps its transforms.
    for d in sorted(early):
        labels = groups[d][0]
        if not estimated or not labels:
            continue
        at = float(np.mean(labels))
        source = min(estimated, key=lambda e: abs(float(np.mean(estimated[e][0])) - at))
        borrowed[d] = (labels, source)
    return estimated, borrowed


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
    controls, which shows up at the edges of the brain after subtraction. For each delay, the
    mean of the corrected label volumes is registered rigidly to the mean of the corrected control
    volumes (AFNI ``3dvolreg``; the two images differ only by the perfusion signal), and the
    correction is composed into the transforms of every label volume of that delay. Delays whose
    labels still hold arterial signal take the estimate of the nearest other delay. Control,
    M0 and derived volumes keep their transforms. Series without both labels and controls pass
    through unchanged.
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
        estimated, borrowed = plan_label_offset_groups(volume_types, plds, asl_type)
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
                nb.Nifti1Image(np.mean(moved, axis=0), ref_img.affine, ref_img.header).to_filename(
                    means[name]
                )
            aff = os.path.join(runtime.cwd, f'offset_{d:g}.aff12.1D')
            prefix = os.path.join(runtime.cwd, f'offset_{d:g}.nii.gz')
            cmd = ['3dvolreg', '-base', means['control'], '-prefix', prefix]
            cmd += ['-1Dmatrix_save', aff, means['label']]
            subprocess.run(cmd, check=True, capture_output=True)
            # Pull from a control-mean point to the label-mean point (AFNI deobliques implicitly;
            # nitransforms reads its matrices in the plumb frames of the two headers).
            pull = AFNILinearTransformArray.from_filename(aff).to_ras(
                moving=nb.load(means['label']), reference=nb.load(means['control'])
            )
            corrections[d] = np.asarray(pull).reshape(-1, 4, 4)[0]
            record['delays'][f'{d:g}'] = {
                'n_labels': len(labels),
                'n_controls': len(controls),
                **_magnitude(corrections[d], ref_img),
            }
        for d, (labels, source) in borrowed.items():
            corrections[d] = corrections[source]
            record['delays'][f'{d:g}'] = {'n_labels': len(labels), 'borrowed_from': f'{source:g}'}
        groups = {d: g[0] for d, g in estimated.items()} | {d: b[0] for d, b in borrowed.items()}
        for d, labels in groups.items():
            for i in labels:
                # The corrected label shows at y what the label mean shows at P y.
                new[i] = pulls[i] @ corrections[d]
        nt.linear.LinearTransformsMapping(new, reference=ref_img).to_filename(out_file, fmt='itk')
        Path(out_report).write_text(json.dumps(record, indent=2) + '\n')
        self._results['out_file'] = out_file
        self._results['out_report'] = out_report
        return runtime


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
