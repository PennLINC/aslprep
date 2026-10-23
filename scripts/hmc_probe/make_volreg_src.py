#!/usr/bin/env python3
"""Copy an ASLPrep source snapshot with AFNI 3dvolreg in place of its HMC's MCFLIRT steps.

For PennLINC/aslprep#767 (per-type references): its two MCFLIRT steps, the volumes of a type
to their middle volume and then to the type's FLIRT-aligned mean, become 3dvolreg (least
squares), and everything else is the PR's own code: the FLIRT alignment, the combination of the
transforms and, with ``--hmc-method typewise_offset``, its label offset correction. Probe-only:
the interface imports the harness's ``candidates._volreg`` (on the path inside the harness
container) and writes MCFLIRT's outputs (one FSL matrix per volume, a ``.par`` file and the
resampled series), so the workflow downstream is unchanged.

Usage::

    python scripts/hmc_probe/make_volreg_src.py --src <dir with aslprep/> --out <new dir>
"""

import argparse
import shutil
from pathlib import Path

INTERFACE = '''"""Probe-only: AFNI 3dvolreg with MCFLIRT's outputs (scripts/hmc_probe/make_volreg_src.py)."""

import os

import nibabel as nb
import numpy as np
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    OutputMultiObject,
    SimpleInterface,
    TraitedSpec,
)


class _VolregMCFLIRTInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='4D series to register.')
    ref_file = File(exists=True, desc='Reference; the middle volume if not given (as MCFLIRT).')


class _VolregMCFLIRTOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='Registered series.')
    mat_file = OutputMultiObject(File(exists=True), desc='FSL matrices, one per volume.')
    par_file = File(exists=True, desc='Rotations (rad) and translations (mm), one row per volume.')


class VolregMCFLIRT(SimpleInterface):
    """3dvolreg of a series to a reference, with MCFLIRT's outputs."""

    input_spec = _VolregMCFLIRTInputSpec
    output_spec = _VolregMCFLIRTOutputSpec

    def _run_interface(self, runtime):
        from candidates import _volreg

        img = nb.load(self.inputs.in_file)
        data = np.asanyarray(img.dataobj, np.float32)
        data = data[..., None] if data.ndim == 3 else data
        series = os.path.join(runtime.cwd, 'series.nii.gz')
        nb.Nifti1Image(data, img.affine, img.header).to_filename(series)
        ref = self.inputs.ref_file if self.inputs.ref_file else None
        if ref is None:
            ref = os.path.join(runtime.cwd, 'middle.nii.gz')
            nb.Nifti1Image(data[..., data.shape[-1] // 2], img.affine).to_filename(ref)
        base = os.path.join(runtime.cwd, 'volreg')
        mats = _volreg(series, ref, base)
        mat_dir = os.path.join(runtime.cwd, 'volreg.mat')
        os.makedirs(mat_dir, exist_ok=True)
        files, rows = [], []
        for k, m in enumerate(mats):
            f = os.path.join(mat_dir, f'MAT_{k:04d}')
            np.savetxt(f, m, fmt='%.10f')
            files.append(f)
            r = m[:3, :3]
            rows.append([np.arctan2(r[1, 2], r[2, 2]), -np.arcsin(np.clip(r[0, 2], -1, 1)),
                         np.arctan2(r[0, 1], r[0, 0]), *m[:3, 3]])  # fmt: skip
        par = os.path.join(runtime.cwd, 'volreg.par')
        np.savetxt(par, np.array(rows), fmt='%.8f')
        self._results['out_file'] = base + '.nii.gz'
        self._results['mat_file'] = files
        self._results['par_file'] = par
        return runtime
'''

SWAPS = (
    (
        "fsl.MCFLIRT(cost='mutualinfo'),\n        name='type_mcflirt',",
        "VolregMCFLIRT(),\n        name='type_mcflirt',",
    ),
    (
        'fsl.MCFLIRT(save_mats=True, save_plots=True, save_rms=False, '
        "cost='mutualinfo'),\n        name='mcflirt',",
        "VolregMCFLIRT(),\n        name='mcflirt',",
    ),
)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--src', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.out.exists():
        raise SystemExit(f'{args.out} exists')
    shutil.copytree(args.src, args.out)
    (args.out / 'aslprep' / 'interfaces' / 'probe_volreg.py').write_text(INTERFACE)
    hmc = args.out / 'aslprep' / 'workflows' / 'asl' / 'hmc.py'
    text = hmc.read_text()
    for old, new in SWAPS:
        if text.count(old) != 1:
            raise SystemExit(f'expected one occurrence of: {old!r}')
        text = text.replace(old, new)
    # The imports are local to init_asl_hmc_wf (the first match), indented by four spaces.
    old = '    from aslprep.interfaces.utility import ('
    if old not in text:
        raise SystemExit(f'expected the import line: {old!r}')
    text = text.replace(
        old, '    from aslprep.interfaces.probe_volreg import VolregMCFLIRT\n' + old, 1
    )
    hmc.write_text(text)
    print(f'{args.out}: type_mcflirt and mcflirt are 3dvolreg')


if __name__ == '__main__':
    main()
