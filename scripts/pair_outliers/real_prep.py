"""Prepare real multi-delay ASL runs from aslprep-zoo for the pair-exclusion evaluation.

For each run:

1. The control and label volumes are motion-corrected together to their mean with FSL
   ``mcflirt`` (``-cost mutualinfo``, as ASLPrep's head-motion correction uses), from the
   ``pennlinc/aslprep:unstable`` image. M0 volumes are left out of the correction.
2. M0 is the mean of the run's M0 volumes (or its separate M0 scan), smoothed with a 5 mm
   FWHM Gaussian, as ASLPrep's ``ExtractCBF`` does. The brain mask comes from ``bet``
   (``-f 0.5``) on the unsmoothed M0.
3. Each pair's delays come from ``PostLabelingDelay``, shifted by ``SliceTiming`` for 2D runs.

The aslprep-zoo is read in place and never written to.

Usage::

    python real_prep.py OUT_DIR

Writes ``OUT_DIR/<run>.npz`` with the corrected control and label volumes, M0, mask and
metadata, and ``OUT_DIR/<run>_motion.par`` with mcflirt's motion parameters.
"""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
from scipy import ndimage

ZOO = Path('/mnt/c/Users/tsalo/Documents/datasets/aslprep-zoo')
IMAGE = 'pennlinc/aslprep:unstable'

#: Runs: name -> (path of the ASL file relative to the zoo, separate M0 file or None)
RUNS = {
    '2dmb_sub-1': ('2d_mb_pcasl/sub-1/perf/sub-1_asl.nii.gz', None),
    'tarricone_sub-Sub002': (
        'tarricone/dset/sub-Sub002/perf/sub-Sub002_asl.nii.gz',
        'tarricone/dset/sub-Sub002/perf/sub-Sub002_m0scan.nii.gz',
    ),
    'tarricone_sub-Sub004': (
        'tarricone/dset/sub-Sub004/perf/sub-Sub004_asl.nii.gz',
        'tarricone/dset/sub-Sub004/perf/sub-Sub004_m0scan.nii.gz',
    ),
    'tarricone_sub-Sub006': (
        'tarricone/dset/sub-Sub006/perf/sub-Sub006_asl.nii.gz',
        'tarricone/dset/sub-Sub006/perf/sub-Sub006_m0scan.nii.gz',
    ),
    'lsa_sub-729': ('lsa-aslprep-debug/dset/sub-729/ses-1/perf/sub-729_ses-1_asl.nii.gz', None),
    'lsa_sub-747': ('lsa-aslprep-debug/dset/sub-747/ses-1/perf/sub-747_ses-1_asl.nii.gz', None),
    'lsa_sub-761': ('lsa-aslprep-debug/dset/sub-761/ses-1/perf/sub-761_ses-1_asl.nii.gz', None),
    'lsa_sub-762': ('lsa-aslprep-debug/dset/sub-762/ses-1/perf/sub-762_ses-1_asl.nii.gz', None),
    'grotzinger_sub-01': (
        'grotzinger_menopausal_asl/dset/sub-01/ses-01/perf/sub-01_ses-01_asl.nii.gz',
        'grotzinger_menopausal_asl/dset/sub-01/ses-01/perf/sub-01_ses-01_m0scan.nii.gz',
    ),
}


def _docker(work, *args):
    """Run an FSL command in the ASLPrep image, with ``work`` mounted at /work."""
    cmd = ['docker', 'run', '--rm', '-u', '1000:1000', '-v', f'{work}:/work', '-w', '/work']
    cmd += ['--entrypoint', args[0], IMAGE, *args[1:]]
    subprocess.run(cmd, check=True, capture_output=True)


def _smooth(data, affine, fwhm=5.0):
    voxel = np.sqrt((affine[:3, :3] ** 2).sum(axis=0))
    sigma = fwhm / np.sqrt(8 * np.log(2)) / voxel
    return ndimage.gaussian_filter(data, sigma)


def prepare(name, asl_rel, m0_rel, out_dir):
    asl_file = ZOO / asl_rel
    metadata = json.loads(
        asl_file.with_name(asl_file.name.replace('.nii.gz', '.json')).read_text()
    )
    context = pd.read_table(str(asl_file).replace('_asl.nii.gz', '_aslcontext.tsv'))
    types = context['volume_type'].to_numpy()
    img = nb.load(asl_file)
    data = img.get_fdata(dtype=np.float32)
    control_idx = np.flatnonzero(types == 'control')
    label_idx = np.flatnonzero(types == 'label')
    if control_idx.size != label_idx.size:
        raise ValueError(f'{name}: {control_idx.size} control and {label_idx.size} label volumes')

    # Pair the k-th control with the k-th label volume.
    pair_idx = np.column_stack([control_idx, label_idx])
    plds_all = np.broadcast_to(np.atleast_1d(metadata['PostLabelingDelay']), (types.size,))
    plds = plds_all[control_idx].astype(float)
    if not np.allclose(plds, plds_all[label_idx]):
        raise ValueError(f'{name}: control and label delays differ')

    if m0_rel is not None:
        m0 = nb.load(ZOO / m0_rel).get_fdata(dtype=np.float32)
        m0 = m0.mean(axis=3) if m0.ndim == 4 else m0
    else:
        m0 = data[..., types == 'm0scan'].mean(axis=3)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        order = pair_idx.ravel()
        nb.Nifti1Image(data[..., order], img.affine, img.header).to_filename(tmp / 'cl.nii.gz')
        _docker(tmp, 'mcflirt', '-in', 'cl.nii.gz', '-out', 'cl_mc', '-meanvol',
                '-cost', 'mutualinfo', '-plots')  # fmt: skip
        corrected = nb.load(tmp / 'cl_mc.nii.gz').get_fdata(dtype=np.float32)
        motion = np.loadtxt(tmp / 'cl_mc.par')
        m0_smooth = _smooth(m0, img.affine)
        nb.Nifti1Image(m0, img.affine).to_filename(tmp / 'm0.nii.gz')
        _docker(tmp, 'bet', 'm0.nii.gz', 'm0_brain', '-m', '-f', '0.5')
        mask = nb.load(tmp / 'm0_brain_mask.nii.gz').get_fdata() > 0

    control = corrected[..., 0::2]
    label = corrected[..., 1::2]
    slice_time = np.zeros(data.shape[2])
    if metadata.get('MRAcquisitionType') == '2D' and 'SliceTiming' in metadata:
        slice_time = np.asarray(metadata['SliceTiming'], dtype=float)

    tau = np.broadcast_to(np.atleast_1d(metadata.get('LabelingDuration', np.nan)), (types.size,))
    np.savez_compressed(
        out_dir / f'{name}.npz',
        control=control,
        label=label,
        m0=m0_smooth.astype(np.float32),
        mask=mask,
        affine=img.affine,
        plds=plds,
        tau=tau[control_idx].astype(float),
        slice_time=slice_time,
        readout=metadata.get('MRAcquisitionType', '3D'),
        multiband=int(metadata.get('MultibandAccelerationFactor', 1)),
        background_suppression=bool(metadata.get('BackgroundSuppression', False)),
    )
    np.savetxt(out_dir / f'{name}_motion.par', motion)
    print(
        name, control.shape, 'mask', int(mask.sum()), 'plds', np.unique(plds).tolist(), flush=True
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('out_dir', type=Path)
    parser.add_argument('--runs', nargs='+', default=list(RUNS))
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name in args.runs:
        if (args.out_dir / f'{name}.npz').exists():
            continue
        prepare(name, *RUNS[name], args.out_dir)


if __name__ == '__main__':
    main()
