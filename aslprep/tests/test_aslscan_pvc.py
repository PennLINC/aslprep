"""Local-regression PVC on aslscan's noiseless truth (PVC design spec Section 7.2(a)).

The estimator alone is checked: no preprocessing, and the fixture's exact partial-volume maps
(box overlap on the acquisition grid, which has no readout blur) as regressors.
The signal is the M0-free perfusion signal, q * delta-M, with q from the independent
quantification in :mod:`aslprep.tests.truth_models`.

Tissue amplitudes vary between the slices of a slice-timed acquisition (tissue T1 and arrival
time enter through each slice's delay), so the reference amplitudes are fitted per slice.
The kinetic spread across slices and the calibration error (amplitude / tissue M0 versus the
phantom's perfusion) are reported, not asserted.
"""

import json

import nibabel as nb
import numpy as np
import pytest

from aslprep.tests import aslscan_fixtures as af
from aslprep.tests import truth_models as tm
from aslprep.utils.pvc import local_regression_pvc

#: The flat phantom: constant pure-tissue perfusion (ml/100 g/min).
FLAT_PERFUSION = {'GM': 60.0, 'WM': 20.0}
#: Local LR must match the per-slice amplitudes this closely (relative).
ORACLE_TOL = 0.01
#: Single-slice fits of the noiseless truth must be this close to exact (relative residual).
LINEARITY_TOL = 1e-5


def _fixture(name, data_dir):
    try:
        return af.fixture_dir(name, data_dir)
    except af.FixtureUnavailable as exc:
        pytest.skip(str(exc))


def _signal(fixture):
    """The M0-free signal q * delta-M on the acquisition grid, and the PV truth."""
    recipe = af.load_recipe('pcasl1pld')
    perf = fixture / 'sub-01' / 'perf'
    stem = f'sub-01_acq-{recipe.acq}'
    metadata = json.loads((perf / f'{stem}_asl.json').read_text())
    context = (perf / f'{stem}_aslcontext.tsv').read_text().split()[1:]
    gt = perf / 'ground-truth'

    deltam = nb.load(gt / 'sub-01_desc-deltam_gt.nii.gz').get_fdata()
    rows = [i for i, kind in enumerate(context) if kind in ('label', 'deltam')]
    assert rows, 'no labeled volumes'
    # delta-M is defined on the labeled rows only; the other rows are zero.
    mean_deltam = deltam[..., rows].mean(axis=-1)

    shape = mean_deltam.shape
    pld = tm.pld_map(metadata, shape, context)
    q = tm.cbf_single_delay(
        np.ones(shape), np.ones(shape), metadata, tm.labeling_efficiency(metadata), pld, context
    )
    pv = {t: nb.load(gt / f'sub-01_desc-pv{t}_gt.nii.gz').get_fdata() for t in ('GM', 'WM', 'CSF')}
    simulation = json.loads((gt / 'simulation.json').read_text())
    return q * mean_deltam, pv, metadata, simulation


def _slice_amplitudes(signal, pv, slice_axis=2):
    """Two-tissue least-squares amplitudes per slice, and the relative residual of each fit."""
    amplitudes, residuals = [], []
    for k in range(signal.shape[slice_axis]):
        index = [slice(None)] * 3
        index[slice_axis] = k
        index = tuple(index)
        sel = (pv['GM'][index] + pv['WM'][index]) > 0
        if sel.sum() < 10:
            amplitudes.append((np.nan, np.nan))
            residuals.append(np.nan)
            continue
        design = np.column_stack([pv['GM'][index][sel], pv['WM'][index][sel]])
        target = signal[index][sel]
        beta, *_ = np.linalg.lstsq(design, target, rcond=None)
        amplitudes.append(tuple(beta))
        residuals.append(np.linalg.norm(target - design @ beta) / np.linalg.norm(target))
    return np.array(amplitudes), np.array(residuals)


def test_lr_pvc_oracle_pcasl1pld(data_dir):
    """Local regression recovers the per-slice tissue amplitudes of the flat phantom."""
    fixture = _fixture('pcasl1pld', data_dir)
    signal, pv, metadata, _ = _signal(fixture)
    slice_axis = 'ijk'.index(metadata.get('SliceEncodingDirection', 'k')[0])

    amplitudes, residuals = _slice_amplitudes(signal, pv, slice_axis)
    fitted = np.isfinite(residuals)
    # The premise: within a slice, the truth is a linear mixture of constant tissue amplitudes.
    assert np.nanmax(residuals) < LINEARITY_TOL, np.nanmax(residuals)

    brain = (pv['GM'] + pv['WM'] + pv['CSF']) >= 0.5
    fit = local_regression_pvc(
        signal, pv['GM'], pv['WM'], mask=brain, valid=brain, kernel_size=(5, 5, 1)
    )
    two = fit['status'] == 1
    assert two.sum() > 1000

    k_index = np.indices(signal.shape)[slice_axis]
    for i_tissue, tissue in enumerate(('gm', 'wm')):
        expected = amplitudes[k_index, i_tissue]
        ok = two & np.isfinite(expected)
        ratio = fit[f'beta_{tissue}'][ok] / expected[ok]
        assert np.abs(ratio - 1).max() < ORACLE_TOL, (tissue, np.abs(ratio - 1).max())

    # Reported, not asserted: how much the amplitudes vary between slices (tissue kinetics),
    # and how they compare with the phantom's perfusion.
    gm = amplitudes[fitted, 0]
    wm = amplitudes[fitted, 1]
    print(
        f'GM amplitude spread across slices: {np.ptp(gm) / np.median(gm):.3f}; '
        f'WM: {np.ptp(wm) / np.median(wm):.3f}'
    )
    print(f'GM/WM amplitude ratio (median over slices): {np.median(gm / wm):.3f}')


def test_lr_pvc_oracle_through_plane_kernel_fails(data_dir):
    """Across slices of a slice-timed acquisition, one amplitude per tissue does not fit.

    This is why slice-timed data are limited to single-slice kernels.
    """
    fixture = _fixture('pcasl1pld', data_dir)
    signal, pv, metadata, _ = _signal(fixture)
    if 'SliceTiming' not in metadata:
        pytest.skip('pcasl1pld is not slice-timed')
    sel = (pv['GM'] + pv['WM']) > 0
    design = np.column_stack([pv['GM'][sel], pv['WM'][sel]])
    beta, *_ = np.linalg.lstsq(design, signal[sel], rcond=None)
    residual = np.linalg.norm(signal[sel] - design @ beta) / np.linalg.norm(signal[sel])
    assert residual > 100 * LINEARITY_TOL
