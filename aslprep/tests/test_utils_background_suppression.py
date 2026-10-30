"""Tests for the background-suppression quality-control measures and figure."""

import json

import nibabel as nb
import numpy as np
import pandas as pd
import pytest

from aslprep.utils.confounds import compute_background_suppression_qc, estimate_noise_floor

SHAPE = (40, 40, 24)
ZOOMS = (3.0, 3.0, 3.0)
# The brain is a box, its CSF fraction rising linearly from 0 to 1 along x.
BRAIN = (slice(12, 28), slice(12, 28), slice(6, 18))
SIGMA = 10.0


def _phantom(tmp_path, tissue, csf_signal, deltam, n_pairs=8, background=True, seed=0):
    """Write a magnitude ASL series from signed static signals and complex Gaussian noise.

    ``tissue`` and ``csf_signal`` are the signed static magnetizations at readout,
    and ``deltam`` is the tissue's label-control difference (labeling lowers the signed
    magnetization). The noise is added to the complex signal before taking the magnitude.
    """
    rng = np.random.default_rng(seed)
    affine = np.diag([*ZOOMS, 1.0])
    mask = np.zeros(SHAPE, dtype=bool)
    mask[BRAIN] = True
    csf = np.zeros(SHAPE)
    x = np.arange(BRAIN[0].start, BRAIN[0].stop)
    csf[BRAIN] = ((x - x[0]) / (x[-1] - x[0]))[:, None, None]
    gm = np.where(mask, 1 - csf, 0)

    control = np.where(mask, gm * tissue + csf * csf_signal, 0)
    label = control - np.where(mask, gm * deltam, 0)
    volumes, context = [], []
    for _ in range(n_pairs):
        for name, signed in (('label', label), ('control', control)):
            noise = rng.normal(0, SIGMA, SHAPE) + 1j * rng.normal(0, SIGMA, SHAPE)
            volumes.append(np.abs(signed + noise))
            context.append(name)
    data = np.stack(volumes, axis=-1)
    if not background:
        data[~mask] = 0

    files = {
        'asl_file': tmp_path / 'sub-01_asl.nii.gz',
        'aslcontext': tmp_path / 'sub-01_aslcontext.tsv',
        'mask': tmp_path / 'mask.nii.gz',
        'gm': tmp_path / 'gm.nii.gz',
        'csf': tmp_path / 'csf.nii.gz',
    }
    nb.Nifti1Image(data.astype(np.float32), affine).to_filename(files['asl_file'])
    pd.DataFrame({'volume_type': context}).to_csv(files['aslcontext'], sep='\t', index=False)
    nb.Nifti1Image(mask.astype(np.uint8), affine).to_filename(files['mask'])
    nb.Nifti1Image(gm.astype(np.float32), affine).to_filename(files['gm'])
    nb.Nifti1Image(csf.astype(np.float32), affine).to_filename(files['csf'])
    return {key: str(value) for key, value in files.items()}


def test_noise_floor_of_rayleigh_background():
    """The floor of single-channel magnitude noise is sqrt(2) sigma, even with bright structures."""
    rng = np.random.default_rng(1)
    shape = (30, 30, 30, 6)
    data = np.abs(rng.normal(0, SIGMA, shape) + 1j * rng.normal(0, SIGMA, shape))
    mask = np.zeros(shape[:3], dtype=bool)
    mask[12:18, 12:18, 12:18] = True
    # A bright structure in a corner of the background, as the face or neck would be
    data[:4, :4, :4] += 1000

    floor, background = estimate_noise_floor(data, mask, ZOOMS)
    assert floor == pytest.approx(np.sqrt(2) * SIGMA, rel=0.05)
    assert not background[mask].any()


def test_noise_floor_needs_a_background():
    """A background set to zero by the scanner gives NaN."""
    data = np.zeros((20, 20, 20, 2))
    mask = np.zeros((20, 20, 20), dtype=bool)
    mask[5:15, 5:15, 5:15] = True
    data[mask] = 100

    floor, _ = estimate_noise_floor(data, mask, ZOOMS)
    assert np.isnan(floor)


def test_bs_qc_signal_well_above_floor(tmp_path):
    """Moderate suppression: far from the floor and delta-M positive along the CSF axis."""
    files = _phantom(tmp_path, tissue=300, csf_signal=200, deltam=30)
    qc = compute_background_suppression_qc(**files)
    metrics = qc['metrics']

    assert metrics['bs_noise_floor'] == pytest.approx(np.sqrt(2) * SIGMA, rel=0.05)
    assert metrics['bs_gm_signal_floor_ratio'] > 15
    assert metrics['bs_gm_deltam_retained'] > 0.99
    assert metrics['bs_perc_voxels_near_noise_floor'] == 0
    assert metrics['bs_csf_profile_signal_floor_ratio_min'] > 10
    assert metrics['bs_csf_profile_deltam_min'] > -0.2
    assert len(qc['profile']) == 10
    assert (qc['profile']['n_voxels'] > 0).all()


def test_bs_qc_signal_near_floor(tmp_path):
    """Strong suppression of the tissue: near the floor, delta-M compressed, but not inverted."""
    files = _phantom(tmp_path, tissue=15, csf_signal=15, deltam=5)
    metrics = compute_background_suppression_qc(**files)['metrics']

    # sqrt(15^2 + 2 sigma^2) / (sqrt(2) sigma) = 1.46
    assert metrics['bs_gm_signal_floor_ratio'] == pytest.approx(1.46, rel=0.1)
    assert metrics['bs_gm_deltam_retained'] < 0.8
    assert metrics['bs_perc_voxels_near_noise_floor'] > 90


def test_bs_qc_inverted_tissue(tmp_path):
    """Inverted tissue and positive CSF: the signal dips to the floor and delta-M is negative."""
    # The CSF fraction steps by 1/15, so the voxels closest to the crossing have |signal| 7
    files = _phantom(tmp_path, tissue=-100, csf_signal=100, deltam=20)
    qc = compute_background_suppression_qc(**files)
    metrics = qc['metrics']

    assert metrics['bs_csf_profile_deltam_min'] == pytest.approx(-1)
    assert metrics['bs_csf_profile_signal_floor_ratio_min'] < 1.5
    # The tissue end of the profile is negative, the CSF end is close to zero
    deltam = qc['profile']['deltam_median'].to_numpy()
    assert deltam[0] < 0
    assert abs(deltam[-1]) < abs(deltam[0]) / 4
    # The dip is where the net magnetization crosses zero, at half CSF
    dip = int(np.argmin(qc['profile']['signal_floor_ratio_median']))
    assert dip in (4, 5)


def test_bs_qc_without_background(tmp_path):
    """Without a background there is no floor, but delta-M along CSF is still measured."""
    files = _phantom(tmp_path, tissue=-300, csf_signal=300, deltam=30, background=False)
    metrics = compute_background_suppression_qc(**files)['metrics']

    assert np.isnan(metrics['bs_noise_floor'])
    assert np.isnan(metrics['bs_gm_signal_floor_ratio'])
    assert np.isnan(metrics['bs_perc_voxels_near_noise_floor'])
    assert np.isnan(metrics['bs_csf_profile_signal_floor_ratio_min'])
    assert metrics['bs_csf_profile_deltam_min'] == pytest.approx(-1)


def test_bs_qc_needs_control_and_label(tmp_path):
    """Series without both control and label volumes (e.g., delta-M only) give NaNs."""
    files = _phantom(tmp_path, tissue=300, csf_signal=200, deltam=30)
    n_volumes = nb.load(files['asl_file']).shape[-1]
    pd.DataFrame({'volume_type': ['deltam'] * n_volumes}).to_csv(
        files['aslcontext'],
        sep='\t',
        index=False,
    )
    qc = compute_background_suppression_qc(**files)

    assert all(np.isnan(value) for value in qc['metrics'].values())
    assert qc['profile'].empty


def test_bs_qc_without_tissue_maps(tmp_path):
    """Without tissue maps, only the floor measures are computed."""
    files = _phantom(tmp_path, tissue=300, csf_signal=200, deltam=30)
    metrics = compute_background_suppression_qc(
        asl_file=files['asl_file'],
        aslcontext=files['aslcontext'],
        mask=files['mask'],
    )['metrics']

    assert np.isfinite(metrics['bs_noise_floor'])
    assert np.isfinite(metrics['bs_perc_voxels_near_noise_floor'])
    assert np.isnan(metrics['bs_gm_signal_floor_ratio'])
    assert np.isnan(metrics['bs_csf_profile_deltam_min'])


@pytest.mark.parametrize('with_aslcontext', [True, False])
def test_cbf_qc_interface_bs_columns(tmp_path, with_aslcontext):
    """ComputeCBFQC reports the background-suppression measures, n/a without aslcontext."""
    from aslprep.interfaces.confounds import BS_QC_METADATA, ComputeCBFQC

    files = _phantom(tmp_path, tissue=-300, csf_signal=300, deltam=30)
    mask = nb.load(files['mask'])
    cbf_file = str(tmp_path / 'mean_cbf.nii.gz')
    nb.Nifti1Image(np.full(SHAPE, 50.0, dtype=np.float32), mask.affine).to_filename(cbf_file)
    confounds_file = str(tmp_path / 'confounds.tsv')
    pd.DataFrame({'framewise_displacement': [0.1], 'rmsd': [0.05]}).to_csv(
        confounds_file,
        sep='\t',
        index=False,
    )
    wm_file = str(tmp_path / 'wm.nii.gz')
    nb.Nifti1Image(np.zeros(SHAPE, dtype=np.float32), mask.affine).to_filename(wm_file)

    kwargs = {'aslcontext': files['aslcontext']} if with_aslcontext else {}
    results = ComputeCBFQC(
        name_source=files['asl_file'],
        asl_file=files['asl_file'],
        mean_cbf=cbf_file,
        gm_tpm=files['gm'],
        wm_tpm=wm_file,
        csf_tpm=files['csf'],
        asl_mask=files['mask'],
        t1w_mask=files['mask'],
        confounds_file=confounds_file,
        **kwargs,
    ).run(cwd=str(tmp_path))

    qc_df = pd.read_csv(results.outputs.qc_file, sep='\t')
    with open(results.outputs.qc_metadata) as fobj:
        metadata = json.load(fobj)
    for key in BS_QC_METADATA:
        assert key in qc_df.columns
        assert metadata[key]['LongName']

    if with_aslcontext:
        assert qc_df['bs_csf_profile_deltam_min'].iloc[0] == pytest.approx(-1)
    else:
        assert qc_df[list(BS_QC_METADATA)].isna().all(axis=None)


def test_background_suppression_plot(tmp_path):
    """The figure renders."""
    from aslprep.interfaces.plotting import BackgroundSuppressionPlot

    files = _phantom(tmp_path, tissue=-300, csf_signal=300, deltam=30)
    results = BackgroundSuppressionPlot(
        asl_file=files['asl_file'],
        aslcontext=files['aslcontext'],
        asl_mask=files['mask'],
        gm_tpm=files['gm'],
        csf_tpm=files['csf'],
    ).run(cwd=str(tmp_path))
    assert results.outputs.out_file.endswith('.svg')
    assert '<svg' in open(results.outputs.out_file).read(2000)


@pytest.mark.parametrize('plot', [True, False])
def test_cbf_confounds_wf_bs_report(plot):
    """The background-suppression reportlet is only added when requested."""
    from aslprep.tests.tests import mock_config
    from aslprep.workflows.asl.confounds import init_cbf_confounds_wf

    with mock_config():
        wf = init_cbf_confounds_wf(plot_background_suppression=plot)
    nodes = {name.split('.')[-1] for name in wf.list_node_names()}
    assert ('ds_report_bs' in nodes) is plot
