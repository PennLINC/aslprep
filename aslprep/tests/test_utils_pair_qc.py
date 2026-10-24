"""Tests for observation-level (pair) quality control."""

import json

import nibabel as nb
import numpy as np
import pandas as pd
import pytest

from aslprep.interfaces.confounds import ComputePairQC
from aslprep.utils import pair_qc as pq


def _series(n_obs=20, n_vox=2000, seed=0):
    """Gaussian delta-M around a smooth positive map, voxels x observations."""
    rng = np.random.default_rng(seed)
    signal = 1.0 + 0.5 * np.sin(np.linspace(0, 6, n_vox))
    return signal[:, None] + rng.normal(0, 1.0, (n_vox, n_obs))


def test_clean_series_has_no_outliers():
    """Gaussian noise alone leaves few voxels beyond 3.5 robust SD and flags nothing.

    With 20 observations the voxelwise MAD is itself noisy, so about 0.5% of voxels pass
    3.5 robust SD (0.05% for exact parameters), well below the 5% flag.
    """
    dm = _series()
    out = pq.observation_measures(dm, np.zeros(dm.shape[1], int))
    assert not out['outlier'].any()
    assert out['deviant_fraction'].max() < 0.02
    # About 13% of Gaussian values lie beyond 1.5 SD
    assert out['deviant_fraction_z1p5'].median() == pytest.approx(0.134, abs=0.03)


def test_corrupted_pairs_are_flagged_despite_contamination():
    """Median/MAD keeps flagging when several of ten pairs are corrupted."""
    dm = _series(n_obs=10)
    bad = [2, 5, 7]
    dm[:600, bad] += 8.0  # a third of the brain, in three of ten pairs
    out = pq.observation_measures(dm, np.zeros(10, int))
    assert out.loc[bad, 'outlier'].all()
    assert not out.drop(index=bad)['outlier'].any()


def test_global_shift_shows_in_mean_deltam_z():
    """A uniform offset (labeling failure) moves the brain-mean z-score."""
    dm = _series(n_obs=12)
    dm[:, 4] -= 1.0
    out = pq.observation_measures(dm, np.zeros(12, int))
    assert out.loc[4, 'mean_deltam_z'] < -5
    assert np.abs(out.drop(index=4)['mean_deltam_z']).max() < 4


def test_groups_are_scored_separately_and_small_groups_are_skipped():
    """Each delay is its own reference; delays with too few observations get NaN."""
    dm = np.concatenate([_series(8, seed=1), 5 * _series(8, seed=2), _series(2, seed=3)], axis=1)
    groups = pq.delay_groups([1.0] * 8 + [2.0] * 8 + [3.0] * 2)
    out = pq.observation_measures(dm, groups)
    # The second delay's larger signal is not an outlier against the first delay
    assert not out['outlier'].any()
    assert out['n_same_delay'].tolist() == [8] * 16 + [2] * 2
    assert out.loc[16:, 'deviant_fraction'].isna().all()
    assert pq.MIN_OBSERVATIONS > 2


def test_delay_groups_use_labeling_duration():
    groups = pq.delay_groups([1.5, 1.5, 1.5], [1.8, 1.8, 0.7])
    assert groups[0] == groups[1] != groups[2]


def test_slice_measure_finds_slice_confined_artifact():
    """An artifact in one slice group raises the slice measure more than the whole-brain one."""
    n_slices, per_slice = 50, 100  # one slice is 2% of the brain, under the 5% flag
    dm = _series(n_obs=16, n_vox=n_slices * per_slice)
    slices = np.repeat(np.arange(n_slices), per_slice)
    dm[slices == 7, 3] += 8.0
    out = pq.observation_measures(dm, np.zeros(16, int), slices)
    assert out.loc[3, 'max_slice_deviant_fraction'] > 0.9
    assert out.loc[3, 'deviant_fraction'] < pq.FRACTION_OUTLIER
    assert not out.loc[3, 'outlier']


def test_slice_groups_from_slice_timing():
    """Multiband slices that share a time form one group; a mismatch falls back to slices."""
    np.testing.assert_array_equal(pq.slice_groups([0, 0.1, 0, 0.1], 4), [0, 1, 0, 1])
    np.testing.assert_array_equal(pq.slice_groups(None, 3), [0, 1, 2])
    np.testing.assert_array_equal(pq.slice_groups([0, 0.1], 3), [0, 1, 2])


def test_observation_volumes_pair_by_rank():
    """The k-th control pairs with the k-th label, as in ExtractCBF, whatever the order."""
    kind, vols = pq.observation_volumes(['m0scan', 'control', 'control', 'label', 'label'])
    assert kind == 'control-label'
    assert vols == [(1, 3), (2, 4)]
    kind, vols = pq.observation_volumes(['m0scan', 'deltam', 'deltam'])
    assert kind == 'deltam'
    assert vols == [(1,), (2,)]
    with pytest.raises(ValueError, match='No control-label'):
        pq.observation_volumes(['m0scan'])


def test_input_checks():
    with pytest.raises(ValueError, match='delay groups'):
        pq.observation_measures(np.zeros((10, 4)), np.zeros(3, int))
    with pytest.raises(ValueError, match='slice labels'):
        pq.observation_measures(np.zeros((10, 4)), np.zeros(4, int), np.zeros(9, int))


def _write(arr, path):
    nb.Nifti1Image(arr.astype(np.float32), np.eye(4)).to_filename(path)
    return str(path)


@pytest.mark.parametrize('multi_delay', [False, True])
def test_compute_pair_qc_interface(tmp_path, multi_delay):
    """The interface writes one row per pair with volumes, delays, FD and the measures."""
    shape, n_pairs = (10, 10, 4), 12
    rng = np.random.default_rng(4)
    deltam = 1.0 + rng.normal(0, 1.0, (*shape, n_pairs))
    deltam[:, :, :2, 5] += 10.0  # half of the brain in pair 5
    mask = np.ones(shape)
    mask[0] = 0

    # m0scan first, then all controls before all labels (as some exports store them)
    context = ['m0scan'] + ['control'] * n_pairs + ['label'] * n_pairs
    (tmp_path / 'aslcontext.tsv').write_text('volume_type\n' + '\n'.join(context) + '\n')
    plds = [1.0] * 6 + [2.0] * 6 if multi_delay else 1.8
    metadata = {
        'PostLabelingDelay': plds,
        'LabelingDuration': 1.5,
        'MRAcquisitionType': '2D',
        'SliceTiming': [0.0, 0.05, 0.1, 0.15],
    }
    fd = np.linspace(0.0, 0.5, len(context))
    fd[0] = np.nan
    pd.DataFrame({'framewise_displacement': fd}).to_csv(
        tmp_path / 'confounds.tsv', sep='\t', index=False, na_rep='n/a'
    )

    node = ComputePairQC(
        deltam=_write(deltam, tmp_path / 'deltam.nii.gz'),
        asl_mask=_write(mask, tmp_path / 'mask.nii.gz'),
        aslcontext=str(tmp_path / 'aslcontext.tsv'),
        metadata=metadata,
        confounds_file=str(tmp_path / 'confounds.tsv'),
    )
    result = node.run(cwd=tmp_path)
    table = pd.read_table(result.outputs.qc_file, na_values='n/a')

    assert len(table) == n_pairs
    assert table['observation_type'].unique().tolist() == ['control-label']
    assert table['control_volume'].tolist() == list(range(1, n_pairs + 1))
    assert table['label_volume'].tolist() == list(range(n_pairs + 1, 2 * n_pairs + 1))
    assert table['volume'].isna().all()
    assert table['labeling_duration'].eq(1.5).all()
    np.testing.assert_allclose(table['framewise_displacement'], fd[table['label_volume']])
    assert table['n_same_delay'].eq(6 if multi_delay else n_pairs).all()
    assert table['outlier'].tolist() == [i == 5 for i in range(n_pairs)]
    if multi_delay:
        assert sorted(table['post_labeling_delay'].unique()) == [1.0, 2.0]

    metadata_out = result.outputs.qc_metadata
    assert set(table.columns) <= set(metadata_out)
    json.dumps(metadata_out)  # the datasink writes it as a JSON sidecar


def test_compute_pair_qc_rejects_mismatched_series(tmp_path):
    shape = (4, 4, 2)
    (tmp_path / 'aslcontext.tsv').write_text('volume_type\ncontrol\nlabel\ncontrol\nlabel\n')
    node = ComputePairQC(
        deltam=_write(np.zeros((*shape, 3)), tmp_path / 'deltam.nii.gz'),
        asl_mask=_write(np.ones(shape), tmp_path / 'mask.nii.gz'),
        aslcontext=str(tmp_path / 'aslcontext.tsv'),
        metadata={'PostLabelingDelay': 1.8},
    )
    with pytest.raises(ValueError, match='implies 2 observations'):
        node.run(cwd=tmp_path)
