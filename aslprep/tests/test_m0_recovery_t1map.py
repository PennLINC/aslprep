"""Tests for correcting each voxel of an M0 scan for T1 recovery with a T1 map on its grid."""

import json

import nibabel as nb
import numpy as np
import pytest

from aslprep.utils.calibration import (
    basil_calibration_inputs,
    calibration_value,
    resolve_calibration,
)
from aslprep.utils.t1map import MAX_M0_RECOVERY_T1, m0_recovery_factors

SHAPE = (6, 6, 4)
AFFINE = np.diag([3.0, 3.0, 5.0, 1.0])
URI = 'bids:raw:sub-01/anat/sub-01_acq-tr1_T1map.nii.gz'
METADATA = {
    'ArterialSpinLabelingType': 'PCASL',
    'MagneticFieldStrength': 3,
    'M0Type': 'Separate',
    'LabelingDuration': 1.8,
    'PostLabelingDelay': 1.8,
    'BackgroundSuppression': False,
}


def _write(data, path, affine=AFFINE):
    nb.Nifti1Image(np.asarray(data, dtype=np.float32), affine).to_filename(path)
    return str(path)


def test_m0_recovery_factors():
    """Each voxel is corrected with its own T1, long T1s are capped, and missing ones fall back."""
    t1 = np.array([0.8, 1.6, 4.0, 20.0, 0.0, -1.0, np.nan])
    factors, fallback = m0_recovery_factors(t1, 2.0, 1.4)

    expected_t1 = np.array([0.8, 1.6, 4.0, MAX_M0_RECOVERY_T1, 1.4, 1.4, 1.4])
    np.testing.assert_allclose(factors, 1 / (1 - np.exp(-2.0 / expected_t1)))
    np.testing.assert_array_equal(fallback, [False, False, False, False, True, True, True])
    # White matter needs less correction than gray matter, and CSF more
    assert factors[0] < factors[1] < factors[2]


def test_resolve_calibration_voxelwise():
    """A T1 map on the M0 grid replaces the single recovery factor."""
    calibration = resolve_calibration(
        METADATA,
        m0tr=2.0,
        m0scan_metadata={'RepetitionTimePreparation': 2.0},
        m0_t1map=URI,
    )
    assert calibration_value(calibration, 'M0RecoveryCorrection') == (
        'voxelwise saturation recovery'
    )
    assert calibration_value(calibration, 'M0RecoveryT1Map') == URI
    assert calibration['Parameters']['M0RecoveryT1Map']['Source'] == 'measured'
    # The gray matter T1 remains, for the voxels without a T1
    assert calibration_value(calibration, 'TissueT1') == pytest.approx(1.607)
    assert 'M0RecoveryFactor' not in calibration['Parameters']
    # oxford_asl must not correct the corrected image again
    basil_calibration = resolve_calibration(
        METADATA,
        m0tr=2.0,
        m0scan_metadata={'RepetitionTimePreparation': 2.0},
        basil=True,
        m0_t1map=URI,
    )
    assert basil_calibration_inputs(basil_calibration)['m0tr'] == 5.0


@pytest.mark.parametrize(
    ('m0type', 'm0tr', 'correction'),
    [
        # Only a separate M0 scan is corrected before it is resampled
        ('Included', 2.0, 'saturation recovery'),
        # A long TR needs no correction
        ('Separate', 6.0, 'none'),
    ],
)
def test_resolve_calibration_not_voxelwise(m0type, m0tr, correction):
    """The T1 map is not used where the voxelwise correction does not apply."""
    calibration = resolve_calibration(
        {**METADATA, 'M0Type': m0type},
        m0tr=m0tr,
        m0scan_metadata={'RepetitionTimePreparation': m0tr} if m0type == 'Separate' else None,
        m0_t1map=URI,
    )
    assert calibration_value(calibration, 'M0RecoveryCorrection') == correction
    assert 'M0RecoveryT1Map' not in calibration['Parameters']


def test_correct_m0_recovery(tmp_path):
    """Every volume of the M0 scan is corrected, with the map's units, on the M0 grid only."""
    from aslprep.interfaces.t1map import CorrectM0Recovery

    t1_ms = np.full(SHAPE, 1600.0)
    t1_ms[0] = 800.0
    t1_ms[1] = 0.0
    m0 = np.stack([np.full(SHAPE, 100.0), np.full(SHAPE, 200.0)], axis=-1)
    m0_file = _write(m0, tmp_path / 'm0scan.nii.gz')
    t1map = _write(t1_ms, tmp_path / 'T1map.nii.gz')

    result = CorrectM0Recovery(
        in_file=m0_file,
        t1map=t1map,
        metadata={'Units': 'ms'},
        repetition_time=2.0,
        fallback_t1=1.4,
    ).run(cwd=str(tmp_path))
    corrected = nb.load(result.outputs.out_file).get_fdata()

    assert corrected.shape == m0.shape
    for volume, value in ((0, 100.0), (1, 200.0)):
        assert corrected[2, 0, 0, volume] == pytest.approx(value / (1 - np.exp(-2.0 / 1.6)))
        assert corrected[0, 0, 0, volume] == pytest.approx(value / (1 - np.exp(-2.0 / 0.8)))
        assert corrected[1, 0, 0, volume] == pytest.approx(value / (1 - np.exp(-2.0 / 1.4)))

    off_grid = _write(t1_ms, tmp_path / 'off_grid.nii.gz', affine=np.eye(4))
    with pytest.raises(ValueError, match='not on the grid'):
        CorrectM0Recovery(
            in_file=m0_file,
            t1map=off_grid,
            repetition_time=2.0,
            fallback_t1=1.4,
        ).run(cwd=str(tmp_path))


def _bids_dataset(root, name='raw'):
    root.mkdir(parents=True)
    description = {'Name': name, 'BIDSVersion': '1.10.0'}
    if name != 'raw':
        description.update(DatasetType='derivative', GeneratedBy=[{'Name': name}])
    (root / 'dataset_description.json').write_text(json.dumps(description))
    session = root / 'sub-01' / 'ses-1'
    (session / 'anat').mkdir(parents=True)
    (session / 'perf').mkdir(parents=True)
    return session


def test_collect_m0scan_t1map(tmp_path):
    """Only maps on the M0 grid are found, derivatives first, and duplicates are refused."""
    from bids.layout import BIDSLayout

    from aslprep.utils.bids import collect_m0scan_t1map

    raw = _bids_dataset(tmp_path / 'raw')
    m0scan = _write(np.ones(SHAPE + (2,)), raw / 'perf' / 'sub-01_ses-1_m0scan.nii.gz')
    # A high-resolution T1 map is not on the M0 grid
    _write(np.ones((12, 12, 8)), raw / 'anat' / 'sub-01_ses-1_acq-mp2rage_T1map.nii.gz')
    layout = BIDSLayout(str(tmp_path / 'raw'), validate=False)
    assert collect_m0scan_t1map(layout, {}, m0scan) is None

    _write(np.ones(SHAPE), raw / 'anat' / 'sub-01_ses-1_acq-tr1_T1map.nii.gz')
    (raw / 'anat' / 'sub-01_ses-1_acq-tr1_T1map.json').write_text(json.dumps({'Units': 'ms'}))
    layout = BIDSLayout(str(tmp_path / 'raw'), validate=False)
    found = collect_m0scan_t1map(layout, {}, m0scan)
    assert found['t1map'].endswith('sub-01_ses-1_acq-tr1_T1map.nii.gz')
    assert found['metadata']['Units'] == 'ms'
    assert found['uri'] == 'bids:raw:sub-01/ses-1/anat/sub-01_ses-1_acq-tr1_T1map.nii.gz'

    # A derivative map, here in perf, is preferred
    deriv = _bids_dataset(tmp_path / 'qmri', name='qmri')
    _write(np.ones(SHAPE), deriv / 'perf' / 'sub-01_ses-1_desc-tdp_T1map.nii.gz')
    found = collect_m0scan_t1map(layout, {'qmri': str(tmp_path / 'qmri')}, m0scan)
    assert found['uri'] == 'bids:qmri:sub-01/ses-1/perf/sub-01_ses-1_desc-tdp_T1map.nii.gz'

    # Two maps on the grid in one dataset are refused, unless a filter selects one
    _write(np.ones(SHAPE), raw / 'anat' / 'sub-01_ses-1_acq-tr2_T1map.nii.gz')
    layout = BIDSLayout(str(tmp_path / 'raw'), validate=False)
    with pytest.raises(ValueError, match='Found 2 T1 maps on the grid'):
        collect_m0scan_t1map(layout, {}, m0scan)

    found = collect_m0scan_t1map(
        layout,
        {},
        m0scan,
        bids_filters={'m0scan_t1map': {'acquisition': 'tr2'}},
    )
    assert found['t1map'].endswith('sub-01_ses-1_acq-tr2_T1map.nii.gz')


def test_collect_m0scan_t1map_sessionless(tmp_path):
    """An M0 scan without a session only uses a map without one, not another session's."""
    from bids.layout import BIDSLayout

    from aslprep.utils.bids import collect_m0scan_t1map

    session = _bids_dataset(tmp_path / 'raw')
    _write(np.ones(SHAPE), session / 'anat' / 'sub-01_ses-1_T1map.nii.gz')
    subject = tmp_path / 'raw' / 'sub-01'
    (subject / 'perf').mkdir()
    (subject / 'anat').mkdir()
    m0scan = _write(np.ones(SHAPE), subject / 'perf' / 'sub-01_m0scan.nii.gz')
    layout = BIDSLayout(str(tmp_path / 'raw'), validate=False)
    assert collect_m0scan_t1map(layout, {}, m0scan) is None

    _write(np.ones(SHAPE), subject / 'anat' / 'sub-01_T1map.nii.gz')
    layout = BIDSLayout(str(tmp_path / 'raw'), validate=False)
    found = collect_m0scan_t1map(layout, {}, m0scan)
    assert found['uri'] == 'bids:raw:sub-01/anat/sub-01_T1map.nii.gz'


@pytest.mark.parametrize(
    ('ignore', 'm0tr', 'expected'),
    [
        ([], 2.0, True),
        (['t1map'], 2.0, False),
        ([], 6.0, False),
    ],
)
def test_init_asl_wf_m0scan_t1map(tmp_path, monkeypatch, ignore, m0tr, expected):
    """The M0 scan is corrected before it is resampled, only when the correction applies."""
    from aslprep.tests.tests import mock_config, reset_config
    from aslprep.utils import bids
    from aslprep.workflows.asl import base

    collect_run_data = base.collect_run_data
    m0scan = _write(np.ones(SHAPE), tmp_path / 'sub-01_m0scan.nii.gz')
    t1map = _write(np.full(SHAPE, 1.4), tmp_path / 'sub-01_T1map.nii.gz')

    def _separate_m0(layout, asl_file):
        run_data = collect_run_data(layout, asl_file)
        run_data['asl_metadata'] = {**run_data['asl_metadata'], 'M0Type': 'Separate'}
        run_data['m0scan'] = m0scan
        run_data['m0scan_metadata'] = {'RepetitionTimePreparation': m0tr}
        return run_data

    monkeypatch.setattr(base, 'collect_run_data', _separate_m0)
    monkeypatch.setattr(
        bids,
        'collect_m0scan_t1map',
        lambda *args, **kwargs: {'t1map': t1map, 'metadata': {'Units': 's'}, 'uri': URI},
    )

    reset_config()
    with mock_config():
        from aslprep import config

        config.workflow.ignore = ignore
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        wf = base.init_asl_wf(asl_file=str(asl_file))

    reset_config()
    node = wf.get_node('correct_m0_recovery')
    assert (node is not None) == expected
    description = wf.get_node('cbf_wf').__desc__
    assert ('each voxel of the calibration volumes' in description) == expected
    if expected:
        assert node.inputs.repetition_time == 2.0
        assert node.inputs.t1map == t1map
        edges = wf._graph.get_edge_data(node, wf.get_node('asl_native_wf'))
        assert ('out_file', 'inputnode.m0scan') in edges['connect']
