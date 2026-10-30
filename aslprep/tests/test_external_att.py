"""Tests for single-delay tissue-T1 CBF with an ATT map from a multi-delay run (gh-505)."""

import json

import nibabel as nb
import numpy as np
import pytest

from aslprep import config
from aslprep.tests.tests import mock_config, reset_config
from aslprep.utils import att as att_utils
from aslprep.utils.cbf import ATT_STATUS, quantify_single_delay_tissue_t1

AFFINE = np.diag([3.0, 3.0, 4.0, 1.0])
SHAPE = (10, 10, 6)
TISSUE_T1_MODEL = {
    'Model': {
        'Description': 'The tissue-T1 model.',
        'Parameters': {
            'GrayMatterT1': 1.6,
            'WhiteMatterT1': 0.84,
            'NominalGrayMatterCBF': 60.0,
            'NominalWhiteMatterCBF': 20.0,
            'BloodT1': 1.65,
        },
    },
}


def _itk(path, translation=(0.0, 0.0, 0.0)):
    """Write an ITK affine transform: a translation, in ITK's (LPS) convention."""
    path.write_text(
        '#Insight Transform File V1.0\n'
        '#Transform 0\n'
        'Transform: MatrixOffsetTransformBase_double_3_3\n'
        f'Parameters: 1 0 0 0 1 0 0 0 1 {" ".join(str(t) for t in translation)}\n'
        'FixedParameters: 0 0 0\n'
    )
    return str(path)


def _write(data, path, affine=AFFINE, dtype=np.float32):
    img = nb.Nifti1Image(np.asarray(data).astype(dtype), affine)
    img.set_data_dtype(dtype)
    img.to_filename(path)
    return str(path)


def _derivatives(root, *, model=TISSUE_T1_MODEL, runs=('01',), skip=()):
    """A derivatives dataset with a multi-delay run's ATT products and the anatomy."""
    root.mkdir()
    description = {'Name': 'aslprep', 'BIDSVersion': '1.9.0', 'DatasetType': 'derivative'}
    (root / 'dataset_description.json').write_text(json.dumps(description))
    perf = root / 'sub-01' / 'perf'
    perf.mkdir(parents=True)
    anat = root / 'sub-01' / 'anat'
    anat.mkdir()
    for run in runs:
        prefix = f'sub-01_acq-mpld_run-{run}'
        files = {
            'att': f'{prefix}_att.nii.gz',
            'att_status': f'{prefix}_attstatus.nii.gz',
            'att_ci_lower': f'{prefix}_desc-ciLower_att.nii.gz',
            'att_ci_upper': f'{prefix}_desc-ciUpper_att.nii.gz',
        }
        for key, name in files.items():
            if key not in skip:
                _write(np.ones(SHAPE), perf / name)

        (perf / f'{prefix}_att.json').write_text(json.dumps(model))
        if 'aslref2anat' not in skip:
            _itk(perf / f'{prefix}_from-aslref_to-T1w_mode-image_desc-coreg_xfm.txt')

    if 't1w_preproc' not in skip:
        _write(np.ones((8, 8, 8)), anat / 'sub-01_desc-preproc_T1w.nii.gz', np.eye(4))

    return str(root)


def test_collect_external_att(tmp_path):
    """A native-space ATT map is found with its companions, and duplicates are refused."""
    from aslprep.utils.bids import _get_layout, collect_external_att

    _get_layout.cache_clear()
    derivatives = {'mpld': _derivatives(tmp_path / 'mpld')}
    found = collect_external_att(derivatives, '01')
    assert found['att'].endswith('sub-01_acq-mpld_run-01_att.nii.gz')
    assert found['att_status'].endswith('_attstatus.nii.gz')
    assert found['att_ci_lower'].endswith('_desc-ciLower_att.nii.gz')
    assert found['att_ci_upper'].endswith('_desc-ciUpper_att.nii.gz')
    assert found['aslref2anat'].endswith('_desc-coreg_xfm.txt')
    assert found['t1w_preproc'].endswith('sub-01_desc-preproc_T1w.nii.gz')
    assert found['uri'] == 'bids:mpld:sub-01/perf/sub-01_acq-mpld_run-01_att.nii.gz'
    assert found['metadata']['Model']['Parameters']['GrayMatterT1'] == 1.6
    assert collect_external_att(derivatives, '02') is None
    assert collect_external_att({}, '01') is None

    # Several runs need a filter
    two = {'mpld': _derivatives(tmp_path / 'two', runs=('01', '02'))}
    with pytest.raises(ValueError, match='Found 2 ATT maps'):
        collect_external_att(two, '01')

    found = collect_external_att(two, '01', bids_filters={'att': {'run': '02'}})
    assert found['att'].endswith('run-02_att.nii.gz')
    assert found['att_status'].endswith('run-02_attstatus.nii.gz')

    # Every companion is needed
    for missing in ('att_status', 'aslref2anat', 't1w_preproc'):
        partial = {'mpld': _derivatives(tmp_path / f'no-{missing}', skip=(missing,))}
        with pytest.raises(ValueError, match=f'lacks the {missing}'):
            collect_external_att(partial, '01')


def test_collect_external_att_source_metadata(tmp_path):
    """The multi-delay run's own sidecar is read from the raw dataset, which has its labeling."""
    from bids.layout import BIDSLayout

    from aslprep.utils.bids import _get_layout, collect_external_att

    _get_layout.cache_clear()
    derivatives = {'mpld': _derivatives(tmp_path / 'mpld')}
    raw = tmp_path / 'raw'
    perf = raw / 'sub-01' / 'perf'
    perf.mkdir(parents=True)
    description = {'Name': 'raw', 'BIDSVersion': '1.9.0'}
    (raw / 'dataset_description.json').write_text(json.dumps(description))
    _write(np.zeros((2, 2, 2, 2)), perf / 'sub-01_acq-mpld_run-01_asl.nii.gz', np.eye(4))
    sidecar = {'ArterialSpinLabelingType': 'PCASL', 'LabelingDistance': 90}
    (perf / 'sub-01_acq-mpld_run-01_asl.json').write_text(json.dumps(sidecar))
    layout = BIDSLayout(raw, validate=False)
    found = collect_external_att(derivatives, '01', layout=layout)
    assert found['source_metadata']['ArterialSpinLabelingType'] == 'PCASL'
    assert collect_external_att(derivatives, '01')['source_metadata'] is None


def test_check_external_att_metadata():
    """A consensus-model ATT map is refused, and labeling that differs warns."""
    with pytest.raises(ValueError, match='not estimated with the tissue-T1 kinetic model'):
        att_utils.check_external_att_metadata(
            {'Model': {'Parameters': {'BloodT1': 1.65}}},
            {},
            {},
            'att',
            'asl',
        )

    source = {'ArterialSpinLabelingType': 'PCASL', 'LabelingDistance': 90}
    assert not att_utils.check_external_att_metadata(TISSUE_T1_MODEL, source, source, 'att', 'asl')
    warnings = att_utils.check_external_att_metadata(
        TISSUE_T1_MODEL,
        source,
        {'ArterialSpinLabelingType': 'PASL', 'LabelingDistance': 90},
        'att',
        'asl',
    )
    assert len(warnings) == 1
    assert 'ArterialSpinLabelingType PCASL' in warnings[0]

    # Without the multi-delay run's sidecar, the labeling cannot be compared
    warnings = att_utils.check_external_att_metadata(TISSUE_T1_MODEL, None, source, 'att', 'asl')
    assert 'could not be compared' in warnings[0]


def test_target_to_source_voxels(tmp_path):
    """Target voxels go to the anatomy through their run's transform, then to the source run's."""
    target_affine = np.diag([2.0, 2.0, 2.0, 1.0])
    target_affine[:3, 3] = [-10, -10, -10]
    source_affine = AFFINE.copy()
    source_affine[:3, 3] = [-12, -12, -12]
    identity = _itk(tmp_path / 'identity.txt')
    shifted = _itk(tmp_path / 'shifted.txt', (3.0, -6.0, 4.0))

    grid = np.indices((4, 4, 4)).reshape(3, -1)
    points = np.vstack([grid, np.ones(grid.shape[1])])
    # With identity transforms, only the grids differ
    coords = att_utils.target_to_source_voxels(
        target_affine,
        (4, 4, 4),
        source_affine,
        identity,
        identity,
    )
    expected = (np.linalg.inv(source_affine) @ target_affine @ points)[:3]
    np.testing.assert_allclose(coords, expected, atol=1e-9)

    # The same transform in both runs cancels
    coords = att_utils.target_to_source_voxels(
        target_affine,
        (4, 4, 4),
        source_affine,
        shifted,
        shifted,
    )
    np.testing.assert_allclose(coords, expected, atol=1e-9)

    # The source run's transform takes anatomical points to its reference. An ITK transform
    # maps points of the fixed (anatomical) image to the moving (ASL reference) image, in LPS,
    # so a translation of (3, -6, 4) in the file is (-3, 6, 4) in RAS.
    ras_shift = np.array([-3.0, 6.0, 4.0, 0.0])[:, None]
    coords = att_utils.target_to_source_voxels(
        target_affine,
        (4, 4, 4),
        source_affine,
        identity,
        shifted,
    )
    expected = (np.linalg.inv(source_affine) @ (target_affine @ points + ras_shift))[:3]
    np.testing.assert_allclose(coords, expected, atol=1e-9)

    # And the target run's transform the other way: its inverse takes target points to the
    # anatomical image
    coords = att_utils.target_to_source_voxels(
        target_affine,
        (4, 4, 4),
        source_affine,
        shifted,
        identity,
    )
    expected = (np.linalg.inv(source_affine) @ (target_affine @ points - ras_shift))[:3]
    np.testing.assert_allclose(coords, expected, atol=1e-9)


def test_weighted_median():
    """The weighted median, which ignores voxels without weight or value."""
    assert att_utils.weighted_median([1.0, 2.0, 3.0], [1, 1, 1]) == 2.0
    assert att_utils.weighted_median([1.0, 2.0, 3.0], [0, 0, 5]) == 3.0
    assert att_utils.weighted_median([1.0, np.nan, 3.0], [1, 10, 1]) == 1.0
    assert np.isnan(att_utils.weighted_median([1.0], [0.0]))


def _source_maps():
    """ATT located in the first half of the image, with an interval elsewhere."""
    x = np.indices(SHAPE)[0].astype(float)
    att = 0.8 + 0.05 * x
    status = np.full(SHAPE, ATT_STATUS['located'], dtype=np.uint8)
    status[5:] = ATT_STATUS['wide_range']
    status[9] = ATT_STATUS['not_fitted']
    lower = att - 0.1
    upper = att + 0.1
    # Where ATT was not located, the interval excludes the median, which it then clips to
    lower[5:9] = 1.6
    upper[5:9] = 2.0
    return att, status, lower, upper


def test_resample_external_att_fill():
    """Located ATT is kept, and the rest is filled and clipped, or left without ATT."""
    att, status, lower, upper = _source_maps()
    coords = np.indices(SHAPE).reshape(3, -1).astype(float)
    n = coords.shape[1]
    gm_pv = np.full(n, 0.6)
    wm_pv = np.full(n, 0.3)
    mask = np.ones(n, bool)
    x = coords[0]
    gm_pv[(x == 6)] = 0.0
    wm_pv[(x == 6)] = 0.0

    target, target_status, summary = att_utils.resample_external_att(
        att,
        status,
        lower,
        upper,
        coords,
        gm_pv,
        wm_pv,
        mask,
        (60.0, 20.0),
        fill='median',
    )
    located = x < 5
    assert (target_status[located] == 1).all()
    np.testing.assert_allclose(target[located], att.reshape(-1)[located])
    # The located ATTs (0.8 to 1.0 s) have a median of 0.9 s, which the interval clips to 1.6 s
    filled = (x >= 5) & (x < 9) & (x != 6)
    assert (target_status[filled] == 2).all()
    np.testing.assert_allclose(target[filled], 1.6)
    np.testing.assert_allclose(summary['GrayMatterMedianATT'], 0.9)
    # No tissue, so no fill; the last slab was not fitted
    assert (target_status[x == 6] == 5).all()
    assert (target_status[x == 9] == 0).all()
    assert np.isnan(target[x == 9]).all()
    # Every status but 6, which quantification sets
    expected_codes = {str(code) for code in att_utils.EXTERNAL_ATT_STATUS if code != 6}
    assert set(summary['ResamplingStatusFractions']) == expected_codes

    target, target_status, _ = att_utils.resample_external_att(
        att,
        status,
        lower,
        upper,
        coords,
        gm_pv,
        wm_pv,
        mask,
        (60.0, 20.0),
        fill='none',
    )
    assert (target_status[filled] == 4).all()
    assert np.isnan(target[filled]).all()

    # Without white matter among the measured voxels, its median is unknown, so no fill
    no_wm = np.where(located, 0.0, wm_pv)
    target, target_status, summary = att_utils.resample_external_att(
        att,
        status,
        lower,
        upper,
        coords,
        gm_pv,
        no_wm,
        mask,
        (60.0, 20.0),
        fill='median',
    )
    # The summary goes to the sidecars, which must be valid JSON
    assert summary['WhiteMatterMedianATT'] is None
    json.dumps(summary, allow_nan=False)
    assert (target_status[filled] == 7).all()
    assert np.isnan(target[filled]).all()


def test_resample_external_att_normalized():
    """Interpolation uses only located voxels, so unlocated ones do not blend in."""
    att = np.full(SHAPE, 1.0)
    att[5:] = 50.0  # values where ATT was not located, which must not leak in
    status = np.full(SHAPE, ATT_STATUS['located'], dtype=np.uint8)
    status[5:] = ATT_STATUS['no_tissue_signal']
    # Halfway between voxels 4 and 5, and three quarters of the way to 5
    coords = np.array([[4.5, 4.75], [2.0, 2.0], [2.0, 2.0]])
    target, target_status, _ = att_utils.resample_external_att(
        att,
        status,
        att,
        att,
        coords,
        np.full(2, 0.5),
        np.full(2, 0.5),
        np.ones(2, bool),
        (60.0, 20.0),
        fill='median',
    )
    # Half the weight is located, so it counts; a quarter does not, and the voxel is filled
    assert target_status[0] == 1
    np.testing.assert_allclose(target[0], 1.0)
    assert target_status[1] == 3
    np.testing.assert_allclose(target[1], 1.0)


def _truth_deltam(cbf, m0a, att, delay, bolus, is_casl, t1prime, weights, labeleff=0.85):
    """Delta-M from the numerically integrated tissue signal, independent of the closed form."""
    from aslprep.tests.truth_models import tissue_t1_signal

    # The time of the readout since labeling began (PCASL) or the inversion (PASL)
    t = delay + bolus if is_casl else delay
    signal = sum(
        weight * tissue_t1_signal(t, att, bolus, 1.65, t1p, is_casl)
        for t1p, weight in zip(t1prime, weights, strict=True)
    )
    # The signal is per unit 2 * alpha * M0a * f, with f in mL/g/s
    return 2 * labeleff * m0a * cbf / 6000 * signal


@pytest.mark.parametrize(('is_casl', 'bolus', 'pld'), [(True, 1.8, 1.8), (False, 0.8, 2.0)])
def test_quantify_single_delay_recovers_cbf(is_casl, bolus, pld):
    """Noise-free delta-M with a known ATT gives back CBF, while arriving and after arrival.

    The data come from numerical integration of the model (truth_models), not from the closed
    form that the quantification uses, and each voxel has its own slice-timed delay.
    """
    t1prime = np.array([1.3, 0.75])
    weights = np.array([0.7, 0.3])
    # Arrived; arriving; past the end of the bolus (no signal)
    end = pld + (bolus if is_casl else 0.0)
    atts = np.array([0.5, end - 0.3, end + 0.2])
    plds = pld + np.array([[0.0, 0.0], [0.04, 0.04], [0.08, 0.08]])
    true_cbf = 55.0
    m0 = np.full(3, 1000.0)
    partition = 0.9
    model = {'bolus': bolus, 'is_casl': is_casl, 't1prime': t1prime, 'weights': weights}
    deltam = np.array(
        [
            [_truth_deltam(true_cbf, m0[i] / partition, atts[i], d, **model) for d in plds[i]]
            for i in range(3)
        ]
    )
    cbf_ts, has_att, has_signal = quantify_single_delay_tissue_t1(
        deltam,
        m0,
        plds,
        atts,
        labeleff=0.85,
        t1blood=1.65,
        partition_coefficient=partition,
        is_casl=is_casl,
        bolus=bolus,
        t1prime=t1prime,
        tissue_weights=weights,
    )
    assert has_att.all()
    np.testing.assert_array_equal(has_signal, [True, True, False])
    np.testing.assert_allclose(cbf_ts[:2], true_cbf, rtol=1e-8)
    assert (cbf_ts[2] == 0).all()

    # Without an ATT, no CBF
    cbf_ts, has_att, has_signal = quantify_single_delay_tissue_t1(
        deltam,
        m0,
        plds,
        np.array([np.nan, -1.0, 0.5]),
        labeleff=0.85,
        t1blood=1.65,
        partition_coefficient=partition,
        is_casl=is_casl,
        bolus=bolus,
        t1prime=t1prime,
        tissue_weights=weights,
    )
    np.testing.assert_array_equal(has_att, [False, False, True])
    assert (cbf_ts[:2] == 0).all()


def test_quantify_single_delay_matches_consensus():
    """With T1' equal to blood T1 and the label arrived, it is the consensus (Alsop) formula."""
    t1blood, labeleff, tau, pld, partition = 1.65, 0.85, 1.8, 1.8, 0.9
    deltam = np.array([[10.0, 12.0]])
    m0 = np.array([1000.0])
    cbf_ts, _, _ = quantify_single_delay_tissue_t1(
        deltam,
        m0,
        np.full((1, 2), pld),
        np.array([1.0]),
        labeleff=labeleff,
        t1blood=t1blood,
        partition_coefficient=partition,
        is_casl=True,
        bolus=tau,
        t1prime=np.array([t1blood, t1blood]),
        tissue_weights=np.array([1.0, 0.0]),
    )
    alsop = (
        6000
        * partition
        * deltam
        * np.exp(pld / t1blood)
        / (2 * labeleff * t1blood * (1 - np.exp(-tau / t1blood)) * m0[:, None])
    )
    np.testing.assert_allclose(cbf_ts, alsop, rtol=1e-12)


def _single_delay_inputs(tmp_path, is_casl=True):
    from aslprep.utils.calibration import resolve_calibration

    metadata = {
        'ArterialSpinLabelingType': 'PCASL' if is_casl else 'PASL',
        'MagneticFieldStrength': 3,
        'MRAcquisitionType': '3D',
        'PostLabelingDelay': 1.8,
        'M0Type': 'Estimate',
        'M0Estimate': 1000.0,
        'BackgroundSuppression': False,
    }
    if is_casl:
        metadata['LabelingDuration'] = 1.8
    else:
        metadata.update(
            {
                'BolusCutOffFlag': True,
                'BolusCutOffTechnique': 'QUIPSSII',
                'BolusCutOffDelayTime': 0.8,
            }
        )

    calibration = resolve_calibration(metadata, kinetic_model='tissue-t1', external_att=True)
    shape = (4, 4, 2)
    affine = np.diag([3.0, 3.0, 5.0, 1.0])
    att = np.full(shape, 1.2)
    att[0] = 3.9  # past the end of the bolus
    status = np.ones(shape, np.uint8)
    status[1] = 4  # not filled
    files = {
        'deltam': _write(np.full((*shape, 3), 8.0), tmp_path / 'deltam.nii.gz', affine),
        'm0_file': _write(np.full(shape, 1000.0), tmp_path / 'm0.nii.gz', affine),
        'gm_pv': _write(np.full(shape, 0.6), tmp_path / 'gm_pv.nii.gz', affine),
        'wm_pv': _write(np.full(shape, 0.3), tmp_path / 'wm_pv.nii.gz', affine),
        'external_att': _write(att, tmp_path / 'att.nii.gz', affine),
        'external_att_status': _write(status, tmp_path / 'status.nii.gz', affine, np.uint8),
    }
    return metadata, calibration, files


def test_compute_cbf_single_delay_tissue_t1(tmp_path):
    """ComputeCBF quantifies single-delay CBF with the external ATT, and marks no-signal voxels."""
    from aslprep.interfaces.cbf import ComputeCBF

    metadata, calibration, files = _single_delay_inputs(tmp_path)
    out = (
        ComputeCBF(
            cbf_only=False,
            metadata=metadata,
            calibration=calibration,
            write_calibration=True,
            **files,
        )
        .run(cwd=tmp_path)
        .outputs
    )
    cbf = nb.load(out.mean_cbf).get_fdata()
    status = nb.load(out.external_att_status).get_fdata()
    valid = nb.load(out.quant_valid).get_fdata() > 0
    assert (cbf[2:] > 0).all()
    assert (status[2:] == 1).all()
    # ATT past the end of the bolus: no signal, no CBF
    assert (status[0] == 6).all()
    assert (cbf[0] == 0).all()
    # Not filled: no CBF
    assert (status[1] == 4).all()
    assert (cbf[1] == 0).all()
    np.testing.assert_array_equal(valid, cbf > 0)
    assert out.att is None

    files.pop('external_att')
    with pytest.raises(ValueError, match='needs an arterial transit time map'):
        ComputeCBF(cbf_only=False, metadata=metadata, calibration=calibration, **files).run(
            cwd=tmp_path
        )


def test_compute_cbf_single_delay_pasl(tmp_path):
    """PASL with QUIPSS II is quantified, and QUIPSS (no TI1) is refused."""
    from aslprep.interfaces.cbf import ComputeCBF

    metadata, calibration, files = _single_delay_inputs(tmp_path, is_casl=False)
    out = (
        ComputeCBF(cbf_only=False, metadata=metadata, calibration=calibration, **files)
        .run(cwd=tmp_path)
        .outputs
    )
    assert (nb.load(out.mean_cbf).get_fdata()[2:] > 0).all()

    quipss = {**metadata, 'BolusCutOffTechnique': 'QUIPSS'}
    with pytest.raises(ValueError, match='only QUIPSS II and Q2TIPS'):
        ComputeCBF(cbf_only=False, metadata=quipss, calibration=calibration, **files).run(
            cwd=tmp_path
        )


def test_external_att_interface(tmp_path):
    """The interface resamples the map onto the run's grid, after checking the anatomy."""
    from aslprep.interfaces.att import ExternalATT

    att, status, lower, upper = _source_maps()
    anat_data = np.random.default_rng(0).uniform(size=(8, 8, 8))
    anat = _write(anat_data, tmp_path / 't1w.nii.gz', np.eye(4))
    other = _write(np.ones((8, 8, 8)), tmp_path / 'other.nii.gz', np.eye(4))
    files = {
        'att': _write(att, tmp_path / 'att.nii.gz'),
        'att_status': _write(status, tmp_path / 'status.nii.gz', dtype=np.uint8),
        'att_ci_lower': _write(lower, tmp_path / 'lower.nii.gz'),
        'att_ci_upper': _write(upper, tmp_path / 'upper.nii.gz'),
        'source_aslref2anat': _itk(tmp_path / 'source.txt'),
        'aslref2anat': _itk(tmp_path / 'target.txt'),
        'reference': _write(np.ones(SHAPE), tmp_path / 'mask.nii.gz', dtype=np.uint8),
        'gm_pv': _write(np.full(SHAPE, 0.6), tmp_path / 'gm.nii.gz'),
        'wm_pv': _write(np.full(SHAPE, 0.3), tmp_path / 'wm.nii.gz'),
    }
    interface = ExternalATT(
        **files,
        source_t1w=anat,
        t1w_preproc=anat,
        nominal_cbf=[60.0, 20.0],
        source_name='bids:mpld:sub-01/perf/sub-01_att.nii.gz',
    )
    out = interface.run(cwd=tmp_path).outputs
    resampled = nb.load(out.att).get_fdata()
    out_status = nb.load(out.att_status).get_fdata()
    np.testing.assert_allclose(resampled[:5], att[:5], rtol=1e-6)
    assert (out_status[:5] == 1).all()
    assert (out_status[5:9] == 2).all()
    assert out.summary['Source'].startswith('bids:mpld:')
    assert out.summary['Fill'] == 'median'

    interface.inputs.t1w_preproc = other
    with pytest.raises(ValueError, match='not the anatomical image of this run'):
        interface.run(cwd=tmp_path)


def test_init_cbf_wf_external_att(tmp_path):
    """Single-delay tissue-T1 CBF resamples the external ATT and gives it to the fit."""
    import pandas as pd

    from aslprep.workflows.asl.cbf import init_cbf_wf

    asl_file = tmp_path / 'sub-01_asl.nii.gz'
    asl_file.touch()
    metadata, calibration, _ = _single_delay_inputs(tmp_path)
    aslcontext = tmp_path / 'sub-01_aslcontext.tsv'
    pd.DataFrame({'volume_type': ['control', 'label'] * 3}).to_csv(
        aslcontext,
        sep='\t',
        index=False,
    )
    external_att = {
        'att': 'att.nii.gz',
        'att_status': 'attstatus.nii.gz',
        'att_ci_lower': 'lower.nii.gz',
        'att_ci_upper': 'upper.nii.gz',
        'aslref2anat': 'xfm.txt',
        't1w_preproc': 'T1w.nii.gz',
        'metadata': TISSUE_T1_MODEL,
        'uri': 'bids:mpld:sub-01/perf/sub-01_att.nii.gz',
    }
    reset_config()
    with mock_config():
        config.workflow.kinetic_model = 'tissue-t1'
        config.workflow.att_fill = 'none'
        wf = init_cbf_wf(
            name_source=str(asl_file),
            processing_target='control',
            aslcontext=str(aslcontext),
            metadata=metadata,
            calibration=calibration,
            external_att=external_att,
        )
        desc = wf.__desc__
        # Without a calibration, the workflow resolves one that allows the external ATT
        resolved = init_cbf_wf(
            name_source=str(asl_file),
            processing_target='control',
            aslcontext=str(aslcontext),
            metadata=metadata,
            external_att=external_att,
        )

    reset_config()
    assert resolved.get_node('resample_external_att') is not None
    resample = wf.get_node('resample_external_att')
    assert resample.inputs.fill == 'none'
    assert resample.inputs.source_name == external_att['uri']
    to_cbf = wf._graph.get_edge_data(resample, wf.get_node('compute_cbf'))['connect']
    assert ('att', 'external_att') in to_cbf
    assert ('att_status', 'external_att_status') in to_cbf
    to_resample = wf._graph.get_edge_data(wf.get_node('inputnode'), resample)['connect']
    assert ('t1w_preproc', 't1w_preproc') in to_resample
    assert wf.get_node('external_att_files').inputs.att == 'att.nii.gz'
    # The reportlet shows the final status (with status 6) and names the source run
    plot = wf.get_node('external_att_plot')
    assert plot.inputs.source_name == external_att['uri']
    to_plot = wf._graph.get_edge_data(wf.get_node('compute_cbf'), plot)['connect']
    assert ('external_att_status', 'att_status') in to_plot
    to_out = wf._graph.get_edge_data(plot, wf.get_node('outputnode'))['connect']
    assert ('out_file', 'external_att_report') in to_out
    assert 'multi-delay ASL run' in desc
    assert 'given no CBF' in desc
    assert 'apparent T1 of the tissue' in desc
    # The boilerplate is a methods description, so it names no files
    assert 'sub-01' not in desc


def test_external_att_sidecars():
    """The sidecars name the source run and describe the status levels."""
    from aslprep.workflows.asl.base import _external_att_metadata

    summary = {
        'Source': 'bids:mpld:sub-01/perf/sub-01_att.nii.gz',
        'StatusLevels': {'1': 'measured'},
        'Fill': 'median',
    }
    att = _external_att_metadata(summary, 'att')
    assert att['Sources'] == [summary['Source']]
    assert summary['Source'] in att['Description']
    status = _external_att_metadata(summary, 'attstatus')
    assert status['Levels'] == {'1': 'measured'}
