"""Tests for the aslprep.interfaces.cbf submodule."""

import os

import nibabel as nb
import numpy as np
import pandas as pd
import pytest
from nipype.interfaces.base import isdefined
from traits.api import TraitError

from aslprep.interfaces import cbf
from aslprep.utils.cbf import ATT_STATUS, ROBUST_TUNING_CONSTANTS


def test_scoreandscrubcbf_wavelet_function():
    """The SCRUB weight function input accepts only implemented weight functions."""
    assert cbf.ScoreAndScrubCBF().inputs.wavelet_function == 'huber'
    for wfun in ROBUST_TUNING_CONSTANTS:
        cbf.ScoreAndScrubCBF(wavelet_function=wfun)

    with pytest.raises(TraitError):
        cbf.ScoreAndScrubCBF(wavelet_function='logistics')


def test_scoreandscrubcbf_writes_float32_with_uint8_mask(tmp_path):
    """Test that the SCORE and SCRUB outputs are float32 when the mask is uint8.

    The outputs used to be written with the mask's header, and nibabel writes an image in its
    header's data type, so the CBF maps were stored as scaled uint8 (at most 256 values).
    """
    from aslprep.utils.cbf import _getcbfscore, _scrubcbf

    rng = np.random.default_rng(0)
    shape = (12, 12, 10)
    mask = np.zeros(shape, dtype=np.uint8)
    mask[2:10, 2:10, 2:8] = 1
    gm = rng.uniform(0, 1, shape).astype(np.float32)
    wm = ((1 - gm) * rng.uniform(0, 1, shape)).astype(np.float32)
    csf = np.clip(1 - gm - wm, 0, 1).astype(np.float32)
    cbf_ts = (60 * gm + 20 * wm)[..., None] + rng.normal(0, 8, shape + (20,))
    cbf_ts = (cbf_ts * mask[..., None]).astype(np.float32)

    interface = cbf.ScoreAndScrubCBF(
        cbf_ts=_save_img(cbf_ts, tmp_path, 'cbf_ts.nii.gz'),
        mask=_save_img(mask, tmp_path, 'mask.nii.gz'),
        gm_tpm=_save_img(gm, tmp_path, 'gm.nii.gz'),
        wm_tpm=_save_img(wm, tmp_path, 'wm.nii.gz'),
        csf_tpm=_save_img(csf, tmp_path, 'csf.nii.gz'),
    )
    results = interface.run(cwd=tmp_path)

    # The same computation as the interface, on the saved inputs
    expected_ts, _ = _getcbfscore(
        cbfts=cbf_ts.astype(float),
        wm=wm.astype(float),
        gm=gm.astype(float),
        csf=csf.astype(float),
        mask=mask.astype(float),
        thresh=0.7,
    )
    expected = {
        'cbf_ts_score': expected_ts,
        'mean_cbf_score': np.mean(expected_ts, axis=3),
        'mean_cbf_scrub': _scrubcbf(
            cbf_ts=expected_ts,
            gm=gm.astype(float),
            wm=wm.astype(float),
            csf=csf.astype(float),
            mask=mask.astype(float),
            wfun='huber',
            thresh=0.7,
        ),
    }
    for name, expected_data in expected.items():
        img = nb.load(getattr(results.outputs, name))
        assert img.get_data_dtype() == np.float32, name
        np.testing.assert_allclose(
            img.get_fdata(), expected_data.astype(np.float32), rtol=1e-6, atol=1e-4
        )


def test_computecbf_casl(datasets, tmp_path_factory):
    """Test aslprep.interfaces.cbf.ComputeCBF with (P)CASL."""
    tmpdir = tmp_path_factory.mktemp('test_computecbf_casl')
    aslcontext_file = os.path.join(
        datasets['nonqtab'],
        'sub-01/perf/sub-01_acq-pcaslsinglepldsiemens3d_run-02_aslcontext.tsv',
    )

    aslcontext = pd.read_table(aslcontext_file)
    n_deltam = aslcontext.loc[aslcontext['volume_type'] == 'label'].shape[0]
    n_volumes = aslcontext.shape[0]

    # Simulate ASL data and a brain mask.
    asl_data = np.random.random((30, 30, 30, n_deltam)).astype(np.float32)
    asl_file = _save_img(asl_data, tmpdir, 'asl.nii.gz')
    asl_mask = np.zeros((30, 30, 30), dtype=np.uint8)
    asl_mask[10:20, 10:20, 10:20] = 1
    m0_file = _save_img(asl_mask, tmpdir, 'm0.nii.gz')

    single_pld = 1.5
    plds = np.zeros(n_volumes)
    temp_plds = np.linspace(0.5, 3.5, n_deltam)
    plds[aslcontext['volume_type'] == 'm0scan'] = 0
    plds[aslcontext['volume_type'] == 'label'] = temp_plds
    plds[aslcontext['volume_type'] == 'control'] = temp_plds
    bad_multiple_plds = plds.tolist()
    good_multiple_plds = plds[aslcontext['volume_type'] == 'control']

    BASE_METADATA = {
        'MagneticFieldStrength': 3,
        'LabelingDuration': 1.6,
    }
    ACQ_DICTS = [
        {'MRAcquisitionType': '3D'},
        {
            'MRAcquisitionType': '2D',
            'SliceTiming': list(np.linspace(0.1, 0.5, 30)),
        },
    ]

    for asltype in ['PCASL', 'CASL']:
        for acq_dict in ACQ_DICTS:
            # Scenario 1: PCASL with a single PostLabelingDelay
            # This should produce CBF time series and mean CBF, but no ATT
            metadata = {
                'ArterialSpinLabelingType': asltype,
                'PostLabelingDelay': single_pld,
                **BASE_METADATA,
                **acq_dict,
            }
            interface = cbf.ComputeCBF(
                cbf_only=False,
                deltam=asl_file,
                metadata=metadata,
                m0_scale=1,
                m0_file=m0_file,
            )
            results = interface.run(cwd=tmpdir)
            assert os.path.isfile(results.outputs.cbf_ts)
            cbf_img = nb.load(results.outputs.cbf_ts)
            assert cbf_img.ndim == 4
            assert cbf_img.shape[3] == n_deltam
            assert os.path.isfile(results.outputs.mean_cbf)
            mean_cbf_img = nb.load(results.outputs.mean_cbf)
            assert mean_cbf_img.ndim == 3
            assert results.outputs.att is None

            # Scenario 2: PCASL with one PostLabelingDelay for each volume (bad)
            metadata = {
                'ArterialSpinLabelingType': asltype,
                'PostLabelingDelay': bad_multiple_plds,
                **BASE_METADATA,
                **acq_dict,
            }

            interface = cbf.ComputeCBF(
                cbf_only=False,
                deltam=asl_file,
                metadata=metadata,
                m0_scale=1,
                m0_file=m0_file,
            )
            with pytest.raises(ValueError, match='Number of PostLabelingDelays'):
                results = interface.run(cwd=tmpdir)

            # Scenario 3: PCASL with one PostLabelingDelay for each deltam volume (good)
            # This should produce ATT and mean CBF volumes, but no CBF time series
            metadata = {
                'ArterialSpinLabelingType': asltype,
                'PostLabelingDelay': good_multiple_plds,
                **BASE_METADATA,
                **acq_dict,
            }

            interface = cbf.ComputeCBF(
                cbf_only=False,
                deltam=asl_file,
                metadata=metadata,
                m0_scale=1,
                m0_file=m0_file,
            )
            results = interface.run(cwd=tmpdir)
            assert results.outputs.cbf_ts is None
            assert os.path.isfile(results.outputs.mean_cbf)
            mean_cbf_img = nb.load(results.outputs.mean_cbf)
            assert mean_cbf_img.ndim == 3
            assert os.path.isfile(results.outputs.att)
            att_img = nb.load(results.outputs.att)
            assert att_img.ndim == 3
            # The bounds enclose the estimates
            for name in ('att', 'abat'):
                estimate = nb.load(getattr(results.outputs, name)).get_fdata()
                lower = nb.load(getattr(results.outputs, f'{name}_lower')).get_fdata()
                upper = nb.load(getattr(results.outputs, f'{name}_upper')).get_fdata()
                # Voxels where the model could not be fit are NaN in every map
                fit = np.isfinite(estimate)
                assert fit.any()
                assert np.all(lower[fit] <= estimate[fit])
                assert np.all(estimate[fit] <= upper[fit])

            # The identifiable arterial signal is no larger than aBV, which it scales down
            abv = nb.load(results.outputs.abv).get_fdata()
            abv_decayed = nb.load(results.outputs.abv_decayed).get_fdata()
            fit = np.isfinite(abv)
            assert np.all(abv_decayed[fit] <= abv[fit])


def test_computecbf_pasl(datasets, tmp_path_factory):
    """Test aslprep.interfaces.cbf.ComputeCBF with PASL."""
    tmpdir = tmp_path_factory.mktemp('test_computecbf_pasl')
    aslcontext_file = os.path.join(
        datasets['nonqtab'],
        'sub-01/perf/sub-01_acq-pcaslsinglepldsiemens3d_run-02_aslcontext.tsv',
    )

    aslcontext = pd.read_table(aslcontext_file)
    n_deltam = aslcontext.loc[aslcontext['volume_type'] == 'label'].shape[0]
    n_volumes = aslcontext.shape[0]

    # Simulate ASL data and a brain mask.
    asl_data = np.random.random((30, 30, 30, n_deltam)).astype(np.float32)
    asl_file = _save_img(asl_data, tmpdir, 'asl.nii.gz')
    asl_mask = np.zeros((30, 30, 30), dtype=np.uint8)
    asl_mask[10:20, 10:20, 10:20] = 1
    m0_file = _save_img(asl_mask, tmpdir, 'm0.nii.gz')

    single_pld = 1.5
    plds = np.zeros(n_volumes)
    temp_plds = np.linspace(0.5, 3.5, n_deltam)
    plds[aslcontext['volume_type'] == 'm0scan'] = 0
    plds[aslcontext['volume_type'] == 'label'] = temp_plds
    plds[aslcontext['volume_type'] == 'control'] = temp_plds
    bad_multiple_plds = plds.tolist()
    good_multiple_plds = plds[aslcontext['volume_type'] == 'control']

    BASE_METADATA = {
        'ArterialSpinLabelingType': 'PASL',
        'MagneticFieldStrength': 3,
    }
    ACQ_DICTS = [
        {'MRAcquisitionType': '3D'},
        {
            'MRAcquisitionType': '2D',
            'SliceTiming': list(np.linspace(0.1, 0.5, 30)),
        },
    ]

    for acq_dict in ACQ_DICTS:
        # Scenario 1: PASL without BolusCutOff (raises ValueError).
        metadata = {
            'BolusCutOffFlag': False,
            'PostLabelingDelay': single_pld,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
        )
        with pytest.raises(ValueError, match='not supported in ASLPrep'):
            results = interface.run(cwd=tmpdir)

        # Scenario 2: QUIPSS PASL with a single PostLabelingDelay
        # This should produce CBF time series and mean CBF, but no ATT
        metadata = {
            'BolusCutOffFlag': True,
            'BolusCutOffTechnique': 'QUIPSS',
            'BolusCutOffDelayTime': 0.5,
            'PostLabelingDelay': single_pld,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
        )
        results = interface.run(cwd=tmpdir)
        assert os.path.isfile(results.outputs.cbf_ts)
        cbf_img = nb.load(results.outputs.cbf_ts)
        assert cbf_img.ndim == 4
        assert cbf_img.shape[3] == n_deltam
        assert os.path.isfile(results.outputs.mean_cbf)
        mean_cbf_img = nb.load(results.outputs.mean_cbf)
        assert mean_cbf_img.ndim == 3

        # Scenario 3: QUIPSS PASL with one PostLabelingDelay for each volume (bad)
        metadata = {
            'BolusCutOffFlag': True,
            'BolusCutOffTechnique': 'QUIPSS',
            'BolusCutOffDelayTime': 0.5,
            'PostLabelingDelay': bad_multiple_plds,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
        )
        with pytest.raises(ValueError, match=r'Unsupported BolusCutOffTechnique \(QUIPSS\)'):
            results = interface.run(cwd=tmpdir)

        # Scenario 4: QUIPSS PASL with one PostLabelingDelay for each deltam volume (good)
        metadata = {
            'BolusCutOffFlag': True,
            'BolusCutOffTechnique': 'QUIPSS',
            'BolusCutOffDelayTime': 0.5,
            'PostLabelingDelay': good_multiple_plds,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
        )
        with pytest.raises(ValueError, match=r'Unsupported BolusCutOffTechnique \(QUIPSS\)'):
            results = interface.run(cwd=tmpdir)

        # Scenario 5: QUIPSSII PASL with one PostLabelingDelay
        # This should produce CBF time series and mean CBF, but no ATT
        metadata = {
            'BolusCutOffFlag': True,
            'BolusCutOffTechnique': 'QUIPSSII',
            'BolusCutOffDelayTime': 0.5,
            'PostLabelingDelay': single_pld,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
        )
        results = interface.run(cwd=tmpdir)
        assert os.path.isfile(results.outputs.cbf_ts)
        cbf_img = nb.load(results.outputs.cbf_ts)
        assert cbf_img.ndim == 4
        assert cbf_img.shape[3] == n_deltam
        assert os.path.isfile(results.outputs.mean_cbf)
        mean_cbf_img = nb.load(results.outputs.mean_cbf)
        assert mean_cbf_img.ndim == 3

        # Scenario 6: QUIPSSII PASL with multiple PostLabelingDelays
        metadata = {
            'BolusCutOffFlag': True,
            'BolusCutOffTechnique': 'QUIPSSII',
            'BolusCutOffDelayTime': 0.5,
            'PostLabelingDelay': good_multiple_plds,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
        )
        results = interface.run(cwd=tmpdir)
        assert os.path.isfile(results.outputs.mean_cbf)
        mean_cbf_img = nb.load(results.outputs.mean_cbf)
        assert mean_cbf_img.ndim == 3
        assert os.path.isfile(results.outputs.att_lower)
        assert os.path.isfile(results.outputs.abat_upper)
        status_img = nb.load(results.outputs.att_status)
        assert status_img.get_data_dtype() == np.uint8
        assert status_img.shape == mean_cbf_img.shape
        assert set(np.unique(status_img.dataobj)) <= set(ATT_STATUS.values())
        # The PASL arterial signal does not depend on aBAT
        assert results.outputs.abv_decayed is None

        # Scenario 7: Q2TIPS PASL with one PostLabelingDelay
        # This should produce CBF time series and mean CBF, but no ATT
        metadata = {
            'BolusCutOffFlag': True,
            'BolusCutOffTechnique': 'Q2TIPS',
            'BolusCutOffDelayTime': [0.7, 1.6],
            'PostLabelingDelay': single_pld,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
        )
        results = interface.run(cwd=tmpdir)
        assert os.path.isfile(results.outputs.cbf_ts)
        cbf_img = nb.load(results.outputs.cbf_ts)
        assert cbf_img.ndim == 4
        assert cbf_img.shape[3] == n_deltam
        assert os.path.isfile(results.outputs.mean_cbf)
        mean_cbf_img = nb.load(results.outputs.mean_cbf)
        assert mean_cbf_img.ndim == 3

        # Scenario 8: Q2TIPS PASL with multiple PostLabelingDelays
        metadata = {
            'BolusCutOffFlag': True,
            'BolusCutOffTechnique': 'Q2TIPS',
            'BolusCutOffDelayTime': [0.7, 1.6],
            'PostLabelingDelay': good_multiple_plds,
            **BASE_METADATA,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
        )
        results = interface.run(cwd=tmpdir)
        assert os.path.isfile(results.outputs.mean_cbf)
        mean_cbf_img = nb.load(results.outputs.mean_cbf)
        assert mean_cbf_img.ndim == 3


@pytest.mark.parametrize(
    'acq_dict',
    [
        {'MRAcquisitionType': '3D'},
        {'MRAcquisitionType': '2D', 'SliceTiming': list(np.linspace(0, 0.25, 30))},
    ],
    ids=['3D', '2D'],
)
def test_computecbf_q2tips_uses_inversion_time(acq_dict, tmp_path):
    """Single-delay Q2TIPS decays over the (slice-shifted) inversion time, like QUIPSS II.

    The second BolusCutOffDelayTime is the end of the saturation pulse train,
    so it must not affect CBF.
    """
    asl_data = np.random.random((30, 30, 30, 4)).astype(np.float32)
    asl_file = _save_img(asl_data, tmp_path, 'asl.nii.gz')
    asl_mask = np.zeros((30, 30, 30), dtype=np.uint8)
    asl_mask[10:20, 10:20, 10:20] = 1
    m0_file = _save_img(asl_mask, tmp_path, 'm0.nii.gz')

    base_metadata = {
        'ArterialSpinLabelingType': 'PASL',
        'MagneticFieldStrength': 3,
        'BolusCutOffFlag': True,
        'PostLabelingDelay': 1.8,
        **acq_dict,
    }
    bolus_cutoffs = {
        'quipss2': {'BolusCutOffTechnique': 'QUIPSSII', 'BolusCutOffDelayTime': 0.7},
        'q2tips': {'BolusCutOffTechnique': 'Q2TIPS', 'BolusCutOffDelayTime': [0.7, 1.6]},
        'q2tips_short': {'BolusCutOffTechnique': 'Q2TIPS', 'BolusCutOffDelayTime': [0.7, 1.2]},
    }
    mean_cbf = {}
    for name, bolus_cutoff in bolus_cutoffs.items():
        out_dir = tmp_path / name
        out_dir.mkdir()
        results = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata={**base_metadata, **bolus_cutoff},
            m0_scale=1,
            m0_file=m0_file,
        ).run(cwd=out_dir)
        mean_cbf[name] = nb.load(results.outputs.mean_cbf).get_fdata()

    assert mean_cbf['quipss2'].max() > 0
    np.testing.assert_allclose(mean_cbf['q2tips'], mean_cbf['quipss2'])
    np.testing.assert_allclose(mean_cbf['q2tips_short'], mean_cbf['quipss2'])


@pytest.mark.parametrize(
    'labeling',
    [
        {'ArterialSpinLabelingType': 'PCASL', 'LabelingDuration': 1.6},
        {
            'ArterialSpinLabelingType': 'PASL',
            'BolusCutOffFlag': True,
            'BolusCutOffTechnique': 'QUIPSSII',
            'BolusCutOffDelayTime': 0.7,
        },
        {
            'ArterialSpinLabelingType': 'PASL',
            'BolusCutOffFlag': True,
            'BolusCutOffTechnique': 'QUIPSS',
            'BolusCutOffDelayTime': 0.7,
        },
        {
            'ArterialSpinLabelingType': 'PASL',
            'BolusCutOffFlag': True,
            'BolusCutOffTechnique': 'Q2TIPS',
            'BolusCutOffDelayTime': [0.7, 1.6],
        },
    ],
    ids=['pcasl', 'quipss2', 'quipss', 'q2tips'],
)
@pytest.mark.parametrize(
    'acq_dict',
    [
        {'MRAcquisitionType': '3D'},
        {'MRAcquisitionType': '2D', 'SliceTiming': list(np.linspace(0, 0.5, 12))},
    ],
    ids=['3D', '2D'],
)
@pytest.mark.parametrize('m0tr', [None, 3.0], ids=['longTR', 'shortTR'])
def test_computecbf_calibration_outputs(labeling, acq_dict, m0tr, tmp_path):
    """The calibration outputs reproduce the CBF quantification (PVC spec test 6).

    Mean CBF times the calibration M0 is the time-averaged q * delta-M (the M0-free signal
    regressed by partial volume correction), with q and the TR correction computed
    independently (:mod:`aslprep.tests.truth_models`). Voxels zeroed for low M0, and voxels
    with a non-finite delta-M sample, are not valid.
    """
    from aslprep.tests import truth_models as tm

    rng = np.random.default_rng(0)
    shape = (10, 11, 12)
    deltam = rng.normal(50, 10, (*shape, 4)).astype(np.float32)
    deltam[5, 5, 5, 1] = np.nan  # a missing sample: the masker would turn it into zero
    deltam_file = _save_img(deltam, tmp_path, 'deltam.nii.gz')
    m0 = rng.uniform(800, 1200, shape).astype(np.float32)
    m0[:2] = 1.0  # almost no signal: zeroed, so not quantifiable
    m0_file = _save_img(m0, tmp_path, 'm0.nii.gz')

    metadata = {'MagneticFieldStrength': 3, 'PostLabelingDelay': 1.8, **labeling, **acq_dict}
    kwargs = {'m0tr': m0tr} if m0tr else {}
    results = cbf.ComputeCBF(
        cbf_only=False,
        deltam=deltam_file,
        metadata=metadata,
        m0_scale=2.0,
        m0_file=m0_file,
        write_calibration=True,
        **kwargs,
    ).run(cwd=tmp_path)

    mean_cbf = nb.load(results.outputs.mean_cbf).get_fdata()
    m0_calib = nb.load(results.outputs.m0_calib).get_fdata()
    quant_valid = nb.load(results.outputs.quant_valid).get_fdata().astype(bool)
    expected_valid = np.ones(shape, bool)
    expected_valid[:2] = False
    expected_valid[5, 5, 5] = False
    np.testing.assert_array_equal(quant_valid, expected_valid)
    assert np.all(m0_calib[:2] == 0)

    # The calibration is the scaled M0, corrected for incomplete recovery when TR < 5 s.
    expected_m0 = tm.m0_tr_correction(2.0 * m0.astype(float), m0tr, 3)
    np.testing.assert_allclose(m0_calib[2:], expected_m0[2:], rtol=1e-5)

    # mean CBF * M0 is q * mean(delta-M), with q from the independent model.
    pld = tm.pld_map(metadata, shape)
    alpha = tm.labeling_efficiency(metadata)
    if labeling.get('BolusCutOffTechnique') == 'QUIPSS':
        # Wong 1998: the bolus duration is TI - TI1, per voxel
        t1b = tm.T1_BLOOD[3]
        q = (
            tm.UNIT_CONVERSION
            * tm.LAMBDA
            * np.exp(pld / t1b)
            / (2 * alpha * (pld - labeling['BolusCutOffDelayTime']))
        )
    else:
        q = tm.cbf_single_delay(np.ones(shape), np.ones(shape), metadata, alpha, pld)
    signal = mean_cbf * m0_calib
    expected = q * deltam.astype(float).mean(axis=3)
    np.testing.assert_allclose(signal[quant_valid], expected[quant_valid], rtol=1e-4)


def test_computecbf_calibration_outputs_off_by_default(tmp_path):
    """Without write_calibration (no partial volume correction), nothing extra is written."""
    rng = np.random.default_rng(0)
    deltam_file = _save_img(
        rng.normal(50, 10, (6, 6, 6, 2)).astype(np.float32), tmp_path, 'd.nii.gz'
    )
    m0_file = _save_img(np.full((6, 6, 6), 1000, np.float32), tmp_path, 'm0.nii.gz')
    results = cbf.ComputeCBF(
        cbf_only=False,
        deltam=deltam_file,
        metadata={
            'ArterialSpinLabelingType': 'PCASL',
            'LabelingDuration': 1.8,
            'MagneticFieldStrength': 3,
            'MRAcquisitionType': '3D',
            'PostLabelingDelay': 1.8,
        },
        m0_scale=1,
        m0_file=m0_file,
    ).run(cwd=tmp_path)
    assert results.outputs.m0_calib is None
    assert results.outputs.quant_valid is None
    assert not list(tmp_path.glob('*_m0calib*'))
    assert not list(tmp_path.glob('*_quantvalid*'))


def test_computecbf_calibration_outputs_multipld(tmp_path):
    """Multi-delay data get the calibration outputs too (PVC spec 5.7.2).

    A non-finite sample in any volume makes the voxel invalid, even though the masker would
    replace it with zero.
    """
    rng = np.random.default_rng(0)
    deltam = rng.normal(50, 10, (6, 6, 6, 3)).astype(np.float32)
    deltam[2, 2, 2, 1] = np.nan
    deltam_file = _save_img(deltam, tmp_path, 'd.nii.gz')
    m0 = np.full((6, 6, 6), 1000, np.float32)
    m0[0] = 1.0  # zeroed as having no signal
    m0_file = _save_img(m0, tmp_path, 'm0.nii.gz')
    results = cbf.ComputeCBF(
        cbf_only=False,
        deltam=deltam_file,
        metadata={
            'ArterialSpinLabelingType': 'PCASL',
            'LabelingDuration': 1.6,
            'MagneticFieldStrength': 3,
            'MRAcquisitionType': '3D',
            'PostLabelingDelay': [1.0, 1.5, 2.0],
        },
        m0_scale=1,
        m0_file=m0_file,
        write_calibration=True,
    ).run(cwd=tmp_path)
    quant_valid = nb.load(results.outputs.quant_valid).get_fdata() > 0
    m0_calib = nb.load(results.outputs.m0_calib).get_fdata()
    expected = np.ones((6, 6, 6), bool)
    expected[0] = False
    expected[2, 2, 2] = False
    np.testing.assert_array_equal(quant_valid, expected)
    np.testing.assert_allclose(m0_calib[1:], 1000.0)
    assert np.all(m0_calib[0] == 0)

    # Voxels without M0 could not be fit, and have their own ATT status
    status_img = nb.load(results.outputs.att_status)
    assert status_img.get_data_dtype() == np.uint8
    status = np.asarray(status_img.dataobj)
    np.testing.assert_array_equal(status == ATT_STATUS['not_fitted'], m0_calib == 0)
    assert set(np.unique(status)) <= set(ATT_STATUS.values())


def test_compare_slicetiming(datasets, tmp_path_factory):
    """CBF estimates should be the same with slice-timing modification as not.

    As long as the slice times are all zero, of course.
    """
    tmpdir = tmp_path_factory.mktemp('test_computecbf_casl')
    aslcontext_file = os.path.join(
        datasets['nonqtab'],
        'sub-01/perf/sub-01_acq-pcaslsinglepldsiemens3d_run-02_aslcontext.tsv',
    )

    aslcontext = pd.read_table(aslcontext_file)
    n_deltam = aslcontext.loc[aslcontext['volume_type'] == 'label'].shape[0]

    # Simulate ASL data and a brain mask.
    asl_data = np.random.random((30, 30, 30, n_deltam)).astype(np.float32)
    asl_file = _save_img(asl_data, tmpdir, 'asl.nii.gz')
    asl_mask = np.zeros((30, 30, 30), dtype=np.uint8)
    asl_mask[10:20, 10:20, 10:20] = 1
    m0_file = _save_img(asl_mask, tmpdir, 'm0.nii.gz')

    ACQ_DICTS = [
        {'MRAcquisitionType': '3D'},
        {
            'MRAcquisitionType': '2D',
            'SliceTiming': list(np.zeros(30)),
        },
    ]

    cbf_data = []
    for acq_dict in ACQ_DICTS:
        metadata = {
            'ArterialSpinLabelingType': 'PCASL',
            'MagneticFieldStrength': 3,
            'LabelingDuration': 1.6,
            'PostLabelingDelay': 1.5,
            **acq_dict,
        }
        interface = cbf.ComputeCBF(
            cbf_only=False,
            deltam=asl_file,
            metadata=metadata,
            m0_scale=1,
            m0_file=m0_file,
        )
        results = interface.run(cwd=tmpdir)
        cbf_data.append(nb.load(results.outputs.cbf_ts).get_fdata())

    assert np.array_equal(cbf_data[0], cbf_data[1])


def test_basilcbf_tmpdir(tmp_path, monkeypatch):
    """oxford_asl writes its temporary files to the working directory, not to /tmp (gh-557)."""
    # A stand-in for oxford_asl that records where FSL's tmpnam would create its files,
    # and writes the outputs that BASILCBF expects.
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    fake_oxford_asl = bin_dir / 'oxford_asl'
    fake_oxford_asl.write_text(
        '#!/bin/sh\n'
        'echo "$TMPDIR" > tmpdir.txt\n'
        'mkdir -p native_space/pvcorr\n'
        'touch native_space/perfusion_calib.nii.gz\n'
        'touch native_space/pvcorr/perfusion_calib.nii.gz\n'
        'touch native_space/pvcorr/perfusion_wm_calib.nii.gz\n'
    )
    fake_oxford_asl.chmod(0o755)
    monkeypatch.setenv('PATH', f'{bin_dir}{os.pathsep}{os.environ["PATH"]}')
    monkeypatch.setenv('TMPDIR', str(tmp_path / 'elsewhere'))

    shape = (4, 4, 4)
    deltam = _save_img(np.ones(shape, dtype=np.float32), tmp_path, 'deltam.nii.gz')
    mask = _save_img(np.ones(shape, dtype=np.uint8), tmp_path, 'mask.nii.gz')
    work_dir = tmp_path / 'work'
    work_dir.mkdir()

    results = cbf.BASILCBF(
        deltam=deltam,
        mask=mask,
        m0_scale=1,
        tis=[3.6],
        bolus=1.8,
        pcasl=True,
    ).run(cwd=str(work_dir))

    assert (work_dir / 'tmpdir.txt').read_text().strip() == str(work_dir)
    assert results.outputs.mean_cbf_basil == str(work_dir / 'native_space/perfusion_calib.nii.gz')


@pytest.mark.parametrize('pvc', [True, False])
def test_basilcbf_pvc_outputs(tmp_path, monkeypatch, pvc):
    """BASILCBF only reports partial-volume corrected maps when oxford_asl writes them."""
    # A stand-in for oxford_asl that writes the outputs that the real command would.
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    fake_oxford_asl = bin_dir / 'oxford_asl'
    fake_oxford_asl.write_text(
        '#!/bin/sh\n'
        'echo "$@" > args.txt\n'
        'mkdir -p native_space\n'
        'touch native_space/perfusion_calib.nii.gz\n'
        'case "$*" in *--pvcorr*)\n'
        '    mkdir -p native_space/pvcorr\n'
        '    touch native_space/pvcorr/perfusion_calib.nii.gz\n'
        '    touch native_space/pvcorr/perfusion_wm_calib.nii.gz\n'
        'esac\n'
    )
    fake_oxford_asl.chmod(0o755)
    monkeypatch.setenv('PATH', f'{bin_dir}{os.pathsep}{os.environ["PATH"]}')

    shape = (4, 4, 4)
    work_dir = tmp_path / 'work'
    work_dir.mkdir()
    results = cbf.BASILCBF(
        deltam=_save_img(np.ones(shape, dtype=np.float32), tmp_path, 'deltam.nii.gz'),
        mask=_save_img(np.ones(shape, dtype=np.uint8), tmp_path, 'mask.nii.gz'),
        m0_scale=1,
        tis=[3.6],
        bolus=1.8,
        pcasl=True,
        pvc=pvc,
    ).run(cwd=str(work_dir))

    assert ('--pvcorr' in (work_dir / 'args.txt').read_text()) is pvc
    assert results.outputs.mean_cbf_basil == str(work_dir / 'native_space/perfusion_calib.nii.gz')
    assert isdefined(results.outputs.mean_cbf_gm_basil) is pvc
    assert isdefined(results.outputs.mean_cbf_wm_basil) is pvc


@pytest.mark.parametrize(('pvc', 'multi_delay'), [(True, True), (True, False), (False, True)])
def test_basilcbf_posterior_sd(tmp_path, monkeypatch, pvc, multi_delay):
    """BASILCBF turns oxford_asl's posterior variance maps into standard deviations (gh-671)."""
    shape = (4, 4, 4)
    mask = np.zeros(shape, dtype=np.uint8)
    mask[1:3, 1:3, 1:3] = 1
    inside = mask.astype(bool)

    # What oxford_asl would write: calibrated CBF variances that are 0 outside the mask,
    # and an uncalibrated ATT variance with its 1e12 filler for voxels it did not fit.
    canned = tmp_path / 'canned'
    variances = {
        'perfusion_var_calib': np.where(inside, 25.0, 0.0),
    }
    if multi_delay:
        variances['arrival_var'] = np.where(inside, 0.04, 1e12)
    if pvc:
        variances['pvcorr/perfusion_var_calib'] = np.where(inside, 9.0, 0.0)
        variances['pvcorr/perfusion_wm_var_calib'] = np.where(inside, 16.0, 0.0)
    (canned / 'native_space' / 'pvcorr').mkdir(parents=True)
    _save_img(mask, str(canned / 'native_space'), 'mask.nii.gz')
    _save_img(np.ones(shape, np.float32), str(canned / 'native_space'), 'perfusion_calib.nii.gz')
    for name, data in variances.items():
        _save_img(data.astype(np.float32), str(canned / 'native_space'), f'{name}.nii.gz')
    if pvc:
        for name in ('perfusion_calib', 'perfusion_wm_calib'):
            _save_img(
                np.ones(shape, np.float32),
                str(canned / 'native_space' / 'pvcorr'),
                f'{name}.nii.gz',
            )

    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    fake_oxford_asl = bin_dir / 'oxford_asl'
    fake_oxford_asl.write_text(f'#!/bin/sh\ncp -r {canned}/native_space .\n')
    fake_oxford_asl.chmod(0o755)
    monkeypatch.setenv('PATH', f'{bin_dir}{os.pathsep}{os.environ["PATH"]}')

    work_dir = tmp_path / 'work'
    work_dir.mkdir()
    out = (
        cbf.BASILCBF(
            deltam=_save_img(np.ones(shape, dtype=np.float32), str(tmp_path), 'deltam.nii.gz'),
            mask=_save_img(mask, str(tmp_path), 'mask.nii.gz'),
            m0_scale=1,
            tis=[3.6],
            bolus=1.8,
            pcasl=True,
            pvc=pvc,
        )
        .run(cwd=str(work_dir))
        .outputs
    )

    expected = {'cbf_sd_basil': 5.0}
    if pvc:
        expected.update({'cbf_gm_sd_basil': 3.0, 'cbf_wm_sd_basil': 4.0})
    if multi_delay:
        expected['att_sd_basil'] = 0.2
    for field, value in expected.items():
        img = nb.load(getattr(out, field))
        assert img.get_data_dtype() == np.float32
        data = img.get_fdata()
        assert np.allclose(data[inside], value)
        # Zero outside the mask, including where oxford_asl wrote its 1e12 filler.
        assert np.all(data[~inside] == 0)

    assert isdefined(out.cbf_gm_sd_basil) is pvc
    assert isdefined(out.cbf_wm_sd_basil) is pvc
    assert isdefined(out.att_sd_basil) is multi_delay


@pytest.mark.parametrize(
    ('direction', 'tis'),
    [('k', [3.6]), ('k-', [2.05, 2.55, 3.05]), ('j', [3.6]), ('i-', [2.05, 2.55])],
)
def test_basil_inflow_time_image(tmp_path, direction, tis):
    """Each voxel's inflow time is shifted by the acquisition time of its slice (gh-743)."""
    shape = (3, 4, 5)
    axis = 'ijk'.index(direction[0])
    n_slices = shape[axis]
    # Interleaved: even slices first, then odd slices
    order = list(range(0, n_slices, 2)) + list(range(1, n_slices, 2))
    slice_times = np.zeros(n_slices)
    slice_times[order] = 0.045 * np.arange(n_slices)
    deltam = _save_img(np.ones(shape + (2,), dtype=np.float32), tmp_path, 'deltam.nii.gz')

    out = (
        cbf.BASILInflowTimeImage(
            deltam=deltam,
            tis=tis,
            metadata={
                'SliceTiming': slice_times.tolist(),
                'SliceEncodingDirection': direction,
            },
        )
        .run(cwd=str(tmp_path))
        .outputs
    )

    img = nb.load(out.out_file)
    assert img.shape == shape + (len(tis),)
    assert img.get_data_dtype() == np.float32
    data = img.get_fdata()
    # Slice index i along the slice axis was acquired at slice_times[i], counting from the
    # far end of the axis for a "-" direction.
    times_by_index = slice_times[::-1] if direction.endswith('-') else slice_times
    for i_slice in range(n_slices):
        slice_data = np.take(data, i_slice, axis=axis)
        for i_ti, ti in enumerate(tis):
            assert np.allclose(slice_data[..., i_ti], ti + times_by_index[i_slice])


def test_basilcbf_tiimg_cmdline(tmp_path):
    """A voxelwise inflow time image replaces --tis and --slicedt in the oxford_asl call."""
    shape = (4, 4, 4)
    deltam = _save_img(np.ones(shape, dtype=np.float32), tmp_path, 'deltam.nii.gz')
    mask = _save_img(np.ones(shape, dtype=np.uint8), tmp_path, 'mask.nii.gz')
    tiimg = _save_img(np.full(shape, 3.6, dtype=np.float32), tmp_path, 'tiimg.nii.gz')
    common = {'deltam': deltam, 'mask': mask, 'm0_scale': 1, 'bolus': 1.8, 'pcasl': True}

    cmdline = cbf.BASILCBF(tiimg=tiimg, **common).cmdline
    assert f'--tiimg {tiimg}' in cmdline
    assert '--tis' not in cmdline

    with pytest.raises(OSError, match='mutually exclusive'):
        cbf.BASILCBF(tiimg=tiimg, tis=[3.6], **common)
    with pytest.raises(OSError, match='mutually exclusive'):
        cbf.BASILCBF(tiimg=tiimg, slice_spacing=0.045, **common)
    with pytest.raises(ValueError, match='requires a value'):
        cbf.BASILCBF(**common).cmdline  # noqa: B018


def _save_img(arr, tmpdir, fname):
    img = nb.Nifti1Image(arr, affine=np.eye(4))
    fname = os.path.join(tmpdir, fname)
    img.to_filename(fname)
    return fname
