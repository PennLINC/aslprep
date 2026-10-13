"""Tests for the aslprep.interfaces.cbf submodule."""

import os

import nibabel as nb
import numpy as np
import pandas as pd
import pytest

from aslprep.interfaces import cbf


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


def _save_img(arr, tmpdir, fname):
    img = nb.Nifti1Image(arr, affine=np.eye(4))
    fname = os.path.join(tmpdir, fname)
    img.to_filename(fname)
    return fname
