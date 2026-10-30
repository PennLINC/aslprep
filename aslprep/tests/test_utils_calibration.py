"""Tests for the resolution of calibration and physiological parameters."""

import json
from pathlib import Path

import nibabel as nb
import numpy as np
import pytest

from aslprep.utils.calibration import (
    BASIL_TISSUE_T1,
    M0_RECOVERY_TR_THRESHOLD,
    MODELS,
    SOURCE_PRECEDENCE,
    TISSUE_T1_PARAMETERS,
    basil_calibration_inputs,
    calibration_kinetic_model,
    calibration_value,
    describe_model,
    m0_partition_coefficient,
    resolve_calibration,
    tissue_t1_kinetics,
)
from aslprep.utils.cbf import T1TISSUE_BY_FIELD_STRENGTH, T1WM_BY_FIELD_STRENGTH

PCASL = {
    'ArterialSpinLabelingType': 'PCASL',
    'M0Type': 'Separate',
    'MagneticFieldStrength': 3,
    'BackgroundSuppression': False,
    'PostLabelingDelay': 1.8,
    'LabelingDuration': 1.8,
}


def _source(calibration, name):
    return calibration['Parameters'][name]['Source']


def test_resolve_calibration_defaults():
    """Without metadata values, every parameter has a documented default or inference."""
    calibration = resolve_calibration(PCASL, m0tr=6.0)

    assert calibration['SourcePrecedence'] == list(SOURCE_PRECEDENCE)
    assert calibration['Warnings'] == []
    # Every value can be written to a JSON file.
    assert json.loads(json.dumps(calibration)) == calibration

    expected = {
        'M0Type': ('Separate', 'metadata'),
        'M0Scale': (1.0, 'default'),
        'M0SmoothingFWHM': (5.0, 'default'),
        'M0RepetitionTime': (6.0, 'metadata'),
        'M0RecoveryCorrection': ('none', 'inferred'),
        'LabelingEfficiency': (0.85, 'default'),
        'BackgroundSuppressionNumberPulses': (0, 'metadata'),
        'EffectiveLabelingEfficiency': (0.85, 'inferred'),
        'BloodT1': (1.65, 'default'),
        'PartitionCoefficient': (0.9, 'default'),
        'CBFUnitConversion': (6000, 'default'),
    }
    for name, (value, source) in expected.items():
        assert calibration_value(calibration, name) == pytest.approx(value), name
        assert _source(calibration, name) == source, name

    # A fully relaxed M0 needs no tissue T1, and BASIL-only parameters are left out.
    for name in ('TissueT1', 'M0RecoveryFactor', 'BASILTissueT1', 'M0Estimate'):
        assert name not in calibration['Parameters']


def test_resolve_calibration_overrides_and_metadata():
    """Command-line values and BIDS metadata are recorded with their sources."""
    metadata = {
        **PCASL,
        'LabelingEfficiency': 0.8,
        'BackgroundSuppression': True,
        'BackgroundSuppressionNumberPulses': 2,
    }
    calibration = resolve_calibration(
        metadata,
        m0tr=3.0,
        m0_scale=0.5,
        m0_scale_source='user',
        smooth_kernel=3.0,
        smooth_kernel_source='user',
    )

    assert calibration_value(calibration, 'M0Scale') == 0.5
    assert _source(calibration, 'M0Scale') == 'user'
    assert calibration['Parameters']['M0Scale']['Origin'] == '--m0-scale'
    assert _source(calibration, 'M0SmoothingFWHM') == 'user'
    assert _source(calibration, 'LabelingEfficiency') == 'metadata'
    assert calibration_value(calibration, 'BackgroundSuppressionEfficiency') == pytest.approx(
        0.95**2
    )
    assert calibration_value(calibration, 'EffectiveLabelingEfficiency') == pytest.approx(
        0.8 * 0.95**2
    )
    # A short M0 TR is corrected for T1 recovery with the gray matter T1.
    assert calibration_value(calibration, 'M0RecoveryCorrection') == 'saturation recovery'
    assert calibration_value(calibration, 'TissueT1') == pytest.approx(1.607)
    assert calibration_value(calibration, 'M0RecoveryFactor') == pytest.approx(
        1 / (1 - np.exp(-3.0 / 1.607))
    )


def test_resolve_calibration_matches_estimators():
    """The resolved values are the ones the existing estimators return."""
    from aslprep.utils.asl import estimate_labeling_efficiency
    from aslprep.utils.cbf import estimate_t1

    metadata = {
        **PCASL,
        'MagneticFieldStrength': 1.5,
        'BackgroundSuppression': True,
        'BackgroundSuppressionNumberPulses': 4,
        'Haematocrit': 0.35,
    }
    calibration = resolve_calibration(metadata, m0tr=2.0)

    t1blood, t1tissue = estimate_t1(metadata)
    assert calibration_value(calibration, 'BloodT1') == pytest.approx(t1blood)
    assert calibration_value(calibration, 'TissueT1') == pytest.approx(t1tissue)
    assert calibration_value(calibration, 'EffectiveLabelingEfficiency') == pytest.approx(
        estimate_labeling_efficiency(metadata, warn=False)
    )


def test_resolve_calibration_pulse_count_from_pulse_times():
    """The number of background suppression pulses is inferred from their times."""
    from aslprep.utils.asl import describe_labeling_efficiency, estimate_labeling_efficiency

    metadata = {
        **PCASL,
        'BackgroundSuppression': True,
        'BackgroundSuppressionPulseTime': [0.5, 1.2, 1.6],
    }
    calibration = resolve_calibration(metadata, m0tr=6.0)

    assert calibration_value(calibration, 'BackgroundSuppressionNumberPulses') == 3
    assert _source(calibration, 'BackgroundSuppressionNumberPulses') == 'inferred'
    assert calibration['Warnings'] == []
    assert estimate_labeling_efficiency(metadata) == pytest.approx(0.85 * 0.95**3)
    assert 'counted from the background suppression pulse times' in (
        describe_labeling_efficiency(metadata)
    )


def test_resolve_calibration_m0tr_rounded():
    """ASLPrep and BASIL decide on T1-recovery correction from the same rounded M0 TR."""
    calibration = resolve_calibration(PCASL, m0tr=4.9999996)
    assert calibration_value(calibration, 'M0RepetitionTime') == 5.0
    assert calibration_value(calibration, 'M0RecoveryCorrection') == 'none'


def test_describe_asl_acquisition_inferred_pulses():
    """The acquisition description gives the pulse count inferred from the pulse times."""
    import pandas as pd

    from aslprep.utils.asl import describe_asl_acquisition

    metadata = {
        **PCASL,
        'BackgroundSuppression': True,
        'BackgroundSuppressionPulseTime': [0.5, 1.2],
    }
    description = describe_asl_acquisition(metadata, pd.DataFrame({'volume_type': ['control']}))
    assert 'with background suppression (2 pulses)' in description


@pytest.mark.parametrize(
    ('extra', 'warning'),
    [
        (
            {'BackgroundSuppression': True},
            'One background suppression pulse is assumed',
        ),
        (
            {'BackgroundSuppression': False, 'BackgroundSuppressionNumberPulses': 2},
            'BackgroundSuppression is false or not defined',
        ),
        (
            {'BackgroundSuppression': True, 'BackgroundSuppressionNumberPulses': 0},
            'BackgroundSuppressionNumberPulses is 0',
        ),
        (
            {
                'BackgroundSuppression': True,
                'BackgroundSuppressionNumberPulses': 2,
                'BackgroundSuppressionPulseTime': [0.5, 1.2, 1.6],
            },
            'does not match the number of BackgroundSuppressionPulseTime values',
        ),
        ({'LabelingEfficiency': 0.3}, 'unusually low'),
        ({'Haematocrit': 0.8}, 'outside the plausible range'),
        ({'MagneticFieldStrength': 2.89}, 'formula of Zhang et al. (2013)'),
        ({'MagneticFieldStrength': 9.4, 'Haematocrit': 0.4}, 'was not used'),
    ],
)
def test_resolve_calibration_warnings(extra, warning):
    """Missing, conflicting and implausible inputs are flagged."""
    calibration = resolve_calibration({**PCASL, **extra}, m0tr=6.0)
    assert any(warning in message for message in calibration['Warnings']), calibration['Warnings']


def test_resolve_calibration_m0scan_field_strength_conflict():
    """An M0 scan from a different field strength is flagged."""
    calibration = resolve_calibration(
        PCASL, m0tr=6.0, m0scan_metadata={'MagneticFieldStrength': 1.5}
    )
    assert any('The M0 scan was acquired at 1.5 T' in w for w in calibration['Warnings'])


@pytest.mark.parametrize(
    ('extra', 'kwargs', 'match'),
    [
        ({'LabelingEfficiency': 1.2}, {}, 'LabelingEfficiency must be a number'),
        ({'LabelingEfficiency': 0}, {}, 'LabelingEfficiency must be a number'),
        ({}, {'m0_scale': 0.0}, '--m0-scale must be'),
        ({'M0Type': 'Estimate', 'M0Estimate': -1}, {}, 'M0Estimate must be'),
        ({'MagneticFieldStrength': 0}, {}, 'MagneticFieldStrength must be'),
        (
            {'BackgroundSuppression': True, 'BackgroundSuppressionNumberPulses': 1.5},
            {},
            'non-negative integer',
        ),
    ],
)
def test_resolve_calibration_invalid(extra, kwargs, match):
    """Values that would make CBF meaningless are rejected."""
    with pytest.raises(ValueError, match=match):
        resolve_calibration({**PCASL, **extra}, m0tr=6.0, **kwargs)


def test_resolve_calibration_tissue_t1_only_when_needed():
    """Tissue T1 is only required, and only defined, for an M0 with a short TR."""
    metadata = {**PCASL, 'MagneticFieldStrength': 2.89}
    calibration = resolve_calibration(metadata, m0tr=6.0)
    assert 'TissueT1' not in calibration['Parameters']
    with pytest.raises(ValueError, match='T1tissue cannot be inferred'):
        resolve_calibration(metadata, m0tr=3.0)


def test_resolve_calibration_haematocrit_source():
    """A haematocrit from a BIDS TSV file is recorded as metadata, with its file."""
    metadata = {**PCASL, 'Haematocrit': 0.35}
    calibration = resolve_calibration(metadata, m0tr=6.0, haematocrit_source='participants.tsv')

    assert _source(calibration, 'Haematocrit') == 'metadata'
    assert calibration['Parameters']['Haematocrit']['Origin'] == 'participants.tsv'
    assert _source(calibration, 'BloodT1') == 'inferred'
    assert 'Hales et al. (2016)' in calibration['Parameters']['BloodT1']['Origin']


@pytest.mark.parametrize('multi_delay', [False, True])
def test_describe_model(multi_delay):
    """Each model's sidecar dictionary has the calibration values that the model used."""
    calibration = resolve_calibration(PCASL, m0tr=6.0, basil=True)

    models = {
        model: describe_model(calibration, model, asl_type='PCASL', multi_delay=multi_delay)
        for model in MODELS
    }
    for model in models.values():
        assert set(model) == {'Description', 'TermURL', 'Parameters'}
        assert model['TermURL'].startswith('https://doi.org/')
        assert model['Parameters']['BloodT1'] == pytest.approx(1.65)
        assert model['Parameters']['M0Type'] == 'Separate'
        assert json.loads(json.dumps(model)) == model

    # Only BASIL uses its kinetic tissue T1, and it does not drop low-M0 voxels.
    assert models['basil']['Parameters']['BASILTissueT1'] == BASIL_TISSUE_T1
    assert 'M0Threshold' not in models['basil']['Parameters']
    for name in ('aslprep', 'lr'):
        assert 'BASILTissueT1' not in models[name]['Parameters']
        assert models[name]['Parameters']['M0Threshold'] == pytest.approx(0.1)
    assert ('multi-delay' in models['aslprep']['Description']) is multi_delay

    with pytest.raises(ValueError, match='Unknown model'):
        describe_model(calibration, 'other', asl_type='PCASL', multi_delay=multi_delay)


@pytest.mark.parametrize(
    ('name', 'model'),
    [
        ('mean_cbf', 'aslprep'),
        ('cbf_ts_score', 'aslprep'),
        ('att', 'aslprep'),
        ('mean_cbf_basil', 'basil'),
        ('cbf_wm_sd_basil', 'basil'),
        ('mean_cbf_gm_lr', 'lr'),
        ('mean_cbf_gm_lr_score', 'lr'),
        ('kinres_wm_lr', 'lr'),
    ],
)
def test_output_model(name, model):
    """Each output's sidecar gets the dictionary of the model that produced it."""
    from aslprep.workflows.asl.outputs import output_fields, output_model

    assert output_model(name) == model
    models = {key: {'Description': key} for key in MODELS}
    assert output_fields(name, models)['Model'] == {'Description': model}
    assert 'Model' not in output_fields(name)


@pytest.mark.parametrize(
    ('m0type', 'm0tr', 'expected_tr', 'expected_t1t'),
    [
        ('Separate', 3.0, 3.0, 1.607),
        ('Separate', 6.0, 6.0, None),
        ('Included', 4.9, 4.9, 1.607),
    ],
)
def test_basil_calibration_inputs(m0type, m0tr, expected_tr, expected_t1t):
    """BASIL is given the parameters ASLPrep's own models use."""
    metadata = {
        **PCASL,
        'M0Type': m0type,
        'LabelingEfficiency': 0.81,
        'BackgroundSuppression': True,
        'BackgroundSuppressionNumberPulses': 4,
    }
    calibration = resolve_calibration(metadata, m0tr=m0tr, m0_scale=0.001, basil=True)
    inputs = basil_calibration_inputs(calibration)

    assert inputs['alpha'] == calibration_value(calibration, 'EffectiveLabelingEfficiency')
    assert inputs['t1blood'] == calibration_value(calibration, 'BloodT1')
    assert inputs['m0_scale'] == 0.001
    assert inputs['t1tissue'] == BASIL_TISSUE_T1
    assert inputs['m0tr'] == pytest.approx(expected_tr)
    # oxford_asl corrects for T1 recovery exactly when ASLPrep does.
    assert (inputs['m0tr'] < M0_RECOVERY_TR_THRESHOLD) is (expected_t1t is not None)
    assert inputs.get('t1tissue_calib') == (pytest.approx(expected_t1t) if expected_t1t else None)
    assert 'm0_value' not in inputs


def test_basil_calibration_inputs_m0_estimate():
    """An M0Estimate goes to oxford_asl as a single blood M0 value, scaled by ASLPrep."""
    from aslprep.interfaces.cbf import BASILCBF

    metadata = {**PCASL, 'M0Type': 'Estimate', 'M0Estimate': 1234.5}
    calibration = resolve_calibration(metadata, m0_scale=0.001, basil=True)
    inputs = basil_calibration_inputs(calibration)

    # oxford_asl applies neither --cgain nor the T1 recovery correction to an --M0 value.
    assert inputs['m0_value'] == pytest.approx(1.2345)
    for name in ('m0_scale', 'm0tr', 't1tissue_calib'):
        assert name not in inputs
    assert inputs['alpha'] == calibration_value(calibration, 'EffectiveLabelingEfficiency')
    assert inputs['t1blood'] == calibration_value(calibration, 'BloodT1')

    args = BASILCBF(tis=[3.6], bolus=1.8, **inputs)._parse_inputs(skip=['deltam', 'mask'])
    assert '--M0 1.2345' in args
    assert not any(arg.startswith(('--cgain', '--tr', '-c ')) for arg in args)

    # bc, which multiplies the value in oxford_asl, does not read exponents.
    args = BASILCBF(tis=[3.6], bolus=1.8, m0_value=1.5e-7)._parse_inputs(skip=['deltam', 'mask'])
    assert '--M0 0.00000015' in args


def test_m0_estimate_is_blood_m0():
    """An M0Estimate (the M0 of blood) is not divided by the partition coefficient."""
    metadata = {**PCASL, 'M0Type': 'Estimate', 'M0Estimate': 1000.0}
    calibration = resolve_calibration(metadata)
    assert 'PartitionCoefficient' not in calibration['Parameters']
    assert m0_partition_coefficient(calibration) == 1.0
    assert 'partition coefficient' in calibration['Parameters']['M0Estimate']['Description']
    for model in MODELS:
        parameters = describe_model(calibration, model, asl_type='PCASL', multi_delay=False)
        assert 'PartitionCoefficient' not in parameters['Parameters']

    separate = resolve_calibration(PCASL, m0tr=6.0)
    assert m0_partition_coefficient(separate) == calibration_value(
        separate, 'PartitionCoefficient'
    )


@pytest.mark.parametrize('multi_delay', [False, True])
def test_computecbf_m0_estimate(tmp_path, multi_delay):
    """CBF from an M0Estimate is that from a uniform M0 image of lambda times the estimate."""
    from aslprep.interfaces.cbf import ComputeCBF
    from aslprep.utils.cbf import PARTITION_COEFFICIENT

    rng = np.random.default_rng(0)
    shape = (4, 4, 4)
    plds = [0.5, 1.0, 1.5, 2.0] if multi_delay else [1.8, 1.8]
    deltam = _write(rng.uniform(5, 15, (*shape, len(plds))), tmp_path / 'deltam.nii.gz')
    m0_estimate = 1000.0
    m0_file = _write(np.full(shape, m0_estimate), tmp_path / 'm0.nii.gz')
    tissue_m0_file = _write(
        np.full(shape, PARTITION_COEFFICIENT * m0_estimate), tmp_path / 'tissue_m0.nii.gz'
    )

    def _mean_cbf(name, m0, m0type):
        metadata = {**PCASL, 'PostLabelingDelay': plds, 'M0Type': m0type}
        if m0type == 'Estimate':
            metadata['M0Estimate'] = m0_estimate
        work_dir = tmp_path / name
        work_dir.mkdir()
        result = ComputeCBF(
            deltam=deltam,
            m0_file=m0,
            metadata=metadata,
            cbf_only=False,
            calibration=resolve_calibration(metadata, m0tr=None if m0type == 'Estimate' else 6.0),
        ).run(cwd=str(work_dir))
        return nb.load(result.outputs.mean_cbf).get_fdata()

    estimate = _mean_cbf('estimate', m0_file, 'Estimate')
    tissue = _mean_cbf('tissue', tissue_m0_file, 'Separate')
    assert np.any(estimate > 0)
    np.testing.assert_allclose(estimate, tissue, rtol=1e-6)


def test_basil_calibration_inputs_keep_precision():
    """oxford_asl gets the effective labeling efficiency and M0 scale without rounding."""
    from aslprep.interfaces.cbf import BASILCBF

    metadata = {
        **PCASL,
        'BackgroundSuppression': True,
        'BackgroundSuppressionNumberPulses': 4,
    }
    calibration = resolve_calibration(metadata, m0tr=3.0, m0_scale=0.00123, basil=True)
    args = BASILCBF(tis=[3.6], bolus=1.8, **basil_calibration_inputs(calibration))._parse_inputs(
        skip=['deltam', 'mask']
    )

    alpha = 0.85 * 0.95**4
    assert f'--alpha {alpha:.6f}' in args
    assert '--cgain 0.00123' in args
    assert '--tr 3.000000' in args
    assert '--t1t 1.607000' in args
    assert '--t1 1.300000' in args


def test_basil_calibration_inputs_partition_coefficient():
    """oxford_asl cannot use any partition coefficient but 0.9."""
    calibration = resolve_calibration(PCASL, m0tr=6.0)
    calibration['Parameters']['PartitionCoefficient']['Value'] = 0.98
    with pytest.raises(ValueError, match='fixed blood-brain partition coefficient'):
        basil_calibration_inputs(calibration)


def _write(data, path):
    nb.Nifti1Image(np.asarray(data, dtype=np.float32), np.eye(4)).to_filename(path)
    return str(path)


@pytest.mark.parametrize('m0tr', [3.0, 6.0])
def test_computecbf_uses_calibration(tmp_path, m0tr):
    """ComputeCBF quantifies with the specification it is given."""
    from aslprep.interfaces.cbf import ComputeCBF

    rng = np.random.default_rng(0)
    shape = (5, 5, 5)
    deltam = _write(rng.uniform(5, 15, (*shape, 2)), tmp_path / 'deltam.nii.gz')
    m0_file = _write(rng.uniform(800, 1200, shape), tmp_path / 'm0.nii.gz')
    metadata = {**PCASL, 'PostLabelingDelay': [1.8, 1.8]}

    def _mean_cbf(name, **inputs):
        work_dir = tmp_path / name
        work_dir.mkdir()
        result = ComputeCBF(
            deltam=deltam, m0_file=m0_file, metadata=metadata, cbf_only=False, **inputs
        ).run(cwd=str(work_dir))
        return nb.load(result.outputs.mean_cbf).get_fdata()

    # The interface resolves the same specification itself when it is not given one.
    calibration = resolve_calibration(metadata, m0tr=m0tr, m0_scale=2.0)
    given = _mean_cbf('given', calibration=calibration)
    resolved = _mean_cbf('resolved', m0_scale=2.0, m0tr=m0tr)
    np.testing.assert_allclose(given, resolved)

    # Halving the labeling efficiency in the specification doubles CBF.
    changed = json.loads(json.dumps(calibration))
    changed['Parameters']['EffectiveLabelingEfficiency']['Value'] /= 2
    np.testing.assert_allclose(_mean_cbf('changed', calibration=changed), 2 * given)


MULTI_DELAY = {
    **PCASL,
    'PostLabelingDelay': [0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
    'LabelingDuration': 1.8,
}


@pytest.mark.parametrize('field_strength', [1.5, 3, 7])
def test_resolve_calibration_tissue_t1(field_strength):
    """The tissue-T1 model adds its gray and white matter T1s and nominal CBFs."""
    metadata = {**MULTI_DELAY, 'MagneticFieldStrength': field_strength}
    consensus = resolve_calibration(metadata, m0tr=6.0)
    assert calibration_kinetic_model(consensus) == 'consensus'
    assert tissue_t1_kinetics(consensus) is None
    for name in TISSUE_T1_PARAMETERS:
        assert name not in consensus['Parameters']

    calibration = resolve_calibration(metadata, m0tr=6.0, kinetic_model='tissue-t1')
    assert calibration_kinetic_model(calibration) == 'tissue-t1'
    t1_gm = T1TISSUE_BY_FIELD_STRENGTH[field_strength]
    t1_wm = T1WM_BY_FIELD_STRENGTH[field_strength]
    assert calibration_value(calibration, 'GrayMatterT1') == t1_gm
    assert calibration_value(calibration, 'WhiteMatterT1') == t1_wm
    for name in TISSUE_T1_PARAMETERS:
        assert _source(calibration, name) == 'default'

    kinetics = tissue_t1_kinetics(calibration)
    np.testing.assert_allclose(kinetics['nominal_cbf'], [60, 20])
    rates = 1 / np.array([t1_gm, t1_wm]) + np.array([60, 20]) / 6000 / 0.9
    np.testing.assert_allclose(kinetics['t1prime'], 1 / rates)


def test_tissue_t1_kinetics_m0_estimate():
    """An M0Estimate has no partition coefficient, but the apparent T1 still uses lambda."""
    metadata = {**MULTI_DELAY, 'M0Type': 'Estimate', 'M0Estimate': 100.0}
    estimate = resolve_calibration(metadata, kinetic_model='tissue-t1')
    separate = resolve_calibration(MULTI_DELAY, m0tr=6.0, kinetic_model='tissue-t1')
    assert 'PartitionCoefficient' not in estimate['Parameters']
    np.testing.assert_allclose(
        tissue_t1_kinetics(estimate)['t1prime'],
        tissue_t1_kinetics(separate)['t1prime'],
    )


def test_resolve_calibration_tissue_t1_refusals():
    """The tissue-T1 model is refused for single-delay data without an ATT map, and unknown models."""
    with pytest.raises(ValueError, match='needs an arterial transit time map'):
        resolve_calibration(PCASL, m0tr=6.0, kinetic_model='tissue-t1')

    # With an ATT map from a multi-delay run, single-delay data can use it
    calibration = resolve_calibration(
        PCASL,
        m0tr=6.0,
        kinetic_model='tissue-t1',
        external_att=True,
    )
    assert 'GrayMatterT1' in calibration['Parameters']

    with pytest.raises(ValueError, match='Unknown kinetic model'):
        resolve_calibration(MULTI_DELAY, m0tr=6.0, kinetic_model='buxton')

    with pytest.raises(ValueError, match=r'cannot be inferred for 2\.89T data'):
        resolve_calibration(
            {**MULTI_DELAY, 'MagneticFieldStrength': 2.89},
            m0tr=6.0,
            kinetic_model='tissue-t1',
        )


def test_describe_model_tissue_t1():
    """The sidecars name the tissue-T1 model, and BASIL's leave out its parameters."""
    calibration = resolve_calibration(MULTI_DELAY, m0tr=6.0, basil=True, kinetic_model='tissue-t1')
    assert 'GrayMatterT1' in calibration['Parameters']['BASILTissueT1']['Description']
    consensus_basil = resolve_calibration(MULTI_DELAY, m0tr=6.0, basil=True)
    assert 'T1 of blood' in consensus_basil['Parameters']['BASILTissueT1']['Description']
    models = {
        model: describe_model(calibration, model, asl_type='PCASL', multi_delay=True)
        for model in MODELS
    }
    for model in ('aslprep', 'lr'):
        assert 'apparent T1 of gray and white matter' in models[model]['Description']
        for name in TISSUE_T1_PARAMETERS:
            assert name in models[model]['Parameters']

    assert 'apparent T1' not in models['basil']['Description']
    for name in TISSUE_T1_PARAMETERS:
        assert name not in models['basil']['Parameters']

    consensus = resolve_calibration(MULTI_DELAY, m0tr=6.0)
    description = describe_model(consensus, 'aslprep', asl_type='PCASL', multi_delay=True)
    assert 'apparent T1' not in description['Description']


def test_computecbf_tissue_t1(tmp_path):
    """ComputeCBF fits the tissue-T1 model with the partial volume maps it is given."""
    from aslprep.interfaces.cbf import ComputeCBF
    from aslprep.tests.truth_models import tissue_t1_signal
    from aslprep.utils.cbf import tissue_mixing_weights

    shape = (3, 3, 2)
    n_voxels = int(np.prod(shape))
    plds = np.array(MULTI_DELAY['PostLabelingDelay'])
    rng = np.random.default_rng(4)
    gm_pv = rng.uniform(0, 1, shape)
    wm_pv = (1 - gm_pv) * rng.uniform(0, 1, shape)
    calibration = resolve_calibration(MULTI_DELAY, m0tr=6.0, kinetic_model='tissue-t1')
    kinetics = tissue_t1_kinetics(calibration)
    weights = tissue_mixing_weights(gm_pv.ravel(), wm_pv.ravel(), *kinetics['nominal_cbf'])
    labeleff = calibration_value(calibration, 'EffectiveLabelingEfficiency')
    t1blood = calibration_value(calibration, 'BloodT1')
    true_cbf, true_att, m0 = 50.0, 1.1, 1000.0
    signal = np.zeros((n_voxels, plds.size))
    for i_tissue, t1prime in enumerate(kinetics['t1prime']):
        tissue = [tissue_t1_signal(t, true_att, 1.8, t1blood, t1prime, True) for t in plds + 1.8]
        signal += weights[:, i_tissue, None] * np.array(tissue)
    deltam = 2 * labeleff * (m0 / 0.9) * true_cbf / 6000 * signal
    deltam_file = _write(deltam.reshape(*shape, plds.size), tmp_path / 'deltam.nii.gz')
    m0_file = _write(np.full(shape, m0), tmp_path / 'm0.nii.gz')
    gm_file = _write(gm_pv, tmp_path / 'gm.nii.gz')
    wm_file = _write(wm_pv, tmp_path / 'wm.nii.gz')
    inputs = {
        'deltam': deltam_file,
        'm0_file': m0_file,
        'metadata': MULTI_DELAY,
        'cbf_only': False,
        'calibration': calibration,
    }

    work_dir = tmp_path / 'fit'
    work_dir.mkdir()
    result = ComputeCBF(**inputs, gm_pv=gm_file, wm_pv=wm_file).run(cwd=str(work_dir))
    np.testing.assert_allclose(nb.load(result.outputs.mean_cbf).get_fdata(), true_cbf, rtol=1e-3)
    np.testing.assert_allclose(nb.load(result.outputs.att).get_fdata(), true_att, atol=1e-3)

    with pytest.raises(ValueError, match='needs gray and white matter partial volume maps'):
        ComputeCBF(**inputs).run(cwd=str(work_dir))

    wrong_shape = _write(np.zeros((3, 3, 3)), tmp_path / 'wrong_shape.nii.gz')
    wrong_affine = tmp_path / 'wrong_affine.nii.gz'
    nb.Nifti1Image(gm_pv.astype(np.float32), np.diag([-1, 1, 1, 1])).to_filename(wrong_affine)
    wrong_ndim = _write(gm_pv[..., None], tmp_path / 'wrong_ndim.nii.gz')
    for wrong in (wrong_shape, str(wrong_affine), wrong_ndim):
        with pytest.raises(ValueError, match='does not match the delta-M grid'):
            ComputeCBF(**inputs, gm_pv=wrong, wm_pv=wm_file).run(cwd=str(work_dir))


def test_update_calibration():
    """Run-time values replace whole parameters, and bad entries are refused."""
    from aslprep.utils.calibration import update_calibration

    calibration = resolve_calibration(MULTI_DELAY, m0tr=6.0, kinetic_model='tissue-t1')
    entry = {
        'Value': 1.4,
        'Units': 's',
        'Source': 'metadata',
        'Origin': 'Median of a T1 map',
        'Description': 'T1 of gray matter.',
    }
    updated = update_calibration(calibration, {'GrayMatterT1': entry})
    assert updated['Parameters']['GrayMatterT1'] == entry
    # The original is unchanged, and the other parameters are kept
    assert calibration_value(calibration, 'GrayMatterT1') == T1TISSUE_BY_FIELD_STRENGTH[3]
    assert updated['Parameters']['BloodT1'] == calibration['Parameters']['BloodT1']
    assert update_calibration(calibration, {}) == calibration

    with pytest.raises(ValueError, match='must have exactly the keys'):
        update_calibration(calibration, {'GrayMatterT1': {'Value': 1.4}})

    with pytest.raises(ValueError, match="Unknown source 'scanner'"):
        update_calibration(calibration, {'GrayMatterT1': {**entry, 'Source': 'scanner'}})

    # A misspelled or new parameter would be recorded but never used
    with pytest.raises(ValueError, match='GreyMatterT1 is not a parameter'):
        update_calibration(calibration, {'GreyMatterT1': entry})

    consensus = resolve_calibration(MULTI_DELAY, m0tr=6.0)
    with pytest.raises(ValueError, match='GrayMatterT1 is not a parameter'):
        update_calibration(consensus, {'GrayMatterT1': entry})


def test_resolve_runtime_calibration(tmp_path):
    """The run-time node gives the updated specification and each model's Model metadata."""
    from aslprep.interfaces.cbf import ResolveRuntimeCalibration

    calibration = resolve_calibration(MULTI_DELAY, m0tr=6.0, basil=True, kinetic_model='tissue-t1')
    inputs = {'calibration': calibration, 'asl_type': 'PCASL', 'multi_delay': True}
    result = ResolveRuntimeCalibration(**inputs).run(cwd=str(tmp_path))
    assert result.outputs.calibration == calibration
    assert result.outputs.model_metadata == {
        model: describe_model(calibration, model, asl_type='PCASL', multi_delay=True)
        for model in MODELS
    }

    updates = {
        'GrayMatterT1': {
            **calibration['Parameters']['GrayMatterT1'],
            'Value': 1.4,
            'Source': 'metadata',
        },
        'BloodT1': {**calibration['Parameters']['BloodT1'], 'Value': 1.7},
    }
    result = ResolveRuntimeCalibration(**inputs, updates=updates).run(cwd=str(tmp_path))
    models = result.outputs.model_metadata
    assert calibration_value(result.outputs.calibration, 'GrayMatterT1') == 1.4
    for model in ('aslprep', 'lr'):
        assert models[model]['Parameters']['GrayMatterT1'] == 1.4
        assert models[model]['Parameters']['BloodT1'] == 1.7

    # BASIL was configured from the specification as built, so its sidecars describe that one
    assert 'GrayMatterT1' not in models['basil']['Parameters']
    assert models['basil']['Parameters']['BloodT1'] == calibration_value(calibration, 'BloodT1')


@pytest.mark.parametrize(
    ('name', 'model', 'base', 'base_from_file', 'map_node'),
    [
        ('mean_cbf', 'aslprep', None, False, False),
        ('mean_cbf', 'aslprep', {'Other': 1}, False, False),
        # A base that already has a Model: the model's own wins, as the fixed argument did
        ('mean_cbf', 'aslprep', {'Model': 'stale', 'Other': 1}, False, False),
        ('mean_cbf_basil', 'basil', None, False, False),
        ('mean_cbf_gm_lr', 'lr', {'Other': 1}, False, False),
        # A base read from a file through a connect function, as for the CIFTI outputs
        ('mean_cbf', 'aslprep', {'Other': 1}, True, False),
        # MapNode sinks, as for the surface outputs
        ('mean_cbf', 'aslprep', None, False, True),
    ],
)
def test_model_metadata_reaches_sidecar(tmp_path, name, model, base, base_from_file, map_node):
    """A sink given the Model metadata at run time writes the sidecar a fixed argument did."""
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe

    from aslprep.interfaces.bids import DerivativesDataSink
    from aslprep.workflows.asl.outputs import _read_json, connect_model_metadata, output_fields

    calibration = resolve_calibration(MULTI_DELAY, m0tr=6.0, basil=True)
    models = {
        each: describe_model(calibration, each, asl_type='PCASL', multi_delay=True)
        for each in MODELS
    }
    in_file = _write(np.ones((2, 2, 2)), tmp_path / 'cbf.nii.gz')
    source_file = tmp_path / 'sub-01_asl.nii.gz'
    source_file.touch()
    base_file = tmp_path / 'base.json'
    base_file.write_text(json.dumps(base))
    base_field = ('base_file', _read_json) if base_from_file else 'base'

    def sidecars(out_dir, run_time):
        workflow = pe.Workflow(name=f'wf_{run_time}', base_dir=str(tmp_path / 'work'))
        inputnode = pe.Node(
            niu.IdentityInterface(fields=['model_metadata', 'base', 'base_file']),
            name='inputnode',
        )
        inputnode.inputs.model_metadata = models
        inputnode.inputs.base = base
        inputnode.inputs.base_file = str(base_file)
        sink_kwargs = {
            'source_file': str(source_file),
            'base_directory': str(out_dir),
            **output_fields(name, None if run_time else models),
        }
        if map_node:
            # Two outputs, told apart by their acq entity, as surface sinks are by hemisphere
            sink = pe.MapNode(
                DerivativesDataSink(**sink_kwargs),
                iterfield=['in_file', 'acquisition'],
                name=f'ds_{name}',
            )
            sink.inputs.in_file = [in_file, in_file]
            sink.inputs.acquisition = ['a', 'b']
        else:
            sink = pe.Node(DerivativesDataSink(in_file=in_file, **sink_kwargs), name=f'ds_{name}')

        if run_time:
            pair = None if base is None else (inputnode, base_field)
            connect_model_metadata(workflow, inputnode, sink, name, base=pair)
        elif base is not None:
            workflow.connect([(inputnode, sink, [(base_field, 'meta_dict')])])
        else:
            workflow.add_nodes([sink])
        workflow.run()
        return [json.loads(path.read_text()) for path in sorted(Path(out_dir).rglob('*.json'))]

    static = sidecars(tmp_path / 'static', run_time=False)
    run_time = sidecars(tmp_path / 'run_time', run_time=True)
    assert len(run_time) == (2 if map_node else 1)
    assert run_time == static
    for sidecar in run_time:
        assert sidecar['Model'] == models[model]
        if base:
            assert sidecar['Other'] == 1


def test_computecbf_calibration_xor():
    """A specification and the values it contains cannot both be given."""
    from aslprep.interfaces.cbf import ComputeCBF

    calibration = resolve_calibration(PCASL, m0tr=6.0)
    with pytest.raises(OSError, match='mutually exclusive'):
        ComputeCBF(calibration=calibration, m0_scale=1.0)
