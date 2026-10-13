"""Tests for the aslprep.utils.asl module."""

import pytest

from aslprep import config
from aslprep.utils.asl import estimate_labeling_efficiency, prepare_basil_kwargs


def ascending(n_slices, spacing=0.09375):
    return [i * spacing for i in range(n_slices)]


@pytest.mark.parametrize(
    ('desc', 'metadata', 'expected'),
    [
        (
            'no SliceTiming key',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'MagneticFieldStrength': 3,
            },
            {},
        ),
        (
            'empty SliceTiming list',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'MagneticFieldStrength': 3,
                'SliceTiming': [],
            },
            {},
        ),
        (
            'single-band ascending (32 slices, TR=3)',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'MagneticFieldStrength': 3,
                'SliceTiming': ascending(32),  # default spacing = 0.09375 = 3/32
            },
            {'slice_spacing': 0.09375},
        ),
        (
            'single-band ascending with different spacing (24 slices, TR=3)',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'MagneticFieldStrength': 3,
                'SliceTiming': ascending(24, spacing=0.125),  # 3/24
            },
            {'slice_spacing': 0.125},
        ),
        (
            'single-band non-monotonic',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'MagneticFieldStrength': 3,
                # duplicate slice-time at index 10
                'SliceTiming': (ascending(32)[:10] + [ascending(32)[9]] + ascending(32)[11:]),
            },
            {},
        ),
        (
            'single-band descending',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'MagneticFieldStrength': 3,
                'SliceTiming': list(reversed(ascending(32))),
            },
            {},
        ),
        (
            'single-band floating-point precision handled',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'MagneticFieldStrength': 3,
                # tiny jitter within rounding tolerance
                'SliceTiming': (
                    [0.0, 0.0937500000001, 0.1875000000002] + [i * 0.09375 for i in range(3, 32)]
                ),
            },
            {'slice_spacing': 0.09375},
        ),
        (
            'multiband factor 2 ascending (16 slices/band, TR=3)',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'MagneticFieldStrength': 3,
                'MultibandAccelerationFactor': 2,
                # spacing = 3.0 / (32 / 2) = 3.0/16
                'SliceTiming': ascending(16, spacing=3.0 / 16) * 2,
            },
            {'sliceband': 16, 'slice_spacing': 3.0 / 16},
        ),
        (
            'multiband factor 2 ascending (16 slices/band, TR=4)',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'MagneticFieldStrength': 3,
                'MultibandAccelerationFactor': 2,
                # spacing = 4.0 / (32 / 2) = 4.0/16
                'SliceTiming': ascending(16, spacing=4.0 / 16) * 2,
            },
            {'sliceband': 16, 'slice_spacing': 4.0 / 16},
        ),
        (
            'multiband factor 2 non-monotonic in first band',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'MagneticFieldStrength': 3,
                'MultibandAccelerationFactor': 2,
                # break monotonicity in the first band only
                'SliceTiming': (
                    ascending(16, spacing=3.0 / 16)[:5]
                    + [ascending(16, spacing=3.0 / 16)[4]]
                    + ascending(16, spacing=3.0 / 16)[6:]
                )
                + ascending(16, spacing=3.0 / 16),
            },
            {},
        ),
        (
            'multiband factor 2 second band descending',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'MagneticFieldStrength': 3,
                'MultibandAccelerationFactor': 2,
                # first band ascending, second band reversed
                'SliceTiming': (
                    ascending(16, spacing=3.0 / 16)
                    + list(reversed(ascending(16, spacing=3.0 / 16)))
                ),
            },
            {},
        ),
        (
            'multiband factor 3 ascending (10 slices/band, TR=3)',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'MagneticFieldStrength': 3,
                'MultibandAccelerationFactor': 3,
                # spacing = 3.0 / (30 / 3) = 3.0/10
                'SliceTiming': ascending(10, spacing=3.0 / 10) * 3,
            },
            {'sliceband': 10, 'slice_spacing': 3.0 / 10},
        ),
    ],
)
def test_prepare_basil_kwargs(desc, metadata, expected):
    """Parametrized tests for prepare_basil_kwargs."""
    result = prepare_basil_kwargs(metadata)
    assert result == expected, f'Failed case: {desc}'


@pytest.mark.parametrize(
    ('desc', 'metadata', 'expected'),
    [
        (
            'default PCASL, no suppression',
            {'ArterialSpinLabelingType': 'PCASL', 'BackgroundSuppression': False},
            0.85,
        ),
        (
            'default PASL, suppression not stated',
            {'ArterialSpinLabelingType': 'PASL'},
            0.98,
        ),
        (
            'default CASL with two suppression pulses',
            {
                'ArterialSpinLabelingType': 'CASL',
                'BackgroundSuppression': True,
                'BackgroundSuppressionNumberPulses': 2,
            },
            0.68 * 0.95**2,
        ),
        (
            'sidecar efficiency, no suppression',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'LabelingEfficiency': 0.7,
                'BackgroundSuppression': False,
            },
            0.7,
        ),
        (
            'sidecar efficiency with four suppression pulses',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'LabelingEfficiency': 0.7,
                'BackgroundSuppression': True,
                'BackgroundSuppressionNumberPulses': 4,
            },
            0.7 * 0.95**4,
        ),
        (
            'pulse count ignored without suppression',
            {
                'ArterialSpinLabelingType': 'PCASL',
                'LabelingEfficiency': 0.7,
                'BackgroundSuppression': False,
                'BackgroundSuppressionNumberPulses': 4,
            },
            0.7,
        ),
    ],
)
def test_estimate_labeling_efficiency(desc, metadata, expected):
    """The suppression loss applies whether or not the sidecar gives LabelingEfficiency."""
    assert estimate_labeling_efficiency(metadata) == pytest.approx(expected), desc


@pytest.mark.parametrize('sidecar', [{}, {'LabelingEfficiency': 0.7}])
def test_estimate_labeling_efficiency_missing_pulse_count(sidecar, monkeypatch):
    """Without BackgroundSuppressionNumberPulses, one pulse is assumed, with a warning."""
    messages = []
    monkeypatch.setattr(config.loggers.utils, 'warning', messages.append)
    metadata = {'ArterialSpinLabelingType': 'PCASL', 'BackgroundSuppression': True, **sidecar}

    labeleff = estimate_labeling_efficiency(metadata)

    assert labeleff == pytest.approx(sidecar.get('LabelingEfficiency', 0.85) * 0.95)
    assert len(messages) == 1
    assert 'BackgroundSuppressionNumberPulses is not defined' in messages[0]
