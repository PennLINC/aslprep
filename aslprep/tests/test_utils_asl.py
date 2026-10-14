"""Tests for the aslprep.utils.asl module."""

import numpy as np
import pytest

from aslprep import config
from aslprep.utils.asl import (
    estimate_labeling_efficiency,
    ignore_volume_types,
    infer_m0tr,
    prepare_basil_kwargs,
    select_processing_target,
)


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


def _write_aslcontext(tmp_path, volume_types):
    aslcontext = tmp_path / 'aslcontext.tsv'
    aslcontext.write_text('\n'.join(['volume_type', *volume_types]) + '\n')
    return str(aslcontext)


@pytest.mark.parametrize(
    ('volume_types', 'metadata', 'm0scan_metadata', 'expected'),
    [
        (
            ['control', 'label'],
            {'M0Type': 'Separate', 'RepetitionTimePreparation': 4.0},
            {'RepetitionTimePreparation': 6.0},
            6.0,
        ),
        (
            ['m0scan', 'control', 'label'],
            {'M0Type': 'Included', 'RepetitionTimePreparation': [8.0, 4.0, 4.0]},
            None,
            [8.0],
        ),
        (
            ['control', 'label'],
            {
                'M0Type': 'Absent',
                'BackgroundSuppression': False,
                'RepetitionTimePreparation': 4.0,
            },
            None,
            4.0,
        ),
        (
            ['control', 'label'],
            {'M0Type': 'Estimate', 'M0Estimate': 100, 'RepetitionTimePreparation': 0},
            None,
            None,
        ),
    ],
)
def test_infer_m0tr(tmp_path, volume_types, metadata, m0scan_metadata, expected):
    """infer_m0tr returns the TR of the volumes used for calibration."""
    m0tr = infer_m0tr(
        aslcontext=_write_aslcontext(tmp_path, volume_types),
        metadata=metadata,
        m0scan_metadata=m0scan_metadata,
    )
    if expected is None:
        assert m0tr is None
    else:
        assert np.array_equal(m0tr, expected)


@pytest.mark.parametrize(
    ('volume_types', 'metadata', 'm0scan_metadata'),
    [
        (
            ['control', 'label'],
            {'M0Type': 'Separate', 'RepetitionTimePreparation': 4.0},
            {'RepetitionTimePreparation': 0},
        ),
        # The ds000240 case from gh-654.
        (
            ['m0scan', 'control', 'label'],
            {'M0Type': 'Included', 'RepetitionTimePreparation': 0},
            None,
        ),
        (
            ['m0scan', 'control', 'label'],
            {'M0Type': 'Included', 'RepetitionTimePreparation': [0, 4.0, 4.0]},
            None,
        ),
        (
            ['control', 'label'],
            {
                'M0Type': 'Absent',
                'BackgroundSuppression': False,
                'RepetitionTimePreparation': 0,
            },
            None,
        ),
    ],
)
def test_infer_m0tr_nonpositive(tmp_path, volume_types, metadata, m0scan_metadata):
    """A zero M0 TR raises an error instead of silently zeroing CBF."""
    with pytest.raises(ValueError, match=r'RepetitionTimePreparation .* must be positive'):
        infer_m0tr(
            aslcontext=_write_aslcontext(tmp_path, volume_types),
            metadata=metadata,
            m0scan_metadata=m0scan_metadata,
        )


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


def _write_volume_types(tmp_path, volume_types):
    aslcontext = tmp_path / 'sub-01_aslcontext.tsv'
    aslcontext.write_text('\n'.join(['volume_type', *volume_types]) + '\n')
    return str(aslcontext)


@pytest.mark.parametrize(
    ('volume_types', 'ignore', 'target', 'reference'),
    [
        # Scanner-calculated CBF next to control-label pairs: the reference is no longer CBF.
        (['m0scan', 'control', 'label', 'cbf'], ['cbf'], 'control', 'm0scan'),
        (['control', 'label', 'cbf'], ['cbf'], 'control', 'control'),
        (['control', 'label', 'deltam', 'cbf'], ['cbf'], 'control', 'deltam'),
        (['control', 'label', 'deltam', 'cbf'], ['cbf', 'deltam'], 'control', 'control'),
        # GE-style series: ignoring one precomputed type falls back on the other.
        (['m0scan', 'deltam', 'cbf'], ['deltam'], 'cbf', 'cbf'),
        (['m0scan', 'deltam', 'cbf'], ['cbf'], 'deltam', 'deltam'),
    ],
)
def test_ignore_volume_types(tmp_path, volume_types, ignore, target, reference):
    """Ignored volume types are not processed and are not used for the reference image."""
    import pandas as pd

    from aslprep.workflows.asl.fit import select_reference_volume_type

    aslcontext = _write_volume_types(tmp_path, volume_types)
    out_file = ignore_volume_types(aslcontext, ignore=ignore, out_dir=tmp_path / 'work')

    assert out_file != aslcontext
    out_types = pd.read_table(out_file)['volume_type'].tolist()
    # The number of rows still matches the number of volumes in the ASL file.
    assert len(out_types) == len(volume_types)
    assert not set(out_types) & set(ignore)
    assert [t for t in out_types if t != 'ignored'] == [t for t in volume_types if t not in ignore]

    assert select_processing_target(out_file) == target
    assert (
        select_reference_volume_type(
            aslcontext=out_file,
            metadata={'M0Type': 'Included' if 'm0scan' in volume_types else 'Absent'},
            prioritize_m0=False,
        )
        == reference
    )


@pytest.mark.parametrize(
    'ignore',
    [[], ['cbf'], ['fieldmaps', 'sbref'], None],
)
def test_ignore_volume_types_unchanged(tmp_path, ignore):
    """The aslcontext file is used as is when it has no volumes to ignore."""
    aslcontext = _write_volume_types(tmp_path, ['m0scan', 'control', 'label'])

    assert ignore_volume_types(aslcontext, ignore=ignore, out_dir=tmp_path / 'work') == aslcontext
    assert not (tmp_path / 'work').exists()


@pytest.mark.parametrize(
    ('volume_types', 'ignore'),
    [
        (['m0scan', 'cbf'], ['cbf']),
        (['m0scan', 'deltam'], ['deltam']),
        (['deltam', 'cbf'], ['cbf', 'deltam']),
        (['control', 'cbf'], ['cbf']),
    ],
)
def test_ignore_volume_types_nothing_left(tmp_path, volume_types, ignore):
    """Ignoring every volume type that can be processed is an error."""
    aslcontext = _write_volume_types(tmp_path, volume_types)

    with pytest.raises(ValueError, match=r'\(--ignore\) leaves no'):
        ignore_volume_types(aslcontext, ignore=ignore, out_dir=tmp_path / 'work')
