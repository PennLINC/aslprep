"""Tests for the aslprep.utils.asl module."""

import numpy as np
import pytest

from aslprep import config
from aslprep.utils.asl import (
    count_reduced_volumes,
    describe_asl_acquisition,
    describe_labeling_efficiency,
    estimate_labeling_efficiency,
    ignore_volume_types,
    infer_m0tr,
    prepare_basil_kwargs,
    select_processing_target,
)


def ascending(n_slices, spacing=0.09375):
    return [i * spacing for i in range(n_slices)]


def rounded(slice_times, step):
    return [round(t / step) * step for t in slice_times]


def interleaved(n_slices, spacing):
    order = list(range(0, n_slices, 2)) + list(range(1, n_slices, 2))
    slice_times = [0.0] * n_slices
    for i_acq, i_slice in enumerate(order):
        slice_times[i_slice] = i_acq * spacing
    return slice_times


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
        (
            'all slices at time zero',
            {'SliceTiming': [0.0] * 20},
            {},
        ),
        (
            'ascending times rounded to 2.5 ms, so steps alternate between 42.5 and 45 ms',
            {'SliceTiming': rounded(ascending(20, spacing=0.04375), 0.0025)},
            {'slice_spacing': 0.04375},
        ),
        (
            'ascending times rounded to 2.5 ms, multiband factor 2',
            {
                'MultibandAccelerationFactor': 2,
                'SliceTiming': rounded(ascending(10, spacing=0.04375), 0.0025) * 2,
            },
            {'sliceband': 10, 'slice_spacing': 0.04375},
        ),
        (
            'descending times, with slices encoded along k-',
            {
                'SliceEncodingDirection': 'k-',
                'SliceTiming': list(reversed(ascending(32))),
            },
            {'slice_spacing': 0.09375},
        ),
    ],
)
def test_prepare_basil_kwargs(desc, metadata, expected):
    """BASIL gets --slicedt (and --sliceband) when it can describe the slice times."""
    basil_kwargs, use_tiimg = prepare_basil_kwargs(metadata)
    assert basil_kwargs == pytest.approx(expected, abs=1e-4), f'Failed case: {desc}'
    assert not use_tiimg, f'Failed case: {desc}'


@pytest.mark.parametrize(
    ('desc', 'metadata'),
    [
        (
            'single-band non-monotonic',
            # duplicate slice time at index 10
            {'SliceTiming': ascending(32)[:10] + [ascending(32)[9]] + ascending(32)[11:]},
        ),
        (
            'multiband factor 2, non-monotonic in the first band',
            {
                'MultibandAccelerationFactor': 2,
                'SliceTiming': (
                    ascending(16)[:5] + [ascending(16)[4]] + ascending(16)[5:15] + ascending(16)
                ),
            },
        ),
        ('single-band descending', {'SliceTiming': list(reversed(ascending(32)))}),
        (
            'interleaved, even slices first',
            {'SliceTiming': interleaved(20, spacing=0.045)},
        ),
        (
            'ascending, but the first slice is not at time zero',
            {'SliceTiming': [0.01 + t for t in ascending(20, spacing=0.045)]},
        ),
        (
            'ascending, but slices are encoded along j',
            {'SliceEncodingDirection': 'j', 'SliceTiming': ascending(20, spacing=0.045)},
        ),
        (
            'multiband factor 2, second band descending',
            {
                'MultibandAccelerationFactor': 2,
                'SliceTiming': ascending(16) + list(reversed(ascending(16))),
            },
        ),
        (
            'multiband factor 2, bands acquired at different times',
            {
                'MultibandAccelerationFactor': 2,
                'SliceTiming': ascending(8, spacing=0.05) + ascending(8, spacing=0.06),
            },
        ),
        (
            'multiband factor that does not divide the slices',
            {'MultibandAccelerationFactor': 3, 'SliceTiming': ascending(16)},
        ),
    ],
)
def test_prepare_basil_kwargs_tiimg(desc, metadata):
    """BASIL needs a voxelwise inflow time image when --slicedt cannot describe the slices."""
    basil_kwargs, use_tiimg = prepare_basil_kwargs(metadata)
    assert basil_kwargs == {}, f'Failed case: {desc}'
    assert use_tiimg, f'Failed case: {desc}'


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


def test_estimate_labeling_efficiency_no_warning(monkeypatch):
    """The boilerplate can ask for the same value without repeating the warning."""
    messages = []
    monkeypatch.setattr(config.loggers.utils, 'warning', messages.append)
    metadata = {'ArterialSpinLabelingType': 'PCASL', 'BackgroundSuppression': True}

    assert estimate_labeling_efficiency(metadata, warn=False) == pytest.approx(0.85 * 0.95)
    assert not messages


@pytest.mark.parametrize(
    ('metadata', 'expected'),
    [
        (
            {'ArterialSpinLabelingType': 'PCASL', 'BackgroundSuppression': False},
            ['a labeling efficiency of 0.85', 'typical labeling efficiency for PCASL'],
        ),
        (
            {'ArterialSpinLabelingType': 'CASL'},
            ['a labeling efficiency of 0.68', '@wang2005amplitude'],
        ),
        (
            {
                'ArterialSpinLabelingType': 'PCASL',
                'LabelingEfficiency': 0.8,
                'BackgroundSuppression': True,
                'BackgroundSuppressionNumberPulses': 2,
            },
            [
                'a labeling efficiency of 0.722',
                'labeling efficiency in the BIDS metadata',
                'each of the 2 background suppression pulses',
            ],
        ),
        (
            {'ArterialSpinLabelingType': 'PASL', 'BackgroundSuppression': True},
            ['a labeling efficiency of 0.931', 'which was assumed'],
        ),
    ],
)
def test_describe_labeling_efficiency(metadata, expected):
    """The description reports the value that estimate_labeling_efficiency returns."""
    description = describe_labeling_efficiency(metadata)
    for text in expected:
        assert text in description


def test_describe_asl_acquisition_pcasl():
    """Multi-delay PCASL reports the OSIPI parameters in the sidecar, ignoring M0 volumes."""
    import pandas as pd

    volume_types = ['m0scan', 'control', 'label', 'control', 'label', 'control', 'label']
    metadata = {
        'ArterialSpinLabelingType': 'PCASL',
        'PCASLType': 'unbalanced',
        'BackgroundSuppression': True,
        'BackgroundSuppressionNumberPulses': 2,
        'BackgroundSuppressionPulseTime': [1.2, 2.1],
        'LabelingDuration': 1.8,
        'PostLabelingDelay': [0, 1.0, 1.0, 1.5, 1.5, 2.0, 2.0],
        'AcquisitionVoxelSize': [3.5, 3.5, 4],
        'VascularCrushing': False,
        'LabelingDistance': -90,
        'LabelingPulseAverageGradient': 0.6,
        'LabelingPulseDuration': 0.5,
    }
    description = describe_asl_acquisition(metadata, pd.DataFrame({'volume_type': volume_types}))

    assert description.startswith(
        'The ASL data were acquired with an unbalanced pseudo-continuous ASL (PCASL) sequence, '
        'with background suppression (2 pulses).'
    )
    assert 'applied 1.2 and 2.1 s after the start of labeling' in description
    assert (
        'The labeling duration was 1.8 s and the post-labeling delays were 1, 1.5, and 2 s.'
    ) in description
    assert 'The ASL file contained 3 control-label pairs.' in description
    assert 'The acquisition voxel size was 3.5 x 3.5 x 4 mm.' in description
    assert 'No vascular crushing was applied.' in description
    assert 'center of the labeling plane was 90 mm from the center' in description
    assert (
        'The labeling pulses had an average gradient of 0.6 mT/m and a duration of 0.5 ms.'
    ) in description


def test_describe_asl_acquisition_pasl(tmp_path):
    """PASL reports inversion times and the bolus cut-off, and nothing it does not know."""
    aslcontext = _write_volume_types(tmp_path, ['label', 'control'] * 4)
    metadata = {
        'ArterialSpinLabelingType': 'PASL',
        'PASLType': 'PICORE',
        'PostLabelingDelay': 1.8,
        'BolusCutOffFlag': True,
        'BolusCutOffTechnique': 'Q2TIPS',
        'BolusCutOffDelayTime': [0.7, 1.6],
        'TotalAcquiredPairs': 40,
        'VascularCrushing': True,
        'VascularCrushingVENC': 4,
        'LabelingSlabThickness': 100,
        'LabelingPulseDuration': 10,
    }
    description = describe_asl_acquisition(metadata, aslcontext)

    assert 'with a pulsed ASL (PASL, PICORE) sequence.' in description
    assert 'background suppression' not in description
    assert 'The inversion time was 1.8 s.' in description
    assert 'post-labeling delay' not in description
    assert (
        'The Q2TIPS bolus cut-off technique was used [@luh1999quipss], with a bolus cut-off '
        'delay time (TI1) of 0.7 s and a last saturation pulse (TI1s) at 1.6 s.'
    ) in description
    assert 'A total of 40 control-label pairs were acquired.' in description
    assert 'The ASL file contained' not in description
    assert 'Vascular crushing was applied (VENC = 4 cm/s).' in description
    assert 'The labeling slab was 100 mm thick.' in description
    # Pulse parameters only apply to (P)CASL
    assert 'labeling pulses' not in description
    assert 'voxel size' not in description


def test_describe_asl_acquisition_scanner_volumes():
    """Files without control-label pairs report the scanner-calculated volumes instead."""
    import pandas as pd

    metadata = {
        'ArterialSpinLabelingType': 'PCASL',
        'BackgroundSuppression': True,
        'BolusCutOffFlag': False,
        'PostLabelingDelay': 2.025,
    }
    aslcontext = pd.DataFrame({'volume_type': ['deltam', 'm0scan']})
    description = describe_asl_acquisition(metadata, aslcontext)

    assert 'a pseudo-continuous ASL (PCASL) sequence, with background suppression.' in description
    assert 'The post-labeling delay was 2.025 s.' in description
    assert 'The ASL file contained 1 delta-M volume calculated on the scanner.' in description
    assert 'bolus cut-off' not in description


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


@pytest.mark.parametrize(
    ('volume_types', 'n_reduced'),
    [
        # GE product output (#728): the CBF volume is dropped, leaving the M0 and delta-M.
        (['m0scan', 'deltam', 'cbf'], 2),
        (['m0scan', 'cbf'], 2),
        (['m0scan', 'm0scan'] + ['label', 'control'] * 3 + ['deltam', 'cbf'], 8),
        (['deltam'] * 4 + ['cbf'], 4),
    ],
)
def test_count_reduced_volumes(tmp_path, volume_types, n_reduced):
    """The count matches the series that ReduceASLFiles writes."""
    import nibabel as nb

    from aslprep.interfaces.utility import ReduceASLFiles

    aslcontext = _write_volume_types(tmp_path, volume_types)
    target = select_processing_target(aslcontext)
    assert count_reduced_volumes(aslcontext, target) == n_reduced

    asl_file = tmp_path / 'asl.nii.gz'
    data = np.zeros((2, 2, 2, len(volume_types)), dtype=np.float32)
    nb.Nifti1Image(data, np.eye(4)).to_filename(asl_file)
    result = ReduceASLFiles(
        asl_file=str(asl_file),
        aslcontext=aslcontext,
        processing_target=target,
        metadata={},
    ).run(cwd=str(tmp_path))
    assert nb.load(result.outputs.asl_file).shape[3] == n_reduced
