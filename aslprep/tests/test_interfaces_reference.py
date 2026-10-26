"""Tests for aslprep.interfaces.reference."""

import nibabel as nb
import numpy as np
import pytest


def _box(shape, lo, hi):
    """A binary mask with a box from ``lo`` to ``hi`` (exclusive) on every axis."""
    data = np.zeros(shape, dtype=np.uint8)
    data[lo[0] : hi[0], lo[1] : hi[1], lo[2] : hi[2]] = 1
    return data


def _run_guard(tmp_path, anat, asl, **kwargs):
    from aslprep.interfaces.reference import GuardAnatomicalMask

    affine = np.diag([3.0, 3.0, 4.0, 1.0])
    paths = {}
    for name, data in (('anat', anat), ('asl', asl)):
        paths[name] = str(tmp_path / f'{name}_mask.nii.gz')
        nb.Nifti1Image(data, affine).to_filename(paths[name])
    result = GuardAnatomicalMask(anat_mask=paths['anat'], asl_mask=paths['asl'], **kwargs).run(
        cwd=str(tmp_path)
    )
    out = nb.load(result.outputs.out_mask)
    assert np.allclose(out.affine, affine)
    return np.asanyarray(out.dataobj), result.outputs


SHAPE = (20, 20, 20)
ANAT = _box(SHAPE, (4, 4, 4), (16, 16, 16))


@pytest.mark.parametrize(
    ('asl', 'used', 'expected_dice'),
    [
        # The masks agree: the anatomical mask is used.
        (_box(SHAPE, (5, 4, 4), (16, 16, 16)), 'anatomical', 2 * 11 / 23),
        # SynthStrip lost the inferior half of the brain (the smaller mask): anatomical.
        (_box(SHAPE, (4, 4, 10), (16, 16, 16)), 'anatomical', 2 * 6 / 18),
        # The anatomical mask is displaced (failed coregistration): SynthStrip.
        (_box(SHAPE, (8, 4, 4), (20, 16, 16)), 'synthstrip', 2 * 8 / 24),
    ],
)
def test_guard_anatomical_mask(tmp_path, asl, used, expected_dice):
    """The anatomical mask is used unless the masks disagree and SynthStrip's is not smaller."""
    out, outputs = _run_guard(tmp_path, ANAT, asl)

    np.testing.assert_array_equal(out, ANAT if used == 'anatomical' else asl)
    assert outputs.qc['anat_mask_synthstrip_dice'] == pytest.approx(expected_dice)
    assert outputs.qc['anat_mask_synthstrip_volume_ratio'] == pytest.approx(asl.sum() / ANAT.sum())
    assert outputs.qc['anat_mask_used'] == int(used == 'anatomical')
    assert outputs.metadata == {
        'SkullStripASL': 'anatomical',
        'BrainMaskSource': used,
        'AnatomicalSynthStripDice': pytest.approx(expected_dice),
    }


def test_guard_anatomical_mask_thresholds(tmp_path):
    """The thresholds are inputs: a stricter Dice threshold rejects a close anatomical mask."""
    asl = _box(SHAPE, (5, 4, 4), (16, 16, 16))
    out, outputs = _run_guard(tmp_path, ANAT, asl, min_dice=0.99, min_volume_ratio=0.9)
    np.testing.assert_array_equal(out, asl)
    assert outputs.qc['anat_mask_used'] == 0


def test_guard_anatomical_mask_empty(tmp_path):
    """An anatomical mask entirely outside the field of view is never used."""
    asl = _box(SHAPE, (4, 4, 4), (16, 16, 16))
    out, outputs = _run_guard(tmp_path, np.zeros(SHAPE, dtype=np.uint8), asl)
    np.testing.assert_array_equal(out, asl)
    assert outputs.qc['anat_mask_synthstrip_dice'] == 0
    assert np.isinf(outputs.qc['anat_mask_synthstrip_volume_ratio'])
    assert outputs.qc['anat_mask_used'] == 0


def test_guard_anatomical_mask_grid(tmp_path):
    """Masks on different grids are an error, not a silent comparison."""
    with pytest.raises(ValueError, match='not on one grid'):
        _run_guard(tmp_path, ANAT, np.ones((20, 20, 21), dtype=np.uint8))
