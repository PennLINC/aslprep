"""F8: head motion between and within the volumes of a segmented 3D GRASE series (spec Section 5).

Each GRASE volume is read in four segments, one per labeling, and the head moves along a known
timed trajectory, so the segments of a volume can be read at different poses (ghosting and blur
that no rigid transform undoes). The anatomical derivatives carry a known rigid offset. Scored:
head-motion correction against the per-volume rigid fit to the segments' poses, coregistration
against the offset, and the T1w-space CBF. Native-space Tier A/B are not asserted, for the reason
given in :mod:`aslprep.tests.test_aslscan_motion`.
"""

import pytest

from aslprep.tests import truth_bounds as tb
from aslprep.tests.aslscan_cli import Known, known_failure, run_recipe, shared_items

RECIPE = 'grasemotion'
pytestmark = [pytest.mark.aslscan, pytest.mark.aslscan_grase_motion]


@pytest.fixture(scope='module')
def f8_run(data_dir, output_dir, working_dir):
    return run_recipe(
        RECIPE,
        'aslscan_grase_motion',
        data_dir,
        output_dir,
        working_dir,
        spaces=('asl', 'T1w'),
        extra_args=('--atlases', '4S156Parcels'),
        score_spaces=('T1w',),
    )


globals().update(
    shared_items(
        'f8_run',
        report_only=[
            ('confounds', 'fd_mean'),
            ('confounds', 'abs_r_trans_z'),
            ('confounds', 'abs_r_rot_x'),
            ('space-T1w', 'GM', 'median'),
            ('space-T1w', 'r_truth'),
            ('space-T1w', 'r_native'),
            ('frames', 'motion', 'rms_error_max_mm'),
            ('frames', 'motion', 'within_volume_rms_max_mm'),
        ],
    )
)


def test_motion_correction_median(f8_run):
    tb.check(f8_run.score, 'motion_median', RECIPE)


def test_motion_correction_worst_volume(f8_run):
    tb.check(f8_run.score, 'motion_max', RECIPE)


def test_coregistration_rotation(f8_run):
    tb.check(f8_run.score, 'coreg_rot', RECIPE)


@known_failure(
    'f8_run',
    Known(
        'PennLINC/aslprep#762: the end-to-end error composes the M0-to-T1w registration '
        '(0.17 voxels alone) with the constant offset left by aligning the background-'
        "suppressed series' mean to the M0 (0.13 voxels); measured 0.2521 voxels",
        ('frames', 'coreg', 'rms_voxels'),
        0.25,
        0.35,
    ),
)
def test_coregistration_displacement(request):
    tb.check(request.getfixturevalue('f8_run').score, 'coreg_rms', RECIPE)


def test_t1w_space_resampling(f8_run):
    """The T1w map is the native map moved through the scored coregistration."""
    tb.check(f8_run.score, 'space_resampling', RECIPE, space='space-T1w')


@pytest.mark.parametrize('tissue', ['GM', 'WM'])
def test_t1w_space_coverage(f8_run, tissue):
    tb.check(f8_run.score, 'space_finite', RECIPE, space='space-T1w', t=tissue)
