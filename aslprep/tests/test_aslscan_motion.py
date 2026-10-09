"""F7: head motion, multiband, and coregistration under a known pose (spec Section 5).

The ASL series moves along a known rigid trajectory, the T1w header carries a known rigid
offset, and ASLPrep preprocesses the anatomy itself. Scored: head-motion correction against
the trajectory, coregistration against the offset, and the alignment of T1w-space CBF.
Native-space Tier A/B and the T1w-space alignment correlation (``r_native``) are not
asserted: the motion-corrected series sits in the aslref frame, whose pose is unknown here, so
sampling it at the phantom's coordinates compares misaligned points. Alignment is asserted
through the coregistration error, which composes the motion-correction transforms, and the
T1w map is checked against the native map moved through that coregistration.
"""

import pytest

from aslprep.tests import truth_bounds as tb
from aslprep.tests.aslscan_cli import Known, known_failure, run_recipe, shared_items

RECIPE = 'headmotion'
pytestmark = [pytest.mark.aslscan, pytest.mark.aslscan_motion]


@pytest.fixture(scope='module')
def f7_run(data_dir, output_dir, working_dir):
    return run_recipe(
        RECIPE,
        'aslscan_motion',
        data_dir,
        output_dir,
        working_dir,
        spaces=('asl', 'T1w'),
        extra_args=('--scorescrub', '--atlases', '4S156Parcels'),
        score_spaces=('T1w',),
        extra_cbf=('score', 'scrub'),
    )


globals().update(
    shared_items(
        'f7_run',
        report_only=[
            ('confounds', 'fd_mean'),
            ('confounds', 'abs_r_trans_x'),
            ('confounds', 'abs_r_rot_z'),
            ('space-T1w', 'GM', 'median'),
            ('space-T1w', 'r_truth'),
            ('space-T1w', 'r_native'),
            ('frames', 'motion', 'rms_error_max_mm'),
            ('desc-score', 'GM', 'median'),
            ('desc-scrub', 'GM', 'median'),
        ],
    )
)


def test_motion_correction_median(f7_run):
    tb.check(f7_run.score, 'motion_median', RECIPE)


def test_motion_correction_worst_volume(f7_run):
    tb.check(f7_run.score, 'motion_max', RECIPE)


def test_coregistration_rotation(f7_run):
    tb.check(f7_run.score, 'coreg_rot', RECIPE)


@known_failure(
    'f7_run',
    Known(
        'RMS end-to-end coregistration error 0.27 voxel (1.37 mm at 5 mm slices) against the '
        '0.25-voxel ceiling; F3 and F6 show the same excess on motion-free data, from the '
        "motion-corrected series' constant offset to its reference; see the plan",
        ('frames', 'coreg', 'rms_voxels'),
        0.25,
        0.4,
    ),
)
def test_coregistration_displacement(request):
    tb.check(request.getfixturevalue('f7_run').score, 'coreg_rms', RECIPE)


def test_t1w_space_resampling(f7_run):
    """The T1w map is the native map moved through the scored coregistration."""
    tb.check(f7_run.score, 'space_resampling', RECIPE, space='space-T1w')


@pytest.mark.parametrize('tissue', ['GM', 'WM'])
def test_t1w_space_coverage(f7_run, tissue):
    tb.check(f7_run.score, 'space_finite', RECIPE, space='space-T1w', t=tissue)
