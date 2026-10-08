"""F7: head motion, multiband, and coregistration under a known pose (spec Section 5).

The ASL series moves along a known rigid trajectory, the T1w header carries a known rigid
offset, and ASLPrep preprocesses the anatomy itself. Scored: head-motion correction against
the trajectory, coregistration against the offset, and the alignment of T1w-space CBF.
Native-space Tier A/B are not scored: the motion-corrected series sits in the aslref frame,
not the phantom's.
"""

import pytest

from aslprep.tests import truth_bounds as tb
from aslprep.tests.aslscan_cli import run_recipe, shared_items

RECIPE = 'motion'
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
            ('desc-score', 'GM', 'median'),
            ('desc-scrub', 'GM', 'median'),
        ],
    )
)


def test_motion_correction_median(f7_run):
    tb.check(f7_run.score, 'motion_median', RECIPE)


def test_motion_correction_worst_volume(f7_run):
    tb.check(f7_run.score, 'motion_max', RECIPE)


def test_coregistration(f7_run):
    tb.check(f7_run.score, 'coreg_rot', RECIPE)
    tb.check(f7_run.score, 'coreg_rms', RECIPE)


def test_t1w_space_alignment(f7_run):
    tb.check(f7_run.score, 'space_alignment', RECIPE, space='space-T1w')


@pytest.mark.parametrize('tissue', ['GM', 'WM'])
def test_t1w_space_coverage(f7_run, tissue):
    tb.check(f7_run.score, 'space_finite', RECIPE, space='space-T1w', t=tissue)
