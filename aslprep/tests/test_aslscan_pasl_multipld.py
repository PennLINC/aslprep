"""F6: multi-delay 2D PASL with Q2TIPS (spec Section 5).

Replaces examples_pasl_multipld: the multi-delay PASL fit, ``--m0_scale=10``, BASIL,
SCORE/SCRUB and atlases.
"""

import pytest

from aslprep.tests.aslscan_cli import Known, run_recipe, scored_items, shared_items

RECIPE = 'pasl_multipld'
pytestmark = [pytest.mark.aslscan, pytest.mark.aslscan_pasl_multipld]


@pytest.fixture(scope='module')
def f6_run(data_dir, output_dir, working_dir):
    return run_recipe(
        RECIPE,
        'aslscan_pasl_multipld',
        data_dir,
        output_dir,
        working_dir,
        extra_args=('--scorescrub', '--basil', '--atlases', '4S156Parcels', '4S1056Parcels'),
        extra_cbf=('basil',),
    )


globals().update(
    shared_items(
        'f6_run',
        multi_delay=True,
        report_only=[
            ('desc-basil', 'GM', 'median'),
            ('native', 'tier_a', 'p95_abs_dev'),
            ('native', 'tier_a_quant', 'p95_abs_dev_multi'),
            ('native', 'tier_b_att', 'GM', 'median_abs_error'),
            ('native', 'tier_a_quant', 'abat_median_abs_diff'),
            ('native', 'tier_a_quant', 'abv_median_abs_diff'),
        ],
        known={
            'test_fit_bound': Known(
                'PennLINC/aslprep#705: motion correction misaligns mid-delay control-label '
                'pairs by 0.12-0.16 mm, pegging about 1 % of edge brain voxels at CBF 300',
                ('native', 'coverage', 'at_fit_bound'),
                0.005,
                0.03,
            ),
        },
    )
)
globals().update(
    scored_items(
        'f6_run',
        RECIPE,
        known={
            'test_coregistration': Known(
                'Motion correction leaves the motion-free series rotated about 1 degree from '
                'its reference (the registration step alone is 0.09 voxel); see the plan',
                ('frames', 'coreg', 'rot_arc_voxels'),
                0.25,
                0.5,
            ),
        },
    )
)
