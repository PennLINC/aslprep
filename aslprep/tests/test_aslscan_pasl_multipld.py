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
                'pairs by 0.12-0.16 mm, pegging 2-3 % of brain voxels at CBF 300 (about 1 % with '
                'the local fit this search replaced)',
                ('native', 'coverage', 'at_fit_bound'),
                0.005,
                0.05,
            ),
        },
    )
)
globals().update(
    scored_items(
        'f6_run',
        RECIPE,
        known={
            'test_tier_a_median': Known(
                'PennLINC/aslprep#705: motion correction changes delta-M voxel by voxel '
                '(preprocessed / raw: 0.2 to 2.1 between the 5th and 95th percentiles in the '
                'interior), so CBF fitted to the preprocessed series is 7 % above CBF fitted '
                'to the raw series in the interior and twice it at the edge; given the same '
                'series, ASLPrep matches the model (tier_a_quant)',
                ('native', 'tier_a', 'median_dev'),
                0.05,
                0.3,
            ),
            'test_coregistration_rotation': Known(
                'PennLINC/aslprep#707: motion correction leaves the motion-free series rotated '
                '1 degree from its reference (the registration step alone is 0.09 voxel)',
                ('frames', 'coreg', 'rot_arc_voxels'),
                0.25,
                0.5,
            ),
        },
    )
)
