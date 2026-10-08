"""F6: multi-delay 2D PASL with Q2TIPS (spec Section 5).

Replaces examples_pasl_multipld: the multi-delay PASL fit, ``--m0_scale=10``, BASIL,
SCORE/SCRUB and atlases.
"""

import pytest

from aslprep.tests.aslscan_cli import run_recipe, scored_items, shared_items

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
            ('native', 'tier_b_att', 'GM', 'median_abs_error'),
        ],
    )
)
globals().update(scored_items('f6_run', RECIPE))
