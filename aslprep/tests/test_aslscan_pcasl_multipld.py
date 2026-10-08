"""F5: multi-delay 2D PCASL, label first, unequal repeats, arterial term (spec Section 5).

Replaces examples_pcasl_multipld: the multi-delay fit (CBF and ATT), label-first ordering,
``--m0_scale=10``, SCORE/SCRUB and atlases. aBV and aBAT are reported, not bounded.
"""

import pytest

from aslprep.tests.aslscan_cli import run_recipe, scored_items, shared_items

RECIPE = 'pcasl_multipld'
pytestmark = [pytest.mark.aslscan, pytest.mark.aslscan_pcasl_multipld]


@pytest.fixture(scope='module')
def f5_run(data_dir, output_dir, working_dir):
    return run_recipe(
        RECIPE,
        'aslscan_pcasl_multipld',
        data_dir,
        output_dir,
        working_dir,
        extra_args=('--scorescrub', '--atlases', '4S156Parcels'),
    )


globals().update(
    shared_items(
        'f5_run',
        multi_delay=True,
        report_only=[
            ('native', 'tier_a', 'p95_abs_dev'),
            ('native', 'tier_b_att', 'GM', 'median_abs_error'),
            ('native', 'tier_b_att', 'WM', 'median_abs_error'),
        ],
    )
)
globals().update(scored_items('f5_run', RECIPE))
