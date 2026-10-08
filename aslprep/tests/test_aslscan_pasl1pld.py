"""F4: single-delay 2D PASL with QUIPSS II and interleaved slices (spec Section 5).

New coverage (no existing real-data test runs single-delay PASL). QUIPSS II rather than
Q2TIPS: single-delay Q2TIPS has a known quantification discrepancy, pinned in the fast tier.
"""

import pytest

from aslprep.tests.aslscan_cli import run_recipe, scored_items, shared_items

RECIPE = 'pasl1pld'
pytestmark = [pytest.mark.aslscan, pytest.mark.aslscan_pasl1pld]


@pytest.fixture(scope='module')
def f4_run(data_dir, output_dir, working_dir):
    return run_recipe(
        RECIPE,
        'aslscan_pasl1pld',
        data_dir,
        output_dir,
        working_dir,
        # without --atlases ASLPrep parcellates with every atlas it ships
        extra_args=('--atlases', '4S156Parcels'),
    )


globals().update(
    shared_items(
        'f4_run',
        report_only=[('native', 'tier_a', 'p95_abs_dev'), ('native', 'tier_a', 'edge', 'median')],
    )
)
globals().update(scored_items('f4_run', RECIPE))
