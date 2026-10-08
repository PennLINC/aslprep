"""F3: Siemens-style 3D GRASE PCASL with a short-TR separate M0 (spec Section 5).

Replaces examples_pcasl_singlepld_siemens: 3D single-delay PCASL, the M0 TR < 5 s correction,
the M0 stored divided by 10 (``--m0_scale=10``), four suppression pulses with
LabelingEfficiency in the sidecar, BASIL, MNI output and atlases.
"""

import pytest

from aslprep.tests.aslscan_cli import run_recipe, scored_items, shared_items

RECIPE = 'pcasl1pld_grase'
pytestmark = [pytest.mark.aslscan, pytest.mark.aslscan_pcasl1pld_grase]


@pytest.fixture(scope='module')
def f3_run(data_dir, output_dir, working_dir):
    return run_recipe(
        RECIPE,
        'aslscan_pcasl1pld_grase',
        data_dir,
        output_dir,
        working_dir,
        spaces=('asl', 'MNI152NLin2009cAsym'),
        extra_args=('--basil', '--atlases', '4S156Parcels', '4S1056Parcels'),
        score_spaces=('MNI152NLin2009cAsym',),
        extra_cbf=('basil',),
    )


globals().update(
    shared_items(
        'f3_run',
        report_only=[
            ('desc-basil', 'GM', 'median'),
            ('native', 'tier_a', 'p95_abs_dev'),
            ('native', 'tier_a', 'edge', 'median'),
        ],
    )
)
globals().update(
    scored_items('f3_run', RECIPE, spaces=('space-MNI152NLin2009cAsym',), suppression_pulses=4)
)
