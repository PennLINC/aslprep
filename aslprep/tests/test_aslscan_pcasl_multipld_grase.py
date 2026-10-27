"""F3: Siemens-style 3D GRASE multi-delay PCASL with a short-TR separate M0 (spec Section 5).

3D multi-delay PCASL (PLDs 0.5-4.0 s), the M0 TR < 5 s correction, the M0 stored divided by
10 (``--m0_scale=10``), four suppression pulses with LabelingEfficiency in the sidecar, BASIL,
MNI output and atlases. Without slice-timing offsets, the delays locate ATT (issue #742).
"""

import pytest

from aslprep.tests.aslscan_cli import run_recipe, scored_items, shared_items

RECIPE = 'pcasl_multipld_grase'
pytestmark = [pytest.mark.aslscan, pytest.mark.aslscan_pcasl_multipld_grase]


@pytest.fixture(scope='module')
def f3_run(data_dir, output_dir, working_dir):
    return run_recipe(
        RECIPE,
        'aslscan_pcasl_multipld_grase',
        data_dir,
        output_dir,
        working_dir,
        spaces=('asl', 'MNI152NLin2009cAsym'),
        extra_args=(
            '--skull-strip-asl',
            'anatomical',
            '--basil',
            '--atlases',
            '4S156Parcels',
            '4S1056Parcels',
        ),
        score_spaces=('MNI152NLin2009cAsym',),
        extra_cbf=('basil',),
    )


globals().update(
    shared_items(
        'f3_run',
        multi_delay=True,
        report_only=[
            ('desc-basil', 'GM', 'median'),
            ('native', 'tier_a', 'p95_abs_dev'),
            ('native', 'tier_a_quant', 'p95_abs_dev_multi'),
            ('native', 'tier_b_att', 'GM', 'median_abs_error'),
            ('native', 'tier_b_att', 'WM', 'median_abs_error'),
            ('native', 'tier_a_quant', 'abat_median_abs_diff'),
            ('native', 'tier_a_quant', 'abv_median_abs_diff'),
        ],
    )
)
globals().update(
    scored_items(
        'f3_run',
        RECIPE,
        spaces=('space-MNI152NLin2009cAsym',),
        suppression_pulses=4,
    )
)
