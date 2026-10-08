"""F2: GE-style 3D spiral PCASL, one delta-M volume with an included M0 (spec Section 5).

Replaces examples_pcasl_singlepld_ge: the precomputed delta-M path, the included M0 stored
divided by 96 (``--m0_scale=96``), four background-suppression pulses with no
LabelingEfficiency in the sidecar, BASIL, SCORE/SCRUB and atlases.
"""

import pytest

from aslprep.tests.aslscan_cli import run_recipe, scored_items, shared_items

RECIPE = 'pcasl1pld_ge3d'
pytestmark = [pytest.mark.aslscan, pytest.mark.aslscan_pcasl1pld_ge3d]


@pytest.fixture(scope='module')
def f2_run(data_dir, output_dir, working_dir):
    return run_recipe(
        RECIPE,
        'aslscan_pcasl1pld_ge3d',
        data_dir,
        output_dir,
        working_dir,
        extra_args=('--scorescrub', '--basil', '--atlases', '4S156Parcels'),
        extra_cbf=('basil',),
    )


globals().update(
    shared_items(
        'f2_run',
        report_only=[
            ('desc-basil', 'GM', 'median'),
            ('native', 'tier_a', 'p95_abs_dev'),
            ('native', 'tier_a', 'edge', 'median'),
        ],
    )
)
globals().update(
    scored_items(
        'f2_run',
        RECIPE,
        suppression_pulses=4,
        known={
            'test_aslref_pose': (
                'Head-motion correction registers the delta-M volume to an M0-like reference '
                '(no shared contrast) and applies the spurious result: on motion-free data the '
                'delta-M moves by about 1.9 degrees and 1.7 mm RMS. Reported; pending a decision.'
            ),
        },
    )
)
