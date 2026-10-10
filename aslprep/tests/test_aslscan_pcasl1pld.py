"""F1: single-delay 2D PCASL with a separate M0, scored against the truth (spec Section 5).

Exercises the slice-time PLD shift, SCORE/SCRUB, BASIL, T1w and MNI outputs and atlases.
The phantom's truth is flat within tissue, so no within-tissue correlation is scored.
"""

import pytest

from aslprep.tests import truth_bounds as tb
from aslprep.tests.aslscan_cli import run_recipe, scored_items, shared_items

RECIPE = 'pcasl1pld'
pytestmark = [pytest.mark.aslscan, pytest.mark.aslscan_pcasl1pld]


@pytest.fixture(scope='module')
def f1_run(data_dir, output_dir, working_dir):
    return run_recipe(
        RECIPE,
        'aslscan_pcasl1pld',
        data_dir,
        output_dir,
        working_dir,
        spaces=('asl', 'T1w', 'MNI152NLin2009cAsym'),
        extra_args=('--scorescrub', '--basil', '--atlases', '4S156Parcels', '4S1056Parcels'),
        score_spaces=('T1w', 'MNI152NLin2009cAsym'),
        extra_cbf=('basil', 'score', 'scrub'),
    )


globals().update(
    shared_items(
        'f1_run',
        report_only=[
            ('desc-basil', 'GM', 'median'),
            ('desc-basil', 'WM', 'median'),
            ('native', 'purity', 'GM'),
            ('native', 'tier_a', 'p95_abs_dev'),
            ('native', 'tier_a', 'interior', 'p95_abs_dev'),
            ('native', 'tier_a', 'edge', 'median'),
            ('native', 'hmc_deltam', 'p95_abs_dev'),
            ('frames', 'motion', 'rms_error_max_mm'),
        ],
    )
)
globals().update(scored_items('f1_run', RECIPE, spaces=('space-T1w', 'space-MNI152NLin2009cAsym')))


@pytest.mark.parametrize('desc', ['score', 'scrub'])
@pytest.mark.parametrize('tissue', ['GM', 'WM'])
def test_scorescrub(f1_run, desc, tissue):
    """Without outliers, SCORE and SCRUB stay at the mean CBF's accuracy."""
    path = (f'desc-{desc}', '{t}', 'median')
    reference = ('native', 'tier_b', '{t}', 'median')
    tb.check_pair(f1_run.score, 'scorescrub', RECIPE, path, reference, t=tissue)
