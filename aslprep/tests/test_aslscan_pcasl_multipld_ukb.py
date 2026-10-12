"""F8: a UK Biobank-like acquisition (spec Section 5; recipe ``pcasl_multipld_ukb``).

Multi-delay PCASL with a 2-shot 3D GRASE readout: coarse voxels (3.4 x 3.4 x 4.5 mm), heavy
through-slice blurring, one label/control pair per delay, and the ASL series stored ten times
larger than the M0 (``--m0_scale=10``), run with BASIL as in
https://neurostars.org/t/aslprep-csf-shows-elevated-cbf-in-basil-output-multi-pld-uk-biobank-data/36453.

The tissue contrast that thread asks about is recorded, not asserted (``tissue_contrast``):
CBF in voxels that are truly CSF is far from zero, and white matter approaches grey matter,
because of partial volume, blurring and noise rather than a defect.

Not asserted, with reasons:

- White-matter Tier B. With one pair per delay the fitted subset holds about 80 white-matter
  voxels at a delta-M SNR near 1, and their median moves by 0.1 between noise draws.
- End-to-end coregistration displacement. The constant motion-correction offset
  (PennLINC/aslprep#707, 0.5 mm here) puts it at the 0.25-voxel ceiling, so the registration
  step alone is asserted instead.
- The correlation of standard-space CBF with native CBF (``r_native``). Noise lowers it (0.77
  in T1w space) below the floor set on the low-noise recipes; alignment is asserted through
  the coregistration items and the T1w resampling margin, as in F7.
"""

import pytest

from aslprep.tests import truth_bounds as tb
from aslprep.tests.aslscan_cli import run_recipe, shared_items

RECIPE = 'pcasl_multipld_ukb'
SPACES = ('space-T1w', 'space-MNI152NLin2009cAsym')
pytestmark = [pytest.mark.aslscan, pytest.mark.aslscan_pcasl_multipld_ukb]


@pytest.fixture(scope='module')
def f8_run(data_dir, output_dir, working_dir):
    return run_recipe(
        RECIPE,
        'aslscan_pcasl_multipld_ukb',
        data_dir,
        output_dir,
        working_dir,
        spaces=('asl', 'T1w', 'MNI152NLin2009cAsym'),
        extra_args=('--basil',),
        score_spaces=('T1w', 'MNI152NLin2009cAsym'),
        extra_cbf=('basil', 'basilGM', 'basilWM'),
    )


globals().update(
    shared_items(
        'f8_run',
        multi_delay=True,
        report_only=[
            ('native', 'tier_a', 'p95_abs_dev'),
            ('native', 'tier_a_quant', 'p95_abs_dev_multi'),
            ('native', 'tier_a_quant', 'abat_median_abs_diff'),
            ('native', 'tier_a_quant', 'abv_median_abs_diff'),
            ('native', 'tier_b_att', 'GM', 'median_abs_error'),
            ('native', 'tier_b_matched', 'WM', 'median'),
            ('frames', 'coreg', 'rms_voxels'),
            ('space-T1w', 'r_native'),
            ('space-MNI152NLin2009cAsym', 'r_native'),
            *(
                ('tissue_contrast', key, tissue)
                for key in ('native', 'desc-basil', 'desc-basilGM', 'desc-basilWM')
                for tissue in ('GM', 'WM', 'CSF')
            ),
        ],
    )
)


def test_tier_a_median(f8_run):
    tb.check(f8_run.score, 'tier_a_median', RECIPE)


def test_tier_a_quantification(f8_run):
    tb.check(f8_run.score, 'tier_a_quant_median', RECIPE)
    tb.check(f8_run.score, 'tier_a_quant_agree', RECIPE)
    tb.check(f8_run.score, 'tier_a_att', RECIPE)


def test_motion_correction(f8_run):
    tb.check(f8_run.score, 'motion_median', RECIPE)
    tb.check(f8_run.score, 'motion_max', RECIPE)


def test_tier_b_grey_matter(f8_run):
    reference = ('native', 'expected_ratio', '{t}')
    tb.check(f8_run.score, 'tier_b_matched', RECIPE, reference=reference, t='GM')


def test_coregistration_rotation(f8_run):
    tb.check(f8_run.score, 'coreg_rot', RECIPE)


def test_coregistration_step(f8_run):
    tb.check(f8_run.score, 'coreg_only_rms', RECIPE)


@pytest.mark.parametrize('space', SPACES)
def test_standard_space(f8_run, space):
    score = f8_run.score
    for tissue in ('GM', 'WM'):
        tb.check(score, 'space_finite', RECIPE, space=space, t=tissue)
    if space == 'space-T1w':
        tb.check(score, 'space_resampling', RECIPE, space=space)
