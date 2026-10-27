"""The tissue-T1 kinetic model on the long-delay 3D GRASE aslscan fixture (F3, gh-744).

aslscan decays the label with the apparent T1 of tissue after it arrives, so on F3, whose delays
(0.5-4.0 s) bracket every voxel's ATT, the consensus model underestimates white matter ATT and
CBF. These tests fit the fixture's delta-M directly, on nearly pure tissue voxels, so they check
the kinetic model against the simulator's physics without the rest of the pipeline.
"""

import numpy as np
import pytest

from aslprep.tests import aslscan_fixtures as af
from aslprep.tests import truth_models as tm
from aslprep.utils import cbf

RECIPE = 'pcasl_multipld_grase'
#: Partial volume above which a voxel counts as pure tissue
PURE = 0.95
#: Voxels fit per tissue
N_VOXELS = 300


@pytest.fixture(scope='module')
def f3(data_dir):
    from aslprep.tests.truth_scoring import Fixture

    try:
        path = af.fixture_dir(RECIPE, data_dir)
    except af.FixtureUnavailable as exc:
        pytest.skip(str(exc))

    fixture = Fixture(path, af.load_recipe(RECIPE).acq)
    metadata = fixture.metadata
    asl = fixture.asl_img.get_fdata()
    deltam, rows = tm.deltam_observations(asl, fixture.context)
    m0, tr = tm.m0_image(
        asl,
        fixture.context,
        metadata,
        fixture.affine,
        fixture.m0scan,
        fixture.m0_metadata,
    )
    m0 = tm.m0_tr_correction(m0, tr, metadata['MagneticFieldStrength']) * fixture.m0_scale
    plds = np.asarray(metadata['PostLabelingDelay'], float)[rows]
    return fixture, deltam, m0, plds


def _fit(f3, tissue, t1prime=None):
    """Fit the voxels where ``tissue`` is pure, returning the median CBF ratio and ATT error."""
    fixture, deltam, m0, plds = f3
    pv = fixture.pv
    voxels = np.argwhere((pv[tissue] >= PURE) & fixture.brain)
    rng = np.random.default_rng(0)
    voxels = voxels[rng.choice(len(voxels), min(N_VOXELS, len(voxels)), replace=False)]
    i, j, k = voxels.T
    kwargs = {}
    if t1prime is not None:
        kwargs = {
            't1prime': t1prime,
            'tissue_weights': cbf.tissue_mixing_weights(
                pv['GM'][i, j, k],
                pv['WM'][i, j, k],
                cbf.NOMINAL_CBF['gm'],
                cbf.NOMINAL_CBF['wm'],
            ),
        }

    fit_cbf, fit_att, *_ = cbf.fit_deltam_multipld(
        deltam_arr=deltam[i, j, k],
        scaled_m0data=m0[i, j, k],
        plds=np.tile(plds, (len(voxels), 1)),
        labeleff=tm.labeling_efficiency(fixture.metadata),
        t1blood=tm.T1_BLOOD[fixture.metadata['MagneticFieldStrength']],
        partition_coefficient=tm.LAMBDA,
        is_casl=True,
        tau=fixture.metadata['LabelingDuration'],
        **kwargs,
    )
    return (
        float(np.median(fit_cbf / fixture.perfusion[i, j, k])),
        float(np.median(fit_att - fixture.att[i, j, k])),
    )


def _t1prime(gm_t1, wm_t1):
    return np.array(
        [
            cbf.apparent_tissue_t1(gm_t1, cbf.NOMINAL_CBF['gm']),
            cbf.apparent_tissue_t1(wm_t1, cbf.NOMINAL_CBF['wm']),
        ]
    )


#: The default tissue T1s at 3 T (Wright et al., 2008)
DEFAULT_T1PRIME = _t1prime(cbf.T1TISSUE_BY_FIELD_STRENGTH[3], cbf.T1WM_BY_FIELD_STRENGTH[3])
#: The tissue T1s aslscan simulated
MATCHED_T1PRIME = _t1prime(af.TISSUES[1]['T1map'], af.TISSUES[2]['T1map'])


def test_tissue_t1_locates_white_matter_att(f3):
    """The tissue-T1 model removes most of the consensus model's white matter ATT bias.

    The default white matter T1 (0.838 s) is close to aslscan's (0.83 s).
    On the generated fixture, the median errors are about -0.37 s with the consensus model and
    -0.06 s with the tissue-T1 model. White matter CBF is not asserted: both T1s give about 1.3
    times the truth, which the kinetic model does not account for.
    """
    _, consensus_att = _fit(f3, 'WM')
    _, tissue_att = _fit(f3, 'WM', DEFAULT_T1PRIME)
    assert consensus_att < -0.25
    assert abs(tissue_att) < 0.15


def test_tissue_t1_matched_gray_matter(f3):
    """With aslscan's gray matter T1, the tissue-T1 model recovers gray matter CBF and ATT.

    With the default gray matter T1 (1.607 s, against aslscan's 1.33 s), CBF stays low
    (about 0.78 times the truth), which shows how much the model depends on tissue T1.
    """
    matched_cbf, matched_att = _fit(f3, 'GM', MATCHED_T1PRIME)
    default_cbf, _ = _fit(f3, 'GM', DEFAULT_T1PRIME)
    assert matched_cbf == pytest.approx(1.0, abs=0.05)
    assert abs(matched_att) < 0.03
    assert default_cbf < matched_cbf - 0.1
