"""Tests for aslprep.interfaces.hmc (label offset correction after head-motion correction)."""

import json
import shutil

import nibabel as nb
import numpy as np
import pytest

from aslprep.interfaces.hmc import CorrectLabelOffset, _resample, plan_label_offset_groups

needs_afni = pytest.mark.skipif(shutil.which('3dvolreg') is None, reason='AFNI 3dvolreg not found')


def test_plan_label_offset_groups_pcasl_borrows_early_delays():
    """A PCASL delay under 0.5 s takes the estimate of the nearest other delay in time."""
    types = ['label', 'control'] * 6
    plds = [0.25, 0.25, 0.25, 0.25, 1.0, 1.0, 1.0, 1.0, 1.5, 1.5, 1.5, 1.5]
    estimated, borrowed = plan_label_offset_groups(types, plds, 'PCASL')
    assert sorted(estimated) == [1.0, 1.5]
    assert estimated[1.0] == ([4, 6], [5, 7])
    assert borrowed == {0.25: ([0, 2], 1.0)}


def test_plan_label_offset_groups_pasl_keeps_its_delays():
    """PASL inversion times at or above 0.8 s are each estimated from their own pair."""
    types = ['label', 'control'] * 3
    plds = [0.9, 0.9, 1.8, 1.8, 2.7, 2.7]
    estimated, borrowed = plan_label_offset_groups(types, plds, 'PASL')
    assert sorted(estimated) == [0.9, 1.8, 2.7]
    assert borrowed == {}


def test_plan_label_offset_groups_late_delay_without_controls_is_left_alone():
    """Only early delays borrow: a later delay with labels but no controls is not corrected."""
    types = ['label', 'control', 'label', 'control', 'label', 'label']
    plds = [1.0, 1.0, 1.0, 1.0, 2.0, 2.0]
    estimated, borrowed = plan_label_offset_groups(types, plds, 'PCASL')
    assert sorted(estimated) == [1.0]
    assert borrowed == {}


def test_plan_label_offset_groups_all_early_and_no_pairs():
    """When every delay is early each is estimated; without labels nothing is."""
    estimated, borrowed = plan_label_offset_groups(
        ['label', 'control'] * 2, [0.2, 0.2, 0.3, 0.3], 'PCASL'
    )
    assert sorted(estimated) == [0.2, 0.3]
    assert borrowed == {}
    estimated, borrowed = plan_label_offset_groups(['deltam'] * 3, [1.8] * 3, 'PCASL')
    assert estimated == {}
    assert borrowed == {}


def _phantom(shape, affine, shift_mm=(0.0, 0.0, 0.0)):
    """A smooth analytic head-like image, sampled with its content displaced by ``shift_mm``."""
    ijk = np.indices(shape).reshape(3, -1)
    xyz = (affine @ np.vstack([ijk, np.ones((1, ijk.shape[1]))]))[:3]
    xyz = xyz - np.asarray(shift_mm, float)[:, None]
    centre = (affine @ np.r_[(np.array(shape) - 1) / 2, 1])[:3]
    r = xyz - centre[:, None]
    head = 1.0 / (1.0 + np.exp(((r[0] / 40) ** 2 + (r[1] / 48) ** 2 + (r[2] / 32) ** 2 - 1) * 12))
    blobs = sum(
        a * np.exp(-((r - np.array(c)[:, None]) ** 2).sum(0) / (2 * s**2))
        for a, c, s in ((0.6, (12, 8, 0), 7), (-0.4, (-14, -6, 6), 9), (0.5, (0, -20, -8), 6))
    )
    return (1000 * (head + 0.3 * head * blobs)).reshape(shape).astype(np.float32)


def _write_series(tmp_path, volumes, types, affine, pulls=None):
    import nitransforms as nt

    asl = tmp_path / 'asl.nii.gz'
    nb.Nifti1Image(np.stack(volumes, -1), affine).to_filename(asl)
    ctx = tmp_path / 'aslcontext.tsv'
    ctx.write_text('\n'.join(['volume_type', *types]) + '\n')
    ref = tmp_path / 'ref.nii.gz'
    nb.Nifti1Image(volumes[0], affine).to_filename(ref)
    xfm = tmp_path / 'xfms.txt'
    if pulls is None:
        pulls = np.repeat(np.eye(4)[None], len(volumes), axis=0)
    nt.linear.LinearTransformsMapping(pulls, reference=nb.load(ref)).to_filename(xfm, fmt='itk')
    return asl, ctx, ref, xfm


def _pose(rot_deg, trans_mm, centre):
    """A rigid RAS transform: rotation about z through ``centre``, then a translation."""
    a = np.radians(rot_deg)
    m = np.eye(4)
    m[:2, :2] = [[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]]
    m[:3, 3] = centre - m[:3, :3] @ centre + np.asarray(trans_mm, float)
    return m


@needs_afni
def test_correct_label_offset_recovers_a_known_shift(tmp_path):
    """Labels displaced by 0.6 mm from the controls are brought back onto them."""
    import nitransforms as nt

    shape, affine = (40, 44, 30), np.diag([3.0, 3.0, 3.0, 1.0])
    control = _phantom(shape, affine)
    label = _phantom(shape, affine, shift_mm=(0.6, 0.0, 0.0)) - 0.01 * control  # perfusion
    types = ['label', 'control'] * 4
    vols = [label if t == 'label' else control for t in types]
    asl, ctx, ref, xfm = _write_series(tmp_path, vols, types, affine)
    res = CorrectLabelOffset(
        asl_file=str(asl),
        aslcontext=str(ctx),
        metadata={'PostLabelingDelay': 1.8, 'ArterialSpinLabelingType': 'PCASL'},
        ref_file=str(ref),
        xforms=str(xfm),
    ).run(cwd=str(tmp_path))
    record = json.loads(open(res.outputs.out_report).read())
    assert record['delays']['1.8']['translation_mm'] == pytest.approx(0.6, abs=0.1)
    out = np.asarray(nt.linear.load(res.outputs.out_file, fmt='itk').matrix).reshape(-1, 4, 4)
    for i, t in enumerate(types):
        if t == 'control':
            np.testing.assert_allclose(out[i], np.eye(4), atol=1e-6)
    # Resampled through the corrected transforms, the label lines up with the control again.
    ref_img = nb.load(ref)
    fixed = _resample(label, affine, out[0], ref_img)
    inner = (slice(5, -5),) * 3
    before = np.abs(label - 0.99 * control)[inner].mean()
    after = np.abs(fixed - 0.99 * control)[inner].mean()
    assert after < 0.3 * before


@needs_afni
def test_correct_label_offset_composes_with_a_moved_head(tmp_path):
    """With non-identity head-motion transforms the correction is composed in the right order.

    Every volume shows the head at a rotated, translated pose ``H`` (the HMC pull from the
    reference to the volume), and the labels are additionally displaced by 0.6 mm in the
    reference frame. Composing the correction the wrong way round (``correction @ H``) would
    misplace the labels, because the rotation and the offset do not commute.
    """
    import nitransforms as nt

    shape, affine = (40, 44, 30), np.diag([3.0, 3.0, 3.0, 1.0])
    ref_img = nb.Nifti1Image(np.zeros(shape, np.float32), affine)
    centre = (affine @ np.r_[(np.array(shape) - 1) / 2, 1])[:3]
    # A large rotation, so that the two composition orders differ measurably (by about
    # 2 * 0.6 * sin(15 degrees) = 0.31 mm); a realistic few degrees would differ by ~0.04 mm.
    pose = _pose(30.0, (1.5, -1.0, 0.5), centre)
    static_c = _phantom(shape, affine)
    static_l = _phantom(shape, affine, shift_mm=(0.6, 0.0, 0.0)) - 0.01 * static_c
    # Volume images V with V(H y) = static(y): V = static resampled through H^-1.
    vol_c = _resample(static_c, affine, np.linalg.inv(pose), ref_img)
    vol_l = _resample(static_l, affine, np.linalg.inv(pose), ref_img)
    types = ['label', 'control'] * 4
    vols = [vol_l if t == 'label' else vol_c for t in types]
    pulls = np.repeat(pose[None], len(types), axis=0)
    asl, ctx, ref, xfm = _write_series(tmp_path, vols, types, affine, pulls=pulls)
    res = CorrectLabelOffset(
        asl_file=str(asl),
        aslcontext=str(ctx),
        metadata={'PostLabelingDelay': 1.8, 'ArterialSpinLabelingType': 'PCASL'},
        ref_file=str(ref),
        xforms=str(xfm),
    ).run(cwd=str(tmp_path))
    out = np.asarray(nt.linear.load(res.outputs.out_file, fmt='itk').matrix).reshape(-1, 4, 4)
    np.testing.assert_allclose(out[1], pose, atol=1e-4)  # controls keep their transforms
    control_ref = _resample(vol_c, affine, pose, ref_img)
    label_before = _resample(vol_l, affine, pose, ref_img)
    label_after = _resample(vol_l, affine, out[0], ref_img)
    correction = np.linalg.inv(pose) @ out[0]
    label_wrong = _resample(vol_l, affine, correction @ pose, ref_img)
    inner = (slice(6, -6),) * 3

    def err(img):
        return np.abs(img - 0.99 * control_ref)[inner].mean()

    assert err(label_after) < 0.3 * err(label_before)
    assert err(label_after) < 0.5 * err(label_wrong)


@needs_afni
def test_correct_label_offset_passes_deltam_through(tmp_path):
    """A series without control-label pairs keeps its transforms."""
    import nitransforms as nt

    shape, affine = (20, 20, 16), np.diag([3.0, 3.0, 3.0, 1.0])
    vol = _phantom(shape, affine)
    asl, ctx, ref, xfm = _write_series(tmp_path, [vol, vol, vol], ['deltam'] * 3, affine)
    res = CorrectLabelOffset(
        asl_file=str(asl),
        aslcontext=str(ctx),
        metadata={'PostLabelingDelay': 1.8},
        ref_file=str(ref),
        xforms=str(xfm),
    ).run(cwd=str(tmp_path))
    out = np.asarray(nt.linear.load(res.outputs.out_file, fmt='itk').matrix).reshape(-1, 4, 4)
    np.testing.assert_allclose(out, np.repeat(np.eye(4)[None], 3, axis=0), atol=1e-6)


def test_magnitude_keeps_small_rotations():
    """A 0.02 degree rotation whose diagonal rounds to 1 is still reported."""
    from aslprep.interfaces.hmc import _magnitude

    ref = nb.Nifti1Image(np.zeros((10, 10, 10), np.float32), np.eye(4))
    a = np.radians(0.02)
    m = np.eye(4)
    m[0, 1], m[1, 0] = -np.round(np.sin(a), 6), np.round(np.sin(a), 6)  # diagonal stays 1.0
    assert _magnitude(m, ref)['rotation_deg'] == pytest.approx(0.02, rel=0.05)


def test_init_asl_hmc_wf_adds_label_offset_only_when_asked():
    """The correction node exists only with --hmc-method typewise_offset."""
    from aslprep import config
    from aslprep.tests.tests import mock_config
    from aslprep.workflows.asl.hmc import init_asl_hmc_wf

    with mock_config():
        config.workflow.hmc_method = 'typewise'
        wf = init_asl_hmc_wf(mem_gb=1, omp_nthreads=1)
        assert 'correct_label_offset' not in wf.list_node_names()
        config.workflow.hmc_method = 'typewise_offset'
        wf = init_asl_hmc_wf(mem_gb=1, omp_nthreads=1)
        assert 'correct_label_offset' in wf.list_node_names()
        config.workflow.hmc_method = 'typewise'


def test_identity_transforms(tmp_path):
    """One identity transform per volume, readable as the HMC transforms are."""
    from nitransforms.linear import load

    from aslprep.interfaces.hmc import IdentityTransforms

    ref = tmp_path / 'ref.nii.gz'
    nb.Nifti1Image(np.zeros((8, 8, 8), np.float32), np.diag([2.0, 2.0, 3.0, 1.0])).to_filename(ref)
    ctx = tmp_path / 'aslcontext.tsv'
    ctx.write_text('volume_type\nm0scan\ndeltam\n')
    res = IdentityTransforms(aslcontext=str(ctx), ref_file=str(ref)).run(cwd=str(tmp_path))
    mats = np.asarray(load(res.outputs.out_file, fmt='itk').matrix).reshape(-1, 4, 4)
    assert mats.shape == (2, 4, 4)
    np.testing.assert_allclose(mats, np.tile(np.eye(4), (2, 1, 1)), atol=1e-12)


@pytest.mark.parametrize(
    ('volume_types', 'use_ge', 'skipped'),
    [
        (['m0scan', 'deltam'], True, True),
        (['m0scan', 'deltam', 'cbf'], True, True),
        (['control', 'label', 'control', 'label'], True, False),
        (['m0scan', 'control', 'label'], True, False),
        (['m0scan', 'deltam'], False, False),
    ],
)
def test_skips_motion_correction(tmp_path, volume_types, use_ge, skipped):
    """Only GE-path series without label-control pairs skip motion correction (#769)."""
    from aslprep.workflows.asl.fit import skips_motion_correction

    ctx = tmp_path / 'aslcontext.tsv'
    ctx.write_text('volume_type\n' + '\n'.join(volume_types) + '\n')
    assert skips_motion_correction(str(ctx), use_ge) is skipped


def test_ge_path_skips_motion_correction(tmp_path):
    """The GE-specific path writes identity transforms instead of running MCFLIRT (#769)."""
    from aslprep import config
    from aslprep.tests.tests import mock_config, reset_config
    from aslprep.workflows.asl.fit import init_asl_fit_wf

    deltam_ctx = tmp_path / 'aslcontext.tsv'
    deltam_ctx.write_text('volume_type\nm0scan\ndeltam\n')
    reset_config()
    with mock_config():
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        paired_ctx = str(asl_file).replace('.nii.gz', 'context.tsv')
        kwargs = {'asl_file': str(asl_file), 'm0scan': None, 'precomputed': {}, 'use_ge': True}
        ge = init_asl_fit_wf(aslcontext=str(deltam_ctx), **kwargs)
        assert ge.get_node('asl_hmc_wf') is None
        assert ge.get_node('identity_xforms') is not None
        sink = ge.get_node('ds_hmc_wf').get_node('ds_xforms')
        assert sink.inputs.meta_dict == {'HMCMethod': 'none'}
        # A short series with label-control pairs is still motion corrected.
        paired = init_asl_fit_wf(aslcontext=paired_ctx, **kwargs)
        assert paired.get_node('asl_hmc_wf') is not None
        assert paired.get_node('identity_xforms') is None
        sink = paired.get_node('ds_hmc_wf').get_node('ds_xforms')
        assert sink.inputs.meta_dict == {'HMCMethod': config.workflow.hmc_method}
