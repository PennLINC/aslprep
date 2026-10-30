"""Truth-free scores of the coregistration probe (plan: "Real data") and their building blocks.

Every transform is a RAS pull map, ASLPrep's convention for ``from-aslref_to-T1w``: a coregistration
``C`` maps a T1w point to an aslref point; a direct registration ``D_ab`` of run a (moving) to run b
(fixed) maps a b point to an a point. Distances are RMS displacements over a set of points (mm),
with ``score_full._transform_error``, the HMC probe's scorer.

- Inconsistency of an arm over runs a and b: ``D_ab^-1 C_a C_b^-1`` over run b's brain.
- Noise floor of a direct registration: ``D_ab D_ba`` over run a's brain.
- Semi-synthetic error: the T1w given the header ``M A`` (world moved by ``M``, no resampling), an
  arm with solution ``C`` on the original should return ``C M^-1``; with ``C'`` its solution on the
  moved T1w, the error is ``C'^-1 C M^-1`` over the brain in the moved world.
- Gross failure by the carried mask: the ASL brain mask carried into T1w space through ``C`` (a T1w
  point is inside if ``C`` maps it into the ASL mask) against the T1w brain mask, as Dice; an arm
  fails if its Dice falls below the run's best arm by more than the drop that a 5 mm shift of the
  best arm's own solution produces (median over the six axis shifts).
"""

import sys
from pathlib import Path

import nibabel as nb
import numpy as np
from nibabel.affines import apply_affine

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'hmc_probe'))

from score_full import _transform_error

#: Gross-failure limits (plan: "Gross failure"), fixed in Task 3; the defaults are the plan's.
GROSS_INCONSISTENCY_MM = 5.0
CALIBRATION_SHIFT_MM = 5.0
CALIBRATION_DRAWS = 12
CALIBRATION_SEED = 20261022


def distance(a, b, points):
    """RMS displacement (mm) of ``a^-1 b`` over ``points``: zero when the two maps agree."""
    return _transform_error(a, b, points)['rms_mm']


def inconsistency(d_ab, c_a, c_b, points_b):
    """RMS displacement of ``D_ab^-1 C_a C_b^-1`` over run b's brain points."""
    return distance(d_ab, c_a @ np.linalg.inv(c_b), points_b)


def noise_floor(d_ab, d_ba, points_a):
    """RMS displacement of ``D_ab D_ba`` (identity for a perfect pair) over run a's brain."""
    return distance(d_ab @ d_ba, np.eye(4), points_a)


def semi_synthetic_expected(c, m):
    """The pull an arm should return after the T1w's world is moved by ``m``."""
    return c @ np.linalg.inv(m)


def semi_synthetic_error(c_moved, c, m, points_moved):
    """RMS displacement of ``C'^-1 C M^-1`` over the brain points in the moved world."""
    return distance(c_moved, semi_synthetic_expected(c, m), points_moved)


def move_header(src, dst, m):
    """Write ``src`` with the affine ``m A``: the same voxels, the world moved by ``m``."""
    img = nb.load(src)
    out = img.__class__(np.asanyarray(img.dataobj), m @ img.affine, img.header)
    out.set_qform(m @ img.affine)
    out.set_sform(m @ img.affine)
    out.to_filename(dst)


def rigid(translation=(0, 0, 0), rotation_deg=(0, 0, 0), centre=(0, 0, 0)):
    """A rigid RAS transform: rotations about x, y, z (in that order) about ``centre``, then a shift."""
    rx, ry, rz = np.deg2rad(rotation_deg)
    cx, sx, cy, sy, cz, sz = np.cos(rx), np.sin(rx), np.cos(ry), np.sin(ry), np.cos(rz), np.sin(rz)
    rot_x = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    rot_y = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    rot_z = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    out = np.eye(4)
    out[:3, :3] = rot_z @ rot_y @ rot_x
    centre = np.asarray(centre, float)
    out[:3, 3] = centre - out[:3, :3] @ centre + np.asarray(translation, float)
    return out


def scaled_rigid(direction, size_mm, points, centre):
    """A rigid transform of RMS displacement ``size_mm`` over ``points``, along ``direction``.

    ``direction`` is a 6-vector (3 shifts in mm, 3 rotations in degrees about ``centre``), scaled
    by bisection; the displacement grows with the scale over the sizes used here.
    """
    direction = np.asarray(direction, float)

    def at(scale):
        return rigid(direction[:3] * scale, direction[3:] * scale, centre)

    def size(scale):
        return distance(np.eye(4), at(scale), points)

    lo, hi = 0.0, 1.0
    while size(hi) < size_mm:
        hi *= 2

    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if size(mid) < size_mm else (lo, mid)
    return at(hi)


def mask_points(mask_file, every=1):
    """World coordinates of every ``every``-th voxel inside a mask."""
    img = nb.load(mask_file)
    ijk = np.argwhere(np.asanyarray(img.dataobj) > 0)[::every]
    return apply_affine(img.affine, ijk)


class CarriedMask:
    """Dice of the ASL brain mask carried into T1w space through a pull, against the T1w mask.

    Evaluated on the T1w grid at a stride (2 voxels by default) inside the T1w mask's bounding box
    grown by ``margin_mm``, so masks that fall outside the T1w brain are counted.
    """

    def __init__(self, asl_mask_file, t1w_mask_file, stride=2, margin_mm=20.0):
        asl = nb.load(asl_mask_file)
        self.asl = np.asanyarray(asl.dataobj) > 0
        self.asl_inv = np.linalg.inv(asl.affine)
        t1 = nb.load(t1w_mask_file)
        t1_data = np.asanyarray(t1.dataobj) > 0
        vox = np.sqrt((t1.affine[:3, :3] ** 2).sum(0))
        grow = np.ceil(margin_mm / vox).astype(int)
        inside = np.argwhere(t1_data)
        lo = np.maximum(inside.min(0) - grow, 0)
        hi = np.minimum(inside.max(0) + grow + 1, t1_data.shape)
        axes = [np.arange(a, b, stride) for a, b in zip(lo, hi, strict=True)]
        ijk = np.stack(np.meshgrid(*axes, indexing='ij'), -1).reshape(-1, 3)
        self.t1 = t1_data[tuple(ijk.T)]
        self.points = apply_affine(t1.affine, ijk)
        self.centre = apply_affine(t1.affine, inside.mean(0))

    def inside(self, pull):
        """Which grid points the pull maps into the ASL mask (nearest voxel)."""
        ijk = np.rint(apply_affine(self.asl_inv @ pull, self.points)).astype(int)
        ok = np.all((ijk >= 0) & (ijk < self.asl.shape), axis=1)
        out = np.zeros(len(ijk), bool)
        out[ok] = self.asl[tuple(ijk[ok].T)]
        return out

    def dice(self, pull):
        carried = self.inside(pull)
        return float(2 * (carried & self.t1).sum() / (carried.sum() + self.t1.sum()))

    def perturbations(self, kind, size_mm):
        """The calibration's perturbations, applied on the T1w side of a pull.

        ``axes``: shifts of ``size_mm`` along each T1w axis, both ways (the plan's wording).
        ``rigid``: ``CALIBRATION_DRAWS`` random rigid transforms (fixed seed) about the brain's
        centre, each of RMS displacement ``size_mm`` over the T1w brain, the kind of error a
        registration makes.
        """
        if kind == 'axes':
            shifts = [sign * size_mm * axis for axis in np.eye(3) for sign in (-1, 1)]
            return [rigid(translation=t) for t in shifts]

        rng = np.random.default_rng(CALIBRATION_SEED)
        brain = self.points[self.t1]
        return [
            scaled_rigid(rng.normal(size=6), size_mm, brain, self.centre)
            for _ in range(CALIBRATION_DRAWS)
        ]

    def calibration(self, pull, shift_mm=CALIBRATION_SHIFT_MM, kind='axes'):
        """Median Dice drop when ``pull`` is perturbed by ``shift_mm`` (see ``perturbations``)."""
        best = self.dice(pull)
        drops = [best - self.dice(pull @ p) for p in self.perturbations(kind, shift_mm)]
        return float(np.median(drops))


def gross_by_mask(carried, arms, kind='axes'):
    """Per arm: Dice, the run's best Dice, the calibrated limit and whether the arm fails.

    ``arms`` is {arm: pull}. The calibration is the best arm's own 5 mm drop (``kind``: see
    ``CarriedMask.perturbations``).
    """
    dice = {arm: carried.dice(pull) for arm, pull in arms.items()}
    best_arm = max(dice, key=dice.get)
    limit = carried.calibration(arms[best_arm], kind=kind)
    return {
        arm: {
            'dice': d,
            'best_dice': dice[best_arm],
            'best_arm': best_arm,
            'limit': limit,
            'gross_mask': bool(dice[best_arm] - d > limit),
        }
        for arm, d in dice.items()
    }


def gross_verdict(partner_values, mask_gross_rigid, limit=GROSS_INCONSISTENCY_MM):
    """The gross-failure rule adopted for Tasks 3-6 (2026-10-29): 'fail', 'ambiguous' or 'pass'.

    ``partner_values`` are the run's inconsistencies with each comparable partner run. With two
    or more partners, the run fails when it exceeds ``limit`` with every partner, so a good run
    beside a failed one is not flagged. With one partner, an inconsistency above ``limit`` cannot
    say which run failed, so both are 'ambiguous'. With none, the carried-mask rule with the rigid
    calibration (``mask_gross_rigid``) decides; it is a backstop, weaker at 6-8 mm.
    """
    values = list(partner_values)
    if len(values) >= 2:
        return 'fail' if all(v > limit for v in values) else 'pass'

    if len(values) == 1:
        return 'ambiguous' if values[0] > limit else 'pass'

    return 'fail' if mask_gross_rigid else 'pass'
