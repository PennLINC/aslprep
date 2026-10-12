"""Geometry for the aslscan ground-truth tests.

Conventions (see Sections 4.4 and 4.6 of
``docs/specs/2026-10-08-aslscan-ground-truth-tests-design.md``):

- Matrices are 4x4 and act on RAS world points in mm, unless a name says LPS.
- A *point mapping* ``M`` sends a point ``p`` to ``M @ p``.
- An ITK *image transform* file named ``from-X_to-Y`` stores a *pull* mapping:
  it sends points of Y's space to X's space (which is what resampling X onto Y needs).
- ITK files store LPS matrices; ``ras_to_lps`` is a change of basis, not an inverse.

Only numpy is required, so the fixture generator can use this module.
"""

from pathlib import Path

import numpy as np

_FLIP = np.diag([-1.0, -1.0, 1.0, 1.0])


def _rot_x(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def _rot_y(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def _rot_z(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def rotation_zyx(rx, ry, rz):
    """``Rz @ Ry @ Rx`` for active right-handed rotations given in radians.

    This is ``rotation_zyx_deg`` in mrsim-acq (``src/mat.rs``), in radians.
    """
    return _rot_z(rz) @ _rot_y(ry) @ _rot_x(rx)


def rigid_ras(tx, ty, tz, rx, ry, rz):
    """A rigid point mapping: rotate about the world origin (degrees), then translate (mm)."""
    m = np.eye(4)
    m[:3, :3] = rotation_zyx(*np.radians([rx, ry, rz]))
    m[:3, 3] = [tx, ty, tz]
    return m


def pose_matrix(trans_mm, rot_rad, center):
    """aslscan's head pose for one volume, as a point mapping of the static phantom.

    mrsim-acq applies ``p' = R (p - c) + c + t`` with ``R = Rz Ry Rx``
    (``Pose::to_matrix``, ``src/motion.rs``), about the simulation grid's field-of-view centre
    ``c`` (``fov_center``). ``desc-motion_gt.tsv`` gives translations in mm and rotations in
    radians.
    """
    r = rotation_zyx(*rot_rad)
    c = np.asarray(center, dtype=float)
    m = np.eye(4)
    m[:3, :3] = r
    m[:3, 3] = c - r @ c + np.asarray(trans_mm, dtype=float)
    return m


def fov_center(affine, shape):
    """World coordinates of the centre of a grid (mrsim-acq ``fov_center``)."""
    ijk = (np.asarray(shape[:3], dtype=float) - 1) / 2
    return affine[:3, :3] @ ijk + affine[:3, 3]


def ras_to_lps(m):
    """Express a RAS point mapping in LPS coordinates (and back: the change is an involution)."""
    return _FLIP @ m @ _FLIP


lps_to_ras = ras_to_lps


def write_itk_affine(path, m_ras):
    """Write a RAS point mapping as an ITK affine transform text file (centre at the origin)."""
    m = ras_to_lps(np.asarray(m_ras, dtype=float))
    params = [*m[:3, :3].ravel(), *m[:3, 3]]
    text = (
        '#Insight Transform File V1.0\n'
        '#Transform 0\n'
        'Transform: AffineTransform_double_3_3\n'
        f'Parameters: {" ".join(repr(float(v)) for v in params)}\n'
        'FixedParameters: 0 0 0\n'
    )
    Path(path).write_text(text)


def rot_angle_deg(m):
    """Rotation angle (degrees) of a rigid matrix's linear part."""
    r = np.asarray(m)[:3, :3]
    # atan2 of (sin, cos) stays accurate near the identity, where arccos(cos) does not
    sin = np.linalg.norm([r[2, 1] - r[1, 2], r[0, 2] - r[2, 0], r[1, 0] - r[0, 1]]) / 2
    cos = (np.trace(r) - 1) / 2
    return float(np.degrees(np.arctan2(sin, cos)))


def rms_displacement(m, points):
    """Root-mean-square distance (mm) that ``m`` moves ``points`` (an N x 3 array)."""
    points = np.asarray(points, dtype=float)
    moved = points @ np.asarray(m)[:3, :3].T + np.asarray(m)[:3, 3]
    return float(np.sqrt(np.mean(np.sum((moved - points) ** 2, axis=1))))


def apply_points(m, points):
    """Apply a 4x4 point mapping to an N x 3 array of points."""
    points = np.asarray(points, dtype=float)
    return points @ np.asarray(m)[:3, :3].T + np.asarray(m)[:3, 3]


def _axis_edges(affine, shape, axis):
    """World coordinates of the cell edges of an axis-aligned grid along one axis."""
    step = affine[axis, axis]
    first = affine[axis, 3] - step / 2
    return first + step * np.arange(shape[axis] + 1)


def _check_axis_aligned(affine, what):
    lin = np.asarray(affine)[:3, :3]
    if np.any(np.abs(lin - np.diag(np.diag(lin))) > 1e-6) or np.any(np.diag(lin) <= 0):
        raise ValueError(
            f'{what} must be axis-aligned RAS with positive voxel sizes; affine is\n{affine}'
        )


def overlap_matrix(src_edges, dst_edges):
    """1D overlap weights: ``W[a, p]`` = length of source cell ``p`` inside target cell ``a``,
    divided by the target cell's length.

    The weights of a target cell sum to the fraction of it that the source covers. This is the
    per-axis factor of aslscan's volume-weighted mean (``src/resample.rs``).
    """
    src_edges = np.asarray(src_edges, dtype=float)
    dst_edges = np.asarray(dst_edges, dtype=float)
    lo = np.maximum(dst_edges[:-1, None], src_edges[None, :-1])
    hi = np.minimum(dst_edges[1:, None], src_edges[None, 1:])
    return np.clip(hi - lo, 0, None) / np.diff(dst_edges)[:, None]


def overlap_mean(data, src_affine, dst_affine, dst_shape):
    """Volume-weighted mean of a 3D array onto an axis-aligned target grid.

    Both grids must be axis-aligned RAS. Applied to a one-hot label map this gives
    partial-volume fractions; applied to a phantom map it reproduces aslscan's ``*_gt`` maps.
    """
    _check_axis_aligned(src_affine, 'source grid')
    _check_axis_aligned(dst_affine, 'target grid')
    out = np.asarray(data, dtype=np.float64)
    for axis in range(3):
        w = overlap_matrix(
            _axis_edges(src_affine, out.shape, axis),
            _axis_edges(dst_affine, dst_shape, axis),
        )
        out = np.moveaxis(np.tensordot(w, out, axes=([1], [axis])), 0, axis)
    return out


def overlap_fractions(labels, src_affine, dst_affine, dst_shape, label_values=(1, 2, 3)):
    """Partial-volume fraction of each label on the target grid, as a dict of arrays."""
    labels = np.asarray(labels)
    return {
        value: overlap_mean(
            (labels == value).astype(np.float32), src_affine, dst_affine, dst_shape
        )
        for value in label_values
    }
