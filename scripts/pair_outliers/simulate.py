"""Simulate multi-delay PCASL control-label pairs with known corrupted pairs.

The phantom is the MNI152NLin2009cAsym tissue probability maps, block-averaged to 4 mm.
The signal model is deliberately richer than the general kinetic model that ASLPrep fits,
so that clean pairs have the systematic misfit real data have:

- the labeled bolus is dispersed by a gamma kernel on its way to the voxel,
- label in tissue relaxes with tissue T1 rather than blood T1,
- every image carries physiological noise proportional to its static tissue signal
  (a global scale and a smooth spatial field), on top of slightly smoothed thermal noise.

Corrupted pairs get one of four artifacts, each with a random severity:

``motion``
    The label image is moved by a small rigid transform relative to the control image,
    as residual motion after head-motion correction leaves it.
``intensity``
    For 2D data, a slab of one to three slices of one image is scaled (dropout or spin
    history). For 3D data, a whole image is scaled (drift or a receiver gain change).
``labeling``
    The perfusion signal of the pair is reduced, as when labeling fails.
``spike``
    The pair's thermal noise is multiplied several-fold.
"""

from dataclasses import dataclass, field

import nibabel as nb
import numpy as np
from scipy import ndimage, special
from templateflow import api as tflow

ARTIFACTS = ('motion', 'intensity', 'labeling', 'spike')

LABELEFF = 0.85
T1BLOOD = 1.65
PARTITION = 0.9
M0_SCALE = 1000.0


@dataclass
class Design:
    """An acquisition protocol."""

    name: str
    plds: list
    repeats: list
    tau: float
    readout: str  # '2d' or '3d'
    background_suppression: bool
    slice_time: float = 0.0  # seconds per slice for 2D readouts
    thermal_sd: float = 1.75  # per image, in M0 units of 1000
    notes: str = ''
    pld_per_pair: np.ndarray = field(init=False)

    def __post_init__(self):
        self.pld_per_pair = np.repeat(np.asarray(self.plds, float), self.repeats)


DESIGNS = {
    'hcpa': Design(
        'hcpa',
        plds=[0.2, 0.7, 1.2, 1.7, 2.2],
        repeats=[6, 6, 6, 10, 15],
        tau=1.5,
        readout='2d',
        background_suppression=False,
        slice_time=0.0125,
        notes='HCP-Aging-like 2D multiband PCASL without background suppression',
    ),
    'rep4': Design(
        'rep4',
        plds=[0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        repeats=[4] * 6,
        tau=1.8,
        readout='3d',
        background_suppression=True,
        notes='3D background-suppressed PCASL, six delays with four pairs each',
    ),
    'rep2': Design(
        'rep2',
        plds=[0.5, 1.0, 1.5, 2.0, 2.5],
        repeats=[2] * 5,
        tau=1.5,
        readout='3d',
        background_suppression=True,
        notes='3D background-suppressed PCASL, five delays with two pairs each',
    ),
    'rep1': Design(
        'rep1',
        plds=[1.25, 1.5, 2.0, 2.5, 3.0],
        repeats=[1] * 5,
        tau=1.8,
        readout='3d',
        background_suppression=True,
        notes='3D background-suppressed PCASL, five delays with one pair each',
    ),
}


def _smooth_field(shape, rng, sigma, mask=None):
    """Make a smooth random field with zero mean and unit standard deviation."""
    noise = ndimage.gaussian_filter(rng.standard_normal(shape), sigma)
    values = noise[mask] if mask is not None else noise
    return (noise - values.mean()) / values.std()


class Phantom:
    """Tissue maps and ground-truth parameters on a 4 mm grid."""

    def __init__(self, seed=0):
        probs = []
        for label in ('GM', 'WM', 'CSF'):
            path = tflow.get('MNI152NLin2009cAsym', resolution=2, label=label, suffix='probseg')
            img = nb.load(path)
            data = img.get_fdata()
            # Block-average 2 mm to 4 mm.
            shape = tuple(s // 2 for s in data.shape)
            data = data[: shape[0] * 2, : shape[1] * 2, : shape[2] * 2]
            data = data.reshape(shape[0], 2, shape[1], 2, shape[2], 2).mean(axis=(1, 3, 5))
            probs.append(data)
        self.gm, self.wm, self.csf = probs
        self.shape = self.gm.shape
        self.affine = img.affine.copy()
        self.affine[:3, :3] *= 2
        tissue = self.gm + self.wm
        self.head = (tissue + self.csf) > 0.3
        self.brain = tissue > 0.5
        self.gm_mask = self.gm > 0.6

        rng = np.random.default_rng(seed)
        head = self.head
        # M0 (proton density) with a smooth receive field
        bias = 1 + 0.15 * _smooth_field(self.shape, rng, 6, head)
        self.m0 = M0_SCALE * (0.85 * self.gm + 0.72 * self.wm + 1.0 * self.csf) * bias
        self.m0[~head] = 0.0
        # Tissue T1 at 3T (seconds) and static fraction left by background suppression
        t1 = 1.45 * self.gm + 0.85 * self.wm + 4.0 * self.csf
        self.t1 = np.where(head, t1 / np.maximum(self.gm + self.wm + self.csf, 1e-3), 1.0)
        self.bs_residual = np.where(head, 0.06 * self.gm + 0.03 * self.wm + 0.12 * self.csf, 0.0)

        # Ground truth: CBF (mL/100 g/min), ATT (s), and macrovascular aBV and aBAT
        cbf = (62 * self.gm + 22 * self.wm) * (1 + 0.15 * _smooth_field(self.shape, rng, 3, head))
        self.cbf = np.where(head, np.maximum(cbf, 0.0), 0.0)
        # ATT is longer in white matter and toward the back of the brain.
        y = np.linspace(1, 0, self.shape[1])[None, :, None]
        att = 1.0 * self.gm + 1.5 * self.wm + 0.35 * y + 0.1 * _smooth_field(self.shape, rng, 3)
        self.att = np.where(head, np.clip(att / np.maximum(self.gm + self.wm, 0.3), 0.6, 2.6), 1.0)
        vessels = (rng.random(self.shape) < 0.02) & head
        self.abv = np.where(vessels, rng.uniform(0.003, 0.012, self.shape), 0.0)
        self.abat = np.clip(self.att - rng.uniform(0.2, 0.5, self.shape), 0.3, None)


def _gamma_cdf(t, shape=2.0, scale=0.1):
    t = np.maximum(t, 0.0)
    return special.gammainc(shape, t / scale)


def true_deltam(phantom, plds, tau, extra_delay=None):
    """Delta-M from a dispersed-bolus, two-T1 kinetic model, for each delay.

    ``plds`` is a sequence of delays; ``extra_delay`` (3D array) shifts them voxelwise,
    as slice timing does. Returns an array of shape (*phantom.shape, n_delays).
    """
    head = phantom.head
    cbf, att, t1t = phantom.cbf[head], phantom.att[head], phantom.t1[head]
    abv, abat, m0 = phantom.abv[head], phantom.abat[head], phantom.m0[head]
    shift = extra_delay[head] if extra_delay is not None else np.zeros_like(cbf)
    dt = 0.005
    out = np.zeros(phantom.shape + (len(plds),))
    for i_pld, pld in enumerate(plds):
        t_read = tau + pld + shift  # time since the start of labeling, voxelwise
        read_step = np.rint(t_read / dt).astype(int)
        tissue = np.zeros_like(cbf)
        result = np.zeros_like(cbf)
        decay = np.exp(-dt / t1t)
        for step in range(read_step.max() + 1):
            # Midpoint rule: label arriving during this step, decayed for half a step
            t = (step - 0.5) * dt
            arriving = np.exp(-att / T1BLOOD) * (_gamma_cdf(t - att) - _gamma_cdf(t - att - tau))
            tissue = tissue * decay + arriving * dt * np.sqrt(decay) * (step > 0)
            result = np.where(step == read_step, tissue, result)
        arterial = np.exp(-abat / T1BLOOD) * (
            _gamma_cdf(t_read - abat) - _gamma_cdf(t_read - abat - tau)
        )
        signal = 2 * LABELEFF * (m0 / PARTITION) * (cbf / 6000 * result + abv * arterial)
        vol = np.zeros(phantom.shape)
        vol[head] = signal
        out[..., i_pld] = vol
    return out


def static_signal(phantom, design, pld, slice_delay):
    """Static tissue signal of a control or label image at one delay."""
    if design.background_suppression:
        return phantom.m0 * phantom.bs_residual
    # Presaturation at the start of labeling, then recovery until the slice is read out.
    recovery = 1 - np.exp(-(design.tau + pld + slice_delay) / phantom.t1)
    return phantom.m0 * recovery


def _rigid(volume, rng, magnitude_mm, voxel_mm=4.0):
    """Move a volume by a random rigid transform of the given magnitude (mm and degrees)."""
    angles = np.deg2rad(rng.normal(0, 1, 3))
    angles *= np.deg2rad(magnitude_mm) / np.linalg.norm(angles)
    shift = rng.normal(0, 1, 3)
    shift *= magnitude_mm / np.linalg.norm(shift) / voxel_mm
    cx, cy, cz = np.cos(angles)
    sx, sy, sz = np.sin(angles)
    rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    rot = rz @ ry @ rx
    center = (np.array(volume.shape) - 1) / 2
    offset = center - rot @ center + shift
    return ndimage.affine_transform(volume, rot, offset=offset, order=1, mode='constant')


def simulate_run(phantom, design, n_bad, seed, artifacts=ARTIFACTS):
    """Simulate one run's delta-M series.

    Returns a dict with the 4D delta-M (artifacts included), the same series without
    artifacts, per-voxel delays, and a table of corrupted pairs.
    """
    rng = np.random.default_rng(seed)
    shape = phantom.shape
    n_pairs = design.pld_per_pair.size
    n_slices = shape[2]
    if design.readout == '2d':
        slice_delay = (np.arange(n_slices) * design.slice_time)[None, None, :] * np.ones(shape)
    else:
        slice_delay = np.zeros(shape)

    signal = true_deltam(phantom, design.plds, design.tau, extra_delay=slice_delay)
    pld_index = np.repeat(np.arange(len(design.plds)), design.repeats)

    bad = np.zeros(n_pairs, dtype=bool)
    bad[rng.choice(n_pairs, size=n_bad, replace=False)] = True
    kinds = np.full(n_pairs, '', dtype=object)
    severity = np.zeros(n_pairs)

    head = phantom.head
    clean = np.zeros(shape + (n_pairs,))
    dirty = np.zeros(shape + (n_pairs,))
    for i_pair in range(n_pairs):
        pld = design.pld_per_pair[i_pair]
        static = static_signal(phantom, design, pld, slice_delay)
        dm = signal[..., pld_index[i_pair]]
        images = []
        for _ in range(2):
            # Physiological noise: a global scale and a smooth field, both proportional to
            # the static signal, and thermal noise smoothed slightly by reconstruction.
            physio = 0.001 * rng.standard_normal() + 0.0015 * _smooth_field(shape, rng, 2)
            thermal = ndimage.gaussian_filter(rng.standard_normal(shape), 0.5)
            thermal *= design.thermal_sd / thermal.std()
            images.append((physio, thermal))

        def make_pair(images, static, dm, thermal_gain=1.0, label_gain=1.0):
            (p_c, n_c), (p_l, n_l) = images
            control = static * (1 + p_c) + thermal_gain * n_c
            label = static * (1 + p_l) - label_gain * dm + thermal_gain * n_l
            control[~head] = thermal_gain * n_c[~head]
            label[~head] = thermal_gain * n_l[~head]
            return control, label

        control, label = make_pair(images, static, dm)
        clean[..., i_pair] = control - label
        if not bad[i_pair]:
            dirty[..., i_pair] = control - label
            continue

        kind = artifacts[rng.integers(len(artifacts))]
        kinds[i_pair] = kind
        if kind == 'motion':
            sev = rng.uniform(0.3, 1.5)  # mm (and degrees)
            label = _rigid(label, rng, sev)
        elif kind == 'intensity':
            sev = rng.uniform(0.02, 0.15) if design.readout == '2d' else rng.uniform(0.01, 0.05)
            target = control if rng.random() < 0.5 else label
            if design.readout == '2d':
                width = rng.integers(1, 4)
                start = rng.integers(5, n_slices - 5 - width)
                target[..., start : start + width] *= 1 - sev
            else:
                target *= 1 + sev * rng.choice([-1, 1])
        elif kind == 'labeling':
            sev = rng.uniform(0.3, 1.0)  # fraction of the perfusion signal lost
            control, label = make_pair(images, static, dm, label_gain=1 - sev)
        else:  # spike
            sev = rng.uniform(2.0, 5.0)  # noise gain
            control, label = make_pair(images, static, dm, thermal_gain=sev)
        severity[i_pair] = sev
        dirty[..., i_pair] = control - label

    plds = design.pld_per_pair[None, :] + slice_delay[head][:, None]
    return {
        'dirty': dirty[head],
        'clean': clean[head],
        'truth': signal[head][:, pld_index],
        'plds': plds,
        'bad': bad,
        'kinds': kinds,
        'severity': severity,
        'pld_per_pair': design.pld_per_pair,
    }
