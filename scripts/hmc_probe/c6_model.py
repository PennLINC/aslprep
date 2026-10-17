"""C6 prototype: a differentiable contrast model fitted jointly with the head motion (PyTorch).

The model is the one that passed the identifiability test (``c6_identifiability.py``; plan,
"C6: differentiable contrast model"). Each control or label volume ``v`` of delay ``d`` and group
``g = (type, d)`` is predicted, in its own contrast, as

    y_v(p) = R_v(T_v p),   R_v = exp(sum_k a_gk B_k) * S  -  w_d * P * 1(label)

where ``p`` is a point of the volume (world mm), ``T_v`` a rigid transform (six parameters, the
mean pose fixed), ``S`` the static image, ``B_k`` the ten polynomials up to order 2, ``P`` one
perfusion image that is zero outside a brain support, and ``w_d`` a weight per delay (the first
fixed at 1). Least squares against the data, with no cross-contrast cost and no fixed reference.

The fit is Adam on the GPU over a smoothing pyramid (data and prediction smoothed alike), after a
first stage with the poses frozen. The result is a RAS pull per volume (reference point to volume
point), the static image (the reference handed on) and the model maps.
"""

import numpy as np
import torch
import torch.nn.functional as F


def poly_basis(shape, device):
    u = [torch.linspace(-1, 1, n, device=device) for n in shape]
    x, y, z = torch.meshgrid(*u, indexing='ij')
    one = torch.ones_like(x)
    return torch.stack([one, x, y, z, x * x, y * y, z * z, x * y, x * z, y * z])


def rigid(params, centre):
    """(n, 6) [tx ty tz mm, rx ry rz degrees] -> (n, 4, 4) world transforms about ``centre``."""
    t = params[:, :3]
    r = torch.deg2rad(params[:, 3:])
    c, s = torch.cos(r), torch.sin(r)
    one, zero = torch.ones_like(c[:, 0]), torch.zeros_like(c[:, 0])
    rx = torch.stack([one, zero, zero, zero, c[:, 0], -s[:, 0], zero, s[:, 0], c[:, 0]], 1)
    ry = torch.stack([c[:, 1], zero, s[:, 1], zero, one, zero, -s[:, 1], zero, c[:, 1]], 1)
    rz = torch.stack([c[:, 2], -s[:, 2], zero, s[:, 2], c[:, 2], zero, zero, zero, one], 1)
    R = rz.view(-1, 3, 3) @ ry.view(-1, 3, 3) @ rx.view(-1, 3, 3)
    T = torch.zeros(params.shape[0], 4, 4, dtype=params.dtype, device=params.device)
    T[:, :3, :3] = R
    T[:, :3, 3] = centre + t - (R @ centre)
    T[:, 3, 3] = 1.0
    return T


def params_from_transforms(T, centre):
    """Inverse of ``rigid``: (n, 4, 4) world transforms -> (n, 6) parameters (R = Rz Ry Rx)."""
    T = np.asarray(T, float)
    out = np.zeros((T.shape[0], 6))
    for v, m in enumerate(T):
        R = m[:3, :3]
        ry = np.arcsin(np.clip(-R[2, 0], -1, 1))
        rx = np.arctan2(R[2, 1], R[2, 2])
        rz = np.arctan2(R[1, 0], R[0, 0])
        out[v, 3:] = np.degrees([rx, ry, rz])
        out[v, :3] = m[:3, 3] - centre + R @ centre
    return out


def kcrop(x, shape, pe_weights=None, pe_axis=1):
    """Acquire an in-plane oversampled image: keep the central ``shape[:2]`` of its 2D k-space.

    ``x`` is (n, 1, oX, oY, Z); returns the magnitude on (n, 1, X, Y, Z), at the same intensity
    scale. Even sizes only (the centred crop is then exact). ``pe_weights`` (one per acquired
    phase-encoding line, centred order) multiply k-space along ``pe_axis`` (0 or 1): the
    readout's T2* decay, the same for every volume.
    """
    X2, Y2 = x.shape[2:4]
    X, Y = shape[:2]
    k = torch.fft.fftshift(torch.fft.fft2(x, dim=(2, 3)), dim=(2, 3))
    sx, sy = (X2 - X) // 2, (Y2 - Y) // 2
    k = k[:, :, sx : sx + X, sy : sy + Y]
    if pe_weights is not None:
        view = [1, 1, 1, 1, 1]
        view[2 + pe_axis] = -1
        k = k * pe_weights.view(*view)
    out = torch.fft.ifft2(torch.fft.ifftshift(k, dim=(2, 3)), dim=(2, 3)) * (X * Y / (X2 * Y2))
    return out.abs()


def kpad(img, o):
    """The inverse: zero-pad an (X, Y, Z) image's 2D k-space by ``o`` in-plane (numpy)."""
    X, Y = img.shape[:2]
    k = np.fft.fftshift(np.fft.fft2(img, axes=(0, 1)), axes=(0, 1))
    pad = (((o - 1) * X // 2,) * 2, ((o - 1) * Y // 2,) * 2, (0, 0))
    k = np.pad(k, pad)
    return np.real(np.fft.ifft2(np.fft.ifftshift(k, axes=(0, 1)), axes=(0, 1))) * o * o


def gaussian(x, sigma_vox):
    """Separable Gaussian smoothing of (n, 1, X, Y, Z) tensors; ``sigma_vox`` per axis."""
    if max(sigma_vox) <= 0:
        return x
    for axis, s in enumerate(sigma_vox):
        if s <= 0:
            continue
        r = int(np.ceil(3 * s))
        k = torch.exp(-0.5 * (torch.arange(-r, r + 1, device=x.device, dtype=x.dtype) / s) ** 2)
        view = [1, 1, 1, 1, 1]
        view[2 + axis] = -1
        k = (k / k.sum()).view(*view)
        pad = [0, 0, 0, 0, 0, 0]
        pad[2 * (2 - axis)] = pad[2 * (2 - axis) + 1] = r
        x = F.conv3d(F.pad(x, pad, mode='replicate'), k)
    return x


def fit(
    data,
    ctx,
    pld,
    affine,
    support,
    init=None,
    stages=((None, 200), (4.0, 300), (2.0, 300), (0.0, 400)),
    device='cuda',
    seed=0,
    freeze_poses=False,
    hires=1,
    pe_axis=None,
    p_bound=None,
):
    """Fit the model; ``data`` (X, Y, Z, n), ``support`` (X, Y, Z) bool, ``init`` (n, 6) or None.

    ``stages``: (FWHM in mm or None, iterations); None freezes the poses (first stage, unsmoothed).
    ``freeze_poses`` keeps the poses at ``init`` throughout (the model's floor at given poses).
    ``hires``: the maps live on a grid ``hires`` times finer in-plane, and each predicted volume
    is acquired by cropping its 2D k-space to the data's matrix (2D readouts: partial volume
    and truncation move with the head, which interpolating a coarse image cannot reproduce).
    ``pe_axis`` (0 or 1, with ``hires`` > 1): also learn a real weight per phase-encoding line,
    shared by all volumes (the readout's T2* filter).
    ``p_bound``: bound the perfusion image to ``±p_bound`` times the initial static image
    (``P = p_bound S0 tanh(q)``). Delta-M is about 1% of the control signal, while absorbing a
    label displacement into ``P`` needs edge dipoles several times larger; an inequality the
    Fisher test cannot use.
    """
    torch.manual_seed(seed)
    dt = torch.float32
    shape = data.shape[:3]
    n = data.shape[-1]
    ctx = np.asarray(ctx)
    pld = np.asarray(pld, float)
    delays = sorted(set(pld))
    groups = [(t, d) for d in delays for t in ('control', 'label')]
    gidx = torch.tensor(
        [groups.index((t, d)) for t, d in zip(ctx, pld, strict=True)], device=device
    )
    d_of = torch.tensor([delays.index(d) for d in pld], device=device)
    is_label = torch.tensor(ctx == 'label', device=device, dtype=dt)

    ctrl = data[..., ctx == 'control']
    scale = float(np.median(ctrl.mean(-1)[support]))
    y = torch.tensor(np.moveaxis(data / scale, -1, 0)[:, None], dtype=dt, device=device)
    head = torch.tensor(data.mean(-1) > 0.05 * data.mean(-1).max(), device=device)

    # Initial maps from the unaligned data.
    S0 = ctrl.mean(-1) / scale
    perf = np.mean(
        [
            data[..., (ctx == 'control') & (pld == d)].mean(-1)
            - data[..., (ctx == 'label') & (pld == d)].mean(-1)
            for d in delays
        ],
        axis=0,
    )
    P0 = np.where(support, perf, 0.0) / scale
    a0 = [
        float(
            np.log(
                np.median(data[..., (ctx == t) & (pld == d)].mean(-1)[support])
                / scale
                / np.median(S0[support])
            )
        )
        for t, d in groups
    ]
    o = int(hires)
    hshape = (shape[0] * o, shape[1] * o, shape[2])
    if o > 1:
        S0, P0 = kpad(S0, o), kpad(P0, o)
        support = support.repeat(o, 0).repeat(o, 1)
    S = torch.tensor(S0, dtype=dt, device=device, requires_grad=True)
    supp = torch.tensor(support, dtype=dt, device=device)
    if p_bound is None:
        P = torch.tensor(P0, dtype=dt, device=device, requires_grad=True)

        def perfusion():
            return P * supp
    else:
        cap = torch.tensor(p_bound * np.abs(S0) + 1e-6, dtype=dt, device=device)
        q0 = np.arctanh(np.clip(P0 / cap.cpu().numpy(), -0.95, 0.95))
        P = torch.tensor(q0, dtype=dt, device=device, requires_grad=True)

        def perfusion():
            return cap * torch.tanh(P) * supp

    a = torch.zeros(len(groups), 10, dtype=dt, device=device)
    a[:, 0] = torch.tensor(a0, dtype=dt)
    a.requires_grad_(True)
    w_free = torch.zeros(max(len(delays) - 1, 0), dtype=dt, device=device, requires_grad=True)
    phi = torch.tensor(np.zeros((n, 6)) if init is None else init, dtype=dt, device=device)
    phi.requires_grad_(True)
    B = poly_basis(hshape, device).to(dt)

    A = torch.tensor(affine, dtype=dt, device=device)
    # Map grid indices <-> acquired-grid indices: in-plane, i = (i_h - (o - 1) / 2) / o.
    to_acq = torch.tensor([1 / o, 1 / o, 1.0], dtype=dt, device=device)
    off = torch.tensor([(o - 1) / (2 * o), (o - 1) / (2 * o), 0.0], dtype=dt, device=device)
    ijk = torch.stack(
        torch.meshgrid(*[torch.arange(s, device=device, dtype=dt) for s in hshape], indexing='ij'),
        -1,
    )
    world = (ijk * to_acq - off) @ A[:3, :3].T + A[:3, 3]  # (oX, oY, Z, 3)
    centre = (A @ torch.tensor([*((np.array(shape) - 1) / 2), 1.0], dtype=dt, device=device))[:3]
    A_inv = torch.linalg.inv(A)
    dims = torch.tensor(hshape, dtype=dt, device=device)
    zooms = np.sqrt((affine[:3, :3] ** 2).sum(0))

    pe = None
    if o > 1 and pe_axis is not None:
        pe = torch.ones(shape[pe_axis], dtype=dt, device=device, requires_grad=True)

    def acquire(x):
        return kcrop(x, shape, pe, pe_axis or 1) if o > 1 else x

    def weights():
        w = torch.cat([torch.ones(1, dtype=dt, device=device), 1.0 + w_free])
        return w[d_of]

    def predict(theta):
        T = rigid(theta, centre)  # volume world -> reference world
        pts = torch.einsum('nij,xyzj->nxyzi', T[:, :3, :3], world) + T[:, None, None, None, :3, 3]
        vox = (pts @ A_inv[:3, :3].T + A_inv[:3, 3] + off) / to_acq  # fine-grid indices
        norm = 2 * vox / (dims - 1) - 1
        grid = norm.flip(-1)  # grid_sample wants (z, y, x) for (D=X, H=Y, W=Z)
        contrast = torch.exp(torch.einsum('gk,kxyz->gxyz', a, B))[gidx]  # (n, X, Y, Z)
        ref = contrast * S - (weights() * is_label)[:, None, None, None] * perfusion()
        moved = F.grid_sample(
            ref[:, None], grid, mode='bilinear', padding_mode='border', align_corners=True
        )
        return acquire(moved)

    history = []
    for fwhm, iters in stages:
        freeze = fwhm is None or freeze_poses
        level = 0.0 if fwhm is None else fwhm
        sig = [level / 2.3548 / z for z in zooms]
        y_s = gaussian(y, sig)
        opt = torch.optim.Adam(
            [{'params': [S, P], 'lr': 5e-3}, {'params': [a, w_free], 'lr': 5e-3}]
            + ([{'params': [pe], 'lr': 2e-3}] if pe is not None else [])
            + ([] if freeze else [{'params': [phi], 'lr': 2e-2}])
        )
        # Cosine decay to a tenth within each stage, so the last steps do not jitter at the
        # learning rate (Adam's steps stay about that size near an optimum).
        sched = torch.optim.lr_scheduler.LambdaLR(
            opt, lambda t, n_it=iters: 0.1 + 0.45 * (1 + np.cos(np.pi * min(t, n_it) / n_it))
        )
        first = None
        for _ in range(iters):
            opt.zero_grad()
            theta = phi - phi.mean(0)
            pred = gaussian(predict(theta), sig)
            loss = ((pred - y_s)[:, 0][:, head] ** 2).mean()
            loss.backward()
            opt.step()
            sched.step()
            first = float(loss) if first is None else first
        history.append({'fwhm_mm': fwhm, 'iters': iters, 'loss': float(loss), 'first': first})

    with torch.no_grad():
        theta = (phi - phi.mean(0)).cpu().double()
        T = rigid(theta, centre.cpu().double()).numpy()
        pulls = np.linalg.inv(T)  # reference point -> volume point
        w = torch.cat([torch.ones(1, dtype=dt, device=device), 1.0 + w_free]).cpu().numpy()
        return {
            'pulls': pulls,
            'theta': theta.numpy(),
            'S': acquire(S[None, None]).cpu().numpy()[0, 0] * scale,
            'P': acquire(perfusion()[None, None]).cpu().numpy()[0, 0] * scale,
            'w': dict(zip(map(float, delays), map(float, w), strict=True)),
            'history': history,
        }
