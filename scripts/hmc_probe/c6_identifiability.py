#!/usr/bin/env python3
"""C6 identifiability test (plan: Task 6), before any optimizer is written.

The C6 model predicts each volume ``v`` of group ``g = (type, delay)`` as

    y_v = W_{theta_v} [ c_g * S  -  w_d(v) * P ]

with ``S`` a static-tissue image, ``c_g`` a contrast field per group, ``P`` a perfusion image
(entering label volumes only, scaled by a weight per delay), and ``W_theta`` a rigid warp with
six parameters per volume. The global frame is anchored by fixing the mean pose.

The test is local and exact: the Fisher information of all parameters at the truth, with the
nuisance parameters (S, P, contrast, delay weights) eliminated by the Schur complement, gives the
Cramer-Rao bound on any pattern of poses. A pattern that the nuisance parameters can absorb has
no information left (an infinite bound). The patterns are those the plan names: one volume
(reference), one delay block moved, every label moved (motion confounded with type) and the labels
of one delay moved.

Models (contrast / perfusion):

- ``free``: an unrestricted contrast map per group and a free perfusion map per delay. The plan
  states this is not identifiable; the test must find that, or it is not testing anything.
- ``smooth``: a smooth contrast field per group (polynomials up to order 2, ten terms, applied as
  ``exp(sum a_k B_k)``), one perfusion map with a weight per delay. The proposed model.
- ``smooth-pd``: the smooth contrast with a free perfusion map per delay.
- ``gain``: a scalar gain per group, one perfusion map.
- ``+supp-s`` / ``+supp-t``: P fixed to zero outside a support, the static-image mask or the true
  grey and white matter (dilated by a voxel; the best case). Added after the first results: a
  free P absorbs label shifts proportional to the delay weights whenever the tissue contrast is
  the same at every delay (an equality constraint can enter the bound; positivity cannot).

Criterion, fixed before the test was run: for the proposed model, every confounded pattern has a
standard deviation of at most 0.1 mm (translations) and 0.1 degrees (rotations) at white noise of
2% of the median brain control signal, a typical raw-image level. (The first statement said "at
the fixture's noise"; the fixtures turned out to be nearly noiseless, 0.01%, which would make the
criterion empty, so the noise is set explicitly. That change was made after the first run had
shown only that the free model is not identifiable, before any noise-level result.) Positivity of
P and other inequality priors cannot enter a Fisher bound; they can only help, so passing without
them is sufficient and failing is reported as such. White noise is optimistic for real data,
whose physiological noise is structured; 5% is reported as a stress level.
"""

import argparse
import json
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
from scipy import ndimage

AXES = ('tx', 'ty', 'tz', 'rx', 'ry', 'rz')


def load(fixture):
    perf = fixture / 'sub-01' / 'perf'
    img = nb.load(next(perf.glob('*_asl.nii.gz')))
    data = np.asanyarray(img.dataobj, np.float64)
    ctx = pd.read_table(next(perf.glob('*_aslcontext.tsv')))['volume_type'].to_numpy()
    meta = json.loads(next(perf.glob('*_asl.json')).read_text())
    pld = np.atleast_1d(np.asarray(meta.get('PostLabelingDelay', 0.0), float))
    pld = np.repeat(pld, len(ctx)) if pld.size == 1 else pld
    keep = np.isin(ctx, ('control', 'label'))
    return img, data[..., keep], ctx[keep], pld[keep]


def tissue_mask(fixture):
    """The true grey plus white matter (partial volume at least 0.5), dilated by one voxel."""
    gt = fixture / 'sub-01' / 'perf' / 'ground-truth'
    pv = sum(
        np.asanyarray(nb.load(gt / f'sub-01_desc-pv{t}_gt.nii.gz').dataobj, np.float64)
        for t in ('GM', 'WM')
    )
    return ndimage.binary_dilation(pv >= 0.5)


def build_state(img, data, ctx, pld):
    """The linearization point: group means, S, contrast, P and delay weights."""
    delays = sorted(set(pld))
    groups = [(t, d) for d in delays for t in ('control', 'label')]
    gidx = np.array([groups.index((t, d)) for t, d in zip(ctx, pld, strict=True)])
    ctrl_mean = {d: data[..., (ctx == 'control') & (pld == d)].mean(-1) for d in delays}
    lab_mean = {d: data[..., (ctx == 'label') & (pld == d)].mean(-1) for d in delays}
    S = np.mean([ctrl_mean[d] for d in delays], axis=0)
    mask = S > 0.2 * np.percentile(S, 99)
    # Tissue contrast is shared by the control and label of one delay; labels add -w_d P.
    # Outside the brain the contrast is 1 (no artificial edge at the mask), where S is ~0 anyway.
    head = S > 1e-3 * S.max()
    c = {d: np.where(head, ctrl_mean[d] / np.where(head, S, 1.0), 1.0) for d in delays}
    P_d = {d: np.where(mask, ctrl_mean[d] - lab_mean[d], 0.0) for d in delays}
    P = np.mean([P_d[d] for d in delays], axis=0)
    w = {d: float((P_d[d] * P).sum() / (P * P).sum()) for d in delays}
    # Pooled within-group noise from repeats of each group.
    resid = []
    for g in range(len(groups)):
        sel = data[..., gidx == g]
        if sel.shape[-1] > 1:
            resid.append(
                (sel - sel.mean(-1, keepdims=True))[mask].ravel()
                * np.sqrt(sel.shape[-1] / (sel.shape[-1] - 1))
            )
    sigma = float(np.std(np.concatenate(resid))) if resid else np.nan
    return {
        'delays': delays,
        'groups': groups,
        'gidx': gidx,
        'S': S,
        'mask': mask,
        'c': c,
        'P': P,
        'P_d': P_d,
        'w': w,
        'sigma': sigma,
        'img': img,
        'ctx': ctx,
        'pld': pld,
    }


def predicted(st):
    """Each volume's image at the truth (no motion in the frame of the linearization)."""
    out = []
    for t, d in zip(st['ctx'], st['pld'], strict=True):
        y = st['c'][d] * st['S']
        if t == 'label':
            y = y - st['w'][d] * st['P']
        out.append(y)
    return np.stack(out, -1)


def pose_jacobian(img, pred):
    """d y_v / d theta_v for each volume: (n_vol, 6, N) in mm and degrees."""
    A = img.affine[:3, :3]
    ijk = np.indices(img.shape[:3]).reshape(3, -1)
    xyz = A @ ijk + img.affine[:3, 3:4]
    r = xyz - xyz.mean(1, keepdims=True)
    inv_t = np.linalg.inv(A).T
    out = np.empty((pred.shape[-1], 6, ijk.shape[1]))
    for v in range(pred.shape[-1]):
        g_vox = np.stack(np.gradient(pred[..., v]), 0).reshape(3, -1)
        g = inv_t @ g_vox  # world gradient, per mm
        out[v, :3] = g
        for k in range(3):
            e = np.zeros(3)
            e[k] = 1.0
            disp = np.cross(e[:, None], r, axis=0) * (np.pi / 180)  # mm per degree
            out[v, 3 + k] = (g * disp).sum(0)
    return out


def poly_basis(shape):
    """Polynomials up to order 2 in normalized voxel coordinates: (10, N)."""
    u = [np.linspace(-1, 1, n) for n in shape]
    x, y, z = (a.ravel() for a in np.meshgrid(*u, indexing='ij'))
    return np.stack([np.ones_like(x), x, y, z, x * x, y * y, z * z, x * y, x * z, y * z])


def information(st, model, chunk=4000):
    """Pose information matrix (6 n_vol square) with every nuisance parameter eliminated.

    ``model`` is a contrast/perfusion model, optionally with ``+supp-s`` (P fixed to zero outside
    the static-image mask) or ``+supp-t`` (outside the true tissue mask, the best case).
    """
    model, _, supp = model.partition('+')
    support = {'': None, 'supp-s': st['mask'], 'supp-t': st.get('tissue')}[supp]
    p_on = np.ones(st['S'].size) if support is None else support.ravel().astype(float)
    img, ctx, pld, delays = st['img'], st['ctx'], st['pld'], st['delays']
    n = len(ctx)
    N = int(np.prod(img.shape[:3]))
    pred = predicted(st)
    Jp = pose_jacobian(img, pred)  # (n, 6, N)
    flat = {k: st[k].ravel() for k in ('S', 'P')}
    c = {d: st['c'][d].ravel() for d in delays}
    gidx = st['gidx']
    is_label = ctx == 'label'
    d_of = np.array([delays.index(p) for p in pld])
    n_groups = len(st['groups'])

    # Global small nuisance parameters, built per chunk: (n, q, nc).
    B = None
    if model in ('smooth', 'smooth-pd'):
        B = poly_basis(img.shape[:3])
    elif model == 'gain':
        B = np.ones((1, N))
    n_w = len(delays) - 1 if model in ('smooth', 'gain') else 0  # the first w_d is fixed
    q = (n_groups * B.shape[0] if B is not None else 0) + n_w

    def small(sl):
        out = np.zeros((n, q, len(sl)))
        if B is not None:
            nb_ = B.shape[0]
            for v in range(n):
                g = gidx[v]
                out[v, g * nb_ : (g + 1) * nb_] = B[:, sl] * (c[pld[v]][sl] * flat['S'][sl])
        for j in range(1, n_w + 1):
            for v in np.where(is_label & (d_of == j))[0]:
                out[v, q - n_w + j - 1] = -flat['P'][sl]
        return out

    # Per-voxel nuisance parameters: their derivatives (n, m) at each voxel.
    def local(sl):
        cols = [np.stack([c[pld[v]][sl] for v in range(n)])]  # S
        if model == 'free':
            cols.extend(
                np.stack([flat['S'][sl] * (gidx[v] == g) for v in range(n)])
                for g in range(n_groups)
            )
        if model in ('free', 'smooth-pd'):
            cols.extend(
                np.stack([-float(is_label[v] and d_of[v] == j) * p_on[sl] for v in range(n)])
                for j in range(len(delays))
            )
        else:
            cols.append(
                np.stack([-(st['w'][pld[v]] if is_label[v] else 0.0) * p_on[sl] for v in range(n)])
            )
        return np.stack(cols, -1)  # (n, nc, m)

    # Parameters: theta (6n), then small (q). F over [theta, small], minus per-voxel Schur terms.
    dim = 6 * n + q
    F = np.zeros((dim, dim))
    support = np.flatnonzero(ndimage.binary_dilation(st['S'] > 1e-3 * st['S'].max(), iterations=2))
    for start in range(0, support.size, chunk):
        sl = support[start : start + chunk]
        nc = len(sl)
        # Global Jacobian on this chunk: (n, dim, nc); a pose column is non-zero in its volume.
        G = np.zeros((n, dim, nc))
        for v in range(n):
            G[v, 6 * v : 6 * v + 6] = Jp[v][:, sl]
        G[:, 6 * n :] = small(sl)
        F += np.einsum('vix,vjx->ij', G, G, optimize=True)
        K = local(sl)  # (n, nc, m)
        Fll = np.einsum('vxa,vxb->xab', K, K, optimize=True)
        Fgl = np.einsum('vix,vxa->xia', G, K, optimize=True)
        Fll_pinv = np.linalg.pinv(Fll, rcond=1e-10, hermitian=True)
        D = np.einsum('xia,xab->xib', Fgl, Fll_pinv, optimize=True)
        F -= np.einsum('xib,xjb->ij', D, Fgl, optimize=True)
    F /= st['sigma'] ** 2
    # Eliminate the small nuisance parameters, then anchor the mean pose.
    Ftt, Fts, Fss = F[: 6 * n, : 6 * n], F[: 6 * n, 6 * n :], F[6 * n :, 6 * n :]
    if q:
        Ftt = Ftt - Fts @ np.linalg.pinv(Fss, rcond=1e-10, hermitian=True) @ Fts.T
    return Ftt


def anchored(Ftt, n):
    """Eigen-decomposition of the pose information with the mean pose fixed.

    Returns the eigenvalues and eigenvectors of the information restricted to the zero-mean
    subspace, and the tolerance below which an eigenvalue counts as no information.
    """
    # An orthonormal basis of the zero-mean subspace (per axis, the last n - 1 columns of a
    # centred QR), so the anchored directions are removed exactly rather than found numerically.
    q, _ = np.linalg.qr(np.eye(n) - np.ones((n, n)) / n)
    Q = np.kron(q[:, : n - 1], np.eye(6))
    lam, W = np.linalg.eigh(Q.T @ Ftt @ Q)
    return lam, Q @ W, 1e-10 * lam.max()


def patterns(ctx, pld):
    delays = sorted(set(pld))
    n = len(ctx)
    pats = {'one volume': np.arange(n) == n // 2}
    if len(delays) > 1:
        mid = delays[len(delays) // 2]
        pats['one delay block'] = pld == mid
        pats['labels of one delay'] = (pld == mid) & (ctx == 'label')
    pats['every label'] = ctx == 'label'
    return pats


def bound(lam, V, tol, sel, axis, n):
    """SD of the displacement amplitude of the volumes in ``sel`` relative to the rest.

    Infinite when the pattern reaches a direction with no information. The eigenvectors span the
    zero-mean subspace only, so the anchored directions never enter.
    """
    u = np.zeros((n, 6))
    u[sel, axis] = 1.0
    u = (u - u.mean(0)).ravel()  # relative to the anchored mean
    lvec = u / (u @ u)
    proj = V.T @ lvec
    null = lam <= tol
    if np.linalg.norm(proj[null]) > 1e-3 * np.linalg.norm(lvec):
        return np.inf
    return float(np.sqrt((proj[~null] ** 2 / lam[~null]).sum()))


def run(fixture, models, noise_fracs):
    img, data, ctx, pld = load(fixture)
    st = build_state(img, data, ctx, pld)
    st['tissue'] = tissue_mask(fixture)
    s_med = float(np.median(st['S'][st['mask']]))
    n = len(ctx)
    rows = []
    for model in models:
        st['sigma'] = 1.0  # information per unit variance; scaled below
        Ftt = information(st, model)
        for frac in noise_fracs:
            sigma = frac / 100 * s_med
            lam, V, tol = anchored(Ftt / sigma**2, n)
            for name, sel in patterns(ctx, pld).items():
                sds = [bound(lam, V, tol, sel, a, n) for a in range(6)]
                rows.append(
                    {
                        'fixture': fixture.name,
                        'model': model,
                        'noise_pct': frac,
                        'pattern': name,
                        'n_volumes': int(sel.sum()),
                        'sigma': sigma,
                        **{f'sd_{ax}': s for ax, s in zip(AXES, sds, strict=True)},
                        'max_sd_trans_mm': max(sds[:3]),
                        'max_sd_rot_deg': max(sds[3:]),
                        'n_null': int((lam <= tol).sum()),
                    }
                )
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--fixtures', type=Path, nargs='+', required=True)
    parser.add_argument(
        '--models',
        nargs='+',
        default=[
            'free',
            'smooth',
            'smooth-pd',
            'gain',
            'smooth+supp-s',
            'smooth+supp-t',
            'gain+supp-t',
        ],
    )
    parser.add_argument(
        '--noise-pct',
        type=float,
        nargs='+',
        default=[0.5, 1.0, 2.0, 5.0],
        help='White noise SD, as a percentage of the median brain control signal.',
    )
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    rows = []
    for fx in args.fixtures:
        rows += run(fx, args.models, args.noise_pct)
        print(fx.name, 'done', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    with pd.option_context('display.width', 200, 'display.max_rows', 300):
        print(
            df[
                [
                    'fixture',
                    'model',
                    'noise_pct',
                    'pattern',
                    'max_sd_trans_mm',
                    'max_sd_rot_deg',
                    'n_null',
                ]
            ]
            .round(4)
            .to_string(index=False)
        )


if __name__ == '__main__':
    main()
