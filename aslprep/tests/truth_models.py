"""Independent ASL quantification for the Tier A conformance checks.

These functions recompute, from raw simulated images, what ASLPrep's documented model should
produce (Section 4.7 of ``docs/specs/2026-10-08-aslscan-ground-truth-tests-design.md``).
They are written from the cited papers, not from ASLPrep's code, and deliberately do not
import from ``aslprep``. Where ASLPrep's behaviour is a convention rather than physics
(for example its labeling-efficiency rule), the function says so.

References
----------
Alsop et al. (2015), Magn Reson Med 73:102-116 (the ASL white paper), equations 1 and 2.
Wong et al. (1998), Magn Reson Med 39:702-708 (QUIPSS II).
Luh et al. (1999), Magn Reson Med 41:1246-1254 (Q2TIPS).
"""

import numpy as np

UNIT_CONVERSION = 6000.0  # ml/g/s -> ml/100 g/min
LAMBDA = 0.9  # blood-brain partition coefficient, ml/g (Alsop 2015)
#: Arterial blood T1 by field strength (Alsop 2015; Zhang et al. 2013 for 7 T).
T1_BLOOD = {1.5: 1.35, 3: 1.65, 7: 2.087}
#: Grey-matter T1 for the short-TR M0 correction (Wright et al. 2008), ASLPrep's constants.
T1_TISSUE = {1.5: 1.197, 3: 1.607, 7: 1.939}
#: Labeling efficiency by labeling type (Alsop 2015; Wang et al. 2005 for CASL).
BASE_EFFICIENCY = {'PCASL': 0.85, 'PASL': 0.98, 'CASL': 0.68}
#: Fraction of the ASL signal each background-suppression pulse keeps (Alsop et al. 2015).
BS_PULSE_EFFICIENCY = 0.95


def volume_indices(context):
    """Indices of each volume type in an aslcontext list."""
    return {
        kind: [i for i, t in enumerate(context) if t == kind]
        for kind in ('control', 'label', 'deltam', 'm0scan', 'cbf')
    }


def deltam_mean(asl, context, swap=False):
    """Mean delta-M (control - label, paired by order; or deltam rows).

    ``swap`` subtracts the other way round (a mutation for the sensitivity tests).
    """
    idx = volume_indices(context)
    if idx['deltam']:
        return asl[..., idx['deltam']].mean(axis=-1)
    control, label = asl[..., idx['control']], asl[..., idx['label']]
    diff = label - control if swap else control - label
    return diff.mean(axis=-1)


def smooth(data, affine, fwhm):
    """Gaussian smoothing in mm (nibabel's implementation, as ASLPrep uses)."""
    if not fwhm:
        return data
    import nibabel as nb
    from nibabel.processing import smooth_image

    return smooth_image(nb.Nifti1Image(data, affine), fwhm=fwhm).get_fdata()


def m0_image(asl, context, metadata, affine, m0scan=None, m0_metadata=None, fwhm=0.0):
    """The calibration image and its repetition time, by ``M0Type`` (BIDS ASL).

    Returns ``(m0, tr)``; ``tr`` is None when no TR correction applies (``Estimate``).
    ``m0`` is the M0 of tissue, which the equations below divide by lambda. BIDS defines
    ``M0Estimate`` as the M0 of arterial blood, so the tissue M0 is lambda times it.
    """
    idx = volume_indices(context)
    m0_type = metadata['M0Type']
    if m0_type == 'Separate':
        data = m0scan if m0scan.ndim == 3 else m0scan.mean(axis=-1)
        return smooth(data, affine, fwhm), np.mean(m0_metadata['RepetitionTimePreparation'])
    if m0_type == 'Included':
        data = smooth(asl[..., idx['m0scan']], affine, fwhm).mean(axis=-1)
        return data, _tr_of(metadata, idx['m0scan'])
    if m0_type == 'Estimate':
        return np.full(asl.shape[:3], LAMBDA * float(metadata['M0Estimate'])), None
    if m0_type == 'Absent':
        if metadata.get('BackgroundSuppression'):
            raise ValueError('background-suppressed control volumes cannot calibrate CBF')
        data = smooth(asl[..., idx['control']], affine, fwhm).mean(axis=-1)
        return data, _tr_of(metadata, idx['control'])
    raise ValueError(f'unknown M0Type {m0_type!r}')


def _tr_of(metadata, rows):
    tr = metadata['RepetitionTimePreparation']
    return float(np.mean(np.asarray(tr)[rows])) if isinstance(tr, list) else float(tr)


def m0_tr_correction(m0, tr, field_strength):
    """Correct an M0 acquired with TR < 5 s for incomplete recovery (Alsop 2015, p. 113)."""
    if tr is None or tr >= 5:
        return m0
    return m0 / (1 - np.exp(-tr / T1_TISSUE[field_strength]))


def labeling_efficiency(metadata):
    """ASLPrep's rule: the sidecar value, else the base value, x 0.95^n for n suppression pulses.

    Each pulse keeps about 95 % of the ASL signal (Alsop et al. 2015). The sidecar's
    LabelingEfficiency describes the labeling alone, so the loss applies to it as well.
    """
    if 'LabelingEfficiency' in metadata:
        alpha = float(metadata['LabelingEfficiency'])
    else:
        alpha = BASE_EFFICIENCY[metadata['ArterialSpinLabelingType']]
    if metadata.get('BackgroundSuppression'):
        alpha *= BS_PULSE_EFFICIENCY ** metadata.get('BackgroundSuppressionNumberPulses', 1)
    return alpha


def labeling_efficiency_physical(simulation, context=None):
    """The efficiency the simulator applied: labeling efficiency x suppression attenuation.

    ``simulation`` is the sidecar's ``AslscanSimulation`` block. The suppression factor is a
    scalar, or one value per volume (M0 volumes are unsuppressed); per-volume factors are taken
    from the label and delta-M volumes of ``context``, which must agree.
    """
    alpha = simulation['Resolved']['LabelingEfficiency']['Value']
    factor = simulation.get('BackgroundSuppressionLabelFactor')
    if factor is None:
        return alpha
    if isinstance(factor, list):
        if context is None or len(context) != len(factor):
            raise ValueError('a per-volume suppression factor needs the matching aslcontext')
        values = {
            round(abs(f), 9)
            for f, t in zip(factor, context, strict=True)
            if t in ('label', 'deltam')
        }
        if len(values) != 1:
            raise ValueError(f'suppression factors differ across labeled volumes: {values}')
        factor = values.pop()
    return alpha * abs(factor)


def labeled_value(metadata, key, context=None):
    """A single-delay timing field (``PostLabelingDelay``, ``LabelingDuration``).

    A scalar is returned as is. A per-volume list is read on the control, label and deltam
    rows only (M0 rows carry placeholders such as 0), and must have one value there.
    """
    values = np.atleast_1d(np.asarray(metadata[key], dtype=float))
    if values.size == 1:
        return float(values[0])
    if context is None or len(context) != values.size:
        raise ValueError(f'a per-volume {key} needs the matching aslcontext')
    rows = [i for i, t in enumerate(context) if t in ('control', 'label', 'deltam')]
    unique = np.unique(values[rows])
    if unique.size != 1:
        raise ValueError(f'{key} differs across labeled volumes: {unique.tolist()}')
    return float(unique[0])


def pld_map(metadata, shape, context=None, reverse_slices=False, pld_shift=0.0):
    """Post-labeling delay (PCASL) or inversion time (PASL) per voxel, with slice timing.

    2D acquisitions add each slice's acquisition time (BIDS ``SliceTiming``, along the
    ``SliceEncodingDirection`` axis, reversed for ``-``). Mutations: ``reverse_slices`` and
    ``pld_shift``.
    """
    pld = labeled_value(metadata, 'PostLabelingDelay', context) + pld_shift
    out = np.full(shape, pld)
    timing = metadata.get('SliceTiming')
    if timing is None:
        return out
    direction = metadata.get('SliceEncodingDirection', 'k')
    axis = 'ijk'.index(direction[0])
    timing = np.asarray(timing, dtype=float)
    if direction.endswith('-') != reverse_slices:
        timing = timing[::-1]
    expand = [None, None, None]
    expand[axis] = slice(None)
    return out + timing[tuple(expand)]


def cbf_single_delay(deltam, m0, metadata, alpha, pld, context=None):
    """Single-delay CBF (ml/100 g/min) from delta-M and M0, per Alsop 2015.

    PCASL/CASL (eq. 1): ``6000 lambda dM e^(PLD/T1b) / (2 alpha T1b M0 (1 - e^(-tau/T1b)))``.
    PASL with QUIPSS II or Q2TIPS (eq. 2): ``6000 lambda dM e^(TI/T1b) / (2 alpha TI1 M0)``,
    where ``TI`` is the inversion time (BIDS ``PostLabelingDelay`` for PASL) and ``TI1`` the
    bolus duration (the first ``BolusCutOffDelayTime``).
    """
    t1b = T1_BLOOD[metadata['MagneticFieldStrength']]
    kind = metadata['ArterialSpinLabelingType']
    if kind in ('PCASL', 'CASL'):
        tau = labeled_value(metadata, 'LabelingDuration', context)
        denominator = 2 * alpha * t1b * (1 - np.exp(-tau / t1b)) * m0
    elif kind == 'PASL':
        technique = metadata.get('BolusCutOffTechnique')
        if technique not in ('QUIPSSII', 'Q2TIPS'):
            raise ValueError(f'PASL bolus cut-off {technique!r} is not modelled here')
        ti1 = np.atleast_1d(metadata['BolusCutOffDelayTime'])[0]
        denominator = 2 * alpha * ti1 * m0
    else:
        raise ValueError(f'unknown labeling type {kind!r}')
    with np.errstate(divide='ignore', invalid='ignore'):
        return UNIT_CONVERSION * LAMBDA * deltam * np.exp(pld / t1b) / denominator


def expected_cbf(
    asl,
    context,
    metadata,
    affine,
    m0scan=None,
    m0_metadata=None,
    m0_scale=1.0,
    fwhm=0.0,
    alpha=None,
    mutation=None,
):
    """Tier A expectation for single-delay data, from raw images.

    ``mutation`` applies one deliberate error (for the sensitivity tests): ``'swap'``,
    ``'reverse_slices'``, ``'pld_shift'``, ``'m0_scale'`` or ``'efficiency_ignored'``.
    """
    deltam = deltam_mean(asl, context, swap=mutation == 'swap')
    m0, tr = m0_image(asl, context, metadata, affine, m0scan, m0_metadata, fwhm)
    m0 = m0_tr_correction(m0, tr, metadata['MagneticFieldStrength']) * m0_scale
    if mutation == 'm0_scale':
        m0 = m0 * 1.1
    if alpha is None:
        alpha = 1.0 if mutation == 'efficiency_ignored' else labeling_efficiency(metadata)
    pld = pld_map(
        metadata,
        asl.shape[:3],
        context,
        reverse_slices=mutation == 'reverse_slices',
        pld_shift=0.1 if mutation == 'pld_shift' else 0.0,
    )
    return cbf_single_delay(deltam, m0, metadata, alpha, pld, context)


# ---------------------------------------------------------------------------------------------
# Multi-delay
# ---------------------------------------------------------------------------------------------
#: Bounds ASLPrep documents for its general-kinetic-model fit (utils/cbf.py,
#: fit_deltam_multipld), as (lower, upper) for (CBF, ATT, aBAT, aBV).
FIT_BOUNDS = ((0.0, 0.0, 0.0, 0.0), (300.0, 5.0, 5.0, 0.1))
#: ATT spacing of the reference fit's exhaustive search, in seconds.
FIT_ATT_STEP = 0.002
#: Fits whose sums of squared errors differ by less than this fraction of the data's sum of
#: squares are treated as equally good.
FIT_TIE_TOLERANCE = 1e-9


def gkm_pcasl(pld, tau, cbf, att, abat, abv, alpha, t1b, m0a, m0b):
    """(P)CASL delta-M: Buxton tissue term plus a plug-flow arterial term.

    Tissue (Buxton et al. 1998; Woods et al. 2023, eq. 2), with ``f = cbf / 6000``:
    0 before arrival; ``2 a M0a f T1b e^(-ATT/T1b) (1 - e^(-(tau + PLD - ATT)/T1b))`` while
    the bolus arrives; ``2 a M0a f T1b e^(-PLD/T1b) (1 - e^(-tau/T1b))`` after.
    Arterial (Chappell et al. 2010): ``2 a M0b aBV e^(-aBAT/T1b)`` while the bolus is in the
    vessel (``aBAT <= tau + PLD < aBAT + tau``), else 0.
    """
    end = tau + pld
    f = cbf / UNIT_CONVERSION
    arriving = 2 * alpha * m0a * f * t1b * np.exp(-att / t1b) * (1 - np.exp(-(end - att) / t1b))
    arrived = 2 * alpha * m0a * f * t1b * np.exp(-pld / t1b) * (1 - np.exp(-tau / t1b))
    tissue = np.where(end < att, 0.0, np.where(end < att + tau, arriving, arrived))
    in_vessel = (abat <= end) & (end < abat + tau)
    arterial = np.where(in_vessel, 2 * alpha * m0b * abv * np.exp(-abat / t1b), 0.0)
    return tissue + arterial


def gkm_pasl(ti, ti1, cbf, att, abat, abv, alpha, t1b, m0a, m0b):
    """PASL delta-M with a bolus of duration ``ti1`` (Buxton et al. 1998; Woods et al. 2023).

    Tissue: 0 before ``ATT``; ``2 a M0a f e^(-TI/T1b) (TI - ATT)`` while arriving;
    ``2 a M0a f e^(-TI/T1b) TI1`` after. Arterial: ``2 a M0b aBV e^(-TI/T1b)`` for
    ``aBAT < TI < aBAT + TI1``.
    """
    f = cbf / UNIT_CONVERSION
    decay = np.exp(-ti / t1b)
    tissue = np.where(
        ti < att,
        0.0,
        np.where(
            ti < att + ti1,
            2 * alpha * m0a * f * decay * (ti - att),
            2 * alpha * m0a * f * decay * ti1,
        ),
    )
    in_vessel = (abat < ti) & (ti < abat + ti1)
    arterial = np.where(in_vessel, 2 * alpha * m0b * abv * decay, 0.0)
    return tissue + arterial


def tissue_t1_signal(t, att, bolus, t1b, t1prime, casl):
    """Tissue delta-M per unit ``2 a M0a f`` with tissue relaxation, by numerical integration.

    The label delivered to the tissue at time ``s`` (from the start of labeling) has decayed with
    blood T1 on its way, and then decays with the tissue's apparent T1 until ``t``
    (Buxton et al. 1998), so the signal is the integral over ``s`` of
    ``c(s) exp(-(t - s) / T1')``. For (P)CASL, label leaves the labeling plane at a steady rate
    for ``bolus`` (tau) seconds and arrives ``att`` later, so ``c(s) = exp(-att / T1b)`` for
    ``att <= s < att + tau``; ``t`` is ``tau + PLD``. For PASL, the whole bolus is labeled at
    once, so blood arriving at ``s`` has decayed for ``s`` seconds:
    ``c(s) = exp(-s / T1b)`` for ``att <= s < att + TI1``; ``t`` is TI.
    This does not use the closed-form expressions, so it checks them independently.
    """
    from scipy.integrate import quad

    end = min(t, att + bolus)
    if end <= att:
        return 0.0

    def integrand(s):
        delivered = np.exp(-att / t1b) if casl else np.exp(-s / t1b)
        return delivered * np.exp(-(t - s) / t1prime)

    value, _ = quad(integrand, att, end, epsabs=0, epsrel=1e-12)
    return value


def deltam_observations(asl, context, swap=False):
    """Every delta-M observation (pairs by order, or deltam rows) and its sidecar row index."""
    idx = volume_indices(context)
    if idx['deltam']:
        return asl[..., idx['deltam']], idx['deltam']
    control, label = asl[..., idx['control']], asl[..., idx['label']]
    diff = label - control if swap else control - label
    return diff, idx['control']


def delay_edges(plds, bolus, casl):
    """The ATT and aBAT values at which a delay's signal changes form, within the fit bounds.

    A delay's tissue term changes form where ATT, and its arterial term where aBAT, crosses the
    delay or the delay plus the labeling duration ((P)CASL), or the inversion time minus TI1 or
    the inversion time (PASL). The bounds are included, so consecutive edges tile them.
    """
    lower, upper = FIT_BOUNDS[0][1], FIT_BOUNDS[1][1]
    edges = np.concatenate([plds, plds + bolus] if casl else [plds - bolus, plds])
    return np.unique(np.clip(np.concatenate([edges, [lower, upper]]), lower, upper))


def _bounded_pair(y, tissue, arterial):
    """Bounded least squares for (CBF, aBV) at every ATT, by enumerating the candidates.

    ``tissue`` is (n_att, n_obs) and ``arterial`` is (n_obs,), both per unit of their
    parameter. The objective is convex, so its bounded minimum is the unconstrained solution
    when that is in bounds, and otherwise lies on an edge of the box, where the other parameter
    is its clipped one-parameter solution. Returns the CBF, aBV and sum of squared errors of the
    best candidate at each ATT; among equally good candidates, the one with the smallest aBV.
    """
    (c_lo, _, _, v_lo), (c_hi, _, _, v_hi) = FIT_BOUNDS
    tt = np.einsum('ao,ao->a', tissue, tissue)
    aa = arterial @ arterial
    ta = tissue @ arterial
    ty = tissue @ y
    ay = arterial @ y
    with np.errstate(divide='ignore', invalid='ignore'):
        det = tt * aa - ta**2
        free_c = (aa * ty - ta * ay) / det
        free_v = (tt * ay - ta * ty) / det
        free_ok = (
            (det > 1e-10 * tt * aa)
            & (free_c >= c_lo)
            & (free_c <= c_hi)
            & (free_v >= v_lo)
            & (free_v <= v_hi)
        )

        def best_c(v):
            return np.where(tt > 0, np.clip((ty - v * ta) / tt, c_lo, c_hi), c_lo)

        def best_v(c):
            return np.clip((ay - c * ta) / aa, v_lo, v_hi) if aa > 0 else np.full(c.shape, v_lo)

    n_att = tissue.shape[0]
    candidates = [
        (best_c(v_lo), np.full(n_att, v_lo)),
        (np.full(n_att, c_lo), best_v(np.full(n_att, c_lo))),
        (np.full(n_att, c_hi), best_v(np.full(n_att, c_hi))),
        (np.where(free_ok, free_c, c_lo), np.where(free_ok, free_v, v_lo)),
        (best_c(v_hi), np.full(n_att, v_hi)),
    ]
    cbf = np.stack([c for c, _ in candidates])
    abv = np.stack([v for _, v in candidates])
    residual = y - cbf[..., None] * tissue - abv[..., None] * arterial
    sse = (residual**2).sum(axis=-1)
    # Smallest aBV among the candidates within the tie tolerance of the best one
    tied = sse <= sse.min(axis=0) + FIT_TIE_TOLERANCE * (y @ y)
    pick = np.argmin(np.where(tied, abv, np.inf), axis=0)
    cols = np.arange(n_att)
    return cbf[pick, cols], abv[pick, cols], sse[pick, cols]


def fit_multi_delay(deltam, m0, plds, metadata, alpha, rows=None):
    """Fit CBF, ATT, aBAT and aBV per voxel by exhaustive bounded least squares.

    ``deltam`` and ``plds`` are (n_voxels, n_observations); ``m0`` is the calibration image
    (already TR-corrected and scaled), (n_voxels,). ``rows`` are the observations' sidecar
    rows, needed when ``LabelingDuration`` is per volume.

    This is the model ASLPrep documents, minimized without ASLPrep's search strategy: every
    observation is used as is, ATT is searched on a :data:`FIT_ATT_STEP` grid over its bounds
    (plus the edges of :func:`delay_edges`), and CBF and aBV are solved at each candidate.
    Four conventions follow ASLPrep's documentation, because the data cannot settle them:
    ``M0b`` equals ``M0a = M0 / lambda``; aBAT is the midpoint of its interval between
    consecutive edges; among fits within :data:`FIT_TIE_TOLERANCE` of the best, the one with
    the smallest ATT, then the smallest aBAT, then the smallest aBV is selected; and where the
    ATTs that follow it within its interval fit equally well, the middle of that run is reported.
    A voxel without a positive, finite M0 or with non-finite data gives NaN.
    """
    t1b = T1_BLOOD[metadata['MagneticFieldStrength']]
    casl = metadata['ArterialSpinLabelingType'] in ('PCASL', 'CASL')
    n_obs = deltam.shape[1]
    if casl:
        bolus = np.atleast_1d(np.asarray(metadata['LabelingDuration'], float))
        if bolus.size > 1:
            if rows is None:
                raise ValueError('a per-volume LabelingDuration needs the observation rows')
            bolus = bolus[rows]
        bolus = np.broadcast_to(bolus, (n_obs,))
        model = gkm_pcasl
    else:
        bolus = float(np.atleast_1d(metadata['BolusCutOffDelayTime'])[0])
        model = gkm_pasl
    lower, upper = FIT_BOUNDS[0][1], FIT_BOUNDS[1][1]
    att_grid = np.arange(lower, upper + FIT_ATT_STEP / 2, FIT_ATT_STEP)

    out = np.full((deltam.shape[0], 4), np.nan)
    for i in range(deltam.shape[0]):
        m0a = m0[i] / LAMBDA
        y = deltam[i]
        if not (np.isfinite(m0a) and m0a > 0 and np.all(np.isfinite(y))):
            continue
        edges = delay_edges(plds[i], bolus, casl)
        att = np.unique(np.concatenate([att_grid, edges]))
        abat = (edges[:-1] + edges[1:]) / 2
        # Signal per unit CBF at each ATT, and per unit aBV at each aBAT
        tissue = model(plds[i], bolus, 1.0, att[:, None], 0.0, 0.0, alpha, t1b, m0a, m0a)
        arterial = model(plds[i], bolus, 0.0, 0.0, abat[:, None], 1.0, alpha, t1b, m0a, m0a)
        # (n_att, n_abat) arrays, so that raveling orders the candidates by ATT, then aBAT
        cbf, abv, sse = (
            np.stack(columns, axis=1)
            for columns in zip(*(_bounded_pair(y, tissue, a) for a in arterial), strict=True)
        )
        tied = sse <= sse.min() + FIT_TIE_TOLERANCE * (y @ y)
        a, b = np.unravel_index(np.argmax(tied), tied.shape)
        # The run of equally good ATTs from there to the end of its interval
        next_edge = min(np.searchsorted(edges, att[a], side='right'), edges.size - 1)
        run = tied[a:, b] & (att[a:] <= edges[next_edge])
        run_length = run.size if run.all() else int(np.argmin(run))
        a += (run_length - 1) // 2
        out[i] = cbf[a, b], att[a], abat[b], abv[a, b]
    return out


def expected_multi_delay(
    asl,
    context,
    metadata,
    affine,
    voxels,
    m0scan=None,
    m0_metadata=None,
    m0_scale=1.0,
    mutation=None,
    fwhm=0.0,
):
    """Tier A expectation for multi-delay data on selected voxels (an (n, 3) index array).

    Returns an (n, 4) array of CBF, ATT, aBAT and aBV. ``mutation`` as in
    :func:`expected_cbf` (``'swap'``, ``'pld_shift'``, ``'m0_scale'``). ``fwhm`` smooths M0 as
    ASLPrep does.
    """
    obs, rows = deltam_observations(asl, context, swap=mutation == 'swap')
    m0, tr = m0_image(asl, context, metadata, affine, m0scan, m0_metadata, fwhm=fwhm)
    m0 = m0_tr_correction(m0, tr, metadata['MagneticFieldStrength']) * m0_scale
    if mutation == 'm0_scale':
        m0 = m0 * 1.1
    plds = np.asarray(metadata['PostLabelingDelay'], float)[rows]
    if mutation == 'pld_shift':
        plds = plds + 0.1
    offset = pld_map({**metadata, 'PostLabelingDelay': 0.0}, asl.shape[:3])
    i, j, k = voxels.T
    voxel_plds = plds[None, :] + offset[i, j, k][:, None]
    return fit_multi_delay(
        obs[i, j, k], m0[i, j, k], voxel_plds, metadata, labeling_efficiency(metadata), rows
    )
