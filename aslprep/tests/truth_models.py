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
#: ASLPrep's assumed efficiency per background-suppression pulse (a convention, not physics).
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
        return np.full(asl.shape[:3], float(metadata['M0Estimate'])), None
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
    """ASLPrep's convention: the sidecar value, else base x 0.95^n for n suppression pulses."""
    if 'LabelingEfficiency' in metadata:
        return float(metadata['LabelingEfficiency'])
    alpha = BASE_EFFICIENCY[metadata['ArterialSpinLabelingType']]
    if metadata.get('BackgroundSuppression'):
        alpha *= BS_PULSE_EFFICIENCY ** metadata.get('BackgroundSuppressionNumberPulses', 1)
    return alpha


def labeling_efficiency_physical(simulation):
    """The efficiency the simulator applied: labeling efficiency x suppression attenuation.

    ``simulation`` is the sidecar's ``AslscanSimulation`` block.
    """
    alpha = simulation['Resolved']['LabelingEfficiency']['Value']
    factor = simulation.get('BackgroundSuppressionLabelFactor')
    return alpha * abs(factor) if factor is not None else alpha


def pld_map(metadata, shape, reverse_slices=False, pld_shift=0.0):
    """Post-labeling delay (PCASL) or inversion time (PASL) per voxel, with slice timing.

    2D acquisitions add each slice's acquisition time (BIDS ``SliceTiming``, along the
    ``SliceEncodingDirection`` axis, reversed for ``-``). Mutations: ``reverse_slices`` and
    ``pld_shift``.
    """
    pld = float(np.mean(metadata['PostLabelingDelay'])) + pld_shift
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


def cbf_single_delay(deltam, m0, metadata, alpha, pld):
    """Single-delay CBF (ml/100 g/min) from delta-M and M0, per Alsop 2015.

    PCASL/CASL (eq. 1): ``6000 lambda dM e^(PLD/T1b) / (2 alpha T1b M0 (1 - e^(-tau/T1b)))``.
    PASL with QUIPSS II or Q2TIPS (eq. 2): ``6000 lambda dM e^(TI/T1b) / (2 alpha TI1 M0)``,
    where ``TI`` is the inversion time (BIDS ``PostLabelingDelay`` for PASL) and ``TI1`` the
    bolus duration (the first ``BolusCutOffDelayTime``).
    """
    t1b = T1_BLOOD[metadata['MagneticFieldStrength']]
    kind = metadata['ArterialSpinLabelingType']
    if kind in ('PCASL', 'CASL'):
        tau = float(np.mean(metadata['LabelingDuration']))
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
        reverse_slices=mutation == 'reverse_slices',
        pld_shift=0.1 if mutation == 'pld_shift' else 0.0,
    )
    return cbf_single_delay(deltam, m0, metadata, alpha, pld)
