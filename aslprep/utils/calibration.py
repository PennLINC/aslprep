"""Resolve the calibration and physiological parameters used to quantify perfusion.

CBF depends on how the calibration (M0) image was obtained and corrected and on physiological
parameters (labeling efficiency, blood T1, the blood-brain partition coefficient, ...).
Each of these may come from the command line, the BIDS metadata, an inference from the
acquisition, or a literature default.
:func:`resolve_calibration` resolves all of them once, records where each value came from,
and checks the inputs for missing, conflicting or implausible values.
Every quantification backend (ASLPrep's kinetic models, BASIL and partial volume correction)
reads its parameters from the resolved specification, so they all use the same values.

The specification is a JSON-serializable dictionary::

    {
        "Description": "...",
        "SourcePrecedence": ["user", "metadata", "inferred", "default"],
        "Parameters": {
            "<name>": {
                "Value": ...,
                "Units": "...",          # None for unitless values
                "Source": "user" | "metadata" | "inferred" | "default",
                "Origin": "...",         # what the value was taken or calculated from
                "Description": "...",
            },
            ...
        },
        "Warnings": ["...", ...],
    }

For each parameter, the first available source in ``SourcePrecedence`` wins:
a command-line value, then the BIDS metadata (including a haematocrit from a BIDS TSV file),
then a documented inference from other acquisition parameters, then a documented default.
"""

from __future__ import annotations

from numbers import Number

import numpy as np

#: The order in which the sources of a parameter's value take precedence.
SOURCE_PRECEDENCE = ('user', 'metadata', 'inferred', 'default')

#: Default of ``--m0-scale``.
DEFAULT_M0_SCALE = 1.0

#: Default of ``--smooth-kernel``, the FWHM in mm of the kernel the M0 image is smoothed with.
DEFAULT_SMOOTH_KERNEL = 5.0

#: M0 images acquired with a shorter repetition time than this, in seconds, are corrected for
#: incomplete T1 recovery (Alsop et al., 2015). oxford_asl uses the same threshold.
M0_RECOVERY_TR_THRESHOLD = 5.0

#: Factor converting CBF from mL/g/s to mL/100 g/min.
CBF_UNIT_CONVERSION = 6000

#: T1 of tissue, in seconds, that BASIL's kinetic model uses (the oxford_asl default).
#: ASLPrep's own models describe the decay of the label with the T1 of blood alone.
BASIL_TISSUE_T1 = 1.3

#: ASLPrep's kinetic models (``--kinetic-model``): the consensus models, in which the label
#: decays with the T1 of blood throughout, and the multi-delay model in which it decays with the
#: apparent T1 of tissue once it reaches the tissue.
KINETIC_MODELS = ('consensus', 'tissue-t1')

#: Calibration parameters that only the tissue-T1 kinetic model uses.
TISSUE_T1_PARAMETERS = (
    'GrayMatterT1',
    'WhiteMatterT1',
    'NominalGrayMatterCBF',
    'NominalWhiteMatterCBF',
)

#: Haematocrit values outside this range (as fractions) are flagged as implausible.
PLAUSIBLE_HAEMATOCRIT = (0.2, 0.65)

#: Labeling efficiencies below this value are flagged as implausible.
MIN_PLAUSIBLE_LABELING_EFFICIENCY = 0.5

_M0_SOURCES = {
    'Separate': 'the separate M0 scan',
    'Included': 'the M0 volumes in the ASL series',
    'Estimate': 'the M0Estimate value in the ASL sidecar',
    'Absent': 'the mean of the control volumes, as no M0 was available',
}


def _parameter(value, source, origin, description, units=None):
    """Build one entry of the specification's ``Parameters``."""
    if source not in SOURCE_PRECEDENCE:
        raise ValueError(f'Unknown source {source!r}; expected one of {SOURCE_PRECEDENCE}.')
    if isinstance(value, np.generic):
        value = value.item()
    return {
        'Value': value,
        'Units': units,
        'Source': source,
        'Origin': origin,
        'Description': description,
    }


def _positive(name, value):
    """Raise an error unless a value is a finite, positive number."""
    if (
        isinstance(value, bool)
        or not isinstance(value, Number)
        or not np.isfinite(value)
        or value <= 0
    ):
        raise ValueError(f'{name} must be a finite, positive number, but it is {value}.')


def resolve_calibration(
    metadata,
    *,
    m0tr=None,
    m0scan_metadata=None,
    m0_scale=DEFAULT_M0_SCALE,
    m0_scale_source='default',
    smooth_kernel=DEFAULT_SMOOTH_KERNEL,
    smooth_kernel_source='default',
    haematocrit_source=None,
    basil=False,
    kinetic_model='consensus',
):
    """Resolve the calibration and physiological parameters for an ASL run.

    Parameters
    ----------
    metadata : :obj:`dict`
        BIDS metadata of the ASL file. ``Haematocrit`` may have been added from a BIDS TSV
        file (see :func:`~aslprep.utils.bids.collect_haematocrit`).
    m0tr : :obj:`float` or None, optional
        Repetition time of the M0 volumes, in seconds
        (see :func:`~aslprep.utils.asl.infer_m0tr`), or None if there are none.
    m0scan_metadata : :obj:`dict` or None, optional
        Metadata of the separate M0 scan, if there is one.
    m0_scale : :obj:`float`, optional
        Factor the M0 image is multiplied by (``--m0-scale``).
    m0_scale_source, smooth_kernel_source : {"user", "default"}, optional
        Whether ``m0_scale`` and ``smooth_kernel`` were set on the command line.
    smooth_kernel : :obj:`float`, optional
        FWHM, in mm, of the Gaussian kernel the M0 image is smoothed with (``--smooth-kernel``).
    haematocrit_source : :obj:`str` or None, optional
        Where ``metadata['Haematocrit']`` came from, if not the ASL sidecar.
    basil : :obj:`bool`, optional
        Whether BASIL also quantifies CBF. Adds the parameters that only BASIL uses.
    kinetic_model : {"consensus", "tissue-t1"}, optional
        ASLPrep's kinetic model (``--kinetic-model``). ``tissue-t1`` adds the gray and white
        matter T1s and nominal CBFs it uses, and is only available for multi-delay data.

    Returns
    -------
    :obj:`dict`
        The calibration specification (see the module docstring).

    Raises
    ------
    ValueError
        If a value cannot be used (for example, a labeling efficiency above one).
    """
    from aslprep.utils.asl import (
        BACKGROUND_SUPPRESSION_PULSE_EFFICIENCY,
        TYPICAL_LABELING_EFFICIENCY,
        background_suppression_pulses,
        determine_multi_pld,
    )
    from aslprep.utils.cbf import (
        M0_THRESHOLD,
        NOMINAL_CBF,
        PARTITION_COEFFICIENT,
        REFERENCE_HAEMATOCRIT,
        estimate_t1blood,
        estimate_t1tissue,
        estimate_t1wm,
    )

    if kinetic_model not in KINETIC_MODELS:
        raise ValueError(
            f'Unknown kinetic model {kinetic_model!r}; expected one of {KINETIC_MODELS}.'
        )

    if kinetic_model == 'tissue-t1' and not determine_multi_pld(metadata=metadata):
        raise ValueError(
            'The tissue-T1 kinetic model (--kinetic-model tissue-t1) is only available for '
            'multi-delay data. Single-delay data need an arterial transit time to use it.'
        )

    params = {}
    warnings = []
    asl_type = metadata['ArterialSpinLabelingType']
    # Interfaces may be run on their own, with metadata that does not say where M0 came from.
    m0type = metadata.get('M0Type')
    field_strength = metadata['MagneticFieldStrength']
    _positive('MagneticFieldStrength', field_strength)
    params['MagneticFieldStrength'] = _parameter(
        field_strength,
        'metadata',
        'MagneticFieldStrength in the ASL sidecar',
        'Magnetic field strength, which sets the default relaxation times.',
        units='T',
    )
    if m0scan_metadata and 'MagneticFieldStrength' in m0scan_metadata:
        m0_field_strength = m0scan_metadata['MagneticFieldStrength']
        if not np.isclose(m0_field_strength, field_strength):
            warnings.append(
                f'The M0 scan was acquired at {m0_field_strength} T, but the ASL series at '
                f'{field_strength} T. The ASL field strength is used.'
            )

    # Calibration image
    if m0type is not None:
        m0_source = _M0_SOURCES[m0type]
        if m0type == 'Absent' and m0tr is None:
            # Without control volumes (only delta-M and CBF volumes), ExtractCBF uses an M0 of 1.
            m0_source = (
                'none, so a uniform M0 of 1 was used, as no M0 or control volumes were available'
            )
        params['M0Type'] = _parameter(
            m0type,
            'metadata',
            'M0Type in the ASL sidecar',
            f'Source of the calibration (M0) image: {m0_source}.',
        )
    if m0type == 'Estimate':
        _positive('M0Estimate', metadata.get('M0Estimate'))
        params['M0Estimate'] = _parameter(
            metadata['M0Estimate'],
            'metadata',
            'M0Estimate in the ASL sidecar',
            (
                'Single M0 value of arterial blood, as BIDS defines M0Estimate, used as a '
                'uniform calibration image. It is not divided by the blood-brain partition '
                'coefficient, which converts the M0 of tissue to that of blood, so the '
                'specification has no PartitionCoefficient.'
            ),
            units='arbitrary',
        )
    else:
        params['M0SmoothingFWHM'] = _parameter(
            smooth_kernel,
            smooth_kernel_source,
            '--smooth-kernel' if smooth_kernel_source == 'user' else 'ASLPrep default',
            'FWHM of the Gaussian kernel the calibration image was smoothed with.',
            units='mm',
        )

    _positive('--m0-scale', m0_scale)
    params['M0Scale'] = _parameter(
        m0_scale,
        m0_scale_source,
        '--m0-scale' if m0_scale_source == 'user' else 'ASLPrep default',
        (
            'Factor the calibration image was multiplied by, for the relative scaling of the '
            'ASL and M0 data. BIDS expects data to be scaled already, so this is normally 1.'
        ),
    )

    if m0tr is not None:
        # Rounded to the precision BASIL is given it, so both backends apply the same threshold.
        m0tr = round(float(np.mean(m0tr)), 6)
        _positive('RepetitionTimePreparation of the M0 volumes', m0tr)
        if m0type == 'Separate':
            origin = 'RepetitionTimePreparation in the M0 scan sidecar'
        elif m0type == 'Included':
            origin = 'RepetitionTimePreparation of the M0 volumes in the ASL sidecar'
        elif m0type == 'Absent':
            origin = 'RepetitionTimePreparation of the control volumes in the ASL sidecar'
        else:
            origin = 'RepetitionTimePreparation of the M0 volumes'
        params['M0RepetitionTime'] = _parameter(
            m0tr,
            'metadata',
            origin,
            'Repetition time of the calibration image.',
            units='s',
        )

    recovery = m0tr is not None and m0tr < M0_RECOVERY_TR_THRESHOLD
    params['M0RecoveryCorrection'] = _parameter(
        'saturation recovery' if recovery else 'none',
        'inferred',
        (
            f'M0 repetition time below {M0_RECOVERY_TR_THRESHOLD:g} s'
            if recovery
            else (
                'No M0 repetition time'
                if m0tr is None
                else f'M0 repetition time of at least {M0_RECOVERY_TR_THRESHOLD:g} s'
            )
        ),
        (
            'Correction of the calibration image for incomplete T1 recovery: '
            'division by 1 - exp(-TR / T1 of tissue) (Alsop et al., 2015).'
        ),
    )
    if recovery:
        t1tissue = estimate_t1tissue(field_strength)
        params['TissueT1'] = _parameter(
            t1tissue,
            'default',
            f'Gray matter T1 at {field_strength:g} T (Wright et al., 2008)',
            'T1 of tissue used to correct the calibration image for incomplete T1 recovery.',
            units='s',
        )
        params['M0RecoveryFactor'] = _parameter(
            1 / (1 - np.exp(-m0tr / t1tissue)),
            'inferred',
            '1 / (1 - exp(-M0RepetitionTime / TissueT1))',
            'Factor the calibration image was multiplied by to correct for T1 recovery.',
        )

    if m0type != 'Estimate':
        params['M0Threshold'] = _parameter(
            M0_THRESHOLD,
            'default',
            'ASLPrep default',
            (
                'Voxels whose calibration value is below this fraction of the 98th percentile '
                'of the calibration image get no CBF estimate.'
            ),
            units='fraction',
        )

    # Labeling efficiency
    if 'LabelingEfficiency' in metadata:
        labeleff = metadata['LabelingEfficiency']
        if isinstance(labeleff, bool) or not isinstance(labeleff, Number) or not 0 < labeleff <= 1:
            raise ValueError(
                f'LabelingEfficiency must be a number between zero and one, but it is {labeleff}.'
            )
        if labeleff < MIN_PLAUSIBLE_LABELING_EFFICIENCY:
            warnings.append(
                f'The LabelingEfficiency of {labeleff:g} is unusually low for {asl_type} data.'
            )
        params['LabelingEfficiency'] = _parameter(
            labeleff,
            'metadata',
            'LabelingEfficiency in the ASL sidecar',
            'Efficiency of the labeling itself, before background suppression.',
            units='fraction',
        )
    else:
        citation = 'Wang et al. (2005)' if asl_type == 'CASL' else 'Alsop et al. (2015)'
        params['LabelingEfficiency'] = _parameter(
            TYPICAL_LABELING_EFFICIENCY[asl_type],
            'default',
            f'Typical labeling efficiency for {asl_type} ({citation})',
            'Efficiency of the labeling itself, before background suppression.',
            units='fraction',
        )

    # Background suppression
    n_pulses, pulse_source = background_suppression_pulses(metadata)
    background_suppression = bool(metadata.get('BackgroundSuppression', False))
    declared_pulses = metadata.get('BackgroundSuppressionNumberPulses')
    pulse_times = metadata.get('BackgroundSuppressionPulseTime')
    if not background_suppression and declared_pulses:
        warnings.append(
            'BackgroundSuppression is false or not defined, but '
            f'BackgroundSuppressionNumberPulses is {declared_pulses}. '
            'No loss of label from background suppression is modeled.'
        )
    if background_suppression and pulse_source == 'metadata':
        if n_pulses == 0:
            warnings.append(
                'BackgroundSuppression is true, but BackgroundSuppressionNumberPulses is 0. '
                'No loss of label from background suppression is modeled.'
            )
        if pulse_times is not None and np.size(pulse_times) != n_pulses:
            warnings.append(
                f'BackgroundSuppressionNumberPulses ({n_pulses}) does not match the number of '
                f'BackgroundSuppressionPulseTime values ({np.size(pulse_times)}). '
                'BackgroundSuppressionNumberPulses is used.'
            )
    if pulse_source == 'default':
        warnings.append(
            'BackgroundSuppression is true, but neither BackgroundSuppressionNumberPulses nor '
            'BackgroundSuppressionPulseTime is defined. One background suppression pulse is '
            'assumed to estimate the labeling efficiency.'
        )
    pulse_origin = {
        'metadata': (
            'BackgroundSuppressionNumberPulses in the ASL sidecar'
            if background_suppression
            else 'BackgroundSuppression in the ASL sidecar'
        ),
        'inferred': 'Number of BackgroundSuppressionPulseTime values in the ASL sidecar',
        'default': 'Assumed, as the number of pulses was not specified',
    }[pulse_source]
    params['BackgroundSuppressionNumberPulses'] = _parameter(
        n_pulses,
        pulse_source,
        pulse_origin,
        'Number of background suppression pulses that reduced the label.',
    )
    if n_pulses:
        params['BackgroundSuppressionPulseEfficiency'] = _parameter(
            BACKGROUND_SUPPRESSION_PULSE_EFFICIENCY,
            'default',
            'Alsop et al. (2015)',
            'Fraction of the label left by each background suppression pulse.',
            units='fraction',
        )
    bs_efficiency = BACKGROUND_SUPPRESSION_PULSE_EFFICIENCY**n_pulses
    params['BackgroundSuppressionEfficiency'] = _parameter(
        bs_efficiency,
        'inferred',
        'BackgroundSuppressionPulseEfficiency ** BackgroundSuppressionNumberPulses',
        'Fraction of the label left after all background suppression pulses.',
        units='fraction',
    )
    params['EffectiveLabelingEfficiency'] = _parameter(
        params['LabelingEfficiency']['Value'] * bs_efficiency,
        'inferred',
        'LabelingEfficiency * BackgroundSuppressionEfficiency',
        'Labeling efficiency used to quantify CBF (alpha in Alsop et al., 2015).',
        units='fraction',
    )

    # Blood T1
    t1blood, t1blood_method = estimate_t1blood(metadata, warn=False)
    if t1blood_method == 'field strength formula':
        warnings.append(
            f'There is no default blood T1 for {field_strength} T data, so it was calculated '
            'from the field strength with the formula of Zhang et al. (2013).'
        )
    if 'Haematocrit' in metadata:
        haematocrit = metadata['Haematocrit']
        low, high = PLAUSIBLE_HAEMATOCRIT
        if not low <= haematocrit <= high:
            warnings.append(
                f'The haematocrit of {haematocrit:g} is outside the plausible range of '
                f'{low:g}-{high:g}.'
            )
        params['Haematocrit'] = _parameter(
            haematocrit,
            'metadata',
            haematocrit_source or 'Haematocrit in the ASL sidecar',
            'Measured haematocrit of the participant.',
            units='fraction',
        )
        if t1blood_method != 'haematocrit':
            warnings.append(
                f'The haematocrit of {haematocrit:g} was not used, because the model that '
                f'converts it to blood T1 is not valid at {field_strength} T.'
            )
    blood_t1_origin = {
        'field strength': (
            f'Blood T1 at {field_strength:g} T '
            + ('(Alsop et al., 2015)' if field_strength in (1.5, 3) else '(Zhang et al., 2013)')
        ),
        'field strength formula': (
            f'(110 * {field_strength:g} + 1316) / 1000, from the field strength '
            '(Zhang et al., 2013)'
        ),
        'haematocrit': (
            'Default blood T1 for the field strength, multiplied by the ratio of the blood T1 '
            'at the measured haematocrit to that at a haematocrit of '
            f'{REFERENCE_HAEMATOCRIT:g}, with the model of Hales et al. (2016)'
        ),
    }[t1blood_method]
    params['BloodT1'] = _parameter(
        t1blood,
        'inferred' if t1blood_method == 'haematocrit' else 'default',
        blood_t1_origin,
        'T1 of arterial blood.',
        units='s',
    )

    if m0type != 'Estimate':
        # An M0Estimate is already the M0 of blood (see m0_partition_coefficient).
        params['PartitionCoefficient'] = _parameter(
            PARTITION_COEFFICIENT,
            'default',
            'Alsop et al. (2015)',
            (
                'Blood-brain partition coefficient (lambda). The calibration image divided by '
                'it is the equilibrium magnetization of arterial blood.'
            ),
            units='mL/g',
        )
    params['CBFUnitConversion'] = _parameter(
        CBF_UNIT_CONVERSION,
        'default',
        'Unit conversion',
        'Factor converting CBF from mL/g/s to mL/100 g/min, the units of the CBF outputs.',
    )

    if kinetic_model == 'tissue-t1':
        for tissue, name, t1 in (
            ('gray', 'GrayMatterT1', estimate_t1tissue(field_strength)),
            ('white', 'WhiteMatterT1', estimate_t1wm(field_strength)),
        ):
            params[name] = _parameter(
                t1,
                'default',
                f'{tissue.capitalize()} matter T1 at {field_strength:g} T (Wright et al., 2008)',
                (
                    f'T1 of {tissue} matter in the tissue-T1 kinetic model, '
                    'in which the label decays with it once it reaches the tissue.'
                ),
                units='s',
            )
        for tissue, name, key in (
            ('gray', 'NominalGrayMatterCBF', 'gm'),
            ('white', 'NominalWhiteMatterCBF', 'wm'),
        ):
            params[name] = _parameter(
                NOMINAL_CBF[key],
                'default',
                'ASLPrep default',
                (
                    f'Nominal CBF of {tissue} matter in the tissue-T1 kinetic model. It sets the '
                    'flow term of the apparent tissue T1, 1 / (1 / T1 + CBF / lambda), and the '
                    "weight of the tissue's signal in voxels that mix gray and white matter."
                ),
                units='mL/100 g/min',
            )

    if basil:
        params['BASILTissueT1'] = _parameter(
            BASIL_TISSUE_T1,
            'default',
            'oxford_asl default',
            (
                "T1 of tissue in BASIL's kinetic model. "
                "ASLPrep's own models use the T1 of blood for the decay of the label."
            ),
            units='s',
        )

    return {
        'Description': (
            'Calibration and physiological parameters used to quantify perfusion, '
            'with the source of each value.'
        ),
        'SourcePrecedence': list(SOURCE_PRECEDENCE),
        'Parameters': params,
        'Warnings': warnings,
    }


def calibration_value(calibration, name, default=None):
    """Get a parameter's value from a calibration specification.

    Parameters
    ----------
    calibration : :obj:`dict`
        The specification, from :func:`resolve_calibration`.
    name : :obj:`str`
        The parameter's name.
    default : optional
        Returned if the specification does not include the parameter.
    """
    param = calibration['Parameters'].get(name)
    return default if param is None else param['Value']


def update_calibration(calibration, updates):
    """Replace parameters of a calibration specification with values found at run time.

    The specification is resolved when the workflow is built, but some values, such as tissue
    T1s measured from a T1 map, are only known when it runs. This replaces those parameters,
    so every output that records the specification records the values that were used.

    Parameters
    ----------
    calibration : :obj:`dict`
        The specification, from :func:`resolve_calibration`.
    updates : :obj:`dict`
        New entries of ``Parameters``, keyed by parameter name. Each is a complete entry, with the
        same keys as those :func:`resolve_calibration` writes.

    Returns
    -------
    :obj:`dict`
        A copy of the specification, with the updated parameters.

    Raises
    ------
    ValueError
        If a parameter is not in the specification, or an entry is incomplete or has an
        unknown source.
    """
    import copy

    keys = {'Value', 'Units', 'Source', 'Origin', 'Description'}
    calibration = copy.deepcopy(calibration)
    for name, entry in updates.items():
        # The workflow was built for the parameters the specification has, so an update can
        # replace them but not add new ones (which could, for example, change the kinetic model).
        if name not in calibration['Parameters']:
            raise ValueError(
                f'{name} is not a parameter of the calibration specification, so it cannot be '
                'updated at run time.'
            )

        if set(entry) != keys:
            raise ValueError(
                f'The update of {name} must have exactly the keys {sorted(keys)}, '
                f'but it has {sorted(entry)}.'
            )

        if entry['Source'] not in calibration['SourcePrecedence']:
            raise ValueError(
                f'Unknown source {entry["Source"]!r} for {name}; '
                f'expected one of {calibration["SourcePrecedence"]}.'
            )

        calibration['Parameters'][name] = _parameter(
            entry['Value'],
            entry['Source'],
            entry['Origin'],
            entry['Description'],
            units=entry['Units'],
        )
    return calibration


def calibration_kinetic_model(calibration):
    """Get the kinetic model a calibration specification was resolved for.

    Parameters
    ----------
    calibration : :obj:`dict`
        The specification, from :func:`resolve_calibration`.

    Returns
    -------
    {"consensus", "tissue-t1"}
        ``tissue-t1`` if the specification has the parameters of the tissue-T1 model.
    """
    if all(name in calibration['Parameters'] for name in TISSUE_T1_PARAMETERS):
        return 'tissue-t1'
    return 'consensus'


def tissue_t1_kinetics(calibration):
    """Get the apparent T1 and nominal CBF of gray and white matter for the tissue-T1 model.

    The apparent T1 of each tissue is calculated with its nominal CBF
    (see :func:`~aslprep.utils.cbf.apparent_tissue_t1`), and the blood-brain partition
    coefficient of the specification, or the default one for an ``M0Estimate``,
    whose specification has none.

    Parameters
    ----------
    calibration : :obj:`dict`
        The specification, from :func:`resolve_calibration`.

    Returns
    -------
    :obj:`dict` or None
        ``t1prime`` (gray, then white matter, in seconds) and ``nominal_cbf``
        (gray, then white matter, in mL/100 g/min), as arrays,
        or None if the specification is for the consensus model.
    """
    from aslprep.utils.cbf import PARTITION_COEFFICIENT, apparent_tissue_t1

    if calibration_kinetic_model(calibration) != 'tissue-t1':
        return None

    partition_coefficient = calibration_value(
        calibration,
        'PartitionCoefficient',
        PARTITION_COEFFICIENT,
    )
    t1tissue = np.array(
        [
            calibration_value(calibration, 'GrayMatterT1'),
            calibration_value(calibration, 'WhiteMatterT1'),
        ]
    )
    nominal_cbf = np.array(
        [
            calibration_value(calibration, 'NominalGrayMatterCBF'),
            calibration_value(calibration, 'NominalWhiteMatterCBF'),
        ]
    )
    return {
        't1prime': apparent_tissue_t1(t1tissue, nominal_cbf, partition_coefficient),
        'nominal_cbf': nominal_cbf,
    }


def m0_partition_coefficient(calibration):
    """Get the factor the calibration image is divided by to give the M0 of arterial blood.

    This is the blood-brain partition coefficient, which converts the M0 of tissue to that of
    blood (Alsop et al., 2015), except for an ``M0Estimate``: BIDS defines it as the M0 of
    blood already, so it is divided by 1.

    Parameters
    ----------
    calibration : :obj:`dict`
        The specification, from :func:`resolve_calibration`.
    """
    if calibration_value(calibration, 'M0Type') == 'Estimate':
        return 1.0
    return calibration_value(calibration, 'PartitionCoefficient')


#: The models that quantify ASLPrep's outputs, as named in :func:`describe_model`.
MODELS = ('aslprep', 'basil', 'lr')

#: Calibration parameters that a model does not use, and so leaves out of its metadata.
_UNUSED_PARAMETERS = {
    'aslprep': ('BASILTissueT1',),
    # oxford_asl does not drop voxels with a low M0, and has its own kinetic model.
    'basil': ('M0Threshold',) + TISSUE_T1_PARAMETERS,
    # Local regression unmixes the signal calibrated by ASLPrep's model, with its parameters.
    'lr': ('BASILTissueT1',),
}


def describe_model(calibration, model, *, asl_type, multi_delay):
    """Describe a model, and the calibration parameters it used, for the sidecars of its maps.

    The metadata follow the organization of model fits proposed for BIDS diffusion derivatives
    (BEP 016, bids-standard/bids-specification#2258): a ``Model`` dictionary,
    identical in the sidecars of all maps from that model, with a ``Description``,
    a ``TermURL`` and the ``Parameters`` of the fit.
    The parameters are the values of the calibration specification that the model used,
    under the same names; their provenance is in ``*_desc-calibration_parameters.json``.

    Parameters
    ----------
    calibration : :obj:`dict`
        The specification, from :func:`resolve_calibration`.
    model : {"aslprep", "basil", "lr"}
        ASLPrep's general kinetic model, BASIL, or local-regression partial volume correction.
    asl_type : :obj:`str`
        The ``ArterialSpinLabelingType``.
    multi_delay : :obj:`bool`
        Whether the data have several delays.

    Returns
    -------
    :obj:`dict`
        The ``Model`` dictionary.
    """
    if model not in MODELS:
        raise ValueError(f'Unknown model {model!r}; expected one of {MODELS}.')

    delays = 'multi-delay' if multi_delay else 'single-delay'
    if calibration_kinetic_model(calibration) == 'tissue-t1':
        kinetic_model = (
            'two-compartment general kinetic model, with tissue and macrovascular terms '
            '(Buxton et al., 1998; Woods et al., 2023), in which the label decays with the T1 of '
            'blood until it reaches the tissue and with the apparent T1 of gray and white matter '
            'afterwards, weighted by their partial volumes and nominal CBFs in mixed voxels'
        )
    elif multi_delay:
        kinetic_model = (
            'two-compartment general kinetic model, with tissue and macrovascular terms '
            '(Buxton et al., 1998; Woods et al., 2023)'
        )
    else:
        kinetic_model = (
            'single-compartment general kinetic model (Buxton et al., 1998; Alsop et al., 2015)'
        )

    if model == 'aslprep':
        description = f'The {kinetic_model} for {delays} {asl_type} data, ' + (
            'fit voxelwise by bounded least squares in ASLPrep.'
            if multi_delay
            else 'solved voxelwise in ASLPrep.'
        )
        term_url = (
            'https://doi.org/10.1002/mrm.1910400308'
            if multi_delay
            else 'https://doi.org/10.1002/mrm.25197'
        )
    elif model == 'basil':
        description = (
            f'Bayesian Inference for Arterial Spin Labeling (BASIL) of {delays} {asl_type} data, '
            "with spatial regularization (Chappell et al., 2009), run with FSL's oxford_asl. "
            'The desc-basilGM and desc-basilWM maps also use its partial volume correction '
            '(Chappell et al., 2011).'
        )
        term_url = 'https://doi.org/10.1109/TSP.2008.2005752'
    else:
        description = (
            'Local-regression partial volume correction (Asllani et al., 2008) of the CBF '
            f'from ASLPrep, which quantified it with the {kinetic_model}'
            + (
                '. Each tissue curve was then fit with the same kinetic model.'
                if multi_delay
                else '.'
            )
        )
        term_url = 'https://doi.org/10.1002/mrm.21670'

    unused = _UNUSED_PARAMETERS[model]
    parameters = {
        name: param['Value']
        for name, param in calibration['Parameters'].items()
        if name not in unused
    }
    if model == 'basil' and 'BASILTissueT1' not in parameters:
        parameters['BASILTissueT1'] = BASIL_TISSUE_T1

    return {'Description': description, 'TermURL': term_url, 'Parameters': parameters}


def basil_calibration_inputs(calibration):
    """Translate a calibration specification into inputs of the BASILCBF interface.

    oxford_asl calibrates voxelwise with a fixed partition coefficient of 0.9, and corrects
    the calibration image for incomplete T1 recovery whenever its ``--tr`` is below 5 s,
    assuming 3.2 s when ``--tr`` is not given.
    ASLPrep therefore always gives ``--tr``, so that BASIL corrects the same calibration
    images as ASLPrep, with the same tissue T1.

    An ``M0Estimate`` is the M0 of blood, so it is given as a single value (``--M0``),
    which oxford_asl does not divide by the partition coefficient.
    oxford_asl does not apply ``--cgain`` to it, so it is multiplied by ``M0Scale`` here.

    Parameters
    ----------
    calibration : :obj:`dict`
        The specification, from :func:`resolve_calibration` with ``basil=True``.

    Returns
    -------
    :obj:`dict`
        Keyword arguments for :class:`~aslprep.interfaces.cbf.BASILCBF`.

    Raises
    ------
    ValueError
        If the specification has a partition coefficient that oxford_asl cannot use.
    """
    inputs = {
        'alpha': calibration_value(calibration, 'EffectiveLabelingEfficiency'),
        't1blood': calibration_value(calibration, 'BloodT1'),
        't1tissue': calibration_value(calibration, 'BASILTissueT1', BASIL_TISSUE_T1),
    }
    m0_scale = calibration_value(calibration, 'M0Scale')
    if calibration_value(calibration, 'M0Type') == 'Estimate':
        inputs['m0_value'] = calibration_value(calibration, 'M0Estimate') * m0_scale
        return inputs

    partition_coefficient = calibration_value(calibration, 'PartitionCoefficient')
    if not np.isclose(partition_coefficient, 0.9):
        raise ValueError(
            'oxford_asl uses a fixed blood-brain partition coefficient of 0.9, but the '
            f'calibration uses {partition_coefficient}.'
        )

    inputs['m0_scale'] = m0_scale
    if calibration_value(calibration, 'M0RecoveryCorrection') == 'saturation recovery':
        inputs['m0tr'] = calibration_value(calibration, 'M0RepetitionTime')
        inputs['t1tissue_calib'] = calibration_value(calibration, 'TissueT1')
    else:
        # A TR at the threshold tells oxford_asl not to correct the calibration image.
        inputs['m0tr'] = max(
            calibration_value(calibration, 'M0RepetitionTime', M0_RECOVERY_TR_THRESHOLD),
            M0_RECOVERY_TR_THRESHOLD,
        )
    return inputs
