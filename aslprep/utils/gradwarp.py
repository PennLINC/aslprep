"""Utilities for gradient nonlinearity correction.

Adapted from QSIPrep (https://github.com/PennLINC/qsiprep; BSD 3-Clause License,
Copyright (c) 2015-2018, the BBL developers team), which runs the same TORTOISE tool.
"""

import dataclasses
from pathlib import Path

#: Extensions TORTOISE's coefficient reader recognizes, plus ready-made ITK displacement fields.
COEFFICIENT_EXTENSIONS = ('.grad', '.dat', '.gc')
FIELD_EXTENSIONS = ('.nii', '.nii.gz')

#: ``--force`` values that set the warp dimensionality, overriding ImageType.
FORCED_WARP_DIMS = {'gradwarp1D': '1D', 'gradwarp3D': '3D'}


@dataclasses.dataclass(frozen=True)
class GradwarpPlan:
    """What gradient nonlinearity correction to apply to one ASL run.

    Attributes
    ----------
    gradient_file
        Coefficient file or ITK displacement field.
    warp_dim
        ``'3D'`` (full correction), ``'1D'`` (through-plane only), or ``None``
        (the scanner already corrected the images in 3D, so nothing is applied).
    is_ge
        Whether the data come from a GE scanner.
    basis
        ``'metadata'`` if ``warp_dim`` came from ImageType, ``'forced'`` if from ``--force``.
    slice_axis
        Voxel axis (``'i'``, ``'j'``, or ``'k'``) along which slices are stacked,
        which defines the through-plane direction for ``'1D'`` correction.
    """

    gradient_file: str
    warp_dim: str | None
    is_ge: bool
    basis: str
    slice_axis: str = 'k'


def is_displacement_field(gradient_file):
    """Check whether ``--gradient-file`` is a ready-made ITK field rather than coefficients.

    The standalone ``CreateNonlinearityDisplacementMap`` tool only expands coefficients.
    Handing it a NIfTI file yields an error or an all-zero field.
    """
    return str(gradient_file).endswith(FIELD_EXTENSIONS)


def image_type_tags(metadata):
    """Normalize ImageType, which may be a list or a backslash-joined string."""
    image_type = metadata.get('ImageType') or ()
    if isinstance(image_type, str):
        image_type = image_type.split('\\')
    return {str(tag).strip().upper() for tag in image_type}


def warp_dim_from_metadata(metadata):
    """Determine the remaining gradient distortion from a run's ImageType.

    ``DIS3D`` means the scanner corrected the image in 3D, so nothing remains.
    ``DIS2D`` means it corrected in-plane, so only the through-plane component remains.
    Otherwise, the image is uncorrected. DIS3D wins if both tags are present.
    """
    tags = image_type_tags(metadata)
    if 'DIS3D' in tags:
        return None
    if 'DIS2D' in tags:
        return '1D'
    return '3D'


def is_ge(metadata):
    """Check whether the Manufacturer field names GE."""
    return str(metadata.get('Manufacturer', '')).strip().upper().startswith('GE')


def forced_warp_dim(force):
    """Return the warp dimensionality set by ``--force``, or None."""
    forced = sorted({value for value in (force or []) if value in FORCED_WARP_DIMS})
    if len(forced) > 1:
        raise ValueError(
            f'"--force {forced[0]}" and "--force {forced[1]}" are mutually exclusive: '
            'a run is corrected in one dimension or in three, not both.'
        )
    return FORCED_WARP_DIMS[forced[0]] if forced else None


_GE_GUARD = (
    'Gradient nonlinearity correction from a coefficient file is not supported for GE data '
    '({}).\n\n'
    "TORTOISE's own pipeline shifts the displacement field's z origin after expanding GE "
    'coefficients, but the standalone CreateNonlinearityDisplacementMap tool used by ASLPrep '
    'does not, so the field would be misplaced.\n\n'
    'Either pass --gradient-file a ready-made ITK displacement field (.nii/.nii.gz), '
    'or pass --ignore gradwarp to skip gradient nonlinearity correction.'
)


def resolve_gradwarp_plan(metadata, asl_file, gradient_file=None, force=None, ignore=None):
    """Decide what gradient nonlinearity correction to apply to one ASL run.

    Parameters
    ----------
    metadata
        The ASL run's metadata.
    asl_file
        The ASL file, for messages.
    gradient_file, force, ignore
        Default to ``config.workflow.gradient_file``, ``.force``, and ``.ignore``.

    Returns
    -------
    GradwarpPlan or None
        None if no correction was requested.
    """
    from aslprep import config

    if gradient_file is None:
        gradient_file = config.workflow.gradient_file
    if force is None:
        force = config.workflow.force or []
    if ignore is None:
        ignore = config.workflow.ignore or []

    if not gradient_file or 'gradwarp' in ignore:
        return None

    warp_dim = forced_warp_dim(force)
    basis = 'forced'
    if warp_dim is None:
        warp_dim = warp_dim_from_metadata(metadata)
        basis = 'metadata'

    plan = GradwarpPlan(
        gradient_file=str(gradient_file),
        warp_dim=warp_dim,
        is_ge=is_ge(metadata),
        basis=basis,
        slice_axis=metadata.get('SliceEncodingDirection', 'k')[0],
    )
    if plan.is_ge and plan.warp_dim is not None and not is_displacement_field(gradient_file):
        raise ValueError(_GE_GUARD.format(Path(asl_file).name))

    return plan


def gradwarp_metadata(plan, jacobian):
    """Describe a run's gradient nonlinearity correction for derivative sidecars."""
    metadata = {
        'GradientNonlinearityCorrection': plan.warp_dim is not None,
        'GradientCoefficientFile': Path(plan.gradient_file).name,
    }
    if plan.warp_dim is not None:
        metadata['GradientWarpDimensions'] = plan.warp_dim
        metadata['GradientWarpJacobian'] = jacobian
    return metadata


def validate_gradient_flags(gradient_file, force, ignore):
    """Validate the ``--gradient-file``/``--force``/``--ignore`` combination.

    An unrecognized extension is rejected, rather than warned about as TORTOISE does,
    because TORTOISE then silently skips the correction.

    Raises
    ------
    ValueError
        If the flags are contradictory or the gradient file has an unrecognized extension.
    """
    forced = forced_warp_dim(force)
    flag = next((value for value in force if value in FORCED_WARP_DIMS), None)

    if forced and 'gradwarp' in ignore:
        raise ValueError(f'"--force {flag}" and "--ignore gradwarp" are contradictory.')

    if forced and not gradient_file:
        raise ValueError(f'"--force {flag}" requires --gradient-file.')

    if 'gradwarp-jacobian' in ignore and not gradient_file:
        raise ValueError('"--ignore gradwarp-jacobian" requires --gradient-file.')

    if gradient_file and not str(gradient_file).endswith(
        COEFFICIENT_EXTENSIONS + FIELD_EXTENSIONS
    ):
        raise ValueError(
            f'--gradient-file {gradient_file} has an unrecognized extension. Expected a '
            "Siemens '.grad', GE '.dat', or TORTOISE '.gc' coefficient file, or an ITK "
            "displacement field ('.nii' or '.nii.gz')."
        )


# Siemens .grad sanitizing
#
# TORTOISE's ``GRADCAL::read_Siemens_format`` parses every line holding a "(" at index 3-9
# together with a "," and a ")". It has no comment handling, so a header comment written in
# the same notation as the data, such as
#
#     #  A(1,1) = 1.1547 (2/Sqrt[3])
#
# is read as a coefficient. Here, "std::stof" throws and the tool aborts. A comment that does
# parse is worse: it is silently added to the expansion. Real Siemens files carry such comments.
# The functions below emulate the reader to drop such comments into a copy of the file, and to
# fail early on data lines the reader cannot parse.

_WHITESPACE = ' \t\n\r\f\v'


def _stof(text):
    """Emulate ``std::stof``: read a leading float, ignore trailing junk, raise if none."""
    stripped = text.lstrip(_WHITESPACE)
    index = 0
    if index < len(stripped) and stripped[index] in '+-':
        index += 1
    digits = 0
    while index < len(stripped) and (stripped[index].isdigit() or stripped[index] == '.'):
        digits += stripped[index].isdigit()
        index += 1
    if not digits:
        raise ValueError('stof')
    return float(stripped[:index])


def _stoi(text):
    """Emulate ``std::stoi``: read a leading integer, ignore trailing junk, raise if none."""
    stripped = text.lstrip(_WHITESPACE)
    index = 0
    if index < len(stripped) and stripped[index] in '+-':
        index += 1
    digits = 0
    while index < len(stripped) and stripped[index].isdigit():
        index += 1
        digits += 1
    if not digits:
        raise ValueError('stoi')
    return int(stripped[:index])


def _substr(text, pos, count):
    """Emulate ``std::string::substr``, where a negative count wraps to "until the end"."""
    if pos > len(text):
        raise IndexError('out_of_range')
    return text[pos:] if count < 0 else text[pos : pos + count]


def _is_comment(line):
    return line.lstrip().startswith('#')


def siemens_reader_verdict(line):
    """Predict what TORTOISE's Siemens reader does with one line.

    Returns
    -------
    status : {'skip', 'term', 'abort'}
        Whether the reader ignores the line, reads a coefficient from it, or throws.
    detail : str or None
        For 'abort', which call failed and on what text.
    """
    pos_open = line.find('(', 3)
    if pos_open == -1 or pos_open >= 10:
        return 'skip', None
    # Both are searched from the start of the line, not from the "(".
    pos_comma = line.find(',')
    pos_close = line.find(')')
    if pos_comma == -1 or pos_close == -1:
        return 'skip', None

    degree = _substr(line, pos_open + 1, pos_comma - pos_open - 1)
    order = _substr(line, pos_comma + 1, pos_close - pos_comma - 1)
    try:
        _stoi(degree)
        _stoi(order)
    except ValueError:
        return 'abort', f'std::stoi on {degree!r}/{order!r}'

    # Everything after ")" except the final character, which is taken as the axis letter.
    coefficient = _substr(line, pos_close + 1, len(line) - pos_close - 2)
    try:
        _stof(coefficient)
    except ValueError:
        return 'abort', f'std::stof on {coefficient!r}'

    return 'term', None


def _lines(path):
    r"""Split like ``std::getline(f, s, '\n')``, which keeps a trailing ``\r``."""
    lines = Path(path).read_bytes().decode('latin-1').split('\n')
    if lines and lines[-1] == '':
        lines.pop()
    return lines


def sanitize_siemens_coefficients(gradient_file, dest_dir, logger=None):
    """Return a coefficient file TORTOISE's Siemens reader can read correctly.

    Comment lines the reader would parse are dropped into a copy in ``dest_dir``,
    with the original file name. Other files, and files with nothing to drop,
    are returned unchanged. The user's file is never modified.

    Raises
    ------
    ValueError
        If a line that is not a comment would make the reader throw.
    """
    gradient_file = Path(gradient_file)
    if gradient_file.suffix != '.grad':
        return gradient_file

    lines = _lines(gradient_file)
    verdicts = [
        (number, line, _is_comment(line), *siemens_reader_verdict(line))
        for number, line in enumerate(lines, start=1)
    ]

    fatal = [v for v in verdicts if v[3] == 'abort' and not v[2]]
    if fatal:
        number, line, _, _, detail = fatal[0]
        raise ValueError(
            f'{gradient_file} cannot be read by TORTOISE: line {number} makes its Siemens '
            f'coefficient reader fail ({detail}). Offending line: {line!r}. '
            'ASLPrep drops comment lines automatically, but this line holds data, '
            'so the file itself must be corrected.'
        )

    dropped = [v for v in verdicts if v[2] and v[3] != 'skip']
    if not dropped:
        return gradient_file

    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    out_file = dest_dir / gradient_file.name
    kept = [line for line in lines if not _is_comment(line)]
    out_file.write_bytes('\n'.join(kept).encode('latin-1') + b'\n')

    if logger is not None:
        logger.warning(
            'Gradient coefficient file %s has %d comment line(s) that TORTOISE would read as '
            'coefficients. Using a copy without them: %s',
            gradient_file,
            len(dropped),
            out_file,
        )
        for number, line, _, status, _ in dropped:
            effect = 'would abort the tool' if status == 'abort' else 'would add a term'
            logger.warning('  dropped line %d (%s): %r', number, effect, line)

    return out_file
