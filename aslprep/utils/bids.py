# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Utilities to handle BIDS inputs."""

from __future__ import annotations

import json
import os
import re
from collections import defaultdict
from functools import cache
from pathlib import Path

import filelock
import yaml
from bids.layout import BIDSLayout
from bids.utils import listify
from nipype import logging
from nipype.interfaces.base import isdefined
from nipype.interfaces.utility.base import _ravel
from packaging.version import Version

from aslprep import config
from aslprep.data import load as load_data
from aslprep.utils.asl import IGNORABLE_VOLUME_TYPES, ignore_volume_types

LOGGER = logging.getLogger('nipype.utils')


@cache
def _get_layout(derivatives_dir: Path) -> BIDSLayout:
    return BIDSLayout(derivatives_dir, config=['bids', 'derivatives'], validate=False)


def collect_run_data(layout, asl_file):
    """Use pybids to retrieve the input data for a given participant."""
    queries = {
        'aslcontext': {'suffix': 'aslcontext', 'extension': '.tsv'},
        'sbref': {'suffix': 'sbref', 'extension': ['.nii', '.nii.gz']},
    }

    bids_file = layout.get_file(asl_file)

    run_data = {
        dtype: layout.get_nearest(
            bids_file.path,
            return_type='file',
            strict=False,  # aslcontext files aren't grabbed when strict=True, for some reason
            **query,
        )
        for dtype, query in queries.items()
    }

    if 'sbref' in config.workflow.ignore:
        config.loggers.workflow.info('Single-band reference files ignored.')
        run_data['sbref'] = None

    # The aslcontext file is required
    if not run_data['aslcontext']:
        raise FileNotFoundError(f'aslcontext file for {asl_file} not found.')

    # Hide the volume types the user asked to ignore (e.g., CBF maps calculated on the scanner)
    ignored_types = sorted(set(config.workflow.ignore or []) & set(IGNORABLE_VOLUME_TYPES))
    if ignored_types:
        original_aslcontext = run_data['aslcontext']
        run_data['aslcontext'] = ignore_volume_types(
            aslcontext=original_aslcontext,
            ignore=ignored_types,
            out_dir=Path(config.execution.work_dir) / 'aslcontext',
        )
        if run_data['aslcontext'] != original_aslcontext:
            config.loggers.workflow.info(
                f'Ignoring the {" and ".join(ignored_types)} volumes in {asl_file}, '
                'as requested with --ignore.'
            )

    # Now let's look for an m0scan
    m0scan_candidates = [
        f
        for f in bids_file.get_associations(kind='InformedBy')
        if f.entities['suffix'] == 'm0scan'
    ]
    if m0scan_candidates:
        if len(m0scan_candidates) > 1:
            config.loggers.workflow.warning(
                f'More than one M0 file found for {asl_file}. '
                f'Using the first one ({m0scan_candidates[0].path})'
            )
        run_data['m0scan'] = m0scan_candidates[0].path
    else:
        run_data['m0scan'] = None

    m0scan_metadata = None
    asl_metadata = layout.get_metadata(asl_file)
    if (asl_metadata['M0Type'] == 'Separate') and not run_data['m0scan']:
        raise FileNotFoundError(f'M0 file for {asl_file} not found.')
    elif asl_metadata['M0Type'] == 'Separate':
        m0scan_metadata = layout.get_file(run_data['m0scan']).get_metadata()
        if not m0scan_metadata:
            raise Exception(f'No metadata for m0scan: {run_data["m0scan"]}')
    elif run_data['m0scan']:
        raise ValueError(
            f'M0Type is {asl_metadata["M0Type"]}, but an M0 scan was found at {run_data["m0scan"]}'
        )

    config.loggers.workflow.info(
        (
            f'Collected run data for {asl_file}:\n'
            f'{yaml.dump(run_data, default_flow_style=False, indent=4)}'
        ),
    )

    # A measured haematocrit sets T1blood (see estimate_t1), so carry it with the ASL metadata.
    haematocrit, source = collect_haematocrit(layout, asl_file, asl_metadata)
    run_data['haematocrit_source'] = None
    if haematocrit is not None:
        config.loggers.workflow.info(
            f'Using a haematocrit of {haematocrit:g} for {asl_file}, from {source}.'
        )
        asl_metadata['Haematocrit'] = haematocrit
        # Recorded in the calibration specification (see aslprep.utils.calibration)
        run_data['haematocrit_source'] = source

    # Add metadata to dictionary now (we don't want to print these with the logger).
    run_data['asl_metadata'] = asl_metadata
    run_data['m0scan_metadata'] = m0scan_metadata

    return run_data


# Names accepted for a measured haematocrit, compared without regard to case.
# Haematocrit is the BIDS metadata field (defined for PET blood data).
HAEMATOCRIT_NAMES = ('haematocrit', 'hematocrit')


def _parse_haematocrit(value, source):
    """Convert a haematocrit value to a fraction, or return None if it cannot be used.

    BIDS defines haematocrit as a fraction of whole-blood volume.
    Values above 1 are read as percentages.
    """
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None

    if value != value:  # NaN, as pandas reads "n/a"
        return None

    if 1 < value <= 100:
        config.loggers.workflow.info(
            f'Reading the haematocrit of {value:g} in {source} as a percentage.'
        )
        value /= 100

    if not 0 < value < 1:
        config.loggers.workflow.warning(
            f'Ignoring the haematocrit of {value:g} in {source}: '
            'it must be a fraction between 0 and 1 (or a percentage).'
        )
        return None

    return value


def _haematocrit_from_tsv(tsv_file, key_column, key):
    """Return the haematocrit in the row of a BIDS TSV file whose ``key_column`` is ``key``."""
    import pandas as pd

    tsv_file = Path(tsv_file)
    if not tsv_file.is_file():
        return None

    df = pd.read_table(tsv_file, dtype=str, na_values='n/a')
    columns = [col for col in df.columns if col.lower() in HAEMATOCRIT_NAMES]
    if not columns or key_column not in df.columns:
        return None

    rows = df.loc[df[key_column] == key, columns[0]]
    if rows.empty:
        return None

    return _parse_haematocrit(rows.iloc[0], f'{tsv_file.name} ({key})')


def collect_haematocrit(layout, asl_file, metadata):
    """Find a measured haematocrit for an ASL file, from the most specific source available.

    The sources are checked in this order, and the first usable value wins:

    1. The ASL file's JSON sidecar, as a ``Haematocrit`` (or ``Hematocrit``) field.
    2. The ``scans.tsv`` file of the subject or session, in the ASL file's row.
    3. The subject's ``sessions.tsv`` file, in the ASL file's session's row.
    4. The dataset's ``participants.tsv`` file, in the subject's row.

    In the TSV files, the value is read from a ``haematocrit`` (or ``hematocrit``) column.

    Parameters
    ----------
    layout : :obj:`bids.layout.BIDSLayout`
        The layout of the input dataset.
    asl_file : :obj:`str`
        Path to the ASL file.
    metadata : :obj:`dict`
        The ASL file's metadata.

    Returns
    -------
    haematocrit : :obj:`float` or None
        The haematocrit, as a fraction of whole-blood volume.
    source : :obj:`str` or None
        A description of where the value was found.
    """
    for key, value in metadata.items():
        if key.lower() in HAEMATOCRIT_NAMES:
            haematocrit = _parse_haematocrit(value, f'the {key} field of the ASL sidecar')
            if haematocrit is not None:
                return haematocrit, f'the {key} field of the ASL sidecar'

    entities = layout.parse_file_entities(asl_file)
    subject = f'sub-{entities["subject"]}'
    session = f'ses-{entities["session"]}' if entities.get('session') else None
    root = Path(layout.root)
    data_dir = root / subject / session if session else root / subject

    scans_file = data_dir / '_'.join(filter(None, [subject, session, 'scans.tsv']))
    # scans.tsv lists files relative to the subject or session directory
    try:
        scans_key = Path(asl_file).resolve().relative_to(data_dir.resolve()).as_posix()
    except ValueError:
        scans_key = None
    candidates = [
        (scans_file, 'filename', scans_key),
        (root / subject / f'{subject}_sessions.tsv', 'session_id', session),
        (root / 'participants.tsv', 'participant_id', subject),
    ]
    for tsv_file, key_column, key in candidates:
        if key is None:
            continue
        haematocrit = _haematocrit_from_tsv(tsv_file, key_column, key)
        if haematocrit is not None:
            return haematocrit, tsv_file.name

    return None, None


def explain_missing_asl(layout, subject_id):
    """Explain why no ASL runs were found, when the dataset's location is the cause.

    PyBIDS versions without bids-standard/pybids#1273 parse the ``datatype`` entity from the
    whole path of each file.
    A dataset stored under a directory named after a BIDS datatype
    (for example, ``/data/projects/motion/``) then gets that datatype for every file,
    so a query for ``datatype='perf'`` finds nothing.

    Parameters
    ----------
    layout : :obj:`bids.layout.BIDSLayout`
        The layout of the input dataset.
    subject_id : :obj:`str`
        The participant label, without the ``sub-`` prefix.

    Returns
    -------
    :obj:`str`
        A hint to append to the error message.
        Empty if the participant has no ASL files with a datatype other than ``perf``.
    """
    asl_files = layout.get(subject=subject_id, suffix='asl', extension=['.nii', '.nii.gz'])
    datatypes = sorted(
        {f.entities.get('datatype') for f in asl_files} - {'perf', None},
    )
    if not datatypes:
        return ''

    root_parts = Path(layout.root).absolute().parts
    culprits = [datatype for datatype in datatypes if datatype in root_parts]
    hint = (
        f'ASLPrep found {len(asl_files)} ASL file(s) for this participant, '
        f'but PyBIDS assigned them the datatype {", ".join(repr(d) for d in datatypes)} '
        "instead of 'perf'. "
    )
    if culprits:
        hint += (
            f'The path to the dataset ({layout.root}) contains a directory named '
            f'{", ".join(repr(d) for d in culprits)}, which PyBIDS mistakes for the datatype '
            'of every file in the dataset (bids-standard/pybids#1270). '
            'Please move the dataset, or rename that directory, so that no directory above '
            'the dataset is named after a BIDS datatype.'
        )
    else:
        hint += 'Please check that the ASL files are in a perf/ directory.'

    return hint


def collect_derivatives(
    derivatives_dir: Path,
    entities: dict,
    fieldmap_id: str | None = None,
    spec: dict | None = None,
    patterns: list[str] | None = None,
):
    """Gather existing derivatives and compose a cache."""
    if spec is None or patterns is None:
        _spec, _patterns = tuple(
            json.loads(load_data.readable('io_spec.json').read_text()).values()
        )

        if spec is None:
            spec = _spec
        if patterns is None:
            patterns = _patterns

    derivs_cache = defaultdict(list, {})
    layout = _get_layout(derivatives_dir)

    # search for both aslrefs
    for k, q in spec['baseline'].items():
        item = _lookup(layout, f'baseline/{k}', q, entities)
        if item:
            derivs_cache[f'{k}_aslref'] = item

    transforms_cache = {}
    for xfm, q in spec['transforms'].items():
        # Transform extension will often not match provided entities
        #   (e.g., ".nii.gz" vs ".txt").
        # And transform suffixes will be "xfm",
        #   whereas relevant src file will be "bold".
        override = {}
        if xfm == 'aslref2fmap':
            if not fieldmap_id:
                # Without a fieldmap, the query would match any transform out of the reference.
                continue
            # fieldmaps have ids like auto_00000
            override['to'] = re.sub(r'[^a-zA-Z0-9]', '', fieldmap_id)
        item = _lookup(layout, f'transforms/{xfm}', q, entities, override)
        if item:
            transforms_cache[xfm] = item
    derivs_cache['transforms'] = transforms_cache

    for k, q in spec['masks'].items():
        # ``aslref_mask``, the key the fit workflow reads (it was stored as ``aslref``, so a
        # precomputed mask was never reused).
        item = _lookup(layout, f'masks/{k}', q, entities)
        if item:
            derivs_cache[f'{k}_mask'] = item

    return derivs_cache


def _lookup(layout, location, alternatives, entities, override=None):
    """Find the one file a query names, as nipost's collector does.

    ``alternatives`` is one entity dict, or a list of them describing the same product under
    different names (current naming first, then legacy naming). The first alternative that
    matches anything is used alone. More than one match is an error, so that a duplicated
    product is reported rather than passed on as a list or dropped.
    """
    if isinstance(alternatives, dict):
        alternatives = [alternatives]
    for i, alternative in enumerate(alternatives):
        query = {**entities, **alternative, **(override or {})}
        found = layout.get(return_type='filename', **query)
        if len(found) > 1:
            raise ValueError(
                f'{location} in {layout.root} (alternative {i}): expected at most one file, '
                f'found {len(found)}: {sorted(found)}. Remove all but one.'
            )
        if found:
            return found[0]
    return None


#: Products defined in the frame of the HMC reference. None needs another product's file, only
#: that frame: the motion and M0 transforms map into it, the fieldmap transform maps out of it,
#: the coregistration reference is that reference unwarped, the brain mask lies on the
#: coregistration reference's grid, and the ASL-to-anatomical transform maps out of it. The HMC
#: reference is built deterministically, so a recomputed one lands in the same frame unless the
#: settings that build it changed.
FRAME_PRODUCTS = (
    ('coreg_aslref',),
    ('aslref_mask',),
    ('transforms', 'hmc'),
    ('transforms', 'm0scan2aslref'),
    ('transforms', 'aslref2fmap'),
    ('transforms', 'aslref2anat'),
)

#: Products whose method ASLPrep records in their sidecar: product -> (sidecar field, option).
DERIVATIVE_METHODS = (
    (('transforms', 'hmc'), 'HMCMethod', '--hmc-method'),
    (('coreg_aslref',), 'SkullStripASL', '--skull-strip-asl'),
    (('aslref_mask',), 'SkullStripASL', '--skull-strip-asl'),
)


def _sidecar(path):
    """The JSON sidecar of a derivative (``.nii.gz``, ``.nii`` or ``.txt``)."""
    path = Path(path if isinstance(path, str | Path) else path[0])
    name = path.name
    for ext in ('.nii.gz', '.nii', '.txt'):
        if name.endswith(ext):
            return path.with_name(name[: -len(ext)] + '.json')
    return path.with_suffix('.json')


def merge_functional_derivatives(
    caches,
    hmc_method='typewise_offset',
    skull_strip_asl='synthstrip',
):
    """Merge the caches of several derivatives datasets into one set of products to reuse.

    Each dataset's ``collect_derivatives`` cache has a ``transforms`` entry, even when empty, so
    updating one dict with the next replaced every transform found before (in particular with
    the internal atlas datasets, which are collected last; PennLINC/aslprep#738). Here every
    product is merged on its own, a later dataset winning.

    Every product found is reused, so that a dataset can be seeded with single files, such as a
    hand-edited brain mask. Two things are logged or acted on:

    - A product whose sidecar records a method (:data:`DERIVATIVE_METHODS`: ``HMCMethod`` for the
      head-motion transforms, ``SkullStripASL`` for the coregistration reference and the brain
      mask) other than the requested one was made by ASLPrep under other settings, and is
      recomputed. A product without the record (made by hand, or before the record existed) is
      reused, and the log says its method is unknown.
    - Products in the frame of the HMC reference (:data:`FRAME_PRODUCTS`) reused while the HMC
      reference is recomputed, or taken from another derivatives dataset than the HMC reference,
      are reused, and the log says so: they must be in the frame of the reference the run uses.
      Their dependence is on that frame, not on each other's files, so they are not tied together
      (a hand-edited mask goes with a recomputed coregistration reference, and so on).

    Parameters
    ----------
    caches : list of dict
        ``collect_derivatives`` outputs, in the order of the derivatives datasets.
    hmc_method : str
        The requested head-motion correction method.
    skull_strip_asl : str
        The requested ASL brain extraction method.

    Returns
    -------
    precomputed : dict
        The products to reuse, in ``collect_derivatives``' layout.
    notes : list of str
        What was recomputed or reused with a caveat, and why, for the log.
    """
    found = {}
    for i, derivs_cache in enumerate(caches):
        for key, value in derivs_cache.items():
            if key == 'transforms':
                for name, path in (value or {}).items():
                    found['transforms', name] = (path, i)
            elif value:
                found[(key,)] = (value, i)

    notes = []
    requested = {'HMCMethod': hmc_method, 'SkullStripASL': skull_strip_asl}
    for key, field, option in DERIVATIVE_METHODS:
        if key not in found:
            continue
        sidecar = _sidecar(found[key][0])
        method = json.loads(sidecar.read_text()).get(field) if sidecar.exists() else None
        name = '/'.join(key)
        if method is None:
            notes.append(f'reusing {name}, whose {option} method is not recorded')
        elif method != requested[field]:
            del found[key]
            notes.append(
                f'recomputing {name}: it was made with {option} {method}, not {requested[field]}'
            )

    in_frame = [key for key in FRAME_PRODUCTS if key in found]
    if in_frame and ('hmc_aslref',) not in found:
        notes.append(
            f'reusing {", ".join("/".join(k) for k in in_frame)} with a recomputed HMC '
            'reference; they must be in its frame (as they are if they came from a run whose '
            'reference was built the same way)'
        )
    elif in_frame:
        source = found['hmc_aslref',][1]
        other = [key for key in in_frame if found[key][1] != source]
        if other:
            notes.append(
                f'reusing {", ".join("/".join(k) for k in other)} from another derivatives '
                'dataset than the HMC reference; they must be in its frame'
            )

    precomputed = {'transforms': {}}
    for key, (value, _) in found.items():
        if key[0] == 'transforms':
            precomputed['transforms'][key[1]] = value
        else:
            precomputed[key[0]] = value
    return precomputed, notes


def write_bidsignore(deriv_dir):
    """Write .bidsignore file."""
    bids_ignore = (
        '*.html',
        'logs/',
        'figures/',  # Reports
        '*_xfm.*',  # Unspecified transform files
        '*.surf.gii',  # Unspecified structural outputs
        # Unspecified functional outputs
        '*_aslref.nii.gz',
        '*_asl.func.gii',
        '*_mixing.tsv',
        '*_timeseries.tsv',
    )
    ignore_file = Path(deriv_dir) / '.bidsignore'

    ignore_file.write_text('\n'.join(bids_ignore) + '\n')


def write_derivative_description(bids_dir, deriv_dir, dataset_links=None):
    """Write derivative dataset_description file.

    Adapted from fMRIPrep, including the BIDS-URI DatasetLinks
    (nipreps/fmriprep#3255, nipreps/fmriprep#3267).
    """
    from aslprep import __version__

    DOWNLOAD_URL = f'https://github.com/PennLINC/aslprep/archive/{__version__}.tar.gz'

    bids_dir = Path(bids_dir)
    deriv_dir = Path(deriv_dir)
    desc = {
        'Name': 'ASLPrep - ASL PREProcessing workflow',
        'BIDSVersion': '1.9.0',
        'DatasetType': 'derivative',
        'GeneratedBy': [
            {
                'Name': 'ASLPrep',
                'Version': __version__,
                'CodeURL': DOWNLOAD_URL,
            },
        ],
        'HowToAcknowledge': (
            'Please cite our paper '
            'and include the generated citation boilerplate within the Methods '
            'section of the text.'
        ),
    }

    # Keys that can only be set by environment
    if 'ASLPREP_DOCKER_TAG' in os.environ:
        desc['GeneratedBy'][0]['Container'] = {
            'Type': 'docker',
            'Tag': f'pennlinc/aslprep:{os.environ["ASLPREP_DOCKER_TAG"]}',
        }
    if 'ASLPREP_SINGULARITY_URL' in os.environ:
        desc['GeneratedBy'][0]['Container'] = {
            'Type': 'singularity',
            'URI': os.getenv('ASLPREP_SINGULARITY_URL'),
        }

    # Keys deriving from source dataset
    orig_desc = {}
    fname = bids_dir / 'dataset_description.json'
    if fname.exists():
        orig_desc = json.loads(fname.read_text())

    if 'DatasetDOI' in orig_desc:
        desc['SourceDatasets'] = [
            {'URL': f'https://doi.org/{orig_desc["DatasetDOI"]}', 'DOI': orig_desc['DatasetDOI']}
        ]
    if 'License' in orig_desc:
        desc['License'] = orig_desc['License']

    # Declare the datasets referenced by BIDS URIs in Sources fields
    if dataset_links:
        desc['DatasetLinks'] = {k: str(v) for k, v in dataset_links.items()}
        if 'templateflow' in dataset_links:
            desc['DatasetLinks']['templateflow'] = 'https://github.com/templateflow/templateflow'

    (deriv_dir / 'dataset_description.json').write_text(json.dumps(desc, indent=4))


def find_atlas_entities(filename):
    """Extract atlas entities from filename."""
    import os

    fname = os.path.basename(filename)
    elements = fname.split('_')

    out = []
    for ent in ('tpl', 'atlas', 'res'):
        ent_parts = [el for el in elements if el.startswith(f'{ent}-')]
        ent_value = None
        if ent_parts:
            ent_value = ent_parts[0].split('-')[1]

        out.append(ent_value)

    suffix = elements[-1].split('.')[0]
    extension = '.' + '.'.join(elements[-1].split('.')[1:])
    out += [suffix, extension]

    return tuple(out)


def collect_anat_derivatives(
    derivatives_dir,
    subject_id,
    std_spaces,
    spec=None,
    patterns=None,
    session_id=None,
):
    """Gather existing derivatives and compose a cache."""
    from niworkflows.data import load as nwf_load

    if spec is None or patterns is None:
        _spec, _patterns = tuple(json.loads(load_data('smriprep.json').read_text()).values())

        if spec is None:
            spec = _spec
        if patterns is None:
            patterns = _patterns

    deriv_config = nwf_load('nipreps.json')
    layout = BIDSLayout(derivatives_dir, config=deriv_config, validate=False)

    derivs_cache = {}

    # Subject and session (if available) will be added to all queries
    qry_base = {'subject': subject_id}
    if session_id:
        qry_base['session'] = session_id

    for key, qry in spec['baseline'].items():
        qry = {**qry, **qry_base}
        item = layout.get(**qry)
        if not item:
            continue

        # Respect label order in queries
        if 'label' in qry:
            item = sorted(item, key=lambda x: qry['label'].index(x.entities['label']))

        paths = [item.path for item in item]

        if key.startswith('t2w_'):
            derivs_cache[key] = paths[0] if len(paths) == 1 else paths
        else:
            derivs_cache[f't1w_{key}'] = paths[0] if len(paths) == 1 else paths

    transforms = derivs_cache.setdefault('transforms', {})
    for _space in std_spaces:
        space = _space.replace(':cohort-', '+')
        for key, qry in spec['transforms'].items():
            qry = {**qry, **qry_base}
            qry['from'] = qry['from'] or space
            qry['to'] = qry['to'] or space
            item = layout.get(return_type='filename', **qry)
            if not item:
                continue
            transforms.setdefault(_space, {})[key] = item[0] if len(item) == 1 else item

    for key, qry in spec['surfaces'].items():
        qry = {**qry, **qry_base}
        item = layout.get(return_type='filename', **qry)
        if not item or len(item) != 2:
            continue

        derivs_cache[key] = sorted(item)

    for key, qry in spec['masks'].items():
        qry = {**qry, **qry_base}
        item = layout.get(return_type='filename', **qry)
        if not item or len(item) != 1:
            continue

        derivs_cache[key] = item[0]

    return derivs_cache


def _find_nearest_path(path_dict, input_path):
    """Find the nearest relative path from an input path to a dictionary of paths.

    If ``input_path`` is not relative to any of the paths in ``path_dict``,
    the absolute path string is returned.
    If ``input_path`` is already a BIDS-URI, then it will be returned unmodified.

    Parameters
    ----------
    path_dict : dict of (str, Path)
        A dictionary of paths.
    input_path : Path
        The input path to match.

    Returns
    -------
    matching_path : str
        The nearest relative path from the input path to a path in the dictionary.
        This is either the concatenation of the associated key from ``path_dict``
        and the relative path from the associated value from ``path_dict`` to ``input_path``,
        or the absolute path to ``input_path`` if no matching path is found from ``path_dict``.

    Examples
    --------
    >>> from pathlib import Path
    >>> path_dict = {
    ...     'bids::': Path('/data/derivatives/fmriprep'),
    ...     'bids:raw:': Path('/data'),
    ...     'bids:deriv-0:': Path('/data/derivatives/source-1'),
    ... }
    >>> input_path = Path('/data/derivatives/source-1/sub-01/func/sub-01_task-rest_bold.nii.gz')
    >>> _find_nearest_path(path_dict, input_path)  # match to 'bids:deriv-0:'
    'bids:deriv-0:sub-01/func/sub-01_task-rest_bold.nii.gz'
    >>> input_path = Path('/out/sub-01/func/sub-01_task-rest_bold.nii.gz')
    >>> _find_nearest_path(path_dict, input_path)  # no match- absolute path
    '/out/sub-01/func/sub-01_task-rest_bold.nii.gz'
    >>> input_path = Path('/data/sub-01/func/sub-01_task-rest_bold.nii.gz')
    >>> _find_nearest_path(path_dict, input_path)  # match to 'bids:raw:'
    'bids:raw:sub-01/func/sub-01_task-rest_bold.nii.gz'
    >>> input_path = 'bids::sub-01/func/sub-01_task-rest_bold.nii.gz'
    >>> _find_nearest_path(path_dict, input_path)  # already a BIDS-URI
    'bids::sub-01/func/sub-01_task-rest_bold.nii.gz'
    """
    # Don't modify BIDS-URIs
    if isinstance(input_path, str) and input_path.startswith('bids:'):
        return input_path

    input_path = Path(input_path)
    matching_path = None
    for key, path in path_dict.items():
        if input_path.is_relative_to(path):
            relative_path = input_path.relative_to(path)
            if (matching_path is None) or (len(relative_path.parts) < len(matching_path.parts)):
                matching_key = key
                matching_path = relative_path

    if matching_path is None:
        matching_path = str(input_path.absolute())
    else:
        matching_path = f'{matching_key}{matching_path}'

    return matching_path


def _get_bidsuris(in_files, dataset_links, out_dir):
    """Convert input paths to BIDS-URIs using a dictionary of dataset links."""
    in_files = listify(in_files)
    in_files = _ravel(in_files)
    # Remove undefined inputs
    in_files = [f for f in in_files if isdefined(f)]
    # Convert the dataset links to BIDS URI prefixes
    updated_keys = {f'bids:{k}:': Path(v) for k, v in dataset_links.items()}
    updated_keys['bids::'] = Path(out_dir)
    # Convert the paths to BIDS URIs
    out = [_find_nearest_path(updated_keys, f) for f in in_files]
    return out


def get_entity(filename, entity):
    """Extract a given entity from a BIDS filename via string manipulation.

    Parameters
    ----------
    filename : :obj:`str`
        Path to the BIDS file.
    entity : :obj:`str`
        The entity to extract from the filename.

    Returns
    -------
    entity_value : :obj:`str` or None
        The BOLD file's entity value associated with the requested entity.
    """
    import os
    import re

    folder, file_base = os.path.split(filename)

    # Allow + sign, which is not allowed in BIDS,
    # but is used by templateflow for the MNIInfant template.
    entity_values = re.findall(f'{entity}-([a-zA-Z0-9+]+)', file_base)
    entity_value = None if len(entity_values) < 1 else entity_values[0]
    if entity == 'space' and entity_value is None:
        foldername = os.path.basename(folder)
        if foldername == 'anat':
            entity_value = 'T1w'
        elif foldername == 'func':
            entity_value = 'native'
        else:
            raise ValueError(f'Unknown space for {filename}')

    return entity_value


def write_atlas_dataset_description(atlas_dir):
    """Write dataset_description.json file for Atlas derivatives.

    Parameters
    ----------
    atlas_dir : :obj:`str`
        Path to the output XCP-D Atlases dataset.
    """
    import json
    import os

    from aslprep import __version__

    DOWNLOAD_URL = f'https://github.com/PennLINC/aslprep/archive/{__version__}.tar.gz'

    desc = {
        'Name': 'ASLPrep Atlases',
        'DatasetType': 'derivative',
        'GeneratedBy': [
            {
                'Name': 'ASLPrep',
                'Version': __version__,
                'CodeURL': DOWNLOAD_URL,
            },
        ],
        'HowToAcknowledge': 'Include the generated boilerplate in the methods section.',
    }
    os.makedirs(atlas_dir, exist_ok=True)

    atlas_dset_description = os.path.join(atlas_dir, 'dataset_description.json')
    if os.path.isfile(atlas_dset_description):
        with open(atlas_dset_description) as fobj:
            old_desc = json.load(fobj)

        old_version = old_desc['GeneratedBy'][0]['Version']
        if Version(__version__).public != Version(old_version).public:
            LOGGER.warning(f'Previous output generated by version {old_version} found.')

    else:
        lock_atlas_dset_description = os.path.join(atlas_dir, 'dataset_description.json.lock')
        with filelock.SoftFileLock(lock_atlas_dset_description, timeout=60):
            with open(atlas_dset_description, 'w') as fobj:
                json.dump(desc, fobj, indent=4, sort_keys=True)
