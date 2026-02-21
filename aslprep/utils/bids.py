# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Utilities to handle BIDS inputs."""

from __future__ import annotations

import json
import os
from collections import defaultdict
from functools import cache
from pathlib import Path

import yaml
from bids.layout import BIDSLayout
from bids.utils import listify
from nipype.interfaces.base import isdefined
from nipype.interfaces.utility.base import _ravel

from aslprep import config
from aslprep.data import load as load_data


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
            f'M0Type is {run_data["asl_metadata"]["M0Type"]}, '
            f'but an M0 scan was found at {run_data["m0scan"]}'
        )

    config.loggers.workflow.info(
        (
            f'Collected run data for {asl_file}:\n'
            f'{yaml.dump(run_data, default_flow_style=False, indent=4)}'
        ),
    )

    # Add metadata to dictionary now (we don't want to print these with the logger).
    run_data['asl_metadata'] = asl_metadata
    run_data['m0scan_metadata'] = m0scan_metadata

    return run_data


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
    derivatives_dir = Path(derivatives_dir)

    # search for both aslrefs
    for k, q in spec['baseline'].items():
        query = {**entities, **q}
        item = layout.get(return_type='filename', **query)
        if not item:
            continue
        derivs_cache[f'{k}_aslref'] = item[0] if len(item) == 1 else item

    transforms_cache = {}
    for xfm, q in spec['transforms'].items():
        # Transform extension will often not match provided entities
        #   (e.g., ".nii.gz" vs ".txt").
        # And transform suffixes will be "xfm",
        #   whereas relevant src file will be "bold".
        query = {**entities, **q}
        if xfm == 'aslref2fmap' and fieldmap_id:
            # fieldmaps have ids like auto_00000
            query['to'] = fieldmap_id.replace('_', '')
        item = layout.get(return_type='filename', **query)
        if not item:
            continue
        transforms_cache[xfm] = item[0] if len(item) == 1 else item
    derivs_cache['transforms'] = transforms_cache

    for k, q in spec['masks'].items():
        query = {**entities, **q}
        item = layout.get(return_type='filename', **query)
        if not item or len(item) != 1:
            continue

        derivs_cache[k] = item[0]

    return derivs_cache


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


def write_derivative_description(bids_dir, deriv_dir):
    """Write derivative dataset_description file."""
    from aslprep.__about__ import DOWNLOAD_URL, __url__, __version__

    bids_dir = Path(bids_dir)
    deriv_dir = Path(deriv_dir)
    desc = {
        'Name': 'ASLPrep - ASL PREProcessing workflow',
        'BIDSVersion': '1.9.0',
        'PipelineDescription': {
            'Name': 'ASLPrep',
            'Version': __version__,
            'CodeURL': DOWNLOAD_URL,
        },
        'CodeURL': __url__,
        'HowToAcknowledge': 'Please cite our paper '
        'and include the generated citation boilerplate within the Methods '
        'section of the text.',
    }

    # Keys that can only be set by environment
    if 'ASLPREP_DOCKER_TAG' in os.environ:
        desc['DockerHubContainerTag'] = os.environ['ASLPREP_DOCKER_TAG']

    if 'ASLPREP_SINGULARITY_URL' in os.environ:
        singularity_url = os.environ['ASLPREP_SINGULARITY_URL']
        desc['SingularityContainerURL'] = singularity_url

        singularity_md5 = _get_shub_version(singularity_url)
        if singularity_md5 and singularity_md5 is not NotImplemented:
            desc['SingularityContainerMD5'] = _get_shub_version(singularity_url)

    # Keys deriving from source dataset
    orig_desc = {}
    fname = bids_dir / 'dataset_description.json'
    if fname.exists():
        with fname.open() as fobj:
            orig_desc = json.load(fobj)

    if 'DatasetDOI' in orig_desc:
        desc['SourceDatasetsURLs'] = [f'https://doi.org/{orig_desc["DatasetDOI"]}']

    if 'License' in orig_desc:
        desc['License'] = orig_desc['License']

    with (deriv_dir / 'dataset_description.json').open('w') as fobj:
        json.dump(desc, fobj, indent=4)


def _get_shub_version(singularity_url):  # noqa: U100
    return NotImplemented


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
