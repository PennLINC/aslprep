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

from aslprep import config
from aslprep.data import load as load_data


@cache
def _get_layout(derivatives_dir: Path) -> BIDSLayout:
    return BIDSLayout(derivatives_dir, config=['bids', 'derivatives'], validate=False)


def collect_data(
    layout,
    participant_label,
    bids_filters=None,
):
    """Use pybids to retrieve the input data for a given participant."""
    queries = {
        'fmap': {'datatype': 'fmap'},
        'flair': {'datatype': 'anat', 'suffix': 'FLAIR'},
        't2w': {'datatype': 'anat', 'suffix': 'T2w'},
        't1w': {'datatype': 'anat', 'suffix': 'T1w'},
        'roi': {'datatype': 'anat', 'suffix': 'roi'},
        'sbref': {'datatype': 'perf', 'suffix': 'sbref'},
        'asl': {'datatype': 'perf', 'suffix': 'asl'},
    }

    bids_filters = bids_filters or {}
    for acq, entities in bids_filters.items():
        queries[acq].update(entities)

    subj_data = {
        dtype: sorted(
            layout.get(
                return_type='file',
                subject=participant_label,
                extension=['nii', 'nii.gz'],
                **query,
            )
        )
        for dtype, query in queries.items()
    }

    return subj_data


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
        if xfm == 'boldref2fmap' and fieldmap_id:
            # fieldmaps have ids like auto_00000
            query['to'] = fieldmap_id.replace('_', '')
        item = layout.get(return_type='filename', **query)
        if not item:
            continue
        transforms_cache[xfm] = item[0] if len(item) == 1 else item
    derivs_cache['transforms'] = transforms_cache
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
