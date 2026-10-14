"""Tests for aslprep.utils.bids."""

import json
from pathlib import Path

import pytest

from aslprep.utils.bids import collect_derivatives, explain_missing_asl

XFM_NAMES = {
    'desc': {
        'hmc': 'sub-01_from-orig_to-aslref_mode-image_desc-hmc_xfm.txt',
        'aslref2anat': 'sub-01_from-aslref_to-T1w_mode-image_desc-coreg_xfm.txt',
        'aslref2fmap': 'sub-01_from-aslref_to-auto00000_mode-image_desc-fmap_xfm.txt',
    },
    'legacy': {
        'hmc': 'sub-01_from-orig_to-aslref_mode-image_xfm.txt',
        'aslref2anat': 'sub-01_from-aslref_to-T1w_mode-image_xfm.txt',
        'aslref2fmap': 'sub-01_from-aslref_to-auto00000_mode-image_xfm.txt',
    },
}


@pytest.mark.parametrize('naming', ['desc', 'legacy'])
def test_collect_derivatives_transforms(tmp_path, naming):
    """Precomputed ASL transforms are found with or without desc entities.

    nipreps/fmriprep#3532 added desc entities (hmc, coreg, fmap) to these transforms,
    but derivatives from earlier versions must still be reusable.
    """
    deriv_dir = tmp_path / 'derivatives'
    perf_dir = deriv_dir / 'sub-01' / 'perf'
    perf_dir.mkdir(parents=True)
    (deriv_dir / 'dataset_description.json').write_text(
        json.dumps({'Name': 'test', 'BIDSVersion': '1.9.0', 'DatasetType': 'derivative'})
    )
    for fname in XFM_NAMES[naming].values():
        (perf_dir / fname).write_text('')

    derivs = collect_derivatives(
        derivatives_dir=deriv_dir,
        entities={'subject': '01'},
        fieldmap_id='auto_00000',
    )
    transforms = derivs['transforms']
    for key, fname in XFM_NAMES[naming].items():
        assert transforms[key] == str(perf_dir / fname)


def test_collect_derivatives_sanitizes_fieldmap_id(tmp_path):
    """Non-alphanumeric characters are removed from the fieldmap ID in transform queries.

    Reproduces nipreps/fmriprep#3490.
    """
    deriv_dir = tmp_path / 'derivatives'
    perf_dir = deriv_dir / 'sub-01' / 'perf'
    perf_dir.mkdir(parents=True)
    (deriv_dir / 'dataset_description.json').write_text(
        json.dumps({'Name': 'test', 'BIDSVersion': '1.9.0', 'DatasetType': 'derivative'})
    )
    xfm = perf_dir / 'sub-01_from-aslref_to-myfmap01_mode-image_desc-fmap_xfm.txt'
    xfm.write_text('')

    derivs = collect_derivatives(
        derivatives_dir=deriv_dir,
        entities={'subject': '01'},
        fieldmap_id='my-fmap_01',
    )
    assert derivs['transforms']['aslref2fmap'] == str(xfm)


def test_write_derivative_description(tmp_path, monkeypatch):
    """DatasetLinks declare the datasets referenced by BIDS URIs in Sources."""
    from aslprep.utils.bids import write_derivative_description

    bids_dir = tmp_path / 'bids'
    bids_dir.mkdir()
    (bids_dir / 'dataset_description.json').write_text(
        json.dumps({'Name': 'raw', 'BIDSVersion': '1.9.0', 'DatasetDOI': '10.1/abc'})
    )
    deriv_dir = tmp_path / 'derivatives'
    deriv_dir.mkdir()
    monkeypatch.setenv('ASLPREP_DOCKER_TAG', 'unstable')

    write_derivative_description(
        bids_dir,
        deriv_dir,
        dataset_links={'raw': bids_dir, 'templateflow': tmp_path / 'tf'},
    )
    desc = json.loads((deriv_dir / 'dataset_description.json').read_text())
    assert desc['DatasetLinks'] == {
        'raw': str(bids_dir),
        'templateflow': 'https://github.com/templateflow/templateflow',
    }
    assert desc['SourceDatasets'] == [{'URL': 'https://doi.org/10.1/abc', 'DOI': '10.1/abc'}]
    assert desc['GeneratedBy'][0]['Container'] == {
        'Type': 'docker',
        'Tag': 'pennlinc/aslprep:unstable',
    }


def test_collect_run_data_rejects_unexpected_m0scan(tmp_path, monkeypatch):
    """An m0scan for a run whose M0Type is not Separate is reported, not a KeyError."""
    import nibabel as nb
    import numpy as np
    from bids.layout import BIDSLayout

    from aslprep import config
    from aslprep.utils.bids import collect_run_data

    monkeypatch.setattr(config.workflow, 'ignore', [])
    (tmp_path / 'dataset_description.json').write_text(
        json.dumps({'Name': 'test', 'BIDSVersion': '1.10.0'})
    )
    perf = tmp_path / 'sub-01' / 'perf'
    perf.mkdir(parents=True)
    img = nb.Nifti1Image(np.zeros((2, 2, 2, 3), np.float32), np.eye(4))
    img.to_filename(perf / 'sub-01_asl.nii.gz')
    (perf / 'sub-01_asl.json').write_text(json.dumps({'M0Type': 'Included'}))
    (perf / 'sub-01_aslcontext.tsv').write_text('volume_type\nm0scan\ncontrol\nlabel\n')
    img.slicer[..., 0].to_filename(perf / 'sub-01_m0scan.nii.gz')
    (perf / 'sub-01_m0scan.json').write_text(
        json.dumps({'IntendedFor': ['bids::sub-01/perf/sub-01_asl.nii.gz']})
    )
    layout = BIDSLayout(tmp_path, validate=False)
    with pytest.raises(ValueError, match='M0Type is Included, but an M0 scan was found'):
        collect_run_data(layout, str(perf / 'sub-01_asl.nii.gz'))


@pytest.mark.parametrize(
    ('ignore', 'expected'),
    [
        ([], ['m0scan', 'control', 'label', 'cbf']),
        (['cbf'], ['m0scan', 'control', 'label', 'ignored']),
    ],
)
def test_collect_run_data_ignores_volume_types(tmp_path, monkeypatch, ignore, expected):
    """--ignore cbf hides CBF volumes from the workflow without touching the dataset."""
    import nibabel as nb
    import numpy as np
    import pandas as pd
    from bids.layout import BIDSLayout

    from aslprep import config
    from aslprep.utils.bids import collect_run_data

    bids_dir = tmp_path / 'bids'
    work_dir = tmp_path / 'work'
    monkeypatch.setattr(config.workflow, 'ignore', ignore)
    monkeypatch.setattr(config.execution, 'work_dir', work_dir)
    perf = bids_dir / 'sub-01' / 'perf'
    perf.mkdir(parents=True)
    (bids_dir / 'dataset_description.json').write_text(
        json.dumps({'Name': 'test', 'BIDSVersion': '1.10.0'})
    )
    img = nb.Nifti1Image(np.zeros((2, 2, 2, 4), np.float32), np.eye(4))
    img.to_filename(perf / 'sub-01_asl.nii.gz')
    (perf / 'sub-01_asl.json').write_text(json.dumps({'M0Type': 'Included'}))
    raw_aslcontext = perf / 'sub-01_aslcontext.tsv'
    raw_text = 'volume_type\nm0scan\ncontrol\nlabel\ncbf\n'
    raw_aslcontext.write_text(raw_text)
    layout = BIDSLayout(bids_dir, validate=False)

    run_data = collect_run_data(layout, str(perf / 'sub-01_asl.nii.gz'))

    assert pd.read_table(run_data['aslcontext'])['volume_type'].tolist() == expected
    assert raw_aslcontext.read_text() == raw_text
    if ignore:
        assert Path(run_data['aslcontext']).parent == work_dir / 'aslcontext'
    else:
        assert Path(run_data['aslcontext']) == raw_aslcontext


def _make_asl_dataset(root):
    """Write a minimal dataset with one ASL run and return its layout."""
    from bids.layout import BIDSLayout

    perf_dir = root / 'sub-01' / 'perf'
    perf_dir.mkdir(parents=True)
    (root / 'dataset_description.json').write_text(
        json.dumps({'Name': 'test', 'BIDSVersion': '1.8.0'})
    )
    (perf_dir / 'sub-01_asl.nii.gz').write_bytes(b'')
    (perf_dir / 'sub-01_asl.json').write_text('{}')
    return BIDSLayout(str(root), validate=False)


def test_explain_missing_asl_no_hint(tmp_path):
    """There is no hint when the ASL files have the right datatype, or do not exist."""
    layout = _make_asl_dataset(tmp_path / 'dset')

    assert explain_missing_asl(layout, '01') == ''
    assert explain_missing_asl(layout, '02') == ''


def test_explain_missing_asl_datatype_in_path(tmp_path):
    """A dataset under a directory named after a BIDS datatype gets a hint (gh-704)."""
    layout = _make_asl_dataset(tmp_path / 'motion' / 'dset')
    asl_files = layout.get(subject='01', suffix='asl', extension='.nii.gz')
    if asl_files[0].entities.get('datatype') == 'perf':
        pytest.skip('This version of PyBIDS only parses entities within the dataset root.')

    hint = explain_missing_asl(layout, '01')
    assert 'found 1 ASL file(s)' in hint
    assert "datatype 'motion' instead of 'perf'" in hint
    assert "contains a directory named 'motion'" in hint


def test_explain_missing_asl_wrong_directory(tmp_path):
    """ASL files outside of perf/ get a hint that does not blame the dataset path."""
    from types import SimpleNamespace

    asl_file = SimpleNamespace(entities={'datatype': 'func', 'suffix': 'asl'})
    layout = SimpleNamespace(root=str(tmp_path / 'dset'), get=lambda **kwargs: [asl_file])

    hint = explain_missing_asl(layout, '01')
    assert "datatype 'func' instead of 'perf'" in hint
    assert 'perf/ directory' in hint
    assert 'contains a directory named' not in hint


def _haematocrit_dataset(
    root, sessions=True, sidecar=None, scans=None, session_rows=None, participant_rows=None
):
    """Write a minimal BIDS dataset with one ASL run and optional haematocrit sources."""
    from bids.layout import BIDSLayout

    root.mkdir(parents=True)
    (root / 'dataset_description.json').write_text(
        json.dumps({'Name': 'test', 'BIDSVersion': '1.10.0'})
    )
    data_dir = root / 'sub-01' / 'ses-1' if sessions else root / 'sub-01'
    prefix = 'sub-01_ses-1' if sessions else 'sub-01'
    perf = data_dir / 'perf'
    perf.mkdir(parents=True)
    asl_file = perf / f'{prefix}_asl.nii.gz'
    asl_file.write_bytes(b'')
    (perf / f'{prefix}_asl.json').write_text(json.dumps({'M0Type': 'Absent', **(sidecar or {})}))
    (perf / f'{prefix}_aslcontext.tsv').write_text('volume_type\ncontrol\nlabel\n')
    if scans is not None:
        (data_dir / f'{prefix}_scans.tsv').write_text(
            f'filename\t{scans[0]}\nperf/{prefix}_asl.nii.gz\t{scans[1]}\n'
        )
    if session_rows is not None:
        (root / 'sub-01' / 'sub-01_sessions.tsv').write_text(session_rows)
    if participant_rows is not None:
        (root / 'participants.tsv').write_text(participant_rows)
    layout = BIDSLayout(root, validate=False)
    return layout, str(asl_file)


@pytest.mark.parametrize(
    ('sources', 'expected', 'source'),
    [
        # Each source on its own
        (
            {'participant_rows': 'participant_id\thaematocrit\nsub-01\t0.41\n'},
            0.41,
            'participants',
        ),
        ({'session_rows': 'session_id\thaematocrit\nses-1\t0.39\n'}, 0.39, 'sessions'),
        ({'scans': ('haematocrit', '0.37')}, 0.37, 'scans'),
        ({'sidecar': {'Haematocrit': 0.35}}, 0.35, 'Haematocrit field'),
        # The most specific source wins
        (
            {
                'sidecar': {'Haematocrit': 0.35},
                'scans': ('haematocrit', '0.37'),
                'session_rows': 'session_id\thaematocrit\nses-1\t0.39\n',
                'participant_rows': 'participant_id\thaematocrit\nsub-01\t0.41\n',
            },
            0.35,
            'Haematocrit field',
        ),
        (
            {
                'scans': ('haematocrit', '0.37'),
                'participant_rows': 'participant_id\thaematocrit\nsub-01\t0.41\n',
            },
            0.37,
            'scans',
        ),
        # A missing value falls through to the next source
        (
            {
                'session_rows': 'session_id\thaematocrit\nses-1\tn/a\n',
                'participant_rows': 'participant_id\thaematocrit\nsub-01\t0.41\n',
            },
            0.41,
            'participants',
        ),
        # Percentages, the American spelling, and other participants' rows
        ({'participant_rows': 'participant_id\tHematocrit\nsub-01\t43\n'}, 0.43, 'participants'),
        ({'sidecar': {'Hematocrit': 0.44}}, 0.44, 'Hematocrit field'),
        ({'participant_rows': 'participant_id\thaematocrit\nsub-02\t0.41\n'}, None, None),
        # Values that are not a haematocrit are ignored
        ({'participant_rows': 'participant_id\thaematocrit\nsub-01\t250\n'}, None, None),
        ({'participant_rows': 'participant_id\thaematocrit\nsub-01\tlow\n'}, None, None),
        ({'participant_rows': 'participant_id\tage\nsub-01\t30\n'}, None, None),
        ({}, None, None),
    ],
)
def test_collect_haematocrit(tmp_path, sources, expected, source):
    """A measured haematocrit is taken from the most specific source that has one (gh-295)."""
    from aslprep.utils.bids import collect_haematocrit

    layout, asl_file = _haematocrit_dataset(tmp_path / 'bids', **sources)
    metadata = layout.get_metadata(asl_file)

    haematocrit, found_in = collect_haematocrit(layout, asl_file, metadata)

    if expected is None:
        assert haematocrit is None
        assert found_in is None
    else:
        assert haematocrit == pytest.approx(expected)
        assert source in found_in


def test_collect_haematocrit_without_sessions(tmp_path):
    """Without sessions, scans.tsv sits in the subject directory."""
    from aslprep.utils.bids import collect_haematocrit

    layout, asl_file = _haematocrit_dataset(
        tmp_path / 'bids',
        sessions=False,
        scans=('haematocrit', '0.38'),
        participant_rows='participant_id\thaematocrit\nsub-01\t0.41\n',
    )

    haematocrit, found_in = collect_haematocrit(layout, asl_file, layout.get_metadata(asl_file))

    assert haematocrit == pytest.approx(0.38)
    assert found_in == 'sub-01_scans.tsv'


def test_collect_run_data_haematocrit(tmp_path, monkeypatch):
    """collect_run_data adds the haematocrit to the ASL metadata, where estimate_t1 reads it."""
    from aslprep import config
    from aslprep.utils.bids import collect_run_data

    monkeypatch.setattr(config.workflow, 'ignore', [])
    layout, asl_file = _haematocrit_dataset(
        tmp_path / 'bids',
        participant_rows='participant_id\tsex\thaematocrit\nsub-01\tF\t38.5\n',
    )

    run_data = collect_run_data(layout, asl_file)

    assert run_data['asl_metadata']['Haematocrit'] == pytest.approx(0.385)
