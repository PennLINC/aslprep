"""Tests for aslprep.utils.bids."""

import json

import pytest

from aslprep.utils.bids import collect_derivatives

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
