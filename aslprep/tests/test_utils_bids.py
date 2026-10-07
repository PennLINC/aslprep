"""Tests for aslprep.utils.bids."""

import json

import pytest

from aslprep.utils.bids import collect_derivatives

XFM_NAMES = {
    'desc': {
        'aslref2anat': 'sub-01_from-aslref_to-T1w_mode-image_desc-coreg_xfm.txt',
        'aslref2fmap': 'sub-01_from-aslref_to-auto00000_mode-image_desc-fmap_xfm.txt',
    },
    'legacy': {
        'aslref2anat': 'sub-01_from-aslref_to-T1w_mode-image_xfm.txt',
        'aslref2fmap': 'sub-01_from-aslref_to-auto00000_mode-image_xfm.txt',
    },
}


@pytest.mark.parametrize('naming', ['desc', 'legacy'])
def test_collect_derivatives_transforms(tmp_path, naming):
    """Precomputed ASL transforms are found with or without desc entities.

    nipreps/fmriprep#3532 added desc-coreg and desc-fmap to these transforms,
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
