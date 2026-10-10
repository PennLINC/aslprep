"""ASLPrep's BIDS collection on the aslscan fixtures (spec Section 4.8, fast tier).

Each fixture must present exactly one ASL run, with the aslcontext and M0 ASLPrep expects.
"""

import re

import pytest
from bids.layout import BIDSLayout, BIDSLayoutIndexer
from niworkflows.utils.bids import collect_data

from aslprep.tests import aslscan_fixtures as af
from aslprep.utils.bids import collect_run_data

#: The queries of aslprep.workflows.base.init_single_subject_wf
QUERIES = {
    'fmap': {'datatype': 'fmap'},
    'flair': {'datatype': 'anat', 'suffix': 'FLAIR'},
    't2w': {'datatype': 'anat', 'suffix': 'T2w'},
    't1w': {'datatype': 'anat', 'suffix': 'T1w'},
    'roi': {'datatype': 'anat', 'suffix': 'roi'},
    'sbref': {'datatype': 'perf', 'suffix': 'sbref'},
    'asl': {'datatype': 'perf', 'suffix': 'asl'},
}


def _layout(path):
    """A layout indexed as aslprep.config.execution builds it."""
    indexer = BIDSLayoutIndexer(
        validate=False,
        ignore=[
            'code',
            'stimuli',
            'sourcedata',
            'models',
            re.compile(r'\/\.\w+|^\.\w+'),
            re.compile(r'sub-[a-zA-Z0-9]+(/ses-[a-zA-Z0-9]+)?/(beh|dwi|eeg|ieeg|meg|func)'),
        ],
    )
    return BIDSLayout(str(path), indexer=indexer)


def _fixture(name, data_dir):
    try:
        return af.fixture_dir(name, data_dir)
    except af.FixtureUnavailable as exc:
        pytest.skip(str(exc))


@pytest.mark.parametrize('name', sorted(af.RECIPES))
def test_fixture_presents_one_asl_run(name, data_dir, monkeypatch):
    from aslprep import config

    monkeypatch.setattr(config.workflow, 'ignore', [])  # as parse_args leaves it by default
    fixture = _fixture(name, data_dir)
    recipe = af.load_recipe(name)
    layout = _layout(fixture)
    subject_data, _ = collect_data(layout, '01', queries=QUERIES, bids_validate=False)

    assert [p.split('/')[-1] for p in subject_data['asl']] == [
        f'sub-01_acq-{recipe.acq}_asl.nii.gz'
    ]
    expected_t1w = [] if recipe.anat == 'none' else ['sub-01_T1w.nii.gz']
    assert [p.split('/')[-1] for p in subject_data['t1w']] == expected_t1w

    run_data = collect_run_data(layout, subject_data['asl'][0])
    assert run_data['aslcontext'].endswith(f'sub-01_acq-{recipe.acq}_aslcontext.tsv')
    metadata = layout.get_metadata(subject_data['asl'][0])
    assert metadata['M0Type'] == recipe.sidecar['M0Type']
    if metadata['M0Type'] == 'Separate':
        assert run_data['m0scan'].endswith(f'sub-01_acq-{recipe.acq}_m0scan.nii.gz')
        assert run_data['m0scan_metadata']['RepetitionTimePreparation'] > 0
    else:
        assert not run_data.get('m0scan')
    if not recipe.keep_labeling_efficiency:
        assert 'LabelingEfficiency' not in metadata
