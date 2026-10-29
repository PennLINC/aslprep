"""Tests for measuring tissue T1s from a T1 map (gh-283)."""

import json
import shutil
from pathlib import Path

import nibabel as nb
import numpy as np
import pytest

from aslprep.utils import t1map as t1m

SHAPE = (20, 20, 20)
AFFINE = np.diag([2.0, 2.0, 2.0, 1.0])
GM_T1, WM_T1 = 1.4, 0.8


def _anatomy():
    """Gray matter at low x, white matter at high x, and a mixed band between them."""
    x = np.arange(SHAPE[0])[:, None, None] * np.ones(SHAPE)
    gm = np.clip((12 - x) / 4, 0, 1)
    wm = 1 - gm
    csf = np.zeros(SHAPE)
    # The voxel's T1 follows the mixture of rates, as in fast exchange
    t1 = 1 / (gm / GM_T1 + wm / WM_T1)
    return t1, gm, wm, csf


def _write(data, path, affine=AFFINE):
    nb.Nifti1Image(np.asarray(data, dtype=np.float32), affine).to_filename(path)
    return str(path)


@pytest.mark.parametrize(
    ('metadata', 'values', 'scale'),
    [
        ({'Units': 's'}, [1.2, 0.8], 1.0),
        ({'Units': 'ms'}, [1200.0, 800.0], 1e-3),
        ({'Units': 'Milliseconds'}, [1200.0], 1e-3),
        ({}, [1200.0, 800.0, 0.0], 1e-3),
        ({}, [1.2, 0.8, np.nan], 1.0),
    ],
)
def test_t1_scale(metadata, values, scale):
    """Units come from the sidecar, or from the magnitude of the values."""
    found, how = t1m.t1_scale(np.array(values), metadata)
    assert found == scale
    assert how


def test_t1_scale_refusals():
    with pytest.raises(ValueError, match='expected seconds or milliseconds'):
        t1m.t1_scale(np.array([1.0]), {'Units': 'Hz'})

    with pytest.raises(ValueError, match='no positive values'):
        t1m.t1_scale(np.zeros(3), {})


def test_tissue_t1_medians():
    """Pure voxels, eroded away from the mixed band, give each tissue's own T1."""
    t1, gm, wm, csf = _anatomy()
    medians = t1m.tissue_t1_medians(t1, gm, wm, csf)
    assert medians['gray']['t1'] == pytest.approx(GM_T1)
    assert medians['white']['t1'] == pytest.approx(WM_T1)
    assert medians['gray']['iqr'] == pytest.approx(0)
    assert medians['gray']['n'] >= t1m.MIN_PURE_VOXELS

    # CSF in a voxel excludes it, even where the tissue probability is high
    csf[:4] = 0.1
    fewer = t1m.tissue_t1_medians(t1, gm, wm, csf)
    assert fewer['gray']['n'] < medians['gray']['n']

    # Without enough pure voxels with a usable T1, the measurement is refused
    no_gm = np.where(gm >= t1m.PURE_TISSUE_PROBABILITY, 0.0, t1)
    with pytest.raises(ValueError, match='pure gray matter voxels'):
        t1m.tissue_t1_medians(no_gm, gm, wm, csf)


def test_implausible_t1_warnings():
    medians = {'gray': {'t1': 1.6}, 'white': {'t1': 0.3}}
    warnings = t1m.implausible_t1_warnings(medians, 3)
    assert len(warnings) == 1
    assert 'white matter' in warnings[0]
    # No defaults for this field strength, or no field strength: no warnings
    assert t1m.implausible_t1_warnings(medians, 2.89) == []
    assert t1m.implausible_t1_warnings(medians, None) == []


def test_check_t1map_field_strength():
    t1map = {'t1map': 'sub-01_T1map.nii.gz', 'metadata': {'MagneticFieldStrength': 7}}
    with pytest.raises(ValueError, match='acquired at 7 T'):
        t1m.check_t1map_field_strength(t1map, 3, 'sub-01_asl.nii.gz')

    t1m.check_t1map_field_strength(t1map, 7.0, 'sub-01_asl.nii.gz')
    t1m.check_t1map_field_strength({**t1map, 'metadata': {}}, 3, 'sub-01_asl.nii.gz')


def test_interface_defaults_can_be_documented():
    """Every ASLPrep interface input with a default can be documented by nipype's Sphinx extension.

    nipype's apidoc calls ``spec.default_value()`` on each such input. A ``default_value``
    argument to a trait is stored as metadata that hides that method, which breaks the docs build.
    """
    import importlib
    import pkgutil

    from nipype.interfaces.base import BaseInterface

    import aslprep.interfaces

    checked = 0
    for module_info in pkgutil.iter_modules(aslprep.interfaces.__path__):
        module = importlib.import_module(f'aslprep.interfaces.{module_info.name}')
        for obj in vars(module).values():
            if not (
                isinstance(obj, type)
                and issubclass(obj, BaseInterface)
                and obj.__module__ == module.__name__
                and obj.input_spec
            ):
                continue

            inputs = obj.input_spec()
            for name, spec in inputs.traits(transient=None).items():
                if spec.usedefault:
                    assert callable(spec.default_value), f'{obj.__name__}.{name}'
                    spec.default_value()
                    checked += 1

    assert checked


def test_prepare_t1map(tmp_path):
    """A map in milliseconds is converted to seconds, with R1 and zeros where T1 is unusable."""
    from aslprep.interfaces.t1map import PrepareT1Map

    values = np.full(SHAPE, 1200.0)
    values[0, 0, 0] = np.nan
    values[1, 0, 0] = -5
    in_file = _write(values[..., None], tmp_path / 'T1map.nii.gz')
    result = PrepareT1Map(in_file=in_file).run(cwd=str(tmp_path))
    t1 = nb.load(result.outputs.t1).get_fdata()
    r1 = nb.load(result.outputs.r1).get_fdata()
    assert t1.shape == SHAPE
    np.testing.assert_allclose(t1[2:], 1.2, rtol=1e-6)
    np.testing.assert_allclose(r1[2:], 1 / 1.2, rtol=1e-6)
    assert t1[0, 0, 0] == t1[1, 0, 0] == r1[0, 0, 0] == 0
    assert 'milliseconds' in result.outputs.units

    with pytest.raises(ValueError, match='must be a 3D image'):
        PrepareT1Map(in_file=_write(np.ones((4, 4, 4, 2)), tmp_path / 'bad.nii.gz')).run(
            cwd=str(tmp_path)
        )


def test_measure_tissue_t1(tmp_path):
    """The measured T1s become complete calibration updates, with their provenance."""
    from aslprep.interfaces.t1map import MeasureTissueT1
    from aslprep.utils.calibration import resolve_calibration, update_calibration

    t1, gm, wm, csf = _anatomy()
    tpms = [
        _write(data, tmp_path / f'{name}.nii.gz')
        for name, data in (('gm', gm), ('wm', wm), ('csf', csf))
    ]
    result = MeasureTissueT1(
        t1=_write(t1, tmp_path / 't1.nii.gz'),
        tpms=tpms,
        t1map_name='bids:raw:sub-01/anat/sub-01_T1map.nii.gz',
        units='Units "s" in its sidecar',
        registration='registered through UNIT1',
        field_strength=3,
    ).run(cwd=str(tmp_path))
    updates = result.outputs.updates
    assert updates['GrayMatterT1']['Value'] == pytest.approx(GM_T1)
    assert updates['WhiteMatterT1']['Value'] == pytest.approx(WM_T1)
    for entry in updates.values():
        assert entry['Source'] == 'measured'
        assert 'bids:raw:sub-01/anat/sub-01_T1map.nii.gz' in entry['Origin']
        assert 'registered through UNIT1' in entry['Origin']

    assert result.outputs.warnings == []

    # The updates apply to a tissue-T1 specification
    metadata = {
        'ArterialSpinLabelingType': 'PCASL',
        'M0Type': 'Separate',
        'MagneticFieldStrength': 3,
        'BackgroundSuppression': False,
        'PostLabelingDelay': [0.5, 1.0, 1.5, 2.0],
        'LabelingDuration': 1.8,
    }
    calibration = resolve_calibration(metadata, m0tr=6.0, kinetic_model='tissue-t1')
    updated = update_calibration(calibration, updates)
    assert updated['Parameters']['GrayMatterT1']['Value'] == pytest.approx(GM_T1)

    off_grid = _write(gm, tmp_path / 'off_grid.nii.gz', affine=np.eye(4))
    with pytest.raises(ValueError, match='is not on the T1 grid'):
        MeasureTissueT1(
            t1=_write(t1, tmp_path / 't1.nii.gz'),
            tpms=[off_grid, *tpms[1:]],
            t1map_name='x',
            units='s',
            registration='x',
        ).run(cwd=str(tmp_path))


def test_runtime_calibration_measured_t1(tmp_path):
    """Measured T1s and their warnings reach the specification, and the kinetic model's T1'."""
    from aslprep.interfaces.cbf import ResolveRuntimeCalibration
    from aslprep.utils.calibration import (
        measured_tissue_t1_updates,
        resolve_calibration,
        tissue_t1_kinetics,
    )
    from aslprep.utils.cbf import apparent_tissue_t1

    metadata = {
        'ArterialSpinLabelingType': 'PCASL',
        'M0Type': 'Separate',
        'MagneticFieldStrength': 3,
        'BackgroundSuppression': False,
        'PostLabelingDelay': [0.5, 1.0, 1.5, 2.0],
        'LabelingDuration': 1.8,
    }
    calibration = resolve_calibration(metadata, m0tr=6.0, kinetic_model='tissue-t1')
    updates = measured_tissue_t1_updates(1.3, 0.8, {'gray': 'GM map', 'white': 'WM map'})
    result = ResolveRuntimeCalibration(
        calibration=calibration,
        updates=updates,
        warnings=['A warning'],
        asl_type='PCASL',
        multi_delay=True,
    ).run(cwd=str(tmp_path))
    runtime = result.outputs.calibration
    assert runtime['Warnings'][-1] == 'A warning'
    np.testing.assert_allclose(
        tissue_t1_kinetics(runtime)['t1prime'],
        apparent_tissue_t1(np.array([1.3, 0.8]), np.array([60.0, 20.0])),
    )
    assert result.outputs.model_metadata['aslprep']['Parameters']['GrayMatterT1'] == 1.3


def _bids_dataset(root, name='raw'):
    root.mkdir(parents=True)
    description = {'Name': name, 'BIDSVersion': '1.10.0'}
    if name != 'raw':
        description.update(DatasetType='derivative', GeneratedBy=[{'Name': name}])
    (root / 'dataset_description.json').write_text(json.dumps(description))
    anat = root / 'sub-01' / 'anat'
    anat.mkdir(parents=True)
    return anat


def test_collect_t1map(tmp_path):
    """Derivative maps come first, a companion must share the grid, and duplicates are refused."""
    from bids.layout import BIDSLayout

    from aslprep.utils.bids import collect_t1map

    anat = _bids_dataset(tmp_path / 'raw')
    _write(np.ones(SHAPE), anat / 'sub-01_T1map.nii.gz')
    (anat / 'sub-01_T1map.json').write_text(json.dumps({'Units': 'ms'}))
    _write(np.ones(SHAPE), anat / 'sub-01_UNIT1.nii.gz')
    # A T1w on another grid is not a companion
    _write(np.ones((10, 10, 10)), anat / 'sub-01_T1w.nii.gz')
    layout = BIDSLayout(str(tmp_path / 'raw'), validate=False)

    found = collect_t1map(layout, {}, '01')
    assert found['t1map'].endswith('sub-01_T1map.nii.gz')
    assert found['metadata']['Units'] == 'ms'
    assert found['uri'] == 'bids:raw:sub-01/anat/sub-01_T1map.nii.gz'
    assert found['companion'].endswith('sub-01_UNIT1.nii.gz')

    # A derivative map is preferred, and the raw UNIT1 on its grid is still its companion
    deriv_anat = _bids_dataset(tmp_path / 'qmri', name='qmri')
    _write(np.ones(SHAPE), deriv_anat / 'sub-01_desc-b1corrected_T1map.nii.gz')
    found = collect_t1map(layout, {'qmri': str(tmp_path / 'qmri')}, '01')
    assert found['uri'] == 'bids:qmri:sub-01/anat/sub-01_desc-b1corrected_T1map.nii.gz'
    assert found['companion'].endswith('sub-01_UNIT1.nii.gz')

    # Without a companion on the grid, there is none
    (anat / 'sub-01_UNIT1.nii.gz').unlink()
    layout = BIDSLayout(str(tmp_path / 'raw'), validate=False)
    assert collect_t1map(layout, {}, '01')['companion'] is None

    # Two maps in one dataset are refused, unless a filter selects one
    _write(np.ones(SHAPE), anat / 'sub-01_acq-other_T1map.nii.gz')
    layout = BIDSLayout(str(tmp_path / 'raw'), validate=False)
    with pytest.raises(ValueError, match='Found 2 T1 maps'):
        collect_t1map(layout, {}, '01')

    found = collect_t1map(layout, {}, '01', bids_filters={'t1map': {'acquisition': 'other'}})
    assert found['t1map'].endswith('sub-01_acq-other_T1map.nii.gz')

    assert collect_t1map(layout, {}, '02') is None


def test_collect_t1map_space_and_symlinks(tmp_path):
    """Maps in a template space are not used, and symlinked maps keep their dataset path."""
    from bids.layout import BIDSLayout

    from aslprep.utils.bids import collect_t1map

    anat = _bids_dataset(tmp_path / 'raw')
    deriv_anat = _bids_dataset(tmp_path / 'qmri', name='qmri')
    # A map resampled to a template cannot be registered rigidly to the participant
    _write(np.ones(SHAPE), deriv_anat / 'sub-01_space-MNI152NLin2009cAsym_T1map.nii.gz')
    # A DataLad dataset stores its files as symlinks into an annex outside the dataset tree
    annex = tmp_path / 'annex'
    annex.mkdir()
    target = _write(np.ones(SHAPE), annex / 'object.nii.gz')
    (anat / 'sub-01_T1map.nii.gz').symlink_to(target)
    layout = BIDSLayout(str(tmp_path / 'raw'), validate=False)

    found = collect_t1map(layout, {'qmri': str(tmp_path / 'qmri')}, '01')
    assert found['uri'] == 'bids:raw:sub-01/anat/sub-01_T1map.nii.gz'


def test_collect_t1map_pymp2rage(tmp_path, monkeypatch):
    """pymp2rage-style maps, all with a space entity, are listed in a warning or selected by filter.

    pymp2rage writes each map in the MP2RAGE, T1w and MNI spaces, with and without B1
    correction, in each session, with a UNIT1 in each space.
    """
    from bids.layout import BIDSLayout

    from aslprep import config
    from aslprep.utils.bids import collect_t1map

    warnings = []
    monkeypatch.setattr(config.loggers.workflow, 'warning', warnings.append)
    _bids_dataset(tmp_path / 'raw')
    layout = BIDSLayout(str(tmp_path / 'raw'), validate=False)
    root = tmp_path / 'pymp2rage'
    _bids_dataset(root, name='pymp2rage')
    grids = {
        'MP2RAGE': AFFINE,
        'T1w': np.diag([2.0, 2.0, 2.0, 1.0]) + np.array([[0, 0, 0, 5.0]] * 3 + [[0] * 4]),
        'MNI152NLin2009cAsym': np.eye(4),
    }
    for session in ('01', '02'):
        anat = root / 'sub-01' / f'ses-{session}' / 'anat'
        anat.mkdir(parents=True)
        for space, affine in grids.items():
            prefix = f'sub-01_ses-{session}_run-01_space-{space}'
            for desc in ('', '_desc-B1corrected'):
                _write(np.ones(SHAPE), anat / f'{prefix}{desc}_T1map.nii.gz', affine=affine)

            _write(np.ones(SHAPE), anat / f'{prefix}_desc-B1correctedbrain_UNIT1.nii.gz', affine)

    derivatives = {'pymp2rage': str(root)}
    # By default, none is used, and the warning names them
    assert collect_t1map(layout, derivatives, '01') is None
    assert len(warnings) == 1
    assert 'has 12 T1 map(s) that were not used' in warnings[0]
    assert 'bids:pymp2rage:sub-01/ses-01/anat/' in warnings[0]

    # A filter selects one, and its companion is the UNIT1 in the same space
    t1w_filter = {'t1map': {'space': 'T1w', 'desc': 'B1corrected', 'session': '01'}}
    found = collect_t1map(layout, derivatives, '01', bids_filters=t1w_filter)
    assert found['uri'] == (
        'bids:pymp2rage:sub-01/ses-01/anat/'
        'sub-01_ses-01_run-01_space-T1w_desc-B1corrected_T1map.nii.gz'
    )
    assert found['companion'].endswith(
        'sub-01_ses-01_run-01_space-T1w_desc-B1correctedbrain_UNIT1.nii.gz'
    )

    # With a session-level anatomical reference, the session comes from the subject workflow
    session_filter = {'t1map': {'space': 'MP2RAGE', 'desc': 'B1corrected'}}
    found = collect_t1map(layout, derivatives, '01', session_id='02', bids_filters=session_filter)
    assert 'ses-02' in found['uri']
    assert 'space-MP2RAGE_desc-B1correctedbrain_UNIT1' in found['companion']

    # Without a session, both sessions match
    with pytest.raises(ValueError, match='Found 2 T1 maps'):
        collect_t1map(layout, derivatives, '01', bids_filters=session_filter)

    # A template-space map cannot be used, even when a filter selects it
    mni_filter = {
        't1map': {'space': 'MNI152NLin2009cAsym', 'desc': 'B1corrected', 'session': '01'}
    }
    with pytest.raises(ValueError, match='template space MNI152NLin2009cAsym'):
        collect_t1map(layout, derivatives, '01', bids_filters=mni_filter)


def test_init_asl_wf_t1map_unused():
    """A run that cannot use the T1 map (here, with the consensus model) does not ask for it."""
    from aslprep.tests.tests import mock_config, reset_config
    from aslprep.workflows.asl.base import init_asl_wf

    reset_config()
    with mock_config():
        from aslprep import config

        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        wf = init_asl_wf(
            asl_file=str(asl_file),
            t1map={
                't1map': 'sub-01_T1map.nii.gz',
                'metadata': {'MagneticFieldStrength': 7},
                'uri': 'bids:raw:sub-01/anat/sub-01_T1map.nii.gz',
                'companion': None,
            },
        )

    reset_config()
    assert wf.t1map_field_strength is None


@pytest.mark.parametrize(('identity', 'companion'), [(True, True), (False, True), (False, False)])
def test_init_t1map_wf(identity, companion):
    """Only a map outside the anatomical image's space is registered, through its companion."""
    from aslprep.tests.tests import mock_config
    from aslprep.workflows.t1map import init_t1map_wf

    t1map = {
        't1map': 'sub-01_T1map.nii.gz',
        'metadata': {'Units': 's'},
        'uri': 'bids:raw:sub-01/anat/sub-01_T1map.nii.gz',
        'companion': 'sub-01_UNIT1.nii.gz' if companion else None,
    }
    with mock_config():
        wf = init_t1map_wf(t1map=t1map, identity=identity, output_dir='.')

    register = wf.get_node('register_t1map')
    assert (register is None) == identity
    assert 'bids:raw:sub-01/anat/sub-01_T1map.nii.gz' in wf.__desc__
    if register is not None:
        # The companion if there is one, otherwise the map's R1
        source = wf.get_node('inputnode') if companion else wf.get_node('prepare_t1map')
        edge = wf._graph.get_edge_data(source, register)['connect']
        assert (('companion' if companion else 'r1'), 'moving_image') in edge


@pytest.mark.parametrize(
    ('name', 'prefix'),
    [
        ('sub-01/anat/sub-01_T1map.nii.gz', 'sub-01'),
        ('sub-01/anat/sub-01_desc-b1corrected_T1map.nii.gz', 'sub-01'),
        # A pymp2rage map selected with --bids-filter-file
        (
            'sub-01/ses-01/anat/sub-01_ses-01_run-01_space-T1w_desc-B1corrected_T1map.nii.gz',
            'sub-01_ses-01_run-01',
        ),
    ],
)
def test_t1map_sinks(tmp_path, name, prefix):
    """The transform and reportlet sinks, as the T1 map workflow configures them, name outputs."""
    from aslprep.tests.tests import mock_config
    from aslprep.workflows.t1map import init_t1map_wf

    source = tmp_path / name
    source.parent.mkdir(parents=True)
    source = _write(np.ones(SHAPE), source)
    xfm = tmp_path / 'identity.txt'
    xfm.write_text('#Insight Transform File V1.0\n')
    svg = tmp_path / 'report.svg'
    svg.write_text('<svg></svg>')
    with mock_config():
        wf = init_t1map_wf(
            t1map={'t1map': source, 'metadata': {}, 'uri': 'x', 'companion': None},
            identity=True,
            output_dir=str(tmp_path / 'out'),
        )
        out_files = {}
        for node_name, in_file in (('ds_t1map2anat_xfm', xfm), ('ds_t1map_report', svg)):
            sink = wf.get_node(node_name)
            sink.inputs.source_file = source
            sink.inputs.in_file = str(in_file)
            sink.base_dir = str(tmp_path / 'work')
            out_files[node_name] = Path(sink.run().outputs.out_file)

    xfm_out = out_files['ds_t1map2anat_xfm']
    assert xfm_out.name == f'{prefix}_from-T1map_to-T1w_mode-image_xfm.txt'
    assert xfm_out.parent.name == 'anat'
    report_out = out_files['ds_t1map_report']
    assert report_out.name == f'{prefix}_desc-coreg_T1map.svg'
    assert report_out.parent.name == 'figures'


def test_t1map_report(tmp_path):
    """The registration reportlet compares the anatomical reference with the map's R1."""
    from aslprep.tests.tests import mock_config
    from aslprep.workflows.t1map import init_t1map_wf

    t1, _, wm, _ = _anatomy()
    with mock_config():
        wf = init_t1map_wf(
            t1map={'t1map': 'x.nii.gz', 'metadata': {}, 'uri': 'x', 'companion': None},
            identity=True,
            output_dir=str(tmp_path / 'out'),
        )
        report = wf.get_node('t1map_report')
        report.inputs.before = _write(1 / t1, tmp_path / 't1w.nii.gz')
        report.inputs.after = _write(1 / t1, tmp_path / 'r1.nii.gz')
        report.inputs.wm_seg = _write(wm, tmp_path / 'wm.nii.gz')
        report.base_dir = str(tmp_path / 'work')
        out_report = Path(report.run().outputs.out_report)

    assert out_report.suffix == '.svg'
    assert 'T1 map (R1)' in out_report.read_text()


@pytest.mark.skipif(
    shutil.which('antsApplyTransforms') is None,
    reason='ANTs is not installed',
)
def test_t1map_wf_measures_identity(tmp_path):
    """Run end to end: a map in the anatomical image's space gives each tissue's T1."""
    from aslprep.tests.tests import mock_config
    from aslprep.workflows.t1map import init_t1map_wf

    t1, gm, wm, csf = _anatomy()
    anat = tmp_path / 'sub-01' / 'anat'
    anat.mkdir(parents=True)
    t1map_file = _write(t1 * 1000, anat / 'sub-01_T1map.nii.gz')
    t1w = _write(1 / t1, anat / 'sub-01_T1w.nii.gz')
    with mock_config():
        wf = init_t1map_wf(
            t1map={
                't1map': t1map_file,
                'metadata': {},
                'uri': 'bids:raw:sub-01/anat/sub-01_T1map.nii.gz',
                'companion': t1w,
            },
            identity=True,
            output_dir=str(tmp_path / 'out'),
        )
        wf.base_dir = str(tmp_path / 'work')
        wf.inputs.inputnode.t1w_preproc = t1w
        wf.inputs.inputnode.t1w_mask = _write(np.ones(SHAPE), tmp_path / 'mask.nii.gz')
        wf.inputs.inputnode.t1w_tpms = [
            _write(data, tmp_path / f'{name}.nii.gz')
            for name, data in (('gm', gm), ('wm', wm), ('csf', csf))
        ]
        result = wf.run()

    (measure,) = [node for node in result.nodes() if node.name == 'measure_tissue_t1']
    updates = measure.result.outputs.updates
    assert updates['GrayMatterT1']['Value'] == pytest.approx(GM_T1, rel=1e-4)
    assert updates['WhiteMatterT1']['Value'] == pytest.approx(WM_T1, rel=1e-4)
    assert list(Path(tmp_path / 'out').rglob('*from-T1map_to-T1w_mode-image_xfm.txt'))
