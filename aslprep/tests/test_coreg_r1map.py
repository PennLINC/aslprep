"""Tests for registering R1, from a T1 or R1 map on the grid of the M0 scan, to the anatomy."""

import json

import nibabel as nb
import numpy as np
import pytest

from aslprep.utils.t1map import coregistration_r1, r1_scale

SHAPE = (6, 6, 4)
AFFINE = np.diag([3.0, 3.0, 5.0, 1.0])
URI = 'bids:raw:sub-01/perf/sub-01_acq-tr1_T1map.nii.gz'


def _write(data, path, affine=AFFINE):
    nb.Nifti1Image(np.asarray(data, dtype=np.float32), affine).to_filename(path)
    return str(path)


@pytest.mark.parametrize(
    ('values', 'metadata', 'suffix', 'expected'),
    [
        # T1 in milliseconds and seconds, from the sidecar or the values
        ([1000.0, 1600.0], {'Units': 'ms'}, 'T1map', [1.0, 0.625]),
        ([1000.0, 1600.0], {}, 'T1map', [1.0, 0.625]),
        ([1.0, 1.6], {'Units': 's'}, 'T1map', [1.0, 0.625]),
        ([1.0, 1.6], {}, 'T1map', [1.0, 0.625]),
        # R1 per second and per millisecond, from the sidecar or the values
        ([1.0, 0.625], {'Units': 's-1'}, 'R1map', [1.0, 0.625]),
        ([1.0, 0.625], {'Units': '1/s'}, 'R1map', [1.0, 0.625]),
        ([1.0, 0.625], {}, 'R1map', [1.0, 0.625]),
        ([0.001, 0.000625], {'Units': 'ms-1'}, 'R1map', [1.0, 0.625]),
        ([0.001, 0.000625], {}, 'R1map', [1.0, 0.625]),
    ],
)
def test_coregistration_r1_units(values, metadata, suffix, expected):
    """T1 and R1 maps, in any of their units, give R1 per second."""
    r1, _ = coregistration_r1(np.array(values), metadata, suffix)
    np.testing.assert_allclose(r1, expected)


@pytest.mark.parametrize('suffix', ['T1map', 'R1map'])
def test_coregistration_r1_invalid(suffix):
    """Voxels without a finite, positive value have an R1 of zero."""
    values = np.array([1.6, 0.0, -1.0, np.nan, np.inf])
    r1, _ = coregistration_r1(values, {'Units': 's' if suffix == 'T1map' else 's-1'}, suffix)
    expected = 1 / 1.6 if suffix == 'T1map' else 1.6
    np.testing.assert_allclose(r1, [expected, 0.0, 0.0, 0.0, 0.0])


def test_coregistration_r1_errors():
    """Unknown units, suffixes and empty maps are refused."""
    with pytest.raises(ValueError, match='expected per second'):
        r1_scale(np.ones(3), {'Units': 'ms'})

    with pytest.raises(ValueError, match='no positive values'):
        r1_scale(np.zeros(3), {})

    with pytest.raises(ValueError, match='expected seconds'):
        coregistration_r1(np.ones(3), {'Units': 's-1'}, 'T1map')

    with pytest.raises(ValueError, match='Expected a T1map or an R1map'):
        coregistration_r1(np.ones(3), {}, 'M0map')


def test_coregistration_r1_interface(tmp_path):
    """The interface writes R1 on the map's grid, from a 3D map or one with one volume."""
    from aslprep.interfaces.t1map import CoregistrationR1

    t1_ms = np.full(SHAPE, 1600.0)
    t1_ms[0] = 800.0
    t1_ms[1] = 0.0
    t1_ms[2] = np.nan
    t1map = _write(t1_ms[..., None], tmp_path / 'T1map.nii.gz')

    result = CoregistrationR1(
        in_file=t1map,
        metadata={'Units': 'ms'},
        suffix='T1map',
    ).run(cwd=str(tmp_path))
    img = nb.load(result.outputs.out_file)
    r1 = img.get_fdata()

    assert r1.shape == SHAPE
    np.testing.assert_allclose(img.affine, AFFINE)
    assert r1[0, 0, 0] == pytest.approx(1.25)
    assert r1[3, 0, 0] == pytest.approx(0.625)
    assert r1[1, 0, 0] == 0
    assert r1[2, 0, 0] == 0
    # White matter (short T1) is brighter than gray matter in R1
    assert r1[0, 0, 0] > r1[3, 0, 0]

    series = _write(np.ones(SHAPE + (2,)), tmp_path / 'series.nii.gz')
    with pytest.raises(ValueError, match='must be a 3D image'):
        CoregistrationR1(in_file=series, suffix='R1map').run(cwd=str(tmp_path))


def _bids_dataset(root, name='raw'):
    root.mkdir(parents=True)
    description = {'Name': name, 'BIDSVersion': '1.10.0'}
    if name != 'raw':
        description.update(DatasetType='derivative', GeneratedBy=[{'Name': name}])
    (root / 'dataset_description.json').write_text(json.dumps(description))
    subject = root / 'sub-01'
    (subject / 'anat').mkdir(parents=True)
    (subject / 'perf').mkdir(parents=True)
    return subject


def test_collect_m0scan_r1map(tmp_path):
    """R1 maps are found only when asked for, and T1 maps win within a dataset."""
    from bids.layout import BIDSLayout

    from aslprep.utils.bids import collect_m0scan_t1map

    suffixes = ('T1map', 'R1map')
    raw = _bids_dataset(tmp_path / 'raw')
    m0scan = _write(np.ones(SHAPE), raw / 'perf' / 'sub-01_m0scan.nii.gz')
    _write(np.ones(SHAPE), raw / 'perf' / 'sub-01_acq-tdp_R1map.nii.gz')
    (raw / 'perf' / 'sub-01_acq-tdp_R1map.json').write_text(json.dumps({'Units': 's-1'}))
    layout = BIDSLayout(str(tmp_path / 'raw'), validate=False)

    # The M0 recovery correction only looks for T1 maps
    assert collect_m0scan_t1map(layout, {}, m0scan) is None
    found = collect_m0scan_t1map(layout, {}, m0scan, suffixes=suffixes)
    assert found['suffix'] == 'R1map'
    assert found['t1map'].endswith('sub-01_acq-tdp_R1map.nii.gz')
    assert found['metadata']['Units'] == 's-1'
    assert found['uri'] == 'bids:raw:sub-01/perf/sub-01_acq-tdp_R1map.nii.gz'

    # A T1 map in the same dataset wins
    _write(np.ones(SHAPE), raw / 'anat' / 'sub-01_acq-tdp_T1map.nii.gz')
    layout = BIDSLayout(str(tmp_path / 'raw'), validate=False)
    found = collect_m0scan_t1map(layout, {}, m0scan, suffixes=suffixes)
    assert found['suffix'] == 'T1map'
    assert found['uri'] == 'bids:raw:sub-01/anat/sub-01_acq-tdp_T1map.nii.gz'

    # A derivative R1 map wins over a raw T1 map, as derivatives are searched first
    deriv = _bids_dataset(tmp_path / 'qmri', name='qmri')
    _write(np.ones(SHAPE), deriv / 'perf' / 'sub-01_desc-tdp_R1map.nii.gz')
    derivatives = {'qmri': str(tmp_path / 'qmri')}
    found = collect_m0scan_t1map(layout, derivatives, m0scan, suffixes=suffixes)
    assert found['suffix'] == 'R1map'
    assert found['uri'] == 'bids:qmri:sub-01/perf/sub-01_desc-tdp_R1map.nii.gz'

    # Two R1 maps on the grid are refused, unless a filter selects one
    # (in a new dataset, as derivative layouts are cached)
    deriv = _bids_dataset(tmp_path / 'qmri2', name='qmri2')
    _write(np.ones(SHAPE), deriv / 'perf' / 'sub-01_desc-tdp_R1map.nii.gz')
    _write(np.ones(SHAPE), deriv / 'perf' / 'sub-01_desc-other_R1map.nii.gz')
    derivatives = {'qmri2': str(tmp_path / 'qmri2')}
    with pytest.raises(ValueError, match='Found 2 R1 maps on the grid'):
        collect_m0scan_t1map(layout, derivatives, m0scan, suffixes=suffixes)

    found = collect_m0scan_t1map(
        layout,
        derivatives,
        m0scan,
        bids_filters={'m0scan_t1map': {'desc': 'other'}},
        suffixes=suffixes,
    )
    assert found['t1map'].endswith('sub-01_desc-other_R1map.nii.gz')


@pytest.mark.parametrize(
    ('kwargs', 'expected'),
    [
        ({}, True),
        # Another reference volume type: the M0 transform is between different contrasts
        ({'reference_volume_type': 'control'}, False),
        # The transform is reused, so there is no registration
        ({'aslref2anat_xform': 'from-aslref_to-T1w_mode-image_xfm.txt'}, False),
        ({'coreg_r1map': None}, False),
    ],
)
def test_use_coregistration_r1(tmp_path, kwargs, expected):
    """R1 is registered only when its transform into the reference's frame is sound."""
    from aslprep.workflows.asl.fit import use_coregistration_r1

    arguments = {
        'coreg_r1map': {'uri': URI},
        'm0scan': _write(np.ones(SHAPE), tmp_path / 'm0scan.nii.gz'),
        'reference_volume_type': 'separate_m0scan',
        'aslref2anat_xform': None,
        **kwargs,
    }
    coreg_r1map = arguments.pop('coreg_r1map')
    assert use_coregistration_r1(coreg_r1map, **arguments) is expected


def test_use_coregistration_r1_multivolume(tmp_path):
    """An M0 scan with several volumes has several transforms, so R1 is not registered."""
    from aslprep.workflows.asl.fit import use_coregistration_r1

    m0scan = _write(np.ones(SHAPE + (2,)), tmp_path / 'm0scan.nii.gz')
    assert not use_coregistration_r1(
        {'uri': URI},
        m0scan=m0scan,
        reference_volume_type='separate_m0scan',
        aslref2anat_xform=None,
    )
    single = _write(np.ones(SHAPE + (1,)), tmp_path / 'm0scan_single.nii.gz')
    assert use_coregistration_r1(
        {'uri': URI},
        m0scan=single,
        reference_volume_type='separate_m0scan',
        aslref2anat_xform=None,
    )


@pytest.mark.parametrize('freesurfer', [False, True])
def test_reverse_bbr_contrast(freesurfer):
    """BBR is set for white matter brighter than gray matter, in FLIRT or bbregister."""
    from fmriprep.workflows.bold.registration import init_bold_reg_wf

    from aslprep.workflows.asl.fit import reverse_bbr_contrast

    asl_reg_wf = init_bold_reg_wf(
        freesurfer=freesurfer,
        use_bbr=None,
        bold2anat_dof=6,
        bold2anat_init='t1w',
        mem_gb=1,
        omp_nthreads=1,
    )
    reverse_bbr_contrast(asl_reg_wf)
    if freesurfer:
        assert asl_reg_wf.get_node('bbreg_wf').get_node('bbregister').inputs.contrast_type == 't1'
    else:
        assert asl_reg_wf.get_node('fsl_bbr_wf').get_node('flt_bbr').inputs.bbrslope == 0.5

    # Without BBR there is nothing to reverse
    reverse_bbr_contrast(
        init_bold_reg_wf(
            freesurfer=False,
            use_bbr=False,
            bold2anat_dof=6,
            bold2anat_init='t1w',
            mem_gb=1,
            omp_nthreads=1,
        )
    )


@pytest.mark.parametrize('use_r1', [True, False])
def test_init_asl_fit_wf_r1(tmp_path, monkeypatch, use_r1):
    """R1 replaces the coregistration reference as the registered image, and only then."""
    from aslprep.tests.tests import mock_config, reset_config
    from aslprep.workflows.asl import fit

    m0scan = _write(np.ones(SHAPE), tmp_path / 'sub-01_m0scan.nii.gz')
    t1map = _write(np.full(SHAPE, 1400.0), tmp_path / 'sub-01_T1map.nii.gz')
    coreg_r1map = {'t1map': t1map, 'suffix': 'T1map', 'metadata': {'Units': 'ms'}, 'uri': URI}
    monkeypatch.setattr(
        fit,
        'select_reference_volume_type',
        lambda **kwargs: 'separate_m0scan' if use_r1 else 'm0scan',
    )

    reset_config()
    with mock_config():
        from aslprep import config

        config.workflow.force = []
        # Production's FSL chain: mri_coreg, then FLIRT BBR
        config.workflow.run_reconall = False
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        wf = fit.init_asl_fit_wf(
            asl_file=str(asl_file),
            aslcontext=str(asl_file).replace('_asl.nii.gz', '_aslcontext.tsv'),
            m0scan=m0scan,
            use_ge=False,
            coreg_r1map=coreg_r1map,
        )

    reset_config()
    asl_reg_wf = wf.get_node('bold_reg_wf')
    asl_r1_wf = wf.get_node('asl_r1_wf')
    bbr_wf = asl_reg_wf.get_node('fsl_bbr_wf')
    flt_bbr = bbr_wf.get_node('flt_bbr')
    regref_edges = wf._graph.get_edge_data(wf.get_node('regref_buffer'), asl_reg_wf)
    ds_edges = wf._graph.get_edge_data(wf.get_node('regref_buffer'), wf.get_node('ds_aslreg_wf'))
    # The ASL reference still names the transform
    assert ('aslref', 'inputnode.source_files') in ds_edges['connect']
    assert ('R1 (1 / T1)' in bbr_wf.__desc__) == use_r1
    if not use_r1:
        assert asl_r1_wf is None
        assert ('aslref', 'inputnode.ref_bold_brain') in regref_edges['connect']
        assert flt_bbr.inputs.bbrslope != 0.5
        return

    assert regref_edges is None
    edges = wf._graph.get_edge_data(asl_r1_wf, asl_reg_wf)
    assert ('outputnode.r1_brain', 'inputnode.ref_bold_brain') in edges['connect']
    assert flt_bbr.inputs.bbrslope == 0.5
    assert URI in bbr_wf.__desc__
    assert asl_r1_wf.get_node('coreg_r1').inputs.in_file == t1map
    # R1 is moved with the M0 scan's transform, and masked with the ASL-estimated mask
    edges = wf._graph.get_edge_data(wf.get_node('m0scanreg_buffer'), asl_r1_wf)
    assert ('m0scan2aslref_xfm', 'inputnode.m0scan2aslref_xfm') in edges['connect']
    edges = wf._graph.get_edge_data(wf.get_node('aslstrip_buffer'), asl_r1_wf)
    assert ('aslmask', 'inputnode.asl_mask') in edges['connect']
    edges = asl_r1_wf._graph.get_edge_data(
        asl_r1_wf.get_node('aslref_r1'),
        asl_r1_wf.get_node('mask_r1'),
    )
    assert ('out_file', 'in_file') in edges['connect']
    # The transform's sidecar names the registered image
    r1_metadata = wf.get_node('r1_metadata')
    assert r1_metadata.inputs.update['RegistrationImage'] == 'R1map'
    assert r1_metadata.inputs.update['RegistrationImageSources'] == [URI]
    edges = wf._graph.get_edge_data(r1_metadata, wf.get_node('ds_aslreg_wf'))
    assert ('out_dict', 'inputnode.metadata') in edges['connect']


@pytest.mark.parametrize(
    ('ignore', 'm0type', 'expected'),
    [
        ([], 'Separate', True),
        (['t1map'], 'Separate', False),
    ],
)
def test_init_asl_wf_coreg_r1map(tmp_path, monkeypatch, ignore, m0type, expected):
    """A map for R1 is looked for, whatever the M0 TR, unless --ignore t1map is given."""
    from aslprep.tests.tests import mock_config, reset_config
    from aslprep.utils import bids
    from aslprep.workflows.asl import base

    collect_run_data = base.collect_run_data
    m0scan = _write(np.ones(SHAPE), tmp_path / 'sub-01_m0scan.nii.gz')
    calls = []

    def _separate_m0(layout, asl_file):
        run_data = collect_run_data(layout, asl_file)
        run_data['asl_metadata'] = {**run_data['asl_metadata'], 'M0Type': m0type}
        run_data['m0scan'] = m0scan
        # A TR long enough that the M0 scan needs no recovery correction
        run_data['m0scan_metadata'] = {'RepetitionTimePreparation': 6.0}
        return run_data

    def _collect(*args, suffixes=('T1map',), **kwargs):
        calls.append(suffixes)

    monkeypatch.setattr(base, 'collect_run_data', _separate_m0)
    monkeypatch.setattr(bids, 'collect_m0scan_t1map', _collect)

    reset_config()
    with mock_config():
        from aslprep import config

        config.workflow.ignore = ignore
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        base.init_asl_wf(asl_file=str(asl_file))

    reset_config()
    assert (('T1map', 'R1map') in calls) == expected
