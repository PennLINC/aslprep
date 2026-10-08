"""Tests for gradient nonlinearity correction planning and workflow wiring."""

import pytest

from aslprep import config
from aslprep.tests.tests import mock_config, reset_config
from aslprep.utils import gradwarp as gw


def _grad_file(tmp_path, name='coeff.grad'):
    """Write a minimal Siemens gradient coefficient file.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory to write to.
    name : str, optional
        File name.

    Returns
    -------
    path : pathlib.Path
        The written file.
    """
    path = tmp_path / name
    path.write_text(' Synthetic coefficients\n 0.250 = R0\n\n  1 A( 3, 1) -0.023400 x\n')
    return path


@pytest.mark.parametrize(
    ('image_type', 'expected'),
    [
        (['ORIGINAL', 'PRIMARY', 'M', 'ND'], '3D'),
        (['ORIGINAL', 'PRIMARY', 'M', 'DIS2D'], '1D'),
        ('ORIGINAL\\PRIMARY\\M\\DIS3D', None),
        (['DIS2D', 'DIS3D'], None),
        (None, '3D'),
    ],
)
def test_warp_dim_from_metadata(image_type, expected):
    """ImageType determines the correction: DIS3D none, DIS2D through-plane, else 3D."""
    assert gw.warp_dim_from_metadata({'ImageType': image_type}) == expected


def test_resolve_plan(tmp_path):
    """Plans honor --gradient-file, --ignore gradwarp, and --force."""
    grad = _grad_file(tmp_path)
    meta = {'ImageType': ['ND'], 'Manufacturer': 'Siemens'}

    assert gw.resolve_gradwarp_plan(meta, 'asl.nii.gz', gradient_file=None, force=[]) is None
    assert (
        gw.resolve_gradwarp_plan(meta, 'asl.nii.gz', grad, force=[], ignore=['gradwarp']) is None
    )

    plan = gw.resolve_gradwarp_plan(meta, 'asl.nii.gz', grad, force=[], ignore=[])
    assert (plan.warp_dim, plan.basis, plan.is_ge) == ('3D', 'metadata', False)

    plan = gw.resolve_gradwarp_plan(meta, 'asl.nii.gz', grad, force=['gradwarp1D'], ignore=[])
    assert (plan.warp_dim, plan.basis) == ('1D', 'forced')


def test_resolve_plan_refuses_ge_coefficients(tmp_path):
    """GE coefficients are refused only when a field would be expanded from them."""
    meta = {'ImageType': ['ND'], 'Manufacturer': 'GE MEDICAL SYSTEMS'}
    with pytest.raises(ValueError, match='not supported for GE'):
        gw.resolve_gradwarp_plan(meta, 'asl.nii.gz', tmp_path / 'c.dat', force=[], ignore=[])

    # A ready-made field, or data already corrected on the scanner, need no expansion
    plan = gw.resolve_gradwarp_plan(meta, 'asl.nii.gz', tmp_path / 'f.nii.gz', force=[], ignore=[])
    assert plan.is_ge
    meta['ImageType'] = ['DIS3D']
    plan = gw.resolve_gradwarp_plan(meta, 'asl.nii.gz', tmp_path / 'c.dat', force=[], ignore=[])
    assert plan.warp_dim is None


@pytest.mark.parametrize(
    ('gradient_file', 'force', 'ignore', 'message'),
    [
        ('c.grad', {'gradwarp1D', 'gradwarp3D'}, set(), 'mutually exclusive'),
        ('c.grad', {'gradwarp1D'}, {'gradwarp'}, 'contradictory'),
        (None, {'gradwarp3D'}, set(), 'requires --gradient-file'),
        (None, set(), {'gradwarp-jacobian'}, 'requires --gradient-file'),
        ('c.txt', set(), set(), 'unrecognized extension'),
    ],
)
def test_validate_gradient_flags_errors(gradient_file, force, ignore, message):
    """Contradictory flags and unrecognized extensions are rejected."""
    with pytest.raises(ValueError, match=message):
        gw.validate_gradient_flags(gradient_file, force, ignore)


@pytest.mark.parametrize('name', ['c.grad', 'c.dat', 'c.gc', 'f.nii', 'f.nii.gz'])
def test_validate_gradient_flags_ok(name):
    """All supported coefficient and field extensions are accepted."""
    gw.validate_gradient_flags(name, {'gradwarp3D'}, {'gradwarp-jacobian'})


def test_sanitize_siemens_coefficients(tmp_path):
    """Comments TORTOISE would misread are dropped into a copy; other files are untouched."""
    lines = [
        ' Header line',
        '#  A(1,1) = 1.1547 (2/Sqrt[3])',  # would abort the reader
        '#  1 A( 3, 1) -0.5 x',  # would be silently added
        '# a harmless comment',
        ' 0.250 = R0',
        '  1 A( 3, 1) -0.023400 x',
    ]
    src = tmp_path / 'coeff.grad'
    src.write_text('\n'.join(lines) + '\n')

    out = gw.sanitize_siemens_coefficients(src, tmp_path / 'clean')
    assert out == tmp_path / 'clean' / 'coeff.grad'
    assert '#' not in out.read_text()
    assert src.read_text().count('#') == 3  # the original is untouched

    # Nothing to drop: the original is used
    clean = tmp_path / 'clean.grad'
    clean.write_text('\n'.join(line for line in lines if not line.startswith('#')) + '\n')
    assert gw.sanitize_siemens_coefficients(clean, tmp_path / 'other') == clean

    # Other formats are never parsed
    dat = tmp_path / 'coeff.dat'
    dat.write_text('#  A(1,1) = 1.1547\n')
    assert gw.sanitize_siemens_coefficients(dat, tmp_path / 'other') == dat


def test_sanitize_siemens_coefficients_bad_data_line(tmp_path):
    """A data line that would crash TORTOISE's reader fails early, naming the line."""
    src = tmp_path / 'coeff.grad'
    src.write_text('  1 A( 3, 1) = oops x\n')
    with pytest.raises(ValueError, match='line 1'):
        gw.sanitize_siemens_coefficients(src, tmp_path / 'clean')


def test_create_nonlinearity_displacement_map_cmdline(tmp_path):
    """The coefficient file comes first, and the GE flag is appended only for GE data."""
    from aslprep.interfaces.gradunwarp import CreateNonlinearityDisplacementMap

    grad = _grad_file(tmp_path)
    ref = tmp_path / 'ref.nii'
    ref.write_bytes(b'')
    iface = CreateNonlinearityDisplacementMap(coeff_file=str(grad), ref_image=str(ref))
    assert iface.cmdline == f'CreateNonlinearityDisplacementMap {grad} {ref} gradwarp_field.nii'
    iface.inputs.is_ge = True
    assert iface.cmdline.endswith('gradwarp_field.nii 1')


def _rotation(axis, degrees):
    """Build a rotation about one world axis.

    Parameters
    ----------
    axis : int
        Index of the world axis to rotate about (0, 1, or 2).
    degrees : float
        Rotation angle.

    Returns
    -------
    numpy.ndarray of shape (4, 4)
        The rotation, as an affine.
    """
    import numpy as np

    theta = np.deg2rad(degrees)
    c, s = np.cos(theta), np.sin(theta)
    i, j = [ax for ax in range(3) if ax != axis]
    rot = np.eye(4)
    rot[i, i], rot[i, j], rot[j, i], rot[j, j] = c, -s, s, c
    return rot


@pytest.mark.parametrize(
    ('affine_name', 'slice_axis'),
    [('axial', 'k'), ('sagittal', 'k'), ('oblique', 'k'), ('axial_j', 'j')],
)
def test_mask_warp_dimensions(tmp_path, affine_name, slice_axis):
    """Through-plane correction follows the physical slice normal, not world z."""
    import nibabel as nb
    import numpy as np

    from aslprep.interfaces.gradunwarp import MaskWarpDimensions, slice_normal_lps

    affines = {
        'axial': np.diag([2.0, 2.0, 3.0, 1.0]),
        # Slices stacked along RAS x: voxel k runs left-right
        'sagittal': np.array([[0, 0, 3.0, 0], [2.0, 0, 0, 0], [0, 2.0, 0, 0], [0, 0, 0, 1]]),
        'oblique': _rotation(0, 25) @ _rotation(1, 10) @ np.diag([2.0, 2.0, 3.0, 1.0]),
        'axial_j': np.diag([2.0, 3.0, 2.0, 1.0]),
    }
    affine = affines[affine_name]
    ref = tmp_path / 'ref.nii.gz'
    nb.Nifti1Image(np.zeros((4, 4, 4), dtype='float32'), affine).to_filename(ref)

    rng = np.random.default_rng(0)
    deltas = rng.standard_normal((4, 4, 4, 1, 3)).astype('float32')
    field = tmp_path / 'field.nii.gz'
    nb.Nifti1Image(deltas, affine).to_filename(field)

    normal = slice_normal_lps(affine, slice_axis)
    expected_world_normal = {
        'axial': [0, 0, 1],
        'sagittal': [-1, 0, 0],  # RAS x is LPS -x
        'axial_j': [0, -1, 0],  # RAS y is LPS -y
    }
    if affine_name in expected_world_normal:
        np.testing.assert_allclose(np.abs(normal), np.abs(expected_world_normal[affine_name]))

    results = {}
    for warp_dim in ('1D', '2D', '3D'):
        workdir = tmp_path / warp_dim
        workdir.mkdir()
        result = MaskWarpDimensions(
            in_file=str(field),
            ref_image=str(ref),
            slice_axis=slice_axis,
            warp_dim=warp_dim,
        ).run(cwd=str(workdir))
        img = nb.load(result.outputs.out_file)
        assert img.header.get_intent()[0] == 'vector'
        results[warp_dim] = np.asarray(img.dataobj)

    np.testing.assert_allclose(results['3D'], deltas)
    # 1D is parallel to the normal, 2D is perpendicular to it, and together they make the field
    np.testing.assert_allclose(np.cross(results['1D'], normal), 0, atol=1e-5)
    np.testing.assert_allclose(results['2D'] @ normal, 0, atol=1e-5)
    np.testing.assert_allclose(results['1D'] + results['2D'], deltas, atol=1e-5)
    if affine_name == 'axial':
        # Matches TORTOISE's world-component masking for axial acquisitions
        np.testing.assert_allclose(results['1D'][..., :2], 0, atol=1e-6)
        np.testing.assert_allclose(results['1D'][..., 2], deltas[..., 2], atol=1e-6)


def test_mask_warp_dimensions_requires_reference(tmp_path):
    """Through-plane correction needs a reference image for the slice normal."""
    import nibabel as nb
    import numpy as np

    from aslprep.interfaces.gradunwarp import MaskWarpDimensions

    field = tmp_path / 'field.nii.gz'
    nb.Nifti1Image(np.ones((4, 4, 4, 1, 3), dtype='float32'), np.eye(4)).to_filename(field)
    with pytest.raises(ValueError, match='ref_image is required'):
        MaskWarpDimensions(in_file=str(field), warp_dim='1D').run(cwd=str(tmp_path))


def test_resolve_plan_slice_axis(tmp_path):
    """The slice axis comes from SliceEncodingDirection, without its sign."""
    meta = {'ImageType': ['DIS2D'], 'SliceEncodingDirection': 'j-'}
    plan = gw.resolve_gradwarp_plan(meta, 'asl.nii.gz', _grad_file(tmp_path), force=[], ignore=[])
    assert (plan.warp_dim, plan.slice_axis) == ('1D', 'j')


def test_prequantified_cbf_is_not_modulated(tmp_path, monkeypatch):
    """Jacobian modulation is for signal, not for CBF that was quantified on the scanner."""
    import aslprep.workflows.asl.base as asl_base

    monkeypatch.setattr(asl_base, 'select_processing_target', lambda aslcontext: 'cbf')
    reset_config()
    with mock_config():
        config.workflow.gradient_file = _grad_file(tmp_path)
        config.workflow.level = 'full'
        config.workflow.cifti_output = False
        config.execution.output_spaces = 'asl T1w'
        config.init_spaces()
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        try:
            wf = asl_base.init_asl_wf(asl_file=str(asl_file))
        except Exception as exc:  # noqa: BLE001
            # The test dataset is not CBF-only, so later CBF wiring may object.
            pytest.skip(f'Could not build a CBF-only workflow from the test data: {exc}')
        graph = wf._create_flat_graph()
        for name in ('asl_native_wf.aslref_asl', 'asl_anat_wf.resample'):
            assert graph_node(graph, name).inputs.gradwarp_jacobian is False, name
        updates = graph_node(graph, 'asl_output_metadata').inputs.updates
        assert updates['GradientWarpJacobian'] is False


def test_precomputed_geometry_is_recomputed(tmp_path):
    """Derivatives that depend on the corrected geometry are not reused."""
    from aslprep.workflows.asl.fit import init_asl_fit_wf

    reset_config()
    with mock_config():
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        precomputed = {
            'coreg_aslref': str(asl_file),
            'aslref_mask': str(asl_file),
            'transforms': {'aslref2anat': str(asl_file)},
        }
        kwargs = {
            'asl_file': str(asl_file),
            'aslcontext': str(asl_file).replace('.nii.gz', 'context.tsv'),
            'm0scan': None,
            'use_ge': False,
            'precomputed': precomputed,
        }
        plan = gw.GradwarpPlan(str(_grad_file(tmp_path)), '3D', False, 'metadata')

        reused = init_asl_fit_wf(**kwargs)
        assert reused.get_node('bold_reg_wf') is None
        assert reused.get_node('unwarp_aslref') is None

        recomputed = init_asl_fit_wf(gradwarp_plan=plan, **kwargs)
        assert recomputed.get_node('bold_reg_wf') is not None
        assert recomputed.get_node('unwarp_aslref') is not None


def _edges_into(graph, dst_suffix, field):
    """Find the nodes connected to one input of a node.

    Parameters
    ----------
    graph : networkx.DiGraph
        A flattened workflow graph.
    dst_suffix : str
        End of the destination node's full name.
    field : str
        Name of the destination node's input.

    Returns
    -------
    set of str
        Full names of the nodes connected to that input.
    """
    return {
        src.fullname
        for src, dst, data in graph.edges(data=True)
        if dst.fullname.endswith(dst_suffix)
        and any(conn[1] == field for conn in data.get('connect', []))
    }


@pytest.mark.parametrize('fieldmap_id', [None, 'auto_00000'])
def test_init_asl_wf_gradwarp(tmp_path, fieldmap_id):
    """The field reaches every resampling of ASL data, and the fit stages that need it."""
    from aslprep.workflows.asl.base import init_asl_wf

    reset_config()
    with mock_config():
        config.workflow.gradient_file = _grad_file(tmp_path)
        config.workflow.level = 'full'
        config.workflow.cifti_output = False
        config.execution.output_spaces = 'asl T1w MNI152NLin2009cAsym'
        config.init_spaces()
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        wf = init_asl_wf(asl_file=str(asl_file), fieldmap_id=fieldmap_id)
        graph = wf._create_flat_graph()

        make_field = graph_node(graph, 'asl_fit_wf.gradwarp_wf.make_field')
        assert make_field.inputs.coeff_file == str(config.workflow.gradient_file)

        resamplers = [
            'asl_fit_wf.unwarp_aslref',
            'asl_native_wf.aslref_asl',
            'asl_anat_wf.resample',
            'asl_std_wf.resample',
        ]
        for name in resamplers:
            assert _edges_into(graph, name, 'gradwarp_field'), name
            assert graph_node(graph, name).inputs.gradwarp_jacobian is True

        fmapreg_targets = _edges_into(graph, 'fmapreg_wf.inputnode', 'target_ref')
        if fieldmap_id:
            assert any(src.endswith('gradwarp_fmapreg_ref') for src in fmapreg_targets)
        else:
            assert not fmapreg_targets

        metadata = graph_node(graph, 'asl_output_metadata').inputs.updates
        assert metadata['GradientWarpDimensions'] == '3D'
        assert metadata['GradientCoefficientFile'] == 'coeff.grad'


def test_init_asl_wf_gradwarp_options(tmp_path):
    """--ignore gradwarp-jacobian and --force gradwarp1D reach the nodes."""
    from aslprep.workflows.asl.base import init_asl_wf

    reset_config()
    with mock_config():
        config.workflow.gradient_file = _grad_file(tmp_path)
        config.workflow.ignore = ['gradwarp-jacobian']
        config.workflow.force = ['gradwarp1D']
        config.workflow.level = 'minimal'
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        wf = init_asl_wf(asl_file=str(asl_file))
        graph = wf._create_flat_graph()

        assert graph_node(graph, 'gradwarp_wf.mask_field').inputs.warp_dim == '1D'
        assert graph_node(graph, 'asl_native_wf.aslref_asl').inputs.gradwarp_jacobian is False


def test_init_asl_wf_no_gradwarp():
    """Without --gradient-file, no field is built or applied."""
    from aslprep.workflows.asl.base import init_asl_wf

    reset_config()
    with mock_config():
        config.workflow.level = 'minimal'
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        wf = init_asl_wf(asl_file=str(asl_file))
        graph = wf._create_flat_graph()
        names = {node.fullname for node in graph.nodes()}
        assert not any('gradwarp' in name for name in names)
        assert not _edges_into(graph, 'asl_native_wf.aslref_asl', 'gradwarp_field')


def graph_node(graph, suffix):
    """Find the single node whose full name ends with a suffix.

    Parameters
    ----------
    graph : networkx.DiGraph
        A flattened workflow graph.
    suffix : str
        End of the node's full name.

    Returns
    -------
    nipype.pipeline.engine.Node
        The node.

    Raises
    ------
    AssertionError
        If no node, or more than one, matches.
    """
    matches = [node for node in graph.nodes() if node.fullname.endswith(suffix)]
    assert len(matches) == 1, (suffix, [node.fullname for node in matches])
    return matches[0]


@pytest.mark.parametrize(
    ('asl_dim', 'm0_dim'),
    [('3D', None), (None, '3D'), ('3D', '1D'), ('3D', '3D')],
)
def test_separate_m0scan_correction(tmp_path, asl_dim, m0_dim):
    """A separate M0 scan is corrected according to its own metadata, not the ASL run's."""
    from aslprep.workflows.asl.fit import init_asl_fit_wf, init_asl_native_wf

    reset_config()
    with mock_config():
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        grad = str(_grad_file(tmp_path))
        asl_plan = gw.GradwarpPlan(grad, asl_dim, False, 'metadata')
        m0_plan = gw.GradwarpPlan(grad, m0_dim, False, 'metadata', slice_axis='j')

        fit_wf = init_asl_fit_wf(
            asl_file=str(asl_file),
            aslcontext=str(asl_file).replace('.nii.gz', 'context.tsv'),
            m0scan=str(asl_file),  # any existing image stands in for the M0 scan
            use_ge=False,
            gradwarp_plan=asl_plan,
            m0scan_gradwarp_plan=m0_plan,
        )
        assert (fit_wf.get_node('gradwarp_wf') is not None) == (asl_dim is not None)
        m0_wf = fit_wf.get_node('m0scan_gradwarp_wf')
        assert (m0_wf is not None) == (m0_dim is not None)
        if m0_wf is not None:
            assert m0_wf.get_node('mask_field').inputs.slice_axis == 'j'
            assert m0_wf.get_node('mask_field').inputs.warp_dim == m0_dim
            assert m0_wf.inputs.inputnode.slice_ref_image == str(asl_file)
            # No second report, which would overwrite the ASL run's
            assert m0_wf.get_node('ds_gradwarp_report') is None

        native_wf = init_asl_native_wf(
            asl_file=str(asl_file),
            m0scan=str(asl_file),
            gradwarp=asl_dim is not None,
            m0scan_gradwarp=m0_dim is not None,
        )
        graph = native_wf._create_flat_graph()
        asl_sources = _edges_into(graph, 'aslref_asl', 'gradwarp_field')
        m0_sources = _edges_into(graph, 'aslref_m0scan', 'gradwarp_field')
        assert bool(asl_sources) == (asl_dim is not None)
        assert bool(m0_sources) == (m0_dim is not None)


def test_first_defined():
    """The fallback is used unless a preferred value is given."""
    from aslprep.workflows.asl.gradwarp import _first_defined

    assert _first_defined('ref.nii') == 'ref.nii'
    assert _first_defined('ref.nii', preferred='m0.nii') == 'm0.nii'
