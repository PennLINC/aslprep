"""Tests for ASL workflow wiring."""

import pytest

from aslprep import config
from aslprep.tests.tests import mock_config, reset_config


def test_remove_keys():
    """Patched metadata keys are removed from preprocessed ASL sidecars."""
    from aslprep.workflows.asl.outputs import _remove_keys

    metadata = {'RepetitionTime': 4.0, 'PostLabelingDelay': [1.8, 1.8]}
    assert _remove_keys(metadata, ['RepetitionTime']) == {'PostLabelingDelay': [1.8, 1.8]}
    assert _remove_keys(metadata, []) == metadata


@pytest.mark.parametrize('level', ['minimal', 'full'])
def test_init_asl_wf_output_metadata(level):
    """Preprocessed ASL sidecars use the reduced series' metadata."""
    reset_config()
    with mock_config():
        config.workflow.level = level
        config.workflow.cifti_output = False
        config.execution.output_spaces = 'asl T1w MNI152NLin2009cAsym'
        config.init_spaces()
        from aslprep.workflows.asl.base import init_asl_wf

        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        wf = init_asl_wf(asl_file=str(asl_file))
        graph = wf._create_flat_graph()
        targets = {
            dst.fullname
            for src, dst in graph.edges()
            if src.fullname.endswith('asl_output_metadata')
        }
        assert any('ds_asl_native_wf' in name for name in targets)
        if level == 'full':
            assert any('ds_asl_t1_wf' in name for name in targets)
            assert any('ds_asl_std_wf' in name for name in targets)


def test_init_asl_fit_wf_fallback_total_readout_time():
    """The TotalReadoutTime fallback reaches the ASL distortion parameters.

    Reproduces the ASL side of nipreps/fmriprep#3423.
    """
    from aslprep.workflows.asl.fit import init_asl_fit_wf, init_asl_native_wf

    reset_config()
    with mock_config():
        config.workflow.fallback_total_readout_time = 0.05
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        fit_wf = init_asl_fit_wf(
            asl_file=str(asl_file),
            aslcontext=str(asl_file).replace('.nii.gz', 'context.tsv'),
            m0scan=None,
            use_ge=False,
            fieldmap_id='auto_00000',
        )
        assert fit_wf.get_node('distortion_params').inputs.fallback == 0.05
        native_wf = init_asl_native_wf(asl_file=str(asl_file), fieldmap_id='auto_00000')
        assert native_wf.get_node('distortion_params').inputs.fallback == 0.05


def test_init_asl_wf_t1w_outputs_on_asl_grid():
    """T1w-space outputs share the ASL-resolution resampling reference, as in fMRIPrep."""
    reset_config()
    with mock_config():
        config.workflow.level = 'full'
        config.workflow.cifti_output = False
        config.execution.output_spaces = 'asl T1w MNI152NLin2009cAsym'
        config.init_spaces()
        from aslprep.workflows.asl.base import init_asl_wf

        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        wf = init_asl_wf(asl_file=str(asl_file))
        ds_inputnode = wf.get_node('ds_asl_t1_wf.inputnode')
        edges = [
            (src, data)
            for src, dst, data in wf._graph.edges(data=True)
            if dst.name == 'ds_asl_t1_wf'
        ]
        ref_sources = [
            src.name
            for src, data in edges
            if ('inputnode.ref_file' in [c[1] for c in data['connect']])
        ]
        assert ref_sources == ['asl_anat_wf']
        assert ds_inputnode is not None


@pytest.mark.parametrize(
    ('delays', 'use_gpu', 'n_procs', 'node_use_gpu'),
    [
        ([1.0, 1.0, 1.5, 1.5, 2.0, 2.0], False, 3, False),
        ([1.0, 1.0, 1.5, 1.5, 2.0, 2.0], True, 1, True),
        (1.8, True, 1, False),
    ],
)
def test_init_cbf_wf_fit_resources(tmp_path, delays, use_gpu, n_procs, node_use_gpu):
    """The multi-delay fit gets threads on the CPU, or one process and a GPU slot on a GPU."""
    import pandas as pd

    from aslprep.workflows.asl.cbf import init_cbf_wf

    asl_file = tmp_path / 'sub-01_asl.nii.gz'
    asl_file.touch()
    metadata = {
        'ArterialSpinLabelingType': 'PCASL',
        'PostLabelingDelay': delays,
        'LabelingDuration': 1.8,
        'M0Type': 'Estimate',
        'M0Estimate': 1000.0,
        'MagneticFieldStrength': 3,
        'BackgroundSuppression': False,
        'RepetitionTimePreparation': 4.0,
    }
    with mock_config():
        config.nipype.omp_nthreads = 8
        config.nipype.use_gpu = use_gpu
        wf = init_cbf_wf(
            name_source=str(asl_file),
            processing_target='control',
            aslcontext=pd.DataFrame({'volume_type': ['control', 'label'] * 3}),
            metadata=metadata,
        )
        node = wf.get_node('compute_cbf')

    assert node.n_procs == n_procs
    assert node.inputs.num_threads == n_procs
    assert node.inputs.use_gpu is node_use_gpu
    assert node.is_gpu_node() is node_use_gpu


@pytest.mark.parametrize(
    ('ignore', 'reference_volume_type'),
    [([], 'cbf'), (['cbf'], 'm0scan')],
)
def test_init_asl_wf_ignore_volume_types(tmp_path, ignore, reference_volume_type):
    """--ignore cbf keeps scanner-calculated CBF volumes out of the reference image."""
    import shutil
    from pathlib import Path

    reset_config()
    with mock_config():
        bids_dir = tmp_path / 'ds000240'
        shutil.copytree(config.execution.bids_dir, bids_dir)
        aslcontext = bids_dir / 'sub-01' / 'perf' / 'sub-01_aslcontext.tsv'
        volume_types = aslcontext.read_text().split()
        # Stand in for a series with scanner-calculated CBF volumes at the end.
        aslcontext.write_text('\n'.join([*volume_types[:-2], 'cbf', 'cbf']) + '\n')

        config.execution.bids_dir = bids_dir
        config.execution.bids_database_dir = None
        config.execution.layout = None
        config.execution._layout = None
        config.execution.init()
        config.workflow.ignore = ignore
        from aslprep.workflows.asl.base import init_asl_wf

        asl_file = bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        wf = init_asl_wf(asl_file=str(asl_file))

        hmc_aslref_wf = wf.get_node('asl_fit_wf').get_node('hmc_aslref_wf')
        assert f'generated from the {reference_volume_type} volumes' in hmc_aslref_wf.__desc__
        used_aslcontext = wf.get_node('asl_fit_wf.inputnode').inputs.aslcontext
        assert (used_aslcontext == str(aslcontext)) is (not ignore)
        assert ('ignored' in Path(used_aslcontext).read_text()) is bool(ignore)


@pytest.mark.parametrize('sloppy', [False, True])
def test_init_cbf_wf_basil_sloppy(tmp_path, sloppy):
    """Sloppy mode skips BASIL's partial volume correction, which takes most of its run time."""
    import pandas as pd

    from aslprep.workflows.asl.cbf import init_cbf_wf

    asl_file = tmp_path / 'sub-01_asl.nii.gz'
    asl_file.touch()
    metadata = {
        'ArterialSpinLabelingType': 'PCASL',
        'PostLabelingDelay': 1.8,
        'LabelingDuration': 1.8,
        'M0Type': 'Estimate',
        'M0Estimate': 1000.0,
        'MagneticFieldStrength': 3,
        'BackgroundSuppression': False,
        'RepetitionTimePreparation': 4.0,
    }
    reset_config()
    with mock_config():
        config.execution.sloppy = sloppy
        wf = init_cbf_wf(
            name_source=str(asl_file),
            processing_target='control',
            aslcontext=pd.DataFrame({'volume_type': ['control', 'label'] * 3}),
            metadata=metadata,
            basil=True,
        )
        basilcbf = wf.get_node('basilcbf')
        outputs = {
            dst_field: (src.name, src_field)
            for src, dst, data in wf._graph.edges(data=True)
            if dst.name == 'outputnode'
            for src_field, dst_field in data['connect']
        }
        basil_inputs = {
            dst_field
            for _, dst, data in wf._graph.edges(data=True)
            if dst.name == 'basilcbf'
            for _, dst_field in data['connect']
        }

    assert basilcbf.inputs.pvc is (not sloppy)
    assert ('--pvcorr' in basilcbf.interface._parse_inputs(skip=['deltam', 'mask'])) is (
        not sloppy
    )
    assert outputs['mean_cbf_basil'] == ('basilcbf', 'mean_cbf_basil')
    if sloppy:
        # The corrected maps are stand-ins, so the rest of the workflow still runs.
        assert outputs['mean_cbf_gm_basil'] == ('basilcbf', 'mean_cbf_basil')
        assert outputs['mean_cbf_wm_basil'] == ('basilcbf', 'mean_cbf_basil')
        assert not {'gm_tpm', 'wm_tpm'} & basil_inputs
        assert 'Partial volume correction was skipped' in wf.__desc__
    else:
        assert outputs['mean_cbf_gm_basil'] == ('basilcbf', 'mean_cbf_gm_basil')
        assert outputs['mean_cbf_wm_basil'] == ('basilcbf', 'mean_cbf_wm_basil')
        assert {'gm_tpm', 'wm_tpm'} <= basil_inputs
        assert 'partial-volume corrected CBF image' in wf.__desc__
