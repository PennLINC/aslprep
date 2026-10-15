"""Tests for ASL workflow wiring."""

import pytest
from nipype.interfaces.base import isdefined

from aslprep import config
from aslprep.tests.tests import mock_config, reset_config


def test_remove_keys():
    """Patched metadata keys are removed from preprocessed ASL sidecars."""
    from aslprep.workflows.asl.outputs import _remove_keys

    metadata = {'RepetitionTime': 4.0, 'PostLabelingDelay': [1.8, 1.8]}
    assert _remove_keys(metadata, ['RepetitionTime']) == {'PostLabelingDelay': [1.8, 1.8]}
    assert _remove_keys(metadata, []) == metadata


def test_init_asl_wf_boilerplate_order():
    """The run's boilerplate opens with the acquisition and gives each section a heading."""
    reset_config()
    with mock_config():
        config.workflow.level = 'minimal'
        config.execution.output_spaces = 'asl'
        config.init_spaces()
        from aslprep.workflows.asl.base import init_asl_wf

        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        wf = init_asl_wf(asl_file=str(asl_file))
        desc = wf.visit_desc()

    assert desc.startswith('\nASL acquisition\n\n: The ASL data were acquired with')
    for heading in (
        'Cerebral blood flow computation and denoising',
        'Confounds',
        'CBF quality control',
        'Resampling',
    ):
        assert f'\n\n{heading}\n\n: ' in desc
    assert desc.index('ASL acquisition') < desc.index('Cerebral blood flow computation')
    assert '####' not in desc


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


def test_init_asl_wf_clipping_qc_inputs():
    """The acquired ASL file reaches the QC node, which checks it for clipped values."""
    reset_config()
    with mock_config():
        from aslprep.workflows.asl.base import init_asl_wf

        asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
        wf = init_asl_wf(asl_file=str(asl_file))
        graph = wf._create_flat_graph()
        qc_inputs = {
            dst_field: (src.fullname, src_field)
            for src, dst, data in graph.edges(data=True)
            if dst.fullname.endswith('cbf_confounds_wf.compute_qc_metrics')
            for src_field, dst_field in data['connect']
        }
        assert qc_inputs['asl_file'][1] == 'asl_file'
        assert qc_inputs['m0scan'][1] == 'm0scan'
        confounds_inputnode = wf.get_node('cbf_confounds_wf.inputnode')
        # The bundled dataset has its M0 volumes in the ASL file, so there is no separate scan.
        assert not isdefined(confounds_inputnode.inputs.m0scan)


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


@pytest.mark.parametrize(
    ('volume_types', 'force', 'expect_relative'),
    [
        # GE product output (#728): the reduced series is the M0 and delta-M, two volumes.
        (['m0scan', 'deltam', 'cbf'], [], False),
        (['m0scan', 'deltam', 'cbf'], ['no-ge'], False),
        (['m0scan'] * 2 + ['deltam'] * 4 + ['cbf'], ['no-ge'], True),
    ],
)
def test_init_asl_wf_relative_confounds_use_reduced_series(
    tmp_path, volume_types, force, expect_relative
):
    """DVARS, FD and the carpet plot depend on the reduced series, not the whole ASL file."""
    import shutil

    import nibabel as nb

    reset_config()
    with mock_config():
        bids_dir = tmp_path / 'ds000240'
        shutil.copytree(config.execution.bids_dir, bids_dir)
        perf_dir = bids_dir / 'sub-01' / 'perf'
        asl_file = perf_dir / 'sub-01_asl.nii.gz'
        img = nb.load(asl_file)
        img.slicer[..., : len(volume_types)].to_filename(asl_file)
        (perf_dir / 'sub-01_aslcontext.tsv').write_text(
            '\n'.join(['volume_type', *volume_types]) + '\n'
        )

        config.execution.bids_dir = bids_dir
        config.execution.bids_database_dir = None
        config.execution.layout = None
        config.execution._layout = None
        config.execution.init()
        config.workflow.force = force
        from aslprep.workflows.asl.base import init_asl_wf

        wf = init_asl_wf(asl_file=str(asl_file))

        confounds_wf = wf.get_node('asl_confounds_wf')
        assert (confounds_wf.get_node('dvars') is not None) is expect_relative
        assert (confounds_wf.get_node('fdisp') is not None) is expect_relative
        assert (wf.get_node('carpetplot_wf') is not None) is expect_relative


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


@pytest.mark.parametrize('disable_n4', [False, True])
def test_init_synthstrip_aslref_wf_nodes(disable_n4):
    """The ASL reference workflow has no 3dUnifize branch, whose output nothing used."""
    from aslprep.workflows.asl.reference import init_synthstrip_aslref_wf

    reset_config()
    with mock_config():
        wf = init_synthstrip_aslref_wf(disable_n4=disable_n4)
        graph = wf._create_flat_graph()

    names = {node.name for node in graph.nodes()}
    assert not names & {'unifize', 'fixhdr_unifize', 'apply_mask'}
    assert ('n4_correct' in names) is not disable_n4
    outputs = {
        dst_field
        for _, dst, data in graph.edges(data=True)
        if dst.fullname == 'synthstrip_aslref_wf.outputnode'
        for _, dst_field in data['connect']
    }
    assert outputs == {'bias_corrected_file', 'mask_file'}


@pytest.mark.parametrize(
    ('sloppy', 'delays'),
    [(False, 1.8), (True, 1.8), (False, [1.0, 1.0, 1.5, 1.5, 2.0, 2.0])],
)
def test_init_cbf_wf_basil_sd(tmp_path, sloppy, delays):
    """BASIL's posterior standard deviations follow the estimates they belong to (gh-671)."""
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
    n_pairs = len(delays) if isinstance(delays, list) else 3
    reset_config()
    with mock_config():
        config.execution.sloppy = sloppy
        wf = init_cbf_wf(
            name_source=str(asl_file),
            processing_target='control',
            aslcontext=pd.DataFrame({'volume_type': ['control', 'label'] * n_pairs}),
            metadata=metadata,
            basil=True,
        )
        outputs = {
            dst_field: (src.name, src_field)
            for src, dst, data in wf._graph.edges(data=True)
            if dst.name == 'outputnode'
            for src_field, dst_field in data['connect']
        }

    assert outputs['cbf_sd_basil'] == ('basilcbf', 'cbf_sd_basil')
    for tissue in ('gm', 'wm'):
        # In sloppy mode, the corrected maps and their SDs are stand-ins for the uncorrected ones.
        source = 'cbf_sd_basil' if sloppy else f'cbf_{tissue}_sd_basil'
        assert outputs[f'cbf_{tissue}_sd_basil'] == ('basilcbf', source)
    if isinstance(delays, list):
        assert outputs['att_sd_basil'] == ('basilcbf', 'att_sd_basil')
    else:
        assert 'att_sd_basil' not in outputs


def test_basil_sd_output_metadata():
    """The posterior SD derivatives have names that do not collide, and units like the estimates."""
    from aslprep.workflows.asl.outputs import BASE_INPUT_FIELDS

    names = {
        field: (spec['desc'], spec['suffix'])
        for field, spec in BASE_INPUT_FIELDS.items()
        if 'desc' in spec
    }
    assert len(set(names.values())) == len(names)
    for field, units in [
        ('cbf_sd_basil', 'mL/100 g/min'),
        ('cbf_gm_sd_basil', 'mL/100 g/min'),
        ('cbf_wm_sd_basil', 'mL/100 g/min'),
        ('att_sd_basil', 's'),
    ]:
        assert BASE_INPUT_FIELDS[field]['Units'] == units
        estimate = field.replace('_sd_basil', '_basil').replace('cbf_', 'mean_cbf_')
        estimate = 'att_basil' if field == 'att_sd_basil' else estimate
        assert BASE_INPUT_FIELDS[estimate]['Units'] == units


@pytest.mark.parametrize('haematocrit', [None, 0.20])
def test_init_cbf_wf_haematocrit(tmp_path, haematocrit):
    """BASIL gets the same T1blood as ASLPrep's model, and the boilerplate reports a haematocrit."""
    import pandas as pd

    from aslprep.utils.cbf import estimate_t1
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
    if haematocrit is not None:
        metadata['Haematocrit'] = haematocrit
    reset_config()
    with mock_config():
        wf = init_cbf_wf(
            name_source=str(asl_file),
            processing_target='control',
            aslcontext=pd.DataFrame({'volume_type': ['control', 'label'] * 3}),
            metadata=metadata,
            basil=True,
        )

    t1blood = wf.get_node('basilcbf').inputs.t1blood
    assert t1blood == pytest.approx(estimate_t1(metadata)[0])
    if haematocrit is None:
        assert t1blood == pytest.approx(1.65)
        assert 'haematocrit' not in wf.__desc__
    else:
        assert t1blood == pytest.approx(1.948, abs=0.0005)
        assert 'measured haematocrit (0.2)' in wf.__desc__
        assert '@hales2016general' in wf.__desc__
        assert 'the value for 3 T [' not in wf.__desc__


@pytest.mark.parametrize('asl_type', ['PCASL', 'PASL'])
def test_init_cbf_wf_boilerplate(tmp_path, asl_type):
    """The boilerplate reports the acquisition and the constants used to calculate CBF."""
    import pandas as pd

    from aslprep.workflows.asl.cbf import init_cbf_wf

    asl_file = tmp_path / 'sub-01_asl.nii.gz'
    asl_file.touch()
    metadata = {
        'ArterialSpinLabelingType': asl_type,
        'M0Type': 'Estimate',
        'M0Estimate': 1000.0,
        'MagneticFieldStrength': 3,
        'BackgroundSuppression': True,
        'BackgroundSuppressionNumberPulses': 4,
        'RepetitionTimePreparation': 4.0,
        'SliceTiming': [0.0, 0.05, 0.1],
    }
    if asl_type == 'PCASL':
        metadata['LabelingDuration'] = 1.8
        metadata['PostLabelingDelay'] = [1.0, 1.0, 2.0, 2.0]
    else:
        metadata['PostLabelingDelay'] = 1.8
        metadata['BolusCutOffFlag'] = True
        metadata['BolusCutOffTechnique'] = 'QUIPSSII'
        metadata['BolusCutOffDelayTime'] = 0.7

    reset_config()
    with mock_config():
        wf = init_cbf_wf(
            name_source=str(asl_file),
            processing_target='control',
            aslcontext=pd.DataFrame({'volume_type': ['control', 'label'] * 2}),
            metadata=metadata,
        )
    desc = wf.__desc__

    # The acquisition is described by the run's workflow, before the preprocessing.
    assert 'ASL acquisition' not in desc
    # The next workflow's description starts on a new line.
    assert desc.endswith('\n')
    assert 'The T1 of arterial blood was set to 1.65 s, the value for 3 T' in desc
    assert 'blood-brain partition coefficient of 0.9 mL/g' in desc
    assert 'each of the 4 background suppression pulses' in desc
    if asl_type == 'PCASL':
        assert 'the post-labeling delays were shifted' in desc
        assert 'M0a and M0b' not in desc
        assert 'arterial blood (M0b) was estimated voxel-wise' in desc
    else:
        assert 'the inversion times were shifted' in desc
        assert 'post-labeling delay' not in desc
        assert 'QUIPSS II bolus cut-off technique [@wong1998quantitative]' in desc
        assert 'as described in @alsop_recommended_2015' in desc
