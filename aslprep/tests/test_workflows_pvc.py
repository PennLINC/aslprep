"""Workflow-construction tests for local-regression partial volume correction."""

import nibabel as nb
import numpy as np
import pandas as pd
import pytest

from aslprep import config
from aslprep.tests.tests import mock_config, reset_config

SINGLE_DELAY = {
    'ArterialSpinLabelingType': 'PCASL',
    'PostLabelingDelay': 1.8,
    'LabelingDuration': 1.8,
    'M0Type': 'Estimate',
    'M0Estimate': 1000.0,
    'MagneticFieldStrength': 3,
    'BackgroundSuppression': False,
    'RepetitionTimePreparation': 4.0,
    'MRAcquisitionType': '2D',
    'SliceTiming': [0.0, 0.05, 0.1],
}


def _edges(wf, dst_name=None, src_name=None):
    return [
        (src.name, src_field, dst.name, dst_field)
        for src, dst, data in wf._graph.edges(data=True)
        for src_field, dst_field in data['connect']
        if (dst_name is None or dst.name == dst_name)
        and (src_name is None or src.name == src_name)
    ]


def _build_cbf_wf(tmp_path, **kwargs):
    from aslprep.workflows.asl.cbf import init_cbf_wf

    asl_file = tmp_path / 'sub-01_asl.nii.gz'
    asl_file.touch()
    defaults = {
        'name_source': str(asl_file),
        'processing_target': 'control',
        'aslcontext': pd.DataFrame({'volume_type': ['control', 'label'] * 3}),
        'metadata': SINGLE_DELAY,
        'pvc_method': 'local-regression',
    }
    defaults.update(kwargs)
    return init_cbf_wf(**defaults)


@pytest.mark.parametrize('scorescrub', [False, True])
def test_init_cbf_wf_lr_pvc(tmp_path, scorescrub):
    """LR PVC corrects mean CBF (and SCORE), from ComputeCBF's calibration outputs."""
    reset_config()
    with mock_config():
        wf = _build_cbf_wf(tmp_path, scorescrub=scorescrub)
        to_pvc = _edges(wf, dst_name='lr_pvc_wf')
        to_output = {
            dst_field: src for src, _, dst, dst_field in _edges(wf, dst_name='outputnode')
        }
        assert wf.get_node('asl_pv_wf') is not None
        lr_wf = wf.get_node('lr_pvc_wf')
        assert lr_wf.get_node('lr_pvc') is not None
        assert (lr_wf.get_node('lr_pvc_score') is not None) is scorescrub
        assert ('compute_cbf', 'm0_calib', 'lr_pvc_wf', 'inputnode.m0_calib') in to_pvc
        assert ('compute_cbf', 'quant_valid', 'lr_pvc_wf', 'inputnode.quant_valid') in to_pvc
        assert to_output['mean_cbf_gm_lr'] == 'lr_pvc_wf'
        assert to_output['gm_pv'] == 'asl_pv_wf'
        assert to_output['csf_pv'] == 'asl_pv_wf'
        assert to_output['mean_cbf_csf_lr'] == 'lr_pvc_wf'
        assert to_output['pvc_qc'] == 'lr_pvc_wf'
        assert ('mean_cbf_gm_lr_score' in to_output) is scorescrub
        assert 'mean_cbf_scrub_lr' not in str(to_output)
        desc = wf.visit_desc()
        assert 'local linear regression' in desc
        assert 'effective resolution' in desc


def test_init_cbf_wf_no_pvc(tmp_path):
    """Without --pvc-method, no PVC nodes or outputs are added."""
    reset_config()
    with mock_config():
        wf = _build_cbf_wf(tmp_path, pvc_method='none')
        assert wf.get_node('asl_pv_wf') is None
        assert wf.get_node('lr_pvc_wf') is None
        assert 'mean_cbf_gm_lr' not in wf.get_node('outputnode').outputs.copyable_trait_names()


@pytest.mark.parametrize(
    ('kwargs', 'message'),
    [
        (
            {'metadata': {**SINGLE_DELAY, 'PostLabelingDelay': [1.0, 1.5, 2.0]}},
            'single-delay',
        ),
        ({'processing_target': 'cbf'}, 'scanner-provided CBF'),
        ({'pvc_kernel_size': (5, 5, 3)}, 'single slice'),
    ],
    ids=['multi-delay', 'cbf-only', 'through-plane-slice-timing'],
)
def test_init_cbf_wf_lr_pvc_unsupported(tmp_path, kwargs, message):
    """LR PVC refuses data it cannot model (spec Sections 5.1 and 5.5)."""
    reset_config()
    with mock_config(), pytest.raises(ValueError, match=message):
        _build_cbf_wf(tmp_path, **kwargs)


def test_init_cbf_wf_lr_pvc_through_plane_3d(tmp_path):
    """Without slice timing, a through-plane kernel is allowed."""
    reset_config()
    with mock_config():
        metadata = {k: v for k, v in SINGLE_DELAY.items() if k != 'SliceTiming'}
        metadata['MRAcquisitionType'] = '3D'
        wf = _build_cbf_wf(tmp_path, metadata=metadata, pvc_kernel_size=(5, 5, 3))
        lr_node = wf.get_node('lr_pvc_wf').get_node('lr_pvc')
        assert lr_node.inputs.kernel_size == [5, 5, 3]


def test_parcellate_cbf_wf_uses_validity_masks():
    """LR maps are parcellated within their own validity masks, with coverage files."""
    from aslprep.workflows.asl.cbf import init_parcellate_cbf_wf

    reset_config()
    with mock_config():
        wf = init_parcellate_cbf_wf(
            cbf_3d=['mean_cbf'],
            sparse_3d={'mean_cbf_gm_lr': 'gm_valid_lr'},
        )
        masks = {
            dst: src_field
            for _, src_field, dst, dst_field in _edges(wf, src_name='inputnode')
            if dst_field == 'mask'
        }
        assert masks['parcellate_mean_cbf'] == 'asl_mask'
        assert masks['parcellate_mean_cbf_gm_lr'] == 'gm_valid_lr'
        assert wf.get_node('ds_coverage_mean_cbf_gm_lr') is not None
        assert wf.get_node('ds_mean_cbf_gm_lr').inputs.desc == 'lrGM'


def _init_asl_wf(pvc_method, level='full'):
    from aslprep.workflows.asl.base import init_asl_wf

    config.workflow.level = level
    config.workflow.cifti_output = False
    config.workflow.pvc_method = pvc_method
    config.execution.output_spaces = 'asl T1w MNI152NLin2009cAsym'
    config.init_spaces()
    asl_file = config.execution.bids_dir / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
    return init_asl_wf(asl_file=str(asl_file))


def test_init_asl_wf_default_has_no_pvc():
    """With the default options, the workflow has no PVC nodes or sinks (spec test 12)."""
    reset_config()
    with mock_config():
        wf = _init_asl_wf('none')
        names = wf.list_node_names()
        assert not [n for n in names if 'lr_pvc' in n or 'asl_pv_wf' in n or 'pvc_native' in n]
        assert not [n for n in names if '_lr' in n]
        assert 'Resampling' in wf.__postdesc__
        assert 'normalized by the resampled mask' not in wf.__postdesc__


def test_init_asl_wf_lr_pvc():
    """LR PVC outputs are written natively, resampled with normalization, and reported."""
    reset_config()
    with mock_config():
        wf = _init_asl_wf('local-regression')
        names = wf.list_node_names()
        assert 'cbf_wf.lr_pvc_wf.lr_pvc' in names
        assert 'cbf_wf.asl_pv_wf.pv_maps' in names
        assert 'ds_pvc_native_wf.ds_status_lr' in names
        assert 'ds_pvc_native_wf.ds_gm_pv' in names
        assert 'ds_pvc_native_wf.ds_csf_pv' in names
        assert 'ds_pvc_native_wf.ds_mean_cbf_csf_lr' in names
        # The CSF QC map is native only, and its scalar summary reaches the QC file.
        assert not [n for n in names if 'csf_lr' in n and not n.startswith('ds_pvc_native_wf')]
        qc_edges = [
            e for e in _edges(wf, dst_name='cbf_confounds_wf') if e[3] == 'inputnode.pvc_qc'
        ]
        assert qc_edges == [
            ('cbf_wf', 'outputnode.pvc_qc', 'cbf_confounds_wf', 'inputnode.pvc_qc')
        ]
        qc_node = wf.get_node('cbf_confounds_wf')
        assert ('inputnode', 'pvc_qc', 'compute_qc_metrics', 'pvc_qc') in _edges(
            qc_node, dst_name='compute_qc_metrics'
        )
        assert 'ds_asl_t1_wf.normalize_mean_cbf_gm_lr' in names
        assert 'ds_asl_std_wf.normalize_mean_cbf_wm_lr' in names
        assert 'cbf_reporting_wf.lr_diagnostics' in names
        # Standard errors and diagnostics are native only.
        assert not [n for n in names if n.startswith('ds_asl_t1_wf') and 'se_lr' in n]
        # LR maps are not scored by the CBF QC workflow.
        assert not [e for e in _edges(wf, dst_name='cbf_confounds_wf') if '_lr' in e[1]]
        assert 'normalized by the resampled mask' in wf.__postdesc__


def test_lr_pvc_derivative_filenames(tmp_path):
    """Every new output builds a valid path through ASLPrep's sink configuration."""
    from aslprep.interfaces.bids import DerivativesDataSink
    from aslprep.workflows.asl.outputs import BASE_INPUT_FIELDS, LR_PVC_SOURCES, lr_pvc_names

    names = ['gm_pv', 'wm_pv', 'csf_pv']
    for source in LR_PVC_SOURCES:
        lr_names = lr_pvc_names(source)
        names += lr_names['cbf'] + lr_names['se'] + lr_names['diagnostics'] + lr_names['csf']

    in_file = tmp_path / 'in.nii.gz'
    nb.Nifti1Image(np.zeros((2, 2, 2), np.float32), np.eye(4)).to_filename(in_file)
    source_file = tmp_path / 'sub-01' / 'perf' / 'sub-01_asl.nii.gz'
    out_names = set()
    for name in names:
        fields = BASE_INPUT_FIELDS[name]
        results = DerivativesDataSink(
            base_directory=str(tmp_path / 'out'),
            source_file=str(source_file),
            in_file=str(in_file),
            compress=True,
            **fields,
        ).run(cwd=tmp_path)
        out_file = results.outputs.out_file
        out_name = out_file.split('/')[-1] if '/' in out_file else out_file.split('\\')[-1]
        out_names.add(out_name)
        assert '/perf/' in out_file.replace('\\', '/')
        assert out_name.startswith('sub-01_')
        assert out_name.endswith(f'_{fields["suffix"]}.nii.gz')
    assert 'sub-01_label-GM_desc-psf_probseg.nii.gz' in out_names
    assert 'sub-01_desc-lrGM_cbf.nii.gz' in out_names
    assert 'sub-01_desc-scoreLRWMSE_cbf.nii.gz' in out_names
    assert 'sub-01_desc-lr_pvstatus.nii.gz' in out_names
    assert 'sub-01_desc-lrCSF_cbf.nii.gz' in out_names
    assert 'sub-01_desc-lrCSF_pvstatus.nii.gz' in out_names
    assert 'sub-01_label-CSF_desc-psf_probseg.nii.gz' in out_names
    assert len(out_names) == len(names)
