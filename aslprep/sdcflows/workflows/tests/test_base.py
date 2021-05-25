"""Test base workflows."""
import pytest

from ..base import init_sdc_estimate_wf

EPI_METADATA = {
    'MultibandAccelerationFactor': 8,
    'PhaseEncodingDirection': 'i',
    'RepetitionTime': 0.72,
    'TaskName': 'Resting-state fMRI',
}
EPI_FMAP_METADATA_1 = {
    'BandwidthPerPixelPhaseEncode': 2290,
    'EPIFactor': 90,
    'EffectiveEchoSpacing': 0.00058,
    'IntendedFor': ['func/sub-HCP101006_task-rest_dir-LR_bold.nii.gz',
                    'func/sub-HCP101006_task-rest_dir-LR_sbref.nii.gz'],
    'MultibandAccelerationFactor': 1,
    'PhaseEncodingDirection': 'i-',
}
EPI_FMAP_METADATA_2 = EPI_FMAP_METADATA_1.copy()
EPI_FMAP_METADATA_2['PhaseEncodingDirection'] = 'i'

PHDIFF_METADATA = {
    'EchoTime1': 0.00492,
    'EchoTime2': 0.00738,
}
PHASE1_METADATA = {
    'EchoTime': 0.00492,
}
PHASE2_METADATA = {
    'EchoTime': 0.00738,
}


FMAP_DICT_ELEMENTS = {
    'epi1': [('sub-HCP101006/fmap/sub-HCP101006_dir-RL_epi.nii.gz',
              EPI_FMAP_METADATA_1.copy())],
    'epi2': [('sub-HCP101006/fmap/sub-HCP101006_dir-RL_epi.nii.gz',
              EPI_FMAP_METADATA_1.copy()),
             ('sub-HCP101006/fmap/sub-HCP101006_dir-LR_epi.nii.gz',
              EPI_FMAP_METADATA_2.copy())],
    'phdiff1': {
        'magnitude': [('sub-HCP101006/fmap/sub-HCP101006_magnitude1.nii.gz', {}),
                      ('sub-HCP101006/fmap/sub-HCP101006_magnitude2.nii.gz', {})],
        'phases': [('sub-HCP101006/fmap/sub-HCP101006_phasediff.nii.gz', PHDIFF_METADATA)],
    },
    'phdiff2': {
        'magnitude': [('sub-HCP101006/fmap/sub-HCP101006_magnitude1.nii.gz', {}),
                      ('sub-HCP101006/fmap/sub-HCP101006_magnitude2.nii.gz', {})],
        'phases': [('sub-HCP101006/fmap/sub-HCP101006_phase1.nii.gz',
                    PHASE1_METADATA.copy()),
                   ('sub-HCP101006/fmap/sub-HCP101006_phase2.nii.gz',
                    PHASE2_METADATA.copy())],
    },
    'fmap1': {
        'magnitude': [('sub-HCP101006/fmap/sub-HCP101006_magnitude.nii.gz', {})],
        'fieldmap': [('sub-HCP101006/fmap/sub-HCP101006_fieldmap.nii.gz', {})],
    },
    'syn': {
        't1w': [('sub-HCP101006/fmap/sub-HCP101006_T1w.nii.gz', {})]
    }
}


@pytest.mark.parametrize('method', ['skip', 'phasediff', 'pepolar', 'fieldmap', 'syn'])
def test_base(method):
    """Check the heuristics are correctly applied."""
    fieldmaps = {
        'epi': FMAP_DICT_ELEMENTS['epi1'].copy(),
        'fieldmap': [FMAP_DICT_ELEMENTS['fmap1'].copy()],
        'phasediff': [FMAP_DICT_ELEMENTS['phdiff1'].copy()],
    }
    epi_meta = EPI_METADATA.copy()

    if method == 'skip':
        wf = init_sdc_estimate_wf(fmaps=None, epi_meta=epi_meta)
        assert wf.inputs.outputnode.method == 'None'

        with pytest.raises(ValueError):
            wf = init_sdc_estimate_wf(fmaps={'unsupported': None}, epi_meta=epi_meta)
    elif method == 'pepolar':
        wf = init_sdc_estimate_wf(fmaps=fieldmaps, epi_meta=epi_meta)
        assert 'PEPOLAR' in wf.inputs.outputnode.method

        fieldmaps['epi'][0][1].pop('PhaseEncodingDirection')
        with pytest.raises(ValueError):
            wf = init_sdc_estimate_wf(fmaps=fieldmaps, epi_meta=epi_meta)

        fieldmaps['epi'][0][1]['PhaseEncodingDirection'] = \
            EPI_FMAP_METADATA_1['PhaseEncodingDirection']
        epi_meta.pop('PhaseEncodingDirection')
        with pytest.raises(ValueError):
            wf = init_sdc_estimate_wf(fmaps=fieldmaps, epi_meta=epi_meta)
        epi_meta = EPI_METADATA.copy()

        fieldmaps['epi'] = FMAP_DICT_ELEMENTS['epi2']
        wf = init_sdc_estimate_wf(fmaps=fieldmaps, epi_meta=epi_meta)
        assert 'PEPOLAR' in wf.inputs.outputnode.method

    elif method == 'fieldmap':
        fieldmaps = {
            'fieldmap': [FMAP_DICT_ELEMENTS['fmap1'].copy()],
            'phasediff': [FMAP_DICT_ELEMENTS['phdiff1'].copy()],
        }
        wf = init_sdc_estimate_wf(fmaps=fieldmaps, epi_meta=epi_meta)
        assert 'directly measured B0 map' in wf.inputs.outputnode.method

    elif method == 'phasediff':
        fieldmaps = {
            'phasediff': [FMAP_DICT_ELEMENTS['phdiff1'].copy()],
        }

        wf = init_sdc_estimate_wf(fmaps=fieldmaps, epi_meta=epi_meta)
        assert 'phase-difference' in wf.inputs.outputnode.method

        epi_meta.pop('PhaseEncodingDirection')
        with pytest.raises(ValueError):
            wf = init_sdc_estimate_wf(fmaps=fieldmaps, epi_meta=epi_meta)

    elif method == 'syn':
        fmaps_onlysyn = {'syn': FMAP_DICT_ELEMENTS['syn']}
        wf = init_sdc_estimate_wf(fmaps=fmaps_onlysyn, epi_meta=epi_meta)
        assert 'SyN' in wf.inputs.outputnode.method

        epi_pe = epi_meta.pop('PhaseEncodingDirection')
        wf = init_sdc_estimate_wf(fmaps=fmaps_onlysyn, epi_meta=epi_meta)
        assert 'SyN' in wf.inputs.outputnode.method

        # Emulate --force-syn
        fieldmaps.update(fmaps_onlysyn)
        with pytest.raises(ValueError):
            wf = init_sdc_estimate_wf(fmaps=fieldmaps, epi_meta=epi_meta)

        epi_meta['PhaseEncodingDirection'] = epi_pe
        wf = init_sdc_estimate_wf(fmaps=fieldmaps, epi_meta=epi_meta)
        assert 'PEPOLAR' in wf.inputs.outputnode.method
