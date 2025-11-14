"""Command-line interface tests."""

import os

import pytest
from nipype import config as nipype_config

from aslprep.cli.parser import parse_args
from aslprep.cli.workflow import build_boilerplate, build_workflow
from aslprep.data import load as load_data
from aslprep.reports.core import generate_reports
from aslprep.tests.utils import (
    check_generated_files,
    download_test_data,
    get_test_data_path,
    update_resources,
)

nipype_config.enable_debug_mode()


@pytest.mark.examples_pasl_multipld
def test_examples_pasl_multipld(data_dir, output_dir, working_dir):
    """Run aslprep on the asl_003 ASL-BIDS examples dataset.

    This dataset has 10 control-label pairs at 10 different PLDs, along with a separate M0 scan.
    The BolusCutOffTechnique is Q2TIPS.
    The manufacturer is Siemens.

    PASL multi-delay data is not yet supported.
    """
    TEST_NAME = 'examples_pasl_multipld'
    PARTICIPANT_LABEL = '01'

    dataset_dir = download_test_data(TEST_NAME, data_dir)
    download_test_data('anatomical', data_dir)
    out_dir = os.path.join(output_dir, TEST_NAME, 'aslprep')
    work_dir = os.path.join(working_dir, TEST_NAME)

    parameters = [
        dataset_dir,
        out_dir,
        'participant',
        f'--participant-label={PARTICIPANT_LABEL}',
        f'-w={work_dir}',
        '--output-spaces=asl',
        '--scorescrub',
        '--basil',
        '--m0_scale=10',
        '--fs-no-resume',
        f'--fs-subjects-dir={os.path.join(data_dir, "anatomical/freesurfer")}',
        '--derivatives',
        f'{os.path.join(data_dir, "anatomical/smriprep")}',
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.examples_pcasl_multipld
def test_examples_pcasl_multipld(data_dir, output_dir, working_dir):
    """Run aslprep on the asl_004 ASL-BIDS examples dataset.

    This dataset has 48 control-label pairs at 6 different PLDs, along with a separate M0 scan.
    The manufacturer is Siemens.
    """
    TEST_NAME = 'examples_pcasl_multipld'
    PARTICIPANT_LABEL = '01'

    dataset_dir = download_test_data(TEST_NAME, data_dir)
    download_test_data('anatomical', data_dir)
    out_dir = os.path.join(output_dir, TEST_NAME, 'aslprep')
    work_dir = os.path.join(working_dir, TEST_NAME)

    parameters = [
        dataset_dir,
        out_dir,
        'participant',
        f'--participant-label={PARTICIPANT_LABEL}',
        f'-w={work_dir}',
        '--output-spaces=asl',
        '--scorescrub',
        '--basil',
        '--m0_scale=10',
        '--fs-no-resume',
        f'--fs-subjects-dir={os.path.join(data_dir, "anatomical/freesurfer")}',
        '--derivatives',
        f'{os.path.join(data_dir, "anatomical/smriprep")}',
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.examples_pcasl_singlepld_ge
def test_examples_pcasl_singlepld_ge(data_dir, output_dir, working_dir):
    """Run aslprep on the asl_001 ASL-BIDS examples dataset.

    This test uses a GE session with two volumes: one deltam and one M0.
    """
    TEST_NAME = 'examples_pcasl_singlepld_ge'
    PARTICIPANT_LABEL = '01'

    dataset_dir = download_test_data('examples_pcasl_singlepld', data_dir)
    download_test_data('anatomical', data_dir)
    # Symlink the Freesurfer derivatives to a sub-01_ses-anat folder
    os.symlink(
        os.path.join(data_dir, 'anatomical/freesurfer/sub-01'),
        os.path.join(data_dir, 'anatomical/freesurfer/sub-01_ses-anat'),
    )
    out_dir = os.path.join(output_dir, TEST_NAME, 'aslprep')
    work_dir = os.path.join(working_dir, TEST_NAME)

    parameters = [
        dataset_dir,
        out_dir,
        'participant',
        f'--participant-label={PARTICIPANT_LABEL}',
        '--session-label',
        'ge3d',
        'anat',
        f'-w={work_dir}',
        '--output-spaces=asl',
        '--scorescrub',
        '--basil',
        '--m0_scale=96',
        '--fs-no-resume',
        f'--fs-subjects-dir={os.path.join(data_dir, "anatomical/freesurfer")}',
        '--derivatives',
        f'{os.path.join(data_dir, "anatomical/smriprep")}',
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.examples_pcasl_singlepld_philips
def test_examples_pcasl_singlepld_philips(data_dir, output_dir, working_dir):
    """Run aslprep on the asl_002 ASL-BIDS examples dataset.

    This test uses a Philips session.
    The appropriate M0 scale is unknown for this dataset, so CBF values will be inflated.
    """
    TEST_NAME = 'examples_pcasl_singlepld_philips'
    PARTICIPANT_LABEL = '01'

    dataset_dir = download_test_data('examples_pcasl_singlepld', data_dir)
    download_test_data('anatomical', data_dir)
    # Symlink the Freesurfer derivatives to a sub-01_ses-anat folder
    os.symlink(
        os.path.join(data_dir, 'anatomical/freesurfer/sub-01'),
        os.path.join(data_dir, 'anatomical/freesurfer/sub-01_ses-anat'),
    )
    out_dir = os.path.join(output_dir, TEST_NAME, 'aslprep')
    work_dir = os.path.join(working_dir, TEST_NAME)

    parameters = [
        dataset_dir,
        out_dir,
        'participant',
        f'--participant-label={PARTICIPANT_LABEL}',
        '--session-label',
        'philips2d',
        'anat',
        f'-w={work_dir}',
        '--output-spaces',
        'asl',
        'fsaverage:den-10k',
        '--scorescrub',
        '--basil',
        '--fs-no-resume',
        f'--fs-subjects-dir={os.path.join(data_dir, "anatomical/freesurfer")}',
        '--derivatives',
        f'{os.path.join(data_dir, "anatomical/smriprep")}',
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.examples_pcasl_singlepld_siemens
def test_examples_pcasl_singlepld_siemens(data_dir, output_dir, working_dir):
    """Run aslprep on the asl_005 ASL-BIDS examples dataset.

    This test uses a Siemens session.
    """
    TEST_NAME = 'examples_pcasl_singlepld_siemens'
    PARTICIPANT_LABEL = '01'

    dataset_dir = download_test_data('examples_pcasl_singlepld', data_dir)
    download_test_data('anatomical', data_dir)
    # Symlink the Freesurfer derivatives to a sub-01_ses-anat folder
    os.symlink(
        os.path.join(data_dir, 'anatomical/freesurfer/sub-01'),
        os.path.join(data_dir, 'anatomical/freesurfer/sub-01_ses-anat'),
    )
    out_dir = os.path.join(output_dir, TEST_NAME, 'aslprep')
    work_dir = os.path.join(working_dir, TEST_NAME)

    parameters = [
        dataset_dir,
        out_dir,
        'participant',
        f'--participant-label={PARTICIPANT_LABEL}',
        '--session-label',
        'siemens3d',
        'anat',
        f'-w={work_dir}',
        '--output-spaces',
        'MNI152NLin2009cAsym',
        '--basil',
        '--m0_scale=10',
        '--fs-no-resume',
        f'--fs-subjects-dir={os.path.join(data_dir, "anatomical/freesurfer")}',
        '--derivatives',
        f'{os.path.join(data_dir, "anatomical/smriprep")}',
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.qtab
def test_qtab(data_dir, output_dir, working_dir):
    """Run aslprep on QTAB data.

    This dataset is Siemens.
    """
    TEST_NAME = 'qtab'
    PARTICIPANT_LABEL = '01'

    dataset_dir = download_test_data(TEST_NAME, data_dir)
    download_test_data('anatomical', data_dir)
    # Symlink the Freesurfer derivatives to a sub-01_ses-anat folder
    os.symlink(
        os.path.join(data_dir, 'anatomical/freesurfer/sub-01'),
        os.path.join(data_dir, 'anatomical/freesurfer/sub-01_ses-01'),
    )
    out_dir = os.path.join(output_dir, TEST_NAME, 'aslprep')
    work_dir = os.path.join(working_dir, TEST_NAME)

    parameters = [
        dataset_dir,
        out_dir,
        'participant',
        f'--participant-label={PARTICIPANT_LABEL}',
        f'-w={work_dir}',
        '--output-spaces',
        'asl',
        'T1w',
        'MNI152NLin2009cAsym',
        '--scorescrub',
        '--use-syn-sdc',
        '--force',
        'no-ge',
        '--subject-anatomical-reference',
        'first-lex',
        '--fs-no-resume',
        f'--fs-subjects-dir={os.path.join(data_dir, "anatomical/freesurfer")}',
        '--derivatives',
        f'{os.path.join(data_dir, "anatomical/smriprep")}',
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.test_001
def test_test_001(data_dir, output_dir, working_dir):
    """Run aslprep on sub-01 data.

    This dataset is Siemens.
    """
    TEST_NAME = 'test_001'
    PARTICIPANT_LABEL = '01'

    dataset_dir = download_test_data(TEST_NAME, data_dir)
    download_test_data('anatomical', data_dir)
    out_dir = os.path.join(output_dir, TEST_NAME, 'aslprep')
    work_dir = os.path.join(working_dir, TEST_NAME)

    parameters = [
        dataset_dir,
        out_dir,
        'participant',
        f'--participant-label={PARTICIPANT_LABEL}',
        f'-w={work_dir}',
        '--output-spaces',
        'asl',
        'T1w',
        'MNI152NLin2009cAsym',
        '--scorescrub',
        '--force',
        'no-ge',
        '--fs-no-resume',
        f'--fs-subjects-dir={os.path.join(data_dir, "anatomical/freesurfer")}',
        '--derivatives',
        f'{os.path.join(data_dir, "anatomical/smriprep")}',
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.test_002
def test_test_002(data_dir, output_dir, working_dir):
    """Run aslprep on sub-01.

    This dataset contains PCASL data from a GE scanner.
    There are two ASL volumes (both deltam) and separate M0 scan.
    """
    TEST_NAME = 'test_002'
    PARTICIPANT_LABEL = '01'

    dataset_dir = download_test_data(TEST_NAME, data_dir)
    download_test_data('anatomical', data_dir)
    out_dir = os.path.join(output_dir, TEST_NAME, 'aslprep')
    work_dir = os.path.join(working_dir, TEST_NAME)

    parameters = [
        dataset_dir,
        out_dir,
        'participant',
        f'--participant-label={PARTICIPANT_LABEL}',
        f'-w={work_dir}',
        '--output-spaces',
        'asl',
        'MNI152NLin2009cAsym',
        '--scorescrub',
        '--use-syn-sdc',
        '--m0_scale=96',
        '--force',
        'ge',
        '--fs-no-resume',
        f'--fs-subjects-dir={os.path.join(data_dir, "anatomical/freesurfer")}',
        '--derivatives',
        f'{os.path.join(data_dir, "anatomical/smriprep")}',
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.test_003_minimal
def test_test_003_minimal(data_dir, output_dir, working_dir):
    """Run ASLPrep minimal workflow on test_003 dataset."""
    base_test_003(
        data_dir,
        output_dir,
        working_dir,
        level='minimal',
        extra_params=['--fs-no-resume'],
    )


@pytest.mark.test_003_resampling
def test_test_003_resampling(data_dir, output_dir, working_dir):
    """Run ASLPrep resampling workflow on test_003 dataset."""
    base_test_003(
        data_dir,
        output_dir,
        working_dir,
        level='resampling',
        extra_params=['--fs-no-resume'],
    )


@pytest.mark.test_003_full
def test_test_003_full(data_dir, output_dir, working_dir):
    """Run ASLPrep full workflow on test_003 dataset."""
    base_test_003(
        data_dir,
        output_dir,
        working_dir,
        level='full',
        extra_params=[
            '--fs-no-resume',
            '--cifti-output',
            '91k',
            '--project-goodvoxels',
            '--no-msm',
        ],
    )


def base_test_003(data_dir, output_dir, working_dir, level, extra_params):
    """Run aslprep on sub-01.

    This dataset is Siemens.
    """
    TEST_NAME = 'test_003'
    PARTICIPANT_LABEL = '01'

    dataset_dir = download_test_data(TEST_NAME, data_dir)
    download_test_data('anatomical', data_dir)
    # Symlink the Freesurfer derivatives to a sub-01_ses-1 folder
    os.symlink(
        os.path.join(data_dir, 'anatomical/freesurfer/sub-01'),
        os.path.join(data_dir, 'anatomical/freesurfer/sub-01_ses-1'),
    )
    level_test_name = f'{TEST_NAME}_{level}'
    out_dir = os.path.join(output_dir, level_test_name, 'aslprep')
    work_dir = os.path.join(working_dir, level_test_name)

    parameters = [
        dataset_dir,
        out_dir,
        'participant',
        f'--participant-label={PARTICIPANT_LABEL}',
        f'-w={work_dir}',
        '--output-spaces',
        'asl',
        '--use-syn-sdc',
        '--m0_scale=10',
        f'--fs-subjects-dir={os.path.join(data_dir, "anatomical/freesurfer")}',
        '--derivatives',
        f'{os.path.join(data_dir, "anatomical/smriprep")}',
        f'--level={level}',
    ]
    parameters += extra_params

    _run_and_generate(level_test_name, PARTICIPANT_LABEL, parameters, out_dir)


def _run_and_generate(test_name, participant_label, parameters, out_dir):
    from aslprep import config

    parameters.append('--sloppy')
    parameters.append('--stop-on-first-crash')
    parameters.append('--clean-workdir')
    parameters.append('-vv')

    # Add concurrency options if they're not already specified
    parameters = update_resources(parameters)

    parse_args(parameters)
    config_file = config.execution.work_dir / f'config-{config.execution.run_uuid}.toml'
    config.loggers.cli.warning(f'Saving config file to {config_file}')
    config.to_filename(config_file)

    retval = build_workflow(config_file, retval={})
    aslprep_wf = retval['workflow']
    aslprep_wf.run(**config.nipype.get_plugin())
    build_boilerplate(str(config_file), aslprep_wf)
    session_list = (
        config.execution.bids_filters.get('asl', {}).get('session')
        if config.execution.bids_filters
        else None
    )
    generate_reports(
        [participant_label],
        out_dir,
        config.execution.run_uuid,
        session_list=session_list,
        bootstrap_file=load_data('reports-spec.yml'),
    )

    output_list_file = os.path.join(get_test_data_path(), f'expected_outputs_{test_name}.txt')
    check_generated_files(out_dir, output_list_file)


def _run_and_fail(parameters):
    from aslprep import config

    parameters.append('--sloppy')
    parameters.append('--stop-on-first-crash')
    parameters.append('-vv')
    parse_args(parameters)
    config_file = config.execution.work_dir / f'config-{config.execution.run_uuid}.toml'
    config.to_filename(config_file)

    with pytest.raises(ValueError, match='Multi-delay data are not supported for PASL'):
        build_workflow(config_file, retval={})
