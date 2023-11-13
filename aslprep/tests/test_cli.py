"""Command-line interface tests."""
import os

import pytest
from nipype import config as nipype_config
from niworkflows.reports.core import generate_reports
from pkg_resources import resource_filename as pkgrf

from aslprep.cli.parser import parse_args
from aslprep.cli.workflow import build_boilerplate, build_workflow
from aslprep.tests.utils import (
    check_generated_files,
    download_test_data,
    get_test_data_path,
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
    TEST_NAME = "examples_pasl_multipld"
    PARTICIPANT_LABEL = "01"

    dataset_dir = download_test_data(TEST_NAME, data_dir)
    out_dir = os.path.join(output_dir, TEST_NAME, "aslprep")
    work_dir = os.path.join(working_dir, TEST_NAME)

    parameters = [
        dataset_dir,
        out_dir,
        "participant",
        f"--participant-label={PARTICIPANT_LABEL}",
        f"-w={work_dir}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces=asl",
        "--scorescrub",
        "--basil",
        "--use-syn-sdc",
        "--m0_scale=10",
        "--fs-no-reconall",
        f"--anat-derivatives={os.path.join(dataset_dir, 'derivatives/smriprep')}",
    ]

    _run_and_fail(parameters)


@pytest.mark.examples_pcasl_multipld
def test_examples_pcasl_multipld(data_dir, output_dir, working_dir):
    """Run aslprep on the asl_004 ASL-BIDS examples dataset.

    This dataset has 48 control-label pairs at 6 different PLDs, along with a separate M0 scan.
    The manufacturer is Siemens.
    """
    TEST_NAME = "examples_pcasl_multipld"
    PARTICIPANT_LABEL = "01"

    dataset_dir = download_test_data(TEST_NAME, data_dir)
    out_dir = os.path.join(output_dir, TEST_NAME, "aslprep")
    work_dir = os.path.join(working_dir, TEST_NAME)

    parameters = [
        dataset_dir,
        out_dir,
        "participant",
        f"--participant-label={PARTICIPANT_LABEL}",
        f"-w={work_dir}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces=asl",
        "--scorescrub",
        "--basil",
        "--m0_scale=10",
        "--fs-no-reconall",
        f"--anat-derivatives={os.path.join(dataset_dir, 'derivatives/smriprep')}",
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.examples_pcasl_singlepld_ge
def test_examples_pcasl_singlepld_ge(data_dir, output_dir, working_dir):
    """Run aslprep on the asl_001 ASL-BIDS examples dataset.

    This test uses a GE session with two volumes: one deltam and one M0.
    """
    TEST_NAME = "examples_pcasl_singlepld_ge"
    PARTICIPANT_LABEL = "103"

    dataset_dir = download_test_data("examples_pcasl_singlepld", data_dir)
    out_dir = os.path.join(output_dir, TEST_NAME, "aslprep")
    work_dir = os.path.join(working_dir, TEST_NAME)
    test_data_dir = get_test_data_path()
    filter_file = os.path.join(test_data_dir, f"{TEST_NAME}_filter.json")

    parameters = [
        dataset_dir,
        out_dir,
        "participant",
        f"--participant-label={PARTICIPANT_LABEL}",
        f"-w={work_dir}",
        f"--bids-filter-file={filter_file}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces=asl",
        "--scorescrub",
        "--basil",
        "--m0_scale=96",
        "--fs-no-reconall",
        f"--anat-derivatives={os.path.join(dataset_dir, 'derivatives/smriprep')}",
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.examples_pcasl_singlepld_philips
def test_examples_pcasl_singlepld_philips(data_dir, output_dir, working_dir):
    """Run aslprep on the asl_002 ASL-BIDS examples dataset.

    This test uses a Philips session.
    The appropriate M0 scale is unknown for this dataset, so CBF values will be inflated.
    """
    TEST_NAME = "examples_pcasl_singlepld_philips"
    PARTICIPANT_LABEL = "103"

    dataset_dir = download_test_data("examples_pcasl_singlepld", data_dir)
    out_dir = os.path.join(output_dir, TEST_NAME, "aslprep")
    work_dir = os.path.join(working_dir, TEST_NAME)
    test_data_dir = get_test_data_path()
    filter_file = os.path.join(test_data_dir, f"{TEST_NAME}_filter.json")

    parameters = [
        dataset_dir,
        out_dir,
        "participant",
        f"--participant-label={PARTICIPANT_LABEL}",
        f"-w={work_dir}",
        f"--bids-filter-file={filter_file}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces=asl",
        "--scorescrub",
        "--basil",
        "--fs-no-reconall",
        f"--anat-derivatives={os.path.join(dataset_dir, 'derivatives/smriprep')}",
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.examples_pcasl_singlepld_siemens
def test_examples_pcasl_singlepld_siemens(data_dir, output_dir, working_dir):
    """Run aslprep on the asl_005 ASL-BIDS examples dataset.

    This test uses a Siemens session.
    """
    TEST_NAME = "examples_pcasl_singlepld_siemens"
    PARTICIPANT_LABEL = "103"

    dataset_dir = download_test_data("examples_pcasl_singlepld", data_dir)
    out_dir = os.path.join(output_dir, TEST_NAME, "aslprep")
    work_dir = os.path.join(working_dir, TEST_NAME)
    test_data_dir = get_test_data_path()
    filter_file = os.path.join(test_data_dir, f"{TEST_NAME}_filter.json")

    parameters = [
        dataset_dir,
        out_dir,
        "participant",
        f"--participant-label={PARTICIPANT_LABEL}",
        f"-w={work_dir}",
        f"--bids-filter-file={filter_file}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces",
        "asl",
        "T1w",
        "--basil",
        "--m0_scale=10",
        "--fs-no-reconall",
        f"--anat-derivatives={os.path.join(dataset_dir, 'derivatives/smriprep')}",
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.test_001
def test_test_001(data_dir, output_dir, working_dir):
    """Run aslprep on sub-01 data.

    This dataset is Siemens.
    """
    TEST_NAME = "test_001"
    PARTICIPANT_LABEL = "01"

    dataset_dir = download_test_data(TEST_NAME, data_dir)
    out_dir = os.path.join(output_dir, TEST_NAME, "aslprep")
    work_dir = os.path.join(working_dir, TEST_NAME)

    parameters = [
        dataset_dir,
        out_dir,
        "participant",
        f"--participant-label={PARTICIPANT_LABEL}",
        f"-w={work_dir}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces",
        "asl",
        "T1w",
        "MNI152NLin2009cAsym",
        "--scorescrub",
        "--use-syn-sdc",
        "--force-no-ge",
        "--fs-no-reconall",
        f"--anat-derivatives={os.path.join(dataset_dir, 'derivatives/smriprep')}",
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.test_002
def test_test_002(data_dir, output_dir, working_dir):
    """Run aslprep on sub-10R01383.

    This dataset contains PCASL data from a GE scanner.
    There are two ASL volumes (both deltam) and separate M0 scan.
    """
    TEST_NAME = "test_002"
    PARTICIPANT_LABEL = "10R01383"

    dataset_dir = download_test_data(TEST_NAME, data_dir)
    out_dir = os.path.join(output_dir, TEST_NAME, "aslprep")
    work_dir = os.path.join(working_dir, TEST_NAME)

    parameters = [
        dataset_dir,
        out_dir,
        "participant",
        f"--participant-label={PARTICIPANT_LABEL}",
        f"-w={work_dir}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces",
        "asl",
        "MNI152NLin2009cAsym",
        "--scorescrub",
        "--use-syn-sdc",
        "--m0_scale=96",
        "--force-ge",
        "--fs-no-reconall",
        f"--anat-derivatives={os.path.join(dataset_dir, 'derivatives/smriprep')}",
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


@pytest.mark.test_003
def test_test_003(data_dir, output_dir, working_dir):
    """Run aslprep on sub-A00086748.

    This dataset is Siemens.
    """
    TEST_NAME = "test_003"
    PARTICIPANT_LABEL = "A00086748"

    dataset_dir = download_test_data(TEST_NAME, data_dir)
    out_dir = os.path.join(output_dir, TEST_NAME, "aslprep")
    work_dir = os.path.join(working_dir, TEST_NAME)

    parameters = [
        dataset_dir,
        out_dir,
        "participant",
        f"--participant-label={PARTICIPANT_LABEL}",
        f"-w={work_dir}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces=asl",
        "--use-syn-sdc",
        "--m0_scale=10",
        "--fs-no-reconall",
        f"--anat-derivatives={os.path.join(dataset_dir, 'derivatives/smriprep')}",
    ]

    _run_and_generate(TEST_NAME, PARTICIPANT_LABEL, parameters, out_dir)


def _run_and_generate(test_name, participant_label, parameters, out_dir):
    from aslprep import config

    parameters.append("--stop-on-first-crash")
    parameters.append("--clean-workdir")
    parse_args(parameters)
    config_file = config.execution.work_dir / f"config-{config.execution.run_uuid}.toml"
    config.loggers.cli.warning(f"Saving config file to {config_file}")
    config.to_filename(config_file)

    retval = build_workflow(config_file, retval={})
    aslprep_wf = retval["workflow"]
    aslprep_wf.run()
    build_boilerplate(str(config_file), aslprep_wf)
    generate_reports(
        [participant_label],
        out_dir,
        config.execution.run_uuid,
        config=pkgrf("aslprep", "data/reports-spec.yml"),
        packagename="aslprep",
    )

    output_list_file = os.path.join(get_test_data_path(), f"expected_outputs_{test_name}.txt")
    check_generated_files(out_dir, output_list_file)


def _run_and_fail(parameters):
    from aslprep import config

    parameters.append("--stop-on-first-crash")
    parse_args(parameters)
    config_file = config.execution.work_dir / f"config-{config.execution.run_uuid}.toml"
    config.to_filename(config_file)

    with pytest.raises(ValueError, match="Multi-delay data are not supported for PASL"):
        build_workflow(config_file, retval={})
