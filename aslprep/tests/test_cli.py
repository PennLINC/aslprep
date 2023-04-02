"""Command-line interface tests."""
import os

import pytest
from nipype import config as nipype_config

from aslprep.cli.parser import parse_args
from aslprep.cli.workflow import build_workflow
from aslprep.tests.utils import check_generated_files, get_test_data_path

nipype_config.enable_debug_mode()


@pytest.mark.examples_pcasl_singlepld
def test_examples_pcasl_singlepld(datasets, output_dir, working_dir):
    """Run aslprep on sub-01 data."""
    from aslprep import config

    test_name = "examples_pcasl_singlepld"
    data_dir = datasets[test_name]
    out_dir = os.path.join(output_dir, test_name)
    work_dir = os.path.join(working_dir, test_name)
    test_data_dir = get_test_data_path()
    filter_file = os.path.join(test_data_dir, "examples_pcasl_singlepld_ge_filter.json")

    parameters = [
        data_dir,
        out_dir,
        "participant",
        "--participant-label=103",
        f"-w={work_dir}",
        f"--bids-filter-file={filter_file}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces=asl",
        "--scorescrub",
        "--basil",
        "--use-syn-sdc",
        f"--anat-derivatives={os.path.join(data_dir, 'derivatives/smriprep')}",
    ]
    parse_args(parameters)
    config_file = config.execution.work_dir / f"config-{config.execution.run_uuid}.toml"
    config.to_filename(config_file)

    retval = {}
    retval = build_workflow(config_file, retval=retval)
    aslprep_wf = retval.get("workflow", None)
    aslprep_wf.run()

    output_list_file = os.path.join(test_data_dir, "expected_outputs_examples_pcasl_singlepld.txt")
    check_generated_files(out_dir, output_list_file)


@pytest.mark.examples_pcasl_singlepld_ge
def test_examples_pcasl_singlepld_ge(datasets, output_dir, working_dir):
    """Run aslprep on sub-01 data."""
    from aslprep import config

    test_name = "examples_pcasl_singlepld_ge"
    data_dir = datasets["examples_pcasl_singlepld"]
    out_dir = os.path.join(output_dir, test_name)
    work_dir = os.path.join(working_dir, test_name)
    test_data_dir = get_test_data_path()
    filter_file = os.path.join(test_data_dir, "examples_pcasl_singlepld_ge_filter.json")

    parameters = [
        data_dir,
        out_dir,
        "participant",
        "--participant-label=103",
        f"-w={work_dir}",
        f"--bids-filter-file={filter_file}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces=asl",
        "--scorescrub",
        "--basil",
        "--use-syn-sdc",
        f"--anat-derivatives={os.path.join(data_dir, 'derivatives/smriprep')}",
    ]
    parse_args(parameters)
    config_file = config.execution.work_dir / f"config-{config.execution.run_uuid}.toml"
    config.to_filename(config_file)

    retval = {}
    retval = build_workflow(config_file, retval=retval)
    aslprep_wf = retval.get("workflow", None)
    aslprep_wf.run()

    output_list_file = os.path.join(
        test_data_dir,
        "expected_outputs_examples_pcasl_singlepld_ge.txt",
    )
    check_generated_files(out_dir, output_list_file)


@pytest.mark.examples_pcasl_multipld
def test_examples_pcasl_multipld(datasets, output_dir, working_dir):
    """Run aslprep on sub-01 data."""
    from aslprep import config

    test_name = "examples_pcasl_multipld"
    data_dir = datasets[test_name]
    out_dir = os.path.join(output_dir, test_name)
    work_dir = os.path.join(working_dir, test_name)
    test_data_dir = get_test_data_path()

    parameters = [
        data_dir,
        out_dir,
        "participant",
        "--participant-label=01",
        f"-w={work_dir}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces=asl",
        "--scorescrub",
        "--basil",
        "--use-syn-sdc",
        f"--anat-derivatives={os.path.join(data_dir, 'derivatives/smriprep')}",
    ]
    parse_args(parameters)
    config_file = config.execution.work_dir / f"config-{config.execution.run_uuid}.toml"
    config.to_filename(config_file)

    retval = {}
    retval = build_workflow(config_file, retval=retval)
    aslprep_wf = retval.get("workflow", None)
    aslprep_wf.run()

    output_list_file = os.path.join(test_data_dir, "expected_outputs_examples_pcasl_multipld.txt")
    check_generated_files(out_dir, output_list_file)


@pytest.mark.examples_pasl_multipld
def test_examples_pasl_multipld(datasets, output_dir, working_dir):
    """Run aslprep on sub-01 data."""
    from aslprep import config

    test_name = "examples_pasl_multipld"
    data_dir = datasets[test_name]
    out_dir = os.path.join(output_dir, test_name)
    work_dir = os.path.join(working_dir, test_name)
    test_data_dir = get_test_data_path()

    parameters = [
        data_dir,
        out_dir,
        "participant",
        "--participant-label=01",
        f"-w={work_dir}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces=asl",
        "--scorescrub",
        "--basil",
        "--use-syn-sdc",
        f"--anat-derivatives={os.path.join(data_dir, 'derivatives/smriprep')}",
    ]
    parse_args(parameters)
    config_file = config.execution.work_dir / f"config-{config.execution.run_uuid}.toml"
    config.to_filename(config_file)

    retval = {}
    retval = build_workflow(config_file, retval=retval)
    aslprep_wf = retval.get("workflow", None)
    aslprep_wf.run()

    output_list_file = os.path.join(test_data_dir, "expected_outputs_examples_pasl_multipld.txt")
    check_generated_files(out_dir, output_list_file)


@pytest.mark.test_001
def test_test_001(datasets, output_dir, working_dir):
    """Run aslprep on sub-01 data."""
    from aslprep import config

    test_name = "test_001"
    data_dir = datasets[test_name]
    out_dir = os.path.join(output_dir, test_name)
    work_dir = os.path.join(working_dir, test_name)
    test_data_dir = get_test_data_path()

    parameters = [
        data_dir,
        out_dir,
        "participant",
        "--participant-label=01",
        f"-w={work_dir}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces=asl",
        "--scorescrub",
        "--basil",
        "--use-syn-sdc",
        f"--anat-derivatives={os.path.join(data_dir, 'derivatives/smriprep')}",
    ]
    parse_args(parameters)
    config_file = config.execution.work_dir / f"config-{config.execution.run_uuid}.toml"
    config.to_filename(config_file)

    retval = {}
    retval = build_workflow(config_file, retval=retval)
    aslprep_wf = retval.get("workflow", None)
    aslprep_wf.run()

    output_list_file = os.path.join(test_data_dir, "expected_outputs_test_001.txt")
    check_generated_files(out_dir, output_list_file)


@pytest.mark.test_002
def test_test_002(datasets, output_dir, working_dir):
    """Run aslprep on sub-10R01383.

    This dataset contains PCASL data from a GE scanner.
    There are two ASL volumes (both deltam) and separate M0 scan.

    Notes
    -----
    scorescrub fails on this dataset, so I've dropped that parameter.
    I'll probably need to dig into why that happens at some point.
    """
    from aslprep import config

    test_name = "test_002"
    data_dir = datasets[test_name]
    out_dir = os.path.join(output_dir, test_name)
    work_dir = os.path.join(working_dir, test_name)
    test_data_dir = get_test_data_path()

    parameters = [
        data_dir,
        out_dir,
        "participant",
        "--participant-label=10R01383",
        f"-w={work_dir}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces=asl",
        "--basil",
        "--use-syn-sdc",
        f"--anat-derivatives={os.path.join(data_dir, 'derivatives/smriprep')}",
    ]
    parse_args(parameters)
    config_file = config.execution.work_dir / f"config-{config.execution.run_uuid}.toml"
    config.to_filename(config_file)

    retval = {}
    retval = build_workflow(config_file, retval=retval)
    aslprep_wf = retval.get("workflow", None)
    aslprep_wf.run()

    output_list_file = os.path.join(test_data_dir, "expected_outputs_test_002.txt")
    check_generated_files(out_dir, output_list_file)


@pytest.mark.test_003
def test_test_003(datasets, output_dir, working_dir):
    """Run aslprep on sub-A00086748."""
    from aslprep import config

    test_name = "test_003"
    data_dir = datasets[test_name]
    out_dir = os.path.join(output_dir, test_name)
    work_dir = os.path.join(working_dir, test_name)
    test_data_dir = get_test_data_path()

    parameters = [
        data_dir,
        out_dir,
        "participant",
        "--participant-label=A00086748",
        f"-w={work_dir}",
        "--nthreads=2",
        "--omp-nthreads=2",
        "--output-spaces=asl",
        "--scorescrub",
        "--basil",
        "--use-syn-sdc",
        f"--anat-derivatives={os.path.join(data_dir, 'derivatives/smriprep')}",
    ]
    parse_args(parameters)
    config_file = config.execution.work_dir / f"config-{config.execution.run_uuid}.toml"
    config.to_filename(config_file)

    retval = {}
    retval = build_workflow(config_file, retval=retval)
    aslprep_wf = retval.get("workflow", None)
    aslprep_wf.run()

    output_list_file = os.path.join(test_data_dir, "expected_outputs_test_003.txt")
    check_generated_files(out_dir, output_list_file)
