"""Command-line interface tests."""
import os

import pytest

from aslprep.cli.parser import parse_args
from aslprep.cli.workflow import build_workflow
from aslprep.tests.utils import check_generated_files, get_test_data_path


@pytest.mark.sub01
def test_sub01(datasets, output_dir, working_dir):
    """Run aslprep on sub-01 data."""
    return
    from aslprep import config

    test_name = "test_sub01"

    data_dir = datasets["dset"]
    smriprep_dir = datasets["smriprep"]
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
        f"--anat-derivatives={smriprep_dir}",
    ]
    parse_args(parameters)
    config_file = config.execution.work_dir / f"config-{config.execution.run_uuid}.toml"
    config.to_filename(config_file)

    retval = {}
    retval = build_workflow(config_file, retval=retval)
    aslprep_wf = retval.get("workflow", None)
    aslprep_wf.run()

    output_list_file = os.path.join(test_data_dir, "test_outputs_sub01.txt")
    check_generated_files(out_dir, output_list_file)


@pytest.mark.subA00086748
def test_subA00086748(datasets, output_dir, working_dir):
    """Run aslprep on sub-A00086748."""
    return
    from aslprep import config

    test_name = "test_subA00086748"

    data_dir = datasets["dset"]
    smriprep_dir = datasets["smriprep"]
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
        f"--anat-derivatives={smriprep_dir}",
    ]
    parse_args(parameters)
    config_file = config.execution.work_dir / f"config-{config.execution.run_uuid}.toml"
    config.to_filename(config_file)

    retval = {}
    retval = build_workflow(config_file, retval=retval)
    aslprep_wf = retval.get("workflow", None)
    aslprep_wf.run()

    output_list_file = os.path.join(test_data_dir, "test_outputs_subA00086748.txt")
    check_generated_files(out_dir, output_list_file)
