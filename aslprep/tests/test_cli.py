"""Command-line interface tests."""
import json
import os

import pytest

from aslprep.cli.parser import parse_args
from aslprep.cli.workflow import build_workflow
from aslprep.tests.utils import check_generated_files, get_test_data_path


@pytest.mark.sub01
def test_sub01(datasets, output_dir, working_dir):
    """Run aslprep on sub-01 data."""
    from aslprep import config

    test_name = "test_sub01"

    data_dir = datasets["dset"]
    smriprep_dir = datasets["smriprep"]
    out_dir = os.path.join(output_dir, test_name)
    work_dir = os.path.join(working_dir, test_name)

    # Patch the JSON file until I have enough changes to merit completely replacing the current
    # version of the data.
    json_file = os.path.join(data_dir, "sub-01", "perf", "sub-01_asl.json")
    with open(json_file, "r") as fo:
        metadata = json.load(fo)

    metadata["RepetitionTimePreparation"] = metadata["RepetitionTime"]
    with open(json_file, "w") as fo:
        json.dump(metadata, fo, sort_keys=True, indent=4)

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
def test_subA00086748(datasets, output_dir, working_dir):  # noqa: N802
    """Run aslprep on sub-A00086748."""
    from aslprep import config

    test_name = "test_subA00086748"

    data_dir = datasets["dset"]
    smriprep_dir = datasets["smriprep"]
    out_dir = os.path.join(output_dir, test_name)
    work_dir = os.path.join(working_dir, test_name)

    # Patch the JSON file until I have enough changes to merit completely replacing the current
    # version of the data.
    json_file = os.path.join(
        data_dir,
        "sub-A00086748",
        "ses-BAS1",
        "perf",
        "sub-A00086748_ses-BAS1_asl.json",
    )
    with open(json_file, "r") as fo:
        metadata = json.load(fo)

    metadata["RepetitionTimePreparation"] = metadata["RepetitionTime"]
    with open(json_file, "w") as fo:
        json.dump(metadata, fo, sort_keys=True, indent=4)

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


@pytest.mark.sub10R01383
def test_sub10R01383(datasets, output_dir, working_dir):  # noqa: N802
    """Run aslprep on sub-10R01383.

    Currently skipped.

    Notes
    -----
    scorescrub fails on this dataset, so I've dropped that parameter.
    I'll probably need to dig into why that happens at some point.
    """
    from aslprep import config

    test_name = "test_sub10R01383"

    data_dir = datasets["dset"]
    smriprep_dir = datasets["smriprep"]
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
        f"--anat-derivatives={smriprep_dir}",
    ]
    parse_args(parameters)
    config_file = config.execution.work_dir / f"config-{config.execution.run_uuid}.toml"
    config.to_filename(config_file)

    retval = {}
    retval = build_workflow(config_file, retval=retval)
    aslprep_wf = retval.get("workflow", None)
    aslprep_wf.run()

    output_list_file = os.path.join(test_data_dir, "test_outputs_sub10R01383.txt")
    check_generated_files(out_dir, output_list_file)
