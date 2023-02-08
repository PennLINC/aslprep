"""Command-line interface tests."""
import os

import pytest

from aslprep.cli.run import build_workflow, get_parser
from aslprep.tests.utils import (check_affines, check_generated_files,
                                 get_test_data_path)


@pytest.mark.sub01
def test_sub01(datasets, output_dir, working_dir):
    """Run aslprep on sub-01 data."""
    test_name = "test_sub01"

    data_dir = datasets["sub01"]
    smriprep_dir = datasets["smriprep"]
    out_dir = os.path.join(output_dir, test_name)
    work_dir = os.path.join(working_dir, test_name)

    test_data_dir = get_test_data_path()
    os.environ["FS_LICENSE"] = os.path.join(test_data_dir, "license.txt")

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
    opts = get_parser().parse_args(parameters)

    retval = {}
    retval = build_workflow(opts, retval=retval)
    # run_uuid = retval.get("run_uuid", None)
    aslprep_wf = retval.get("workflow", None)
    plugin_settings = retval["plugin_settings"]
    aslprep_wf.run(**plugin_settings)

    output_list_file = os.path.join(test_data_dir, "test_outputs_sub01.txt")
    check_generated_files(out_dir, output_list_file)

    check_affines(data_dir, out_dir, input_type="nifti")


@pytest.mark.subA00086748
def test_subA00086748(datasets, output_dir, working_dir):
    """Run aslprep on sub-A00086748."""
    test_name = "test_subA00086748"

    data_dir = datasets["dset"]
    smriprep_dir = datasets["smriprep"]
    out_dir = os.path.join(output_dir, test_name)
    work_dir = os.path.join(working_dir, test_name)

    test_data_dir = get_test_data_path()
    os.environ["FS_LICENSE"] = os.path.join(test_data_dir, "license.txt")

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
    opts = get_parser().parse_args(parameters)

    retval = {}
    retval = build_workflow(opts, retval=retval)
    # run_uuid = retval.get("run_uuid", None)
    aslprep_wf = retval.get("workflow", None)
    plugin_settings = retval["plugin_settings"]
    aslprep_wf.run(**plugin_settings)

    output_list_file = os.path.join(test_data_dir, "test_outputs_subA00086748.txt")
    check_generated_files(out_dir, output_list_file)

    check_affines(data_dir, out_dir, input_type="nifti")
