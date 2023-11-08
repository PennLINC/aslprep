# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Utilities to handle BIDS inputs."""
import json
import os
import sys
from pathlib import Path

import yaml

from aslprep import config


def collect_data(
    layout,
    participant_label,
    bids_filters=None,
):
    """Use pybids to retrieve the input data for a given participant."""
    queries = {
        "fmap": {"datatype": "fmap"},
        "flair": {"datatype": "anat", "suffix": "FLAIR"},
        "t2w": {"datatype": "anat", "suffix": "T2w"},
        "t1w": {"datatype": "anat", "suffix": "T1w"},
        "roi": {"datatype": "anat", "suffix": "roi"},
        "sbref": {"datatype": "perf", "suffix": "sbref"},
        "asl": {"datatype": "perf", "suffix": "asl"},
    }

    bids_filters = bids_filters or {}
    for acq, entities in bids_filters.items():
        queries[acq].update(entities)

    subj_data = {
        dtype: sorted(
            layout.get(
                return_type="file",
                subject=participant_label,
                extension=["nii", "nii.gz"],
                **query,
            )
        )
        for dtype, query in queries.items()
    }

    return subj_data


def collect_run_data(layout, asl_file):
    """Use pybids to retrieve the input data for a given participant."""
    queries = {
        "aslcontext": {"suffix": "aslcontext", "extension": ".tsv"},
        "sbref": {"suffix": "sbref", "extension": [".nii", ".nii.gz"]},
    }

    bids_file = layout.get_file(asl_file)

    run_data = {
        dtype: layout.get_nearest(
            bids_file.path,
            return_type="file",
            strict=False,  # aslcontext files aren't grabbed when strict=True, for some reason
            **query,
        )
        for dtype, query in queries.items()
    }

    if "sbref" in config.workflow.ignore:
        config.loggers.workflow.info("Single-band reference files ignored.")
        run_data["sbref"] = None

    # The aslcontext file is required
    if not run_data["aslcontext"]:
        raise FileNotFoundError(f"aslcontext file for {asl_file} not found.")

    # Now let's look for an m0scan
    m0scan_candidates = [
        f
        for f in bids_file.get_associations(kind="InformedBy")
        if f.entities["suffix"] == "m0scan"
    ]
    if m0scan_candidates:
        if len(m0scan_candidates) > 1:
            config.loggers.workflow.warning(
                f"More than one M0 file found for {asl_file}. "
                f"Using the first one ({m0scan_candidates[0].path})"
            )
        run_data["m0scan"] = m0scan_candidates[0].path
    else:
        run_data["m0scan"] = None

    m0scan_metadata = None
    asl_metadata = layout.get_metadata(asl_file)
    if (asl_metadata["M0Type"] == "Separate") and not run_data["m0scan"]:
        raise FileNotFoundError(f"M0 file for {asl_file} not found.")
    elif asl_metadata["M0Type"] == "Separate":
        m0scan_metadata = layout.get_file(run_data["m0scan"]).get_metadata()
        if not m0scan_metadata:
            raise Exception(f"No metadata for m0scan: {run_data['m0scan']}")
    elif run_data["m0scan"]:
        raise ValueError(
            f"M0Type is {run_data['asl_metadata']['M0Type']}, "
            f"but an M0 scan was found at {run_data['m0scan']}"
        )

    config.loggers.workflow.info(
        (
            f"Collected run data for {asl_file}:\n"
            f"{yaml.dump(run_data, default_flow_style=False, indent=4)}"
        ),
    )

    # Add metadata to dictionary now (we don't want to print these with the logger).
    run_data["asl_metadata"] = asl_metadata
    run_data["m0scan_metadata"] = m0scan_metadata

    return run_data


def write_bidsignore(deriv_dir):
    """Write .bidsignore file."""
    bids_ignore = (
        "*.html",
        "logs/",
        "figures/",  # Reports
        "*_xfm.*",  # Unspecified transform files
        "*.surf.gii",  # Unspecified structural outputs
        # Unspecified functional outputs
        "*_aslref.nii.gz",
        "*_asl.func.gii",
        "*_mixing.tsv",
        "*_timeseries.tsv",
    )
    ignore_file = Path(deriv_dir) / ".bidsignore"

    ignore_file.write_text("\n".join(bids_ignore) + "\n")


def write_derivative_description(bids_dir, deriv_dir):
    """Write derivative dataset_description file."""
    from aslprep.__about__ import DOWNLOAD_URL, __url__, __version__

    bids_dir = Path(bids_dir)
    deriv_dir = Path(deriv_dir)
    desc = {
        "Name": "ASLPrep - ASL PREProcessing workflow",
        "BIDSVersion": "1.1.1",
        "PipelineDescription": {
            "Name": "ASLPrep",
            "Version": __version__,
            "CodeURL": DOWNLOAD_URL,
        },
        "CodeURL": __url__,
        "HowToAcknowledge": "Please cite our paper coming :), "
        "and include the generated citation boilerplate within the Methods "
        "section of the text.",
    }

    # Keys that can only be set by environment
    if "ASLPREP_DOCKER_TAG" in os.environ:
        desc["DockerHubContainerTag"] = os.environ["ASLPREP_DOCKER_TAG"]

    if "ASLPREP_SINGULARITY_URL" in os.environ:
        singularity_url = os.environ["ASLPREP_SINGULARITY_URL"]
        desc["SingularityContainerURL"] = singularity_url

        singularity_md5 = _get_shub_version(singularity_url)
        if singularity_md5 and singularity_md5 is not NotImplemented:
            desc["SingularityContainerMD5"] = _get_shub_version(singularity_url)

    # Keys deriving from source dataset
    orig_desc = {}
    fname = bids_dir / "dataset_description.json"
    if fname.exists():
        with fname.open() as fobj:
            orig_desc = json.load(fobj)

    if "DatasetDOI" in orig_desc:
        desc["SourceDatasetsURLs"] = [f"https://doi.org/{orig_desc['DatasetDOI']}"]

    if "License" in orig_desc:
        desc["License"] = orig_desc["License"]

    with (deriv_dir / "dataset_description.json").open("w") as fobj:
        json.dump(desc, fobj, indent=4)


def validate_input_dir(exec_env, bids_dir, participant_label):
    """Validate input directory.

    This function checks necessary validator requirements,
    but skips ones that should not influence the pipeline.
    """
    import subprocess
    import tempfile

    validator_config_dict = {
        "ignore": [
            "EVENTS_COLUMN_ONSET",
            "EVENTS_COLUMN_DURATION",
            "TSV_EQUAL_ROWS",
            "TSV_EMPTY_CELL",
            "TSV_IMPROPER_NA",
            "VOLUME_COUNT_MISMATCH",
            "BVAL_MULTIPLE_ROWS",
            "BVEC_NUMBER_ROWS",
            "DWI_MISSING_BVAL",
            "INCONSISTENT_SUBJECTS",
            "INCONSISTENT_PARAMETERS",
            "BVEC_ROW_LENGTH",
            "B_FILE",
            "PARTICIPANT_ID_COLUMN",
            "PARTICIPANT_ID_MISMATCH",
            "TASK_NAME_MUST_DEFINE",
            "PHENOTYPE_SUBJECTS_MISSING",
            "STIMULUS_FILE_MISSING",
            "DWI_MISSING_BVEC",
            "EVENTS_TSV_MISSING",
            "TSV_IMPROPER_NA",
            "ACQTIME_FMT",
            "Participants age 89 or higher",
            "DATASET_DESCRIPTION_JSON_MISSING",
            "FILENAME_COLUMN",
            "WRONG_NEW_LINE",
            "MISSING_TSV_COLUMN_CHANNELS",
            "MISSING_TSV_COLUMN_IEEG_CHANNELS",
            "MISSING_TSV_COLUMN_IEEG_ELECTRODES",
            "UNUSED_STIMULUS",
            "CHANNELS_COLUMN_SFREQ",
            "CHANNELS_COLUMN_LOWCUT",
            "CHANNELS_COLUMN_HIGHCUT",
            "CHANNELS_COLUMN_NOTCH",
            "CUSTOM_COLUMN_WITHOUT_DESCRIPTION",
            "ACQTIME_FMT",
            "SUSPICIOUSLY_LONG_EVENT_DESIGN",
            "SUSPICIOUSLY_SHORT_EVENT_DESIGN",
            "MALFORMED_BVEC",
            "MALFORMED_BVAL",
            "MISSING_TSV_COLUMN_EEG_ELECTRODES",
            "MISSING_SESSION",
        ],
        "error": ["NO_T1W"],
        "ignoredFiles": ["/dataset_description.json", "/participants.tsv"],
    }
    # Limit validation only to data from requested participants
    if participant_label:
        all_subs = {s.name[4:] for s in bids_dir.glob("sub-*")}
        selected_subs = {s[4:] if s.startswith("sub-") else s for s in participant_label}
        if bad_labels := selected_subs.difference(all_subs):
            error_msg = (
                "Data for requested participant(s) label(s) not found. "
                f"Could not find data for participant(s): {','.join(bad_labels)}. "
                "Please verify the requested participant labels."
            )
            if exec_env == "docker":
                error_msg += (
                    " This error can be caused by the input data not being "
                    "accessible inside the docker container. "
                    "Please make sure all volumes are mounted properly "
                    "(see https://docs.docker.com/engine/reference/commandline/"
                    "run/#mount-volume--v---read-only)"
                )

            elif exec_env == "singularity":
                error_msg += (
                    " This error can be caused by the input data not being "
                    "accessible inside the singularity container. "
                    "Please make sure all paths are mapped properly "
                    "(see https://www.sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html)"
                )

            raise RuntimeError(error_msg)

        if ignored_subs := all_subs.difference(selected_subs):
            for sub in ignored_subs:
                validator_config_dict["ignoredFiles"].append(f"/sub-{sub}/**")

    with tempfile.NamedTemporaryFile("w+") as temp:
        temp.write(json.dumps(validator_config_dict))
        temp.flush()
        try:
            subprocess.check_call(["bids-validator", bids_dir, "-c", temp.name])
        except FileNotFoundError:
            print("bids-validator does not appear to be installed", file=sys.stderr)


def _get_shub_version(singularity_url):  # noqa: U100
    return NotImplemented


def find_atlas_entities(filename):
    """Extract atlas entities from filename."""
    import os

    fname = os.path.basename(filename)
    elements = fname.split("_")

    out = []
    for ent in ("tpl", "atlas", "res"):
        ent_parts = [el for el in elements if el.startswith(f"{ent}-")]
        ent_value = None
        if ent_parts:
            ent_value = ent_parts[0].split("-")[1]

        out.append(ent_value)

    suffix = elements[-1].split(".")[0]
    extension = "." + ".".join(elements[-1].split(".")[1:])
    out += [suffix, extension]

    return tuple(out)


def get_estimator(layout, fname):
    """Identify estimator for a file."""
    field_source = layout.get_metadata(fname).get("B0FieldSource")
    if isinstance(field_source, str):
        field_source = (field_source,)

    if field_source is None:
        import re
        from pathlib import Path

        from sdcflows.fieldmaps import get_identifier

        # Fallback to IntendedFor
        intended_rel = re.sub(r"^sub-[a-zA-Z0-9]*/", "", str(Path(fname).relative_to(layout.root)))
        field_source = get_identifier(intended_rel)

    return field_source
