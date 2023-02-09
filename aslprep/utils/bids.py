# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Utilities to handle BIDS inputs."""
import json
import os
import re
import sys
from pathlib import Path

from bids.layout import BIDSLayout


def collect_data(
    bids_dir,
    participant_label,
    echo=None,
    bids_validate=False,
    bids_filters=None,
):
    """Uses pybids to retrieve the input data for a given participant"""
    if isinstance(bids_dir, BIDSLayout):
        layout = bids_dir
    else:
        layout = BIDSLayout(str(bids_dir), validate=bids_validate)

    queries = {
        "fmap": {"datatype": "fmap"},
        "bold": {"datatype": "func", "suffix": "bold"},
        "sbref": {"datatype": "func", "suffix": "sbref"},
        "flair": {"datatype": "anat", "suffix": "FLAIR"},
        "t2w": {"datatype": "anat", "suffix": "T2w"},
        "t1w": {"datatype": "anat", "suffix": "T1w"},
        "roi": {"datatype": "anat", "suffix": "roi"},
        "asl": {"datatype": "perf", "suffix": "asl"},
        "cbf": {"datatype": "perf", "suffix": "DELTAM"},
        "m0z": {"datatype": "perf", "suffix": "MoScan"},
    }

    bids_filters = bids_filters or {}
    for acq, entities in bids_filters.items():
        queries[acq].update(entities)

    if echo:
        queries["asl"]["echo"] = echo

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

    # Special case: multi-echo BOLD, grouping echos
    if any(["_echo-" in bold for bold in subj_data["asl"]]):
        subj_data["asl"] = group_multiecho(subj_data["asl"])

    return subj_data, layout


def group_multiecho(bold_sess):
    """Multiplexes multi-echo EPIs into arrays.

    Dual-echo is a special case of multi-echo, which is treated as single-echo data.

    >>> bold_sess = ["sub-01_task-rest_echo-1_run-01_bold.nii.gz",
    ...              "sub-01_task-rest_echo-2_run-01_bold.nii.gz",
    ...              "sub-01_task-rest_echo-1_run-02_bold.nii.gz",
    ...              "sub-01_task-rest_echo-2_run-02_bold.nii.gz",
    ...              "sub-01_task-rest_echo-3_run-02_bold.nii.gz",
    ...              "sub-01_task-rest_run-03_bold.nii.gz"]
    >>> group_multiecho(bold_sess)
    ['sub-01_task-rest_echo-1_run-01_bold.nii.gz',
     'sub-01_task-rest_echo-2_run-01_bold.nii.gz',
    ['sub-01_task-rest_echo-1_run-02_bold.nii.gz',
     'sub-01_task-rest_echo-2_run-02_bold.nii.gz',
     'sub-01_task-rest_echo-3_run-02_bold.nii.gz'],
     'sub-01_task-rest_run-03_bold.nii.gz']

    >>> bold_sess.insert(2, "sub-01_task-rest_echo-3_run-01_bold.nii.gz")
    >>> group_multiecho(bold_sess)
    [['sub-01_task-rest_echo-1_run-01_bold.nii.gz',
      'sub-01_task-rest_echo-2_run-01_bold.nii.gz',
      'sub-01_task-rest_echo-3_run-01_bold.nii.gz'],
     ['sub-01_task-rest_echo-1_run-02_bold.nii.gz',
      'sub-01_task-rest_echo-2_run-02_bold.nii.gz',
      'sub-01_task-rest_echo-3_run-02_bold.nii.gz'],
      'sub-01_task-rest_run-03_bold.nii.gz']

    >>> bold_sess += ["sub-01_task-beh_echo-1_run-01_bold.nii.gz",
    ...               "sub-01_task-beh_echo-2_run-01_bold.nii.gz",
    ...               "sub-01_task-beh_echo-1_run-02_bold.nii.gz",
    ...               "sub-01_task-beh_echo-2_run-02_bold.nii.gz",
    ...               "sub-01_task-beh_echo-3_run-02_bold.nii.gz",
    ...               "sub-01_task-beh_run-03_bold.nii.gz"]
    >>> group_multiecho(bold_sess)
    [['sub-01_task-rest_echo-1_run-01_bold.nii.gz',
      'sub-01_task-rest_echo-2_run-01_bold.nii.gz',
      'sub-01_task-rest_echo-3_run-01_bold.nii.gz'],
     ['sub-01_task-rest_echo-1_run-02_bold.nii.gz',
      'sub-01_task-rest_echo-2_run-02_bold.nii.gz',
      'sub-01_task-rest_echo-3_run-02_bold.nii.gz'],
      'sub-01_task-rest_run-03_bold.nii.gz',
      'sub-01_task-beh_echo-1_run-01_bold.nii.gz',
      'sub-01_task-beh_echo-2_run-01_bold.nii.gz',
     ['sub-01_task-beh_echo-1_run-02_bold.nii.gz',
      'sub-01_task-beh_echo-2_run-02_bold.nii.gz',
      'sub-01_task-beh_echo-3_run-02_bold.nii.gz'],
      'sub-01_task-beh_run-03_bold.nii.gz']

    Some tests from https://neurostars.org/t/fmriprep-from\
-singularity-unboundlocalerror/3299/7

    >>> bold_sess = ['sub-01_task-AudLoc_echo-1_bold.nii',
    ...              'sub-01_task-AudLoc_echo-2_bold.nii',
    ...              'sub-01_task-FJT_echo-1_bold.nii',
    ...              'sub-01_task-FJT_echo-2_bold.nii',
    ...              'sub-01_task-LDT_echo-1_bold.nii',
    ...              'sub-01_task-LDT_echo-2_bold.nii',
    ...              'sub-01_task-MotLoc_echo-1_bold.nii',
    ...              'sub-01_task-MotLoc_echo-2_bold.nii']
    >>> group_multiecho(bold_sess) == bold_sess
    True

    >>> bold_sess += ['sub-01_task-MotLoc_echo-3_bold.nii']
    >>> groups = group_multiecho(bold_sess)
    >>> len(groups[:-1])
    6
    >>> [isinstance(g, list) for g in groups]
    [False, False, False, False, False, False, True]
    >>> len(groups[-1])
    3
    """
    from itertools import groupby

    def _grp_echos(x):
        if "_echo-" not in x:
            return x
        echo = re.search("_echo-\\d*", x).group(0)
        return x.replace(echo, "_echo-?")

    ses_uids = []
    for _, bold in groupby(bold_sess, key=_grp_echos):
        bold = list(bold)
        # If single- or dual-echo, flatten list; keep list otherwise.
        action = getattr(ses_uids, "append" if len(bold) > 2 else "extend")
        action(bold)
    return ses_uids


def write_derivative_description(bids_dir, deriv_dir):
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
        desc["SourceDatasetsURLs"] = ["https://doi.org/{}".format(orig_desc["DatasetDOI"])]
    if "License" in orig_desc:
        desc["License"] = orig_desc["License"]

    with (deriv_dir / "dataset_description.json").open("w") as fobj:
        json.dump(desc, fobj, indent=4)


def validate_input_dir(exec_env, bids_dir, participant_label):
    # Ignore issues and warnings that should not influence
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
        all_subs = set([s.name[4:] for s in bids_dir.glob("sub-*")])
        selected_subs = set([s[4:] if s.startswith("sub-") else s for s in participant_label])
        bad_labels = selected_subs.difference(all_subs)
        if bad_labels:
            error_msg = (
                "Data for requested participant(s) label(s) not found. Could "
                "not find data for participant(s): %s. Please verify the requested "
                "participant labels."
            )
            if exec_env == "docker":
                error_msg += (
                    " This error can be caused by the input data not being "
                    "accessible inside the docker container. Please make sure all "
                    "volumes are mounted properly (see https://docs.docker.com/"
                    "engine/reference/commandline/run/#mount-volume--v---read-only)"
                )
            if exec_env == "singularity":
                error_msg += (
                    " This error can be caused by the input data not being "
                    "accessible inside the singularity container. Please make sure "
                    "all paths are mapped properly (see https://www.sylabs.io/"
                    "guides/3.0/user-guide/bind_paths_and_mounts.html)"
                )
            raise RuntimeError(error_msg % ",".join(bad_labels))

        ignored_subs = all_subs.difference(selected_subs)
        if ignored_subs:
            for sub in ignored_subs:
                validator_config_dict["ignoredFiles"].append("/sub-%s/**" % sub)
    with tempfile.NamedTemporaryFile("w+") as temp:
        temp.write(json.dumps(validator_config_dict))
        temp.flush()
        try:
            subprocess.check_call(["bids-validator", bids_dir, "-c", temp.name])
        except FileNotFoundError:
            print("bids-validator does not appear to be installed", file=sys.stderr)


def _get_shub_version(singularity_url):
    return NotImplemented


def fieldmap_wrangler(layout, target_image, use_syn=False, force_syn=False):
    """Query the BIDSLayout for fieldmaps, and arrange them for the orchestration workflow."""
    from collections import defaultdict

    fmap_bids = layout.get_fieldmap(target_image, return_list=True)
    fieldmaps = defaultdict(list)
    for fmap in fmap_bids:
        if fmap["suffix"] == "epi":
            fieldmaps["epi"].append((fmap["epi"], layout.get_metadata(fmap["epi"])))

        if fmap["suffix"] == "fieldmap":
            fieldmaps["fieldmap"].append(
                {
                    "magnitude": [(fmap["magnitude"], layout.get_metadata(fmap["magnitude"]))],
                    "fieldmap": [(fmap["fieldmap"], layout.get_metadata(fmap["fieldmap"]))],
                }
            )

        if fmap["suffix"] == "phasediff":
            fieldmaps["phasediff"].append(
                {
                    "magnitude": [
                        (fmap[k], layout.get_metadata(fmap[k]))
                        for k in sorted(fmap.keys())
                        if k.startswith("magnitude")
                    ],
                    "phases": [(fmap["phasediff"], layout.get_metadata(fmap["phasediff"]))],
                }
            )

        if fmap["suffix"] == "phase":
            fieldmaps["phasediff"].append(
                {
                    "magnitude": [
                        (fmap[k], layout.get_metadata(fmap[k]))
                        for k in sorted(fmap.keys())
                        if k.startswith("magnitude")
                    ],
                    "phases": [
                        (fmap[k], layout.get_metadata(fmap[k]))
                        for k in sorted(fmap.keys())
                        if k.startswith("phase")
                    ],
                }
            )

    if fieldmaps and force_syn:
        # syn: True -> Run SyN in addition to fieldmap-based SDC
        fieldmaps["syn"] = True
    elif not fieldmaps and (force_syn or use_syn):
        # syn: False -> Run SyN as only SDC
        fieldmaps["syn"] = False
    return fieldmaps
