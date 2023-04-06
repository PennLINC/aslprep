"""Functions for calculating and collecting confounds."""
import os
import re

import pandas as pd
from nipype.interfaces.base import isdefined

from aslprep import config


def _less_breakable(a_string):
    """Harden the string to different environments, whatever that means."""
    return "".join(a_string.split()).strip("#")


def _camel_to_snake(name):
    """Convert camelCase string to snake_case.

    Taken from https://stackoverflow.com/questions/1175208/.
    If we end up using it more than just here, probably worth pulling in a well-tested package.
    """
    s1 = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub("([a-z0-9])([A-Z])", r"\1_\2", s1).lower()


def _adjust_indices(left_df, right_df):
    """Force missing values to appear at the beginning of the DataFrame instead of the end.

    Taylor: I wonder if this needs to be a nested function to deal with namespace stuff.
    """
    index_diff = len(left_df.index) - len(right_df.index)
    if index_diff > 0:
        # right_df is shorter
        config.loggers.utils.warning(f"Mismatch A between {left_df} and {right_df}")
        right_df = right_df.set_index(
            range(index_diff, len(right_df.index) + index_diff),
            drop=True,
            inplace=False,
        )
    elif index_diff < 0:
        # left_df is shorter
        left_df = left_df.set_index(
            range(-index_diff, len(left_df.index) - index_diff),
            drop=True,
            inplace=False,
        )

    return left_df, right_df


def _gather_confounds(
    signals=None,
    dvars=None,
    std_dvars=None,
    fdisp=None,
    rmsd=None,
    motion=None,
    newpath=None,
):
    """Load confounds from the filenames, concatenate together horizontally, and save new file.

    For some confounds (e.g., FD), the number of rows in the file will be one less than the
    number of volumes. This will be adjusted automatically in this function.
    """
    all_files = []
    confounds_list = []
    for confound, name in (
        (signals, "Global signals"),
        (std_dvars, "Standardized DVARS"),
        (dvars, "DVARS"),
        (fdisp, "Framewise displacement"),
        (rmsd, "RMSD"),
        (motion, "Motion parameters"),
    ):
        if confound is not None and isdefined(confound):
            confounds_list.append(name)
            if os.path.exists(confound) and os.stat(confound).st_size > 0:
                all_files.append(confound)

    confounds_data = pd.DataFrame()
    for file_name in all_files:  # assumes they all have headings already
        new = pd.read_csv(file_name, sep="\t")
        for column_name in new.columns:
            new.rename(
                columns={column_name: _camel_to_snake(_less_breakable(column_name))}, inplace=True
            )

        confounds_data, new = _adjust_indices(confounds_data, new)
        confounds_data = pd.concat((confounds_data, new), axis=1)

    if newpath is None:
        newpath = os.getcwd()

    combined_out = os.path.join(newpath, "confounds.tsv")
    confounds_data.to_csv(combined_out, sep="\t", index=False, na_rep="n/a")

    return combined_out, confounds_list
