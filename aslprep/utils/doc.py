# -*- coding: utf-8 -*-
"""Functions related to the documentation.

docdict contains the standard documentation entries used across xcp_d.

source: Eric Larson and MNE-python team.
https://github.com/mne-tools/mne-python/blob/main/mne/utils/docs.py
"""
import sys

###################################
# Standard documentation entries
#
docdict = dict()

docdict[
    "omp_nthreads"
] = """
omp_nthreads : :obj:`int`
    Maximum number of threads an individual process may use.
"""

docdict[
    "mem_gb"
] = """
mem_gb : :obj:`float`
    Memory limit, in gigabytes.
"""

docdict[
    "output_dir"
] = """
output_dir : :obj:`str`
    Path to the output directory for ``aslprep`` derivatives.
    This should not include the ``aslprep`` folder.
    For example, "/path/to/dset/derivatives/".
"""

docdict[
    "work_dir"
] = """
work_dir : :obj:`str`
    Directory in which to store workflow execution state and temporary files.
"""

docdict[
    "analysis_level"
] = """
analysis_level : {"participant"}
    The analysis level for ``aslprep``. Must be specified as "participant" since ASLPrep
    performs analyses at the participant level.
"""

docdict[
    "basil"
] = """
basil : :obj:`bool`
    Run BASIL, FSL utils to compute CBF with spatial regularization and partial volume correction.
    BASIL will not be run if the ASL file only contains pre-calculated CBF images.
"""

docdict[
    "scorescrub"
] = """
scorescrub : :obj:`bool`
    Run SCORE and SCRUB, Sudipto Dolui's algorithms for denoising CBF.
    SCORE and SCRUB will not be run if the ASL file is short (i.e., if the GE workflow is used).
"""

docdict[
    "m0_scale"
] = """
m0_scale : :obj:`float`, optional
    Relative scale between ASL (delta-M) and M0 volumes.
    The M0 volumes will be multiplied by ``m0_scale`` when CBF is calculated.
    The default is 1 (no scaling).
"""

docdict[
    "smooth_kernel"
] = """
smooth_kernel : :obj:`float`
    Kernel size for smoothing M0.
"""

docdict[
    "processing_target"
] = """
processing_target : {"controllabel", "deltam", "cbf"}
    The target image types from the ASL file to process.
"""

docdict[
    "dummy_vols"
] = """
dummy_vols : :obj:`int`
    Number of label-control volume pairs to delete before CBF computation.
"""

docdict[
    "name"
] = """
name : :obj:`str`, optional
    Name of the workflow. This is used for working directories and workflow graphs.
"""

docdict[
    "aslcontext"
] = """
aslcontext : :obj:`str`
    Path to the ASL context file.
"""

docdict[
    "name_source"
] = """
name_source : :obj:`str`
    Path to the raw ASL file. Used as the base name for derivatives.
"""

docdict_indented = {}


def _indentcount_lines(lines):
    """Minimum indent for all lines in line list.

    >>> lines = [' one', '  two', '   three']
    >>> _indentcount_lines(lines)
    1
    >>> lines = []
    >>> _indentcount_lines(lines)
    0
    >>> lines = [' one']
    >>> _indentcount_lines(lines)
    1
    >>> _indentcount_lines(['    '])
    0

    """
    indentno = sys.maxsize
    for line in lines:
        stripped = line.lstrip()
        if stripped:
            indentno = min(indentno, len(line) - len(stripped))
    if indentno == sys.maxsize:
        return 0
    return indentno


def fill_doc(f):
    """Fill a docstring with docdict entries.

    Parameters
    ----------
    f : callable
        The function to fill the docstring of. Will be modified in place.

    Returns
    -------
    f : callable
        The function, potentially with an updated ``__doc__``.

    """
    docstring = f.__doc__
    if not docstring:
        return f
    lines = docstring.splitlines()
    # Find the minimum indent of the main docstring, after first line
    if len(lines) < 2:
        icount = 0
    else:
        icount = _indentcount_lines(lines[1:])
    # Insert this indent to dictionary docstrings
    try:
        indented = docdict_indented[icount]
    except KeyError:
        indent = " " * icount
        docdict_indented[icount] = indented = {}
        for name, dstr in docdict.items():
            lines = dstr.splitlines()
            try:
                newlines = [lines[0]]
                for line in lines[1:]:
                    newlines.append(indent + line)
                indented[name] = "\n".join(newlines)
            except IndexError:
                indented[name] = dstr
    try:
        f.__doc__ = docstring % indented
    except (TypeError, ValueError, KeyError) as exp:
        funcname = f.__name__
        funcname = docstring.split("\n")[0] if funcname is None else funcname
        raise RuntimeError(f"Error documenting {funcname}:\n{str(exp)}")
    return f
