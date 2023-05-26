"""Interfaces from FSL."""
from nipype.interfaces.base import File, TraitedSpec, traits
from nipype.interfaces.fsl import Split as FSLSplit


class _SplitOutputSpec(TraitedSpec):
    out_files = traits.List(File(exists=True))


class Split(FSLSplit):
    """Adapted version of FSLSplit interface that enforces a list output.

    In the official Nipype version of Split, out_files is an OutputMultiPath,
    which will return a string if the output is a single file.
    """

    output_spec = _SplitOutputSpec
