"""Adapted interfaces from Niworkflows."""
from aslprep.niworkflows.interfaces.bids import (
    DerivativesDataSink as BaseDerivativesDataSink,
)


class DerivativesDataSink(BaseDerivativesDataSink):
    """Store derivative files.

    A child class of the niworkflows DerivativesDataSink, using aslprep's output directory.
    """

    out_path_base = "aslprep"
