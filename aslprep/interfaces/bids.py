"""Adapted interfaces from Niworkflows."""
from niworkflows.interfaces.bids import DerivativesDataSink as BaseDerivativesDataSink


class DerivativesDataSink(BaseDerivativesDataSink):
    """Store derivative files.

    A child class of the niworkflows DerivativesDataSink, using aslprep's configuration files.
    """

    out_path_base = "aslprep"