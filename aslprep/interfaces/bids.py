"""Adapted interfaces from Niworkflows."""
from json import loads
from pathlib import Path

from bids.layout import Config
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    OutputMultiObject,
    SimpleInterface,
    Str,
    TraitedSpec,
    traits,
)
from niworkflows.interfaces.bids import DerivativesDataSink as BaseDerivativesDataSink
from pkg_resources import resource_filename as pkgrf

from aslprep import config

# NOTE: Modified for aslprep's purposes
aslprep_spec = loads(Path(pkgrf("aslprep", "data/aslprep_bids_config.json")).read_text())
bids_config = Config.load("bids")
deriv_config = Config.load("derivatives")

aslprep_entities = {v["name"]: v["pattern"] for v in aslprep_spec["entities"]}
merged_entities = {**bids_config.entities, **deriv_config.entities}
merged_entities = {k: v.pattern for k, v in merged_entities.items()}
merged_entities = {**merged_entities, **aslprep_entities}
merged_entities = [{"name": k, "pattern": v} for k, v in merged_entities.items()]
config_entities = frozenset({e["name"] for e in merged_entities})


class _BIDSDataGrabberInputSpec(BaseInterfaceInputSpec):
    subject_data = traits.Dict(Str, traits.Any)
    subject_id = Str()


class _BIDSDataGrabberOutputSpec(TraitedSpec):
    out_dict = traits.Dict(desc="output data structure")
    fmap = OutputMultiObject(desc="output fieldmaps")
    asl = OutputMultiObject(desc="output ASL images")
    sbref = OutputMultiObject(desc="output sbrefs")
    t1w = OutputMultiObject(desc="output T1w images")
    roi = OutputMultiObject(desc="output ROI images")
    t2w = OutputMultiObject(desc="output T2w images")
    flair = OutputMultiObject(desc="output FLAIR images")


class BIDSDataGrabber(SimpleInterface):
    """Collect files from a BIDS directory structure."""

    input_spec = _BIDSDataGrabberInputSpec
    output_spec = _BIDSDataGrabberOutputSpec
    _require_funcs = True

    def __init__(self, *args, **kwargs):
        anat_only = kwargs.pop("anat_only")
        super(BIDSDataGrabber, self).__init__(*args, **kwargs)
        if anat_only is not None:
            self._require_funcs = not anat_only

    def _run_interface(self, runtime):
        bids_dict = self.inputs.subject_data

        self._results["out_dict"] = bids_dict
        self._results.update(bids_dict)

        if not bids_dict["t1w"]:
            raise FileNotFoundError(
                f"No T1w images found for subject sub-{self.inputs.subject_id}"
            )

        if self._require_funcs and not bids_dict["asl"]:
            raise FileNotFoundError(
                f"No ASL images found for subject sub-{self.inputs.subject_id}"
            )

        for imtype in ["t2w", "flair", "fmap", "sbref", "roi", "asl"]:
            if not bids_dict[imtype]:
                config.loggers.interface.info(
                    'No "%s" images found for sub-%s',
                    imtype,
                    self.inputs.subject_id,
                )

        return runtime


class DerivativesDataSink(BaseDerivativesDataSink):
    """Store derivative files.

    A child class of the niworkflows DerivativesDataSink, using aslprep's configuration files.
    """

    out_path_base = ""
    _allowed_entities = set(config_entities)
    _config_entities = config_entities
    _config_entities_dict = merged_entities
    _file_patterns = aslprep_spec["default_path_patterns"]


class OverrideDerivativesDataSink:
    """A context manager for temporarily overriding the definition of DerivativesDataSink.

    Parameters
    ----------
    None

    Attributes
    ----------
    original_class (type): The original class that is replaced during the override.

    Methods
    -------
    __init__()
        Initialize the context manager.
    __enter__()
        Enters the context manager and performs the class override.
    __exit__(exc_type, exc_value, traceback)
        Exits the context manager and restores the original class definition.
    """

    def __init__(self, module):
        """Initialize the context manager with the target module.

        Parameters
        -----------
        module
            The module where SomeClass should be overridden.
        """
        self.module = module

    def __enter__(self):
        """Enter the context manager and perform the class override.

        Returns
        -------
        OverrideConfoundsDerivativesDataSink
            The instance of the context manager.
        """
        # Save the original class
        self.original_class = self.module.DerivativesDataSink
        # Replace SomeClass with YourOwnClass
        self.module.DerivativesDataSink = DerivativesDataSink
        return self

    def __exit__(self, exc_type, exc_value, traceback):  # noqa: U100
        """Exit the context manager and restore the original class definition.

        Parameters
        ----------
        exc_type : type
            The type of the exception (if an exception occurred).
        exc_value : Exception
            The exception instance (if an exception occurred).
        traceback : traceback
            The traceback information (if an exception occurred).

        Returns
        -------
        None
        """
        # Restore the original class
        self.module.DerivativesDataSink = self.original_class
