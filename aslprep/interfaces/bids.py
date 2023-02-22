"""Adapted interfaces from Niworkflows."""
from nipype import logging
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    OutputMultiObject,
    SimpleInterface,
    Str,
    TraitedSpec,
    traits,
)
from niworkflows.interfaces.bids import DerivativesDataSink as BaseDerivativesDataSink

LOGGER = logging.getLogger("nipype.interface")


class DerivativesDataSink(BaseDerivativesDataSink):
    """Store derivative files.

    A child class of the niworkflows DerivativesDataSink, using aslprep's output directory.
    """

    out_path_base = "aslprep"


class _BIDSDataGrabberInputSpec(BaseInterfaceInputSpec):
    subject_data = traits.Dict(Str, traits.Any)
    subject_id = Str()


class _BIDSDataGrabberOutputSpec(TraitedSpec):
    out_dict = traits.Dict(desc="output data structure")
    fmap = OutputMultiObject(desc="output fieldmaps")
    bold = OutputMultiObject(desc="output functional images")
    sbref = OutputMultiObject(desc="output sbrefs")
    t1w = OutputMultiObject(desc="output T1w images")
    roi = OutputMultiObject(desc="output ROI images")
    t2w = OutputMultiObject(desc="output T2w images")
    flair = OutputMultiObject(desc="output FLAIR images")
    asl = OutputMultiObject(desc="output ASL images")
    m0z = OutputMultiObject(desc="output MZeros images")
    cbf = OutputMultiObject(desc="output CBF images")


class BIDSDataGrabber(SimpleInterface):
    """
    Collect files from a BIDS directory structure.

    #>>> bids_src = BIDSDataGrabber(anat_only=False)
    #>>> bids_src.inputs.subject_data = bids_collect_data(
    #...     str(datadir / 'ds114'), '01', bids_validate=False)[0]
    #>>> bids_src.inputs.subject_id = '01'
    #>>> res = bids_src.run()
    #>>> res.outputs.t1w  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    ['.../ds114/sub-01/ses-retest/anat/sub-01_ses-retest_T1w.nii.gz',
     '.../ds114/sub-01/ses-test/anat/sub-01_ses-test_T1w.nii.gz']

    """

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
                f"No functional images found for subject sub-{self.inputs.subject_id}"
            )

        for imtype in ["bold", "t2w", "flair", "fmap", "sbref", "roi", "asl", "m0z", "cbf"]:
            if not bids_dict[imtype]:
                LOGGER.info(f'No "{imtype}" images found for sub-{self.inputs.subject_id}')

        return runtime
