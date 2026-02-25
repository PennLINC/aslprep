"""Adapted interfaces from Niworkflows."""

import os
import shutil
from json import dump, loads

import filelock
import nibabel as nb
import numpy as np
from bids.layout import Config
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    Directory,
    DynamicTraitedSpec,
    File,
    OutputMultiObject,
    SimpleInterface,
    Str,
    TraitedSpec,
    traits,
)
from nipype.interfaces.io import add_traits
from niworkflows.interfaces.bids import DerivativesDataSink as BaseDerivativesDataSink

from aslprep import config
from aslprep.data import load as load_data
from aslprep.utils.bids import _get_bidsuris, get_entity

# NOTE: Modified for aslprep's purposes
aslprep_spec = loads(load_data.readable('aslprep_bids_config.json').read_text())
bids_config = Config.load('bids')
deriv_config = Config.load('derivatives')

aslprep_entities = {v['name']: v['pattern'] for v in aslprep_spec['entities']}
merged_entities = {**bids_config.entities, **deriv_config.entities}
merged_entities = {k: v.pattern for k, v in merged_entities.items()}
merged_entities = {**merged_entities, **aslprep_entities}
merged_entities = [{'name': k, 'pattern': v} for k, v in merged_entities.items()]
config_entities = frozenset({e['name'] for e in merged_entities})


class _BIDSDataGrabberInputSpec(BaseInterfaceInputSpec):
    subject_data = traits.Dict(Str, traits.Any)
    subject_id = Str()


class _BIDSDataGrabberOutputSpec(TraitedSpec):
    out_dict = traits.Dict(desc='output data structure')
    fmap = OutputMultiObject(desc='output fieldmaps')
    bold = OutputMultiObject(desc='output functional images')
    sbref = OutputMultiObject(desc='output sbrefs')
    t1w = OutputMultiObject(desc='output T1w images')
    roi = OutputMultiObject(desc='output ROI images')
    t2w = OutputMultiObject(desc='output T2w images')
    flair = OutputMultiObject(desc='output FLAIR images')
    pet = OutputMultiObject(desc='output PET images')
    dwi = OutputMultiObject(desc='output DWI images')
    asl = OutputMultiObject(desc='output ASL images')


class BIDSDataGrabber(SimpleInterface):
    """Collect files from a BIDS directory structure.

    .. testsetup::

        >>> data_dir_canary()

    >>> bids_src = BIDSDataGrabber(anat_only=False)
    >>> bids_src.inputs.subject_data = bids_collect_data(
    ...     str(datadir / 'ds114'), '01', bids_validate=False)[0]
    >>> bids_src.inputs.subject_id = '01'
    >>> res = bids_src.run()
    >>> res.outputs.t1w  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    ['.../ds114/sub-01/ses-retest/anat/sub-01_ses-retest_T1w.nii.gz',
     '.../ds114/sub-01/ses-test/anat/sub-01_ses-test_T1w.nii.gz']

    """

    input_spec = _BIDSDataGrabberInputSpec
    output_spec = _BIDSDataGrabberOutputSpec
    _require_funcs = True

    def __init__(self, *args, **kwargs):
        anat_only = kwargs.pop('anat_only')
        anat_derivatives = kwargs.pop('anat_derivatives', None)
        super().__init__(*args, **kwargs)
        if anat_only is not None:
            self._require_funcs = not anat_only
        self._require_t1w = anat_derivatives is None

    def _run_interface(self, runtime):
        bids_dict = self.inputs.subject_data

        self._results['out_dict'] = bids_dict
        self._results.update(bids_dict)

        if self._require_t1w and not bids_dict['t1w']:
            raise FileNotFoundError(
                f'No T1w images found for subject sub-{self.inputs.subject_id}'
            )

        if self._require_funcs and not bids_dict['asl']:
            raise FileNotFoundError(
                f'No ASL images found for subject sub-{self.inputs.subject_id}'
            )

        for imtype in ['bold', 't2w', 'flair', 'fmap', 'sbref', 'roi', 'pet', 'asl']:
            if not bids_dict.get(imtype):
                config.loggers.interface.info(
                    'No "%s" images found for sub-%s', imtype, self.inputs.subject_id
                )

        return runtime


class DerivativesDataSink(BaseDerivativesDataSink):
    """Store derivative files.

    A child class of the niworkflows DerivativesDataSink, using aslprep's configuration files.
    """

    out_path_base = ''
    _allowed_entities = set(config_entities)
    _config_entities = config_entities
    _config_entities_dict = merged_entities
    _file_patterns = aslprep_spec['default_path_patterns']


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


class FunctionOverrideContext:
    """Override a function in imported code with a context manager.

    Even though this class is *currently* unused, I'm keeping it around for when I need to override
    prepare_timing_parameters once fMRIPrep's init_bold_surf_wf is usable
    (i.e., once the DerivativesDataSink import is moved out of the function).

    Here's how it worked before:

    def _fake_params(metadata):  # noqa: U100
        return {"SliceTimingCorrected": False}

    # init_bold_surf_wf uses prepare_timing_parameters, which uses the config object.
    # The uninitialized fMRIPrep config will have config.workflow.ignore set to None
    # instead of a list, which will raise an error.
    with FunctionOverrideContext(resampling, "prepare_timing_parameters", _fake_params):
        asl_surf_wf = resampling.init_bold_surf_wf(
            mem_gb=mem_gb["resampled"],
            metadata=metadata,
            surface_spaces=freesurfer_spaces,
            medial_surface_nan=config.workflow.medial_surface_nan,
            output_dir=config.execution.aslprep_dir,
            name="asl_surf_wf",
        )
    """

    def __init__(self, module, function_name, new_function):
        self.module = module
        self.function_name = function_name
        self.new_function = new_function
        self.original_function = None

    def __enter__(self):
        """Enter the context manager and perform the function override.

        Returns
        -------
        FunctionOverrideContext
            The instance of the context manager.
        """
        self.original_function = getattr(self.module, self.function_name)
        setattr(self.module, self.function_name, self.new_function)

    def __exit__(self, exc_type, exc_value, traceback):  # noqa: U100
        """Exit the context manager and restore the original function definition.

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
        setattr(self.module, self.function_name, self.original_function)


class _BIDSURIInputSpec(DynamicTraitedSpec):
    dataset_links = traits.Dict(mandatory=True, desc='Dataset links')
    out_dir = traits.Str(mandatory=True, desc='Output directory')
    metadata = traits.Dict(desc='Metadata dictionary')
    field = traits.Str(
        'Sources',
        usedefault=True,
        desc='Field to use for BIDS URIs in metadata dict',
    )


class _BIDSURIOutputSpec(TraitedSpec):
    out = traits.List(
        traits.Str,
        desc='BIDS URI(s) for file',
    )
    metadata = traits.Dict(
        desc="Dictionary with 'Sources' field.",
    )


class BIDSURI(SimpleInterface):
    """Convert input filenames to BIDS URIs, based on links in the dataset.

    This interface can combine multiple lists of inputs.
    """

    input_spec = _BIDSURIInputSpec
    output_spec = _BIDSURIOutputSpec

    def __init__(self, numinputs=0, **inputs):
        super().__init__(**inputs)
        self._numinputs = numinputs
        if numinputs >= 1:
            input_names = [f'in{i + 1}' for i in range(numinputs)]
        else:
            input_names = []
        add_traits(self.inputs, input_names)

    def _run_interface(self, runtime):
        inputs = [getattr(self.inputs, f'in{i + 1}') for i in range(self._numinputs)]
        uris = _get_bidsuris(inputs, self.inputs.dataset_links, self.inputs.out_dir)
        self._results['out'] = uris

        # Add the URIs to the metadata dictionary.
        metadata = self.inputs.metadata or {}
        metadata = metadata.copy()
        metadata[self.inputs.field] = metadata.get(self.inputs.field, []) + uris
        self._results['metadata'] = metadata

        return runtime


class _CopyAtlasInputSpec(BaseInterfaceInputSpec):
    name_source = traits.Str(
        desc="The source file's name.",
        mandatory=False,
    )
    in_file = File(
        exists=True,
        desc='The atlas file to copy.',
        mandatory=True,
    )
    meta_dict = traits.Either(
        traits.Dict(),
        None,
        desc='The atlas metadata dictionary.',
        mandatory=False,
    )
    output_dir = Directory(
        exists=True,
        desc='The output directory.',
        mandatory=True,
    )
    atlas = traits.Str(
        desc='The atlas name.',
        mandatory=True,
    )
    Sources = traits.List(
        traits.Str,
        desc='List of sources for the atlas.',
        mandatory=False,
    )


class _CopyAtlasOutputSpec(TraitedSpec):
    out_file = File(
        exists=True,
        desc='The copied atlas file.',
    )


class CopyAtlas(SimpleInterface):
    """Copy atlas file to output directory.

    Parameters
    ----------
    name_source : :obj:`str`
        The source name of the atlas file.
    in_file : :obj:`str`
        The atlas file to copy.
    output_dir : :obj:`str`
        The output directory.
    atlas : :obj:`str`
        The name of the atlas.

    Returns
    -------
    out_file : :obj:`str`
        The path to the copied atlas file.

    Notes
    -----
    I can't use DerivativesDataSink because it has a problem with dlabel CIFTI files.
    It gives the following error:
    "AttributeError: 'Cifti2Header' object has no attribute 'set_data_dtype'"

    I can't override the CIFTI atlas's data dtype ahead of time because setting it to int8 or int16
    somehow converts all of the values in the data array to weird floats.
    This could be a version-specific nibabel issue.

    I've also updated this function to handle JSON and TSV files as well.
    """

    input_spec = _CopyAtlasInputSpec
    output_spec = _CopyAtlasOutputSpec

    def _run_interface(self, runtime):
        output_dir = self.inputs.output_dir
        in_file = self.inputs.in_file
        meta_dict = self.inputs.meta_dict
        name_source = self.inputs.name_source
        atlas = self.inputs.atlas
        Sources = self.inputs.Sources

        tpl = get_entity(name_source, 'tpl')
        if not tpl:
            tpl = get_entity(name_source, 'space')

        if not tpl:
            raise ValueError(f'Could not determine template from {name_source}')

        cohort, cohort_str = None, ''
        if '+' in tpl:
            # Split the template and cohort
            tpl, cohort = tpl.split('+')
            cohort_str = f'cohort-{cohort}_'

        if not cohort:
            cohort = get_entity(name_source, 'cohort')
            cohort_str = f'cohort-{cohort}_' if cohort else ''

        tpl_str = f'tpl-{tpl}_'

        output_dir = os.path.join(output_dir, f'tpl-{tpl}')
        if cohort:
            output_dir = os.path.join(output_dir, f'cohort-{cohort}')

        res = get_entity(name_source, 'res')
        res_str = f'_res-{res}' if res else ''

        den = get_entity(name_source, 'den')
        den_str = f'_den-{den}' if den else ''

        out_basename = f'{tpl_str}{cohort_str}atlas-{atlas}{res_str}{den_str}_dseg'
        if in_file.endswith('.tsv'):
            extension = '.tsv'
        elif in_file.endswith('.dlabel.nii'):
            extension = '.dlabel.nii'
        else:
            extension = '.nii.gz'

        os.makedirs(output_dir, exist_ok=True)
        out_file = os.path.join(output_dir, f'{out_basename}{extension}')

        if out_file.endswith('.nii.gz') and os.path.isfile(out_file):
            # Check that native-resolution atlas doesn't have a different resolution from the last
            # run's atlas.
            old_img = nb.load(out_file)
            new_img = nb.load(in_file)
            if not np.allclose(old_img.affine, new_img.affine):
                raise ValueError(
                    f"Existing '{atlas}' atlas affine ({out_file}) is different from the input "
                    f'file affine ({in_file}).'
                )

        # Don't copy the file if it exists, to prevent any race conditions between parallel
        # processes.
        if not os.path.isfile(out_file):
            lock_file = os.path.join(output_dir, f'{out_basename}{extension}.lock')
            with filelock.SoftFileLock(lock_file, timeout=60):
                shutil.copyfile(in_file, out_file)

        # Only write out a sidecar if metadata are provided
        if meta_dict or Sources:
            meta_file = os.path.join(output_dir, f'{out_basename}.json')
            lock_meta = os.path.join(output_dir, f'{out_basename}.json.lock')
            meta_dict = meta_dict or {}
            meta_dict = meta_dict.copy()
            if Sources:
                meta_dict['Sources'] = meta_dict.get('Sources', []) + Sources

            with filelock.SoftFileLock(lock_meta, timeout=60):
                with open(meta_file, 'w') as fobj:
                    dump(meta_dict, fobj, sort_keys=True, indent=4)

        self._results['out_file'] = out_file

        return runtime


class _CopyAtlasDescriptionInputSpec(BaseInterfaceInputSpec):
    in_dir = Directory(
        exists=True,
        desc='The atlas directory to copy the description file from.',
        mandatory=True,
    )
    atlas_name = traits.Str(
        desc='The name of the atlas.',
        mandatory=True,
    )
    output_dir = Directory(
        exists=True,
        desc='The output directory.',
        mandatory=True,
    )


class _CopyAtlasDescriptionOutputSpec(TraitedSpec):
    out_file = File(
        exists=True,
        desc='The copied atlas file.',
    )


class CopyAtlasDescription(SimpleInterface):
    """Copy atlas description file to output directory.

    Parameters
    ----------
    in_file : :obj:`str`
        The atlas file to copy.
    output_dir : :obj:`str`
        The output directory.

    Returns
    -------
    out_file : :obj:`str`
        The path to the copied atlas description file.
    """

    input_spec = _CopyAtlasDescriptionInputSpec
    output_spec = _CopyAtlasDescriptionOutputSpec

    def _run_interface(self, runtime):
        output_dir = self.inputs.output_dir
        in_dir = self.inputs.in_dir
        atlas_name = self.inputs.atlas_name

        os.makedirs(output_dir, exist_ok=True)

        in_file = os.path.join(in_dir, f'atlas-{atlas_name}_description.json')
        if not os.path.isfile(in_file):
            raise FileNotFoundError(f'Atlas description file not found: {in_file}')

        out_file = os.path.join(output_dir, f'atlas-{atlas_name}_description.json')
        if not os.path.isfile(out_file):
            shutil.copyfile(in_file, out_file)

        self._results['out_file'] = out_file

        return runtime
