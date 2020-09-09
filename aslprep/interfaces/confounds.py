# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Handling confounds.

    .. testsetup::

    >>> import os
    >>> import pandas as pd

"""
import os
import re
import shutil
import numpy as np
import pandas as pd
from nipype import logging
from nipype.utils.filemanip import fname_presuffix
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec, File, Directory, isdefined,
    SimpleInterface
)

LOGGER = logging.getLogger('nipype.interface')


class GatherConfoundsInputSpec(BaseInterfaceInputSpec):
    signals = File(exists=True, desc='input signals')
    dvars = File(exists=True, desc='file containing DVARS')
    std_dvars = File(exists=True, desc='file containing standardized DVARS')
    fd = File(exists=True, desc='input framewise displacement')
    rmsd = File(exists=True, desc='input RMS framewise displacement')
    motion = File(exists=True, desc='input motion parameters')


class GatherConfoundsOutputSpec(TraitedSpec):
    confounds_file = File(exists=True, desc='output confounds file')
    confounds_list = traits.List(traits.Str, desc='list of headers')


class GatherConfounds(SimpleInterface):
    r"""
    Combine various sources of confounds in one TSV file

    .. testsetup::

    >>> from tempfile import TemporaryDirectory
    >>> tmpdir = TemporaryDirectory()
    >>> os.chdir(tmpdir.name)

    .. doctest::

    >>> pd.DataFrame({'a': [0.1]}).to_csv('signals.tsv', index=False, na_rep='n/a')
    >>> pd.DataFrame({'b': [0.2]}).to_csv('dvars.tsv', index=False, na_rep='n/a')

    >>> gather = GatherConfounds()
    >>> gather.inputs.signals = 'signals.tsv'
    >>> gather.inputs.dvars = 'dvars.tsv'
    >>> res = gather.run()
    >>> res.outputs.confounds_list
    ['Global signals', 'DVARS']

    >>> pd.read_csv(res.outputs.confounds_file, sep='\s+', index_col=None,
    ...             engine='python')  # doctest: +NORMALIZE_WHITESPACE
         a    b
    0  0.1  0.2

    .. testcleanup::

    >>> tmpdir.cleanup()

    """
    input_spec = GatherConfoundsInputSpec
    output_spec = GatherConfoundsOutputSpec

    def _run_interface(self, runtime):
        combined_out, confounds_list = _gather_confounds(
            signals=self.inputs.signals,
            dvars=self.inputs.dvars,
            std_dvars=self.inputs.std_dvars,
            fdisp=self.inputs.fd,
            rmsd=self.inputs.rmsd,
            motion=self.inputs.motion,
            newpath=runtime.cwd,
        )
        self._results['confounds_file'] = combined_out
        self._results['confounds_list'] = confounds_list
        return runtime




def _gather_confounds(signals=None, dvars=None, std_dvars=None, fdisp=None,
                      rmsd=None, motion=None, newpath=None):
    r"""
    Load confounds from the filenames, concatenate together horizontally
    and save new file.

    >>> from tempfile import TemporaryDirectory
    >>> tmpdir = TemporaryDirectory()
    >>> os.chdir(tmpdir.name)
    >>> pd.DataFrame({'Global Signal': [0.1]}).to_csv('signals.tsv', index=False, na_rep='n/a')
    >>> pd.DataFrame({'stdDVARS': [0.2]}).to_csv('dvars.tsv', index=False, na_rep='n/a')
    >>> out_file, confound_list = _gather_confounds('signals.tsv', 'dvars.tsv')
    >>> confound_list
    ['Global signals', 'DVARS']

    >>> pd.read_csv(out_file, sep='\s+', index_col=None,
    ...             engine='python')  # doctest: +NORMALIZE_WHITESPACE
       global_signal  std_dvars
    0            0.1        0.2
    >>> tmpdir.cleanup()


    """

    def less_breakable(a_string):
        ''' hardens the string to different envs (i.e., case insensitive, no whitespace, '#' '''
        return ''.join(a_string.split()).strip('#')

    # Taken from https://stackoverflow.com/questions/1175208/
    # If we end up using it more than just here, probably worth pulling in a well-tested package
    def camel_to_snake(name):
        s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
        return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()

    def _adjust_indices(left_df, right_df):
        # This forces missing values to appear at the beggining of the DataFrame
        # instead of the end
        index_diff = len(left_df.index) - len(right_df.index)
        if index_diff > 0:
            right_df.index = range(index_diff,
                                   len(right_df.index) + index_diff)
        elif index_diff < 0:
            left_df.index = range(-index_diff,
                                  len(left_df.index) - index_diff)

    all_files = []
    confounds_list = []
    for confound, name in ((signals, 'Global signals'),
                           (std_dvars, 'Standardized DVARS'),
                           (dvars, 'DVARS'),
                           (fdisp, 'Framewise displacement'),
                           (rmsd, 'Framewise displacement (RMS)'),
                           (motion, 'Motion parameters')):
        if confound is not None and isdefined(confound):
            confounds_list.append(name)
            if os.path.exists(confound) and os.stat(confound).st_size > 0:
                all_files.append(confound)

    confounds_data = pd.DataFrame()
    for file_name in all_files:  # assumes they all have headings already
        new = pd.read_csv(file_name, sep="\t")
        for column_name in new.columns:
            new.rename(columns={column_name: camel_to_snake(less_breakable(column_name))},
                       inplace=True)

        _adjust_indices(confounds_data, new)
        confounds_data = pd.concat((confounds_data, new), axis=1)

    if newpath is None:
        newpath = os.getcwd()

    combined_out = os.path.join(newpath, 'confounds.tsv')
    confounds_data.to_csv(combined_out, sep='\t', index=False,
                          na_rep='n/a')

    return combined_out, confounds_list


class ASLSummaryInputSpec(BaseInterfaceInputSpec):
    in_func = File(exists=True, mandatory=True,
                   desc='input ASL time-series (4D file)')
    in_mask = File(exists=True,
                   desc='3D brain mask')
    in_segm = File(exists=True, desc='resampled segmentation')
    confounds_file = File(exists=True,
                          desc="BIDS' _confounds.tsv file")

    str_or_tuple = traits.Either(
        traits.Str,
        traits.Tuple(traits.Str, traits.Either(None, traits.Str)),
        traits.Tuple(traits.Str, traits.Either(None, traits.Str), traits.Either(None, traits.Str)))
    confounds_list = traits.List(
        str_or_tuple, minlen=1,
        desc='list of headers to extract from the confounds_file')
    tr = traits.Either(None, traits.Float, usedefault=True,
                       desc='the repetition time')


class ASLSummaryOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='written file path')


class ASLSummary(SimpleInterface):
    """
    Copy the x-form matrices from `hdr_file` to `out_file`.
    """
    input_spec = ASLSummaryInputSpec
    output_spec = ASLSummaryOutputSpec

    def _run_interface(self, runtime):
        from ..niworkflows.viz.plots import ASLPlot

        self._results['out_file'] = fname_presuffix(
            self.inputs.in_func,
            suffix='_aslplot.svg',
            use_ext=False,
            newpath=runtime.cwd)

        dataframe = pd.read_csv(
            self.inputs.confounds_file,
            sep="\t", index_col=None, dtype='float32',
            na_filter=True, na_values='n/a')

        headers = []
        units = {}
        names = {}

        for conf_el in self.inputs.confounds_list:
            if isinstance(conf_el, (list, tuple)):
                headers.append(conf_el[0])
                if conf_el[1] is not None:
                    units[conf_el[0]] = conf_el[1]

                if len(conf_el) > 2 and conf_el[2] is not None:
                    names[conf_el[0]] = conf_el[2]
            else:
                headers.append(conf_el)

        if not headers:
            data = None
            units = None
        else:
            data = dataframe[headers]

        colnames = data.columns.ravel().tolist()

        for name, newname in list(names.items()):
            colnames[colnames.index(name)] = newname

        data.columns = colnames

        fig = ASLPlot(
            self.inputs.in_func,
            mask_file=self.inputs.in_mask if isdefined(self.inputs.in_mask) else None,
            seg_file=(self.inputs.in_segm
                      if isdefined(self.inputs.in_segm) else None),
            tr=self.inputs.tr,
            data=data,
            units=units,
        ).plot()
        fig.savefig(self._results['out_file'], bbox_inches='tight')
        return runtime
