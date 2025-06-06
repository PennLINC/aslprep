"""Utility interfaces for ASLPrep."""

import os

import nibabel as nb
import pandas as pd
from nilearn import image
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec
from nipype.interfaces.nilearn import NilearnBaseInterface
from nipype.utils.filemanip import fname_presuffix, load_json, save_json

from aslprep.utils.asl import reduce_metadata_lists


class _ReduceASLFilesInputSpec(BaseInterfaceInputSpec):
    asl_file = File(exists=True, mandatory=True, desc='ASL file to split.')
    aslcontext = File(exists=True, mandatory=True, desc='aslcontext TSV.')
    processing_target = traits.Str()
    metadata = traits.Dict()


class _ReduceASLFilesOutputSpec(TraitedSpec):
    asl_file = File(exists=True, desc='Modified ASL file.')
    aslcontext = File(exists=True, desc='Modified aslcontext file.')
    metadata = traits.Dict()


class ReduceASLFiles(SimpleInterface):
    """Split ASL data into different files."""

    input_spec = _ReduceASLFilesInputSpec
    output_spec = _ReduceASLFilesOutputSpec

    def _run_interface(self, runtime):
        aslcontext = pd.read_table(self.inputs.aslcontext)
        asl_img = nb.load(self.inputs.asl_file)
        if asl_img.shape[3] != aslcontext.shape[0]:
            raise ValueError(
                f"Number of volumes in {self.inputs.asl_file} ({asl_img.shape[3]}) doesn't equal "
                f'number of rows in {self.inputs.aslcontext} ({aslcontext.shape[0]}).'
            )

        if self.inputs.processing_target == 'control':
            files_to_keep = ['control', 'label', 'm0scan']
        elif self.inputs.processing_target == 'deltam':
            files_to_keep = ['deltam', 'm0scan']
        else:
            files_to_keep = ['cbf', 'm0scan']

        n_volumes = aslcontext.shape[0]
        asl_idx = aslcontext.loc[aslcontext['volume_type'].isin(files_to_keep)].index.values
        asl_idx = asl_idx.astype(int)
        self._results['metadata'] = reduce_metadata_lists(
            metadata=self.inputs.metadata,
            n_volumes=n_volumes,
            keep_idx=asl_idx,
        )

        asl_img = image.index_img(asl_img, asl_idx)

        self._results['asl_file'] = fname_presuffix(
            self.inputs.asl_file,
            suffix='_reduced',
            newpath=runtime.cwd,
            use_ext=True,
        )
        asl_img.to_filename(self._results['asl_file'])

        aslcontext = aslcontext.loc[asl_idx]
        self._results['aslcontext'] = fname_presuffix(
            self.inputs.aslcontext,
            suffix='_reduced',
            newpath=runtime.cwd,
            use_ext=True,
        )
        aslcontext.to_csv(self._results['aslcontext'], sep='\t', index=False)

        return runtime


class _RMSDiffInputSpec(FSLCommandInputSpec):
    matrixfile1 = File(
        exists=True,
        position=0,
        argstr='%s',
        desc='First matrix file.',
        mandatory=True,
    )
    matrixfile2 = File(
        exists=True,
        position=1,
        argstr='%s',
        desc='Second matrix file.',
        mandatory=True,
    )
    ref_vol = File(
        exists=True,
        position=2,
        argstr='%s',
        desc='Reference volume.',
        mandatory=True,
    )


class _RMSDiffOutputSpec(TraitedSpec):
    rmsd = traits.Float()


class RMSDiff(FSLCommand):
    """Run rmsdiff."""

    _cmd = 'rmsdiff'
    input_spec = _RMSDiffInputSpec
    output_spec = _RMSDiffOutputSpec

    def aggregate_outputs(self, runtime=None, needed_outputs=None):  # noqa: U100
        """Taken from nipype.interfaces.afni.preprocess.ClipLevel."""
        outputs = self._outputs()

        outfile = os.path.join(os.getcwd(), 'stat_result.json')

        if runtime is None:
            try:
                rmsd = load_json(outfile)['stat']
            except OSError:
                return self.run().outputs
        else:
            rmsd = []
            for line in runtime.stdout.split('\n'):
                if line:
                    values = line.split()
                    if len(values) > 1:
                        rmsd.append([float(val) for val in values])
                    else:
                        rmsd.extend([float(val) for val in values])

            if len(rmsd) == 1:
                rmsd = rmsd[0]

            save_json(outfile, {'stat': rmsd})

        outputs.rmsd = rmsd

        return outputs


class _PairwiseRMSDiffInputSpec(BaseInterfaceInputSpec):
    in_files = traits.List(
        File(exists=True),
        desc='Matrix files to compare with each other.',
        mandatory=True,
    )
    ref_file = File(
        exists=True,
        desc='Reference volume.',
        mandatory=True,
    )


class _PairwiseRMSDiffOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='Output txt file.')


class PairwiseRMSDiff(SimpleInterface):
    """Run rmsdiff on each contiguous pair of transform files to build a txt file of rmsd values.

    This interface uses :class:`~aslprep.interfaces.utility.RMSDiff` internally, which may not be
    a proper nipype pattern.
    """

    input_spec = _PairwiseRMSDiffInputSpec
    output_spec = _PairwiseRMSDiffOutputSpec

    def _run_interface(self, runtime):
        rmsd = []
        for i_file in range(len(self.inputs.in_files) - 1):
            j_file = i_file + 1
            file1 = self.inputs.in_files[i_file]
            file2 = self.inputs.in_files[j_file]
            # run rmsdiff
            rmsdiff = RMSDiff(matrixfile1=file1, matrixfile2=file2, ref_vol=self.inputs.ref_file)
            res = rmsdiff.run()
            if not isinstance(res.outputs.rmsd, float):
                raise ValueError(f'Expected a float, got {res.outputs.rmsd}')
            rmsd.append(str(res.outputs.rmsd))

        self._results['out_file'] = fname_presuffix(
            self.inputs.ref_file,
            suffix='_rmsd.txt',
            newpath=runtime.cwd,
            use_ext=False,
        )
        with open(self._results['out_file'], 'w') as fobj:
            fobj.write('\n'.join(rmsd))

        return runtime


class _CombineMotionParametersInputSpec(BaseInterfaceInputSpec):
    aslcontext = File(exists=True)
    volume_types = traits.List(traits.Str())
    mat_files = traits.List(traits.List(File(exists=True)))
    par_files = traits.List(File(exists=True))


class _CombineMotionParametersOutputSpec(TraitedSpec):
    mat_file_list = traits.List(File(exists=True))
    combined_par_file = File(exists=True)


class CombineMotionParameters(SimpleInterface):
    """Combine motion parameter files from MCFLIRT across image types."""

    input_spec = _CombineMotionParametersInputSpec
    output_spec = _CombineMotionParametersOutputSpec

    def _run_interface(self, runtime):
        aslcontext = pd.read_table(self.inputs.aslcontext)
        out_par = [None] * aslcontext.shape[0]
        out_mat_files = [None] * aslcontext.shape[0]

        if len(self.inputs.volume_types) != len(self.inputs.mat_files):
            raise ValueError('Number of volume types and number of mat files must be the same.')
        if len(self.inputs.volume_types) != len(self.inputs.par_files):
            raise ValueError('Number of volume types and number of par files must be the same.')

        for i_type, volume_type in enumerate(self.inputs.volume_types):
            type_mat_files = self.inputs.mat_files[i_type]
            type_par_file = self.inputs.par_files[i_type]

            type_idx = aslcontext.loc[aslcontext['volume_type'] == volume_type].index.values

            with open(type_par_file) as fobj:
                par = fobj.readlines()

            for i_vol, vol_idx in enumerate(type_idx):
                out_par[vol_idx] = par[i_vol]
                out_mat_files[vol_idx] = type_mat_files[i_vol]

        self._results['combined_par_file'] = fname_presuffix(
            type_par_file,
            suffix='_combined',
            newpath=runtime.cwd,
            use_ext=True,
        )
        with open(self._results['combined_par_file'], 'w') as fobj:
            fobj.write(''.join(out_par))

        self._results['mat_file_list'] = out_mat_files

        return runtime


class _SplitByVolumeTypeInputSpec(BaseInterfaceInputSpec):
    aslcontext = File(exists=True)
    asl_file = File(exists=True)


class _SplitByVolumeTypeOutputSpec(TraitedSpec):
    out_files = traits.List(File(exists=True))
    volume_types = traits.List(traits.Str())


class SplitByVolumeType(SimpleInterface):
    """Split out a specific volume type from the ASL file."""

    input_spec = _SplitByVolumeTypeInputSpec
    output_spec = _SplitByVolumeTypeOutputSpec

    def _run_interface(self, runtime):
        aslcontext = pd.read_table(self.inputs.aslcontext)
        volume_types = sorted(aslcontext['volume_type'].unique())
        out_files = []
        for volume_type in volume_types:
            volumetype_df = aslcontext.loc[aslcontext['volume_type'] == volume_type]
            volumetype_idx = volumetype_df.index.tolist()

            out_img = image.index_img(self.inputs.asl_file, volumetype_idx)
            out_file = fname_presuffix(
                self.inputs.asl_file,
                suffix=f'_{volume_type}',
                newpath=runtime.cwd,
                use_ext=True,
            )
            out_img.to_filename(out_file)
            out_files.append(out_file)

        self._results['out_files'] = out_files
        self._results['volume_types'] = volume_types

        return runtime


class _SmoothInputSpec(BaseInterfaceInputSpec):
    in_file = File(
        exists=True,
        mandatory=True,
        desc='An image to smooth.',
    )
    fwhm = traits.Either(
        traits.Float(),
        traits.List(
            traits.Float(),
            minlen=3,
            maxlen=3,
        ),
        desc=(
            'Full width at half maximum. '
            'Smoothing strength, as a full-width at half maximum, in millimeters.'
        ),
    )
    out_file = File(
        'smooth_img.nii.gz',
        usedefault=True,
        exists=False,
        desc='The name of the smoothed file to write out. smooth_img.nii.gz by default.',
    )


class _SmoothOutputSpec(TraitedSpec):
    out_file = File(
        exists=True,
        desc='Smoothed output file.',
    )


class Smooth(NilearnBaseInterface, SimpleInterface):
    """Smooth image."""

    input_spec = _SmoothInputSpec
    output_spec = _SmoothOutputSpec

    def _run_interface(self, runtime):
        from nilearn.image import smooth_img

        img_smoothed = smooth_img(self.inputs.in_file, fwhm=self.inputs.fwhm)
        self._results['out_file'] = fname_presuffix(
            self.inputs.in_file,
            suffix='_sm.nii.gz',
            newpath=runtime.cwd,
            use_ext=False,
        )
        img_smoothed.to_filename(self._results['out_file'])

        return runtime


class _IndexImageInputSpec(BaseInterfaceInputSpec):
    in_file = File(
        exists=True,
        mandatory=True,
        desc='A 4D image to index.',
    )
    index = traits.Int(
        0,
        usedefault=True,
        desc='Volume index to select from in_file.',
    )
    out_file = File(
        'img_3d.nii.gz',
        usedefault=True,
        exists=False,
        desc='The name of the indexed file.',
    )


class _IndexImageOutputSpec(TraitedSpec):
    out_file = File(
        exists=True,
        desc='Concatenated output file.',
    )


class IndexImage(NilearnBaseInterface, SimpleInterface):
    """Select a specific volume from a 4D image."""

    input_spec = _IndexImageInputSpec
    output_spec = _IndexImageOutputSpec

    def _run_interface(self, runtime):
        from nilearn.image import index_img

        img_3d = index_img(self.inputs.in_file, self.inputs.index)
        self._results['out_file'] = os.path.join(runtime.cwd, self.inputs.out_file)
        img_3d.to_filename(self._results['out_file'])

        return runtime


class _MeanImageInputSpec(BaseInterfaceInputSpec):
    in_file = File(
        exists=True,
        mandatory=True,
        desc='A 4D image to index.',
    )
    out_file = File(
        'img_3d.nii.gz',
        usedefault=True,
        exists=False,
        desc='The name of the indexed file.',
    )


class _MeanImageOutputSpec(TraitedSpec):
    out_file = File(
        exists=True,
        desc='Concatenated output file.',
    )


class MeanImage(NilearnBaseInterface, SimpleInterface):
    """Calculate the mean image of a 4D image."""

    input_spec = _MeanImageInputSpec
    output_spec = _MeanImageOutputSpec

    def _run_interface(self, runtime):
        from nilearn.image import mean_img

        img_mean = mean_img(self.inputs.in_file)
        self._results['out_file'] = os.path.join(runtime.cwd, self.inputs.out_file)
        img_mean.to_filename(self._results['out_file'])

        return runtime


class _Ensure4DInputSpec(BaseInterfaceInputSpec):
    in_file = File(
        exists=True,
        mandatory=True,
        desc='A 4D image to index.',
    )
    out_file = File(
        'img_4d.nii.gz',
        usedefault=True,
        exists=False,
        desc='The name of the 4D file.',
    )


class _Ensure4DOutputSpec(TraitedSpec):
    out_file = File(
        exists=True,
        desc='The name of the 4D file.',
    )


class Ensure4D(NilearnBaseInterface, SimpleInterface):
    """Ensure a 4D image."""

    input_spec = _Ensure4DInputSpec
    output_spec = _Ensure4DOutputSpec

    def _run_interface(self, runtime):
        import nibabel as nb

        img = nb.load(self.inputs.in_file)
        if img.ndim == 3:
            img = img.slicer[..., None]
            self._results['out_file'] = os.path.join(runtime.cwd, self.inputs.out_file)
            img.to_filename(self._results['out_file'])
        else:
            self._results['out_file'] = self.inputs.in_file

        return runtime


class _GetImageTypeInputSpec(BaseInterfaceInputSpec):
    image = File(exists=True, mandatory=True)


class _GetImageTypeOutputSpec(TraitedSpec):
    image_type = traits.Enum(0, 1, 2, 3)


class GetImageType(SimpleInterface):
    """Use to determine what to send to --input-image-type."""

    input_spec = _GetImageTypeInputSpec
    output_spec = _GetImageTypeOutputSpec

    def _run_interface(self, runtime):
        img = nb.load(self.inputs.image)
        if img.ndim == 4:
            self._results['image_type'] = 3
        else:
            self._results['image_type'] = 0
        return runtime
