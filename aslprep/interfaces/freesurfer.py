# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
SynthStrip interfaces
~~~~~~~~~~~~~~~~~~~~~

"""

import os.path as op

import nibabel as nb
from nipype.interfaces.afni import Zeropad
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.interfaces.freesurfer.base import FSCommandOpenMP, FSTraitedSpec
from nipype.utils.filemanip import fname_presuffix
from niworkflows.utils.images import _copyxform


class FSTraitedSpecOpenMP(FSTraitedSpec):
    num_threads = traits.Int(desc='allows for specifying more threads', nohash=True)


class _PrepareSynthStripGridInputSpec(BaseInterfaceInputSpec):
    input_image = File(exists=True, mandatory=True)


class _PrepareSynthStripGridOutputSpec(TraitedSpec):
    prepared_image = File(exists=True)


class PrepareSynthStripGrid(SimpleInterface):
    input_spec = _PrepareSynthStripGridInputSpec
    output_spec = _PrepareSynthStripGridOutputSpec

    def _run_interface(self, runtime):
        out_fname = fname_presuffix(
            self.inputs.input_image,
            newpath=runtime.cwd,
            suffix='_SynthStripGrid.nii',
            use_ext=False,
        )
        self._results['prepared_image'] = out_fname

        # possibly downsample the image for sloppy mode. Always ensure float32
        img = nb.load(self.inputs.input_image)
        if not img.ndim == 3:
            raise Exception('3D inputs are required for Synthstrip')
        xvoxels, yvoxels, zvoxels = img.shape

        def get_padding(nvoxels):
            extra_slices = nvoxels % 64
            if extra_slices == 0:
                return 0
            complete_64s = nvoxels // 64
            return 64 * (complete_64s + 1) - nvoxels

        def split_padding(padding):
            halfpad = padding // 2
            return halfpad, halfpad + halfpad % 2

        spad = get_padding(zvoxels)
        rpad, lpad = split_padding(get_padding(xvoxels))
        apad, ppad = split_padding(get_padding(yvoxels))

        zeropad = Zeropad(
            S=spad,
            R=rpad,
            L=lpad,
            A=apad,
            P=ppad,
            in_files=self.inputs.input_image,
            out_file=out_fname,
        )

        _ = zeropad.run()
        assert op.exists(out_fname)
        return runtime


class _SynthStripInputSpec(FSTraitedSpecOpenMP):
    input_image = File(argstr='-i %s', exists=True, mandatory=True)
    no_csf = traits.Bool(argstr='--no-csf', desc='Exclude CSF from brain border.')
    border = traits.Int(argstr='-b %d', desc='Mask border threshold in mm. Default is 1.')
    gpu = traits.Bool(argstr='-g')
    out_brain = File(
        argstr='-o %s',
        name_template='%s_brain.nii.gz',
        name_source=['input_image'],
        keep_extension=False,
        desc='skull stripped image with corrupt sform',
    )
    out_brain_mask = File(
        argstr='-m %s',
        name_template='%s_mask.nii.gz',
        name_source=['input_image'],
        keep_extension=False,
        desc='mask image with corrupt sform',
    )


class _SynthStripOutputSpec(TraitedSpec):
    out_brain = File(exists=True)
    out_brain_mask = File(exists=True)


class SynthStrip(FSCommandOpenMP):
    input_spec = _SynthStripInputSpec
    output_spec = _SynthStripOutputSpec
    _cmd = 'mri_synthstrip'

    def _num_threads_update(self):
        if self.inputs.num_threads:
            self.inputs.environ.update({'OMP_NUM_THREADS': '1'})


class FixHeaderSynthStrip(SynthStrip):
    def _run_interface(self, runtime, correct_return_codes=(0,)):
        # Run normally
        runtime = super()._run_interface(runtime, correct_return_codes)

        outputs = self._list_outputs()
        if not op.exists(outputs['out_brain']):
            raise Exception('mri_synthstrip failed!')

        if outputs.get('out_brain_mask'):
            _copyxform(self.inputs.input_image, outputs['out_brain_mask'])

        _copyxform(self.inputs.input_image, outputs['out_brain'])

        return runtime


class MockSynthStrip(SimpleInterface):
    input_spec = _SynthStripInputSpec
    output_spec = _SynthStripOutputSpec

    def _run_interface(self, runtime):
        from nipype.interfaces.fsl import BET

        this_bet = BET(
            mask=True,
            in_file=self.inputs.input_image,
            output_type='NIFTI_GZ',
        )
        result = this_bet.run()
        self._results['out_brain'] = result.outputs.out_file
        self._results['out_brain_mask'] = result.outputs.mask_file

        return runtime
