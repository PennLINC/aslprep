# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
SynthStrip interfaces
~~~~~~~~~~~~~~~~~~~~~

"""

import os.path as op

import nibabel as nb
from nipype.interfaces import afni
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.interfaces.freesurfer.base import FSCommandOpenMP, FSTraitedSpec
from niworkflows.utils.images import _copyxform


class FSTraitedSpecOpenMP(FSTraitedSpec):
    num_threads = traits.Int(desc='allows for specifying more threads', nohash=True)


def _get_padding(nvoxels):
    """Return the number of voxels that brings ``nvoxels`` up to a multiple of 64."""
    extra_slices = nvoxels % 64
    if extra_slices == 0:
        return 0
    complete_64s = nvoxels // 64
    return 64 * (complete_64s + 1) - nvoxels


def _split_padding(padding):
    halfpad = padding // 2
    return halfpad, halfpad + halfpad % 2


class _SynthStripPaddingInputSpec(BaseInterfaceInputSpec):
    input_image = File(exists=True, mandatory=True, desc='3D image to pad.')


class _SynthStripPaddingOutputSpec(TraitedSpec):
    S = traits.Int(desc='Number of slices to add at the superior edge.')
    R = traits.Int(desc='Number of slices to add at the right edge.')
    L = traits.Int(desc='Number of slices to add at the left edge.')
    A = traits.Int(desc='Number of slices to add at the anterior edge.')
    P = traits.Int(desc='Number of slices to add at the posterior edge.')


class SynthStripPadding(SimpleInterface):
    """Determine the zero-padding that puts an image on a SynthStrip-friendly grid.

    SynthStrip works on grids whose dimensions are multiples of 64 voxels.
    The padding is applied by a separate ``afni.Zeropad`` node.
    """

    input_spec = _SynthStripPaddingInputSpec
    output_spec = _SynthStripPaddingOutputSpec

    def _run_interface(self, runtime):
        img = nb.load(self.inputs.input_image)
        if not img.ndim == 3:
            raise ValueError('3D inputs are required for SynthStrip')
        xvoxels, yvoxels, zvoxels = img.shape

        self._results['S'] = _get_padding(zvoxels)
        self._results['R'], self._results['L'] = _split_padding(_get_padding(xvoxels))
        self._results['A'], self._results['P'] = _split_padding(_get_padding(yvoxels))
        return runtime


class Zeropad(afni.Zeropad):
    """AFNI's 3dZeropad, reporting the file it writes.

    Nipype's interface leaves ``out_file`` undefined in its outputs,
    so ``out_file`` must be set, and is returned as given.
    """

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = op.abspath(self.inputs.out_file)
        return outputs


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
