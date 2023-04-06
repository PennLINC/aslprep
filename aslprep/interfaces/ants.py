"""ANTS interfaces."""
from nipype.interfaces.base import traits
from nipype.utils.filemanip import fname_presuffix
from niworkflows.interfaces.fixes import (
    FixHeaderApplyTransforms,
    _FixTraitApplyTransformsInputSpec,
)


class _ApplyTransformsInputSpec(_FixTraitApplyTransformsInputSpec):
    # Nipype's version doesn't have GenericLabel
    interpolation = traits.Enum(
        "Linear",
        "NearestNeighbor",
        "CosineWindowedSinc",
        "WelchWindowedSinc",
        "HammingWindowedSinc",
        "LanczosWindowedSinc",
        "MultiLabel",
        "Gaussian",
        "BSpline",
        "GenericLabel",
        argstr="%s",
        usedefault=True,
    )


class ApplyTransforms(FixHeaderApplyTransforms):
    """A modified version of FixHeaderApplyTransforms from niworkflows.

    The niworkflows version of ApplyTransforms "fixes the resampled image header
    to match the xform of the reference image".
    This modification overrides the allowed interpolation values,
    since FixHeaderApplyTransforms doesn't support GenericLabel,
    which is preferred over MultiLabel.
    """

    input_spec = _ApplyTransformsInputSpec

    def _run_interface(self, runtime):
        # Run normally
        self.inputs.output_image = fname_presuffix(
            self.inputs.input_image,
            suffix="_trans.nii.gz",
            newpath=runtime.cwd,
            use_ext=False,
        )
        runtime = super(ApplyTransforms, self)._run_interface(runtime)
        return runtime
