# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Interfaces to generate speciality reportlets."""
import numpy as np
from nilearn.image import threshold_img, load_img
from niworkflows import NIWORKFLOWS_LOG
from niworkflows.viz.utils import cuts_from_bbox, compose_view
from nipype.interfaces.base import File, isdefined, traits
from nipype.interfaces.mixins import reporting

from ..viz.utils import plot_registration, coolwarm_transparent


class _FieldmapReportletInputSpec(reporting.ReportCapableInputSpec):
    reference = File(exists=True, mandatory=True, desc='input reference')
    moving = File(exists=True, desc='input moving')
    fieldmap = File(exists=True, mandatory=True, desc='input fieldmap')
    max_alpha = traits.Float(0.7, usedefault=True, desc='maximum alpha channel')
    mask = File(exists=True, desc='brain mask')
    out_report = File('report.svg', usedefault=True,
                      desc='filename for the visual report')
    show = traits.Enum(1, 0, 'both', usedefault=True,
                       desc='where the fieldmap should be shown')
    reference_label = traits.Str('Reference', usedefault=True,
                                 desc='a label name for the reference mosaic')
    moving_label = traits.Str('Fieldmap (Hz)', usedefault=True,
                              desc='a label name for the reference mosaic')


class FieldmapReportlet(reporting.ReportCapableInterface):
    """An abstract mixin to registration nipype interfaces."""

    _n_cuts = 7
    input_spec = _FieldmapReportletInputSpec
    output_spec = reporting.ReportCapableOutputSpec

    def __init__(self, **kwargs):
        """Instantiate FieldmapReportlet."""
        self._n_cuts = kwargs.pop('n_cuts', self._n_cuts)
        super(FieldmapReportlet, self).__init__(generate_report=True, **kwargs)

    def _run_interface(self, runtime):
        return runtime

    def _generate_report(self):
        """Generate a reportlet."""
        NIWORKFLOWS_LOG.info('Generating visual report')

        movnii = refnii = load_img(self.inputs.reference)
        fmapnii = load_img(self.inputs.fieldmap)

        if isdefined(self.inputs.moving):
            movnii = load_img(self.inputs.moving)

        contour_nii = mask_nii = None
        if isdefined(self.inputs.mask):
            contour_nii = load_img(self.inputs.mask)
            maskdata = contour_nii.get_fdata() > 0
        else:
            mask_nii = threshold_img(refnii, 1e-3)
            maskdata = mask_nii.get_fdata() > 0
        cuts = cuts_from_bbox(contour_nii or mask_nii, cuts=self._n_cuts)
        fmapdata = fmapnii.get_fdata()
        vmax = max(abs(np.percentile(fmapdata[maskdata], 99.8)),
                   abs(np.percentile(fmapdata[maskdata], 0.2)))

        fmap_overlay = [{
            'overlay': fmapnii,
            'overlay_params': {
                'cmap': coolwarm_transparent(max_alpha=self.inputs.max_alpha),
                'vmax': vmax,
                'vmin': -vmax,
            }
        }] * 2

        if self.inputs.show != 'both':
            fmap_overlay[not self.inputs.show] = {}

        # Call composer
        compose_view(
            plot_registration(movnii, 'moving-image',
                              estimate_brightness=True,
                              cuts=cuts,
                              label=self.inputs.moving_label,
                              contour=contour_nii,
                              compress=False,
                              **fmap_overlay[1]),
            plot_registration(refnii, 'fixed-image',
                              estimate_brightness=True,
                              cuts=cuts,
                              label=self.inputs.reference_label,
                              contour=contour_nii,
                              compress=False,
                              **fmap_overlay[0]),
            out_file=self._out_report
        )
