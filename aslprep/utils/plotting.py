"""Plotting functions and classes."""

import nibabel as nb
from lxml import etree
from nilearn import image, plotting
from nilearn._utils.niimg import load_niimg
from nireports._vendored.svgutils.transform import SVGFigure
from nireports.reportlets.utils import (
    compose_view,
    cuts_from_bbox,
    extract_svg,
    robust_set_limits,
)
from niworkflows import NIWORKFLOWS_LOG


class CBFPlot:
    """Generate the CBF Summary Plot.

    This plot restricts CBF values to -20 (if there are negative values) or 0 (if not) to 100.
    """

    __slots__ = ['cbf', 'ref_vol', 'label', 'outfile', 'vmax']

    def __init__(self, cbf, ref_vol, label, outfile, vmax):
        self.cbf = cbf
        self.ref_vol = ref_vol
        self.label = label
        self.outfile = outfile
        self.vmax = vmax

    def plot(self):
        """Generate the plot.

        This plot restricts CBF values to -20 (if there are negative values) or 0 (if not) to 100.
        """
        cbf_img = nb.load(self.cbf)
        cbf_data = cbf_img.get_fdata()
        cbf_data[cbf_data < -20] = -20
        cbf_data[cbf_data > 100] = 100
        cbf_img = nb.Nifti1Image(cbf_data, affine=cbf_img.affine, header=cbf_img.header)
        statfile = plot_stat_map(
            cbf=cbf_img,
            ref_vol=self.ref_vol,
            vmax=self.vmax,
            label=self.label,
        )
        compose_view(bg_svgs=statfile, fg_svgs=None, out_file=self.outfile)


def plot_stat_map(
    cbf,
    ref_vol,
    plot_params=None,
    order=('z', 'x', 'y'),
    vmax=100,
    estimate_brightness=False,
    label=None,
    compress='auto',
):
    """Plot statistical map."""
    plot_params = {} if plot_params is None else plot_params

    image_nii = load_niimg(cbf)
    if image_nii.ndim > 3:
        image_nii = image.mean_img(image_nii)

    data = image_nii.get_fdata()

    bbox_nii = image.threshold_img(image_nii, 1)

    cuts = cuts_from_bbox(bbox_nii, cuts=7)

    out_files = []
    if estimate_brightness:
        plot_params = robust_set_limits(data, plot_params)

    # Plot each cut axis
    for i, mode in enumerate(list(order)):
        display = plotting.plot_stat_map(
            stat_map_img=image_nii,
            bg_img=ref_vol,
            resampling_interpolation='nearest',
            display_mode=mode,
            cut_coords=cuts[mode],
            vmax=vmax,
            threshold=0.02,
            draw_cross=False,
            colorbar=True,
            symmetric_cbar=False,
            cmap='coolwarm',
            title=label if i == 0 else None,
        )
        svg = extract_svg(display, compress=compress)
        display.close()

        # Find and replace the figure_1 id.
        try:
            xml_data = etree.fromstring(svg)  # noqa: S320
        except etree.XMLSyntaxError as e:
            NIWORKFLOWS_LOG.info(e)
            return

        svg_fig = SVGFigure()
        svg_fig.root = xml_data
        out_files.append(svg_fig)

    return out_files
