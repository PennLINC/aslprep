"""Plotting functions and classes."""

import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
import pandas as pd
from lxml import etree
from matplotlib import gridspec as mgs
from nilearn import image, plotting
from nilearn._utils.niimg import load_niimg
from niworkflows import NIWORKFLOWS_LOG
from niworkflows.viz.plots import confoundplot, spikesplot
from niworkflows.viz.utils import (
    compose_view,
    cuts_from_bbox,
    extract_svg,
    robust_set_limits,
)
from svgutils.transform import SVGFigure


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
            xml_data = etree.fromstring(svg)
        except etree.XMLSyntaxError as e:
            NIWORKFLOWS_LOG.info(e)
            return

        svg_fig = SVGFigure()
        svg_fig.root = xml_data
        out_files.append(svg_fig)

    return out_files


class fMRIPlot:  # noqa:N801
    """Generates the fMRI Summary Plot."""

    __slots__ = (
        'timeseries',
        'segments',
        'tr',
        'confounds',
        'spikes',
        'nskip',
        'sort_carpet',
        'paired_carpet',
    )

    def __init__(
        self,
        timeseries,
        segments,
        confounds=None,
        conf_file=None,
        tr=None,
        usecols=None,
        units=None,
        vlines=None,
        spikes_files=None,
        nskip=0,
        sort_carpet=True,
        paired_carpet=False,
    ):
        self.timeseries = timeseries
        self.segments = segments
        self.tr = tr
        self.nskip = nskip
        self.sort_carpet = sort_carpet
        self.paired_carpet = paired_carpet

        if units is None:
            units = {}
        if vlines is None:
            vlines = {}
        self.confounds = {}
        if confounds is None and conf_file:
            confounds = pd.read_csv(conf_file, sep=r'[\t\s]+', usecols=usecols, index_col=False)

        if confounds is not None:
            for name in confounds.columns:
                self.confounds[name] = {
                    'values': confounds[[name]].values.squeeze().tolist(),
                    'units': units.get(name),
                    'cutoff': vlines.get(name),
                }

        self.spikes = []
        if spikes_files:
            for sp_file in spikes_files:
                self.spikes.append((np.loadtxt(sp_file), None, False))

    def plot(self, figure=None):
        """Generate fMRI plot."""
        import seaborn as sns
        from niworkflows.viz.plots import plot_carpet as plt_carpet

        sns.set_style('whitegrid')
        sns.set_context('paper', font_scale=0.8)

        if figure is None:
            figure = plt.gcf()

        nconfounds = len(self.confounds)
        nspikes = len(self.spikes)
        nrows = 1 + nconfounds + nspikes

        # Create grid
        grid = mgs.GridSpec(
            nrows, 1, wspace=0.0, hspace=0.05, height_ratios=[1] * (nrows - 1) + [5]
        )

        grid_id = 0
        for tsz, name, iszs in self.spikes:
            spikesplot(tsz, title=name, outer_gs=grid[grid_id], tr=self.tr, zscored=iszs)
            grid_id += 1

        if self.confounds:
            from seaborn import color_palette

            palette = color_palette('husl', nconfounds)

        for i, (name, kwargs) in enumerate(self.confounds.items()):
            tseries = kwargs.pop('values')
            try:
                confoundplot(
                    tseries, grid[grid_id], tr=self.tr, color=palette[i], name=name, **kwargs
                )
            except Exception:
                raise ValueError(name) from None
            grid_id += 1

        plt_carpet(
            self.timeseries,
            segments=self.segments,
            subplot=grid[-1],
            tr=self.tr,
            sort_rows=self.sort_carpet,
            drop_trs=self.nskip,
            cmap='paired' if self.paired_carpet else None,
            # This is the only modification we need for ASLPrep
            detrend=False,
        )
        return figure
