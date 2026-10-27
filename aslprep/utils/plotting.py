"""Plotting functions and classes."""

from numbers import Number

import nibabel as nb
import numpy as np
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

    __slots__ = ['in_file', 'label', 'outfile', 'ref_vol', 'vmax', 'vmin']

    def __init__(self, in_file, ref_vol, label, outfile, vmin, vmax):
        self.in_file = in_file
        self.ref_vol = ref_vol
        self.label = label
        self.outfile = outfile
        self.vmin = vmin
        self.vmax = vmax

    def plot(self):
        """Generate the plot.

        This plot restricts CBF values to -20 (if there are negative values) or 0 (if not) to 100.
        """
        img = nb.load(self.in_file)
        data = img.get_fdata()
        data[data < self.vmin] = self.vmin
        data[data > self.vmax] = self.vmax
        if np.any(data < 0):
            colormap = 'coolwarm'
            vmin = -20
        else:
            colormap = 'Reds'
            vmin = 0

        img = nb.Nifti1Image(data, affine=img.affine, header=img.header)
        statfile = plot_stat_map(
            img=img,
            ref_vol=self.ref_vol,
            vmin=vmin,
            vmax=self.vmax,
            label=self.label,
            cmap=colormap,
        )
        compose_view(bg_svgs=statfile, fg_svgs=None, out_file=self.outfile)


def plot_stat_map(
    img,
    ref_vol,
    plot_params=None,
    order=('z', 'x', 'y'),
    vmin=None,
    vmax=100,
    estimate_brightness=False,
    label=None,
    compress='auto',
    cmap='coolwarm',
):
    """Plot statistical map."""
    plot_params = {} if plot_params is None else plot_params

    image_nii = load_niimg(img)
    if image_nii.ndim > 3:
        image_nii = image.mean_img(image_nii)

    data = image_nii.get_fdata()

    bbox_nii = image.threshold_img(image_nii, 1)

    cuts = cuts_from_bbox(bbox_nii, cuts=7)

    out_files = []
    if estimate_brightness:
        plot_params = robust_set_limits(data, plot_params)

    if isinstance(vmin, Number) and (vmin >= 0):
        # Scale from vmin (0) to vmax (100)
        symmetric_cbar = False
        arg_vmin = vmin
    else:
        # Scale from -vmax (-100) to vmax (100)
        symmetric_cbar = True
        arg_vmin = None

    # Plot each cut axis
    for i, mode in enumerate(list(order)):
        display = plotting.plot_stat_map(
            stat_map_img=image_nii,
            bg_img=ref_vol,
            resampling_interpolation='nearest',
            display_mode=mode,
            cut_coords=cuts[mode],
            vmin=arg_vmin,
            vmax=vmax,
            threshold=0.00001,
            draw_cross=False,
            colorbar=True,
            symmetric_cbar=symmetric_cbar,
            cmap=cmap,
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


TISSUE_COLORS = {'GM': '#1b60a5', 'WM': '#2da467', 'CSF': '#9d8f25'}


def plot_tissue_ridgeplot(panels, out_file):
    """Plot overlapping density ridges of one or more parameters by map and tissue type.

    Each panel is one parameter (e.g., CBF or ATT), with its own x-axis. Each row of a panel
    is one map of that parameter (e.g., basic, SCORE, or BASIL CBF), with one filled density
    per tissue type. The rows overlap, as in seaborn's ``kde_ridgeplot`` example.
    Each density integrates to one on its own, so a tissue with few voxels (CSF) is as visible
    as one with many.

    Parameters
    ----------
    panels : :obj:`dict`
        Maps each panel's x-axis label to a dict of rows. That dict maps each row label to a
        dict that maps each tissue type in :data:`TISSUE_COLORS` to a 1D array of voxel values.
        Panels and rows are drawn from top to bottom in this order.
    out_file : :obj:`str`
        Path to the output figure.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Layout in inches: the step between the baselines of consecutive rows,
    # how much of a row's height reaches into the row above,
    # the space for the legend, and the space below each panel for its x-axis.
    width, left, right = 9, 1.8, 0.3
    row_step, overlap = 0.75, 0.3
    row_height = row_step / (1 - overlap)
    legend_space, xaxis_space, panel_gap = 0.8, 0.6, 0.35
    panel_heights = [row_height + (len(rows) - 1) * row_step for rows in panels.values()]
    height = (
        legend_space
        + sum(panel_heights)
        + len(panels) * xaxis_space
        + (len(panels) - 1) * panel_gap
    )

    with sns.axes_style('white', rc={'axes.facecolor': (0, 0, 0, 0)}):
        fig = plt.figure(figsize=(width, height))
        top = height - legend_space
        for (x_label, rows), panel_height in zip(panels.items(), panel_heights, strict=True):
            xlim = _robust_xlim(rows)
            n_rows = len(rows)
            for i_row, (row_label, row) in enumerate(rows.items()):
                bottom = top - row_height - i_row * row_step
                ax = fig.add_axes(
                    [
                        left / width,
                        bottom / height,
                        (width - left - right) / width,
                        row_height / height,
                    ],
                )
                # Stack each row above the one below it, so the tails of a lower row's
                # densities don't hide this row's baseline and median ticks.
                ax.set_zorder(n_rows - i_row)
                _plot_ridge_row(ax, row, row_label, xlim)
                if i_row == n_rows - 1:
                    ax.set_xlabel(x_label)
                else:
                    ax.tick_params(bottom=False, labelbottom=False)

            top -= panel_height + xaxis_space + panel_gap

        handles = [
            plt.Rectangle((0, 0), 1, 1, color=color, alpha=0.6) for color in TISSUE_COLORS.values()
        ]
        fig.legend(
            handles,
            TISSUE_COLORS.keys(),
            loc='upper center',
            ncol=len(TISSUE_COLORS),
            frameon=False,
            title='Tissue type (dotted lines: medians)',
        )
        fig.savefig(out_file)
        plt.close(fig)


def _robust_xlim(rows):
    """Clip a panel's x-axis to the bulk of its values, and include zero.

    A few extreme voxels from one map would otherwise flatten the densities of the others.
    """
    values = [v for row in rows.values() for v in row.values() if v.size]
    if not values:
        return np.array([0.0, 1.0])

    xlim = np.percentile(np.concatenate(values), [0.5, 99.5])
    xlim[0] = min(xlim[0], 0)
    if xlim[0] == xlim[1]:
        xlim = xlim + np.array([-1, 1])

    return xlim


def _plot_ridge_row(ax, row, row_label, xlim):
    """Plot one row of a ridge plot: a filled density and median tick per tissue type."""
    import seaborn as sns

    for tissue_type, color in TISSUE_COLORS.items():
        tissue_vals = row.get(tissue_type, np.empty(0))
        # A KDE needs at least two distinct values.
        if np.unique(tissue_vals).size < 2:
            continue

        kde_kws = {'x': tissue_vals, 'clip': xlim, 'bw_adjust': 0.75, 'ax': ax}
        sns.kdeplot(fill=True, alpha=0.45, linewidth=0, color=color, **kde_kws)
        sns.kdeplot(linewidth=1.5, color=color, **kde_kws)
        ax.axvline(
            np.median(tissue_vals),
            ymax=0.25,
            color=color,
            linewidth=1.5,
            linestyle=':',
        )

    ax.axhline(0, color='0.3', linewidth=1, clip_on=False)
    ax.text(
        -0.01,
        0.1,
        row_label,
        fontweight='bold',
        ha='right',
        va='center',
        transform=ax.transAxes,
    )
    ax.set_yticks([])
    ax.set_ylabel('')
    ax.set_xlim(xlim)
    sns.despine(ax=ax, left=True, bottom=True)
