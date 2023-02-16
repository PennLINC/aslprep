"""Plotting functions and classes."""
import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
import pandas as pd
import seaborn as sns
from lxml import etree
from matplotlib import gridspec as mgs
from nilearn.image import threshold_img
from nilearn.plotting import plot_stat_map
from seaborn import color_palette
from svgutils.transform import SVGFigure

from aslprep.niworkflows import NIWORKFLOWS_LOG
from aslprep.niworkflows.interfaces.plotting import _get_tr
from aslprep.niworkflows.viz.plots import confoundplot, plot_carpet, spikesplot
from aslprep.niworkflows.viz.utils import (
    _3d_in_file,
    compose_view,
    cuts_from_bbox,
    extract_svg,
    robust_set_limits,
)


class ASLPlot:
    """Generates the ASL Summary Plot."""

    __slots__ = ("func_file", "mask_data", "tr", "seg_data", "confounds", "spikes")

    def __init__(
        self,
        func_file,
        mask_file=None,
        data=None,
        conf_file=None,
        seg_file=None,
        tr=None,
        usecols=None,
        units=None,
        vlines=None,
        spikes_files=None,
    ):
        func_img = nb.load(func_file)
        self.func_file = func_file
        self.tr = tr or _get_tr(func_img)
        self.mask_data = None
        self.seg_data = None

        if not isinstance(func_img, nb.Cifti2Image):
            self.mask_data = nb.fileslice.strided_scalar(func_img.shape[:3], np.uint8(1))
            if mask_file:
                self.mask_data = np.asanyarray(nb.load(mask_file).dataobj).astype("uint8")
            if seg_file:
                self.seg_data = np.asanyarray(nb.load(seg_file).dataobj)

        if units is None:
            units = {}
        if vlines is None:
            vlines = {}
        self.confounds = {}
        if data is None and conf_file:
            data = pd.read_csv(conf_file, sep=r"[\t\s]+", usecols=usecols, index_col=False)

        if data is not None:
            for name in data.columns.ravel():
                self.confounds[name] = {
                    "values": data[[name]].values.ravel().tolist(),
                    "units": units.get(name),
                    "cutoff": vlines.get(name),
                }

        self.spikes = []
        if spikes_files:
            for sp_file in spikes_files:
                self.spikes.append((np.loadtxt(sp_file), None, False))

    def plot(self, figure=None):
        """Generate the plot."""
        sns.set_style("whitegrid")
        sns.set_context("paper", font_scale=0.8)

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
            palette = color_palette("husl", nconfounds)

        for i, (name, kwargs) in enumerate(self.confounds.items()):
            tseries = kwargs.pop("values")
            confoundplot(tseries, grid[grid_id], tr=self.tr, color=palette[i], name=name, **kwargs)
            grid_id += 1

        plot_carpet(self.func_file, atlaslabels=self.seg_data, subplot=grid[-1], tr=self.tr)
        # spikesplot_cb([0.7, 0.78, 0.2, 0.008])
        return figure


class CBFtsPlot(object):
    """Generate the CBF time series Summary Plot."""

    __slots__ = ["cbf_file", "tr", "seg_data", "fd_file"]

    def __init__(
        self,
        cbf_file,
        conf_file=None,
        seg_file=None,
        scoreindex=None,
        tr=None,
        usecols=None,
        units=None,
        vlines=None,
    ):
        self.cbf_file = cbf_file
        volindex = np.loadtxt(scoreindex)
        func_nii = nb.load(cbf_file)
        self.tr = tr if tr is not None else func_nii.header.get_zooms()[-1]

        self.seg_data = None
        if seg_file:
            self.seg_data = np.asanyarray(nb.load(seg_file).dataobj)

        if units is None:
            units = {}

        if vlines is None:
            vlines = {}

        self.fd_file = {}
        if conf_file:
            data = pd.read_csv(conf_file, sep=r"[\t\s]+", index_col=False)
            fdlist = data["framewise_displacement"].tolist()
            fdlist = fdlist[::2]
            # fdlist=np.nan_to_num(fdlist)
            self.fd_file["FD"] = {"values": fdlist}
        if scoreindex:
            self.fd_file["The SCORE index"] = {"values": volindex}

    def plot(self, figure=None):
        """Generate the plot."""
        sns.set_style("whitegrid")
        sns.set_context("paper", font_scale=0.8)

        if figure is None:
            figure = plt.gcf()

        nconfounds = len(self.fd_file)
        nrows = 1 + nconfounds

        # Create grid
        grid = mgs.GridSpec(
            nrows, 1, wspace=0.0, hspace=0.05, height_ratios=[1] * (nrows - 1) + [5]
        )

        grid_id = 0

        if self.fd_file:
            palette = color_palette("husl", nconfounds)

        for i, (name, kwargs) in enumerate(self.fd_file.items()):
            tseries = kwargs.pop("values")
            confoundplotx(
                tseries, grid[grid_id], tr=self.tr, color=palette[i], name=name, **kwargs
            )
            grid_id += 1

        plot_carpet(self.cbf_file, self.seg_data, subplot=grid[-1], tr=self.tr)
        # spikesplot_cb([0.7, 0.78, 0.2, 0.008])
        return figure


class CBFPlot(object):
    """Generate the CBF Summary Plot."""

    __slots__ = ["cbf", "ref_vol", "label", "outfile", "vmax"]

    def __init__(self, cbf, ref_vol, label, outfile, vmax):
        self.cbf = cbf
        self.ref_vol = ref_vol
        self.label = label
        self.outfile = outfile
        self.vmax = vmax

    def plot(self, figure=None):
        """Generate the plot."""
        statfile = plotstatsimg(
            cbf=self.cbf, ref_vol=self.ref_vol, vmax=self.vmax, label=self.label
        )
        compose_view(bg_svgs=statfile, fg_svgs=None, out_file=self.outfile)


def plotstatsimg(
    cbf,
    ref_vol,
    plot_params=None,
    order=("z", "x", "y"),
    vmax=90,
    estimate_brightness=False,
    label=None,
    compress="auto",
):
    """Plot statistical map."""
    plot_params = {} if plot_params is None else plot_params

    image_nii = _3d_in_file(cbf)
    data = image_nii.get_fdata()

    bbox_nii = threshold_img(nb.load(cbf), 1)

    cuts = cuts_from_bbox(bbox_nii, cuts=7)

    out_files = []
    if estimate_brightness:
        plot_params = robust_set_limits(data, plot_params)

    # Plot each cut axis
    for i, mode in enumerate(list(order)):
        plot_params["display_mode"] = mode
        plot_params["cut_coords"] = cuts[mode]
        plot_params["draw_cross"] = False
        plot_params["symmetric_cbar"] = True
        plot_params["vmax"] = vmax
        plot_params["threshold"] = 0.02
        plot_params["colorbar"] = True
        plot_params["cmap"] = "coolwarm"
        if i == 0:
            plot_params["title"] = label
            plot_params["colorbar"] = True
        else:
            plot_params["title"] = None

        display = plot_stat_map(
            stat_map_img=cbf, bg_img=ref_vol, resampling_interpolation="nearest", **plot_params
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


def confoundplotx(
    tseries,
    gs_ts,
    gs_dist=None,
    name=None,
    units=None,
    tr=None,
    hide_x=True,
    color="b",
    nskip=0,
    cutoff=None,
    ylims=None,
):
    """Generate a modified version of confoundplot from niworkflows.

    TODO: Figure out what was modified and why.
    TODO: See if I can use the regular confoundplot instead.
    """
    # Define TR and number of frames
    notr = False
    if tr is None:
        notr = True
        tr = 1.0
    ntsteps = len(tseries)
    tseries = np.array(tseries)

    # Define nested GridSpec
    gs = mgs.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_ts, width_ratios=[1, 100], wspace=0.0)

    ax_ts = plt.subplot(gs[1])
    ax_ts.grid(False)

    # Set 10 frame markers in X axis
    interval = max((ntsteps // 10, ntsteps // 5, 1))
    xticks = list(range(0, ntsteps)[::interval])
    ax_ts.set_xticks(xticks)

    if not hide_x:
        if notr:
            ax_ts.set_xlabel("time (frame #)")
        else:
            ax_ts.set_xlabel("time (s)")
            labels = tr * np.array(xticks)
            ax_ts.set_xticklabels([f"{t:.02f}" for t in labels.tolist()])
    else:
        ax_ts.set_xticklabels([])

    if name is not None:
        if units is not None:
            name += f" [{units}]"

        ax_ts.annotate(
            name,
            xy=(0.0, 0.7),
            xytext=(0, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            va="center",
            ha="left",
            color=color,
            size=8,
            bbox={
                "boxstyle": "round",
                "fc": "w",
                "ec": "none",
                "color": "none",
                "lw": 0,
                "alpha": 0.8,
            },
        )

    for side in ["top", "right"]:
        ax_ts.spines[side].set_color("none")
        ax_ts.spines[side].set_visible(False)

    if not hide_x:
        ax_ts.spines["bottom"].set_position(("outward", 20))
        ax_ts.xaxis.set_ticks_position("bottom")
    else:
        ax_ts.spines["bottom"].set_color("none")
        ax_ts.spines["bottom"].set_visible(False)

    # ax_ts.spines["left"].set_position(('outward', 30))
    ax_ts.spines["left"].set_color("none")
    ax_ts.spines["left"].set_visible(False)
    # ax_ts.yaxis.set_ticks_position('left')

    ax_ts.set_yticks([])
    ax_ts.set_yticklabels([])

    nonnan = tseries[~np.isnan(tseries)]
    if nonnan.size > 0:
        # Calculate Y limits
        valrange = nonnan.max() - nonnan.min()
        def_ylims = [nonnan.min() - 0.1 * valrange, nonnan.max() + 0.1 * valrange]
        if ylims is not None:
            if ylims[0] is not None:
                def_ylims[0] = min([def_ylims[0], ylims[0]])
            if ylims[1] is not None:
                def_ylims[1] = max([def_ylims[1], ylims[1]])

        # Add space for plot title and mean/SD annotation
        def_ylims[0] -= 0.1 * (def_ylims[1] - def_ylims[0])

        ax_ts.set_ylim(def_ylims)

        # Annotate stats
        maxv = nonnan.max()
        mean = nonnan.mean()
        stdv = nonnan.std()
        p95 = np.percentile(nonnan, 95.0)
    else:
        maxv = 0
        mean = 0
        stdv = 0
        p95 = 0

    stats_label = (
        r"max: {max:.3f}{units} $\bullet$ mean: {mean:.3f}{units} "
        r"$\bullet$ $\sigma$: {sigma:.3f}"
    ).format(max=maxv, mean=mean, units=units or "", sigma=stdv)
    ax_ts.annotate(
        stats_label,
        xy=(0.98, 0.7),
        xycoords="axes fraction",
        xytext=(0, 0),
        textcoords="offset points",
        va="center",
        ha="right",
        color=color,
        size=4,
        bbox={
            "boxstyle": "round",
            "fc": "w",
            "ec": "none",
            "color": "none",
            "lw": 0,
            "alpha": 0.8,
        },
    )

    # Annotate percentile 95
    ax_ts.plot((0, ntsteps - 1), [p95] * 2, linewidth=0.1, color="lightgray")
    ax_ts.annotate(
        f"{p95:.2f}",
        xy=(0, p95),
        xytext=(-1, 0),
        textcoords="offset points",
        va="center",
        ha="right",
        color="lightgray",
        size=3,
    )

    if cutoff is None:
        cutoff = []

    for i, thr in enumerate(cutoff):
        ax_ts.plot((0, ntsteps - 1), [thr] * 2, linewidth=0.2, color="dimgray")

        ax_ts.annotate(
            f"{thr:.2f}",
            xy=(0, thr),
            xytext=(-1, 0),
            textcoords="offset points",
            va="center",
            ha="right",
            color="dimgray",
            size=3,
        )

    # ax_ts.plot(tseries, color=color, linewidth=.8)
    # ax_ts.set_xlim((0, ntsteps - 1))
    ax_ts.step(range(0, ntsteps), tseries, color=color)
    ax_ts.set_xlim((0, ntsteps - 1))
    if gs_dist is not None:
        ax_dist = plt.subplot(gs_dist)
        sns.displot(tseries, vertical=True, ax=ax_dist)
        ax_dist.set_xlabel("Timesteps")
        ax_dist.set_ylim(ax_ts.get_ylim())
        ax_dist.set_yticklabels([])

        return [ax_ts, ax_dist], gs
    return ax_ts, gs
