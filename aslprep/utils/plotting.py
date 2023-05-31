"""Plotting functions and classes."""
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
import pandas as pd
import seaborn as sns
from lxml import etree
from matplotlib import gridspec as mgs
from matplotlib.colors import ListedColormap
from nilearn import image, plotting
from nilearn._utils.niimg import _safe_get_data, check_niimg_4d, load_niimg
from niworkflows import NIWORKFLOWS_LOG
from niworkflows.interfaces.plotting import _get_tr
from niworkflows.viz.plots import confoundplot, spikesplot
from niworkflows.viz.utils import (
    compose_view,
    cuts_from_bbox,
    extract_svg,
    robust_set_limits,
)
from seaborn import color_palette
from svgutils.transform import SVGFigure


class ASLPlot:
    """Generates the ASL Summary Plot."""

    __slots__ = ("func_file", "mask_data", "tr", "seg_data", "confounds", "spikes")

    def __init__(
        self,
        func_file,
        mask_file=None,
        data=None,
        confounds_file=None,
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
        if data is None and confounds_file:
            data = pd.read_csv(confounds_file, sep=r"[\t\s]+", usecols=usecols, index_col=False)

        if data is not None:
            for name in data.columns.ravel():
                self.confounds[name] = {
                    "values": data[[name]].values.ravel().tolist(),
                    "units": units.get(name),
                    "cutoff": vlines.get(name),
                }

        self.spikes = []
        if spikes_files:
            self.spikes.extend((np.loadtxt(sp_file), None, False) for sp_file in spikes_files)

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
        confounds_file=None,
        seg_file=None,
        score_outlier_index=None,
        tr=None,
        units=None,
        vlines=None,
    ):
        self.cbf_file = cbf_file
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
        if confounds_file:
            data = pd.read_csv(confounds_file, sep=r"[\t\s]+", index_col=False)
            fdlist = data["framewise_displacement"].tolist()
            fdlist = fdlist[::2]
            # fdlist=np.nan_to_num(fdlist)
            self.fd_file["FD"] = {"values": fdlist}

        if score_outlier_index:
            volindex = np.loadtxt(score_outlier_index)
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
    """Generate the CBF Summary Plot.

    This plot restricts CBF values to -20 (if there are negative values) or 0 (if not) to 100.
    """

    __slots__ = ["cbf", "ref_vol", "label", "outfile", "vmax"]

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
    order=("z", "x", "y"),
    vmax=100,
    estimate_brightness=False,
    label=None,
    compress="auto",
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
            resampling_interpolation="nearest",
            display_mode=mode,
            cut_coords=cuts[mode],
            vmax=vmax,
            threshold=0.02,
            draw_cross=False,
            colorbar=True,
            symmetric_cbar=False,
            cmap="coolwarm",
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


def confoundplotx(
    tseries,
    gs_ts,
    gs_dist=None,
    name=None,
    units=None,
    tr=None,
    hide_x=True,
    color="b",
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
    xticks = list(range(ntsteps)[::interval])
    ax_ts.set_xticks(xticks)

    if hide_x:
        ax_ts.set_xticklabels([])

    elif notr:
        ax_ts.set_xlabel("time (frame #)")
    else:
        ax_ts.set_xlabel("time (s)")
        labels = tr * np.array(xticks)
        ax_ts.set_xticklabels([f"{t:.02f}" for t in labels.tolist()])
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

    for thr in cutoff:
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
    ax_ts.step(range(ntsteps), tseries, color=color)
    ax_ts.set_xlim((0, ntsteps - 1))
    if gs_dist is not None:
        ax_dist = plt.subplot(gs_dist)
        sns.displot(tseries, vertical=True, ax=ax_dist)
        ax_dist.set_xlabel("Timesteps")
        ax_dist.set_ylim(ax_ts.get_ylim())
        ax_dist.set_yticklabels([])

        return [ax_ts, ax_dist], gs
    return ax_ts, gs


def plot_carpet(
    func,
    atlaslabels=None,
    detrend=True,
    size=(950, 800),
    subplot=None,
    output_file=None,
    legend=False,
    tr=None,
    lut=None,
):
    """Generate a carpet plot.

    Plot an image representation of voxel intensities across time also know
    as the "carpet plot" or "Power plot". See Jonathan Power Neuroimage
    2017 Jul 1; 154:150-158.

    Parameters
    ----------
    func : string
        Path to NIfTI or CIFTI ASL image
    atlaslabels: ndarray, optional
        A 3D array of integer labels from an atlas, resampled into ``img`` space.
        Required if ``func`` is a NIfTI image.
    detrend : boolean, optional
        Detrend and standardize the data prior to plotting.
    nskip : int, optional
        Number of volumes at the beginning of the scan marked as nonsteady state.
        Not used.
    size : tuple, optional
        Size of figure.
    subplot : matplotlib Subplot, optional
        Subplot to plot figure on.
    title : string, optional
        The title displayed on the figure.
    output_file : string, or None, optional
        The name of an image file to export the plot to. Valid extensions
        are .png, .pdf, .svg. If output_file is not None, the plot
        is saved to a file, and the display is closed.
    legend : bool
        Whether to render the average functional series with ``atlaslabels`` as
        overlay.
    tr : float , optional
        Specify the TR, if specified it uses this value. If left as None,
        # of frames is plotted instead of time.
    lut : ndarray, optional
        Look up table for segmentations
    """
    epinii = None
    segnii = None
    nslices = None
    img = nb.load(func)

    if isinstance(img, nb.Cifti2Image):
        assert img.nifti_header.get_intent()[0] == "ConnDenseSeries", "Not a dense timeseries"

        data = img.get_fdata().T
        matrix = img.header.matrix
        struct_map = {
            "LEFT_CORTEX": 1,
            "RIGHT_CORTEX": 2,
            "SUBCORTICAL": 3,
            "CEREBELLUM": 4,
        }
        seg = np.zeros((data.shape[0],), dtype="uint32")
        for bm in matrix.get_index_map(1).brain_models:
            if "CORTEX" in bm.brain_structure:
                lidx = (1, 2)["RIGHT" in bm.brain_structure]
            elif "CEREBELLUM" in bm.brain_structure:
                lidx = 4
            else:
                lidx = 3
            index_final = bm.index_offset + bm.index_count
            seg[bm.index_offset : index_final] = lidx
        assert len(seg[seg < 1]) == 0, "Unassigned labels"

        # Decimate data
        data, seg = _decimate_data(data, seg, size)
        # preserve as much continuity as possible
        order = seg.argsort(kind="stable")

        cmap = ListedColormap([cm.get_cmap("Paired").colors[i] for i in (1, 0, 7, 3)])
        assert len(cmap.colors) == len(
            struct_map
        ), "Mismatch between expected # of structures and colors"

        # ensure no legend for CIFTI
        legend = False

    else:  # Volumetric NIfTI
        img_nii = check_niimg_4d(
            img,
            dtype="auto",
        )
        func_data = _safe_get_data(img_nii, ensure_finite=True)
        ntsteps = func_data.shape[-1]
        data = func_data[atlaslabels > 0].reshape(-1, ntsteps)
        oseg = atlaslabels[atlaslabels > 0].reshape(-1)

        # Map segmentation
        if lut is None:
            lut = np.zeros((256,), dtype="int")
            lut[1:11] = 1
            lut[255] = 2
            lut[30:99] = 3
            lut[100:201] = 4
        # Apply lookup table
        seg = lut[oseg.astype(int)]

        # Decimate data
        data, seg = _decimate_data(data, seg, size)
        # Order following segmentation labels
        order = np.argsort(seg)[::-1]
        # Set colormap
        cmap = ListedColormap(cm.get_cmap("tab10").colors[:4][::-1])

        if legend:
            epiavg = func_data.mean(3)
            epinii = nb.Nifti1Image(epiavg, img_nii.affine, img_nii.header)
            segnii = nb.Nifti1Image(lut[atlaslabels.astype(int)], epinii.affine, epinii.header)
            segnii.set_data_dtype("uint8")
            nslices = epiavg.shape[-1]

    return _carpet(
        data,
        seg,
        order,
        cmap,
        epinii=epinii,
        segnii=segnii,
        nslices=nslices,
        tr=tr,
        detrend=detrend,
        subplot=subplot,
        output_file=output_file,
    )


def _carpet(
    data,
    seg,
    order,
    cmap,
    tr=None,
    detrend=True,
    subplot=None,
    legend=False,
    output_file=None,
    epinii=None,
    segnii=None,
    nslices=None,
):
    notr = False
    if tr is None:
        notr = True
        tr = 1.0

    # Detrend data
    v = (None, None)
    if detrend:
        from nilearn.signal import clean

        data = clean(data.T, t_r=tr).T
        v = (-2, 2)

    # If subplot is not defined
    if subplot is None:
        subplot = mgs.GridSpec(1, 1)[0]

    # Define nested GridSpec
    wratios = [1, 100, 20]
    gs = mgs.GridSpecFromSubplotSpec(
        1,
        2 + int(legend),
        subplot_spec=subplot,
        width_ratios=wratios[: 2 + int(legend)],
        wspace=0.0,
    )

    # Segmentation colorbar
    ax0 = plt.subplot(gs[0])
    ax0.set_yticks([])
    ax0.set_xticks([])
    ax0.imshow(seg[order, np.newaxis], interpolation="none", aspect="auto", cmap=cmap)

    ax0.grid(False)
    ax0.spines["left"].set_visible(False)
    ax0.spines["bottom"].set_color("none")
    ax0.spines["bottom"].set_visible(False)

    # Carpet plot
    ax1 = plt.subplot(gs[1])
    ax1.imshow(
        data[order],
        interpolation="nearest",
        aspect="auto",
        cmap="gray",
        vmin=v[0],
        vmax=v[1],
    )

    ax1.grid(False)
    ax1.set_yticks([])
    ax1.set_yticklabels([])

    # Set 10 frame markers in X axis
    interval = max((int(data.shape[-1] + 1) // 10, int(data.shape[-1] + 1) // 5, 1))
    xticks = list(range(data.shape[-1])[::interval])
    ax1.set_xticks(xticks)
    ax1.set_xlabel("time (frame #)" if notr else "time (s)")
    labels = tr * (np.array(xticks))
    ax1.set_xticklabels([f"{t:.02f}" for t in labels.tolist()], fontsize=5)

    # Remove and redefine spines
    for side in ["top", "right"]:
        # Toggle the spine objects
        ax0.spines[side].set_color("none")
        ax0.spines[side].set_visible(False)
        ax1.spines[side].set_color("none")
        ax1.spines[side].set_visible(False)

    ax1.yaxis.set_ticks_position("left")
    ax1.xaxis.set_ticks_position("bottom")
    ax1.spines["bottom"].set_visible(False)
    ax1.spines["left"].set_color("none")
    ax1.spines["left"].set_visible(False)

    ax2 = None
    if legend:
        gslegend = mgs.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs[2], wspace=0.0, hspace=0.0)
        coords = np.linspace(int(0.10 * nslices), int(0.95 * nslices), 5).astype(np.uint8)
        for i, c in enumerate(coords.tolist()):
            ax2 = plt.subplot(gslegend[i])
            plotting.plot_img(
                segnii,
                bg_img=epinii,
                axes=ax2,
                display_mode="z",
                annotate=False,
                cut_coords=[c],
                threshold=0.1,
                cmap=cmap,
                interpolation="nearest",
            )

    if output_file is not None:
        figure = plt.gcf()
        figure.savefig(output_file, bbox_inches="tight")
        plt.close(figure)
        figure = None
        return output_file

    return (ax0, ax1, ax2), gs


def _decimate_data(data, seg, size):
    """Decimate timeseries data.

    Parameters
    ----------
    data : ndarray
        2 element array of timepoints and samples
    seg : ndarray
        1 element array of samples
    size : tuple
        2 element for P/T decimation
    """
    p_dec = 1 + data.shape[0] // size[0]
    if p_dec:
        data = data[::p_dec, :]
        seg = seg[::p_dec]
    t_dec = 1 + data.shape[1] // size[1]
    if t_dec:
        data = data[:, ::t_dec]
    return data, seg
