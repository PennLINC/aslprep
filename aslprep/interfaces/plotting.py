"""Plotting interfaces."""

import nibabel as nb
import numpy as np
import pandas as pd
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    InputMultiObject,
    SimpleInterface,
    TraitedSpec,
    isdefined,
    traits,
)
from nipype.utils.filemanip import fname_presuffix
from nireports.reportlets.modality.func import fMRIPlot
from niworkflows.utils.timeseries import _cifti_timeseries, _nifti_timeseries

from aslprep.utils.plotting import CBFPlot


class _ASLCarpetPlotInputSpec(BaseInterfaceInputSpec):
    in_nifti = File(exists=True, mandatory=True, desc='input ASL (4D NIfTI file)')
    in_cifti = File(exists=True, desc='input ASL (CIFTI dense timeseries)')
    in_segm = File(exists=True, desc='volumetric segmentation corresponding to in_nifti')
    confounds_file = File(exists=True, desc="BIDS' _confounds.tsv file")

    str_or_tuple = traits.Either(
        traits.Str,
        traits.Tuple(traits.Str, traits.Either(None, traits.Str)),
        traits.Tuple(traits.Str, traits.Either(None, traits.Str), traits.Either(None, traits.Str)),
    )
    confounds_list = traits.Either(
        traits.List(str_or_tuple, minlen=1),
        None,
        desc='list of headers to extract from the confounds_file',
    )
    tr = traits.Either(None, traits.Float, usedefault=True, desc='the repetition time')
    drop_trs = traits.Int(0, usedefault=True, desc='dummy scans')


class _ASLCarpetPlotOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='written file path')


class ASLCarpetPlot(SimpleInterface):
    """Create a combined carpet/line plot for ASL time series data.

    Copy the x-form matrices from `hdr_file` to `out_file`.
    """

    input_spec = _ASLCarpetPlotInputSpec
    output_spec = _ASLCarpetPlotOutputSpec

    def _run_interface(self, runtime):
        self._results['out_file'] = fname_presuffix(
            self.inputs.in_nifti, suffix='_fmriplot.svg', use_ext=False, newpath=runtime.cwd
        )

        has_cifti = isdefined(self.inputs.in_cifti)

        # Read input object and create timeseries + segments object
        seg_file = self.inputs.in_segm if isdefined(self.inputs.in_segm) else None
        dataset, segments = _nifti_timeseries(
            nb.load(self.inputs.in_nifti),
            nb.load(seg_file),
            remap_rois=False,
            labels=(
                ('WM+CSF', 'Edge')
                if has_cifti
                else ('Ctx GM', 'dGM', 'sWM+sCSF', 'dWM+dCSF', 'Cb', 'Edge')
            ),
        )

        # Process CIFTI
        if has_cifti:
            cifti_data, cifti_segments = _cifti_timeseries(nb.load(self.inputs.in_cifti))

            if seg_file is not None:
                # Append WM+CSF and Edge masks
                cifti_length = cifti_data.shape[0]
                dataset = np.vstack((cifti_data, dataset))
                segments = {k: np.array(v) + cifti_length for k, v in segments.items()}
                cifti_segments.update(segments)
                segments = cifti_segments
            else:
                dataset, segments = cifti_data, cifti_segments

        dataframe = pd.read_table(
            self.inputs.confounds_file,
            index_col=None,
            dtype='float32',
            na_filter=True,
            na_values='n/a',
        )

        headers = []
        units = {}
        names = {}

        for conf_el in self.inputs.confounds_list:
            if isinstance(conf_el, list | tuple):
                headers.append(conf_el[0])
                if conf_el[1] is not None:
                    units[conf_el[0]] = conf_el[1]

                if len(conf_el) > 2 and conf_el[2] is not None:
                    names[conf_el[0]] = conf_el[2]
            else:
                headers.append(conf_el)

        if not headers:
            data = None
            units = None
        else:
            data = dataframe[headers]
            data = data.rename(columns=names)

        fig = fMRIPlot(
            dataset,
            segments=segments,
            tr=self.inputs.tr,
            confounds=data,
            units=units,
            nskip=self.inputs.drop_trs,
            paired_carpet=has_cifti,
            detrend=False,
        ).plot()
        fig.savefig(self._results['out_file'], bbox_inches='tight')
        fig.clf()
        return runtime


class _CBFSummaryPlotInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='Image to plot')
    label = traits.Str(mandatory=True, desc='label')
    vmin = traits.Float(-20, usedefault=True, mandatory=False, desc='Minimum value for figure')
    vmax = traits.Float(100, usedefault=True, mandatory=True, desc='Maximum value for figure')
    ref_vol = File(exists=True, mandatory=True, desc='Reference volume to use as underlay')


class _CBFSummaryPlotOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='written file path')


class CBFSummaryPlot(SimpleInterface):
    """Prepare an CBF summary plot for the report.

    This plot restricts CBF values to -20 (if there are negative values) or 0 (if not) to 100.
    """

    input_spec = _CBFSummaryPlotInputSpec
    output_spec = _CBFSummaryPlotOutputSpec

    def _run_interface(self, runtime):
        self._results['out_file'] = fname_presuffix(
            self.inputs.in_file,
            suffix=f'_{self.inputs.label}plot.svg',
            use_ext=False,
            newpath=runtime.cwd,
        )
        CBFPlot(
            in_file=self.inputs.in_file,
            label=self.inputs.label,
            ref_vol=self.inputs.ref_vol,
            vmin=self.inputs.vmin,
            vmax=self.inputs.vmax,
            outfile=self._results['out_file'],
        ).plot()
        return runtime


class _CBFByTissueTypePlotInputSpec(BaseInterfaceInputSpec):
    in_files = InputMultiObject(
        File(exists=True),
        mandatory=True,
        desc='Maps to plot, one row each, all in the space of seg_file',
    )
    labels = traits.List(
        traits.Str,
        mandatory=True,
        desc='Row label for each map in in_files',
    )
    seg_file = File(
        exists=True,
        mandatory=True,
        desc='Segmentation file, with 1 for GM, 2 for WM and 3 for CSF',
    )


class _ByTissueTypePlotOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='written file path')


def _values_by_tissue(in_file, seg):
    """Split a map's values by tissue type.

    Voxels that are exactly zero or not finite are left out, since the gray matter partial
    volume-corrected maps use zero for voxels where gray matter CBF was not estimated.
    """
    from aslprep.utils.plotting import TISSUE_COLORS

    data = nb.load(in_file).get_fdata(dtype=np.float32)
    if data.shape[:3] != seg.shape[:3]:
        raise ValueError(
            f'{in_file} has shape {data.shape}, but the segmentation has shape {seg.shape}.'
        )

    keep = np.isfinite(data) & (data != 0)
    return {
        tissue_type: data[keep & (seg == i_tissue_type + 1)]
        for i_tissue_type, tissue_type in enumerate(TISSUE_COLORS)
    }


class CBFByTissueTypePlot(SimpleInterface):
    """Plot the distributions of one or more CBF maps by tissue type, as a ridge plot.

    Voxels that are exactly zero or not finite are left out of the distributions.
    """

    input_spec = _CBFByTissueTypePlotInputSpec
    output_spec = _ByTissueTypePlotOutputSpec

    def _run_interface(self, runtime):
        from aslprep.utils.plotting import plot_tissue_ridgeplot

        if len(self.inputs.in_files) != len(self.inputs.labels):
            raise ValueError(
                f'Got {len(self.inputs.in_files)} files but {len(self.inputs.labels)} labels.'
            )

        self._results['out_file'] = fname_presuffix(
            self.inputs.in_files[0],
            suffix='_cbfbytissue.svg',
            use_ext=False,
            newpath=runtime.cwd,
        )

        seg = np.asanyarray(nb.load(self.inputs.seg_file).dataobj)
        rows = {
            label: _values_by_tissue(in_file, seg)
            for in_file, label in zip(self.inputs.in_files, self.inputs.labels, strict=True)
        }
        plot_tissue_ridgeplot(
            {'Cerebral Blood Flow (mL/100 g/min)': rows},
            out_file=self._results['out_file'],
        )
        return runtime


class _MultiPLDByTissueTypePlotInputSpec(BaseInterfaceInputSpec):
    att = File(exists=True, desc='Arterial transit time from the basic model')
    abat = File(exists=True, desc='Arterial bolus arrival time from the basic model')
    abv = File(exists=True, desc='Arterial blood volume from the basic model')
    att_basil = File(exists=True, desc='Arterial transit time from BASIL')
    abv_basil = File(exists=True, desc='Arterial blood volume from BASIL')
    seg_file = File(
        exists=True,
        mandatory=True,
        desc='Segmentation file, with 1 for GM, 2 for WM and 3 for CSF',
    )


class MultiPLDByTissueTypePlot(SimpleInterface):
    """Plot the distributions of the multi-delay parameters by tissue type, as ridge plots.

    There is one panel per parameter (ATT, aBAT and aBV), with one row per map of that
    parameter (the basic model and, for ATT and aBV, BASIL). Maps that are not defined are left
    out, and so is a panel with no maps.
    Voxels that are exactly zero or not finite are left out of the distributions.
    """

    input_spec = _MultiPLDByTissueTypePlotInputSpec
    output_spec = _ByTissueTypePlotOutputSpec

    def _run_interface(self, runtime):
        from aslprep.utils.plotting import plot_tissue_ridgeplot

        layout = {
            'Arterial Transit Time (s)': {'Basic': 'att', 'BASIL': 'att_basil'},
            'Arterial Bolus Arrival Time (s)': {'Basic': 'abat'},
            'Arterial Blood Volume (fraction)': {'Basic': 'abv', 'BASIL': 'abv_basil'},
        }
        seg = np.asanyarray(nb.load(self.inputs.seg_file).dataobj)
        panels = {}
        for x_label, row_fields in layout.items():
            rows = {
                row_label: _values_by_tissue(getattr(self.inputs, field), seg)
                for row_label, field in row_fields.items()
                if isdefined(getattr(self.inputs, field))
            }
            if rows:
                panels[x_label] = rows

        if not panels:
            raise ValueError('None of the multi-delay parameter maps are defined.')

        self._results['out_file'] = fname_presuffix(
            self.inputs.seg_file,
            suffix='_multipldbytissue.svg',
            use_ext=False,
            newpath=runtime.cwd,
        )
        plot_tissue_ridgeplot(panels, out_file=self._results['out_file'])
        return runtime


class _ATTStatusPlotInputSpec(BaseInterfaceInputSpec):
    att = File(exists=True, mandatory=True, desc='Arterial transit time map.')
    earliest_delay = File(
        exists=True,
        mandatory=True,
        desc='Earliest delay of each voxel, including slice timing.',
    )
    att_status = File(exists=True, mandatory=True, desc='ATT status map.')
    ref_vol = File(exists=True, mandatory=True, desc='ASL reference image.')
    mask = File(exists=True, mandatory=True, desc='Brain mask.')
    seg_file = File(
        exists=True,
        desc='Tissue segmentation in the same space (1: GM, 2: WM, 3: CSF).',
    )
    slice_axis = traits.Int(2, usedefault=True, desc='Axis along which slices were acquired.')
    n_cuts = traits.Int(5, usedefault=True, desc='Number of cuts to show.')


class _ATTStatusPlotOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='The figure.')


class ATTStatusPlot(SimpleInterface):
    """Plot where the multi-delay fit located the arterial transit time.

    The image rows are cut perpendicular to the acquired slices (sagittal for axial slices),
    so that a slice-timing effect appears as bands along the slice direction.
    The top row shows the ATT status (:data:`~aslprep.utils.cbf.ATT_STATUS`),
    and the bottom row shows ATT only where it was located, with other brain voxels in gray.
    The right-hand panel shows the fraction of brain voxels with each status in every slice,
    together with the earliest delay of that slice, including its slice-timing offset.
    """

    input_spec = _ATTStatusPlotInputSpec
    output_spec = _ATTStatusPlotOutputSpec

    def _run_interface(self, runtime):
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
        from matplotlib.patches import Patch
        from nibabel.affines import apply_affine
        from nilearn import plotting

        from aslprep.utils.cbf import ATT_STATUS

        axis = self.inputs.slice_axis
        status_img = nb.load(self.inputs.att_status)
        ref_img = nb.load(self.inputs.ref_vol)
        mask_img = nb.load(self.inputs.mask)
        affine = status_img.affine
        att, earliest_delay, status, mask = (
            nb.load(fname).get_fdata()
            for fname in (
                self.inputs.att,
                self.inputs.earliest_delay,
                self.inputs.att_status,
                self.inputs.mask,
            )
        )
        status = np.round(status).astype(int)
        mask = mask > 0
        located = mask & (status == ATT_STATUS['located'])

        labels = {
            ATT_STATUS['located']: ('located', '#2da467'),
            ATT_STATUS['bounded_above']: ('bounded above only', '#1b60a5'),
            ATT_STATUS['wide_range']: ('other wide range', '#c0392b'),
            ATT_STATUS['no_tissue_signal']: ('no tissue signal', '#9d8f25'),
            ATT_STATUS['not_fitted']: ('not fitted', '#555555'),
        }
        n_levels = max(labels) + 1
        status_cmap = ListedColormap([labels[level][1] for level in range(n_levels)])

        # Cut perpendicular to the slices: sagittal, unless the slices are themselves sagittal
        slice_code = nb.aff2axcodes(affine)[axis]
        display_mode = 'y' if slice_code in ('L', 'R') else 'x'
        world = apply_affine(affine, np.argwhere(mask))
        coord = world[:, 'xyz'.index(display_mode)]
        cuts = np.percentile(coord, np.linspace(15, 85, self.inputs.n_cuts))

        fig = plt.figure(figsize=(20, 8))
        grid = fig.add_gridspec(
            2,
            2,
            width_ratios=[3.2, 1],
            wspace=0.15,
            hspace=0.25,
            left=0.03,
            right=0.98,
            top=0.88,
            bottom=0.2,
        )
        common = {
            'bg_img': ref_img,
            'display_mode': display_mode,
            'cut_coords': cuts,
            'figure': fig,
            'black_bg': False,
            'draw_cross': False,
            'annotate': True,
        }

        # Status, shifted up by one so that background voxels (zero) are not drawn
        shifted = nb.Nifti1Image(np.where(mask, status + 1, 0).astype(np.uint8), affine)
        plotting.plot_roi(
            shifted,
            axes=fig.add_subplot(grid[0, 0]),
            title='ATT status',
            cmap=status_cmap,
            vmin=0.5,
            vmax=n_levels + 0.5,
            alpha=0.9,
            # The legend next to the per-slice fractions names the colors
            colorbar=False,
            **common,
        )

        located_att = nb.Nifti1Image(np.where(located, att, 0).astype(np.float32), affine)
        display = plotting.plot_stat_map(
            located_att,
            axes=fig.add_subplot(grid[1, 0]),
            title='Located ATT (s)',
            cmap='viridis',
            vmin=0,
            vmax=3,
            symmetric_cbar=False,
            threshold=1e-6,
            colorbar=True,
            resampling_interpolation='nearest',
            **common,
        )
        display.add_overlay(
            nb.Nifti1Image((mask & ~located).astype(np.uint8), affine),
            cmap=ListedColormap(['#bdbdbd']),
            threshold=0.5,
        )

        # Fraction of brain voxels with each status, per slice
        status, earliest_delay, mask = (
            np.moveaxis(arr, axis, -1) for arr in (status, earliest_delay, mask)
        )
        in_brain = np.flatnonzero(mask.any(axis=(0, 1)))
        ax = fig.add_subplot(grid[:, 1])
        n_brain = mask.sum(axis=(0, 1))
        left = np.zeros(in_brain.size)
        for level, (label, color) in labels.items():
            fraction = (mask & (status == level)).sum(axis=(0, 1))[in_brain] / n_brain[in_brain]
            ax.barh(in_brain, fraction, left=left, height=1, color=color, label=label)
            left += fraction
        ax.set_xlim(0, 1)
        ax.set_ylim(in_brain[0] - 0.5, in_brain[-1] + 0.5)
        ax.set_xlabel('Fraction of brain voxels')
        ax.set_ylabel('Slice')

        # ATTs shorter than a slice's earliest delay cannot be told apart
        earliest = np.array([np.median(earliest_delay[..., k][mask[..., k]]) for k in in_brain])
        delay_ax = ax.twiny()
        delay_ax.plot(earliest, in_brain, color='black', marker='.', label='earliest delay')
        delay_ax.set_xlabel('Earliest delay, including slice timing (s)')
        handles = [Patch(color=color, label=label) for label, color in labels.values()]
        handles += delay_ax.get_legend_handles_labels()[0]
        ax.legend(
            handles=handles,
            loc='upper center',
            bbox_to_anchor=(0.5, -0.1),
            ncol=2,
            fontsize=8,
            frameon=False,
        )

        # Fractions located over the whole brain and by tissue
        title = f'ATT located in {100 * located.sum() / max(mask.sum(), 1):.1f}% of brain voxels'
        if isdefined(self.inputs.seg_file):
            seg = nb.load(self.inputs.seg_file).get_fdata()
            brain = mask_img.get_fdata() > 0
            tissue_str = []
            for name, value in (('GM', 1), ('WM', 2)):
                tissue = brain & (np.round(seg) == value)
                if tissue.any():
                    tissue_str.append(f'{name} {100 * located[tissue].mean():.1f}%')

            if tissue_str:
                title += f' ({", ".join(tissue_str)})'

        fig.suptitle(title)

        self._results['out_file'] = fname_presuffix(
            self.inputs.att_status, suffix='_attstatusplot.svg', use_ext=False, newpath=runtime.cwd
        )
        fig.savefig(self._results['out_file'])
        plt.close(fig)
        return runtime


class _ExternalATTPlotInputSpec(BaseInterfaceInputSpec):
    att = File(
        exists=True,
        mandatory=True,
        desc="External ATT map on this run's ASL reference grid, including filled voxels.",
    )
    att_status = File(
        exists=True,
        mandatory=True,
        desc='External ATT status map (see aslprep.utils.att.EXTERNAL_ATT_STATUS).',
    )
    ref_vol = File(exists=True, mandatory=True, desc='ASL reference image.')
    mask = File(exists=True, mandatory=True, desc='Brain mask.')
    gm_pv = File(exists=True, mandatory=True, desc='Gray matter partial volume map.')
    wm_pv = File(exists=True, mandatory=True, desc='White matter partial volume map.')
    source_name = traits.Str(mandatory=True, desc='The multi-delay ATT map the ATT came from.')
    n_cuts = traits.Int(5, usedefault=True, desc='Number of cuts to show.')


class _ExternalATTPlotOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='The figure.')


class ExternalATTPlot(SimpleInterface):
    """Plot the arterial transit time a single-delay run took from a multi-delay run.

    The ATT is not estimated from the single-delay data: it is resampled from the multi-delay
    run's map, and filled where that run did not locate it.
    The top row shows the status of each brain voxel
    (:data:`~aslprep.utils.att.EXTERNAL_ATT_STATUS`),
    the middle row the ATT the multi-delay run measured,
    and the bottom row the filled ATT, with other brain voxels in gray.
    The right-hand panel shows the fraction of brain, gray matter and white matter voxels with
    each status.
    """

    input_spec = _ExternalATTPlotInputSpec
    output_spec = _ExternalATTPlotOutputSpec

    def _run_interface(self, runtime):
        import os

        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
        from matplotlib.patches import Patch
        from nibabel.affines import apply_affine
        from nilearn import plotting

        status_img = nb.load(self.inputs.att_status)
        ref_img = nb.load(self.inputs.ref_vol)
        affine = status_img.affine
        att, status, mask, gm_pv, wm_pv = (
            np.squeeze(nb.load(fname).get_fdata())
            for fname in (
                self.inputs.att,
                self.inputs.att_status,
                self.inputs.mask,
                self.inputs.gm_pv,
                self.inputs.wm_pv,
            )
        )
        status = np.round(status).astype(int)
        mask = mask > 0

        # Short labels of aslprep.utils.att.EXTERNAL_ATT_STATUS
        labels = {
            0: ('no ATT', '#555555'),
            1: ('measured by the multi-delay run', '#2da467'),
            2: ('filled, clipped into its interval', '#1b60a5'),
            3: ('filled, without an interval', '#7fb3e0'),
            4: ('not filled (--att-fill none)', '#bdbdbd'),
            5: ('not filled: too little GM and WM', '#9d8f25'),
            6: ("no tissue signal at this run's delay", '#c0392b'),
            7: ('not filled: no GM or WM median', '#8e44ad'),
        }
        n_levels = max(labels) + 1
        status_cmap = ListedColormap([labels[level][1] for level in range(n_levels)])

        # Sagittal cuts, which show how the ATT changes along the superior-inferior axis
        world = apply_affine(affine, np.argwhere(mask))
        cuts = np.percentile(world[:, 0], np.linspace(15, 85, self.inputs.n_cuts))

        fig = plt.figure(figsize=(20, 11))
        grid = fig.add_gridspec(
            3,
            2,
            width_ratios=[3.2, 1],
            wspace=0.15,
            hspace=0.25,
            left=0.03,
            right=0.98,
            top=0.9,
            bottom=0.06,
        )
        common = {
            'bg_img': ref_img,
            'display_mode': 'x',
            'cut_coords': cuts,
            'figure': fig,
            'black_bg': False,
            'draw_cross': False,
            'annotate': True,
        }

        # Status, shifted up by one so that background voxels (zero) are not drawn
        shifted = nb.Nifti1Image(np.where(mask, status + 1, 0).astype(np.uint8), affine)
        plotting.plot_roi(
            shifted,
            axes=fig.add_subplot(grid[0, 0]),
            title='External ATT status',
            cmap=status_cmap,
            vmin=0.5,
            vmax=n_levels + 0.5,
            alpha=0.9,
            # The legend next to the fractions names the colors
            colorbar=False,
            **common,
        )

        rows = (
            ('ATT measured by the multi-delay run (s)', (1,)),
            ('ATT filled by tissue median (s)', (2, 3)),
        )
        for row, (title, codes) in enumerate(rows, start=1):
            shown = mask & np.isin(status, codes)
            display = plotting.plot_stat_map(
                nb.Nifti1Image(np.where(shown, att, 0).astype(np.float32), affine),
                axes=fig.add_subplot(grid[row, 0]),
                title=title,
                cmap='viridis',
                vmin=0,
                vmax=3,
                symmetric_cbar=False,
                threshold=1e-6,
                colorbar=True,
                resampling_interpolation='nearest',
                **common,
            )
            if (mask & ~shown).any():
                display.add_overlay(
                    nb.Nifti1Image((mask & ~shown).astype(np.uint8), affine),
                    cmap=ListedColormap(['#bdbdbd']),
                    threshold=0.5,
                )

        # Fraction of the voxels of the brain and of each tissue with each status
        groups = {
            'Brain': mask,
            'GM (PV >= 0.5)': mask & (gm_pv >= 0.5),
            'WM (PV >= 0.5)': mask & (wm_pv >= 0.5),
        }
        # The legend goes below the bars, beside the bottom row
        ax = fig.add_subplot(grid[:2, 1])
        positions = np.arange(len(groups))[::-1]
        left = np.zeros(len(groups))
        for level, (label, color) in labels.items():
            fraction = np.array(
                [(sel & (status == level)).sum() / max(sel.sum(), 1) for sel in groups.values()]
            )
            ax.barh(positions, fraction, left=left, height=0.6, color=color, label=label)
            left += fraction
        ax.set_yticks(positions, labels=list(groups))
        ax.set_xlim(0, 1)
        ax.set_xlabel('Fraction of voxels')
        ax.legend(
            handles=[Patch(color=color, label=label) for label, color in labels.values()],
            loc='upper center',
            bbox_to_anchor=(0.5, -0.12),
            ncol=1,
            fontsize=9,
            frameon=False,
        )

        source = os.path.basename(self.inputs.source_name.split(':')[-1])
        measured = mask & (status == 1)
        fig.suptitle(
            f'Arterial transit time taken from the multi-delay run ({source}), not estimated '
            f'from this run.\nMeasured there in {100 * measured.sum() / max(mask.sum(), 1):.1f}% '
            'of brain voxels; filled values are imputed.'
        )

        self._results['out_file'] = fname_presuffix(
            self.inputs.att_status,
            suffix='_externalattplot.svg',
            use_ext=False,
            newpath=runtime.cwd,
        )
        fig.savefig(self._results['out_file'])
        plt.close(fig)
        return runtime


class _PVCDiagnosticsPlotInputSpec(BaseInterfaceInputSpec):
    mean_cbf = File(exists=True, mandatory=True, desc='Uncorrected mean CBF.')
    gm_cbf = File(exists=True, mandatory=True, desc='Gray matter CBF after LR PVC.')
    gm_valid = File(exists=True, mandatory=True, desc='Where gray matter CBF was estimated.')
    gm_pv = File(exists=True, mandatory=True, desc='Gray matter PV map.')
    status = File(exists=True, mandatory=True, desc='LR PVC fit status.')
    gm_noise = File(exists=True, mandatory=True, desc='Gray matter noise amplification.')
    summary = traits.Dict(mandatory=True, desc='Fraction of brain voxels with each status.')
    pv_metadata = traits.Dict(mandatory=True, desc='How the PV maps were built.')
    qc = traits.Dict(mandatory=False, desc='Scalar QC measures, including the CSF estimate.')


class _PVCDiagnosticsPlotOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='The figure.')


class PVCDiagnosticsPlot(SimpleInterface):
    """Plot local-regression partial volume correction diagnostics for the report.

    The panels show the gray matter PV map, the fit status, the gray matter noise
    amplification, and the distributions of uncorrected and corrected CBF where gray matter
    CBF was estimated.
    """

    input_spec = _PVCDiagnosticsPlotInputSpec
    output_spec = _PVCDiagnosticsPlotOutputSpec

    def _run_interface(self, runtime):
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap

        mean_cbf = nb.load(self.inputs.mean_cbf).get_fdata()
        gm_cbf = nb.load(self.inputs.gm_cbf).get_fdata()
        gm_valid = nb.load(self.inputs.gm_valid).get_fdata() > 0
        gm_pv = nb.load(self.inputs.gm_pv).get_fdata()
        status = nb.load(self.inputs.status).get_fdata()
        gm_noise = nb.load(self.inputs.gm_noise).get_fdata()

        # Show the axial slice with the most brain voxels
        k = int(np.argmax((status > 0).sum(axis=(0, 1))))

        summary = self.inputs.summary
        meta = self.inputs.pv_metadata
        fwhm = ' x '.join(f'{v:.2g}' for v in meta.get('EffectiveResolution', []))
        status_str = ', '.join(f'{name}: {100 * frac:.1f}%' for name, frac in summary.items())
        title = (
            f'Effective resolution {fwhm} mm FWHM ({meta.get("EffectiveResolutionSource")}). '
            f'Brain voxels: {status_str}'
        )
        if isdefined(self.inputs.qc) and self.inputs.qc.get('lr_csf_n_voxels'):
            qc = self.inputs.qc
            title += (
                f'\nApparent CSF CBF (QC, expected near 0): median {qc["lr_csf_cbf_median"]:.1f}'
                f', IQR {qc["lr_csf_cbf_iqr"]:.1f} mL/100 g/min '
                f'in {qc["lr_csf_n_voxels"]} voxels with at least 50% CSF'
            )

        self._results['out_file'] = fname_presuffix(
            self.inputs.gm_cbf, suffix='_pvcplot.svg', use_ext=False, newpath=runtime.cwd
        )

        # Other plots change the global matplotlib settings (nireports' carpet plot sets
        # font.size to 24), and nipype workers run many nodes in one process, so the layout
        # is only predictable with the default style.
        with plt.style.context('default'):
            fig, axes = plt.subplots(1, 4, figsize=(20, 5.5))
            im = axes[0].imshow(np.rot90(gm_pv[:, :, k]), cmap='gray', vmin=0, vmax=1)
            axes[0].set_title('Gray matter fraction (effective resolution)')
            fig.colorbar(im, ax=axes[0], fraction=0.046)

            status_cmap = ListedColormap(['black', '#1b60a5', '#2da467', '#9d8f25', '#c0392b'])
            im = axes[1].imshow(
                np.rot90(status[:, :, k]),
                cmap=status_cmap,
                vmin=-0.5,
                vmax=4.5,
                interpolation='nearest',
            )
            axes[1].set_title('Fit status')
            cbar = fig.colorbar(im, ax=axes[1], fraction=0.046, ticks=range(5))
            cbar.ax.set_yticklabels(['outside', 'GM + WM', 'GM only', 'WM only', 'none'])

            noise = np.where(gm_valid, gm_noise, np.nan)
            finite = noise[np.isfinite(noise)]
            vmax = float(np.percentile(finite, 99)) if finite.size else 1.0
            im = axes[2].imshow(np.rot90(noise[:, :, k]), cmap='magma', vmin=0, vmax=vmax)
            axes[2].set_title('GM noise amplification')
            fig.colorbar(im, ax=axes[2], fraction=0.046)
            for ax in axes[:3]:
                ax.axis('off')

            if gm_valid.any():
                values = np.concatenate([mean_cbf[gm_valid], gm_cbf[gm_valid]])
                low, high = np.percentile(values, [1, 99])
                bins = np.linspace(low, high, 60)
                axes[3].hist(mean_cbf[gm_valid], bins=bins, alpha=0.6, label='uncorrected')
                axes[3].hist(gm_cbf[gm_valid], bins=bins, alpha=0.6, label='GM corrected')
                axes[3].legend()

            axes[3].set_xlabel('CBF (mL/100 g/min)')
            axes[3].set_title('Where GM CBF was estimated')
            fig.suptitle(title, fontsize=11)
            fig.tight_layout()
            fig.savefig(self._results['out_file'])
            plt.close(fig)

        return runtime
