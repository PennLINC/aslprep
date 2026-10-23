"""Plotting interfaces."""

import nibabel as nb
import numpy as np
import pandas as pd
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
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
    in_file = File(exists=True, mandatory=True, desc='')
    seg_file = File(exists=True, mandatory=True, desc='Segmentation file')
    img_type = traits.Enum(
        'cbf',
        'att',
        'abat',
        'abv',
        mandatory=True,
        desc='Image type',
    )


class _CBFByTissueTypePlotOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='written file path')


class CBFByTissueTypePlot(SimpleInterface):
    """Prepare an CBF summary plot for the report."""

    input_spec = _CBFByTissueTypePlotInputSpec
    output_spec = _CBFByTissueTypePlotOutputSpec

    def _run_interface(self, runtime):
        import matplotlib.pyplot as plt
        import seaborn as sns
        from nilearn import image, masking

        self._results['out_file'] = fname_presuffix(
            self.inputs.in_file,
            suffix=f'_{self.inputs.img_type}plot.svg',
            use_ext=False,
            newpath=runtime.cwd,
        )

        column_names = {
            'cbf': 'Cerebral Blood Flow\n(mL/100 g/min)',
            'att': 'Arterial Transit Time\n(s)',
            'abat': 'Arterial Bolus Arrival Time\n(s)',
            'abv': 'Arterial Blood Volume\n(fraction)',
        }
        unit_str = column_names[self.inputs.img_type]

        dfs = []
        for i_tissue_type, tissue_type in enumerate(['GM', 'WM', 'CSF']):
            tissue_type_val = i_tissue_type + 1
            mask_img = image.math_img(
                f'(img == {tissue_type_val}).astype(int)',
                img=self.inputs.seg_file,
            )
            tissue_type_vals = masking.apply_mask(self.inputs.in_file, mask_img)
            df = pd.DataFrame(
                columns=[unit_str, 'Tissue Type'],
                data=list(
                    map(
                        list,
                        zip(
                            *[tissue_type_vals, [tissue_type] * tissue_type_vals.size],
                            strict=False,
                        ),
                    ),
                ),
            )
            dfs.append(df)

        df = pd.concat(dfs, axis=0)

        # Create the plot
        with sns.axes_style('whitegrid'), sns.plotting_context(font_scale=3):
            fig, ax = plt.subplots(figsize=(16, 8))
            sns.despine(ax=ax, bottom=True, left=True)
            sns.boxenplot(
                y=unit_str,
                data=df,
                width=0.6,
                showfliers=True,
                palette={'GM': '#1b60a5', 'WM': '#2da467', 'CSF': '#9d8f25'},
                hue='Tissue Type',
                legend=True,
                ax=ax,
            )
            fig.tight_layout()
            fig.savefig(self._results['out_file'])
            fig.clf()

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
        fig.suptitle(title, fontsize=11)
        fig.tight_layout()

        self._results['out_file'] = fname_presuffix(
            self.inputs.gm_cbf, suffix='_pvcplot.svg', use_ext=False, newpath=runtime.cwd
        )
        fig.savefig(self._results['out_file'])
        plt.close(fig)
        return runtime


class _BackgroundSuppressionPlotInputSpec(BaseInterfaceInputSpec):
    asl_file = File(exists=True, mandatory=True, desc='ASL file, as acquired.')
    aslcontext = File(exists=True, mandatory=True, desc='aslcontext file of asl_file.')
    asl_mask = File(exists=True, mandatory=True, desc='Brain mask in ASL space.')
    gm_tpm = File(exists=True, mandatory=True, desc='Gray matter probability map in ASL space.')
    csf_tpm = File(exists=True, mandatory=True, desc='CSF probability map in ASL space.')


class _BackgroundSuppressionPlotOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='The figure.')


class BackgroundSuppressionPlot(SimpleInterface):
    """Plot how close the control and label magnitudes are to the noise floor.

    The maps show the control signal over the noise floor and delta-M (control minus label)
    in three axial slices. The profiles show both against the CSF fraction of brain voxels:
    a dip of the signal to the floor at intermediate CSF fractions, with negative delta-M on
    the tissue side, means background suppression inverted the tissue (or the CSF)
    magnetization and the magnitude images flipped its sign.
    See :func:`~aslprep.utils.confounds.compute_background_suppression_qc`.
    """

    input_spec = _BackgroundSuppressionPlotInputSpec
    output_spec = _BackgroundSuppressionPlotOutputSpec

    def _run_interface(self, runtime):
        import matplotlib.pyplot as plt

        from aslprep.utils.confounds import (
            BS_NEAR_FLOOR_RATIO,
            compute_background_suppression_qc,
        )

        qc = compute_background_suppression_qc(
            asl_file=self.inputs.asl_file,
            aslcontext=self.inputs.aslcontext,
            mask=self.inputs.asl_mask,
            gm=self.inputs.gm_tpm,
            csf=self.inputs.csf_tpm,
        )
        metrics, profile, mask = qc['metrics'], qc['profile'], qc['mask']
        ratio, deltam = qc['signal_floor_ratio'], qc['deltam']

        fig = plt.figure(figsize=(18, 8))
        grid = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 1.6])
        slices = []
        if mask.any():
            z_extent = np.flatnonzero(mask.any(axis=(0, 1)))
            slices = [int(np.percentile(z_extent, q)) for q in (30, 50, 70)]

        ratio_vmax = 10
        dm_in_mask = deltam[mask & np.isfinite(deltam)]
        dm_vmax = float(np.percentile(np.abs(dm_in_mask), 98)) if dm_in_mask.size else 1.0
        for i_slice, k in enumerate(slices):
            ax_map = fig.add_subplot(grid[0, i_slice])
            im_ratio = ax_map.imshow(
                np.rot90(ratio[:, :, k]), cmap='magma', vmin=0, vmax=ratio_vmax
            )
            near_floor = mask[:, :, k] & (ratio[:, :, k] < BS_NEAR_FLOOR_RATIO)
            ax_map.imshow(
                np.rot90(np.where(near_floor, 1.0, np.nan)),
                cmap='cool',
                vmin=0,
                vmax=1,
                alpha=0.6,
            )
            ax_map.contour(np.rot90(mask[:, :, k]), levels=[0.5], colors='white', linewidths=0.5)
            ax_map.set_title(f'Control / noise floor, slice {k}')
            ax_map.axis('off')

            ax_dm_map = fig.add_subplot(grid[1, i_slice])
            im_dm = ax_dm_map.imshow(
                np.rot90(np.where(mask[:, :, k], deltam[:, :, k], np.nan)),
                cmap='RdBu_r',
                vmin=-dm_vmax,
                vmax=dm_vmax,
            )
            ax_dm_map.set_title(f'Delta-M (control - label), slice {k}')
            ax_dm_map.axis('off')

        if slices:
            fig.colorbar(im_ratio, ax=ax_map, fraction=0.046, extend='max')
            fig.colorbar(im_dm, ax=ax_dm_map, fraction=0.046)

        ax_ratio = fig.add_subplot(grid[0, 3])
        ax_dm = fig.add_subplot(grid[1, 3], sharex=ax_ratio)
        if not profile.empty:
            centers = (profile['csf_low'] + profile['csf_high']) / 2
            ax_ratio.fill_between(
                centers,
                profile['signal_floor_ratio_q25'],
                profile['signal_floor_ratio_q75'],
                alpha=0.3,
                color='#1b60a5',
            )
            ax_ratio.plot(centers, profile['signal_floor_ratio_median'], 'o-', color='#1b60a5')
            scale = profile['deltam_median'].abs().max()
            scale = scale if scale > 0 else 1.0
            ax_dm.fill_between(
                centers,
                profile['deltam_q25'] / scale,
                profile['deltam_q75'] / scale,
                alpha=0.3,
                color='#2da467',
            )
            ax_dm.plot(centers, profile['deltam_median'] / scale, 'o-', color='#2da467')
        ax_ratio.axhline(1, color='#c0392b', linestyle='--', label='noise floor')
        ax_ratio.axhline(
            BS_NEAR_FLOOR_RATIO, color='#9d8f25', linestyle=':', label='twice the noise floor'
        )
        ax_ratio.set_ylim(bottom=0)
        ax_ratio.set_ylabel('Control / noise floor\n(median and IQR)')
        ax_ratio.legend(loc='upper right')
        ax_ratio.set_title('Along the tissue-to-CSF partial volume axis')
        ax_dm.axhline(0, color='black', linewidth=0.8)
        ax_dm.set_ylabel('Delta-M / largest |median|\n(median and IQR)')
        ax_dm.set_xlabel('CSF fraction')
        ax_dm.set_xlim(0, 1)

        def _fmt(value, spec='.2f'):
            return 'n/a' if not np.isfinite(value) else format(value, spec)

        title = (
            f'Noise floor {_fmt(metrics["bs_noise_floor"], ".1f")}; '
            f'gray matter control / floor {_fmt(metrics["bs_gm_signal_floor_ratio"], ".1f")} '
            f'(delta-M retained {_fmt(metrics["bs_gm_deltam_retained"])}); '
            f'{_fmt(metrics["bs_perc_voxels_near_noise_floor"], ".1f")}% of brain voxels '
            f'below twice the floor.\nAlong CSF fraction: lowest control / floor '
            f'{_fmt(metrics["bs_csf_profile_signal_floor_ratio_min"])}, lowest delta-M '
            f'{_fmt(metrics["bs_csf_profile_deltam_min"])} (negative: label brighter than '
            'control, as when the magnetization is inverted at readout).'
        )
        fig.suptitle(title, fontsize=12)
        fig.tight_layout()

        self._results['out_file'] = fname_presuffix(
            self.inputs.asl_file, suffix='_bsplot.svg', use_ext=False, newpath=runtime.cwd
        )
        fig.savefig(self._results['out_file'])
        plt.close(fig)
        return runtime
