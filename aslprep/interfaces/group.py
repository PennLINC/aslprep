# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Interfaces for calculating and collecting confounds."""

import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    isdefined,
    traits,
)


class _AggregateCBFQCInputSpec(BaseInterfaceInputSpec):
    in_files = traits.List(File(exists=True), desc='list of QC tsv files')


class _AggregateCBFQCOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='aggregated QC file')


class AggregateCBFQC(SimpleInterface):
    """Aggregate CBF QC metrics across runs, sessions, and subjects."""

    input_spec = _AggregateCBFQCInputSpec
    output_spec = _AggregateCBFQCOutputSpec

    def _run_interface(self, runtime):
        in_files = self.inputs.in_files
        dfs = []
        for in_file in in_files:
            temp_df = pd.read_table(in_file)
            temp_df['source_file'] = os.path.basename(in_file)
            dfs.append(temp_df)

        df = pd.concat(dfs)

        self._results['out_file'] = os.path.join(runtime.cwd, 'aggregated_qc.tsv')
        df.to_csv(self._results['out_file'], index=None, sep='\t')

        return runtime


class _PlotAggregatedCBFQCInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, desc='aggregated QC file')


class _PlotAggregatedCBFQCOutputSpec(TraitedSpec):
    out_files = traits.List(File(exists=True), desc='list of output figures')
    measures = traits.List(traits.String, desc='list of measures that have been plotted')


class PlotAggregatedCBFQC(SimpleInterface):
    """Plot aggregated CBF QC metrics across runs, sessions, and subjects."""

    input_spec = _PlotAggregatedCBFQCInputSpec
    output_spec = _PlotAggregatedCBFQCOutputSpec

    def _run_interface(self, runtime):
        from bids.layout import parse_file_entities

        df = pd.read_table(self.inputs.in_file)
        all_filenames = df['source_file'].unique()
        all_entities = [list(parse_file_entities(f).keys()) for f in all_filenames]
        # flatten list of lists
        all_entities = [item for sublist in all_entities for item in sublist]
        all_entities = list(set(all_entities))

        cols_to_ignore = all_entities + ['source_file']
        measures = [c for c in df.columns if c not in cols_to_ignore]

        entities_to_ignore = ['subject', 'session', 'run']
        entities_to_group = [e for e in all_entities if e not in entities_to_ignore]

        for measure in measures:
            figure_file = self.plot_aggregated_qc_metrics(
                df, measure, entities_to_group, runtime.cwd
            )
            self._results['out_files'].append(figure_file)
            self._results['measures'].append(measure)

        return runtime

    def plot_aggregated_qc_metrics(self, df, measure, entity_columns, cwd):
        """Plot aggregated QC metrics for a given measure and entity columns."""
        # Identify entities with varying values in dataset
        grouping_columns = []
        for col in entity_columns:
            if df[col].unique().size > 1:
                grouping_columns.append(col)
                df[col] = df[col].fillna('NONE')

        # Create a new column that indexes combinations of grouping columns
        df['group'] = ''
        for col in grouping_columns:
            df['group'] += col + '-' + df[col] + '_'

        # Remove trailing underscore
        df['group'] = df['group'].str.rstrip('_')

        kwargs = {}
        if len(grouping_columns):
            kwargs = {'hue': 'group'}

        fig, ax = plt.subplots()
        sns.stripplot(data=df, y=measure, ax=ax, figure=fig, **kwargs)
        figure_file = os.path.join(cwd, f'{measure}.svg')
        fig.savefig(figure_file)
        return figure_file


class _MeanCBFInputSpec(BaseInterfaceInputSpec):
    in_files = traits.List(File(exists=True), desc='list of CBF files')


class _MeanCBFOutputSpec(TraitedSpec):
    mean_map = File(exists=True, desc='mean CBF map')
    sd_map = File(exists=True, desc='SD CBF map')


class MeanCBF(SimpleInterface):
    """Calculate mean and SD CBF maps dataset."""

    input_spec = _MeanCBFInputSpec
    output_spec = _MeanCBFOutputSpec

    def _run_interface(self, runtime):
        from nilearn import image

        mean_map = image.mean_img(self.inputs.in_files)
        sd_map = image.math_img('np.std(img, axis=3)', img=self.inputs.in_files)

        mean_map_file = os.path.join(runtime.cwd, 'mean_cbf.nii.gz')
        sd_map_file = os.path.join(runtime.cwd, 'sd_cbf.nii.gz')

        mean_map.to_filename(mean_map_file)
        sd_map.to_filename(sd_map_file)

        self._results['mean_map'] = mean_map_file
        self._results['sd_map'] = sd_map_file

        return runtime


class _StatisticalMapInputSpecRPT(BaseInterfaceInputSpec):
    overlay = File(
        exists=True,
        mandatory=True,
        desc='FC inflation time series',
    )
    underlay = File(
        exists=True,
        mandatory=True,
        desc='Underlay image',
    )
    mask = File(
        exists=True,
        mandatory=False,
        desc='Mask image',
    )
    cmap = traits.Str(
        'viridis',
        desc='Colormap',
        usedefault=True,
    )
    out_report = File(
        'statistical_map_report.svg',
        usedefault=True,
        desc='Filename for the visual report generated by Nipype.',
    )


class _StatisticalMapOutputSpecRPT(TraitedSpec):
    out_report = File(
        exists=True,
        desc='Filename for the visual report generated by Nipype.',
    )


class StatisticalMapRPT(SimpleInterface):
    """Create a reportlet for Rapidtide outputs."""

    input_spec = _StatisticalMapInputSpecRPT
    output_spec = _StatisticalMapOutputSpecRPT

    def _run_interface(self, runtime):
        from uuid import uuid4

        from nilearn import image, masking, plotting
        from nireports._vendored.svgutils.transform import fromstring
        from nireports.reportlets.utils import compose_view, cuts_from_bbox, extract_svg

        out_file = os.path.abspath(self.inputs.out_report)

        if isdefined(self.inputs.mask):
            mask_img = image.load_img(self.inputs.mask)
            overlay_img = masking.unmask(
                masking.apply_mask(self.inputs.overlay, self.inputs.mask),
                self.inputs.mask,
            )
            # since the moving image is already in the fixed image space we
            # should apply the same mask
            underlay_img = image.load_img(self.inputs.underlay)
        else:
            overlay_img = image.load_img(self.inputs.overlay)
            underlay_img = image.load_img(self.inputs.underlay)
            mask_img = image.threshold_img(overlay_img, 1e-3)

        n_cuts = 7
        cuts = cuts_from_bbox(mask_img, cuts=n_cuts)
        order = ('z', 'x', 'y')
        out_files = []

        # Plot each cut axis
        plot_params = {}
        for mode in list(order):
            plot_params['display_mode'] = mode
            plot_params['cut_coords'] = cuts[mode]
            plot_params['title'] = None
            plot_params['cmap'] = self.inputs.cmap

            # Generate nilearn figure
            display = plotting.plot_stat_map(
                overlay_img,
                bg_img=underlay_img,
                **plot_params,
            )

            svg = extract_svg(display, compress=False)
            display.close()

            # Find and replace the figure_1 id.
            svg = svg.replace('figure_1', f'{mode}-{uuid4()}', 1)
            out_files.append(fromstring(svg))

        compose_view(bg_svgs=out_files, fg_svgs=None, out_file=out_file)
        self._results['out_report'] = out_file
        return runtime
