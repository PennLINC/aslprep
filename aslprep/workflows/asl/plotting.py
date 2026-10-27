# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for plotting ASLPrep derivatives."""

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow

from aslprep import config
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.bids import DerivativesDataSink
from aslprep.interfaces.confounds import GatherCBFConfounds
from aslprep.interfaces.plotting import (
    CBFByTissueTypePlot,
    CBFSummaryPlot,
    MultiPLDByTissueTypePlot,
    PVCDiagnosticsPlot,
)
from aslprep.interfaces.reports import CBFSummary
from aslprep.utils.misc import _select_last_in_list
from aslprep.workflows.asl.confounds import init_carpetplot_wf


def init_cbf_reporting_wf(
    metadata,
    plot_timeseries=True,
    scorescrub=False,
    basil=False,
    is_multi_pld=False,
    lr_pvc=False,
    name='cbf_reporting_wf',
):
    """Generate CBF reports.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.plotting import init_cbf_reporting_wf

            wf = init_cbf_reporting_wf(
                metadata={
                    "RepetitionTime": 4,
                    "RepetitionTimePreparation": 4,
                },
            )

    Parameters
    ----------
    metadata : :obj:`dict`
        ASL metadata.
    plot_timeseries : :obj:`bool`
        Whether to plot the CBF time series carpet plot.
    scorescrub : :obj:`bool`
        Whether SCORE/SCRUB outputs are available to plot.
    basil : :obj:`bool`
        Whether BASIL outputs are available to plot.
    is_multi_pld : :obj:`bool`
        Whether the data have several delays (plots ATT, aBAT and aBV).
    lr_pvc : :obj:`bool`
        Whether local-regression partial volume correction ran. If so, the corrected gray
        matter CBF and a diagnostics figure are added to the report.
    name : :obj:`str`
        Workflow name.

    Inputs
    ------
    aslref, asl_mask, t1w_dseg, aslref2anat_xfm, std2anat_xfm
        Reference image, brain mask, anatomical segmentation and transforms for the plots.
    confounds_file, qc_file
        Confounds and CBF QC tables.
    crown_mask, acompcor_masks
        Masks for the carpet plot (with ``plot_timeseries``).
    mean_cbf, cbf_ts, cifti_cbf_ts
        Mean CBF, volume-wise CBF (single-delay data) and its CIFTI version (with CIFTI
        outputs).
    att, abat, abv, abv_decayed
        Multi-delay outputs (with ``is_multi_pld``).
    cbf_ts_score, mean_cbf_score, mean_cbf_scrub, score_outlier_index
        SCORE/SCRUB outputs (with ``scorescrub``).
    mean_cbf_basil, mean_cbf_gm_basil, mean_cbf_wm_basil, att_basil, abv_basil
        BASIL outputs (with ``basil``). ATT and aBV are only plotted with ``is_multi_pld``.
    cbf_sd_basil, cbf_gm_sd_basil, cbf_wm_sd_basil, att_sd_basil
        BASIL posterior standard deviations (accepted, not plotted).
    mean_cbf_gm_lr, gm_valid_lr, gm_pv, status_lr, gm_noise_lr
        Local-regression gray matter CBF, where it was estimated, the gray matter PV map, the
        fit status and the noise amplification (with ``lr_pvc``).
    pvc_summary, pvc_qc, pv_metadata
        The fit-status fractions, the scalar QC measures and how the PV maps were built, for
        the diagnostics figure's annotations (with ``lr_pvc``).

    Outputs
    -------
    None. The reportlets are written to the ``figures`` directory.
    """
    from niworkflows.interfaces.images import SignalExtraction

    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'aslref',
                'asl_mask',
                't1w_dseg',
                'aslref2anat_xfm',
                'std2anat_xfm',
                'confounds_file',
                'qc_file',
                # If plot_timeseries is True
                'crown_mask',
                'acompcor_masks',
                # CBF outputs
                'mean_cbf',
                # Single-delay outputs
                'cbf_ts',  # only for non-GE
                # If CIFTI is enabled
                'cifti_cbf_ts',
                # Multi-delay outputs
                'att',
                'abat',
                'abv',
                'abv_decayed',  # unused
                # SCORE/SCRUB outputs
                'cbf_ts_score',  # unused
                'mean_cbf_score',
                'mean_cbf_scrub',
                'score_outlier_index',
                # BASIL outputs
                'mean_cbf_basil',
                'mean_cbf_gm_basil',
                'mean_cbf_wm_basil',  # unused
                'att_basil',
                'abv_basil',
                'cbf_sd_basil',  # unused
                'cbf_gm_sd_basil',  # unused
                'cbf_wm_sd_basil',  # unused
                'att_sd_basil',  # unused
                # Local-regression PVC outputs
                'mean_cbf_gm_lr',
                'gm_valid_lr',
                'gm_pv',
                'status_lr',
                'gm_noise_lr',
                'pvc_summary',
                'pvc_qc',
                'pv_metadata',
            ],
        ),
        name='inputnode',
    )

    summary = pe.Node(
        CBFSummary(),
        name='summary',
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )
    workflow.connect([
        (inputnode, summary, [
            ('confounds_file', 'confounds_file'),
            ('qc_file', 'qc_file'),
        ])
    ])  # fmt:skip

    # Warp dseg file from T1w space to ASL reference space
    warp_t1w_dseg_to_aslref = pe.Node(
        ApplyTransforms(
            float=True,
            dimension=3,
            default_value=0,
            interpolation='GenericLabel',
            invert_transform_flags=[True],
            args='-v',
        ),
        name='warp_t1w_dseg_to_aslref',
    )
    workflow.connect([
        (inputnode, warp_t1w_dseg_to_aslref, [
            ('asl_mask', 'reference_image'),
            ('t1w_dseg', 'input_image'),
            ('aslref2anat_xfm', 'transforms'),
        ]),
    ])  # fmt:skip

    if plot_timeseries:
        # Global and segment regressors
        signals_class_labels = [
            'global_signal',
            'csf',
            'white_matter',
            'csf_wm',
        ]
        merge_rois = pe.Node(
            niu.Merge(2, ravel_inputs=True),
            name='merge_rois',
            run_without_submitting=True,
        )
        signals = pe.Node(
            SignalExtraction(class_labels=signals_class_labels),
            name='signals',
            mem_gb=2,
        )
        workflow.connect([
            (inputnode, merge_rois, [
                ('asl_mask', 'in1'),
                ('acompcor_masks', 'in2'),
            ]),
            (inputnode, signals, [('cbf_ts', 'in_file')]),
            (merge_rois, signals, [('out', 'label_files')]),
        ])  # fmt:skip

        # Time series are only available for non-GE data.
        # Create confounds file with SCORE index
        cbf_confounds = pe.Node(
            GatherCBFConfounds(),
            name='cbf_confounds',
        )
        workflow.connect([
            (inputnode, cbf_confounds, [('score_outlier_index', 'score')]),
            (signals, cbf_confounds, [('out_file', 'signals')]),
        ])  # fmt:skip

        carpetplot_wf = init_carpetplot_wf(
            mem_gb=2,
            confounds_list=[
                ('global_signal', None, 'GS'),
                ('csf', None, 'GSCSF'),
                ('white_matter', None, 'GSWM'),
            ]
            + ([('score_outlier_index', None, 'SCORE Index')] if scorescrub else []),
            metadata=metadata,
            cifti_output=False,
            suffix='cbf',
            name='cbf_carpetplot_wf',
        )
        workflow.connect([
            (inputnode, carpetplot_wf, [
                ('std2anat_xfm', 'inputnode.std2anat_xfm'),
                ('cbf_ts', 'inputnode.asl'),
                ('asl_mask', 'inputnode.asl_mask'),
                ('aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
                ('crown_mask', 'inputnode.crown_mask'),
                (('acompcor_masks', _select_last_in_list), 'inputnode.acompcor_mask'),
                ('cifti_cbf_ts', 'inputnode.cifti_asl'),
            ]),
            (cbf_confounds, carpetplot_wf, [('confounds_file', 'inputnode.confounds_file')]),
        ])  # fmt:skip

    cbf_summary = pe.Node(CBFSummaryPlot(label='cbf', vmax=100), name='cbf_summary', mem_gb=1)
    workflow.connect([
        (inputnode, cbf_summary, [
            ('mean_cbf', 'in_file'),
            ('aslref', 'ref_vol'),
        ]),
    ])  # fmt:skip

    ds_report_cbf = pe.Node(
        DerivativesDataSink(datatype='figures', desc='brain', suffix='cbf'),
        name='ds_report_cbf',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([(cbf_summary, ds_report_cbf, [('out_file', 'in_file')])])

    # One ridge plot of every CBF output by tissue type, one row per output.
    cbf_rows = [('Basic', 'mean_cbf')]
    if scorescrub:
        cbf_rows += [('SCORE', 'mean_cbf_score'), ('SCRUB', 'mean_cbf_scrub')]

    if basil:
        cbf_rows += [('BASIL', 'mean_cbf_basil'), ('BASIL GM-PVC', 'mean_cbf_gm_basil')]

    if lr_pvc:
        cbf_rows += [('LR GM-PVC', 'mean_cbf_gm_lr')]

    merge_cbf_maps = pe.Node(
        niu.Merge(len(cbf_rows)),
        name='merge_cbf_maps',
        run_without_submitting=True,
    )
    for i_row, (_, field) in enumerate(cbf_rows):
        workflow.connect([(inputnode, merge_cbf_maps, [(field, f'in{i_row + 1}')])])

    cbf_by_tt_plot = pe.Node(
        CBFByTissueTypePlot(labels=[label for label, _ in cbf_rows]),
        name='cbf_by_tt_plot',
        mem_gb=1,
    )
    workflow.connect([
        (merge_cbf_maps, cbf_by_tt_plot, [('out', 'in_files')]),
        (warp_t1w_dseg_to_aslref, cbf_by_tt_plot, [('output_image', 'seg_file')]),
    ])  # fmt:skip

    ds_report_cbf_by_tt = pe.Node(
        DerivativesDataSink(
            datatype='figures',
            desc='cbfByTissueType',
            suffix='cbf',
        ),
        name='ds_report_cbf_by_tt',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([(cbf_by_tt_plot, ds_report_cbf_by_tt, [('out_file', 'in_file')])])

    if is_multi_pld:
        # Limits for the different figures.
        # Make sure these match the hardcoded limits in the model-fitting function.
        lims = {
            'att': (0, 5),
            'abat': (0, 5),
            'abv': (0, 0.1),
        }
        for img_type in ['att', 'abat', 'abv']:
            img_summary = pe.Node(
                CBFSummaryPlot(
                    label=img_type,
                    vmin=lims[img_type][0],
                    vmax=lims[img_type][1],
                ),
                name=f'{img_type}_summary',
                mem_gb=1,
            )
            workflow.connect([
                (inputnode, img_summary, [
                    (img_type, 'in_file'),
                    ('aslref', 'ref_vol'),
                ]),
            ])  # fmt:skip

            ds_report_img = pe.Node(
                DerivativesDataSink(
                    datatype='figures',
                    desc='brain',
                    suffix=img_type,
                ),
                name=f'ds_report_{img_type}',
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            workflow.connect([(img_summary, ds_report_img, [('out_file', 'in_file')])])

        # One figure of ATT, aBAT and aBV by tissue type, one panel per parameter.
        # BASIL may not write its aBV map, in which case that row is left out.
        multipld_fields = ['att', 'abat', 'abv'] + (['att_basil', 'abv_basil'] if basil else [])
        multipld_by_tt_plot = pe.Node(
            MultiPLDByTissueTypePlot(),
            name='multipld_by_tt_plot',
            mem_gb=1,
        )
        workflow.connect([
            (inputnode, multipld_by_tt_plot, [(field, field) for field in multipld_fields]),
            (warp_t1w_dseg_to_aslref, multipld_by_tt_plot, [('output_image', 'seg_file')]),
        ])  # fmt:skip

        ds_report_multipld_by_tt = pe.Node(
            DerivativesDataSink(
                datatype='figures',
                desc='multiPLDByTissueType',
                suffix='asl',
            ),
            name='ds_report_multipld_by_tt',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([
            (multipld_by_tt_plot, ds_report_multipld_by_tt, [('out_file', 'in_file')]),
        ])  # fmt:skip

    if scorescrub:
        score_summary = pe.Node(
            CBFSummaryPlot(label='score', vmax=100),
            name='score_summary',
            mem_gb=1,
        )
        workflow.connect([
            (inputnode, score_summary, [
                ('mean_cbf_score', 'in_file'),
                ('aslref', 'ref_vol'),
            ]),
        ])  # fmt:skip

        ds_report_score = pe.Node(
            DerivativesDataSink(datatype='figures', desc='score', suffix='cbf', keep_dtype=True),
            name='ds_report_score',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(score_summary, ds_report_score, [('out_file', 'in_file')])])

        scrub_summary = pe.Node(
            CBFSummaryPlot(label='scrub', vmax=100),
            name='scrub_summary',
            mem_gb=1,
        )
        workflow.connect([
            (inputnode, scrub_summary, [
                ('mean_cbf_scrub', 'in_file'),
                ('aslref', 'ref_vol'),
            ]),
        ])  # fmt:skip

        ds_report_scrub = pe.Node(
            DerivativesDataSink(datatype='figures', desc='scrub', suffix='cbf', keep_dtype=True),
            name='ds_report_scrub',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(scrub_summary, ds_report_scrub, [('out_file', 'in_file')])])

    if basil:
        basil_summary = pe.Node(
            CBFSummaryPlot(label='basil', vmax=100),
            name='basil_summary',
            mem_gb=1,
        )
        workflow.connect([
            (inputnode, basil_summary, [
                ('mean_cbf_basil', 'in_file'),
                ('aslref', 'ref_vol'),
            ]),
        ])  # fmt:skip

        ds_report_basil = pe.Node(
            DerivativesDataSink(datatype='figures', desc='basil', suffix='cbf', keep_dtype=True),
            name='ds_report_basil',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(basil_summary, ds_report_basil, [('out_file', 'in_file')])])

        pvc_summary = pe.Node(
            CBFSummaryPlot(label='pvc', vmax=120),
            name='pvc_summary',
            mem_gb=1,
        )
        workflow.connect([
            (inputnode, pvc_summary, [
                ('mean_cbf_gm_basil', 'in_file'),
                ('aslref', 'ref_vol'),
            ]),
        ])  # fmt:skip

        ds_report_pvc = pe.Node(
            DerivativesDataSink(datatype='figures', desc='basilGM', suffix='cbf', keep_dtype=True),
            name='ds_report_pvc',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(pvc_summary, ds_report_pvc, [('out_file', 'in_file')])])

    if lr_pvc:
        lr_summary = pe.Node(
            CBFSummaryPlot(label='pvc', vmax=120),
            name='lr_summary',
            mem_gb=1,
        )
        workflow.connect([
            (inputnode, lr_summary, [
                ('mean_cbf_gm_lr', 'in_file'),
                ('aslref', 'ref_vol'),
            ]),
        ])  # fmt:skip

        ds_report_lr = pe.Node(
            DerivativesDataSink(datatype='figures', desc='lrGM', suffix='cbf', keep_dtype=True),
            name='ds_report_lr',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(lr_summary, ds_report_lr, [('out_file', 'in_file')])])

        lr_diagnostics = pe.Node(PVCDiagnosticsPlot(), name='lr_diagnostics', mem_gb=1)
        workflow.connect([
            (inputnode, lr_diagnostics, [
                ('mean_cbf', 'mean_cbf'),
                ('mean_cbf_gm_lr', 'gm_cbf'),
                ('gm_valid_lr', 'gm_valid'),
                ('gm_pv', 'gm_pv'),
                ('status_lr', 'status'),
                ('gm_noise_lr', 'gm_noise'),
                ('pvc_summary', 'summary'),
                ('pvc_qc', 'qc'),
                ('pv_metadata', 'pv_metadata'),
            ]),
        ])  # fmt:skip

        ds_report_lr_diagnostics = pe.Node(
            DerivativesDataSink(datatype='figures', desc='lr', suffix='pvc', keep_dtype=True),
            name='ds_report_lr_diagnostics',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([
            (lr_diagnostics, ds_report_lr_diagnostics, [('out_file', 'in_file')]),
        ])  # fmt:skip

    return workflow
