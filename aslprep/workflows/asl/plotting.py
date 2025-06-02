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
from aslprep.interfaces.plotting import CBFByTissueTypePlot, CBFSummaryPlot
from aslprep.interfaces.reports import CBFSummary
from aslprep.utils.misc import _select_last_in_list
from aslprep.workflows.asl.confounds import init_carpetplot_wf


def init_cbf_reporting_wf(
    metadata,
    plot_timeseries=True,
    scorescrub=False,
    basil=False,
    is_multi_pld=False,
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
                # SCORE/SCRUB outputs
                'cbf_ts_score',  # unused
                'mean_cbf_score',
                'mean_cbf_scrub',
                'score_outlier_index',
                # BASIL outputs
                'mean_cbf_basil',
                'mean_cbf_gm_basil',
                'mean_cbf_wm_basil',  # unused
                'att_basil',  # unused
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
        DerivativesDataSink(datatype='figures', desc='cbf', suffix='cbf', keep_dtype=True),
        name='ds_report_cbf',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([(cbf_summary, ds_report_cbf, [('out_file', 'in_file')])])

    cbf_by_tt_plot = pe.Node(
        CBFByTissueTypePlot(img_type='cbf'),
        name='cbf_by_tt_plot',
    )
    workflow.connect([
        (inputnode, cbf_by_tt_plot, [('mean_cbf', 'in_file')]),
        (warp_t1w_dseg_to_aslref, cbf_by_tt_plot, [('output_image', 'seg_file')]),
    ])  # fmt:skip

    ds_report_cbf_by_tt = pe.Node(
        DerivativesDataSink(
            datatype='figures',
            desc='cbfByTissueType',
            suffix='cbf',
            keep_dtype=True,
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
                    suffix=img_type,
                    keep_dtype=True,
                ),
                name=f'ds_report_{img_type}',
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            workflow.connect([(img_summary, ds_report_img, [('out_file', 'in_file')])])

            img_by_tt_plot = pe.Node(
                CBFByTissueTypePlot(img_type=img_type),
                name=f'{img_type}_by_tt_plot',
            )
            workflow.connect([
                (inputnode, img_by_tt_plot, [(img_type, 'in_file')]),
                (warp_t1w_dseg_to_aslref, img_by_tt_plot, [('output_image', 'seg_file')]),
            ])  # fmt:skip

            ds_report_img_by_tt = pe.Node(
                DerivativesDataSink(
                    datatype='figures',
                    desc=f'{img_type}ByTissueType',
                    suffix=img_type,
                    keep_dtype=True,
                ),
                name=f'ds_report_{img_type}_by_tt',
                run_without_submitting=True,
                mem_gb=config.DEFAULT_MEMORY_MIN_GB,
            )
            workflow.connect([(img_by_tt_plot, ds_report_img_by_tt, [('out_file', 'in_file')])])

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

        score_by_tt_plot = pe.Node(
            CBFByTissueTypePlot(img_type='cbf'),
            name='score_by_tt_plot',
        )
        workflow.connect([
            (inputnode, score_by_tt_plot, [('mean_cbf_score', 'in_file')]),
            (warp_t1w_dseg_to_aslref, score_by_tt_plot, [('output_image', 'seg_file')]),
        ])  # fmt:skip

        ds_report_score_by_tt = pe.Node(
            DerivativesDataSink(
                datatype='figures',
                desc='scoreByTissueType',
                suffix='cbf',
                keep_dtype=True,
            ),
            name='ds_report_score_by_tt',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(score_by_tt_plot, ds_report_score_by_tt, [('out_file', 'in_file')])])

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

        scrub_by_tt_plot = pe.Node(
            CBFByTissueTypePlot(img_type='cbf'),
            name='scrub_by_tt_plot',
        )
        workflow.connect([
            (inputnode, scrub_by_tt_plot, [('mean_cbf_scrub', 'in_file')]),
            (warp_t1w_dseg_to_aslref, scrub_by_tt_plot, [('output_image', 'seg_file')]),
        ])  # fmt:skip

        ds_report_scrub_by_tt = pe.Node(
            DerivativesDataSink(
                datatype='figures',
                desc='scrubByTissueType',
                suffix='cbf',
                keep_dtype=True,
            ),
            name='ds_report_scrub_by_tt',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(scrub_by_tt_plot, ds_report_scrub_by_tt, [('out_file', 'in_file')])])

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

        basil_by_tt_plot = pe.Node(
            CBFByTissueTypePlot(img_type='cbf'),
            name='basil_by_tt_plot',
        )
        workflow.connect([
            (inputnode, basil_by_tt_plot, [('mean_cbf_basil', 'in_file')]),
            (warp_t1w_dseg_to_aslref, basil_by_tt_plot, [('output_image', 'seg_file')]),
        ])  # fmt:skip

        ds_report_basil_by_tt = pe.Node(
            DerivativesDataSink(
                datatype='figures',
                desc='basilByTissueType',
                suffix='cbf',
                keep_dtype=True,
            ),
            name='ds_report_basil_by_tt',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(basil_by_tt_plot, ds_report_basil_by_tt, [('out_file', 'in_file')])])

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

        pvc_by_tt_plot = pe.Node(
            CBFByTissueTypePlot(img_type='cbf'),
            name='pvc_by_tt_plot',
        )
        workflow.connect([
            (inputnode, pvc_by_tt_plot, [('mean_cbf_gm_basil', 'in_file')]),
            (warp_t1w_dseg_to_aslref, pvc_by_tt_plot, [('output_image', 'seg_file')]),
        ])  # fmt:skip

        ds_report_pvc_by_tt = pe.Node(
            DerivativesDataSink(
                datatype='figures',
                desc='basilGMByTissueType',
                suffix='cbf',
                keep_dtype=True,
            ),
            name='ds_report_pvc_by_tt',
            run_without_submitting=True,
            mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([(pvc_by_tt_plot, ds_report_pvc_by_tt, [('out_file', 'in_file')])])

    return workflow
