# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for calculating CBF."""

from numbers import Number

from nipype.interfaces import utility as niu
from nipype.interfaces.fsl import Info
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.images import RobustAverage
from templateflow.api import get as get_template

from aslprep import config
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.bids import BIDSURI, DerivativesDataSink
from aslprep.interfaces.cbf import (
    BASILCBF,
    ComputeCBF,
    ExtractCBF,
    ScoreAndScrubCBF,
)
from aslprep.interfaces.parcellation import NiftiParcellate
from aslprep.utils.asl import (
    BOLUS_CUTOFF_TECHNIQUES,
    describe_labeling_efficiency,
    determine_multi_pld,
    estimate_labeling_efficiency,
    get_bolus_duration,
    get_inflow_times,
    infer_m0tr,
    pcasl_or_pasl,
    prepare_basil_kwargs,
)
from aslprep.utils.cbf import (
    M0_THRESHOLD,
    PARTITION_COEFFICIENT,
    REFERENCE_HAEMATOCRIT,
    estimate_t1,
    haematocrit_applies,
)
from aslprep.utils.pvc import check_lr_pvc_supported
from aslprep.workflows.asl.outputs import BASE_INPUT_FIELDS
from aslprep.workflows.asl.pv import init_asl_pv_wf, init_lr_pvc_wf, lr_pvc_fields


def _lr_pvc_mem_gb(name_source, num_threads):
    """Memory for multi-delay LR PVC, from the ASL series' shape (2 GB if it cannot be read)."""
    import nibabel as nb

    from aslprep.utils.pvc import estimate_lr_pvc_mem_gb

    try:
        shape = nb.load(name_source).shape
    except Exception:  # noqa: BLE001 - the name source may be a placeholder, as in tests
        return 2.0 + 0.1 * num_threads
    return estimate_lr_pvc_mem_gb(shape, num_threads)


def init_cbf_wf(
    *,
    name_source,
    processing_target,
    aslcontext,
    metadata,
    m0scan_metadata=None,
    scorescrub=False,
    basil=False,
    m0_scale=1,
    smooth_kernel=5,
    pvc_method='none',
    pvc_kernel_size=(5, 5, 1),
    asl_effective_resolution='nominal',
    name='cbf_wf',
):
    """Create a workflow for :abbr:`CCBF (compute cbf)`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            import json

            from aslprep.tests.tests import mock_config
            from aslprep import config
            from aslprep.workflows.asl.cbf import init_cbf_wf

            with mock_config():
                perf_dir = config.execution.bids_dir / "sub-01" / "perf"
                with open(perf_dir / "sub-01_asl.json", "r") as fobj:
                    metadata = json.load(fobj)

                wf = init_cbf_wf(
                    name_source=str(perf_dir / "sub-01_asl.nii.gz"),
                    processing_target="control",
                    metadata=metadata,
                    scorescrub=True,
                    basil=True,
                )

    Parameters
    ----------
    name_source : :obj:`str`
        Path to the raw ASL file.
    metadata : :obj:`dict`
        BIDS metadata for asl file.
        List-type fields in this dict have not been reduced based on selected volumes,
        so this should only be used for metadata that is not affected by volume selection.
    scorescrub
    basil
    m0_scale
    smooth_kernel
    pvc_method : {'none', 'local-regression'}
        Additional partial volume correction of the mean CBF (and SCORE) maps.
    pvc_kernel_size : :obj:`tuple` of :obj:`int`
        Local-regression kernel size, ordered (in-plane, in-plane, slice).
    asl_effective_resolution : :obj:`str` or :obj:`list` of :obj:`float`
        Effective resolution used to build the partial volume maps.
    name : :obj:`str`
        Name of workflow (default: ``cbf_wf``)

    Inputs
    ------
    asl_file
        asl series NIfTI file, after preprocessing
    aslcontext : :obj:`str`
    metadata : :obj:`dict`
        Metadata dictionary associated with asl_file, after preprocessing.
        List-type fields in this dict *have* been reduced based on selected volumes,
        so this is the dictionary that should be used to select volumes.
    m0scan : :obj:`str` or None
        An M0 scan (if available as a separate file) in aslref space.
    m0scan_metadata : :obj:`dict` or None
    asl_mask
        asl mask NIFTI file
    t1w_tpms
        t1w probability maps
    aslref2anat_xfm
        asl to t1w transformation file

    Outputs
    -------
    mean_cbf
        Mean CBF.
    cbf_ts
        Volume-wise CBF (single-delay data only).
    att, abat, abv
        Arterial transit time, arterial bolus arrival time and arterial blood volume
        (multi-delay data only).
    att_lower, att_upper, abat_lower, abat_upper
        The ranges of ATT and aBAT that fit as well as the estimates (multi-delay data only).
    abv_decayed
        aBV times the label decay over aBAT (multi-delay (P)CASL data only).
    plds
        The post-labeling delays used for quantification.
    cbf_ts_score, mean_cbf_score, mean_cbf_scrub, score_outlier_index
        SCORE and SCRUB outputs (with ``scorescrub``).
    mean_cbf_basil, mean_cbf_gm_basil, mean_cbf_wm_basil, att_basil, abv_basil
        BASIL outputs (with ``basil``); ATT and aBV for multi-delay data only.
    cbf_sd_basil, cbf_gm_sd_basil, cbf_wm_sd_basil, att_sd_basil
        BASIL posterior standard deviations (with ``basil``).
    gm_pv, wm_pv, csf_pv, pv_metadata
        Partial volume maps on the ASL reference grid and how they were built
        (with ``pvc_method='local-regression'``; see
        :func:`~aslprep.workflows.asl.pv.init_asl_pv_wf`).
    mean_cbf_gm_lr, ..., csf_noise_lr, pvc_summary, pvc_metadata, pvc_qc
        Local-regression partial volume correction outputs for ``mean_cbf`` (with
        ``pvc_method='local-regression'``; see :func:`~aslprep.workflows.asl.pv.init_lr_pvc_wf`
        for each field).
    mean_cbf_gm_lr_score, ..., csf_noise_lr_score, pvc_summary_score
        The same for ``mean_cbf_score`` (with ``pvc_method='local-regression'`` and
        ``scorescrub``).
    """
    workflow = Workflow(name=name)

    workflow.__desc__ = """

Cerebral blood flow computation and denoising

"""

    m0type = metadata['M0Type']
    m0tr = infer_m0tr(aslcontext=aslcontext, metadata=metadata, m0scan_metadata=m0scan_metadata)
    is_casl = pcasl_or_pasl(metadata=metadata)
    is_multi_pld = determine_multi_pld(metadata=metadata)
    if (processing_target == 'cbf') and not basil:
        config.loggers.workflow.info(f'Only CBF volumes are detected in {name_source}.')
    elif processing_target == 'cbf':
        config.loggers.workflow.warning(
            f'Only CBF volumes are detected in {name_source}. '
            'BASIL will automatically be disabled.'
        )
        basil = False

    if m0type in ('Included', 'Separate'):
        m0_str = (
            'Calibration (M0) volumes associated with the ASL scan were smoothed with a '
            f'Gaussian kernel (FWHM={smooth_kernel} mm) and the average calibration image was '
            f'calculated. '
        )
        if isinstance(m0tr, Number) and m0tr < 5:
            _, t1tissue = estimate_t1(metadata=metadata)
            m0_str += (
                f'As the TR of the calibration volumes was less than 5 seconds ({m0tr}), '
                'the average calibration volume was divided by 1 - exp(-TR / T1) to '
                'compensate for incomplete T1 recovery, as recommended in '
                f'@alsop_recommended_2015, using a gray matter T1 of {t1tissue:g} s '
                '[@wright2008water]. '
            )

        m0_str += f'The calibration volume was then scaled by {m0_scale}.'

    elif m0type == 'Estimate':
        m0_str = (
            f'A single M0 estimate of {metadata["M0Estimate"]} was used to produce a calibration '
            f"'image' and was scaled by {m0_scale}."
        )
    else:
        m0_str = (
            f'As no calibration images or provided M0 estimate was available for the ASL scan, '
            'the control volumes were used as a substitute. '
            'The control volumes in the ASL scans were smoothed with a '
            f'Gaussian kernel (FWHM={smooth_kernel} mm) and the average control image was '
            f'calculated and scaled by {m0_scale}.'
        )

    if m0type != 'Estimate':
        m0_str += (
            f' Voxels where it fell below {100 * M0_THRESHOLD:g}% of its 98th percentile were '
            'treated as having no signal, and were given no CBF estimate.'
        )

    if processing_target == 'cbf':
        workflow.__desc__ += """\
: *ASLPrep* loaded pre-calculated cerebral blood flow (CBF) data from the ASL file.
"""

    elif is_multi_pld:
        abv_str = ''
        if is_casl:
            abv_str = (
                '\nBecause the arterial signal only determines the product of aBV and the decay '
                'of the label over aBAT, the reported aBV depends on the reported aBAT.\n'
                'That product, which does not, was also written out.'
            )

        workflow.__desc__ += f"""\
: *ASLPrep* calculated cerebral blood flow (CBF) from the multi-delay
{metadata['ArterialSpinLabelingType']} data using a two-compartment general kinetic model (GKM)
[@buxton1998general], as recommended and extended in @woods2023recommendations.
This model contains separate terms for tissue and macrovascular delta-M signals.

{m0_str}
The equilibrium magnetization of arterial blood (M0b) was estimated voxel-wise as
M0 divided by the blood-brain partition coefficient,
and was used in both the tissue and the macrovascular terms of the GKM.

CBF, arterial transit time (ATT; the arrival time of the labeled blood in the tissue),
arterial bolus arrival time (aBAT; the arrival time of the labeled blood in the arteries
within each voxel), and arterial blood volume (aBV) were estimated with a bounded least squares
fit.
Because the model is linear in CBF and aBV for fixed ATT and aBAT,
ATT and aBAT were found by a grid search,
with CBF and aBV solved exactly for each candidate.
The data only constrain aBAT to an interval between consecutive labeled-bolus arrival and
departure times, so aBAT was reported as the midpoint of its interval.
Where a range of ATT values fit equally well
(for example, when ATT was shorter than every delay),
ATT was reported as a single value from that range.
The lower and upper bounds of the ATT and aBAT values that fit equally well were also calculated,
and are written out with the native-space outputs.{abv_str}{
            ' The fit was run on a GPU with CuPy.' if config.nipype.use_gpu else ''
        }
"""

    elif is_casl:
        # Single-delay (P)CASL data
        workflow.__desc__ += f"""\
: *ASLPrep* calculated cerebral blood flow (CBF) from the single-delay
{metadata['ArterialSpinLabelingType']} data using a one-compartment general kinetic model
[@buxton1998general].
{m0_str}
"""

    else:
        # Single-delay PASL data, with different bolus cut-off techniques
        bcut = metadata.get('BolusCutOffTechnique')
        bcut_name, bcut_citation = BOLUS_CUTOFF_TECHNIQUES.get(bcut, (bcut, None))
        # Sources of the CBF formula for each technique
        singlepld_pasl_strs = {
            'QUIPSS': '@wong1998quantitative',
            'QUIPSSII': '@alsop_recommended_2015',
            'Q2TIPS': '@alsop_recommended_2015 and @noguchi2015technical',
        }
        bcut_str = f'the {bcut_name} bolus cut-off technique'
        if bcut_citation:
            bcut_str += f' [{bcut_citation}]'

        workflow.__desc__ += f"""\
: *ASLPrep* calculated cerebral blood flow (CBF) from the single-delay PASL data
using a one-compartment general kinetic model [@buxton1998general]
modified for {bcut_str}, as described in {singlepld_pasl_strs.get(bcut, bcut_citation)}.
{m0_str}
"""

    if 'SliceTiming' in metadata:
        delay_name = 'post-labeling delays' if is_casl else 'inversion times'
        workflow.__desc__ += (
            f'Prior to calculating CBF, the {delay_name} were shifted on a slice-wise '
            'basis based on the slice timing.'
        )

    t1blood = None
    if processing_target != 'cbf':
        # Also passed to BASIL, so both models use the same value.
        t1blood, _ = estimate_t1(metadata=metadata)
        if haematocrit_applies(metadata):
            workflow.__desc__ += (
                f' The T1 of arterial blood was set to {t1blood:.3f} s: '
                f'the value for {metadata["MagneticFieldStrength"]:g} T, scaled for the '
                f"participant's measured haematocrit ({metadata['Haematocrit']:g}) relative to "
                f'a typical haematocrit of {REFERENCE_HAEMATOCRIT:g} [@lu2004determining], '
                'using the model of @hales2016general.'
            )
        else:
            field_strength = metadata['MagneticFieldStrength']
            t1blood_citation = (
                '@alsop_recommended_2015' if field_strength in (1.5, 3) else '@zhang2013vivo'
            )
            workflow.__desc__ += (
                f' The T1 of arterial blood was set to {t1blood:g} s, '
                f'the value for {field_strength:g} T [{t1blood_citation}].'
            )

        workflow.__desc__ += (
            ' CBF was calculated with a blood-brain partition coefficient of '
            f'{PARTITION_COEFFICIENT:g} mL/g [@alsop_recommended_2015] and '
            f'{describe_labeling_efficiency(metadata)}.'
        )

    # End the paragraph, so that the next workflow's description does not run into it.
    workflow.__desc__ += '\n'

    lr_pvc = pvc_method == 'local-regression'
    lr_sources = []
    lr_fields = []
    if lr_pvc:
        check_lr_pvc_supported(
            metadata=metadata,
            processing_target=processing_target,
            is_multi_pld=is_multi_pld,
            kernel_size=pvc_kernel_size,
        )
        lr_sources = ['mean_cbf'] + (['mean_cbf_score'] if scorescrub else [])
        lr_fields = ['gm_pv', 'wm_pv', 'csf_pv', 'pv_metadata', 'pvc_metadata', 'pvc_qc']
        lr_fields += lr_pvc_fields(lr_sources, multi_delay=is_multi_pld)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'asl_file',
                'aslcontext',
                'metadata',
                'm0scan',
                'm0scan_metadata',
                'asl_mask',
                't1w_tpms',
                'aslref2anat_xfm',
            ],
        ),
        name='inputnode',
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'mean_cbf',
                'cbf_ts',  # Only calculated for single-delay data
                'att',  # Only calculated for multi-delay data
                'abat',  # Only calculated for multi-delay data
                'abv',  # Only calculated for multi-delay data
                'att_lower',  # Only calculated for multi-delay data
                'att_upper',  # Only calculated for multi-delay data
                'abat_lower',  # Only calculated for multi-delay data
                'abat_upper',  # Only calculated for multi-delay data
                'abv_decayed',  # Only calculated for multi-delay (P)CASL data
                'plds',
                # SCORE/SCRUB outputs
                'cbf_ts_score',
                'mean_cbf_score',
                'mean_cbf_scrub',
                'score_outlier_index',
                # BASIL outputs
                'mean_cbf_basil',
                'mean_cbf_gm_basil',
                'mean_cbf_wm_basil',
                'att_basil',
                'abv_basil',
                # BASIL posterior standard deviations
                'cbf_sd_basil',
                'cbf_gm_sd_basil',
                'cbf_wm_sd_basil',
                'att_sd_basil',  # Only for multi-delay data
            ]
            + lr_fields
        ),
        name='outputnode',
    )

    # Warp tissue probability maps to ASL space
    def _pick_gm(files):
        if not isinstance(files, list):
            raise ValueError(f'_pick_gm: input is not list ({files})')
        return files[0]

    def _pick_wm(files):
        if not isinstance(files, list):
            raise ValueError(f'_pick_wm: input is not list ({files})')
        return files[1]

    def _pick_csf(files):
        if not isinstance(files, list):
            raise ValueError(f'_pick_csf: input is not list ({files})')
        return files[2]

    def _getfiledir(file):
        import os

        return os.path.dirname(file)

    gm_tfm = pe.Node(
        ApplyTransforms(
            interpolation='Gaussian',
            float=True,
            invert_transform_flags=[True],
            args='-v',
        ),
        name='gm_tfm',
        mem_gb=0.1,
    )

    workflow.connect([
        (inputnode, gm_tfm, [
            ('asl_mask', 'reference_image'),
            ('aslref2anat_xfm', 'transforms'),
            (('t1w_tpms', _pick_gm), 'input_image'),
        ]),
    ])  # fmt:skip

    wm_tfm = pe.Node(
        ApplyTransforms(
            interpolation='Gaussian',
            float=True,
            invert_transform_flags=[True],
            args='-v',
        ),
        name='wm_tfm',
        mem_gb=0.1,
    )

    workflow.connect([
        (inputnode, wm_tfm, [
            ('asl_mask', 'reference_image'),
            ('aslref2anat_xfm', 'transforms'),
            (('t1w_tpms', _pick_wm), 'input_image'),
        ]),
    ])  # fmt:skip

    csf_tfm = pe.Node(
        ApplyTransforms(
            interpolation='Gaussian',
            float=True,
            invert_transform_flags=[True],
            args='-v',
        ),
        name='csf_tfm',
        mem_gb=0.1,
    )

    workflow.connect([
        (inputnode, csf_tfm, [
            ('asl_mask', 'reference_image'),
            ('aslref2anat_xfm', 'transforms'),
            (('t1w_tpms', _pick_csf), 'input_image'),
        ]),
    ])  # fmt:skip

    extract_deltam = pe.Node(
        ExtractCBF(
            name_source=name_source,
            fwhm=smooth_kernel,
        ),
        mem_gb=0.2,
        run_without_submitting=True,
        name='extract_deltam',
    )

    workflow.connect([
        (inputnode, extract_deltam, [
            ('asl_file', 'asl_file'),
            ('metadata', 'metadata'),
            ('aslcontext', 'aslcontext'),
            ('m0scan_metadata', 'm0scan_metadata'),
        ]),
    ])  # fmt:skip

    if metadata['M0Type'] == 'Separate':
        mean_m0 = pe.Node(RobustAverage(), name='mean_m0', mem_gb=1)
        workflow.connect([
            (inputnode, mean_m0, [('m0scan', 'in_file')]),
            (mean_m0, extract_deltam, [('out_file', 'm0scan')]),
        ])  # fmt:skip

    # The multi-delay model fit is the only part of ComputeCBF that is slow enough to thread.
    # Its speedup levels off after a few threads, because it is limited by memory bandwidth.
    # On a GPU, chunks are fit one after another, so extra threads would not help.
    fit_on_gpu = config.nipype.use_gpu and is_multi_pld
    fit_nthreads = min(config.nipype.omp_nthreads or 1, 3) if is_multi_pld else 1
    fit_nthreads = 1 if fit_on_gpu else fit_nthreads
    compute_cbf = pe.Node(
        ComputeCBF(
            cbf_only=processing_target == 'cbf',
            m0_scale=m0_scale,
            num_threads=fit_nthreads,
            use_gpu=fit_on_gpu,
            write_calibration=pvc_method != 'none',
        ),
        mem_gb=1 + 0.1 * fit_nthreads if is_multi_pld else 0.2,
        n_procs=fit_nthreads,
        run_without_submitting=not is_multi_pld,
        name='compute_cbf',
    )

    workflow.connect([
        (extract_deltam, compute_cbf, [
            ('out_file', 'deltam'),
            ('m0_file', 'm0_file'),
            ('m0tr', 'm0tr'),
            ('metadata', 'metadata'),
        ]),
        (compute_cbf, outputnode, [
            ('cbf_ts', 'cbf_ts'),
            ('mean_cbf', 'mean_cbf'),
            ('att', 'att'),
            ('plds', 'plds'),
            ('abat', 'abat'),
            ('abv', 'abv'),
            ('att_lower', 'att_lower'),
            ('att_upper', 'att_upper'),
            ('abat_lower', 'abat_lower'),
            ('abat_upper', 'abat_upper'),
            ('abv_decayed', 'abv_decayed'),
        ]),
    ])  # fmt:skip

    if scorescrub:
        score_and_scrub_cbf = pe.Node(
            ScoreAndScrubCBF(tpm_threshold=0.7, wavelet_function='huber'),
            mem_gb=0.2,
            name='score_and_scrub_cbf',
            run_without_submitting=True,
        )

        workflow.__desc__ += """

Structural Correlation based Outlier Rejection (SCORE) algorithm was applied to the CBF timeseries
to discard CBF volumes with outlying values [@dolui2017structural] before computing the mean CBF.
Following SCORE, the Structural Correlation with RobUst Bayesian (SCRUB) algorithm was applied to
the CBF maps using structural tissue probability maps to reweight the mean CBF
[@dolui2017structural;@dolui2016scrub].
"""
        workflow.connect([
            (inputnode, score_and_scrub_cbf, [('asl_mask', 'mask')]),
            (compute_cbf, score_and_scrub_cbf, [('cbf_ts', 'cbf_ts')]),
            (gm_tfm, score_and_scrub_cbf, [('output_image', 'gm_tpm')]),
            (wm_tfm, score_and_scrub_cbf, [('output_image', 'wm_tpm')]),
            (csf_tfm, score_and_scrub_cbf, [('output_image', 'csf_tpm')]),
            (score_and_scrub_cbf, outputnode, [
                ('cbf_ts_score', 'cbf_ts_score'),
                ('score_outlier_index', 'score_outlier_index'),
                ('mean_cbf_score', 'mean_cbf_score'),
                ('mean_cbf_scrub', 'mean_cbf_scrub'),
            ]),
        ])  # fmt:skip

    if basil:
        # Partial volume correction takes most of oxford_asl's run time
        # (its spatial step runs ten times as many iterations as the uncorrected fit),
        # so it is skipped in sloppy mode.
        basil_pvc = not config.execution.sloppy
        if basil_pvc:
            pvc_str = ' and\nadditionally calculates a partial-volume corrected CBF image [@chappell_pvc].'
        else:
            pvc_str = (
                '.\nPartial volume correction was skipped, because *ASLPrep* was run in sloppy '
                'mode.\nThe partial-volume corrected CBF outputs are copies of the uncorrected '
                'BASIL CBF image.'
            )

        workflow.__desc__ += f"""

CBF was also computed with Bayesian Inference for Arterial Spin Labeling (BASIL)
[@chappell2008variational], as implemented in *FSL* {Info.version()}.
BASIL computes CBF using a spatial regularization of the estimated perfusion image{pvc_str}
"""
        if not basil_pvc:
            config.loggers.workflow.warning(
                'Sloppy mode: BASIL partial volume correction will be skipped. '
                'The partial-volume corrected CBF outputs (desc-basilGM and desc-basilWM) '
                'will be copies of the uncorrected BASIL CBF image.'
            )

        # Node to define bolus
        determine_bolus_duration = pe.Node(
            niu.Function(
                function=get_bolus_duration,
                input_names=['metadata', 'is_casl'],
                output_names=['bolus'],
            ),
            name='determine_bolus_duration',
        )
        determine_bolus_duration.inputs.is_casl = is_casl
        workflow.connect([(extract_deltam, determine_bolus_duration, [('metadata', 'metadata')])])

        # Node to define tis
        determine_inflow_times = pe.Node(
            niu.Function(
                function=get_inflow_times,
                input_names=['metadata', 'is_casl'],
                output_names=['tis'],
            ),
            name='determine_inflow_times',
        )
        determine_inflow_times.inputs.is_casl = is_casl

        workflow.connect([(extract_deltam, determine_inflow_times, [('metadata', 'metadata')])])

        # Node to estimate labeling efficiency
        estimate_alpha = pe.Node(
            niu.Function(
                function=estimate_labeling_efficiency,
                input_names=['metadata'],
                output_names=['labeling_efficiency'],
            ),
            name='estimate_alpha',
        )

        workflow.connect([(extract_deltam, estimate_alpha, [('metadata', 'metadata')])])

        basil_kwargs = prepare_basil_kwargs(metadata)

        basilcbf = pe.Node(
            BASILCBF(
                m0_scale=m0_scale,
                pvc=basil_pvc,
                pcasl=is_casl,
                t1blood=t1blood,
                **basil_kwargs,
            ),
            name='basilcbf',
            run_without_submitting=True,
            mem_gb=0.2,
        )

        workflow.connect([
            (inputnode, basilcbf, [('asl_mask', 'mask')]),
            (extract_deltam, basilcbf, [
                ('out_file', 'deltam'),
                ('m0_file', 'mzero'),
            ]),
            (determine_bolus_duration, basilcbf, [('bolus', 'bolus')]),
            (determine_inflow_times, basilcbf, [('tis', 'tis')]),
            (estimate_alpha, basilcbf, [('labeling_efficiency', 'alpha')]),
            (basilcbf, outputnode, [
                ('mean_cbf_basil', 'mean_cbf_basil'),
                ('cbf_sd_basil', 'cbf_sd_basil'),
            ]),
        ])  # fmt:skip

        if basil_pvc:
            workflow.connect([
                (gm_tfm, basilcbf, [('output_image', 'gm_tpm')]),
                (wm_tfm, basilcbf, [('output_image', 'wm_tpm')]),
                (basilcbf, outputnode, [
                    ('mean_cbf_gm_basil', 'mean_cbf_gm_basil'),
                    ('mean_cbf_wm_basil', 'mean_cbf_wm_basil'),
                    ('cbf_gm_sd_basil', 'cbf_gm_sd_basil'),
                    ('cbf_wm_sd_basil', 'cbf_wm_sd_basil'),
                ]),
            ])  # fmt:skip
        else:
            # Stand in for the corrected maps, so the rest of the workflow still runs.
            workflow.connect([
                (basilcbf, outputnode, [
                    ('mean_cbf_basil', 'mean_cbf_gm_basil'),
                    ('mean_cbf_basil', 'mean_cbf_wm_basil'),
                    ('cbf_sd_basil', 'cbf_gm_sd_basil'),
                    ('cbf_sd_basil', 'cbf_wm_sd_basil'),
                ]),
            ])  # fmt:skip

        if is_multi_pld:
            # BASIL's ATT and aBV estimates are only meaningful for multi-PLD data.
            # For single-delay data, ATT is fixed, and its "uncertainty" is only the prior.
            workflow.connect([
                (basilcbf, outputnode, [
                    ('att_basil', 'att_basil'),
                    ('abv_basil', 'abv_basil'),
                    ('att_sd_basil', 'att_sd_basil'),
                ]),
            ])  # fmt:skip

        if metadata['M0Type'] != 'Estimate':
            workflow.connect([(extract_deltam, basilcbf, [('m0tr', 'm0tr')])])

    if lr_pvc:
        asl_pv_wf = init_asl_pv_wf(
            metadata=metadata,
            effective_resolution=asl_effective_resolution,
        )
        lr_pvc_wf = init_lr_pvc_wf(
            metadata=metadata,
            sources=lr_sources,
            kernel_size=tuple(pvc_kernel_size),
            multi_delay=is_multi_pld,
            num_threads=fit_nthreads,
            use_gpu=fit_on_gpu,
            mem_gb=_lr_pvc_mem_gb(name_source, fit_nthreads) if is_multi_pld else None,
        )
        if is_multi_pld:
            # Each delay volume is unmixed, so the fit needs the delta-M series and the
            # metadata of the selected volumes, plus the uncorrected fit for the CSF QC template.
            workflow.connect([
                (extract_deltam, lr_pvc_wf, [
                    ('out_file', 'inputnode.deltam'),
                    ('metadata', 'inputnode.asl_metadata'),
                ]),
                (compute_cbf, lr_pvc_wf, [
                    ('att', 'inputnode.att'),
                    ('abat', 'inputnode.abat'),
                    ('abv', 'inputnode.abv'),
                ]),
            ])  # fmt:skip
        workflow.connect([
            (inputnode, asl_pv_wf, [
                ('t1w_tpms', 'inputnode.t1w_tpms'),
                ('asl_mask', 'inputnode.asl_mask'),
                ('aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
            ]),
            (inputnode, lr_pvc_wf, [('asl_mask', 'inputnode.asl_mask')]),
            (compute_cbf, lr_pvc_wf, [
                ('mean_cbf', 'inputnode.mean_cbf'),
                ('m0_calib', 'inputnode.m0_calib'),
                ('quant_valid', 'inputnode.quant_valid'),
            ]),
            (asl_pv_wf, lr_pvc_wf, [
                ('outputnode.gm_pv', 'inputnode.gm_pv'),
                ('outputnode.wm_pv', 'inputnode.wm_pv'),
                ('outputnode.csf_pv', 'inputnode.csf_pv'),
            ]),
            (asl_pv_wf, outputnode, [
                ('outputnode.gm_pv', 'gm_pv'),
                ('outputnode.wm_pv', 'wm_pv'),
                ('outputnode.csf_pv', 'csf_pv'),
                ('outputnode.pv_metadata', 'pv_metadata'),
            ]),
            (lr_pvc_wf, outputnode, [
                (f'outputnode.{field}', field)
                for field in lr_fields
                if field not in ('gm_pv', 'wm_pv', 'csf_pv', 'pv_metadata')
            ]),
        ])  # fmt:skip
        if scorescrub:
            workflow.connect([
                (score_and_scrub_cbf, lr_pvc_wf, [('mean_cbf_score', 'inputnode.mean_cbf_score')]),
            ])  # fmt:skip

    return workflow


def init_parcellate_cbf_wf(
    cbf_3d,
    sparse_3d=None,
    min_coverage=0.5,
    mem_gb=0.1,
    omp_nthreads=1,
    name='parcellate_cbf_wf',
):
    """Parcellate CBF results using a set of atlases.

    Atlases are in MNI152NLin6Asym space, but ASLPrep is not guaranteed to generate
    MNI152NLin6Asym-space derivatives.
    However, it *is* guaranteed to generate MNI152NLin2009cAsym derivatives,
    so this workflow warps the atlases to asl reference space via the path:
    MNI152NLin6Asym --> MNI152NLin2009cAsym --> anat --> asl.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.cbf import init_parcellate_cbf_wf

            wf = init_parcellate_cbf_wf(cbf_3d=["mean_cbf", "mean_cbf_score"])

    Parameters
    ----------
    cbf_3d : list of str
        CBF maps to parcellate within the ASL brain mask.
    sparse_3d : dict or None
        CBF maps that are only estimated in some voxels, keyed by name, with the name of their
        validity mask as the value. They are parcellated within that mask, so voxels that were
        not estimated (stored as zero) do not enter parcel means, and their coverage files
        report estimation coverage.
    min_coverage : float
        Parcels with less coverage than this are set to zero.
    mem_gb : float
        Memory for the atlas warping and parcellation nodes.
    omp_nthreads : int
        Threads for the atlas warping.
    name : str
        Default is "parcellate_cbf_wf".

    Inputs
    ------
    source_file : str
        The ASL file the parcellated tables are named after.
    asl_mask : str
        ASL brain mask; the grid the atlases are warped to, and the mask for ``cbf_3d`` maps.
    aslref2anat_xfm : str
        Transform from the ASL reference to the anatomical image.
    mni2009c2anat_xfm : str
        The transform from MNI152NLin2009cAsym to the subject's anatomical space.
    atlas_names : list of str
        Names of the atlases, from :func:`~aslprep.workflows.asl.cbf.init_load_atlases_wf`.
    atlas_files : list of str
        The atlases, in MNI152NLin6Asym space.
    atlas_labels_files : list of str
        The atlases' label tables.
    <each name in cbf_3d> : str
        The CBF maps to parcellate within ``asl_mask``.
    <each key and value of sparse_3d> : str
        The sparsely estimated CBF maps, and the validity masks they are parcellated within.

    Outputs
    -------
    None. For each map and atlas, a parcellated CBF table (``*_cbf.tsv``) is written to the
    output directory. ``mean_cbf``, ``mean_cbf_basil`` and the ``sparse_3d`` maps also get a
    coverage table (``*_coverage.tsv``).
    """
    CBF_ENTITIES = {
        'mean_cbf': {},
        'mean_cbf_score': {
            'desc': 'score',
        },
        'mean_cbf_scrub': {
            'desc': 'scrub',
        },
        'mean_cbf_basil': {
            'desc': 'basil',
        },
        'mean_cbf_gm_basil': {
            'desc': 'basilGM',
        },
        'mean_cbf_wm_basil': {
            'desc': 'basilWM',
        },
    }
    sparse_3d = sparse_3d or {}
    for sparse_name in sparse_3d:
        CBF_ENTITIES[sparse_name] = {'desc': BASE_INPUT_FIELDS[sparse_name]['desc']}
    workflow = Workflow(name=name)
    output_dir = config.execution.output_dir

    workflow.__desc__ = f"""\
CBF maps were parcellated in ASL reference space using these atlases.
In cases of partial coverage, either uncovered voxels (values of all zeros or NaNs) were
ignored (when the parcel had >{min_coverage * 100}% coverage)
or the whole parcel was set to zero (when the parcel had <{min_coverage * 100}% coverage).
"""

    input_fields = [
        'source_file',
        'asl_mask',
        'aslref2anat_xfm',
        'mni2009c2anat_xfm',
        'atlas_names',
        'atlas_files',
        'atlas_labels_files',
    ]
    input_fields += cbf_3d
    input_fields += list(sparse_3d) + sorted(set(sparse_3d.values()))
    inputnode = pe.Node(
        niu.IdentityInterface(fields=input_fields),
        name='inputnode',
    )

    # Prepare to warp atlases from MNI152NLin6Asym to ASL reference space.
    # One of the output spaces selected by the user *may* be MNI152NLin6Asym,
    # but MNI152NLin2009cAsym is always used, so it's safer to go:
    # MNI152NLin6Asym --> MNI152NLin2009cAsym --> anat --> asl
    # Atlases are in MNI152NLin6Asym
    MNI152NLin6Asym_to_MNI152NLin2009cAsym = str(
        get_template(
            template='MNI152NLin2009cAsym',
            mode='image',
            suffix='xfm',
            extension='.h5',
            **{'from': 'MNI152NLin6Asym'},
        ),
    )

    merge_xforms = pe.Node(niu.Merge(3), name='merge_xforms')
    merge_xforms.inputs.in1 = MNI152NLin6Asym_to_MNI152NLin2009cAsym
    workflow.connect([
        (inputnode, merge_xforms, [
            ('mni2009c2anat_xfm', 'in2'),
            ('aslref2anat_xfm', 'in3'),
        ]),
    ])  # fmt:skip

    # Using the generated transforms, apply them to get everything in the correct MNI form
    warp_atlases_to_asl_space = pe.MapNode(
        ApplyTransforms(
            interpolation='GenericLabel',
            input_image_type=3,
            dimension=3,
            invert_transform_flags=[False, False, True],
            args='-v',
        ),
        name='warp_atlases_to_asl_space',
        iterfield=['input_image'],
        mem_gb=mem_gb,
        n_procs=omp_nthreads,
    )
    workflow.connect([
        (inputnode, warp_atlases_to_asl_space, [
            ('asl_mask', 'reference_image'),
            ('atlas_files', 'input_image'),
        ]),
        (merge_xforms, warp_atlases_to_asl_space, [('out', 'transforms')]),
    ])  # fmt:skip

    atlas_srcs = pe.MapNode(
        BIDSURI(
            numinputs=5,
            dataset_links=config.execution.dataset_links,
            out_dir=str(output_dir),
        ),
        name='atlas_srcs',
        iterfield=['in1'],
    )
    atlas_srcs.inputs.in5 = MNI152NLin6Asym_to_MNI152NLin2009cAsym
    workflow.connect([
        (inputnode, atlas_srcs, [
            ('atlas_files', 'in1'),
            ('asl_mask', 'in2'),
            ('aslref2anat_xfm', 'in3'),
            ('mni2009c2anat_xfm', 'in4'),
        ]),
    ])  # fmt:skip

    for cbf_type in cbf_3d + list(sparse_3d):
        parcellate_cbf = pe.MapNode(
            NiftiParcellate(min_coverage=min_coverage),
            name=f'parcellate_{cbf_type}',
            iterfield=['atlas', 'atlas_labels'],
            mem_gb=mem_gb,
        )

        workflow.connect([
            (inputnode, parcellate_cbf, [
                (cbf_type, 'in_file'),
                (sparse_3d.get(cbf_type, 'asl_mask'), 'mask'),
                ('atlas_labels_files', 'atlas_labels'),
            ]),
            (warp_atlases_to_asl_space, parcellate_cbf, [('output_image', 'atlas')]),
        ])  # fmt:skip

        ds_cbf = pe.MapNode(
            DerivativesDataSink(
                base_directory=config.execution.aslprep_dir,
                check_hdr=False,
                suffix='cbf',
                **CBF_ENTITIES[cbf_type],
            ),
            name=f'ds_{cbf_type}',
            iterfield=['atlas', 'in_file', 'meta_dict'],
            run_without_submitting=True,
        )
        workflow.connect([
            (inputnode, ds_cbf, [
                ('source_file', 'source_file'),
                ('atlas_names', 'atlas'),
            ]),
            (parcellate_cbf, ds_cbf, [('timeseries', 'in_file')]),
            (atlas_srcs, ds_cbf, [('metadata', 'meta_dict')]),
        ])  # fmt:skip

        # Maps with their own validity masks also get coverage files: estimation coverage.
        if cbf_type in ('mean_cbf', 'mean_cbf_basil') or cbf_type in sparse_3d:
            # I think it is easier to only retain the coverage file for the regular CBF estimates.
            # SCORE/SCRUB CBF should have the same coverage as the regular CBF.
            # BASIL might have different coverage because it drops any voxels with negative CBF,
            # but the different PVC types will have the same coverage as the main BASIL CBF.
            ds_coverage = pe.MapNode(
                DerivativesDataSink(
                    base_directory=config.execution.aslprep_dir,
                    check_hdr=False,
                    suffix='coverage',
                    **CBF_ENTITIES[cbf_type],
                ),
                name=f'ds_coverage_{cbf_type}',
                iterfield=['atlas', 'in_file', 'meta_dict'],
                run_without_submitting=True,
            )
            workflow.connect([
                (inputnode, ds_coverage, [
                    ('source_file', 'source_file'),
                    ('atlas_names', 'atlas'),
                ]),
                (parcellate_cbf, ds_coverage, [('coverage', 'in_file')]),
                (atlas_srcs, ds_coverage, [('metadata', 'meta_dict')]),
            ])  # fmt:skip

    return workflow


def init_load_atlases_wf(
    name='load_atlases_wf',
):
    """Collect atlases and copy them to the sourcedata/atlases directory.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep import config
            from aslprep.data import load as load_data
            from aslprep.workflows.cbf import init_load_atlases_wf

            with mock_config():
                config.execution.derivatives = {
                    "xcpdatlases": str(load_data("atlases")),
                }
                wf = init_load_atlases_wf()

    Parameters
    ----------
    %(name)s
        Default is "load_atlases_wf".

    Outputs
    -------
    atlas_names
    atlas_files
    atlas_labels_files
    """
    import os

    from aslprep.interfaces.bids import CopyAtlas, CopyAtlasDescription
    from aslprep.utils.atlas import collect_atlases
    from aslprep.utils.boilerplate import describe_atlases

    workflow = Workflow(name=name)
    output_dir = config.execution.output_dir
    atlas_output_dir = os.path.join(output_dir, 'sourcedata', 'atlases')
    os.makedirs(atlas_output_dir, exist_ok=True)

    atlases = collect_atlases(
        datasets=config.execution.derivatives,
        atlases=config.execution.atlases,
        bids_filters=config.execution.bids_filters,
    )

    # Reorganize the atlas file information
    atlas_names, atlas_files, atlas_labels_files, atlas_metadata = [], [], [], []
    atlas_datasets, atlas_dataset_paths = [], []
    for atlas, atlas_dict in atlases.items():
        config.loggers.workflow.info(f'Loading atlas: {atlas}')
        atlas_names.append(atlas)
        atlas_datasets.append(atlas_dict['dataset'])
        atlas_dataset_paths.append(atlas_dict['dataset_path'])
        atlas_files.append(atlas_dict['image'])
        atlas_labels_files.append(atlas_dict['labels'])
        atlas_metadata.append(atlas_dict['metadata'])

    # Write a description
    atlas_str = describe_atlases(atlas_names)
    workflow.__desc__ = f"""\
The following atlases were used in the workflow: {atlas_str}.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                # Set within the workflow
                'atlas_names',
                'atlas_datasets',
                'atlas_dataset_paths',
                'atlas_files',
                'atlas_labels_files',
                'atlas_metadata',
            ],
        ),
        name='inputnode',
    )
    inputnode.inputs.atlas_names = atlas_names
    inputnode.inputs.atlas_datasets = atlas_datasets
    inputnode.inputs.atlas_dataset_paths = atlas_dataset_paths
    inputnode.inputs.atlas_files = atlas_files
    inputnode.inputs.atlas_labels_files = atlas_labels_files
    inputnode.inputs.atlas_metadata = atlas_metadata

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'atlas_names',
                'atlas_files',
                'atlas_labels_files',
            ],
        ),
        name='outputnode',
    )
    workflow.connect([(inputnode, outputnode, [('atlas_names', 'atlas_names')])])

    copy_atlas_description = pe.MapNode(
        CopyAtlasDescription(output_dir=atlas_output_dir),
        name='copy_atlas_description',
        iterfield=['in_dir', 'atlas_name'],
        run_without_submitting=True,
    )
    workflow.connect([
        (inputnode, copy_atlas_description, [
            ('atlas_dataset_paths', 'in_dir'),
            ('atlas_names', 'atlas_name'),
        ]),
    ])  # fmt:skip

    atlas_src = pe.MapNode(
        BIDSURI(
            numinputs=1,
            dataset_links=config.execution.dataset_links,
            out_dir=str(output_dir),
        ),
        name='atlas_src',
        iterfield=['in1', 'metadata'],
    )
    workflow.connect([
        (inputnode, atlas_src, [
            ('atlas_files', 'in1'),
            ('atlas_metadata', 'metadata'),
        ]),
    ])  # fmt:skip

    copy_atlas_file = pe.MapNode(
        CopyAtlas(output_dir=atlas_output_dir),
        name='copy_atlas_file',
        iterfield=['name_source', 'in_file', 'atlas', 'meta_dict'],
        run_without_submitting=True,
    )
    workflow.connect([
        (inputnode, copy_atlas_file, [
            ('atlas_names', 'atlas'),
            ('atlas_files', 'name_source'),
            ('atlas_files', 'in_file'),
        ]),
        (atlas_src, copy_atlas_file, [('metadata', 'meta_dict')]),
        (copy_atlas_file, outputnode, [('out_file', 'atlas_files')]),
    ])  # fmt:skip

    copy_atlas_labels_file = pe.MapNode(
        CopyAtlas(output_dir=atlas_output_dir),
        name='copy_atlas_labels_file',
        iterfield=['name_source', 'in_file', 'atlas'],
        run_without_submitting=True,
    )
    workflow.connect([
        (inputnode, copy_atlas_labels_file, [
            ('atlas_names', 'atlas'),
            ('atlas_labels_files', 'name_source'),
            ('atlas_labels_files', 'in_file'),
        ]),
        (copy_atlas_labels_file, outputnode, [('out_file', 'atlas_labels_files')]),
    ])  # fmt:skip

    return workflow
