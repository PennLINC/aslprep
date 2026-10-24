# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for calculating CBF."""

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
    BASILInflowTimeImage,
    ComputeCBF,
    ExtractCBF,
    ScoreAndScrubCBF,
)
from aslprep.interfaces.parcellation import NiftiParcellate
from aslprep.utils.asl import (
    BOLUS_CUTOFF_TECHNIQUES,
    describe_labeling_efficiency,
    determine_multi_pld,
    get_bolus_duration,
    get_inflow_times,
    infer_m0tr,
    pcasl_or_pasl,
    prepare_basil_kwargs,
)
from aslprep.utils.calibration import (
    DEFAULT_M0_SCALE,
    DEFAULT_SMOOTH_KERNEL,
    basil_calibration_inputs,
    calibration_value,
    resolve_calibration,
)
from aslprep.utils.cbf import REFERENCE_HAEMATOCRIT, T1BLOOD_BY_FIELD_STRENGTH
from aslprep.utils.misc import estimate_asl_mem_usage
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


def basil_pvc_enabled(*, basil, processing_target):
    """Whether BASIL runs partial volume correction.

    BASIL does not run on pre-calculated CBF. Partial volume correction takes most of
    oxford_asl's run time (its spatial step runs ten times as many iterations as the uncorrected
    fit), so it is skipped in sloppy mode.
    """
    return basil and processing_target != 'cbf' and not config.execution.sloppy


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
    calibration=None,
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
        Used only to resolve ``calibration`` when it is None.
    smooth_kernel
    pvc_method : {'none', 'local-regression'}
        Additional partial volume correction of the mean CBF (and SCORE) maps.
    pvc_kernel_size : :obj:`tuple` of :obj:`int`
        Local-regression kernel size, ordered (in-plane, in-plane, slice).
    asl_effective_resolution : :obj:`str` or :obj:`list` of :obj:`float`
        Effective resolution used to build the partial volume maps,
        which both BASIL and local-regression PVC use.
    calibration : :obj:`dict` or None
        The calibration specification that every quantification step uses
        (see :func:`~aslprep.utils.calibration.resolve_calibration`).
        If None, it is resolved from ``metadata``, ``m0_scale`` and ``smooth_kernel``.
        Not used for pre-calculated CBF data.
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
        (with ``pvc_method='local-regression'``, or with ``basil`` when BASIL runs partial
        volume correction; see :func:`~aslprep.workflows.asl.pv.init_asl_pv_wf`).
    mean_cbf_gm_lr, ..., csf_noise_lr, pvc_summary, pvc_metadata, pvc_qc
        Local-regression partial volume correction outputs for ``mean_cbf`` (with
        ``pvc_method='local-regression'``; see :func:`~aslprep.workflows.asl.pv.init_lr_pvc_wf`
        for each field).
    mean_cbf_gm_lr_score, ..., csf_noise_lr_score, pvc_summary_score
        The same for ``mean_cbf_score`` (with ``pvc_method='local-regression'`` and
        ``scorescrub``; single-delay data only).
    att_gm_lr, ..., csf_kinres_lr
        For multi-delay data, the local-regression tissue ATT maps, their bounds, the
        kinetic-fit status and residuals, and the CSF curve's residual; the standard errors
        are absent (see :func:`~aslprep.workflows.asl.pv.init_lr_pvc_wf`).
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

    if processing_target != 'cbf' and calibration is None:
        calibration = resolve_calibration(
            metadata,
            m0tr=m0tr,
            m0scan_metadata=m0scan_metadata,
            m0_scale=m0_scale,
            m0_scale_source='user' if m0_scale != DEFAULT_M0_SCALE else 'default',
            smooth_kernel=smooth_kernel,
            smooth_kernel_source=('user' if smooth_kernel != DEFAULT_SMOOTH_KERNEL else 'default'),
            basil=basil,
        )
        for warning in calibration['Warnings']:
            config.loggers.workflow.warning(warning)

    m0_str = ''
    if calibration is not None:
        m0_scale = calibration_value(calibration, 'M0Scale')
        # The M0 image is smoothed as the specification records (absent for an M0Estimate).
        smooth_kernel = calibration_value(calibration, 'M0SmoothingFWHM', smooth_kernel)
        recovery_tr = calibration_value(calibration, 'M0RepetitionTime')
        t1tissue = calibration_value(calibration, 'TissueT1')
        m0_threshold = calibration_value(calibration, 'M0Threshold')

    if processing_target == 'cbf':
        # Pre-calculated CBF is not calibrated.
        pass

    elif m0type in ('Included', 'Separate'):
        m0_str = (
            'Calibration (M0) volumes associated with the ASL scan were smoothed with a '
            f'Gaussian kernel (FWHM={smooth_kernel} mm) and the average calibration image was '
            f'calculated. '
        )
        if t1tissue is not None:
            m0_str += (
                f'As the TR of the calibration volumes was less than 5 seconds ({recovery_tr}), '
                'the average calibration volume was divided by 1 - exp(-TR / T1) to '
                'compensate for incomplete T1 recovery, as recommended in '
                f'@alsop_recommended_2015, using a gray matter T1 of {t1tissue:g} s '
                '[@wright2008water]. '
            )

        m0_str += f'The calibration volume was then scaled by {m0_scale}.'

    elif m0type == 'Estimate':
        m0_str = (
            f'A single M0 estimate of {metadata["M0Estimate"]}, scaled by {m0_scale}, was used '
            'as a uniform calibration image. '
            'As BIDS defines this estimate as the M0 of arterial blood, it was not divided by '
            'the blood-brain partition coefficient.'
        )
    else:
        m0_str = (
            f'As no calibration images or provided M0 estimate was available for the ASL scan, '
            'the control volumes were used as a substitute. '
            'The control volumes in the ASL scans were smoothed with a '
            f'Gaussian kernel (FWHM={smooth_kernel} mm) and the average control image was '
            f'calculated'
        )
        if t1tissue is not None:
            m0_str += (
                f'. As the TR of the control volumes was less than 5 seconds ({recovery_tr}), '
                'the average control image was divided by 1 - exp(-TR / T1) to '
                'compensate for incomplete T1 recovery, as recommended in '
                f'@alsop_recommended_2015, using a gray matter T1 of {t1tissue:g} s '
                '[@wright2008water], and then scaled'
            )
        else:
            m0_str += ' and scaled'
        m0_str += f' by {m0_scale}.'

    if m0_str and m0_threshold is not None:
        m0_str += (
            f' Voxels where it fell below {100 * m0_threshold:g}% of its 98th percentile were '
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

        m0b_str = (
            'the M0 estimate'
            if m0type == 'Estimate'
            else 'M0 divided by the blood-brain partition coefficient'
        )
        workflow.__desc__ += f"""\
: *ASLPrep* calculated cerebral blood flow (CBF) from the multi-delay
{metadata['ArterialSpinLabelingType']} data using a two-compartment general kinetic model (GKM)
[@buxton1998general], as recommended and extended in @woods2023recommendations.
This model contains separate terms for tissue and macrovascular delta-M signals.

{m0_str}
The equilibrium magnetization of arterial blood (M0b) was estimated voxel-wise as
{m0b_str},
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

    if calibration is not None:
        # Every model (including BASIL) uses these values, from the calibration specification.
        t1blood = calibration_value(calibration, 'BloodT1')
        field_strength = calibration_value(calibration, 'MagneticFieldStrength')
        if calibration['Parameters']['BloodT1']['Source'] == 'inferred':
            workflow.__desc__ += (
                f' The T1 of arterial blood was set to {t1blood:.3f} s: '
                f'the value for {field_strength:g} T, scaled for the '
                f"participant's measured haematocrit "
                f'({calibration_value(calibration, "Haematocrit"):g}) relative to '
                f'a typical haematocrit of {REFERENCE_HAEMATOCRIT:g} [@lu2004determining], '
                'using the model of @hales2016general.'
            )
        elif field_strength not in T1BLOOD_BY_FIELD_STRENGTH:
            workflow.__desc__ += (
                f' The T1 of arterial blood was set to {t1blood:g} s, '
                f'calculated from the field strength of {field_strength:g} T [@zhang2013vivo].'
            )
        else:
            t1blood_citation = (
                '@alsop_recommended_2015' if field_strength in (1.5, 3) else '@zhang2013vivo'
            )
            workflow.__desc__ += (
                f' The T1 of arterial blood was set to {t1blood:g} s, '
                f'the value for {field_strength:g} T [{t1blood_citation}].'
            )

        if m0type == 'Estimate':
            workflow.__desc__ += (
                f' CBF was calculated with {describe_labeling_efficiency(metadata)}.'
            )
        else:
            workflow.__desc__ += (
                ' CBF was calculated with a blood-brain partition coefficient of '
                f'{calibration_value(calibration, "PartitionCoefficient"):g} mL/g '
                f'[@alsop_recommended_2015] and {describe_labeling_efficiency(metadata)}.'
            )

    # End the paragraph, so that the next workflow's description does not run into it.
    workflow.__desc__ += '\n'

    lr_pvc = pvc_method == 'local-regression'
    basil_pvc = basil_pvc_enabled(basil=basil, processing_target=processing_target)
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
        lr_fields = ['pvc_metadata', 'pvc_qc']
        lr_fields += lr_pvc_fields(lr_sources, multi_delay=is_multi_pld)
    # Both partial volume corrections use the same voxel-averaged maps (gh-740)
    pv_fields = ['gm_pv', 'wm_pv', 'csf_pv', 'pv_metadata'] if (lr_pvc or basil_pvc) else []

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
            + pv_fields
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

    if scorescrub:
        # SCORE/SCRUB samples the tissue probability maps at each ASL voxel centre
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

    # Nodes that load the ASL series, or a CBF time series of the same size, run in workers,
    # so the scheduler needs their memory to pack them.
    series_mem_gb = estimate_asl_mem_usage(name_source)[1]['resampled']
    extract_deltam = pe.Node(
        ExtractCBF(
            name_source=name_source,
            fwhm=smooth_kernel,
        ),
        mem_gb=series_mem_gb,
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
            num_threads=fit_nthreads,
            use_gpu=fit_on_gpu,
            write_calibration=pvc_method != 'none',
        ),
        mem_gb=1 + 0.1 * fit_nthreads if is_multi_pld else series_mem_gb,
        n_procs=fit_nthreads,
        name='compute_cbf',
    )
    if calibration is not None:
        compute_cbf.inputs.calibration = calibration

    workflow.connect([
        (extract_deltam, compute_cbf, [
            ('out_file', 'deltam'),
            ('m0_file', 'm0_file'),
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

    if pv_fields:
        asl_pv_wf = init_asl_pv_wf(
            metadata=metadata,
            effective_resolution=asl_effective_resolution,
        )
        workflow.connect([
            (inputnode, asl_pv_wf, [
                ('t1w_tpms', 'inputnode.t1w_tpms'),
                ('asl_mask', 'inputnode.asl_mask'),
                ('aslref2anat_xfm', 'inputnode.aslref2anat_xfm'),
            ]),
            (asl_pv_wf, outputnode, [
                ('outputnode.gm_pv', 'gm_pv'),
                ('outputnode.wm_pv', 'wm_pv'),
                ('outputnode.csf_pv', 'csf_pv'),
                ('outputnode.pv_metadata', 'pv_metadata'),
            ]),
        ])  # fmt:skip

    if scorescrub:
        score_and_scrub_cbf = pe.Node(
            ScoreAndScrubCBF(tpm_threshold=0.7, wavelet_function='huber'),
            mem_gb=series_mem_gb,
            name='score_and_scrub_cbf',
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
        if basil_pvc:
            pvc_str = (
                ' and\nadditionally calculates a partial-volume corrected CBF image [@chappell_pvc].'
                '\nBASIL was given the gray and white matter partial volume maps described below.'
            )
        else:
            pvc_str = (
                '.\nPartial volume correction was skipped, because *ASLPrep* was run in sloppy '
                'mode.\nThe partial-volume corrected CBF outputs are copies of the uncorrected '
                'BASIL CBF image.'
            )

        # ASLPrep's own models treat the bolus duration as known, so BASIL fixes it, too.
        # Without --fixbolus, oxford_asl estimates it for multi-delay PASL data (gh-759).
        bolus_str = 'labeling duration' if is_casl else 'bolus cut-off delay time'
        workflow.__desc__ += f"""

CBF was also computed with Bayesian Inference for Arterial Spin Labeling (BASIL)
[@chappell2008variational], as implemented in *FSL* {Info.version()}.
BASIL computes CBF using a spatial regularization of the estimated perfusion image{pvc_str}
As in *ASLPrep*'s own model, the bolus duration was fixed to the {bolus_str},
rather than estimated.
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

        basil_kwargs, use_tiimg = prepare_basil_kwargs(metadata)
        if use_tiimg:
            workflow.__desc__ += (
                'BASIL was given the inflow time of each voxel, shifted by the acquisition '
                'time of its slice.\n'
            )
        elif 'slice_spacing' in basil_kwargs:
            band_str = (
                f' within each multiband band of {basil_kwargs["sliceband"]} slices'
                if 'sliceband' in basil_kwargs
                else ''
            )
            workflow.__desc__ += (
                'BASIL accounted for slice timing by adding '
                f'{1000 * basil_kwargs["slice_spacing"]:.2f} ms to the inflow time for each '
                f'successive slice{band_str}.\n'
            )

        basilcbf = pe.Node(
            BASILCBF(
                pvc=basil_pvc,
                pcasl=is_casl,
                fixbolus=True,
                # The calibration and physiological parameters of ASLPrep's own models
                **basil_calibration_inputs(calibration),
                **basil_kwargs,
            ),
            name='basilcbf',
            # fabber peaked at 1.2 times the series estimate for multi-delay data
            mem_gb=1.5 * series_mem_gb,
        )

        workflow.connect([
            (inputnode, basilcbf, [('asl_mask', 'mask')]),
            (extract_deltam, basilcbf, [('out_file', 'deltam')]),
            (determine_bolus_duration, basilcbf, [('bolus', 'bolus')]),
            (basilcbf, outputnode, [
                ('mean_cbf_basil', 'mean_cbf_basil'),
                ('cbf_sd_basil', 'cbf_sd_basil'),
            ]),
        ])  # fmt:skip
        if m0type == 'Estimate':
            # BASIL gets the M0 estimate as a single value of blood M0 (m0_value)
            workflow.__desc__ += (
                'BASIL was given the scaled M0 estimate as the M0 of arterial blood.\n'
            )
        else:
            workflow.connect([(extract_deltam, basilcbf, [('m0_file', 'mzero')])])

        if use_tiimg:
            make_tiimg = pe.Node(
                BASILInflowTimeImage(),
                name='make_tiimg',
                mem_gb=series_mem_gb,
            )
            workflow.connect([
                (extract_deltam, make_tiimg, [
                    ('out_file', 'deltam'),
                    ('metadata', 'metadata'),
                ]),
                (determine_inflow_times, make_tiimg, [('tis', 'tis')]),
                (make_tiimg, basilcbf, [('out_file', 'tiimg')]),
            ])  # fmt:skip
        else:
            workflow.connect([(determine_inflow_times, basilcbf, [('tis', 'tis')])])

        if basil_pvc:
            # Voxel-averaged maps, not the point-sampled TPMs (gh-740)
            workflow.connect([
                (asl_pv_wf, basilcbf, [
                    ('outputnode.gm_pv', 'gm_tpm'),
                    ('outputnode.wm_pv', 'wm_tpm'),
                ]),
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

    if lr_pvc:
        lr_pvc_wf = init_lr_pvc_wf(
            metadata=metadata,
            sources=lr_sources,
            kernel_size=tuple(pvc_kernel_size),
            multi_delay=is_multi_pld,
            num_threads=fit_nthreads,
            use_gpu=fit_on_gpu,
            mem_gb=_lr_pvc_mem_gb(name_source, fit_nthreads) if is_multi_pld else None,
            calibration=calibration,
        )
        if is_multi_pld:
            # Each delay volume is unmixed, so the fit needs the delta-M series and the
            # metadata of the selected volumes, plus the uncorrected ATT for the CSF QC basis.
            workflow.connect([
                (extract_deltam, lr_pvc_wf, [
                    ('out_file', 'inputnode.deltam'),
                    ('metadata', 'inputnode.asl_metadata'),
                ]),
                (compute_cbf, lr_pvc_wf, [('att', 'inputnode.att')]),
            ])  # fmt:skip
        workflow.connect([
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
            (lr_pvc_wf, outputnode, [(f'outputnode.{field}', field) for field in lr_fields]),
        ])  # fmt:skip
        if scorescrub:
            workflow.connect([
                (score_and_scrub_cbf, lr_pvc_wf, [('mean_cbf_score', 'inputnode.mean_cbf_score')]),
            ])  # fmt:skip

    return workflow


def init_parcellate_cbf_wf(
    maps_3d,
    sparse_3d=None,
    written_maps=(),
    min_coverage=0.5,
    mem_gb=0.1,
    omp_nthreads=1,
    name='parcellate_cbf_wf',
):
    """Parcellate CBF, ATT and other perfusion maps using a set of atlases.

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

            wf = init_parcellate_cbf_wf(maps_3d=["mean_cbf", "mean_cbf_score"])

    Parameters
    ----------
    maps_3d : list of str
        Maps to parcellate within the ASL brain mask.
    sparse_3d : dict or None
        Maps that are only estimated in some voxels, keyed by name, with the name of their
        validity mask as the value. They are parcellated within that mask, so voxels that were
        not estimated (stored as zero) do not enter parcel means, and their coverage columns
        report estimation coverage.
    written_maps : sequence of str
        The maps in ``maps_3d`` or ``sparse_3d`` that are written to the output directory in
        ASL reference space. Their written paths are listed, as BIDS-URIs, in the ``Sources``
        of their ``scalarstats`` columns.
    min_coverage : float
        Parcels with less coverage than this are set to n/a.
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
        ASL brain mask; the grid the atlases are warped to, and the mask for ``maps_3d``.
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
    <each name in maps_3d> : str
        The maps to parcellate within ``asl_mask``.
    <each key and value of sparse_3d> : str
        The sparsely estimated maps, and the validity masks they are parcellated within.
    <each name in written_maps>_file : str
        The written path of that map.

    Outputs
    -------
    None. For each atlas, a table of parcel means (``*_scalarstats.tsv``) and a coverage table
    (``*_coverage.tsv``) are written to the output directory. Each has a ``node`` column,
    with one row per parcel, then one column per map, named ``<suffix>[_<desc>]`` after the
    map's voxel-wise output (for example, ``cbf``, ``cbf_scrub`` or ``att_basil``).
    """
    sparse_3d = sparse_3d or {}
    workflow = Workflow(name=name)
    output_dir = config.execution.output_dir

    workflow.__desc__ = f"""\
CBF and other perfusion maps were parcellated in ASL reference space using these atlases.
In cases of partial coverage, either uncovered voxels were
ignored (when the parcel had >{min_coverage * 100}% coverage)
or the whole parcel was set to n/a (when the parcel had <{min_coverage * 100}% coverage).
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
    input_fields += maps_3d
    input_fields += list(sparse_3d) + sorted(set(sparse_3d.values()))
    map_names = maps_3d + list(sparse_3d)
    if unknown := sorted(set(written_maps) - set(map_names)):
        raise ValueError(f'Written maps {unknown} are not parcellated.')
    input_fields += [f'{map_name}_file' for map_name in map_names if map_name in written_maps]
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

    # Every map is parcellated with each atlas in one node, which writes one table of parcel
    # means (scalarstats) and one coverage table, with a row per parcel and a column per map.
    columns = []
    for map_name in map_names:
        fields = BASE_INPUT_FIELDS[map_name]
        columns.append(
            '_'.join([fields['suffix']] + ([fields['desc']] if 'desc' in fields else []))
        )

    merge_maps = pe.Node(
        niu.Merge(len(map_names)),
        name='merge_maps',
        run_without_submitting=True,
    )
    merge_masks = pe.Node(
        niu.Merge(len(map_names)),
        name='merge_masks',
        run_without_submitting=True,
    )
    for i_map, map_name in enumerate(map_names, start=1):
        workflow.connect([
            (inputnode, merge_maps, [(map_name, f'in{i_map}')]),
            (inputnode, merge_masks, [(sparse_3d.get(map_name, 'asl_mask'), f'in{i_map}')]),
        ])  # fmt:skip

    parcellate = pe.MapNode(
        NiftiParcellate(names=columns, min_coverage=min_coverage),
        name='parcellate',
        iterfield=['atlas', 'atlas_labels'],
        mem_gb=mem_gb,
    )
    workflow.connect([
        (inputnode, parcellate, [('atlas_labels_files', 'atlas_labels')]),
        (merge_maps, parcellate, [('out', 'in_files')]),
        (merge_masks, parcellate, [('out', 'masks')]),
        (warp_atlases_to_asl_space, parcellate, [('output_image', 'atlas')]),
    ])  # fmt:skip

    # Describe each column in the sidecars, as BIDS does for tabular files.
    node_meta = {'Description': 'Name of the parcel, from the atlas labels table.'}
    column_meta = {
        'scalarstats': {'node': node_meta},
        'coverage': {'node': node_meta},
    }
    for column, map_name in zip(columns, map_names, strict=True):
        column_meta['scalarstats'][column] = {
            'Description': (
                f'Mean of {column} in each parcel, '
                f'or n/a for parcels with less than {min_coverage * 100}% coverage. '
                f'{BASE_INPUT_FIELDS[map_name]["Description"]}'
            ),
            'Units': BASE_INPUT_FIELDS[map_name]['Units'],
        }
        column_meta['coverage'][column] = {
            'Description': f'Fraction of each parcel in which {column} was estimated.',
        }

    # Each scalarstats column also lists the written map it came from, when there is one.
    # The column metadata goes in meta_dict, rather than the datasinks' keyword arguments,
    # because the datasinks would replace each column's entry (with its Sources) with theirs.
    sourced_maps = [map_name for map_name in map_names if map_name in written_maps]
    if sourced_maps:
        merge_written = pe.Node(
            niu.Merge(len(sourced_maps)),
            name='merge_written',
            run_without_submitting=True,
        )
        for i_map, map_name in enumerate(sourced_maps, start=1):
            workflow.connect([
                (inputnode, merge_written, [(f'{map_name}_file', f'in{i_map}')]),
            ])  # fmt:skip

    for suffix in ('scalarstats', 'coverage'):
        table_meta = pe.MapNode(
            niu.Function(
                function=_add_column_metadata,
                input_names=[
                    'meta_dict',
                    'column_meta',
                    'sourced_columns',
                    'source_files',
                    'dataset_links',
                    'out_dir',
                ],
                output_names=['metadata'],
            ),
            name=f'{suffix}_meta',
            iterfield=['meta_dict'],
            run_without_submitting=True,
        )
        table_meta.inputs.column_meta = column_meta[suffix]
        table_meta.inputs.dataset_links = config.execution.dataset_links
        table_meta.inputs.out_dir = str(output_dir)
        if suffix == 'scalarstats' and sourced_maps:
            table_meta.inputs.sourced_columns = [
                columns[map_names.index(map_name)] for map_name in sourced_maps
            ]
            workflow.connect([(merge_written, table_meta, [('out', 'source_files')])])
        else:
            table_meta.inputs.sourced_columns = []
            table_meta.inputs.source_files = []

        ds_table = pe.MapNode(
            DerivativesDataSink(
                base_directory=config.execution.aslprep_dir,
                check_hdr=False,
                suffix=suffix,
            ),
            name=f'ds_{suffix}',
            iterfield=['atlas', 'in_file', 'meta_dict'],
            run_without_submitting=True,
        )
        workflow.connect([
            (atlas_srcs, table_meta, [('metadata', 'meta_dict')]),
            (inputnode, ds_table, [
                ('source_file', 'source_file'),
                ('atlas_names', 'atlas'),
            ]),
            (parcellate, ds_table, [(suffix, 'in_file')]),
            (table_meta, ds_table, [('metadata', 'meta_dict')]),
        ])  # fmt:skip

    return workflow


def _add_column_metadata(
    meta_dict,
    column_meta,
    sourced_columns,
    source_files,
    dataset_links,
    out_dir,
):
    """Add column descriptions, and each sourced column's Sources, to a table's metadata."""
    from copy import deepcopy

    from aslprep.utils.bids import _get_bidsuris

    if len(sourced_columns) != len(source_files):
        raise ValueError(
            f'Got {len(sourced_columns)} sourced columns but {len(source_files)} files.'
        )

    metadata = deepcopy(meta_dict)
    metadata.update(deepcopy(column_meta))
    for column, source_file in zip(sourced_columns, source_files, strict=True):
        metadata[column]['Sources'] = _get_bidsuris([source_file], dataset_links, out_dir)
    return metadata


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
