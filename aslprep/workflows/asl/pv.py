# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for partial volume maps and partial volume correction."""

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow

from aslprep.interfaces.pvc import (
    LocalRegressionPVC,
    MultiDelayLocalRegressionPVC,
    PartialVolumeMaps,
)
from aslprep.workflows.asl.outputs import lr_pvc_names


def lr_pvc_fields(sources, multi_delay=False):
    """Output fields of :func:`init_lr_pvc_wf` (besides ``pvc_metadata`` and ``pvc_qc``)."""
    fields = []
    for source in sources:
        names = lr_pvc_names(source, multi_delay)
        fields += (
            names['cbf']
            + names['se']
            + list(names['valid'].values())
            + names['diagnostics']
            + names['csf']
            + names['att']
            + list(names['att_valid'].values())
            + names['att_bounds']
            + names['kinetic']
        )
        fields.append(f'pvc_summary{"" if source == "mean_cbf" else "_score"}')
    return fields


def _pick_gm(files):
    return files[0]


def _pick_wm(files):
    return files[1]


def _pick_csf(files):
    return files[2]


def init_asl_pv_wf(
    *,
    metadata,
    effective_resolution='nominal',
    name='asl_pv_wf',
):
    """Build gray and white matter partial volume maps on the ASL reference grid.

    The tissue probability maps are averaged over each ASL voxel and convolved to the
    acquisition's effective resolution, so they match the tissue mixture each voxel measures
    :footcite:p:`petr2018effects,chappell2021partial`.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.pv import init_asl_pv_wf

            wf = init_asl_pv_wf(metadata={'MRAcquisitionType': '2D'})

    Parameters
    ----------
    metadata : :obj:`dict`
        ASL metadata.
    effective_resolution : :obj:`str` or :obj:`list` of :obj:`float`
        ``'nominal'`` (the default), ``'auto'``, or the effective FWHM in mm
        (in-plane, in-plane, slice).
    name : :obj:`str`
        Workflow name.

    Inputs
    ------
    t1w_tpms
        Anatomical tissue probability maps, as a list ordered [GM, WM, CSF].
    asl_mask
        ASL brain mask. Only its grid is used: the PV maps are built on it.
    aslref2anat_xfm
        Affine transform from the ASL reference to the anatomical image (ITK format).

    Outputs
    -------
    gm_pv
        Gray matter PV map on the ASL reference grid (fraction, 0-1).
    wm_pv
        White matter PV map on the ASL reference grid (fraction, 0-1).
    csf_pv
        CSF PV map on the ASL reference grid (fraction, 0-1).
    pv_metadata
        Dictionary describing how the maps were built (effective resolution and its source,
        slice-encoding axis, supersampling, Gaussian sigma), for the maps' sidecars.

    References
    ----------
    .. footbibliography::
    """
    workflow = Workflow(name=name)
    if effective_resolution == 'nominal':
        resolution_str = (
            'The effective resolution was taken to be the voxel size, '
            'so readout blurring beyond the voxel was not modelled.'
        )
    elif effective_resolution == 'auto':
        resolution_str = (
            'The maps were then convolved to the effective resolution of the acquisition, '
            "inferred from its metadata (see the PV maps' sidecars)."
        )
    else:
        resolution_str = (
            'The maps were then convolved to a user-specified effective resolution of '
            f'{" x ".join(f"{v:g}" for v in effective_resolution)} mm FWHM '
            '(in-plane, in-plane, slice).'
        )
    # Start a new paragraph: niworkflows joins the parent's description and this one directly.
    workflow.__desc__ = f"""
Gray matter, white matter and CSF partial volume maps were built on the ASL reference grid
by averaging the linearly interpolated tissue probability maps over each ASL voxel
[@petr2018effects; @chappell2021partial].
{resolution_str}
Blurring from head motion and distortion correction resampling was not modelled.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['t1w_tpms', 'asl_mask', 'aslref2anat_xfm']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['gm_pv', 'wm_pv', 'csf_pv', 'pv_metadata']),
        name='outputnode',
    )

    pv_maps = pe.Node(
        PartialVolumeMaps(metadata=metadata, effective_resolution=effective_resolution),
        name='pv_maps',
        mem_gb=1.0,
    )
    workflow.connect([
        (inputnode, pv_maps, [
            (('t1w_tpms', _pick_gm), 'gm_tpm'),
            (('t1w_tpms', _pick_wm), 'wm_tpm'),
            (('t1w_tpms', _pick_csf), 'csf_tpm'),
            ('asl_mask', 'reference'),
            ('aslref2anat_xfm', 'aslref2anat_xfm'),
        ]),
        (pv_maps, outputnode, [
            ('gm_pv', 'gm_pv'),
            ('wm_pv', 'wm_pv'),
            ('csf_pv', 'csf_pv'),
            ('pv_metadata', 'pv_metadata'),
        ]),
    ])  # fmt:skip

    return workflow


def init_lr_pvc_wf(
    *,
    metadata,
    sources,
    kernel_size=(5, 5, 1),
    multi_delay=False,
    num_threads=1,
    use_gpu=False,
    name='lr_pvc_wf',
):
    """Correct CBF maps for partial volume effects by local linear regression.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.workflows.asl.pv import init_lr_pvc_wf

            wf = init_lr_pvc_wf(metadata={'MRAcquisitionType': '2D'}, sources=['mean_cbf'])

    Parameters
    ----------
    metadata : :obj:`dict`
        ASL metadata.
    sources : :obj:`list` of :obj:`str`
        The CBF maps to correct: ``mean_cbf``, and ``mean_cbf_score`` with SCORE.
    kernel_size : :obj:`tuple` of :obj:`int`
        Kernel size in voxels, ordered (in-plane, in-plane, slice).
    multi_delay : :obj:`bool`
        Whether the data have several delays. If so, each delay is unmixed and each tissue curve
        is fit with the general kinetic model (spec Section 5.7).
    num_threads, use_gpu
        Threads and GPU use for the multi-delay kinetic fits.
    name : :obj:`str`
        Workflow name.

    Inputs
    ------
    mean_cbf
        Mean CBF map to correct. Present when ``'mean_cbf'`` is in ``sources``.
        For multi-delay data it is the uncorrected CBF, used for the CSF quality-control
        template.
    mean_cbf_score
        SCORE CBF map to correct. Present when ``'mean_cbf_score'`` is in ``sources``
        (single-delay data only).
    deltam, asl_metadata
        Multi-delay data only: the delta-M series that is unmixed, and its metadata (with
        the delays of the selected volumes).
    att, abat, abv
        Multi-delay data only: the uncorrected kinetic maps, for the CSF quality-control
        template.
    m0_calib
        The M0 image the CBF maps were divided by (ComputeCBF's calibration image).
    quant_valid
        Mask of voxels where CBF could be quantified.
    asl_mask
        ASL brain mask.
    gm_pv, wm_pv, csf_pv
        Gray matter, white matter and CSF PV maps on the ASL reference grid.

    Outputs
    -------
    mean_cbf_gm_lr, mean_cbf_wm_lr
        Gray and white matter CBF after correction of ``mean_cbf``.
        Zero where not estimated.
    cbf_gm_se_lr, cbf_wm_se_lr
        Standard errors of the corrected gray and white matter CBF.
    gm_valid_lr, wm_valid_lr
        Masks of voxels where each tissue's CBF was estimated.
    status_lr
        Fit status of each voxel (0: outside the mask; 1: two-tissue fit; 2: GM only;
        3: WM only; 4: not identifiable).
    condition_lr
        Condition number of each voxel's normalized design.
    gm_noise_lr, wm_noise_lr
        Noise amplification of the gray and white matter estimates.
    rmse_lr
        Root mean square residual of each voxel's regression.
    mean_cbf_csf_lr, cbf_csf_se_lr, status_csf_lr, csf_noise_lr
        Apparent CSF CBF from the separate three-tissue quality-control fit, its standard
        error, fit status (0: outside the mask; 1: GM, WM and CSF; 2: GM and CSF; 3: WM and
        CSF; 4: not identifiable; 5: CSF only) and noise amplification.
    pvc_summary
        Dictionary with the fraction of brain voxels in each fit status.
    mean_cbf_gm_lr_score, ..., csf_noise_lr_score, pvc_summary_score
        The same outputs for ``mean_cbf_score``, with the ``_score`` suffix.
        Present when ``'mean_cbf_score'`` is in ``sources``.

    For multi-delay data, the outputs are the single-delay ``mean_cbf`` outputs above
    without the standard errors (``cbf_gm_se_lr``, ``cbf_wm_se_lr``, ``cbf_csf_se_lr``),
    with ``rmse_lr`` pooled over volumes and ``mean_cbf_csf_lr`` the signed projection of the
    CSF curve, plus:

    att_gm_lr, att_wm_lr
        Gray and white matter arterial transit time from each tissue curve's kinetic fit.
        Zero where not estimated.
    gm_att_valid_lr, wm_att_valid_lr
        Masks of voxels where each tissue's ATT was estimated.
    att_gm_lower_lr, att_gm_upper_lr, att_wm_lower_lr, att_wm_upper_lr
        The shortest and longest ATT that fit each tissue curve as well as the estimate.
    kinstatus_gm_lr, kinstatus_wm_lr
        Status of each tissue's kinetic fit (0: not fitted; 1: fitted; 2: CBF at its lower
        bound; 3: CBF at its upper bound; 4: input not finite).
    kinres_gm_lr, kinres_wm_lr
        Relative residual of each tissue's kinetic fit.
    csf_kinres_lr
        Relative residual of the CSF curve against its quality-control template.

    pvc_metadata
        Dictionary of the correction's settings, for the outputs' sidecars.
    pvc_qc
        Dictionary of scalar quality-control measures (from ``mean_cbf``) for the CBF QC
        file: ``lr_csf_cbf_median``, ``lr_csf_cbf_iqr``, ``lr_csf_n_voxels`` and
        ``lr_fraction_two_tissue``.
    """
    workflow = Workflow(name=name)
    kernel_str = ' x '.join(str(k) for k in kernel_size)
    workflow.__desc__ = f"""\
Partial volume correction was additionally performed by local linear regression
[@asllani2008regression].
Within a {kernel_str} voxel kernel (in-plane, in-plane, slice) around each voxel,
the M0-independent perfusion signal (CBF multiplied by the M0 it was quantified with)
was regressed on the gray and white matter partial volume maps, without an intercept,
assuming CSF has no perfusion and tissue perfusion is constant within the kernel.
The tissue signals were then divided by the voxel's M0, the same calibration used for the
uncorrected CBF.
Voxels were fit for both tissues when each had at least half a voxel of support in the kernel
and the normalized design was well conditioned (condition number at most 100),
and for one tissue only when the other tissue's contamination of that estimate was below 5%;
otherwise, no estimate was made.
Tissue estimates were retained where that tissue's partial volume was at least 10%.
The uncorrected CBF maps were retained, as partial volume correction is recommended as an
additional analysis [@chappell2021partial].
For quality control, a separate regression on the gray matter, white matter and CSF partial
volume maps estimated an apparent CSF CBF, which should be near zero;
CSF was left out of the gray and white matter fit.
"""

    out_fields = ['pvc_metadata', 'pvc_qc'] + lr_pvc_fields(sources, multi_delay)

    in_fields = sources + ['m0_calib', 'quant_valid', 'asl_mask', 'gm_pv', 'wm_pv', 'csf_pv']
    if multi_delay:
        in_fields += ['deltam', 'asl_metadata', 'att', 'abat', 'abv']
    inputnode = pe.Node(niu.IdentityInterface(fields=in_fields), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=out_fields), name='outputnode')

    if multi_delay:
        workflow.__desc__ = f"""\
Partial volume correction was additionally performed by local linear regression
[@asllani2008regression], extended to multi-delay data.
Within a {kernel_str} voxel kernel (in-plane, in-plane, slice) around each voxel,
the delta-M of each volume was regressed on the gray and white matter partial volume maps,
without an intercept, assuming CSF has no perfusion and that each tissue's signal is constant
within the kernel. This gave a gray and a white matter curve across delays,
each of which was then fit with the same general kinetic model as the uncorrected data,
calibrated by the voxel's M0 and with a non-negative arterial nuisance term,
giving tissue CBF and arterial transit time.
Voxels were fit for both tissues when each had at least half a voxel of support in the kernel
and the normalized design was well conditioned (condition number at most 100),
and for one tissue only when the other tissue's contamination of that estimate was below 5%.
Tissue estimates were retained where that tissue's partial volume was at least 10%,
and tissue arterial transit times only where the tissue's CBF was positive.
Arterial signal and variation of the kinetics within the kernel are not modelled by the
unmixing; the residual of each tissue's kinetic fit is reported.
The uncorrected CBF maps were retained, as partial volume correction is recommended as an
additional analysis [@chappell2021partial].
For quality control, a separate regression on the gray matter, white matter and CSF partial
volume maps estimated a CSF curve, whose signed projection onto each voxel's uncorrected
kinetic curve gave an apparent CSF CBF, which should be near zero.
"""
        names = lr_pvc_names('mean_cbf', multi_delay=True)
        lr_pvc = pe.Node(
            MultiDelayLocalRegressionPVC(
                kernel_size=list(kernel_size),
                num_threads=num_threads,
                use_gpu=use_gpu,
            ),
            name='lr_pvc',
            mem_gb=1 + 0.1 * num_threads,
            n_procs=num_threads,
        )
        gm_cbf, wm_cbf = names['cbf']
        gm_att, wm_att = names['att']
        gm_lower, gm_upper, wm_lower, wm_upper = names['att_bounds']
        gm_kinstatus, wm_kinstatus, gm_kinres, wm_kinres = names['kinetic']
        status, condition, gm_noise, wm_noise, rmse = names['diagnostics']
        csf_cbf, csf_status, csf_noise, csf_kinres = names['csf']
        workflow.connect([
            (inputnode, lr_pvc, [
                ('deltam', 'deltam'),
                ('asl_metadata', 'metadata'),
                ('m0_calib', 'm0_calib'),
                ('quant_valid', 'quant_valid'),
                ('asl_mask', 'mask'),
                ('gm_pv', 'gm_pv'),
                ('wm_pv', 'wm_pv'),
                ('csf_pv', 'csf_pv'),
                ('mean_cbf', 'mean_cbf'),
                ('att', 'att'),
                ('abat', 'abat'),
                ('abv', 'abv'),
            ]),
            (lr_pvc, outputnode, [
                ('gm_cbf', gm_cbf),
                ('wm_cbf', wm_cbf),
                ('gm_valid', names['valid'][gm_cbf]),
                ('wm_valid', names['valid'][wm_cbf]),
                ('gm_att', gm_att),
                ('wm_att', wm_att),
                ('gm_att_valid', names['att_valid'][gm_att]),
                ('wm_att_valid', names['att_valid'][wm_att]),
                ('gm_att_lower', gm_lower),
                ('gm_att_upper', gm_upper),
                ('wm_att_lower', wm_lower),
                ('wm_att_upper', wm_upper),
                ('gm_kinstatus', gm_kinstatus),
                ('wm_kinstatus', wm_kinstatus),
                ('gm_kinres', gm_kinres),
                ('wm_kinres', wm_kinres),
                ('status', status),
                ('condition', condition),
                ('gm_noise', gm_noise),
                ('wm_noise', wm_noise),
                ('rmse', rmse),
                ('csf_cbf', csf_cbf),
                ('csf_status', csf_status),
                ('csf_noise', csf_noise),
                ('csf_kinres', csf_kinres),
                ('summary', 'pvc_summary'),
                ('pvc_metadata', 'pvc_metadata'),
                ('qc', 'pvc_qc'),
            ]),
        ])  # fmt:skip
        return workflow

    for source in sources:
        names = lr_pvc_names(source)
        suffix = '' if source == 'mean_cbf' else '_score'
        lr_pvc = pe.Node(
            LocalRegressionPVC(metadata=metadata, kernel_size=list(kernel_size)),
            name=f'lr_pvc{suffix}',
            mem_gb=0.5,
        )
        gm_cbf, wm_cbf = names['cbf']
        gm_se, wm_se = names['se']
        gm_valid, wm_valid = names['valid'][gm_cbf], names['valid'][wm_cbf]
        status, condition, gm_noise, wm_noise, rmse = names['diagnostics']
        csf_cbf, csf_se, csf_status, csf_noise = names['csf']
        workflow.connect([
            (inputnode, lr_pvc, [
                (source, 'cbf'),
                ('m0_calib', 'm0_calib'),
                ('quant_valid', 'quant_valid'),
                ('asl_mask', 'mask'),
                ('gm_pv', 'gm_pv'),
                ('wm_pv', 'wm_pv'),
                ('csf_pv', 'csf_pv'),
            ]),
            (lr_pvc, outputnode, [
                ('gm_cbf', gm_cbf),
                ('wm_cbf', wm_cbf),
                ('gm_se', gm_se),
                ('wm_se', wm_se),
                ('gm_valid', gm_valid),
                ('wm_valid', wm_valid),
                ('status', status),
                ('condition', condition),
                ('gm_noise', gm_noise),
                ('wm_noise', wm_noise),
                ('rmse', rmse),
                ('csf_cbf', csf_cbf),
                ('csf_se', csf_se),
                ('csf_status', csf_status),
                ('csf_noise', csf_noise),
                ('summary', f'pvc_summary{suffix}'),
            ]),
        ])  # fmt:skip
        if source == 'mean_cbf':
            workflow.connect([
                (lr_pvc, outputnode, [('pvc_metadata', 'pvc_metadata'), ('qc', 'pvc_qc')]),
            ])  # fmt:skip

    return workflow
