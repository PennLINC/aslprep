# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Workflows for partial volume maps and partial volume correction."""

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow

from aslprep.interfaces.pvc import LocalRegressionPVC, PartialVolumeMaps
from aslprep.workflows.asl.outputs import lr_pvc_names


def _pick_gm(files):
    return files[0]


def _pick_wm(files):
    return files[1]


def _pick_csf(files):
    return files[2]


def init_asl_pv_wf(
    *,
    metadata,
    effective_resolution='auto',
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
        ``'auto'``, ``'nominal'``, or the effective FWHM in mm (in-plane, in-plane, slice).
    name : :obj:`str`
        Workflow name.

    Inputs
    ------
    t1w_tpms
        Anatomical tissue probability maps ([GM, WM, CSF]).
    asl_mask
        ASL brain mask, on the ASL reference grid.
    aslref2anat_xfm
        Affine transform from the ASL reference to the anatomical image.

    Outputs
    -------
    gm_pv, wm_pv
        Gray and white matter PV maps on the ASL reference grid.
    pv_metadata
        How the maps were built.

    References
    ----------
    .. footbibliography::
    """
    workflow = Workflow(name=name)
    if effective_resolution in ('auto', 'nominal'):
        resolution_str = {
            'auto': (
                'the effective resolution of the acquisition, inferred from its metadata '
                "(see the PV maps' sidecars)"
            ),
            'nominal': 'the voxel size, without additional readout blurring',
        }[effective_resolution]
    else:
        resolution_str = (
            'a user-specified effective resolution of '
            f'{" x ".join(f"{v:g}" for v in effective_resolution)} mm FWHM '
            '(in-plane, in-plane, slice)'
        )
    workflow.__desc__ = f"""\
Gray matter, white matter and CSF partial volume maps were built on the ASL reference grid
by averaging the linearly interpolated tissue probability maps over each ASL voxel and
convolving them to {resolution_str} [@petr2018effects; @chappell2021partial].
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
    name : :obj:`str`
        Workflow name.

    Inputs
    ------
    mean_cbf, mean_cbf_score
        CBF maps.
    m0_calib
        The M0 image CBF was divided by.
    quant_valid
        Where CBF could be quantified.
    asl_mask
        ASL brain mask.
    gm_pv, wm_pv, csf_pv
        Gray matter, white matter and CSF PV maps.

    Outputs
    -------
    The corrected maps, standard errors, validity masks, diagnostics and CSF quality-control
    maps named by :func:`~aslprep.workflows.asl.outputs.lr_pvc_names`, plus ``pvc_metadata``,
    ``pvc_summary`` and ``pvc_qc`` (scalar QC measures, from ``mean_cbf``).
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

    out_fields = ['pvc_metadata', 'pvc_qc']
    for source in sources:
        names = lr_pvc_names(source)
        out_fields += (
            names['cbf']
            + names['se']
            + list(names['valid'].values())
            + names['diagnostics']
            + names['csf']
        )
        out_fields.append(f'pvc_summary{"" if source == "mean_cbf" else "_score"}')

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=sources + ['m0_calib', 'quant_valid', 'asl_mask', 'gm_pv', 'wm_pv', 'csf_pv'],
        ),
        name='inputnode',
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=out_fields), name='outputnode')

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
