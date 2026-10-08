"""Workflows for gradient nonlinearity correction."""

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from aslprep import config
from aslprep.interfaces.bids import DerivativesDataSink
from aslprep.interfaces.gradunwarp import CreateNonlinearityDisplacementMap, MaskWarpDimensions
from aslprep.interfaces.resampling import ResampleSeries
from aslprep.utils.gradwarp import GradwarpPlan, is_displacement_field

_CORRECTION_TEXT = {
    '3D': (
        'Gradient nonlinearity distortions were corrected in all three dimensions, '
        'as the ASL images had not been corrected on the scanner.'
    ),
    '1D': (
        'Gradient nonlinearity distortions were corrected in the through-plane direction only, '
        'as the ASL images had already been corrected in-plane on the scanner '
        '(ImageType: DIS2D).'
    ),
}
_FORCED_CORRECTION_TEXT = {
    '3D': 'Gradient nonlinearity distortions were corrected in all three dimensions.',
    '1D': 'Gradient nonlinearity distortions were corrected in the through-plane direction only.',
}


def gradwarp_boilerplate(plan: GradwarpPlan, jacobian: bool) -> str:
    """Describe the gradient nonlinearity correction for the methods section.

    Parameters
    ----------
    plan : GradwarpPlan
        The run's resolved correction.
    jacobian : bool
        Whether intensities are modulated by the field's Jacobian determinant.

    Returns
    -------
    str
        Methods text (Markdown, with a citation key for TORTOISE).
    """
    if plan.warp_dim is None:
        return (
            'Gradient nonlinearity correction was not applied, as the ASL images had already '
            'been corrected in three dimensions on the scanner (ImageType: DIS3D).'
        )

    texts = _FORCED_CORRECTION_TEXT if plan.basis == 'forced' else _CORRECTION_TEXT
    if is_displacement_field(plan.gradient_file):
        source = 'A displacement field computed from the scanner gradient coefficients was used.'
    else:
        source = (
            "A displacement field was computed from the scanner's gradient coefficients with "
            "TORTOISE's *CreateNonlinearityDisplacementMap* [@tortoisev4]."
        )
    applied = (
        'The field was included in the single resampling step that also applies '
        'head-motion and susceptibility distortion correction'
    )
    if jacobian:
        applied += ', and intensities were modulated by its Jacobian determinant.'
    else:
        applied += '.'
    return f'{texts[plan.warp_dim]} {source} {applied}'


def init_gradwarp_wf(
    *,
    plan: GradwarpPlan,
    asl_file: str,
    report: bool = True,
    name: str = 'gradwarp_wf',
) -> pe.Workflow:
    """Build the gradient nonlinearity displacement field for one ASL run.

    The field is generated on the grid of the motion correction reference, which shares the
    raw ASL series' grid. It is applied downstream, in the same resampling step as head-motion
    and susceptibility distortion correction.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from tempfile import NamedTemporaryFile

            from aslprep.tests.tests import mock_config
            from aslprep.utils.gradwarp import GradwarpPlan
            from aslprep.workflows.asl.gradwarp import init_gradwarp_wf

            coeff_file = NamedTemporaryFile(suffix='.grad', delete=False).name
            with mock_config():
                wf = init_gradwarp_wf(
                    plan=GradwarpPlan(
                        gradient_file=coeff_file,
                        warp_dim='3D',
                        is_ge=False,
                        basis='metadata',
                    ),
                    asl_file='sub-01_asl.nii.gz',
                )

    Parameters
    ----------
    plan
        The resolved correction for this run. ``plan.warp_dim`` must not be None.
    asl_file
        The ASL file, for naming the report.
    report
        Whether to write a before/after report of the correction of ``ref_image``.
    name
        Name of the workflow (default: ``gradwarp_wf``).

    Inputs
    ------
    ref_image
        The 3D motion correction reference image, on whose grid the field is generated.
    slice_ref_image
        The image whose slice axis (``plan.slice_axis``) defines the through-plane direction.
        Only used if ``plan.warp_dim`` is ``'1D'``. Defaults to ``ref_image``.

    Outputs
    -------
    gradwarp_field
        ITK displacement field (LPS, mm) on the reference grid.

    Returns
    -------
    workflow : niworkflows.engine.workflows.LiterateWorkflow
        The workflow.

    Raises
    ------
    ValueError
        If ``plan.warp_dim`` is None, as no field is needed.
    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.reportlets.registration import SimpleBeforeAfterRPT

    if plan.warp_dim is None:
        raise ValueError('No displacement field is needed for images corrected on the scanner.')

    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['ref_image', 'slice_ref_image']),
        name='inputnode',
    )
    slice_ref = pe.Node(
        niu.Function(function=_first_defined, output_names=['out']),
        name='slice_ref',
        run_without_submitting=True,
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=['gradwarp_field']), name='outputnode')

    mask_field = pe.Node(
        MaskWarpDimensions(warp_dim=plan.warp_dim, slice_axis=plan.slice_axis),
        name='mask_field',
    )
    if is_displacement_field(plan.gradient_file):
        mask_field.inputs.in_file = plan.gradient_file
    else:
        make_field = pe.Node(
            CreateNonlinearityDisplacementMap(coeff_file=plan.gradient_file, is_ge=plan.is_ge),
            name='make_field',
        )
        workflow.connect([
            (inputnode, make_field, [('ref_image', 'ref_image')]),
            (make_field, mask_field, [('out_field', 'in_file')]),
        ])  # fmt:skip

    workflow.connect([
        (inputnode, slice_ref, [('slice_ref_image', 'preferred'), ('ref_image', 'fallback')]),
        (slice_ref, mask_field, [('out', 'ref_image')]),
        (mask_field, outputnode, [('out_file', 'gradwarp_field')]),
    ])  # fmt:skip

    if not report:
        return workflow

    # The displacements are small (millimetres, largest at the edges of the field of view),
    # so a before/after reportlet, which alternates the images in place, shows them best.
    corrected_ref = pe.Node(
        ResampleSeries(jacobian=False, gradwarp_jacobian=False),
        name='corrected_ref',
        mem_gb=0.5,
    )
    gradwarp_report = pe.Node(
        SimpleBeforeAfterRPT(before_label='Distorted', after_label='Corrected'),
        name='gradwarp_report',
        mem_gb=0.1,
    )
    ds_gradwarp_report = pe.Node(
        DerivativesDataSink(
            source_file=asl_file,
            base_directory=config.execution.aslprep_dir,
            desc='gradwarp',
            suffix='asl',
            datatype='figures',
            dismiss_entities=('echo',),
        ),
        name='ds_gradwarp_report',
        run_without_submitting=True,
    )

    workflow.connect([
        (inputnode, corrected_ref, [
            ('ref_image', 'in_file'),
            ('ref_image', 'ref_file'),
        ]),
        (mask_field, corrected_ref, [('out_file', 'gradwarp_field')]),
        (inputnode, gradwarp_report, [('ref_image', 'before')]),
        (corrected_ref, gradwarp_report, [('out_file', 'after')]),
        (gradwarp_report, ds_gradwarp_report, [('out_report', 'in_file')]),
    ])  # fmt:skip

    return workflow


def _first_defined(fallback, preferred=None):
    """Return ``preferred`` if it is set, else ``fallback``.

    Parameters
    ----------
    fallback : str
        The value to use if ``preferred`` is not set.
    preferred : str or None, optional
        The preferred value. Nipype leaves undefined inputs out of the call,
        so it defaults to None.

    Returns
    -------
    str
        ``preferred`` if it is set, else ``fallback``.
    """
    return preferred or fallback
