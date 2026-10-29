"""Workflows that measure tissue T1s from a quantitative T1 map."""

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.reportlets.registration import (
    SimpleBeforeAfterRPT as SimpleBeforeAfter,
)

from aslprep import config
from aslprep.interfaces.ants import ApplyTransforms
from aslprep.interfaces.bids import DerivativesDataSink
from aslprep.interfaces.t1map import MeasureTissueT1, PrepareT1Map, UnmixTissueR1


def init_t1map_wf(*, t1map, identity, output_dir, field_strength=None, name='t1map_wf'):
    """Measure the T1 of gray and white matter from the participant's T1 map.

    The T1 map is aligned with the anatomical image (``t1w_preproc``), resampled onto it, and the
    T1 of each tissue is the median over its pure voxels
    (see :func:`~aslprep.utils.t1map.tissue_t1_medians`).
    The values replace the defaults of the tissue-T1 kinetic model at run time
    (see :class:`~aslprep.interfaces.cbf.ResolveRuntimeCalibration`).
    With ``--t1map-method local-regression``, the R1 of each tissue is also estimated in each
    voxel (see :func:`~aslprep.utils.t1map.local_regression_r1`), and the medians are used only
    where it is not.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from aslprep.tests.tests import mock_config
            from aslprep.workflows.t1map import init_t1map_wf

            with mock_config():
                wf = init_t1map_wf(
                    t1map={
                        't1map': 'sub-01_T1map.nii.gz',
                        'metadata': {'Units': 's'},
                        'uri': 'bids:raw:sub-01/anat/sub-01_T1map.nii.gz',
                        'companion': 'sub-01_UNIT1.nii.gz',
                    },
                    identity=False,
                    output_dir='.',
                )

    Parameters
    ----------
    t1map : :obj:`dict`
        The T1 map, from :func:`~aslprep.utils.bids.collect_t1map`.
    identity : :obj:`bool`
        Whether the T1 map is already in the anatomical image's space: its companion image is the
        only T1-weighted input of the anatomical workflow, whose reference then shares its scanner
        space. Otherwise, the companion image (or, without one, the map's R1 = 1/T1) is rigidly
        registered to the anatomical image.
    output_dir : :obj:`str`
        The output directory.
    field_strength : :obj:`float` or None
        Field strength, in tesla, of the ASL runs that use the map, to check the measured T1s
        against the defaults when the map's sidecar does not give its own.
    name : :obj:`str`
        Workflow name.

    Inputs
    ------
    t1w_preproc, t1w_mask
        The anatomical image and its brain mask.
    t1w_tpms
        Gray matter, white matter and CSF probability maps of the anatomical image.

    Outputs
    -------
    updates
        ``GrayMatterT1`` and ``WhiteMatterT1`` entries for the calibration specification.
    warnings
        Warnings about the measured values.
    t1map2anat_xfm
        The transform from the T1 map to the anatomical image.
    gm_r1, wm_r1, r1_summary
        With ``--t1map-method local-regression``, the gray and white matter R1 maps on the
        anatomical grid, and how they were made.

    The workflow also writes a reportlet (``desc-coreg_T1map.svg``) comparing the anatomical image
    with the map's R1 (1/T1) after the same resampling as the T1s, with the white matter
    probability map's contour, so the alignment can be checked.
    """
    from aslprep.data import load as load_data

    workflow = Workflow(name=name)
    companion = t1map['companion']
    if identity:
        registration = (
            'which shared the scanner space of the anatomical image through the T1-weighted '
            f'image on its grid ({companion}), the only input of the anatomical workflow'
        )
        registration_desc = (
            'which shared the scanner space of the anatomical reference, as the T1-weighted '
            'image on its grid was the only anatomical input'
        )
    elif companion:
        registration = (
            'rigidly registered to the anatomical image through the T1-weighted image on its grid '
            f'({companion})'
        )
        registration_desc = (
            'which was aligned with the anatomical reference by rigidly registering the '
            'T1-weighted image on its grid to it with `antsRegistration` [@ants]'
        )
    else:
        config.loggers.workflow.warning(
            f'No T1-weighted image shares the grid of the T1 map {t1map["t1map"]}, so the map '
            'itself (as R1 = 1/T1) is registered to the anatomical image. '
            'Check the registration.'
        )
        registration = 'rigidly registered to the anatomical image through its R1 (1/T1) map'
        registration_desc = (
            'which was aligned with the anatomical reference by rigidly registering its R1 (1/T1) '
            'map to it with `antsRegistration` [@ants]'
        )

    workflow.__desc__ = f"""\
The T1s of gray and white matter in the tissue-T1 kinetic model were measured from the
participant's T1 map, {registration_desc}.
The map was resampled onto the anatomical reference with linear interpolation, and each tissue's
T1 was the median over the voxels with a tissue probability of at least 95% and a CSF probability
below 5%, eroded by one voxel.
"""
    if config.workflow.t1map_method == 'local-regression':
        from aslprep.utils.t1map import MAX_CONTAMINATION, MAX_R1_RATIO

        workflow.__desc__ += f"""\
The T1s of gray and white matter were also estimated in each voxel, assuming that a voxel's R1
(1/T1) is the sum of its tissues' R1s weighted by their probabilities (fast exchange): in a
{config.workflow.t1map_kernel_size:g} mm cube around each voxel, R1 was regressed on the gray
matter, white matter and CSF probability maps without an intercept, by local linear regression
[@asllani2008regression].
A tissue was left out of a fit only when its contamination of the gray and white matter estimates
was at most {MAX_CONTAMINATION:g}.
Where a tissue was not in the fit, or its estimate was more than a factor of {MAX_R1_RATIO:g} from
its median, the median was used.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['t1map', 'companion', 't1w_preproc', 't1w_mask', 't1w_tpms'],
        ),
        name='inputnode',
    )
    # Set here, rather than on the nodes, so that files are only checked when the workflow runs
    inputnode.inputs.t1map = t1map['t1map']
    inputnode.inputs.companion = companion
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['updates', 'warnings', 't1map2anat_xfm', 'gm_r1', 'wm_r1', 'r1_summary'],
        ),
        name='outputnode',
    )

    prepare = pe.Node(PrepareT1Map(metadata=t1map['metadata']), name='prepare_t1map')
    workflow.connect([(inputnode, prepare, [('t1map', 'in_file')])])
    xfm = pe.Node(niu.IdentityInterface(fields=['xfm']), name='t1map2anat_xfm')
    if identity:
        from smriprep import load_data as load_smriprep_data

        xfm.inputs.xfm = str(load_smriprep_data('itkIdentityTransform.txt'))
    else:
        from nipype.interfaces.ants import Registration

        # The ANTs random seed is set in the environment (config.seeds).
        omp_nthreads = config.nipype.omp_nthreads or 1
        register = pe.Node(
            Registration(
                from_file=str(load_data('t1map_rigid_registration.json')),
                num_threads=omp_nthreads,
            ),
            name='register_t1map',
            n_procs=omp_nthreads,
        )
        if companion:
            workflow.connect([(inputnode, register, [('companion', 'moving_image')])])
        else:
            workflow.connect([(prepare, register, [('r1', 'moving_image')])])

        workflow.connect([
            (inputnode, register, [
                ('t1w_preproc', 'fixed_image'),
                ('t1w_mask', 'fixed_image_mask'),
            ]),
            (register, xfm, [(('forward_transforms', _first), 'xfm')]),
        ])  # fmt:skip

    resample = pe.Node(
        ApplyTransforms(interpolation='Linear', float=True, input_image_type=0),
        name='resample_t1map',
    )
    measure = pe.Node(
        MeasureTissueT1(
            t1map_name=t1map['uri'],
            registration=registration,
            field_strength=t1map['metadata'].get('MagneticFieldStrength', field_strength),
        ),
        name='measure_tissue_t1',
    )
    ds_xfm = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            to='T1w',
            mode='image',
            suffix='xfm',
            extension='.txt',
            # A selected map may carry the space and desc of its pipeline (e.g., space-T1w)
            dismiss_entities=('desc', 'space'),
            **{'from': 'T1map'},
        ),
        name='ds_t1map2anat_xfm',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # Registration reportlet: the anatomical reference and the map's R1 (1/T1), which has a
    # T1-weighted contrast, after resampling with the same transform as the T1s
    resample_r1 = pe.Node(
        ApplyTransforms(interpolation='Linear', float=True, input_image_type=0),
        name='resample_r1',
    )
    t1map_report = pe.Node(
        SimpleBeforeAfter(before_label='T1w', after_label='T1 map (R1)'),
        name='t1map_report',
        mem_gb=0.1,
    )
    ds_t1map_report = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='coreg',
            suffix='T1map',
            datatype='figures',
            dismiss_entities=('desc', 'space'),
        ),
        name='ds_t1map_report',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([
        (prepare, resample_r1, [('r1', 'input_image')]),
        (inputnode, resample_r1, [('t1w_preproc', 'reference_image')]),
        (xfm, resample_r1, [('xfm', 'transforms')]),
        (inputnode, t1map_report, [
            ('t1w_preproc', 'before'),
            (('t1w_tpms', _second), 'wm_seg'),
        ]),
        (resample_r1, t1map_report, [('output_image', 'after')]),
        (t1map_report, ds_t1map_report, [('out_report', 'in_file')]),
        (inputnode, ds_t1map_report, [('t1map', 'source_file')]),
    ])  # fmt:skip

    workflow.connect([
        (prepare, resample, [('t1', 'input_image')]),
        (inputnode, resample, [('t1w_preproc', 'reference_image')]),
        (xfm, resample, [('xfm', 'transforms')]),
        (resample, measure, [('output_image', 't1')]),
        (prepare, measure, [('units', 'units')]),
        (inputnode, measure, [('t1w_tpms', 'tpms')]),
        (xfm, ds_xfm, [('xfm', 'in_file')]),
        (inputnode, ds_xfm, [('t1map', 'source_file')]),
        (measure, outputnode, [('warnings', 'warnings')]),
        (ds_xfm, outputnode, [('out_file', 't1map2anat_xfm')]),
    ])  # fmt:skip

    if config.workflow.t1map_method == 'local-regression':
        unmix = pe.Node(
            UnmixTissueR1(kernel_mm=config.workflow.t1map_kernel_size),
            name='unmix_tissue_r1',
            # About 1.6 GB at 1 mm, with the inputs, since the box sums are taken in slabs
            mem_gb=3,
        )
        workflow.connect([
            (resample, unmix, [('output_image', 't1')]),
            (inputnode, unmix, [
                ('t1w_tpms', 'tpms'),
                ('t1w_mask', 'mask'),
            ]),
            (measure, unmix, [
                ('medians', 'medians'),
                ('updates', 'updates'),
            ]),
            (unmix, outputnode, [
                ('updates', 'updates'),
                ('gm_r1', 'gm_r1'),
                ('wm_r1', 'wm_r1'),
                ('summary', 'r1_summary'),
            ]),
        ])  # fmt:skip
    else:
        workflow.connect([(measure, outputnode, [('updates', 'updates')])])

    return workflow


def _first(in_list):
    return in_list[0]


def _second(in_list):
    return in_list[1]
