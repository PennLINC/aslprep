# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Apply the estimated fieldmap to perform susceptibility-derived distortion correction."""
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.registration import ANTSApplyTransformsRPT
from niworkflows.func.util import init_enhance_and_skullstrip_bold_wf


def init_sdc_unwarp_wf(omp_nthreads, debug, name='sdc_unwarp_wf'):
    """
    Apply the warping given by a displacements fieldmap.

    This workflow takes in a displacements field through which the
    input reference can be corrected for susceptibility-derived distortion.

    It also calculates a new mask for the input dataset, after the distortions
    have been accounted for.

    .. workflow ::
        :graph2use: orig
        :simple_form: yes

        from sdcflows.workflows.unwarp import init_sdc_unwarp_wf
        wf = init_sdc_unwarp_wf(omp_nthreads=8,
                                debug=False)

    Parameters
    ----------
    omp_nthreads : int
        Maximum number of threads an individual process may use.
    debug : bool
        Run fast configurations of registrations.
    name : str
        Unique name of this workflow.

    Inputs
    ------
    in_warp : os.pathlike
        The :abbr:`DFM (displacements field map)` that corrects for
        susceptibility-derived distortions estimated elsewhere.
    in_reference : os.pathlike
        the reference image to be unwarped.
    in_reference_mask : os.pathlike
        the reference image mask to be unwarped

    Outputs
    -------
    out_reference : str
        the ``in_reference`` after unwarping
    out_reference_brain : str
        the ``in_reference`` after unwarping and skullstripping
    out_warp : str
        the ``in_warp`` field is forwarded for compatibility
    out_mask : str
        mask of the unwarped input file

    """
    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_warp', 'in_reference', 'in_reference_mask']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_reference', 'out_reference_brain', 'out_warp', 'out_mask']),
        name='outputnode')

    unwarp_reference = pe.Node(ANTSApplyTransformsRPT(dimension=3,
                                                      generate_report=False,
                                                      float=True,
                                                      interpolation='LanczosWindowedSinc'),
                               name='unwarp_reference')

    unwarp_mask = pe.Node(ANTSApplyTransformsRPT(
        dimension=3, generate_report=False, float=True,
        interpolation='NearestNeighbor'), name='unwarp_mask')

    enhance_and_skullstrip_bold_wf = init_enhance_and_skullstrip_bold_wf(omp_nthreads=omp_nthreads,
                                                                         pre_mask=True)
    workflow.connect([
        (inputnode, unwarp_reference, [
            ('in_warp', 'transforms'),
            ('in_reference', 'reference_image'),
            ('in_reference', 'input_image')]),
        (inputnode, unwarp_mask, [
            ('in_warp', 'transforms'),
            ('in_reference_mask', 'reference_image'),
            ('in_reference_mask', 'input_image')]),
        (unwarp_reference, enhance_and_skullstrip_bold_wf, [
            ('output_image', 'inputnode.in_file')]),
        (unwarp_mask, enhance_and_skullstrip_bold_wf, [
            ('output_image', 'inputnode.pre_mask')]),
        (inputnode, outputnode, [('in_warp', 'out_warp')]),
        (unwarp_reference, outputnode, [('output_image', 'out_reference')]),
        (enhance_and_skullstrip_bold_wf, outputnode, [
            ('outputnode.mask_file', 'out_mask'),
            ('outputnode.skull_stripped_file', 'out_reference_brain')]),
    ])
    return workflow
