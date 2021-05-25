# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Phase-difference B0 estimation.

.. _sdc_phasediff :

Phase-difference B0 estimation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The field inhomogeneity inside the scanner (fieldmap) is proportional to the
phase drift between two subsequent :abbr:`GRE (gradient recall echo)`
sequence.

This corresponds to `this section of the BIDS specification
<https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#two-phase-images-and-two-magnitude-images>`__.


"""

from nipype.interfaces import fsl, utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow

from ..interfaces.fmap import Phasediff2Fieldmap, PhaseMap2rads, SubtractPhases
from .gre import init_fmap_postproc_wf, init_magnitude_wf


def init_phdiff_wf(omp_nthreads, name='phdiff_wf'):
    r"""
    Distortion correction of EPI sequences using phase-difference maps.

    Estimates the fieldmap using a phase-difference image and one or more
    magnitude images corresponding to two or more :abbr:`GRE (Gradient Echo sequence)`
    acquisitions.
    The most delicate bit of this workflow is the phase-unwrapping process: phase maps
    are clipped in the range :math:`[0 \dotsb 2\pi )`.
    To find the integer number of offsets that make a region continously smooth with
    its neighbour, FSL PRELUDE is run [Jenkinson2003]_.
    FSL PRELUDE takes wrapped maps in the range 0 to 6.28, `as per the user guide
    <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE/Guide#Step_2_-_Getting_.28wrapped.29_phase_in_radians>`__.
    For the phase-difference maps, recentering back to :math:`[-\pi \dotsb \pi )` is necessary.
    After some massaging and the application of the effective echo spacing parameter,
    the phase-difference maps can be converted into a *B0 field map* in Hz units.
    The `original code was taken from nipype
    <https://github.com/nipy/nipype/blob/0.12.1/nipype/workflows/dmri/fsl/artifacts.py#L514>`_.

    Workflow Graph
        .. workflow ::
            :graph2use: orig
            :simple_form: yes

            from sdcflows.workflows.phdiff import init_phdiff_wf
            wf = init_phdiff_wf(omp_nthreads=1)

    Parameters
    ----------
    omp_nthreads : int
        Maximum number of threads an individual process may use

    Inputs
    ------
    magnitude : list of os.pathlike
        List of path(s) the GRE magnitude maps.
    phasediff : list of tuple(os.pathlike, dict)
        List containing one GRE phase-difference map with its corresponding metadata
        (requires ``EchoTime1`` and ``EchoTime2``), or the phase maps for the two
        subsequent echoes, with their metadata (requires ``EchoTime``).

    Outputs
    -------
    fmap_ref : pathlike
        The average magnitude image, skull-stripped
    fmap_mask : pathlike
        The brain mask applied to the fieldmap
    fmap : pathlike
        The estimated fieldmap in Hz

    References
    ----------
    .. [Jenkinson2003] Jenkinson, M. (2003) Fast, automated, N-dimensional phase-unwrapping
        algorithm. MRM 49(1):193-197. doi:`10.1002/mrm.10354 <10.1002/mrm.10354>`__.

    """
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
A B0-nonuniformity map (or *fieldmap*) was estimated based on a phase-difference map
calculated with a dual-echo GRE (gradient-recall echo) sequence, processed with a
custom workflow of *SDCFlows* inspired by the
[`epidewarp.fsl` script](http://www.nmr.mgh.harvard.edu/~greve/fbirn/b0/epidewarp.fsl)
and further improvements in HCP Pipelines [@hcppipelines].
"""

    inputnode = pe.Node(niu.IdentityInterface(fields=['magnitude', 'phasediff']),
                        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['fmap', 'fmap_ref', 'fmap_mask']), name='outputnode')

    split = pe.MapNode(niu.Function(function=_split, output_names=['map_file', 'meta']),
                       iterfield=['phasediff'], run_without_submitting=True, name='split')

    magnitude_wf = init_magnitude_wf(omp_nthreads=omp_nthreads)

    # phase diff -> radians
    phmap2rads = pe.MapNode(PhaseMap2rads(), name='phmap2rads',
                            iterfield=['in_file'], run_without_submitting=True)
    # FSL PRELUDE will perform phase-unwrapping
    prelude = pe.Node(fsl.PRELUDE(), name='prelude')

    calc_phdiff = pe.Node(SubtractPhases(), name='calc_phdiff',
                          run_without_submitting=True)

    fmap_postproc_wf = init_fmap_postproc_wf(omp_nthreads=omp_nthreads,
                                             fmap_bspline=False)
    compfmap = pe.Node(Phasediff2Fieldmap(), name='compfmap')

    workflow.connect([
        (inputnode, split, [('phasediff', 'phasediff')]),
        (inputnode, magnitude_wf, [('magnitude', 'inputnode.magnitude')]),
        (magnitude_wf, prelude, [('outputnode.fmap_ref', 'magnitude_file'),
                                 ('outputnode.fmap_mask', 'mask_file')]),
        (split, phmap2rads, [('map_file', 'in_file')]),
        (phmap2rads, calc_phdiff, [('out_file', 'in_phases')]),
        (split, calc_phdiff, [('meta', 'in_meta')]),
        (calc_phdiff, prelude, [('phase_diff', 'phase_file')]),
        (prelude, fmap_postproc_wf, [('unwrapped_phase_file', 'inputnode.fmap')]),
        (calc_phdiff, fmap_postproc_wf, [('metadata', 'inputnode.metadata')]),
        (magnitude_wf, fmap_postproc_wf, [
            ('outputnode.fmap_mask', 'inputnode.fmap_mask'),
            ('outputnode.fmap_ref', 'inputnode.fmap_ref')]),
        (fmap_postproc_wf, compfmap, [('outputnode.out_fmap', 'in_file'),
                                      ('outputnode.metadata', 'metadata')]),
        (compfmap, outputnode, [('out_file', 'fmap')]),
        (magnitude_wf, outputnode, [('outputnode.fmap_ref', 'fmap_ref'),
                                    ('outputnode.fmap_mask', 'fmap_mask')]),
    ])

    return workflow


def _split(phasediff):
    return phasediff
