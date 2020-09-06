# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Generate T2* map from multi-echo ASL images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_asl_t2s_wf

"""
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from ...interfaces import T2SMap
from ... import config


LOGGER = config.loggers.workflow


# pylint: disable=R0914
def init_asl_t2s_wf(echo_times, mem_gb, omp_nthreads,
                     name='asl_t2s_wf'):
    """
    Combine multiple echos of :abbr:`ME-EPI (multi-echo echo-planar imaging)`.

    This workflow wraps the `tedana`_ `T2* workflow`_ to optimally
    combine multiple echos and derive a T2* map.
    The following steps are performed:

    #. :abbr:`HMC (head motion correction)` on individual echo files.
    #. Compute the T2* map
    #. Create an optimally combined ME-EPI time series

    .. _tedana: https://github.com/me-ica/tedana
    .. _`T2* workflow`: https://tedana.readthedocs.io/en/latest/generated/tedana.workflows.t2smap_workflow.html#tedana.workflows.t2smap_workflow  # noqa

    Parameters
    ----------
    echo_times : :obj:`list`
        list of TEs associated with each echo
    mem_gb : :obj:`float`
        Size of asl file in GB
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    name : :obj:`str`
        Name of workflow (default: ``asl_t2s_wf``)

    Inputs
    ------
    asl_file
        list of individual echo files

    Outputs
    -------
    asl
        the optimally combined time series for all supplied echos

    """
    from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
A T2* map was estimated from the preprocessed asl by fitting to a monoexponential signal
decay model with nonlinear regression, using T2*/S0 estimates from a log-linear
regression fit as initial values.
For each voxel, the maximal number of echoes with reliable signal in that voxel were
used to fit the model.
The calculated T2* map was then used to optimally combine preprocessed asl across
echoes following the method described in [@posse_t2s].
The optimally combined time series was carried forward as the *preprocessed asl*.
"""

    inputnode = pe.Node(niu.IdentityInterface(fields=['asl_file']), name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(fields=['asl']), name='outputnode')

    LOGGER.log(25, 'Generating T2* map and optimally combined ME-EPI time series.')

    t2smap_node = pe.Node(T2SMap(echo_times=echo_times), name='t2smap_node')

    workflow.connect([
        (inputnode, t2smap_node, [('asl_file', 'in_files')]),
        (t2smap_node, outputnode, [('optimal_comb', 'asl')]),
    ])

    return workflow
