# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""

Pre-processing ASL - ASL signal workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: aslprep.workflows.asl.base
.. automodule:: aslprep.workflows.asl.hmc
.. automodule:: aslprep.workflows.asl.stc
.. automodule:: aslprep.workflows.asl.t2s
.. automodule:: aslprep.workflows.asl.registration
.. automodule:: aslprep.workflows.asl.resampling
.. automodule:: aslprep.workflows.asl.confounds
.. automodule:: aslprep.workflows.asl.cbf


"""

from .base import init_asl_preproc_wf
from .gecbf import init_asl_gepreproc_wf
from .hmc import init_asl_hmc_wf
from .stc import init_asl_stc_wf
from .t2s import init_asl_t2s_wf
from .registration import (
    init_asl_t1_trans_wf,
    init_asl_reg_wf,
)
from .resampling import (
    init_asl_std_trans_wf,
    init_asl_surf_wf,
    init_asl_preproc_trans_wf,
)

from .confounds import (
    init_asl_confs_wf
)

from .cbf import (
      init_cbf_compt_wf,
      init_cbfqc_compt_wf,
      init_cbfplot_wf,
      init_gecbfplot_wf,
      init_cbfroiquant_wf,
      init_gecbf_compt_wf,
      init_cbfgeqc_compt_wf)

from .ge_utils import (init_asl_geref_wf, init_asl_gereg_wf,
             init_asl_t1_getrans_wf,init_asl_gestd_trans_wf)
__all__ = [
    'init_asl_confs_wf',
    'init_gecbf_compt_wf',
    'init_asl_t1_getrans_wf',
    'init_asl_geref_wf',
    'init_asl_gereg_wf',
    'init_asl_gestd_trans_wf',
    'init_asl_hmc_wf',
    'init_asl_std_trans_wf',
    'init_asl_preproc_trans_wf',
    'init_asl_reg_wf',
    'init_asl_stc_wf',
    'init_asl_surf_wf',
    'init_asl_t1_trans_wf',
    'init_asl_t2s_wf',
    'init_asl_preproc_wf',
    'init_cbf_compt_wf',
    'init_cbfqc_compt_wf',
    'init_cbfplot_wf',
    'init_cbfroiquant_wf',
    'init_cbfgeqc_compt_wf'
]
