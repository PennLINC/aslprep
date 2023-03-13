.. include:: links.rst

################
Developers - API
################

*****************************
Internal configuration system
*****************************

.. automodule:: aslprep.config
   :members: from_dict, load, get, dumps, to_filename, init_spaces


***********************************
:mod:`aslprep.workflows`: Workflows
***********************************

.. currentmodule:: aslprep

.. autosummary::
   :toctree: generated/
   :template: function.rst

   workflows.base.init_aslprep_wf
   workflows.base.init_single_subject_wf
   workflows.asl.base.init_asl_preproc_wf
   workflows.asl.cbf.init_cbf_compt_wf
   workflows.asl.cbf.init_cbfqc_compt_wf
   workflows.asl.cbf.init_cbfgeqc_compt_wf
   workflows.asl.cbf.init_cbfplot_wf
   workflows.asl.cbf.init_gecbfplot_wf
   workflows.asl.cbf.init_cbfroiquant_wf
   workflows.asl.cbf.init_gecbf_compt_wf
   workflows.asl.confounds.init_asl_confs_wf
   workflows.asl.confounds.init_carpetplot_wf
   workflows.asl.ge_utils.init_asl_geref_wf
   workflows.asl.ge_utils.init_asl_gereg_wf
   workflows.asl.ge_utils.init_asl_t1_getrans_wf
   workflows.asl.ge_utils.init_asl_gestd_trans_wf
   workflows.asl.gecbf.init_asl_gepreproc_wf
   workflows.asl.hmc.init_asl_hmc_wf
   workflows.asl.outputs.init_asl_derivatives_wf
   workflows.asl.outputs.init_geasl_derivatives_wf
   workflows.asl.registration.init_asl_reg_wf
   workflows.asl.registration.init_asl_t1_trans_wf
   workflows.asl.registration.init_fsl_bbr_wf
   workflows.asl.resampling.init_asl_surf_wf
   workflows.asl.resampling.init_asl_std_trans_wf
   workflows.asl.resampling.init_asl_preproc_trans_wf
   workflows.asl.stc.init_asl_stc_wf
   workflows.asl.t2s.init_asl_t2s_wf


*************************************
:mod:`aslprep.interfaces`: Interfaces
*************************************

.. automodule:: aslprep.interfaces
   :no-members:
   :no-inherited-members:

.. currentmodule:: aslprep

.. autosummary::
   :toctree: generated/
   :template: module.rst

   aslprep.interfaces.bids
   aslprep.interfaces.cbf_computation
   aslprep.interfaces.confounds
   aslprep.interfaces.ge
   aslprep.interfaces.multiecho
   aslprep.interfaces.plotting
   aslprep.interfaces.reports


*******************************
:mod:`aslprep.utils`: Utilities
*******************************

.. automodule:: aslprep.utils
   :no-members:
   :no-inherited-members:

.. currentmodule:: aslprep

.. autosummary::
   :toctree: generated/
   :template: module.rst

   aslprep.utils.bids
   aslprep.utils.confounds
   aslprep.utils.meepi
   aslprep.utils.misc
   aslprep.utils.plotting
   aslprep.utils.qc
   aslprep.utils.sentry
