.. include:: links.rst

================
Developers - API
================

Internal configuration system
-----------------------------

.. automodule:: aslprep.config
   :members: from_dict, load, get, dumps, to_filename, init_spaces

:mod:`aslprep.workflows`: Workflows
-----------------------------------

.. automodule:: aslprep.workflows
   :no-members:
   :no-inherited-members:

.. currentmodule:: aslprep

.. autosummary::
   :toctree: generated/
   :template: module.rst

   aslprep.workflows.base
   aslprep.workflows.asl.base
   aslprep.workflows.asl.cbf
   aslprep.workflows.asl.confounds
   aslprep.workflows.asl.ge_utils
   aslprep.workflows.asl.gecbf
   aslprep.workflows.asl.hmc
   aslprep.workflows.asl.outputs
   aslprep.workflows.asl.registration
   aslprep.workflows.asl.resampling
   aslprep.workflows.asl.stc
   aslprep.workflows.asl.t2s

:mod:`aslprep.interfaces`: Interfaces
-------------------------------------

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


:mod:`aslprep.utils`: Utilities
-------------------------------

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
