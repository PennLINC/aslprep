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
   :template: module.rst

   aslprep.workflows.base
   aslprep.workflows.asl.apply
   aslprep.workflows.asl.base
   aslprep.workflows.asl.cbf
   aslprep.workflows.asl.confounds
   aslprep.workflows.asl.fit
   aslprep.workflows.asl.hmc
   aslprep.workflows.asl.outputs
   aslprep.workflows.asl.plotting
   aslprep.workflows.asl.reference
   aslprep.workflows.asl.util


*************************************
:mod:`aslprep.interfaces`: Interfaces
*************************************

.. currentmodule:: aslprep

.. autosummary::
   :toctree: generated/
   :template: module.rst

   aslprep.interfaces.ants
   aslprep.interfaces.bids
   aslprep.interfaces.cbf
   aslprep.interfaces.confounds
   aslprep.interfaces.parcellation
   aslprep.interfaces.plotting
   aslprep.interfaces.reference
   aslprep.interfaces.reports
   aslprep.interfaces.utility


*******************************
:mod:`aslprep.utils`: Utilities
*******************************

.. currentmodule:: aslprep

.. autosummary::
   :toctree: generated/
   :template: module.rst

   aslprep.utils.asl
   aslprep.utils.atlas
   aslprep.utils.bids
   aslprep.utils.cbf
   aslprep.utils.confounds
   aslprep.utils.misc
   aslprep.utils.plotting
   aslprep.utils.sentry

***********************************
:mod:`aslprep.data`: Data resources
***********************************

.. currentmodule:: aslprep

.. autosummary::
   :toctree: generated/
   :template: module.rst

   aslprep.data
