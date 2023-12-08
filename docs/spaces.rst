.. include:: links.rst

###############################
Standard and nonstandard spaces
###############################

The command line interface of *ASLPrep* allows resampling of the pre-processed
and computed CBF maps onto other output spaces.
That is achieved using the ``--output-spaces`` argument, where standard and
nonstandard spaces can be inserted.

This argument is specified in the same manner as *fMRIPrep*.
For more information about this parameter,
please see `fMRIPrep's documentation <https://fmriprep.org/en/stable/spaces.html>`_.

.. important::

  Just note that, instead of ``func``, ``bold``, and ``boldref``,
  *ASLPrep* uses ``asl`` and ``aslref`` for its non-standard spaces.
