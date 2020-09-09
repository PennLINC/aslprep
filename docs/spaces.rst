.. include:: links.rst

.. _output-spaces:

Standard and nonstandard spaces 
================================
The command line interface of *ASLPPrep* allows resampling of the pre-processed  
and computed CBF maps onto other output spaces.
That is achieved using the ``--output-spaces`` argument, where standard and
nonstandard spaces can be inserted.

.. important::
   *ASLPrep* will reduce the amount of output spaces to just spaces listed in ``--output-spaces``,
   even if other options require resampling the preprocessed data into intermediary spaces.

Standard spaces
"""""""""""""""
By default, *ASLPrep* uses ``MNI152NLin2009cAsym`` as spatial-standardization reference and 
as output space. Valid template identifiers can be used from the 
`TemplateFlow repository <https://github.com/templateflow/templateflow>`__.

However, many users will be interested in utilizing a coarse gridding (typically 2mm isotropic)
of the target template.Such a behavior can be achieved applying modifiers to the template identifier,
separated by a ``:`` character.

For instance, ``--output-spaces MNI152NLin6Asym:res-2` will generate
ASL and CBF on standard space (``MNI152NLin6Asym``) with the template's 2mm 
isotropic resolution.


Nonstandard spaces
""""""""""""""""""
Additionally, ``--output-spaces`` accepts identifiers of spatial references
that do not generate *standardized* coordinate spaces:

  * ``T1w`` or ``anat``: data are resampled into the individual's anatomical
    reference generated with the T1w and  available within the BIDS structure.

  * ``func``,``run``, ``aslref`` or ``sbref`` can be used to
    generate ASL/CBF data in their original grid, after slice-timing,
    head-motion, and susceptibility-distortion corrections as well as CBF
    computation. 
    These keywords are experimental.

Modifiers are not allowed when providing nonstandard spaces.

