.. include:: links.rst

.. _output-spaces:

Defining standard and nonstandard spaces where data will be resampled
=====================================================================
The command line interface of *ASLPPrep* allows resampling of the pre-processed data
onto other output spaces.
That is achieved using the ``--output-spaces`` argument, where standard and
nonstandard spaces can be inserted.

.. important::
   *ASLPrep* will reduce the amount of output spaces to just spaces listed in ``--output-spaces``,
   even if other options require resampling the preprocessed data into intermediary spaces.


.. _TemplateFlow:

*TemplateFlow*
""""""""""""""
*TemplateFlow* is a software library and a repository of neuroimaging templates
that allows end-user applications such as *fMRIPrep*, *QSIPrep*, as well as *ASLPrep*,
to flexibly query and pull template and atlas information.
In other words, *TemplateFlow* enables *ASLPrep* to access a wide range
of templates (and also custom templates, see below), and is central to defining 
*ASLPPrep*'s interface regarding templates and atlases.
For more general information about *TemplateFlow*, visit
`TemplateFlow.org <https://www.templateflow.org>`__.


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

Pre-processing blocks depending on standard templates
"""""""""""""""""""""""""""""""""""""""""""""""""""""
Some modules of the pipeline (e.g., the generation of HCP compatible *grayordinates* 
files, or the *fieldmap-less* distortion correction) operate in specific template spaces.
Selecting those modules to be included (using any of the following flags:
``--cifti-outputs``, ``--use-syn-sdc``) will modify the list of
output spaces to include the space identifiers they require, should the
identifier not be found within the ``--output-spaces`` list already.
In other words, running *ASLPrep* with ``--output-spaces MNI152NLin6Asym:res-2
--use-syn-sdc`` will expand the list of output spaces to 
``MNI152NLin6Asym:res-2 MNI152NLin2009cAsym``.
