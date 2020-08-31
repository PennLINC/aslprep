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
of templates (and also custom templates, see below), and is central to defining *ASLPPrep*'s interface regarding templates and atlases.
For more general information about *TemplateFlow*, visit
`TemplateFlow.org <https://www.templateflow.org>`__.


Standard spaces
"""""""""""""""
When using *ASLPrep* in a workflow that will investigate effects that span 
analytical groupings, neuroimagers typically resample their data on to a standard,
stereotactic coordinate system.
The most extended standard space for functaional data analysis is generally MNI.
By default, *ASLPrep* uses ``MNI152NLin2009cAsym`` as spatial-standardization reference.
Valid template identifiers (``MNI152NLin6Asym``, ``MNI152NLin2009cAsym``, etc.) come from
the `TemplateFlow repository <https://github.com/templateflow/templateflow>`__.
For instance, to instruct *ASLPrep* to use the MNI template brain distributed with
FSL as a coordinate reference, the option will read as follows: ``--output-spaces MNI152NLin6Asym``.


Therefore, *ASLPrep* will run nonlinear registration processes against the template
T1w image corresponding to all the standard spaces supplied with the argument
``--output-spaces``.
By default, *ASLPrep* will resample the preprocessed data on those spaces (labeling the
corresponding outputs with the `space-<template-identifier>` BIDS entity) but keeping
the original resolution of the BOLD data to produce smaller files which are more consistent with
the original data gridding.

However, many users will be interested in utilizing a coarse gridding (typically 2mm isotropic)
of the target template.
Such a behavior can be achieved applying modifiers to the template identifier, separated by
a ``:`` character.
For instance, ``--output-spaces MNI152NLin6Asym:res-2 MNI152NLin2009cAsym`` will generate
pre-processed BOLD 4D files on two standard spaces (``MNI152NLin6Asym``,
and ``MNI152NLin2009cAsym``) with the template's 2mm isotropic resolution for
the data on ``MNI152NLin6Asym`` space and the original BOLD resolution
(e.g., 2x2x2.5 [mm]) for the case of ``MNI152NLin2009cAsym``.
This is equivalent to saying
``--output-spaces MNI152NLin6Asym:res-2 MNI152NLin2009cAsym:res-native``.

Other possible modifiers include the ``cohort`` selector.
For instance, ``--output-spaces MNIPediatricAsym:res-1:cohort-2`` selects
the resolution ``1`` of ``cohort-2`` which, for the ``MNIPediatricAsym``
template, corresponds to the `pre-puberty phase
<https://github.com/templateflow/tpl-MNIPediatricAsym/blob/bcf77616f547f327ee53c01dadf689ab6518a097/template_description.json#L22-L26>`__
(4.5--8.5 years old).

Space modifiers such as ``res`` are combinatorial. The argument
``--output-spaces MNIPediatricAsym:cohort-1:cohort-2:res-native:res-1`` will
generate conversions for the following combinations:

  * cohort ``1`` and "native" resolution (meaning, the original BOLD resolution)
  * cohort ``1`` and resolution ``1`` of the template
  * cohort ``2`` and "native" resolution (meaning, the original BOLD resolution)
  * cohort ``2`` and resolution ``1`` of the template

Please note that the selected resolutions specified must exist within TemplateFlow.

When specifying surface spaces (e.g., ``fsaverage``), the legacy identifiers from
FreeSurfer will be supported (e.g., ``fsaverage5``) although the use of the density
modifier would be preferred (e.g., ``fsaverage:den-10k`` for ``fsaverage5``).

Custom standard spaces
""""""""""""""""""""""
To make your custom templates visible by *ASLPrep*, and usable via
the ``--output-spaces`` argument, store your template under
*TemplateFlow*'s home directory.
The default *TemplateFlow* home directory is ``$HOME/.cache/templateflow``
but the path can be arbitrarily changed by setting
the ``$TEMPLATEFLOW_HOME`` environment variable.
A minimal example of the necessary files for a template called
``MyCustom`` (and therefore callable via ``--output-spaces MyCustom``)
follows::

  $TEMPLATEFLOW_HOME/
      tpl-MyCustom/
          template_description.json
          tpl-MyCustom_res-1_T1w.nii.gz
          tpl-MyCustom_res-1_desc-brain_mask.nii.gz
          tpl-MyCustom_res-2_T1w.nii.gz
          tpl-MyCustom_res-2_desc-brain_mask.nii.gz

For further information about how custom templates must be organized and their
corresponding naming conventions, please check `the TemplateFlow tutorials
<https://www.templateflow.org/python-client/tutorials.html>`__.

Nonstandard spaces
""""""""""""""""""
Additionally, ``--output-spaces`` accepts identifiers of spatial references
that do not generate *standardized* coordinate spaces:

  * ``T1w`` or ``anat``: data are resampled into the individual's anatomical
    reference generated with the T1w and T2w images available within the
    BIDS structure.
  * ``fsnative``: similarly to the ``anat`` space for volumetric references,
    including the ``fsnative`` space will instruct *fMRIPrep* to sample the
    original BOLD data onto FreeSurfer's reconstructed surfaces for this
    individual.
  * ``func``, ``bold``, ``run``, ``boldref`` or ``sbref`` can be used to
    generate BOLD data in their original grid, after slice-timing,
    head-motion, and susceptibility-distortion corrections.
    These keywords are experimental.

Modifiers are not allowed when providing nonstandard spaces.

Pre-processing blocks depending on standard templates
""""""""""""""""""""""""""""""""""""""""""""""""""""
Some modules of the pipeline (e.g., the generation of HCP compatible *grayordinates* 
files, or the *fieldmap-less* distortion correction) operate in specific template spaces.
Selecting those modules to be included (using any of the following flags:
``--cifti-outputs``, ``--use-syn-sdc``) will modify the list of
output spaces to include the space identifiers they require, should the
identifier not be found within the ``--output-spaces`` list already.
In other words, running *ASLPrep* with ``--output-spaces MNI152NLin6Asym:res-2
--use-syn-sdc`` will expand the list of output spaces to 
``MNI152NLin6Asym:res-2 MNI152NLin2009cAsym``.
