*ASLPrep*: A Robust Preprocessing Pipeline for ASL Data
=========================================================

This pipeline is developed by the `Sattarwaite lab at the University of Pennysilvania
<https://www.satterthwaitelab.com/>`_ for use at the `The Lifespan Informatics and Neuroimaging Center 
at the University of Pennylvannia <https://www.satterthwaitelab.com/>`_, as well as for
open-source software distribution.

About
-----


.. image:: https://raw.githubusercontent.com/a3sha2/aslprep/4b24ab000d6736a99874b80e04f23fe9a64a0eba/docs/_static/aslprepworkflow.png


*ASLPrep* is a Arterial Spin Labeling  (ASL) data
preprocessing  and Cerebral Blood FLow (CBF) computation pipeline 
that is designed to provide an easily accessible,
state-of-the-art interface that is robust to variations in scan acquisition
protocols and that requires minimal user input, while providing easily
interpretable and comprehensive error and output reporting.
It performs basic processing steps (coregistration, normalization, unwarping,
noise component extraction, segmentation, skullstripping etc.) providing
outputs that can be easily submitted to a variety of group level analyses,
including task-based or resting-state CBF, graph theory measures, surface or
volume-based statistics, etc.


The *ASLPrep* pipeline uses a combination of tools from well-known software
packages, including FSL_, ANTs_, FreeSurfer_ and AFNI_.
This pipeline was designed to provide the best software implementation for each
state of preprocessing, and will be updated as newer and better neuroimaging
software become available.

This tool allows you to easily do the following:

- Take ASL data from raw to fully preprocessed form.
- Implement tools from different software packages.
- Achieve optimal data processing quality by using the best tools available.
- Receive verbose output concerning the stage of preprocessing for each
  subject, including meaningful errors.
- Automate and parallelize processing steps, which provides a significant
  speed-up from typical linear, manual processing.
- compute Cerebral Blood Flow(CBF), denoising and partial volume correction

More information and documentation can be found at
https://aslprep.readthedocs.io/

Principles
----------

*ASLPrep* is built around three principles:

1. **Robustness** - The pipeline adapts the preprocessing steps depending on
   the input dataset and should provide results as good as possible
   independently of scanner make, scanning parameters or presence of additional
   correction scans (such as fieldmaps).
2. **Ease of use** - Thanks to dependence on the BIDS standard, manual
   parameter input is reduced to a minimum, allowing the pipeline to run in an
   automatic fashion.
3. **"Glass box"** philosophy - Automation should not mean that one should not
   visually inspect the results or understand the methods.
   Thus, *ASLPrep* provides visual reports for each subject, detailing the
   accuracy of the most important processing steps.
   This, combined with the documentation, can help researchers to understand
   the process and decide which subjects should be kept for the group level
   analysis.


Limitations and reasons not to use *ASLPrep*
---------------------------------------------

1. Very narrow :abbr:`FoV (field-of-view)` images oftentimes do not contain
   enough information for standard image registration methods to work correctly.
   Also, problems may arise when extracting the brain from these data.
   Supporting these particular images is already a future line of the development
   road-map.
2. *ASLPrep* may also underperform for particular populations (e.g., infants) and
   non-human brains, although appropriate templates can be provided to overcome the
   issue.
3. If you really want unlimited flexibility (which is obviously a double-edged sword).
4. If you want students to suffer through implementing each step for didactic purposes,
   or to learn shell-scripting or Python along the way.
5. If you are trying to reproduce some *in-house* lab pipeline.


Acknowledgements
----------------

Please acknowledge this work using the citation boilerplate that *ASLPrep* includes
in the visual report generated for every subject processed.
