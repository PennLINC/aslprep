*ASLPrep*: A Robust Preprocessing Pipeline for ASL Data
=========================================================

This pipeline is developed by the `Satterthwaite lab at the University of Pennysilvania
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
noise component extraction, segmentation, skullstripping etc.), CBF computation,
denoising CBF, CBF partial volume correction and providing
outputs that can be easily submitted to a variety of group level analyses,
including task-based or resting-state CBF, graph theory measures, surface or
volume-based statistics, etc.


The *ASLPrep* pipeline uses a combination of tools from well-known software
packages, including FSL_, ANTs_, FreeSurfer_ and AFNI_ .
This pipeline was designed to provide the best software implementation for each
state of preprocessing, and will be updated as newer and better neuroimaging
software become available.

This tool allows you to easily do the following:

- Take ASL data from raw to fully preprocessed form.
- Compute Cerebral Blood Flow(CBF), denoising and partial volume correction
- Implement tools from different software packages.
- Achieve optimal data processing quality by using the best tools available.
- Receive verbose output concerning the stage of preprocessing for each
  subject, including meaningful errors.
- Automate and parallelize processing steps, which provides a significant
  speed-up from typical linear, manual processing.
- 

More information and documentation can be found at
https://aslprep.readthedocs.io/

ASLPrep
--------

*ASLPrep* adapts the preprocessing steps depending on the input dataset 
and provide results as good as possible independently of scanner make and scanning parameters 
With the BIDS input, little or no parameter are required allowing ease of operation. 
*ASLPrep*  also provides visual reports for each subject,
detailing the the most important processing steps.



Acknowledgements
----------------

Please acknowledge this work using the citation boilerplate that *ASLPrep* includes
in the visual report generated for every subject processed. `(link) <https://aslprep.readthedocs.io/en/latest/citing.html>`_
