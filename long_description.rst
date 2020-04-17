Preprocessing of arterial spin labeling (ASL)  involves numerous steps to clean and standardize
the data before statistical analysis.
Generally, researchers create ad hoc preprocessing workflows for each dataset,
building upon a large inventory of available tools.
The complexity of these workflows has snowballed with rapid advances in
acquisition and processing.
ASLPrep is an analysis-agnostic tool that addresses the challenge of robust and
reproducible preprocessing for task-based and resting ASL data.
ASLPrep automatically adapts a best-in-breed workflow to the idiosyncrasies of
virtually any dataset, ensuring high-quality preprocessing without manual intervention.
ASLPrep robustly produces high-quality results on diverse ASL data.
Additionally, ASLPrep introduces less uncontrolled spatial smoothness than observed
with commonly used preprocessing tools.
ASLPrep equips neuroscientists with an easy-to-use and transparent preprocessing
workflow, which can help ensure the validity of inference and the interpretability
of results.

The workflow is based on `Nipype <https://nipype.readthedocs.io>`_ and encompases a large
set of tools from well-known neuroimaging packages, including
`FSL <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/>`_,
`ANTs <https://stnava.github.io/ANTs/>`_,
`FreeSurfer <https://surfer.nmr.mgh.harvard.edu/>`_,
`AFNI <https://afni.nimh.nih.gov/>`_,
and `Nilearn <https://nilearn.github.io/>`_.
This pipeline was designed to provide the best software implementation for each state of
preprocessing, and will be updated as newer and better neuroimaging software becomes
available.

ASLPrep performs basic preprocessing steps (coregistration, normalization, unwarping,segmentation, 
skullstripping  and computation of  cerebral blood flow (CBF)) providing outputs that can be
easily submitted to a variety of group level analyses, including task-based or resting-state
CBF, graph theory measures, surface or volume-based statistics, etc.
ASLPrep allows you to easily do the following:

  * Take ASL data from *unprocessed* (only reconstructed) to ready for analysis.
  * Implement tools from different software packages.
  * Achieve optimal data processing quality by using the best tools available.
  * Generate preprocessing-assessment reports, with which the user can easily identify problems.
  * Receive verbose output concerning the stage of preprocessing for each subject, including
    meaningful errors.
  * Automate and parallelize processing steps, which provides a significant speed-up from
    typical linear, manual processing.

[Documentation `aslprep.org <https://aslprep.readthedocs.io>`_]


