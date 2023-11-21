#######################################################
*ASLPrep*: A Robust Preprocessing Pipeline for ASL Data
#######################################################

.. image:: https://img.shields.io/badge/Source%20Code-pennlinc%2Faslprep-purple
   :target: https://github.com/PennLINC/aslprep
   :alt: GitHub Repository

.. image:: https://readthedocs.org/projects/aslprep/badge/?version=latest
   :target: http://aslprep.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation

.. image:: https://img.shields.io/badge/docker-pennlinc/aslprep-brightgreen.svg?logo=docker&style=flat
   :target: https://hub.docker.com/r/pennlinc/aslprep/tags/
   :alt: Docker

.. image:: https://circleci.com/gh/PennLINC/aslprep.svg?style=svg
   :target: https://circleci.com/gh/PennLINC/aslprep
   :alt: Test Status

.. image:: https://codecov.io/gh/PennLINC/aslprep/branch/main/graph/badge.svg
   :target: https://app.codecov.io/gh/PennLINC/aslprep/tree/main
   :alt: Codecov

.. image:: https://img.shields.io/badge/Nature%20Methods-10.1038%2Fs41592--022--01458--7-purple
   :target: https://doi.org/10.1038/s41592-022-01458-7
   :alt: Publication DOI

.. image:: https://zenodo.org/badge/256420694.svg
   :target: https://zenodo.org/badge/latestdoi/256420694
   :alt: Zenodo DOI

.. image:: https://img.shields.io/badge/License-BSD--3--Clause-green
   :target: https://opensource.org/licenses/BSD-3-Clause
   :alt: License

This pipeline is developed by the `Satterthwaite lab at the University of Pennsylvania
<https://www.satterthwaitelab.com/>`_ for use at the `The Lifespan Informatics and Neuroimaging Center
at the University of Pennsylvania <https://www.satterthwaitelab.com/>`_, as well as for
open-source software distribution.

*****
About
*****

.. image:: https://raw.githubusercontent.com/PennLINC/aslprep/main/docs/_static/aslprepworkflow.png

*ASLPrep* is a Arterial Spin Labeling  (ASL) data
preprocessing  and Cerebral Blood Flow (CBF) computation pipeline
that is designed to provide an easily accessible,
state-of-the-art interface that is robust to variations in scan acquisition
protocols and that requires minimal user input, while providing easily
interpretable and comprehensive error and output reporting.
It performs basic processing steps (coregistration, normalization, unwarping,
noise component extraction, segmentation, skullstripping etc.),
CBF computation, denoising CBF, CBF partial volume correction,
and providing outputs that can be easily submitted to a variety of group level analyses,
including task-based or resting-state CBF, graph theory measures, surface or volume-based statistics, etc.

The *ASLPrep* pipeline uses a combination of tools from well-known software
packages, including FSL_, ANTs_, FreeSurfer_ and AFNI_.
This pipeline was designed to provide the best software implementation for each state of preprocessing,
and will be updated as newer and better neuroimaging software become available.

This tool allows you to easily do the following:

- Take ASL data from raw to fully preprocessed form.
- Compute Cerebral Blood Flow (CBF), denoising and partial volume correction.
- Implement tools from different software packages.
- Achieve optimal data processing quality by using the best tools available.
- Receive verbose output concerning the stage of preprocessing for each
  subject, including meaningful errors.
- Automate and parallelize processing steps, which provides a significant
  speed-up from typical linear, manual processing.

More information and documentation can be found at https://aslprep.readthedocs.io/.

*******
ASLPrep
*******

*ASLPrep* adapts the preprocessing steps depending on the input dataset
and provide results as good as possible independently of scanner make and scanning parameters
With the BIDS input, little or no parameter are required allowing ease of operation.
*ASLPrep* also provides visual reports for each subject,
detailing the most important processing steps.

****************
Acknowledgements
****************

Please acknowledge this work using the citation boilerplate that *ASLPrep* includes
in the visual report generated for every subject processed.

************************************************
On the relationship between ASLPrep and fMRIPrep
************************************************

ASLPrep is largely based on fMRIPrep, as ASL processing is very similar to fMRI processing.
fMRIPrep is developed by a larger group and receives more regular feedback on its workflow,
so we (the ASLPrep developers) try to update ASLPrep in line with fMRIPrep.

ASLPrep and fMRIPrep are both part of the NiPreps community.

There are several crucial differences between ASL and fMRI processing,
which must be accounted for in ASLPrep:

1. ASL processing does not include slice timing correction.
   Instead, post-labeling delays are shifted based on the slice timing when CBF is calculated.
2. While ASL motion correction may use the same algorithm as fMRI motion correction,
   different volume types in the ASL time series exhibit different contrasts, which can introduce
   artifacts into the motion-corrected data.
   As such, ASLPrep motion corrects each volume type separately,
   and then concatenated the corrected time series back together.
3. ASLPrep includes extra steps to do the following:
   1. Calculate CBF.
   2. Calculate CBF QC metrics, paralleling the ASL confound calculation.
   3. Plot CBF results, paralleling ASL plots.
   4. Parcellating CBF results with a range of atlases.
4. fMRIPrep contains a lot of code to handle multi-echo fMRI.
   While multi-echo ASl does exist, it is very rare, so we do not include any multi-echo-specific
   elements in ASLPrep.
