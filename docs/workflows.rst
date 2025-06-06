.. include:: links.rst

###########################
Processing pipeline details
###########################

*ASLPrep* :footcite:p:`aslprep_nature_methods,aslprep_zenodo`
adapts its pipeline depending on what data and metadata are available and are used as inputs.
It requires the input data to be BIDS-valid and include necessary ASL parameters.

.. workflow::
   :graph2use: orig
   :simple_form: yes

   from aslprep.workflows.base import init_single_subject_wf

   wf = init_single_subject_wf("01")


************************
Structural Preprocessing
************************

The anatomical sub-workflow is from `sMRIPrep <https://github.com/nipreps/smriprep>`_.
It first constructs an average image by conforming all found T1w images
to a common voxel size, and, in the case of multiple images,
averages them into a single reference template.

.. workflow::
   :graph2use: orig
   :simple_form: yes

   from niworkflows.utils.spaces import Reference, SpatialReferences
   from smriprep.workflows.anatomical import init_anat_fit_wf

   wf = init_anat_fit_wf(
      bids_root=".",
      output_dir=".",
      freesurfer=True,
      hires=True,
      longitudinal=False,
      msm_sulc=True,
      t1w=["sub-01/anat/sub-01_T1w.nii.gz"],
      t2w=["sub-01/anat/sub-01_T2w.nii.gz"],
      flair=["sub-01/anat/sub-01_acq-flair_T2w.nii.gz"],
      skull_strip_mode="force",
      skull_strip_template=Reference("MNI152NLin2009cAsym"),
      spaces=SpatialReferences([
         ("MNI152NLin6Asym", {}),
         ("T1w", {}),
      ]),
      skull_strip_fixed_seed=False,
      omp_nthreads=1,
      sloppy=False,
   )

See also *sMRIPrep*'s
:py:func:`~smriprep.workflows.anatomical.init_anat_preproc_wf`.


Brain extraction, brain tissue segmentation and spatial normalization
=====================================================================

Next, the T1w reference is skull-stripped using a Nipype implementation of
the ``antsBrainExtraction.sh`` tool (ANTs), which is an atlas-based
brain extraction workflow:

.. workflow::
   :graph2use: orig
   :simple_form: yes

   from niworkflows.func.util import init_enhance_and_skullstrip_bold_wf

   wf = init_enhance_and_skullstrip_bold_wf()

An example of brain extraction is shown below:

.. figure:: _static/brainextraction_t1.svg

   Brain extraction

Once the brain mask is computed, FSL ``fast`` is used for brain tissue segmentation.

.. figure:: _static/segmentation.svg

   Brain tissue segmentation

Finally, spatial normalization to standard spaces is performed using ANTs' ``antsRegistration``
in a multiscale, mutual-information based, nonlinear registration scheme.
See :doc:`spaces` for more information on how standard and nonstandard spaces can
be set to resample the preprocessed data onto the final output spaces.

.. figure:: _static/T1MNINormalization.svg

   Animation showing spatial normalization of T1w onto the ``MNI152NLin2009cAsym`` template.


*****************
ASL preprocessing
*****************

:py:func:`~aslprep.workflows.asl.base.init_asl_wf`

.. workflow::
   :graph2use: orig
   :simple_form: yes

   from aslprep.tests.tests import mock_config
   from aslprep import config
   from aslprep.workflows.asl.base import init_asl_wf

   with mock_config():
      asl_file = str(config.execution.bids_dir / "sub-01" / "perf"/ "sub-01_asl.nii.gz")
      wf = init_asl_wf(asl_file=asl_file)

Preprocessing of :abbr:`ASL (Arterial Spin Labelling)` files is
split into multiple sub-workflows described below.

.. topic:: Processing GE data

   ASLPrep can process data from any of the big three manufacturers (Siemens, Philips, GE),
   but the GE ASL product sequence is unique in that it typically only produces a single deltaM or
   CBF volume (optionally along with an M0 volume).

   For these short sequences, ASLPrep operates in basically the same manner as it does with other
   sequences, except that the M0 volume (when available) is used as the reference image instead of
   the deltaM or CBF volume, since GE deltaM volumes tend to have extreme background noise in a
   circle around the brain.
   Additionally, SCORE/SCRUB cannot be used with these short sequences,
   as the denoising method requires a long time series from which to identify outliers.


.. _asl_ref:

ASL reference image estimation
==============================

:py:func:`~aslprep.workflows.asl.util.init_asl_reference_wf`

.. workflow::
   :graph2use: orig
   :simple_form: yes

   from aslprep.workflows.asl.reference import init_raw_aslref_wf

   wf = init_raw_aslref_wf()

This workflow estimates a reference image for an :abbr:`ASL (Arterial Spin Labelling)` series.
The reference image is then used to calculate a brain mask for the
:abbr:`ASL (Arterial Spin Labelling)` signal using *NiWorkflow's*
:py:func:`~niworkflows.func.util.init_enhance_and_skullstrip_bold_wf`.
Subsequently, the reference image is fed to the :ref:`head-motion estimation
workflow <asl_hmc>` and the :ref:`registration workflow <asl_reg>` to map the
ASL series onto the T1w image of the same subject.

.. figure:: _static/brainextraction.svg

   Calculation of a brain mask from the ASL series.


.. _asl_hmc:

Head-motion estimation
======================

:py:func:`~aslprep.workflows.asl.hmc.init_asl_hmc_wf`

.. workflow::
   :graph2use: orig
   :simple_form: yes

   from aslprep.workflows.asl.hmc import init_asl_hmc_wf

   wf = init_asl_hmc_wf(
      name="asl_hmc_wf",
      mem_gb=1,
      omp_nthreads=1,
   )

Using the previously estimated reference scan,
FSL ``mcflirt`` is used to estimate head-motion.
As a result, one rigid-body transform with respect to
the reference image is written for each :abbr:`ASL (Arterial Spin Labelling)`
time-step.
Additionally, a list of 6-parameters (three rotations and
three translations) per time-step is written and fed to the
:ref:`confounds workflow <asl_confounds>`,
for a more accurate estimation of head-motion.


.. _asl_confounds:

Confounds estimation
====================

:py:func:`~aslprep.workflows.asl.confounds.init_asl_confounds_wf`

.. workflow::
   :graph2use: orig
   :simple_form: yes

   from aslprep.workflows.asl.confounds import init_asl_confounds_wf

   wf = init_asl_confounds_wf(
      n_volumes=50,
      mem_gb=1,
      freesurfer=True,
      name="confound_wf",
   )


Calculated confounds include framewise displacement, 6 motion parameters, and DVARS.


Susceptibility Distortion Correction (SDC)
==========================================

One of the major problems that affects ASL data is the spatial distortion caused by the
inhomogeneity of the field inside the scanner.
Please refer to the `SDCFlows documentation <www.nipreps.org/sdcflows>`_
for details on the available workflows.

.. figure:: _static/unwarping.svg

   Applying susceptibility-derived distortion correction, based on fieldmap estimation

See also *SDCFlows*' :py:func:`~sdcflows.workflows.apply.correction.init_unwarp_wf` and
:py:func:`~sdcflows.workflows.apply.registration.init_coeff2epi_wf`


.. _asl_preproc:

Preprocessed ASL in native space
================================

A new *preproc* :abbr:`ASL (Arterial Spin Labelling)` series is generated
from the original data in the original space.
All volumes in the :abbr:`ASL (Arterial Spin Labelling)` series are
resampled in their native space by concatenating the mappings found in previous correction workflows
(:abbr:`HMC (head-motion correction)` and
:abbr:`SDC (susceptibility-derived distortion correction)`, if executed)
for a one-shot interpolation process.
Interpolation uses a Lanczos kernel.

.. figure:: _static/sub-20589_ses-11245_task-rest_desc-carpetplot_asl.svg

   The preprocessed ASL with label and control.
   The signal plots above the carpet plot are framewise displacement (FD) and DVARS.


.. _cbf_preproc:

*******************************
CBF Computation in native space
*******************************

:py:func:`~aslprep.workflows.asl.cbf.init_cbf_wf`

.. workflow::
   :graph2use: orig
   :simple_form: yes

   import json
   from pathlib import Path

   from aslprep.data import load as load_data
   from aslprep.workflows.asl.cbf import init_cbf_wf

   bids_dir = load_data("../tests/data/ds000240").absolute()
   nii_file = bids_dir / "sub-01" / "perf"/ "sub-01_asl.nii.gz"
   metadata_file = bids_dir / "sub-01" / "perf"/ "sub-01_asl.json"
   with open(metadata_file) as f:
      metadata = json.load(f)

   wf = init_cbf_wf(
      name_source=str(nii_file),
      processing_target="control",
      scorescrub=False,
      basil=False,
      metadata=metadata,
      m0_scale=1,
      smooth_kernel=5,
   )

ASL data consist of multiple pairs of labeled and control images.
*ASLPrep* first checks for proton density-weighted volume(s) (M0 scans).
In the absence of M0 images or an M0 estimate provided in the metadata,
the average of control images is used as the reference image.

After :ref:`preprocessing <asl_preproc>`, the pairs of labeled and control images are subtracted:

.. math::
    \Delta{M} = M_{C} - M_{L}

.. topic:: 2D Acquisitions and Slice Timing

   If slice timing is available for the ASL data, then ASLPrep will shift post-labeling delay
   values on a slice-wise basis.

   .. figure:: _static/slice_timing_and_pld.svg

      Slice time-shifted post-labeling delay values plotted next to a corresponding delta-M volume.


Single-Delay ASL
================

The CBF computation of single-delay (post labeling delay) ASL data is done using a one-compartment
model :footcite:p:`alsop_recommended_2015`.

.. sidebar:: Notation

   :math:`\tau` : Labeling duration, in seconds.

   :math:`\lambda` : Brain-blood partition coefficient.
   Set to 0.9 g/mL :footcite:p:`alsop_recommended_2015`.

   :math:`\alpha` : Labeling efficiency.
   This may be directly encoded in the ASL file's sidecar file with the ``LabelingEfficiency`` field.
   If that field is not set, then :math:`\alpha` will be determined based on the ASL type
   (PASL, CASL, or PCASL) and the number of background suppression pulses.

   :math:`w` : Post-labeling delay (PLD) for PCASL data, in seconds.
   In BIDS, this is encoded with the ``PostLabelingDelay`` field.

   :math:`TI` : Inversion time for PASL data, in seconds.
   In BIDS, this is encoded with the ``PostLabelingDelay`` field.

   :math:`T_{1,blood}` : Relaxation time for arterial blood, in seconds.
   ASLPrep infers this automatically based on the magnetic field strength
   :footcite:p:`zhang2013vivo,alsop_recommended_2015`.

   :math:`T_{1,tissue}` : Relaxation time for gray matter,
   used in the two-compartment model for multi-delay data (in seconds).
   ASLPrep infers this automatically based on the magnetic field strength
   :footcite:p:`wright2008water`.

   :math:`M_{0}` : Fully relaxed, equilibrium tissue magnetization.

   :math:`\Delta{TI}` : Post-labeling delay minus bolus cutoff delay time, in seconds.
   Per :footcite:t:`alsop_recommended_2015`, this is QUIPSS II's equivalent to (P)CASL's :math:`w`.

   :math:`TI_{1}` : Bolus cutoff delay time, in seconds.
   For Q2TIPS, this is the *first* bolus cutoff.

   :math:`TI_{2}` : For Q2TIPS, this is the *last* bolus cutoff delay time, in seconds.
   The other methods do not have this variable.


(Pseudo-)Continuous ASL
-----------------------

For (P)CASL ([pseudo-]continuous ASL),
CBF is calculated using a general kinetic model :footcite:p:`buxton1998general`.

.. math::
   CBF = \frac{ 6000 \cdot \lambda \cdot \Delta{M} \cdot e ^ \frac{ w }{ T_{1,blood} } }
   {2 \cdot \alpha \cdot M_{0} \cdot T_{1,blood} \cdot (1 - e^{\frac{ - \tau }{ T_{1,blood} } }) }

:math:`\tau`, :math:`\lambda`, :math:`\alpha`, and :math:`w` are labeling duration,
brain-blood partition coefficient, labeling efficiency, and post-labeling delay (PLD), respectively.

In the absence of any of these parameters, standard values are used based on the scan type and
scanning parameters.

The element which differentiates single-delay PCASL's CBF calculation from the PASL equivalents is
:math:`T1_{blood} \cdot (1 - e^{\frac{ - \tau }{ T1_{blood} } })`.


Pulsed ASL
----------

Currently, ASLPrep does not support PASL data without a bolus cut-off technique applied.


QUIPSS Modification
^^^^^^^^^^^^^^^^^^^

For pulsed ASL (PASL) data with the QUIPSS bolus cut-off technique,
the formula from :footcite:t:`wong1998quantitative` is used.

.. math::
   CBF = \frac{ 6000 \cdot \lambda \cdot \Delta{ M } \cdot e ^ \frac{ TI }{ T1_{blood} } }
   { 2 \cdot \alpha \cdot M_{0} \cdot \Delta{TI} }

where :math:`\Delta{TI}` is the post-labeling delay (PLD) minus the bolus cutoff delay time.

Given that :math:`TI` is equivalent to :math:`w` in BIDS datasets
(i.e., as the ``PostLabelingDelay`` field),
the formula for QUIPSS is the same as PCASL,
except :math:`\Delta{TI}` replaces
:math:`T1_{blood} \cdot (1 - e^{\frac{ - \tau }{ T1_{blood} } })`.


QUIPSS II Modification
^^^^^^^^^^^^^^^^^^^^^^

For PASL data with the QUIPSS II bolus cut-off technique,
the formula from :footcite:t:`alsop_recommended_2015` is used.

.. math::
   CBF = \frac{ 6000 \cdot \lambda \cdot \Delta{M} \cdot e ^ \frac{ TI }{ T1_{blood} } }
   {2 \cdot \alpha \cdot M_{0} \cdot TI_{1} }

where :math:`TI_{1}` is the bolus cutoff delay time.

Note that the formula for QUIPSS II is the same as the one for QUIPSS,
except :math:`TI_{1}` replaces :math:`\Delta{TI}`.


Q2TIPS Modification
^^^^^^^^^^^^^^^^^^^

For PASL data with the Q2TIPS bolus cut-off technique,
the formula from the commercial Q2TIPS CBF calculation is used,
as described in :footcite:t:`noguchi2015technical`.

.. math::
   CBF = \frac{ 6000 \cdot \lambda \cdot \Delta{M} \cdot e ^ { \frac{TI_{2} }{ T1_{blood} } } }
   { 2 \cdot \alpha \cdot M_{0} \cdot TI_{1} }

where :math:`TI_{1}` is the first bolus cutoff delay time and
:math:`TI_{2}` is the last bolus cutoff delay time.

Note that the formula for Q2TIPS is the same as QUIPSS II,
except :math:`TI_{2}` replaces :math:`TI` in the numerator.


Multi-Delay ASL
===============

In multi-delay ASL, control-label pairs are acquired for multiple post-labeling delay values.
This type of acquisition requires more complicated models, but it also results in more accurate
CBF estimates.
Also, multi-delay ASL allows for the estimation of arterial transit time (ATT),
arterial bolus arrival time (aBAT), and arterial blood volume (aBV).

For multi-delay ASL, ASLPrep uses a two-compartment general kinetic model recommended by
:footcite:t:`woods2023recommendations`.
The two-compartment model predicts the :math:`\Delta{M}` signal as a sum of two components:
a tissue compartment and an intravascular compartment.
The tissue compartment reflects labeled blood water that has reached the capillary bed in the
tissue at a given post-labeling delay.
However, some labeled blood water will remain in arteries at that post-labeling delay.
This is the intravascular compartment, which, when unaccounted for, will result in an
overestimation of CBF at voxels near the arteries.


Pseudo-Continuous ASL
---------------------

The tissue compartment is modeled as:

.. math::
   :label: woods_eq2

   \Delta{M}_{tiss} = \begin{cases}
      0 & \text{if } 0 < LD + PLD < ATT \\
      \frac{2 \cdot \alpha \cdot \alpha_{BS} \cdot T_{1b} \cdot M_{0a} \cdot CBF \cdot e^{- \frac{ATT}{T_{1b}}} \cdot \left(1 - e^{- \frac{LD + PLD - ATT}{T_{1b}}}\right)}{6000} & \text{if } ATT < LD + PLD < ATT + LD \\
      \frac{2 \cdot \alpha \cdot \alpha_{BS} \cdot T_{1b} \cdot M_{0a} \cdot CBF \cdot e^{- \frac{PLD}{T_{1b}}} \cdot \left(1 - e^{- \frac{LD}{T_{1b}}}\right)}{6000} & \text{if } ATT < PLD \\
   \end{cases}

The intravascular compartment is modeled as:

.. math::
   :label: woods_eq3

   \Delta{M}_{art} = \begin{cases}
      0 & \text{if } 0 < LD + PLD < aBAT \\
      2 \cdot \alpha \cdot \alpha_{BS} \cdot M_{0b} \cdot aBV \cdot e^{- \frac{aBAT}{T_{1b}}} & \text{if } aBAT < LD + PLD < aBAT + LD \\
      0 & \text{if } aBAT < PLD \\
   \end{cases}

For reference, please see Equations 2 and 4 in :footcite:t:`woods2023recommendations`.


Pulsed ASL
----------

.. math::
   :label: woods_eq3

   \Delta{M}_{tiss} = \begin{cases}
      0 & \text{if } 0 < TI < ATT \\
      \frac{2 \cdot \alpha \cdot \alpha_{BS} \cdot M_{0a} \cdot CBF \cdot e^{- \frac{TI}{T_{1b}}} \cdot (TI - ATT)}{6000} & \text{if } ATT < TI < ATT + TI_1 \\
      \frac{2 \cdot \alpha \cdot \alpha_{BS} \cdot M_{0a} \cdot CBF \cdot e^{- \frac{TI}{T_{1b}}} \cdot TI_1}{6000} & \text{if } ATT + TI_1 < TI \\
   \end{cases}

The intravascular compartment is modeled as:

.. math::
   :label: woods_eq3

   \Delta{M}_{art} = \begin{cases}
      0 & \text{if } 0 < TI < aBAT \\
      2 \cdot \alpha \cdot \alpha_{BS} \cdot M_{0b} \cdot aBV \cdot e^{- \frac{TI}{T_{1b}}} & \text{if } aBAT < TI < aBAT + TI_1 \\
      0 & \text{if } aBAT + TI_1 < TI \\
   \end{cases}

For reference, please see Equations 3 and 5 in :footcite:t:`woods2023recommendations`.


Additional Denoising Options
============================

For cases where data may be especially noisy (e.g., due to motion or a low-SNR protocol),
ASLPrep includes options to additionally denoise CBF estimates.

The two current options are SCORE/SCRUB and BASIL.

SCORE and SCRUB
---------------

ASLPrep includes the ability to denoise CBF with SCORE and SCRUB.

Structural Correlation based Outlier Rejection (SCORE) :footcite:p:`dolui2017structural`
detects and discards extreme outliers in the CBF volume(s) from the CBF time series.
SCORE first discards CBF volumes whose CBF within grey matter (GM)
means are 2.5 standard deviations away from the median of the CBF within GM.
Next, it iteratively removes volumes that are most structurally correlated
to the intermediate mean CBF map unless the variance within each tissue type starts increasing
(which implies an effect of white noise removal as opposed to outlier rejection).

The mean CBF after denoising by SCORE is plotted below

.. figure:: _static/sub-20589_ses-11245_task-rest_desc-scoreplot_asl.svg

   Computed CBF maps denoised by SCORE

After discarding extreme outlier CBF volume(s) (if present) by SCORE,
SCRUB (Structural Correlation with RobUst Bayesian) uses robust Bayesian estimation
of CBF using iterative reweighted least square method :footcite:p:`dolui2016scrub` to denoise CBF.
The SCRUB algorithm is described below:

.. math::
   CBF_{SCRUB} =  \arg\max_{\theta} \sum_{t=1}^N \rho(CBF_{t} -\theta)  + \lambda(\theta -\mu)^2

   \mu =\sum_{i \in Tissue type} p \cdot \mu_{i}

:math:`CBF_{t}`, :math:`\mu`, :math:`\theta`, and :math:`p` equal CBF time series
(after any extreme outliers are discarded by SCORE),
mean CBF, ratio of temporal variance at each voxel to overall variance of all voxels,
and probability tissue maps, respectively.
Other variables include :math:`\lambda` and :math:`\rho` that represent the weighting parameter
and Tukey's bisquare function, respectively.

An example of CBF denoised by SCRUB is shown below.

.. figure:: _static/sub-20589_ses-11245_task-rest_desc-scrubplot_asl.svg

   Computed CBF maps denoised by SCRUB

BASIL
-----

*ASLPrep* also includes the option to compute CBF using BASIL.

Bayesian Inference for Arterial Spin Labeling (`BASIL <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BASIL>`_)
is an FSL tool for CBF estimation.

BASIL implements a simple kinetic model as described above,
but uses Bayesian inference principles :footcite:p:`chappell2008variational`.
BASIL is mostly suitable for multi-delay ASL data.
It includes bolus arrival time estimation with spatial regularization :footcite:p:`groves2009combined`
and the correction of partial volume effects :footcite:p:`chappell2011partial`.

A sample of BASIL CBF with spatial regularization is shown below:

.. figure:: _static/sub-20589_ses-11245_task-rest_desc-basilplot_asl.svg

   Computed CBF maps by BASIL

The CBF map shown below is the result of partial volume corrected CBF computed by BASIL.

.. figure:: _static/sub-20589_ses-11245_task-rest_desc-pvcplot_asl.svg

   Partial volume corrected CBF maps by BASIL


************************
Quality control measures
************************

:py:func:`~aslprep.workflows.asl.confounds.init_cbf_confounds_wf`

.. workflow::
   :graph2use: orig
   :simple_form: yes

   from pathlib import Path

   from aslprep.data import load as load_data
   from aslprep.workflows.asl.qc import init_cbf_confounds_wf

   bids_dir = load_data("../tests/data/ds000240").absolute()
   asl_file = bids_dir / "sub-01" / "perf"/ "sub-01_asl.nii.gz"
   metadata = bids_dir / "sub-01" / "perf"/ "sub-01_asl.json"

   wf = init_cbf_confounds_wf()

Quality control (QC) measures such as FD (framewise displacement), coregistration, normalization index,
and quality evaluation index (QEI) are included for all CBF maps.
The QEI :footcite:p:`dolui2024automated` evaluates the quality of the computed CBF maps considering
three factors:
structural similarity, spatial variability, and percentage of voxels in GM with negative CBF.


.. _asl_reg:

*******************************
ASL and CBF to T1w registration
*******************************

:py:func:`~fmriprep.workflows.bold.registration.init_bold_reg_wf`

.. workflow::
   :graph2use: orig
   :simple_form: yes

   from fmriprep.workflows.bold.registration import init_bold_reg_wf

   wf = init_bold_reg_wf(
      freesurfer=True,
      use_bbr=True,
      bold2anat_dof=6,
      bold2anat_init="auto",
      mem_gb=1,
      omp_nthreads=1,
      name="bold_reg_wf",
      sloppy=False,
   )

*ASLPrep* uses the ``FSL BBR`` routine to calculate the alignment between each run's
:abbr:`ASL (arterial spin labelling)` reference image and the reconstructed subject using the
gray/white matter boundary.


.. figure:: _static/EPIT1Normalization.svg

   Animation showing :abbr:`ASL (arterial spin labelling)` to T1w registration.


FSL ``flirt`` is run with the :abbr:`BBR (boundary-based registration)` cost function, using the
``fast`` segmentation to establish the gray/white matter boundary.
After :abbr:`BBR (boundary-based registration)` is run,
the resulting affine transform will be compared to the initial transform found by ``flirt``.
Excessive deviation will result in rejection of the BBR refinement and acceptance
of the original affine registration.
The computed :ref:`CBF <cbf_preproc>` is registered to T1w using the transformation from ASL-T1w
registration.


Resampling ASL and CBF runs onto standard spaces
================================================

:py:func:`~fmriprep.workflows.bold.apply.init_bold_volumetric_resample_wf`

.. workflow::
   :graph2use: orig
   :simple_form: yes

   from fmriprep.workflows.bold.apply import init_bold_volumetric_resample_wf

   wf = init_bold_volumetric_resample_wf(
      metadata={
         'RepetitionTime': 2.0,
         'PhaseEncodingDirection': 'j-',
         'TotalReadoutTime': 0.03,
      },
      fieldmap_id='my_fieldmap',
   )

This sub-workflow concatenates the transforms calculated upstream
(see `Head-motion estimation`_, `Susceptibility Distortion Correction (SDC)`_)
if fieldmaps are available, and an anatomical-to-standard transform from
`Structural Preprocessing`_ to map the ASL and CBF images to the standard spaces is given by the
``--output-spaces`` argument (see :doc:`spaces`).
It also maps the T1w-based mask to each of those standard spaces.

Transforms are concatenated and applied all at once, with one interpolation (Lanczos) step,
so as little information is lost as possible.


References
==========

.. footbibliography::
