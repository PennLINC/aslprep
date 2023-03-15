.. include:: links.rst

###############################
ASL processing pipeline details
###############################

*ASLPrep* :footcite:p:`aslprep_nature_methods,aslprep_zenodo`
adapts its pipeline depending on what data and metadata are
available and are used as inputs.
It requires the input data to be BIDS-valid and include necessary ASL parameters.

.. workflow::
    :graph2use: orig
    :simple_form: yes

    from aslprep.workflows.base import init_single_subject_wf

    wf = init_single_subject_wf('01')


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

    from smriprep.workflows.anatomical import init_anat_preproc_wf

    from aslprep.utils.niworkflows import Reference, SpatialReferences

    wf = init_anat_preproc_wf(
        bids_root='.',
        freesurfer=False,
        hires=True,
        longitudinal=False,
        omp_nthreads=1,
        output_dir='.',
        skull_strip_mode='force',
        skull_strip_template=Reference('MNI152NLin2009cAsym'),
        spaces=SpatialReferences([
            ('MNI152Lin', {}),

            ('T1w', {}),
            ('fsnative', {})
        ]),
        skull_strip_fixed_seed=False,
        t1w=['sub-01/anat/sub-01_T1w.nii.gz'],
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

    from niworkflows.anat.ants import init_brain_extraction_wf
    wf = init_brain_extraction_wf()


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

    Animation showing spatial normalization of T1w onto the ``MNI152NLin2009cAsym`` template


*****************
ASL preprocessing
*****************

:py:func:`~aslprep.workflows.asl.base.init_asl_preproc_wf`

.. workflow::
    :graph2use: orig
    :simple_form: yes

    from aslprep.tests.tests import mock_config
    from aslprep import config
    from aslprep.workflows.asl.base import init_asl_preproc_wf
    with mock_config():
        asl_file = config.execution.bids_dir / 'sub-01' / 'perf'/ 'sub-01_asl.nii.gz'
        wf = init_asl_preproc_wf(str(asl_file))

Preprocessing of :abbr:`ASL (Arterial Spin Labelling)` files is
split into multiple sub-workflows described below.


.. _asl_ref:

ASL reference image estimation
==============================

:py:func:`~aslrep.niworkflows.func.util.init_asl_reference_wf`

.. workflow::
    :graph2use: orig
    :simple_form: yes

    from niworkflows.func.util import init_bold_reference_wf
    wf = init_bold_reference_wf(omp_nthreads=1, brainmask_thresh=0.1)

This workflow estimates a reference image for an
:abbr:`ASL (Arterial Spin Labelling)` series.
The reference image is then used to calculate a brain mask for the
:abbr:`ASL (Arterial Spin Labelling)` signal using *NiWorkflow's*
:py:func:`~niworkflows.func.util.init_enhance_and_skullstrip_asl_wf`.
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
        mem_gb=1,
        omp_nthreads=1,
    )

Using the previously :ref:`estimated reference scan <asl_ref>`,
FSL ``mcflirt`` or AFNI ``3dvolreg`` is used to estimate head-motion.
As a result, one rigid-body transform with respect to
the reference image is written for each :abbr:`ASL (Arterial Spin Labelling)`
time-step.
Additionally, a list of 6-parameters (three rotations and
three translations) per time-step is written and fed to the
:ref:`confounds workflow <asl_confounds>`,
for a more accurate estimation of head-motion.


Slice time correction
=====================

:py:func:`~aslprep.workflows.asl.stc.init_asl_stc_wf`

.. workflow::
    :graph2use: orig
    :simple_form: yes

    from aslprep.workflows.asl.stc import init_asl_stc_wf

    wf = init_asl_stc_wf(
        metadata={
            'RepetitionTime': 2.0,
            'SliceTiming': [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
        },
    )

If the ``SliceTiming`` field is available within the input dataset metadata,
this workflow performs slice time correction prior to other signal resampling
processes.
Slice time correction is performed using AFNI ``3dTShift``.
All slices are realigned in time to the middle of each TR.

Slice time correction can be disabled with the ``--ignore slicetiming``
command line argument.


.. _asl_confounds:

Confounds estimation
====================

:py:func:`~aslprep.workflows.asl.confounds.init_asl_confs_wf`

.. workflow::
    :graph2use: orig
    :simple_form: yes

    from aslprep.workflows.asl.confounds import init_asl_confs_wf

    wf = init_asl_confs_wf(name="confound_wf", mem_gb=1)


Calculated confounds include framewise displacement, 6 motion parameters, and DVARS.


Susceptibility Distortion Correction (SDC)
==========================================

One of the major problems that affects :abbr:`EPI (echo planar imaging)` data
is the spatial distortion caused by the inhomogeneity of the field inside
the scanner.
Please refer to :doc:`sdc` for details on the available workflows.

.. figure:: _static/unwarping.svg

    Applying susceptibility-derived distortion correction, based on fieldmap estimation

See also *SDCFlows*' :py:func:`~sdcflows.workflows.base.init_sdc_estimate_wf`


.. _asl_preproc:

Preprocessed ASL in native space
================================

:py:func:`~aslprep.workflows.asl.resampling.init_asl_preproc_trans_wf`

.. workflow::
    :graph2use: orig
    :simple_form: yes

    from aslprep.workflows.asl.resampling import init_asl_preproc_trans_wf

    wf = init_asl_preproc_trans_wf(
        mem_gb=3,
        omp_nthreads=1,
    )

A new *preproc* :abbr:`ASL (Arterial Spin Labelling)` series is generated
from either the slice-timing corrected data or the original data (if
:abbr:`STC (slice-timing correction)` was not applied) in the
original space.
All volumes in the :abbr:`ASL (Arterial Spin Labelling)` series are
resampled in their native space by concatenating the mappings found in previous
correction workflows (:abbr:`HMC (head-motion correction)` and
:abbr:`SDC (susceptibility-derived distortion correction)`, if executed)
for a one-shot interpolation process.
Interpolation uses a Lanczos kernel.

.. figure:: _static/sub-20589_ses-11245_task-rest_desc-carpetplot_asl.svg

    The preprocessed ASL with label and control.
    The signal plots above the carpet plot are framewise diplacement (FD) and DVRAS.


.. _cbf_preproc:

*******************************
CBF Computation in native space
*******************************

:py:func:`~aslprep.workflows.asl.cbf.init_cbf_compt_wf`

.. workflow::
    :graph2use: orig
    :simple_form: yes

    import json
    from pathlib import Path
    from pkg_resources import resource_filename as pkgrf

    from aslprep.workflows.asl.cbf import init_cbf_compt_wf

    bids_dir=Path(pkgrf('aslprep', 'tests/data/ds000240')).absolute()
    metadata_file = bids_dir / 'sub-01' / 'perf'/ 'sub-01_asl.json'
    with open(metadata_file) as f:
        metadata = json.load(f)

    wf = init_cbf_compt_wf(
        bids_dir=str(bids_dir),
        scorescrub=False,
        basil=False,
        metadata=metadata,
        M0Scale=1,
        smooth_kernel=5,
        dummy_vols=0,
    )

ASL data consist of multiple pairs of labeled and control images. *ASLPrep* first checks for
the reference ASL volume(s) (M0,``sub-task_xxxx-acq-YYY_m0scan.nii.gz``). In the absence of M0 images,
the average of control images is used as the reference image.

After :ref:`preprocessing <asl_preproc>`, the pairs of labeled and control images are
subtracted:

.. math::
   ASL_{signal}= M_{C} - M_{L}

The CBF computation of either single or multiple PLD (post labelling delay)
is done using a relatively simple model. For P/CASL (pseudo continuous ASL), CBF
is calculated by using a general kinetic model :footcite:p:`buxton1998general`:

.. math::
   CBF = \frac{ 6000 * \lambda * (M_{C} - M_{L})* e ^ {PLD/T1_{blood}}} {2 * \alpha * T1_{blood}  * M_{0} * (1 - e^{- \tau / T1_{blood} }) }


PASL (Pulsed ASL) is also computed by the QUIPSS model :footcite:p:`wong1998quantitative`:


.. math::
   CBF = \frac{ 6000 * \lambda * (M_{C} - M_{L})* e ^ {PLD/T1_{blood}}} {2 * \alpha * TI  * M_{0}}


:math:`\tau,\quad \lambda,\quad and\quad \alpha\quad` are label duration, brain-blood partition
coefficient, and labeling efficiency, respectively.
In the absence of any of these parameters, standard values are used based on the scan type and
scanning parameters.

The computed CBF time series is shown in carpet plot below.

.. figure:: _static/sub-20589_ses-11245_task-rest_desc-cbftsplot_asl.svg

    The carpet plot of computed CBF. The step plot above indicated the volume(s) marked by SCORE
    algorithm to be contaminated by noise.


Mean CBF is computed from the average of CBF timeseries.

.. figure:: _static/sub-20589_ses-11245_task-rest_desc-cbfplot_asl.svg

   Computed CBF maps

For multi-PLDs (Post Labeling Delay) ASL data,
the CBF is first computed for each PLD and the weighted average CBF is computed
over all PLDs at time = t :footcite:p:`dai2012reduced`.

.. math::
   CBF_{t} =\frac {\sum_{i}^{NPLDs} PLD_{i} * CBF_{i}} { \sum_{i}^{NPLDs} PLD_{i} }

ASLPrep includes option of CBF denoising by SCORE and SCRUB.
Structural Correlation based Outlier Rejection (SCORE) :footcite:p:`dolui2017structural`
detects and discards extreme outliers in the CBF volume(s) from the CBF time series.
SCORE first discards CBF volumes whose CBF within grey matter (GM)
means are 2.5 standard deviations away from the median of the CBF within GM.
Next, it iteratively removes volumes that are most structurally correlated
to the intermediate mean CBF map unless the variance within each
tissue type starts increasing
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

   \mu =\sum_{i \in Tissue type} p *\mu_{i}

:math:`CBF_{t},\quad  \mu,\quad  \theta,\quad  and\quad  p` equal CBF time series
(after any extreme outliers are discarded by SCORE),
mean CBF, ratio of temporal variance at each voxel to overall variance of all voxels,
and probability tissue maps, respectively.
Other variables include :math:`\lambda\quad and\quad \rho` that represent the weighting parameter
and Tukey's bisquare function, respectively.

An example of CBF denoised by SCRUB is shown below.

.. figure:: _static/sub-20589_ses-11245_task-rest_desc-scrubplot_asl.svg

   Computed CBF maps denoised by SCRUB

*ASLPrep* also includes option of CBF computation by Bayesian Inference for Arterial Spin Labeling
`(BASIL) <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BASIL>`_.
BASIL also implements a simple kinetic model as described above,
but using Bayesian Inference principles :footcite:p:`chappell2008variational`.
BASIL is mostly suitable for multi-PLD.
It includes bolus arrival time estimation with spatial regularization :footcite:p:`groves2009combined`
and the correction of partial volume effects :footcite:p:`chappell2011partial`.

The sample of BASIL CBF with spatial regularization is shown below:

.. figure:: _static/sub-20589_ses-11245_task-rest_desc-basilplot_asl.svg

   Computed CBF maps by BASIL

The CBF map shown below is the result of partial volume corrected CBF computed by BASIL.

.. figure:: _static/sub-20589_ses-11245_task-rest_desc-pvcplot_asl.svg

   Partial volume corrected CBF maps by BASIL


************************
Quality control measures
************************

:py:func:`~aslprep.workflows.asl.cbf.init_cbfqc_compt_wf`

.. workflow::
    :graph2use: orig
    :simple_form: yes

    from pathlib import Path
    from pkg_resources import resource_filename as pkgrf

    from aslprep.workflows.asl.cbf import init_cbfqc_compt_wf

    bids_dir = Path(pkgrf('aslprep', 'tests/data/ds000240')).absolute()
    asl_file = bids_dir / 'sub-01' / 'perf'/ 'sub-01_asl.nii.gz'
    metadata = bids_dir / 'sub-01' / 'perf'/ 'sub-01_asl.json'

    wf = init_cbfqc_compt_wf(asl_file=str(asl_file))

Quality control (QC) measures such as FD (framewise displacement), coregistration, normalization index,
and quality evaluation index (QEI) are included for all CBF maps.
The QEI :footcite:p:`dolui2017automated` evaluates the quality of the computed CBF maps considering
three factors:
structural similarity, spatial variability, and percentage of voxels in GM with negative CBF.


.. _asl_reg:

*******************************
ASL and CBF to T1w registration
*******************************

:py:func:`~aslprep.workflows.asl.registration.init_asl_reg_wf`

.. workflow::
    :graph2use: orig
    :simple_form: yes

    from aslprep.workflows.asl.registration import init_asl_reg_wf

    wf = init_asl_reg_wf(
        use_bbr=True,
        asl2t1w_dof=6,
        asl2t1w_init='register',
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

:py:func:`~aslprep.workflows.asl.resampling.init_asl_std_trans_wf`

.. workflow::
    :graph2use: orig
    :simple_form: yes

    from aslprep.utils.niworkflows import SpatialReferences
    from aslprep.workflows.asl.resampling import init_asl_std_trans_wf
    wf = init_asl_std_trans_wf(
        mem_gb=3,
        omp_nthreads=1,
        spaces=SpatialReferences(
            spaces=[('MNI152Lin', {})],
            checkpoint=True,
        ),
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
