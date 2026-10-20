.. include:: links.rst

####################
Outputs of *ASLPrep*
####################

*ASLPrep* generates three broad classes of outputs:

1. **Visual QA (quality assessment) reports**:
   one :abbr:`HTML (hypertext markup language)` per subject,
   that allows the user to conduct a thorough visual assessment of the raw and
   processed data.
   This also includes quality control measures.

2. **Derivatives (preprocessed data and computed CBF):** the input ASL data ready for
   analysis, (i.e., after the various preparation procedures have been applied),
   the brain mask, and ASL images after masking has been applied.
   Other outputs include computed CBF maps and post-processed data such as denoised and partial
   volume-corrected CBF.

3. **Confounds and quality control metrics**: a confounds matrix that includes framewise
   displacement, motion parameters, coregistration and registration quality indices,
   and CBF quality control metrics.


******
Layout
******

Assuming ASLPrep is invoked with::

   aslprep <input_dir>/ <output_dir>/ participant [OPTIONS]

The outputs will be a `BIDS Derivatives`_ dataset of the form::

   <output_dir>/
      logs/
      atlases/
      sub-<label>/
      sub-<label>.html
      dataset_description.json
      .bidsignore

For each participant in the dataset,
a directory of derivatives (``sub-<label>/``)
and a visual report (``sub-<label>.html``) are generated.
The log directory contains the citation boilerplate text (``CITATION.md``, ``CITATION.html``,
``CITATION.tex`` and ``CITATION.bib``), which describes the processing that was run and
lists the references to cite.
``dataset_description.json`` is a metadata file in which ASLPrep
records metadata recommended by the BIDS standard.
The atlases directory contains the relevant json, tsv, and nifti files
for each atlas, and the .bidsignore file contains a record of the files that bids-validator ignored.


****************
Processing level
****************

As of version 0.6.0, ASLPrep supports three levels of derivatives:

* ``--level minimal``: This processing mode aims to produce the smallest
  working directory and output dataset possible, while enabling all further
  processing results to be deterministically generated.
  Most components of the `visual reports`_ can be generated at this level,
  so the quality of preprocessing can be assessed.
  ASL-reference-space CBF derivatives will automatically be generated if this level is selected.
* ``--level resampling``: This processing mode aims to produce additional
  derivatives that enable third-party resampling, resampling ASL series
  in the working directory as needed, but these are not saved to the output directory.
  ASL-reference-space CBF derivatives will automatically be generated if this level is selected.
  Currently, 'minimal' and 'resampling' are effectively the same.
* ``--level full``: This processing mode aims to produce all derivatives
  that have previously been a part of the ASLPrep output dataset.
  This is the default processing level.


**************
Visual Reports
**************

*ASLPrep* outputs summary reports are written to ``<output dir>/sub-<label>.html``.
These reports provide a quick way to make visual inspection of the results easy.
`View a sample report. <_static/sub-01.html>`_


************************
Derivatives of *ASLPrep*
************************

Preprocessed, or derivative, data are written to
``<output dir>/sub-<label>/[ses-<label>/]``.
The `BIDS Derivatives`_ specification describes the naming and metadata conventions that we follow.


Anatomical Derivatives
======================

Anatomical derivatives are placed in each subject's ``anat`` subfolder::

   sub-<label>/
      anat/
         <source_entities>[_space-<label>]_desc-preproc_T1w.nii.gz
         <source_entities>[_space-<label>]_desc-preproc_T2w.nii.gz
         <source_entities>[_space-<label>]_desc-brain_mask.nii.gz
         <source_entities>[_space-<label>]_dseg.nii.gz
         <source_entities>[_space-<label>]_label-CSF_probseg.nii.gz
         <source_entities>[_space-<label>]_label-GM_probseg.nii.gz
         <source_entities>[_space-<label>]_label-WM_probseg.nii.gz

Spatially-standardized derivatives are denoted with a space label,
such as ``MNI152NLin2009cAsym``,
while derivatives in the original ``T1w`` space omit the ``space-`` keyword.

T2w images are aligned to the anatomical (``T1w``) space, if found.

.. note::

   T2w derivatives are only generated if FreeSurfer processing is enabled.

Additionally, the following transforms are saved::

   sub-<label>/[ses-<label>/]
      anat/
         sub-<label>_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
         sub-<label>_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5
         sub-<label>_from-<space>_to-T1w_mode-image_xfm.h5
         sub-<label>_from-T1w_to-<space>_mode-image_xfm.h5

If FreeSurfer reconstructions are used, the following surface files are generated::

  sub-<label>/
    anat/
      sub-<label>_hemi-[LR]_white.surf.gii
      sub-<label>_hemi-[LR]_midthickness.surf.gii
      sub-<label>_hemi-[LR]_pial.surf.gii
      sub-<label>_hemi-[LR]_desc-reg_sphere.surf.gii
      sub-<label>_hemi-[LR]_space-fsLR_desc-reg_sphere.surf.gii
      sub-<label>_hemi-[LR]_space-fsLR_desc-msmsulc_sphere.surf.gii
      sub-<label>_hemi-[LR]_space-fsLR_den-32k_desc-preproc_white.surf.gii
      sub-<label>_hemi-[LR]_space-fsLR_den-32k_desc-preproc_midthickness.surf.gii
      sub-<label>_hemi-[LR]_space-fsLR_den-32k_desc-preproc_pial.surf.gii

The registration spheres target ``fsaverage`` and ``fsLR`` spaces. If MSM
is enabled (i.e., the ``--no-msm`` flag is not passed), then the ``msmsulc``
spheres are generated and used for creating ``space-fsLR`` derivatives.

And the affine translation (and inverse) between the original T1w sampling and FreeSurfer's
conformed space for surface reconstruction (``fsnative``) is stored in::

  sub-<label>/
    anat/
      sub-<label>_from-fsnative_to-T1w_mode-image_xfm.txt
      sub-<label>_from-T1w_to-fsnative_mode-image_xfm.txt

Finally, cortical thickness, curvature, and sulcal depth maps are converted to GIFTI
and CIFTI-2::

  sub-<label>/
    anat/
      sub-<label>_hemi-[LR]_thickness.shape.gii
      sub-<label>_hemi-[LR]_curv.shape.gii
      sub-<label>_hemi-[LR]_sulc.shape.gii
      sub-<label>_space-fsLR_den-32k_thickness.dscalar.nii
      sub-<label>_space-fsLR_den-32k_curv.dscalar.nii
      sub-<label>_space-fsLR_den-32k_sulc.dscalar.nii

.. warning::

   GIFTI metric files follow the FreeSurfer conventions and are not modified
   by *ASLPrep* in any way.

   The Human Connectome Project (HCP) inverts the sign of the curvature and
   sulcal depth maps. For consistency with HCP, *ASLPrep* follows these
   conventions and masks the medial wall of CIFTI-2 dscalar files.


.. _fsderivs:

FreeSurfer derivatives
~~~~~~~~~~~~~~~~~~~~~~

If FreeSurfer is run, then a FreeSurfer subjects directory is created in
``<output dir>/sourcedata/freesurfer`` or the directory indicated with the
``--fs-subjects-dir`` flag.
Additionally, FreeSurfer segmentations are resampled into the ASL reference space,
and lookup tables are provided. ::

    <output_dir>/
      sourcedata/
        freesurfer/
          fsaverage{,5,6}/
              mri/
              surf/
              ...
          sub-<label>/
              mri/
              surf/
              ...
          ...
      desc-aparc_dseg.tsv
      desc-aparcaseg_dseg.tsv

Copies of the ``fsaverage`` subjects distributed with the running version of
FreeSurfer are copied into this subjects directory, if any functional data are
sampled to those subject spaces.

Note that the use of ``sourcedata/`` recognizes FreeSurfer derivatives as an input to
the ASLPrep workflow.
This is strictly true when pre-computed FreeSurfer derivatives are provided either in
the ``sourcedata/`` directory or passed via the ``--fs-subjects-dir`` flag;
if ASLPrep runs FreeSurfer, then there is a mutual dependency.


Perfusion Derivatives
=====================

ASL and CBF derivatives are stored in the ``perf/`` subfolder::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>[_space-<label>]_desc-brain_mask.nii.gz  # asl brain mask
         <source_entities>[_space-<label>]_desc-preproc_asl.nii.gz  # preprocessed asl timeseries

.. note::

   The mask file is part of the *minimal* processing level.
   The ASL series is only generated at the *full* processing level.

**Motion correction outputs**.

Head-motion correction (HMC) produces a reference image to which all volumes are aligned,
and a corresponding transform that maps the original ASL series to the reference image::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>_desc-hmc_aslref.nii.gz  # asl reference image for HMC
         <source_entities>_from-orig_to-aslref_mode-image_desc-hmc_xfm.txt  # HMC transforms from raw ASL to aslref

.. note::

   Motion correction outputs are part of the *minimal* processing level.

**Coregistration outputs**.

Registration of the ASL series to the T1w image generates a further reference image and affine
transform::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>_desc-coreg_aslref.nii.gz
         <source_entities>_from-aslref_to-T1w_mode-image_desc-coreg_xfm.txt

.. note::

   Coregistration outputs are part of the *minimal* processing level.

**Fieldmap registration**.

If a fieldmap is used for the correction of an ASL series,
then a registration is calculated between the ASL series and the fieldmap.
If, for example, the fieldmap is identified with ``"B0Identifier": "TOPUP"``,
the generated transform will be named::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>_from-aslref_to-TOPUP_mode-image_desc-fmap_xfm.txt

If the association is discovered through the ``IntendedFor`` field of the
fieldmap metadata, then the transform will be given an auto-generated name::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>_from-aslref_to-auto000XX_mode-image_desc-fmap_xfm.txt

.. note::

   Fieldmap registration outputs are part of the *minimal* processing level.

**Regularly gridded outputs (images)**.
Volumetric output spaces labels (``space-<label>`` above, and in the following) include
``T1w`` and ``MNI152NLin2009cAsym`` (default).
``T1w``-space outputs (the preprocessed ASL series, ASL reference, brain mask,
and CBF derivatives) are aligned to the T1w image but keep the voxel size of the ASL data,
as in fMRIPrep.
Standard-space outputs are written at the resolution requested with the ``res-<label>``
specifier, or at the ASL data's resolution if none is given.

**Surfaces, segmentations and parcellations from FreeSurfer**.
If FreeSurfer reconstructions are used,
the ``(aparc+)aseg`` segmentations are aligned to the subject's T1w space and resampled to the ASL grid,
and the CBF derivatives are resampled to the requested surface meshes::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>_space-T1w_desc-aparcaseg_dseg.nii.gz
         <source_entities>_space-T1w_desc-aseg_dseg.nii.gz
         <source_entities>_hemi-[LR]_space-<label>[_den-<label>][_desc-<label>]_cbf.func.gii

FreeSurfer surface output spaces include ``fsnative`` (full density subject-specific mesh),
``fsaverage`` and the down-sampled meshes ``fsaverage6`` (41k vertices) and
``fsaverage5`` (10k vertices, default).
These are sampled with FreeSurfer's ``mri_vol2surf``.

Other surface templates, such as ``fsLR`` or ``onavg``, may be requested with a density,
for example ``--output-spaces fsLR:den-32k onavg:den-41k``.
These are resampled with the Connectome Workbench, as in fMRIPrep,
using the template's spheres registered to fsLR.
When ``--project-goodvoxels`` is used, the "goodvoxels" mask is applied during
volume-to-surface sampling.

**Grayordinates files**.
`CIFTI <https://www.nitrc.org/forum/attachment.php?attachid=333&group_id=454&forum_id=1955>`_ is
a container format that holds both volumetric (regularly sampled in a grid) and surface
(sampled on a triangular mesh) samples.
Sub-cortical time series are sampled on a regular grid derived from one MNI template, while
cortical time series are sampled on surfaces projected from the Glasser template
:footcite:p:`Glasser_2016`.
If CIFTI outputs are requested (with the ``--cifti-outputs`` argument), the ASL series are also
saved as ``dtseries.nii`` CIFTI2 files::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>_asl.dtseries.nii

CIFTI output resolution can be specified as an optional parameter after ``--cifti-output``.
By default, '91k' outputs are produced and match up to the standard `HCP Pipelines`_ CIFTI
output (91282 grayordinates @ 2mm).
However, '170k' outputs are also possible, and produce higher resolution CIFTI output
(170494 grayordinates @ 1.6mm).

CBF Derivatives
~~~~~~~~~~~~~~~

Cerebral blood flow (CBF) derivatives are generated in each of the volumetric spaces,
surface spaces, and CIFTI2 resolutions requested for the preprocessed ASL data.
For the sake of brevity, we only show the NIfTI outputs below,
but these descriptions and suffixes apply to GIFTI and CIFTI outputs as well.

CBF Outputs::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>[_space-<label>]_cbf.nii.gz  # mean CBF
         <source_entities>[_space-<label>]_desc-timeseries_cbf.nii.gz  # computed CBF timeseries
         <source_entities>[_space-<label>]_att.nii.gz  # arterial transit time (multi-delay data only)
         <source_entities>[_space-<label>]_abv.nii.gz  # arterial blood volume (multi-delay data only)
         <source_entities>[_space-<label>]_abat.nii.gz  # arterial bolus arrival time (multi-delay data only)
         <source_entities>_desc-<lower|upper>_att.nii.gz  # range of equally good arterial transit times (multi-delay data only)
         <source_entities>_desc-<lower|upper>_abat.nii.gz  # range of equally good arterial bolus arrival times (multi-delay data only)
         <source_entities>[_space-<label>]_desc-decayed_abv.nii.gz  # arterial blood volume times exp(-aBAT / T1blood) (multi-delay (P)CASL data only)

For multi-delay data, the model cannot always locate the arterial transit time
or the arterial bolus arrival time, and the ``att`` and ``abat`` maps then hold one value from
the range of values that fit equally well.
The ``desc-lower`` and ``desc-upper`` maps hold the ends of that range.
These maps are only written in the ASL run's native space,
which requires ``--level minimal``, ``--level resampling``, or ``asl`` in ``--output-spaces``,
because the ranges change abruptly between neighboring voxels and slices,
and interpolating them to another space would produce values that are no longer bounds.
See :ref:`multidelay_identifiability` before interpreting these maps.

If ``--scorescrub`` is used::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>[_space-<label>]_desc-scoreTimeseries_cbf.nii.gz  # CBF timeseries denoised with SCORE
         <source_entities>[_space-<label>]_desc-score_cbf.nii.gz  # mean CBF denoised with SCORE
         <source_entities>[_space-<label>]_desc-scrub_cbf.nii.gz  # mean CBF denoised with SCRUB

If ``--basil`` is used::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>[_space-<label>]_desc-basil_cbf.nii.gz  # mean CBF computed with BASIL
         <source_entities>[_space-<label>]_desc-basilGM_cbf.nii.gz  # GM partial volume corrected CBF with BASIL
         <source_entities>[_space-<label>]_desc-basilWM_cbf.nii.gz  # WM partial volume corrected CBF with BASIL
         <source_entities>[_space-<label>]_desc-basil_att.nii.gz  # arterial transit time (multi-delay data only)
         <source_entities>[_space-<label>]_desc-basil_abv.nii.gz  # arterial blood volume (multi-delay data only)
         <source_entities>[_space-<label>]_desc-basilSD_cbf.nii.gz  # posterior SD of desc-basil_cbf
         <source_entities>[_space-<label>]_desc-basilGMSD_cbf.nii.gz  # posterior SD of desc-basilGM_cbf
         <source_entities>[_space-<label>]_desc-basilWMSD_cbf.nii.gz  # posterior SD of desc-basilWM_cbf
         <source_entities>[_space-<label>]_desc-basilSD_att.nii.gz  # posterior SD of desc-basil_att (multi-delay data only)

.. _basil_outputs:

Interpreting the BASIL outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*ASLPrep* runs BASIL through FSL's ``oxford_asl`` command, and
each BASIL derivative is one of the files that ``oxford_asl`` writes.
The table below gives the source file of each derivative, relative to the ``oxford_asl``
output directory.

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Derivative
     - ``oxford_asl`` file
     - Contents
   * - ``desc-basil_cbf``
     - ``native_space/perfusion_calib``
     - CBF in mL/100 g/min, with spatial regularization and without partial volume correction.
       Compare this map with *ASLPrep*'s own ``cbf`` map.
   * - ``desc-basilGM_cbf``
     - ``native_space/pvcorr/perfusion_calib``
     - Gray matter CBF in mL/100 g/min, after partial volume correction.
   * - ``desc-basilWM_cbf``
     - ``native_space/pvcorr/perfusion_wm_calib``
     - White matter CBF in mL/100 g/min, after partial volume correction.
   * - ``desc-basil_att``
     - ``native_space/arrival``
     - Arterial transit time in seconds. Multi-delay data only.
   * - ``desc-basil_abv``
     - ``native_space/aCBV_calib``
     - Arterial blood volume. Multi-delay data only.
   * - ``desc-basilSD_cbf``
     - ``native_space/perfusion_var_calib``
     - Posterior standard deviation of ``desc-basil_cbf``, in mL/100 g/min.
   * - ``desc-basilGMSD_cbf``
     - ``native_space/pvcorr/perfusion_var_calib``
     - Posterior standard deviation of ``desc-basilGM_cbf``, in mL/100 g/min.
   * - ``desc-basilWMSD_cbf``
     - ``native_space/pvcorr/perfusion_wm_var_calib``
     - Posterior standard deviation of ``desc-basilWM_cbf``, in mL/100 g/min.
   * - ``desc-basilSD_att``
     - ``native_space/arrival_var``
     - Posterior standard deviation of ``desc-basil_att``, in seconds. Multi-delay data only.

A voxel in an ASL image usually holds a mixture of gray matter, white matter and CSF,
so the CBF measured in that voxel is a weighted average of the CBF in each tissue.
Partial volume correction :footcite:p:`chappell2011partial` separates that average into
one CBF map for gray matter and one for white matter.
Keep the following in mind when using the corrected maps:

- Each value is the CBF of the *tissue* in that voxel, not of the voxel.
  A voxel that is 30% gray matter and has a ``desc-basilGM_cbf`` value of 60
  holds gray matter with a CBF of 60 mL/100 g/min.
  Do not multiply the map by the gray matter fraction again.
- Each map is only meaningful where its tissue is present.
  Mask ``desc-basilGM_cbf`` with a gray matter mask and ``desc-basilWM_cbf`` with a white matter
  mask before averaging or comparing values.
  Estimates are least reliable in voxels with only a small fraction of the tissue.
- The partial volume estimates are the gray matter and white matter probability maps from the
  anatomical workflow, resampled to the ASL image.
  The corrected maps therefore depend on the quality of the anatomical segmentation and
  of the ASL-to-anatomical registration.
- Corrected and uncorrected values are not interchangeable.
  Gray matter CBF is typically higher after correction, because the uncorrected value
  is diluted by white matter and CSF.
- The arterial transit time and arterial blood volume maps come from the analysis without
  partial volume correction.

BASIL also estimates how uncertain each value is.
*ASLPrep* writes this as a posterior standard deviation map (``SD`` in the ``desc`` entity)
for each CBF map, and for the arterial transit time with multi-delay data.
Each is the square root of the posterior variance that ``oxford_asl`` writes,
so it has the same units as the estimate.
Voxels outside BASIL's analysis mask are set to zero.
Keep the following in mind:

- A small standard deviation means the data constrain the estimate well.
  It does not mean the estimate is correct: it only reflects the noise and the model,
  not errors such as misregistration or a wrong labeling efficiency.
- Relative to the estimate, values are usually larger in white matter and at the edge of
  the brain, where the ASL signal is weaker.
- There is no standard deviation map for the arterial blood volume,
  because ``oxford_asl`` does not write its variance.

If ``--pvc-method local-regression`` is used with single-delay data::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>[_space-<label>]_desc-lrGM_cbf.nii.gz  # GM CBF after local-regression PVC
         <source_entities>[_space-<label>]_desc-lrWM_cbf.nii.gz  # WM CBF after local-regression PVC
         <source_entities>_desc-lrGMSE_cbf.nii.gz  # standard error of desc-lrGM_cbf (aslref only)
         <source_entities>_desc-lrWMSE_cbf.nii.gz  # standard error of desc-lrWM_cbf (aslref only)
         <source_entities>_desc-lr_pvstatus.nii.gz  # fit status of each voxel (aslref only)
         <source_entities>_desc-lr_pvcondition.nii.gz  # condition number of each fit (aslref only)
         <source_entities>_desc-lrGM_pvnoise.nii.gz  # GM noise amplification (aslref only)
         <source_entities>_desc-lrWM_pvnoise.nii.gz  # WM noise amplification (aslref only)
         <source_entities>_desc-lr_pvrmse.nii.gz  # residual of each fit (aslref only)
         <source_entities>_label-GM_desc-psf_probseg.nii.gz  # GM fraction used (aslref only)
         <source_entities>_label-WM_desc-psf_probseg.nii.gz  # WM fraction used (aslref only)
         <source_entities>_label-CSF_desc-psf_probseg.nii.gz  # CSF fraction used (aslref only)
         <source_entities>_desc-lrCSF_cbf.nii.gz  # apparent CSF CBF, for QC (aslref only)
         <source_entities>_desc-lrCSFSE_cbf.nii.gz  # standard error of desc-lrCSF_cbf (aslref only)
         <source_entities>_desc-lrCSF_pvstatus.nii.gz  # status of the CSF QC fit (aslref only)
         <source_entities>_desc-lrCSF_pvnoise.nii.gz  # CSF noise amplification (aslref only)

With ``--scorescrub``, the same correction is applied to the SCORE map,
with ``desc-scoreLR`` in place of ``desc-lr``.

With multi-delay data, the outputs are::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>[_space-<label>]_desc-lrGM_cbf.nii.gz  # GM CBF after local-regression PVC
         <source_entities>[_space-<label>]_desc-lrWM_cbf.nii.gz  # WM CBF after local-regression PVC
         <source_entities>[_space-<label>]_desc-lrGM_att.nii.gz  # GM arterial transit time
         <source_entities>[_space-<label>]_desc-lrWM_att.nii.gz  # WM arterial transit time
         <source_entities>_desc-lrGMLower_att.nii.gz  # shortest GM ATT that fits as well (aslref only)
         <source_entities>_desc-lrGMUpper_att.nii.gz  # longest GM ATT that fits as well (aslref only)
         <source_entities>_desc-lrWMLower_att.nii.gz  # shortest WM ATT that fits as well (aslref only)
         <source_entities>_desc-lrWMUpper_att.nii.gz  # longest WM ATT that fits as well (aslref only)
         <source_entities>_desc-lrGM_pvkinstatus.nii.gz  # GM kinetic fit status (aslref only)
         <source_entities>_desc-lrWM_pvkinstatus.nii.gz  # WM kinetic fit status (aslref only)
         <source_entities>_desc-lrGM_pvkinres.nii.gz  # GM kinetic fit residual (aslref only)
         <source_entities>_desc-lrWM_pvkinres.nii.gz  # WM kinetic fit residual (aslref only)
         <source_entities>_desc-lrCSF_cbf.nii.gz  # signed apparent CSF CBF, for QC (aslref only)
         <source_entities>_desc-lrCSF_pvkinres.nii.gz  # residual of the CSF projection (aslref only)

plus the spatial-regression diagnostics and partial volume maps listed above.
No standard error maps are written for multi-delay data.

.. _lr_pvc_outputs:

Interpreting the local-regression PVC outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Local-regression partial volume correction :footcite:p:`asllani2008regression`
assumes that, within a small kernel around each voxel (by default 5 x 5 voxels in-plane and
one slice), gray and white matter each have one perfusion value, and that CSF has none.
The perfusion signal of the voxels in the kernel is regressed on their gray and white matter
fractions, and the tissue signals are calibrated with the voxel's M0,
the same calibration the uncorrected CBF uses.
The uncorrected CBF maps are always kept; partial volume correction is recommended as an
additional analysis, not a replacement :footcite:p:`chappell2021partial`.

- As for BASIL, each value estimates the CBF of the *tissue* in that voxel, and each map is
  written only where its tissue makes up at least 10% of the voxel.
  Like the uncorrected CBF, it is calibrated with the voxel's own M0, not a tissue-specific
  one. In a mixed voxel that M0 averages the tissues' M0, so where gray and white matter M0
  differ (their water contents differ by roughly 10-20%), each tissue's estimate carries
  part of that difference.
  Elsewhere, and where the fit failed, the value is zero.
  ``desc-lr_pvstatus`` says which: 1 for a gray and white matter fit,
  2 or 3 for a fit of gray or white matter alone (when the other tissue's contribution to that
  estimate was below 5%), and 4 where neither tissue could be estimated.
  Parcellated values average only the voxels with an estimate,
  and their coverage files give the fraction of each parcel that had one.
- The correction smooths CBF within the kernel. Fine detail within a kernel is lost.
- The tissue fractions (``label-*_desc-psf_probseg``) are blurred to the acquisition's
  effective resolution; see :ref:`effective_resolution`.
- ``desc-lrCSF_cbf`` is a quality-control map, not a correction.
  It comes from a separate regression with CSF as a third tissue
  (the gray and white matter maps above come from the two-tissue fit, which assumes CSF has no
  perfusion). Because CSF has no perfusion, this apparent CSF CBF should be close to zero.
  Values far from zero, especially over large CSF spaces, suggest that the partial volume maps
  do not match the ASL data (segmentation or registration errors, or a wrong effective
  resolution), or that the data contain macrovascular signal or model error.
  Real signal from the choroid plexus, which sits in the ventricles and is perfused,
  can also appear here.
  ``desc-lrCSF_pvstatus`` says which tissues each voxel's fit included
  (1: all three; 2: gray matter and CSF; 3: white matter and CSF; 5: CSF alone;
  4: not identifiable), and the CBF QC file (``desc-qualitycontrol_cbf.tsv``) adds its median
  (``lr_csf_cbf_median``) and interquartile range (``lr_csf_cbf_iqr``) over voxels that are at
  least 50% CSF, and the fraction of brain voxels with a two-tissue fit
  (``lr_fraction_two_tissue``).
- ``desc-lrGMSE_cbf`` and ``desc-lrWMSE_cbf`` are standard errors that assume independent
  noise between voxels, conditional on M0. Neighboring ASL voxels are correlated, so these
  are approximate. They are written in the ASL reference space only, because interpolating a
  standard error does not give the standard error of the interpolated estimate.
- The corrected maps are resampled to other spaces with linear interpolation normalized by the
  resampled mask of estimated voxels, so voxels without an estimate do not lower their
  neighbors' values.
- With slice-timed (2D multi-slice) data, the kernel must span a single slice,
  because the tissue signal depends on each slice's delay.

With multi-delay data, the regression is applied to the delta-M of every volume, with one design
for all volumes, which gives a gray matter and a white matter curve across the delays.
Each curve is then fit with the same general kinetic model as the uncorrected data,
calibrated by the voxel's M0, giving tissue CBF and ATT.

- The tissue CBF is non-negative, as for the uncorrected multi-delay fit. Its upper bound is
  raised to 1000 mL/100 g/min, because a tissue's signal calibrated by the voxel's mixed M0 can
  exceed physiological values; ``desc-lr*_pvkinstatus`` flags estimates at either bound.
  The arterial blood volume bound is raised from 0.1 to 1 for the same reason: unmixing assigns
  a voxel's arterial signal to its tissue fraction. A fit with aBV still at that bound is flagged
  (status 5), and neither its CBF nor its ATT is written, because the arterial signal it could
  not absorb would bias both.
- Tissue curves are noisier than the voxel's own curve, and the non-negative fit skews noisy
  estimates upward, so a parcel's mean tissue CBF can be biased up in noisy data; its median is
  more robust.
- Tissue ATT is written only where the tissue's CBF is positive. As for the uncorrected ATT,
  check its bounds (``desc-lr*Lower_att``, ``desc-lr*Upper_att``): where a range of ATT fits
  equally well, for example when ATT is shorter than every delay, the middle of the range is
  reported.
- The unmixing assumes each tissue's signal, including its arrival time, is the same throughout
  the kernel. Arterial (macrovascular) signal does not follow the tissue maps: the kinetic fit's
  arterial term absorbs some of it, but it can bias the tissue estimates, in either direction.
  ``desc-lr*_pvkinres`` is the relative residual of each tissue's kinetic fit, and is the
  diagnostic for these problems; the spatial residual (``desc-lr_pvrmse``) cannot detect them.
- ``desc-lrCSF_cbf`` is the signed least-squares CBF of the CSF curve against a unit-CBF tissue
  curve without an arterial term, calibrated by the voxel's M0, at the voxel's uncorrected ATT
  (or the median uncorrected ATT where it was not estimated, as in most CSF). Like the
  single-delay map, it should be near zero.

.. _effective_resolution:

Choosing the effective resolution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An ASL image is blurrier than its voxel size suggests.
The readout spreads each point of tissue over a region larger than one voxel
(the point-spread function), so each voxel's value is a weighted average of the tissue in and
around it.
The *effective resolution* is the width of that blur: its full width at half maximum (FWHM),
in mm, along each axis of the image.
It is never smaller than the voxel itself.

Partial volume correction needs it because the gray, white matter and CSF fractions it
regresses on must describe the same blurred mixture of tissue that each ASL voxel measured.
If the fractions are sharper than the data, the correction attributes signal to the wrong
tissue near tissue boundaries, which biases gray matter CBF :footcite:p:`petr2018effects`.
The setting only affects the tissue maps used by ``--pvc-method``;
it changes neither the uncorrected CBF nor BASIL.

How blurry an acquisition is depends mostly on its readout:

- **2D EPI** is close to the voxel size :footcite:p:`petr2018effects`.
- **3D GRASE** blurs along the slice direction, because the signal decays over the long echo
  train that encodes the slices. For single-shot GRASE, simulations give 1.38 times the slice
  thickness, and about 1.1 times the voxel size in-plane :footcite:p:`vidorreta2013comparison`.
- **3D (stack-of-)spiral** readouts blur most along the slice direction. For one GE protocol
  with 3.8 x 3.8 x 4 mm voxels, the measured effective resolution was 4.6 x 4.3 x 11.8 mm,
  about three times the slice thickness :footcite:p:`petr2018effects`.

``--asl-effective-resolution`` takes one of:

``nominal`` (default)
   The voxel size: each voxel is treated as measuring exactly the tissue inside it.

``auto``
   Uses the sidecar's ``MRAcquisitionType``, ``PulseSequenceType`` and ``NumberShots``.
   Single-shot 3D GRASE gets 1.1 x 1.1 x 1.38 times the voxel size, and 3D spiral gets
   1.17 x 1.17 x 2.95 times the voxel size (the measurement above, with the two in-plane values
   averaged). 2D acquisitions use the voxel size.
   Anything else, including segmented (multi-shot) GRASE, uses the voxel size with a warning,
   because no validated value is available.
   The spiral factors come from a single protocol, so treat them as an approximation.

Three numbers
   The FWHM in mm, in the order in-plane, in-plane, slice direction.
   The slice direction is the axis named by ``SliceEncodingDirection``
   (the third image axis if it is absent).
   Each value includes the voxel itself, so the smallest allowed value is the voxel size,
   and a value equal to the voxel size means no extra blur along that axis.
   For example, for the spiral protocol above:
   ``--asl-effective-resolution 4.6 4.3 11.8``.

``nominal`` is the default because assuming too much blur does more harm than assuming none.
Extra assumed blur moves signal from white matter to gray matter: in simulations of 2D EPI,
3D spiral and 3D GRASE readouts, each extra half-voxel of slice blur raised corrected gray
matter CBF by a few percent and lowered white matter CBF by 10-40%.
The readouts' true point-spread functions are sharply peaked with long tails, while the model
is a box convolved with a Gaussian. Matched by FWHM, the model blurs the core too much, so even
the true FWHM of a readout overcorrects. In those simulations, the voxel size was the best
setting or close to it for every readout, including 3D GRASE with a through-plane FWHM of
almost five slices. The published factors used by ``auto`` overcorrected the simulated 3D
spiral badly.

Give a wider value only when the voxel is known to be smaller than what was acquired.
The common case is an image reconstructed or interpolated onto a finer grid than it was
acquired on, for example Siemens 2D EPI with ``Interpolation2D`` (a 64 matrix stored as 128).
There, give the acquired voxel size as the in-plane FWHM.
If you know the effective resolution of your protocol, for example from a measurement like
that of :footcite:t:`petr2018effects`, which smooths a high-resolution tissue-based image until
it best matches the ASL image, compare it with ``nominal`` before relying on it.
The value used, and where it came from, are written to the ``EffectiveResolution`` and
``EffectiveResolutionSource`` fields of the ``label-*_desc-psf_probseg.json`` sidecars.
Blurring added by preprocessing, such as interpolation during motion and distortion correction,
is not included.

BASIL does not model the effective resolution.
*ASLPrep* gives it gray and white matter maps that were interpolated to the ASL grid
without additional blurring, and ``oxford_asl`` uses them as they are.
- For single-delay data, the arterial transit time is fixed rather than estimated,
  so neither it nor its standard deviation is written.
- In standard and anatomical spaces, the maps are resampled like the estimates.
  Interpolation can introduce small negative values near sharp edges.

See the `oxford_asl user guide <https://asl-docs.readthedocs.io/en/latest/oxford_asl_userguide.html>`_
for a full description of these files, and :ref:`the BASIL section <basil_workflow>` for how
*ASLPrep* runs BASIL.


*************
ASL Confounds
*************

For each :abbr:`ASL (arterial spin labeling)` run processed with *ASLPrep*, an
accompanying *confounds* file will be generated.
The confounds are saved as a :abbr:`TSV (tab-separated value)` file::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>_desc-confounds_timeseries.tsv
         <source_entities>_desc-confounds_timeseries.json

These :abbr:`TSV (tab-separated values)` tables look like the example below,
where each row of the file corresponds to one time point found in the
corresponding :abbr:`ASL (arterial spin labeling)` time series::

     std_dvars	dvars	framewise_displacement	trans_x	trans_y	trans_z	rot_x	rot_y	rot_z
     n/a	n/a	n/a	0	0	0	-0.00017029	0	0
     1.168398	17.575331	0.0721193	0	0.0207752	0.0463124	-0.000270924	0	0
     1.085204	16.323904	0.0348966	0	0	0.0457372	0	0	0
     1.01591	15.281561	0.0333937	0.010164	-0.0103568	0.0424513	0	0	0.00019174

Confounds include the six head-motion parameters (three rotations and three translations),
which are common outputs from the head-motion correction (also known as *realignment*).
*ASLPrep* also generates framewise displacement, ``DVARS``, and ``std_dvars``.
Confound variables calculated in *ASLprep* are stored separately for each subject,
session and run in :abbr:`TSV (tab-separated value)` files,
with one column for each confound variable.


*******************
CBF Quality Control
*******************

*ASLPrep* produces a quality control (QC) file for each ASL run::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>_desc-qualitycontrol_cbf.tsv

The following QC measures were estimated: framewise displacement and relative root mean square (relRMS).
Other QC measurers include Dice and Jaccard indices, cross-correlation,
coverage estimates of the coregistration quality of :abbr:`ASL (arterial spin labeling)` and T1w images,
and normalization quality of ASL image to the template.
The quality evaluation index (QEI) was also computed for CBF [Sudipto Dolui 2016].
The QEI is included for objective quality evaluation of CBF maps.
It quantifies the quality of the CBF image based on structural similarity, spatial variability, and
percentage of voxels in gray matter with negative CBF values.

The QC file also reports signal clipping in the acquired data.
Some scanners clip the signal at the top of their dynamic range (for example, at 4095),
so that many voxels share exactly the same, maximal value.
CBF is not reliable in those voxels.
``perc_voxels_clipped`` is the percentage of voxels in the brain mask that reach
the maximum intensity of the ASL file in at least one volume, and
``perc_voxels_clipped_m0scan`` is the same measure for a separate M0 scan.
Both are close to zero for data that are not clipped.


***********************
Parcellated CBF Results
***********************

*ASLPrep* produces parcellated CBP outputs using a series of atlases.

The atlases currently used in *ASLPrep* can be separated into three groups: subcortical, cortical,
and combined cortical/subcortical.
The two subcortical atlases are the Tian atlas (``desc-Tian``; :footcite:t:`tian2020topographic`) and
the CIFTI subcortical parcellation (``desc-HCP``).
The cortical atlases are the Glasser :footcite:p:`Glasser_2016` and the
Gordon :footcite:p:`Gordon_2014`.
The combined cortical/subcortical atlases are 10 different resolutions of the
4S (Schaefer Supplemented with Subcortical Structures) atlas (``desc-4S<*>56Parcels``).

The 4S atlas combines the Schaefer 2018 cortical atlas (version v0143) :footcite:p:`Schaefer_2017`
at 10 different resolutions (100, 200, 300, 400, 500, 600, 700, 800, 900, and 1000 parcels) with
the CIT168 subcortical atlas :footcite:p:`pauli2018high`,
the Diedrichson cerebellar atlas :footcite:p:`king2019functional`,
the HCP thalamic atlas :footcite:p:`najdenovska2018vivo`,
and the amygdala and hippocampus parcels from the HCP CIFTI subcortical parcellation
:footcite:p:`glasser2013minimal`.
The 4S atlas is used in the same manner across three PennLINC BIDS Apps:
ASLPrep, QSIPrep, and XCP-D, to produce synchronized outputs across modalities.
For more information about the 4S atlas, please see https://github.com/PennLINC/AtlasPack.

MNI152NLin6Asym-space atlases are warped to the ASL reference image space before parcellation.
*ASLPrep* will output the MNI152NLin6Asym-space atlases to the output directory,
as outputting the reference image-space versions would produce too many extra outputs.

.. code-block::

   aslprep/
      sourcedata/
         atlases/
            dataset_description.json
            atlas-<label>_description.json
            tpl-MNI152NLin6Asym/
               tpl-MNI152NLin6Asym_atlas-<label>_res-01_dseg.json
               tpl-MNI152NLin6Asym_atlas-<label>_res-01_dseg.tsv
               tpl-MNI152NLin6Asym_atlas-<label>_res-01_dseg.nii.gz

      sub-<label>/[ses-<label>/]
         perf/
            <source_entities>_space-<label>_atlas-<label>[_desc-<basil>]_coverage.tsv
            <source_entities>_space-<label>_atlas-<label>[_desc-<basil|basilGM|basilWM|score|scrub>]_cbf.tsv


**********
References
**********

.. footbibliography::
