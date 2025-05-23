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
The log directory contains `citation boilerplate`_ text.
``dataset_description.json`` is a metadata file in which ASLPrep
records metadata recommended by the BIDS standard.


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
         <source_entities>_from-orig_to-aslref_mode-image_xfm.txt  # HMC transforms from raw ASL to aslref

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
         <source_entities>_from-aslref_to-TOPUP_mode-image_xfm.nii.gz

If the association is discovered through the ``IntendedFor`` field of the
fieldmap metadata, then the transform will be given an auto-generated name::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>_from-aslref_to-auto000XX_mode-image_xfm.txt

.. note::

   Fieldmap registration outputs are part of the *minimal* processing level.

**Regularly gridded outputs (images)**.
Volumetric output spaces labels (``space-<label>`` above, and in the following) include
``T1w`` and ``MNI152NLin2009cAsym`` (default).

**Surfaces, segmentations and parcellations from FreeSurfer**.
If FreeSurfer reconstructions are used,
the ``(aparc+)aseg`` segmentations are aligned to the subject's T1w space and resampled to the ASL grid,
and the ASL series are resampled to the mid-thickness surface mesh::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>_space-T1w_desc-aparcaseg_dseg.nii.gz
         <source_entities>_space-T1w_desc-aseg_dseg.nii.gz
         <source_entities>_hemi-[LR]_space-<label>_asl.func.gii

Surface output spaces include ``fsnative`` (full density subject-specific mesh),
``fsaverage`` and the down-sampled meshes ``fsaverage6`` (41k vertices) and
``fsaverage5`` (10k vertices, default).

**Grayordinates files**.
`CIFTI <https://www.nitrc.org/forum/attachment.php?attachid=333&group_id=454&forum_id=1955>`_ is
a container format that holds both volumetric (regularly sampled in a grid) and surface
(sampled on a triangular mesh) samples.
Sub-cortical time series are sampled on a regular grid derived from one MNI template, while
cortical time series are sampled on surfaces projected from the [Glasser2016]_ template.
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
---------------

Cerebral blood flow (CBF) derivatives are generated in each of the volumetric spaces,
surface spaces, and CIFTI2 resolutions requested for the preprocessed ASL data.
For the sake of brevity, we only show the NIfTI outputs below,
but these descriptions and suffixes apply to GIFTI and CIFTI outputs as well.

CBF Outputs::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>[_space-<label>]_cbf.nii.gz  # mean CBF
         <source_entities>[_space-<label>]_desc-timeseries_cbf.nii.gz  # computed CBF timeseries
         <source_entities>[_space-<label>]_att.nii.gz  # arterial transit time (multi-PLD data only)

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
         <source_entities>[_space-<label>]_att.nii.gz  # bolus arrival time/arterial transit time (in seconds)


*************
ASL Confounds
*************

For each :abbr:`ASL (arterial spin labelling)` run processed with *ASLPrep*, an
accompanying *confounds* file will be generated.
`CBF Confounds`_ are saved as a :abbr:`TSV (tab-separated value)` file::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>_desc-confounds_timeseries.tsv
         <source_entities>_desc-confounds_timeseries.json

These :abbr:`TSV (tab-separated values)` tables look like the example below,
where each row of the file corresponds to one time point found in the
corresponding :abbr:`ASL (arterial spin labelling)` time series::

     std_dvars	dvars	framewise_displacement	trans_x	trans_y	trans_z	rot_x	rot_y	rot_z
     n/a	n/a	n/a	0	0	0	-0.00017029	0	0
     1.168398	17.575331	0.0721193	0	0.0207752	0.0463124	-0.000270924	0	0
     1.085204	16.323904	0.0348966	0	0	0.0457372	0	0	0
     1.01591	15.281561	0.0333937	0.010164	-0.0103568	0.0424513	0	0	0.00019174

Confounds include the six head-motion parameters (three rotations and three translations),
which are common outputs from the head-motion correction (also known as *realignment*).
*ASLPrep* also generates framewise displacement, `DVARS`, and `std_dvars`.
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
coverage estimates of the coregistration quality of :abbr:`ASL (arterial spin labelling)` and T1w images,
and normalization quality of ASL image to the template.
The quality evaluation index (QEI) was also computed for CBF [Sudipto Dolui 2016].
The QEI is included for objective quality evaluation of CBF maps.
It quantifies the quality of the CBF image based on structural similarity, spatial variability, and
percentage of voxels in gray matter with negative CBF values.


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
ASLPrep, QSIPrep_, and XCP-D, to produce synchronized outputs across modalities.
For more information about the 4S atlas, please see https://github.com/PennLINC/AtlasPack.

MNI152NLin6Asym-space atlases are warped to the ASL reference image space before parcellation.
*ASLPrep* will output the MNI152NLin6Asym-space atlases to the output directory,
as outputting the reference image-space versions would produce too many extra outputs.

.. code-block::

   aslprep/
      atlases/
         atlas-<label>/
            atlas-<label>_dseg.json
            atlas-<label>_dseg.tsv
            space-<label>_atlas-<label>_dseg.nii.gz

      sub-<label>/[ses-<label>/]
         perf/
            <source_entities>_space-<label>_atlas-<label>[_desc-<basil>]_coverage.tsv
            <source_entities>_space-<label>_atlas-<label>[_desc-<basil|basilGM|basilWM|score|scrub>]_cbf.tsv
