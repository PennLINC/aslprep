.. include:: links.rst

####################
Outputs of *ASLPrep*
####################

*ASLPrep* generates three broad classes of outputs:

1. **Visual QA (quality assessment) reports**:
   one :abbr:`HTML (hypertext markup language)` per subject, per session (if applicable),
   that allows the user to conduct a thorough visual assessment of the raw and
   processed data.
   This also includes quality control measures.

2. **Derivatives (preprocessed data and computed CBF):** the input ASL data ready for
   analysis, (e.g., after the various preparation procedures have been applied),
   the brain mask, and ASL images after masking has been applied. Other outputs
   include computed CBF maps and post-processed data such as denoised and partial
   volume-corrected CBF.

3. **Confounds and quality control metrics**: a confounds matrix that inlcudes framewise
   displacement, motion parameters, coregistration and registration quality indices,
   and CBF quality control metrics.


**************
Visual Reports
**************

*ASLPrep* outputs summary reports are written to ``<output dir>/aslprep/sub-<label>.html``.
These reports provide a quick way to make visual inspection of the results easy.
`View a sample report. <_static/sub-01.html>`_


************************
Derivatives of *ASLPrep*
************************

Preprocessed, or derivative, data are written to
``<output dir>/aslprep/sub-<label>/[ses-<label>/]``.
The `BIDS Derivatives RC1`_ specification describes the naming and metadata conventions
that we follow.


Anatomical derivatives
======================

Anatomical derivatives are placed in each subject's ``anat`` subfolder.
These derivatives are the same as smriprep output::

   sub-<label>/[ses-<label>/]
      anat/
         <source_entities>[_space-<label>]_desc-preproc_T1w.nii.gz
         <source_entities>[_space-<label>]_desc-brain_mask.nii.gz
         <source_entities>[_space-<label>]_dseg.nii.gz
         <source_entities>[_space-<label>]_label-CSF_probseg.nii.gz
         <source_entities>[_space-<label>]_label-GM_probseg.nii.gz
         <source_entities>[_space-<label>]_label-WM_probseg.nii.gz

Spatially-standardized derivatives are denoted with a space label,
such as ``MNI152NLin2009cAsym``, while derivatives in
the original ``T1w`` space omit the ``space-`` keyword.

Additionally, the following transforms are saved::

   sub-<label>/[ses-<label>/]
      anat/
         sub-<label>_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
         sub-<label>_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5


ASL and CBF derivatives
=======================

ASL derivatives are stored in the ``perf/`` subfolder.
All derivatives contain ``task-<task_label>`` (mandatory) and ``run-<run_index>`` (optional), and
these will be indicated with ``[specifiers]``::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>[_space-<label>]_aslref.nii.gz  # asl reference image
         <source_entities>[_space-<label>]_desc-brain_mask.nii.gz  # asl brain mask
         <source_entities>[_space-<label>]_desc-preproc_asl.nii.gz  # preprocessed asl timeseries
         <source_entities>[_space-<label>]_cbf.nii.gz  # mean CBF
         <source_entities>[_space-<label>]_desc-timeseries_cbf.nii.gz  # computed CBF timeseries

SCORE and SCRUB Outputs::

         <source_entities>[_space-<label>]_desc-scoreTimeseries_cbf.nii.gz  # CBF timeseries denoised with SCORE
         <source_entities>[_space-<label>]_desc-score_cbf.nii.gz  # mean CBF denoised with SCORE
         <source_entities>[_space-<label>]_desc-scrub_cbf.nii.gz  # mean CBF denoised with SCRUB

BASIL outputs::

         <source_entities>[_space-<label>]_desc-basil_cbf.nii.gz  # mean CBF computed with BASIL
         <source_entities>[_space-<label>]_desc-basilGM_cbf.nii.gz  # GM partial volume corrected CBF with BASIL
         <source_entities>[_space-<label>]_desc-basilWM_cbf.nii.gz  # WM partial volume corrected CBF with BASIL
         <source_entities>[_space-<label>]_att.nii.gz  # bolus arrival time/arterial transit time (in seconds)


**Regularly gridded outputs (images):**
Volumetric output space labels (``<label>`` above, and in the following) include
``T1w`` and ``MNI152NLin2009cAsym`` (default).


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
CBF quality control
*******************

*ASLPrep* produces a quality control (QC) file for each ASL run::

   sub-<label>/[ses-<label>/]
      perf/
         <source_entities>_desc-qualitycontrol_cbf.csv

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
The two subcortical atlases are the Tian atlas (desc-Tian; :footcite:t:`tian2020topographic`) and
the CIFTI subcortical parcellation (desc-HCP).
The cortical atlases are the Glasser :footcite:p:`Glasser_2016` and the
Gordon :footcite:p:`Gordon_2014`.
The combined cortical/subcortical atlases are 10 different resolutions of the
4S (Schaefer Supplemented with Subcortical Structures) atlas.

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
      atlas-<label>_dseg.json
      atlas-<label>_dseg.tsv
      space-MNI152NLin6Asym_atlas-<label>_dseg.nii.gz

      sub-<label>/[ses-<label>/]
         perf/
            <source_entities>_space-<label>_atlas-<label>[_desc-<basil>]_coverage.tsv
            <source_entities>_space-<label>_atlas-<label>[_desc-<basil|basilGM|basilWM|score|scrub>]_cbf.tsv
