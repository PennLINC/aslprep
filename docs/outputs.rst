.. include:: links.rst

.. _outputs:

---------------------
Outputs of *ASLPrep*
---------------------
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

3. **Confounds and quality controls**: a confounds matrix that inlcudes framewise
   displacement, motion parameters, coregistration and registration quality indices,
   and cbf quality controls.


Visual Reports
--------------
*ASLPrep* outputs summary reports are written to ``<output dir>/aslprep/sub-<subject_label>.html``.
These reports provide a quick way to make visual inspection of the results easy.
`View a sample report. <_static/sub-01.html>`_


Derivatives of *ASLPrep*
------------------------
Preprocessed, or derivative, data are written to
``<output dir>/aslprep/sub-<subject_label>/``.
The `BIDS Derivatives RC1`_ specification describes the naming and metadata conventions
that we follow.



Anatomical derivatives
~~~~~~~~~~~~~~~~~~~~~~
Anatomical derivatives are placed in each subject's ``anat`` subfolder::
These derivatives are the same as smriprep output::

  sub-<subject_label>/
    anat/
      sub-<subject_label>[_space-<space_label>]_desc-preproc_T1w.nii.gz
      sub-<subject_label>[_space-<space_label>]_desc-brain_mask.nii.gz
      sub-<subject_label>[_space-<space_label>]_dseg.nii.gz
      sub-<subject_label>[_space-<space_label>]_label-CSF_probseg.nii.gz
      sub-<subject_label>[_space-<space_label>]_label-GM_probseg.nii.gz
      sub-<subject_label>[_space-<space_label>]_label-WM_probseg.nii.gz

Spatially-standardized derivatives are denoted with a space label,
such as ``MNI152NLin2009cAsym``, while derivatives in
the original ``T1w`` space omit the ``space-`` keyword.

Additionally, the following transforms are saved::

  sub-<subject_label>/
    anat/
      sub-<subject_label>_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
      sub-<subject_label>_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5

.. _fsderivs:


ASL and CBF derivatives
~~~~~~~~~~~~~~~~~~~~~~~~
ASL derivatives are stored in the ``perf/`` subfolder.
All derivatives contain ``task-<task_label>`` (mandatory) and ``run-<run_index>`` (optional), and
these will be indicated with ``[specifiers]``::

  sub-<subject_label>/
    perf/
      sub-<subject_label>_[specifiers]_space-<space_label>_aslref.nii.gz # asl reference image
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-brain_mask.nii.gz # asl brain mask 
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-preproc_asl.nii.gz # preprocessed asl timeseries 
      sub-<subject_label>_[specifiers]_space-<space_label>_cbf.nii.gz  # computed cbf timeseries
      sub-<subject_label>_[specifiers]_space-<space_label>_mean_cbf.nii.gz # mean cbf 

SCORE and SCRUB Outputs::

      sub-<subject_label>_[specifiers]_space-<space_label>_desc-score_cbf.nii.gz # cbf timeseries denoised with  SCORE
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-score_mean_cbf.nii.gz # mean cbf denoised with SCORE 
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-scrub_cbf.nii.gz # mean CBF denoised with SCRUB 

BASIL outputs::

      sub-<subject_label>_[specifiers]_space-<space_label>_desc-basil_cbf.nii.gz # cbf computed with BASIL 
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-pvGM_cbf.nii.gz #  GM partial volume corrected cbf with BASIL
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-pvWM_cbf.nii.gz #  WM partial volume corrected cbf with BASIL
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-bat_cbf.nii.gz # bolus arrivsl time (in seconds)


**Regularly gridded outputs (images):**
Volumetric output space labels (``<space_label>`` above, and in the following) include
``T1w`` and ``MNI152NLin2009cAsym`` (default).

**Extracted confounding time series:**
For each :abbr:`ASL (arterial spin labelling)` run processed with *ASLPrep*, an
accompanying *confounds* file will be generated.
`CBF Confounds`_ are saved as a :abbr:`TSV (tab-separated value)` file::

  sub-<subject_label>/
    perf/
      sub-<subject_label>_[specifiers]_desc-confounds_regressors.tsv
      sub-<subject_label>_[specifiers]_desc-confounds_regressors.json

These :abbr:`TSV (tab-separated values)` tables look like the example below,
where each row of the file corresponds to one time point found in the
corresponding :abbr:`ASL (arterial spin labelling)` time series::

     std_dvars	dvars	framewise_displacement	trans_x	trans_y	trans_z	rot_x	rot_y	rot_z
     n/a	n/a	n/a	0	0	0	-0.00017029	0	0
     1.168398	17.575331	0.0721193	0	0.0207752	0.0463124	-0.000270924	0	0
     1.085204	16.323904	0.0348966	0	0	0.0457372	0	0	0
     1.01591	15.281561	0.0333937	0.010164	-0.0103568	0.0424513	0	0	0.00019174


CBF quality control
--------------------
*ASLPrep* produces a quality control (QC) file for each ASL run::

    sub-<subject_label>/
      perf/
      sub-<subject_label>_[specifiers]_cbfqc.csv

The following QC measures were estimated: framewise displacement
and relative root mean square (relRMS). Other QC measurers include dice and jaccard indices,
cross-correlation, coverage estimates of the coregistration quality of :abbr:`ASL (arterial spin labelling)`
and T1w images, and normalization quality of ASL image to the template. The quality
evaluation index (QEI) was also computed for CBF [Sudipto Dolui 2016]. The QEI is
included for objective quality evaluation of CBF maps. It quantifies the quality
of the CBF image based on structural similarity, spatial variability, and
percentage of voxels in gray matter with negative CBF values.

CBF Confounds
--------------
Confounds include the six head-motion parameters (three rotations and three
translations), which are common outputs from the head-motion correction (also
known as *realignment*). *ASLPrep* also generates framewise displacement, `DVARS`
and `std_dvars`.
Confound variables calculated in *ASLprep* are stored separately for each subject,
session and run in :abbr:`TSV (tab-separated value)` files, with one column for
each confound variable.
