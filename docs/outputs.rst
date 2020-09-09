.. include:: links.rst

.. _outputs:

---------------------
Outputs of *ASLPrep*
---------------------
*ASLPrep* generates three broad classes of outcomes:

1. **Visual QA (quality assessment) reports**:
   one :abbr:`HTML (hypertext markup language)` per subject and per session( if applicable),
   that allows the user a thorough visual assessment of the quality
   of processing and ensures the transparency of *ASLPrep* operation.
   This also includes quality conrol measures. 

2. **Derivatives (preprocessed data and computed CBF):** the input ASL data ready for
   analysis, (e.g., after the various preparation procedures have been applied.
   the brain mask or ASL data  images afte. Other output includes comouted CBF maps
   and post-processing data such as denoised and partial volume-corrected CBF. 
   
3. **Confounds and quality controls**: confounds matrixs inlcuding framewise displacement, 
motion parameters, coregistration and registration quality indexes and cbf quality controls. 


Visual Reports
--------------
*ASLPrep* outputs summary reports, written to ``<output dir>/aslprep/sub-<subject_label>.html``.
These reports provide a quick way to make visual inspection of the results easy.
`View a sample report. <_static/sub-01.html>`_


Derivatives of *ASLPrep* 
------------------------
Preprocessed, or derivative, data are written to
``<output dir>/aslprep/sub-<subject_label>/``.
The `BIDS Derivatives RC1`_ specification describes the naming and metadata conventions we follow



Anatomical derivatives
~~~~~~~~~~~~~~~~~~~~~~
Anatomical derivatives are placed in each subject's ``anat`` subfolder::
These derivatives are the same as smriprep output

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

If FreeSurfer reconstructions are used, the following surface files are generated::

  sub-<subject_label>/
    anat/
      sub-<subject_label>_hemi-[LR]_smoothwm.surf.gii
      sub-<subject_label>_hemi-[LR]_pial.surf.gii
      sub-<subject_label>_hemi-[LR]_midthickness.surf.gii
      sub-<subject_label>_hemi-[LR]_inflated.surf.gii

The affine translation (and inverse) between the original T1w sampling and FreeSurfer's
conformed space for surface reconstruction (``fsnative``) is stored in::

  sub-<subject_label>/
    anat/
      sub-<subject_label>_from-fsnative_to-T1w_mode-image_xfm.txt
      sub-<subject_label>_from-T1w_to-fsnative_mode-image_xfm.txt

.. _fsderivs:

FreeSurfer derivatives
~~~~~~~~~~~~~~~~~~~~~~
A FreeSurfer subjects directory is created in ``<output dir>/freesurfer``, or the
directory indicated with the ``--fs-subjects-dir`` flag. ::

    <output_dir>/
        aslprep/
            ...
        freesurfer/
            fsaverage{,5,6}/
                mri/
                surf/
                ...
            sub-<subject_label>/
                mri/
                surf/
                ...
            ...

Copies of the ``fsaverage`` subjects distributed with the running version of
FreeSurfer are copied into this subjects directory, if any functional data are
sampled to those subject spaces.

ASL and CBF derivatives
~~~~~~~~~~~~~~~~~~~~~~~~
ASL derivatives are stored in the ``perf/`` subfolder.
All derivatives contain ``task-<task_label>`` (mandatory) and ``run-<run_index>`` (optional), and
these will be indicated with ``[specifiers]``::

  sub-<subject_label>/
    perf/
      sub-<subject_label>_[specifiers]_space-<space_label>_aslref.nii.gz
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-brain_mask.nii.gz
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-preproc_asl.nii.gz
      sub-<subject_label>_[specifiers]_space-<space_label>_cbf.nii.gz
      sub-<subject_label>_[specifiers]_space-<space_label>_mean_cbf.nii.gz
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-score_cbf.nii.gz
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-score_mean_cbf.nii.gz
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-scrub_cbf.nii.gz
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-basil_cbf.nii.gz
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-pvc_cbf.nii.gz


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

The following QC measures were  estimated: framewise displacement 
and relative root mean square (relRMS). Other QC measurers include dice and jaccard indices, 
cross-correlation, coverage estimates of the coregistration quality of :abbr:`ASL (arterial spin labelling)` 
and T1w images, and normalization quality of ASL to the template. Quality evaluation index (QEI) was 
also computed for CBF [Sudipto Dolui 2016]. The QEI is automated for objective quality 
evaluation of CBF maps and measured the CBF quality based  on structural similarity,
spatial variability, and the percentatge of voxels with negtaive CBF within grey matter

CBF Confounds
--------------
The six head-motion parameters (three rotations and three translations) - 
the common output of the head-motion correction (also known as *realignment*)
*ASLPrep* also generates a framewise displacement,`DVARS` and `std_dvars`.
Confounding variables calculated in *Aslprep* are stored separately for each subject,
session and run in :abbr:`TSV (tab-separated value)` files with one column for 
each confound variable.