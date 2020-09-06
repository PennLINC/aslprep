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

2. **Derivatives (preprocessed data and computed CBF):** the input ASL data ready for
   analysis, (e.g., after the various preparation procedures
   have been applied.
   For example, :abbr:`INU (intensity non-uniformity)`-corrected versions
   of the T1-weighted image (per subject), the brain mask or ASL data  
   images after head-motion correction, slice-timing correction, and alignment into
   the same-subject's T1w space or other standard space. Other output includes CBF derivatives 
   and post-processing data such as denoised and partial volume-corrected CBF. 
   
3. **Confounds and quality controls**: confounds matrixs inlcuding framewise displacement, 
motion parameters, coregistration and registration quality indexes and cbf quality controls. 


Visual Reports
--------------
*ASLPrep* outputs summary reports, written to ``<output dir>/aslprep/sub-<subject_label>.html``.
These reports provide a quick way to make visual inspection of the results easy.
`View a sample report. <_static/sub-01.html>`_

Derivatives of *ASLPrep* (preprocessed data)
---------------------------------------------
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

ASL derivatives
~~~~~~~~~~~~~~~~~~~~~~
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

**Surfaces, segmentations and parcellations from FreeSurfer:**
If FreeSurfer reconstructions are used, the ``(aparc+)aseg`` segmentations are aligned to the
subject's T1w space and resampled to the BOLD grid, and the BOLD series are resampled to the
mid-thickness surface mesh::

  sub-<subject_label>/
    perf/
      sub-<subject_label>_[specifiers]_space-T1w_desc-aparcaseg_dseg.nii.gz
      sub-<subject_label>_[specifiers]_space-T1w_desc-aseg_dseg.nii.gz
      sub-<subject_label>_[specifiers]_space-<space_label>_hemi-[LR].func.gii

Surface output spaces include ``fsnative`` (full density subject-specific mesh),
``fsaverage``, and the down-sampled meshes ``fsaverage6`` (41k vertices) and
``fsaverage5`` (10k vertices, default).

**Extracted confounding time series:**
For each :abbr:`ASL (arterial spin labelling)` run processed with *ASLPrep*, an
accompanying *confounds* file will be generated.
Confounds_ are saved as a :abbr:`TSV (tab-separated value)` file::

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
-------------------
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

Confounds
---------
The :abbr:`ASL (arterial spin labelling)` signal measures cerebral blood flow (CBF) and noise.
ASL, like BOLD fMRI, is subjected to artifacts induced by head motion (spin history effects). 
These artifacts introduce noise, which corrupts the measurement and cannot be corrected by simple realignment.
However, the calculation of CBF from ASL data involves subtracting subsequently acquired volumes. 
Thus, head motion over two frames contribute artifactual signal which exacerbates the confound.

*Confounds*  are variables (non-CBF fluctuations)  representing fluctuations with a potential 
for quality control or denoising.
Such non-CBF fluctuations may drive spurious results in fMRI data analysis.

The most well established confounding variables in neuroimaging are the six head-motion parameters
(three rotations and three translations) - the common output of the head-motion correction
(also known as *realignment*) of popular fMRI preprocessing software
such as SPM_ or FSL_.
Beyond the standard head-motion parameters, the *ASLPrep* pipeline generates a framewise displacement,`DVARS` and `std_dvars`.
Confounding variables calculated in *Aslprep* are stored separately for each subject,
session and run in :abbr:`TSV (tab-separated value)` files with one column for each confound variable.

Confound regressors description
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Basic confounds:** The most commonly used confounding time series:

- Estimated head-motion parameters:
  ``trans_x``, ``trans_y``, ``trans_z``, ``rot_x``, ``rot_y``, ``rot_z`` 
- the 6 rigid-body motion
  parameters (3 translations and 3 rotation), estimated relative to a reference image;

**Outlier detection:**
These confounds can be used to detect potential outlier time points -
frames with sudden and large motion or intensity spikes.

- ``framewise_displacement`` is a quantification of the estimated bulk-head motion calculated using
  formula proposed by [Power2012]_
 
- ``dvars`` is the derivative of RMS variance over voxels (or :abbr:`DVARS (derivative of
  RMS variance over voxels)`) [Power2012]_
  
- ``std_dvars`` - standardized :abbr:`DVARS (derivative of RMS variance over voxels)`


Confounds and "carpet"-plot on the visual reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The visual reports provide several sections per task and run to aid in designing
a denoising strategy for subsequent analysis.
Some of the estimated confounds are plotted with a "carpet" visualization of the
:abbr:`ASL (Areterial Spin labelling)` time series [Power2016]_.
An example of these plots is in progress.


See implementation on :mod:`~aslprep.workflows.bold.confounds.init_bold_confs_wf`.

.. topic:: References

  

  .. [Friston1996] Friston KJ1, Williams S, Howard R, Frackowiak RS, Turner R,
     Movement‐Related effects in fMRI time‐series. Magnetic Resonance in Medicine. 1996.
     doi:`10.1002/mrm.191035031 <https://doi.org/10.1002/mrm.1910350312>`_

  .. [Glasser2016] Glasser MF, Coalson TS Robinson EC, Hacker CD, Harwell J, Yacoub E, Ugurbil K,
     Andersson J, Beckmann CF, Jenkinson M, Smith SM, Van Essen DC.
     A multi-modal parcellation of human cerebral cortex. Nature. 2016.
     doi:`10.1038/nature18933 <https://doi.org/10.1038/nature18933>`_

  .. [Jenkinson2002] Jenkinson M, Bannister P, Brady M, Smith S. Improved optimization for the
     robust and accurate linear registration and motion correction of brain images. Neuroimage.
     2002. doi:`10.1016/s1053-8119(02)91132-8 <https://doi.org/10.1016/s1053-8119(02)91132-8>`__.


  .. [Power2012] Power JD, Barnes KA, Snyder AZ, Schlaggar BL, Petersen, SA, Spurious but systematic
     correlations in functional connectivity MRI networks arise from subject motion. NeuroImage. 2012.
     doi:`10.1016/j.neuroimage.2011.10.018 <https://doi.org/10.1016/j.neuroimage.2011.10.018>`_

  .. [Power2016] Power JD, A simple but useful way to assess fMRI scan qualities. NeuroImage. 2016.
     doi:`10.1016/j.neuroimage.2016.08.009 <http://doi.org/10.1016/j.neuroimage.2016.08.009>`_
