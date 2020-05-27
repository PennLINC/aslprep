.. include:: links.rst

.. _outputs:

---------------------
Outputs of *ASLPrep*
---------------------
*ASLPrep* generates three broad classes of outcomes:

1. **Visual QA (quality assessment) reports**:
   one :abbr:`HTML (hypertext markup language)` per subject,
   that allows the user a thorough visual assessment of the quality
   of processing and ensures the transparency of *ASLPrep* operation.

2. **Derivatives (preprocessed data)** the input ASL data ready for
   analysis, i.e., after the various preparation procedures
   have been applied.
   For example, :abbr:`INU (intensity non-uniformity)`-corrected versions
   of the T1-weighted image (per subject), the brain mask or ASL data  
   images after head-motion correction, slice-timing correction and aligned into
   the same-subject's T1w space or in some standard space. Other include CBF derivatives 
   and postpostpricessing data such s denoised and partial volume corrected CBF
   

3. **Confounds**: this is a special family of derivatives that can be utilized
   for subsequent quality controls. 


Visual Reports
--------------
*ASLPrep* outputs summary reports, written to ``<output dir>/aslprep/sub-<subject_label>.html``.
These reports provide a quick way to make visual inspection of the results easy.
Each report is self contained and thus can be easily shared with collaborators (for example via email).
`View a sample report. <_static/sample_report.html>`_

Derivatives of *ASLPrep* (preprocessed data)
---------------------------------------------
Preprocessed, or derivative, data are written to
``<output dir>/aslprep/sub-<subject_label>/``.
The `BIDS Derivatives RC1`_ specification describes the naming and metadata conventions we follow.

Anatomical derivatives
~~~~~~~~~~~~~~~~~~~~~~
Anatomical derivatives are placed in each subject's ``anat`` subfolder::

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

And the affine translation (and inverse) between the original T1w sampling and FreeSurfer's
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
        fmriprep/
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
ASL derivatives are stored in the ``asl/`` subfolder.
All derivatives contain ``task-<task_label>`` (mandatory) and ``run-<run_index>`` (optional), and
these will be indicated with ``[specifiers]``::

  sub-<subject_label>/
    asl/
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


**Regularly gridded outputs (images)**.
Volumetric output spaces labels (``<space_label>`` above, and in the following) include
``T1w`` and ``MNI152NLin2009cAsym`` (default).

**Surfaces, segmentations and parcellations from FreeSurfer**.
If FreeSurfer reconstructions are used, the ``(aparc+)aseg`` segmentations are aligned to the
subject's T1w space and resampled to the BOLD grid, and the BOLD series are resampled to the
mid-thickness surface mesh::

  sub-<subject_label>/
    asl/
      sub-<subject_label>_[specifiers]_space-T1w_desc-aparcaseg_dseg.nii.gz
      sub-<subject_label>_[specifiers]_space-T1w_desc-aseg_dseg.nii.gz
      sub-<subject_label>_[specifiers]_space-<space_label>_hemi-[LR].func.gii

Surface output spaces include ``fsnative`` (full density subject-specific mesh),
``fsaverage`` and the down-sampled meshes ``fsaverage6`` (41k vertices) and
``fsaverage5`` (10k vertices, default).

**Grayordinates files**.
`CIFTI <https://www.nitrc.org/forum/attachment.php?attachid=333&group_id=454&forum_id=1955>`_ is
a container format that holds both volumetric (regularly sampled in a grid) and surface
(sampled on a triangular mesh) samples.
Sub-cortical time series are sampled on a regular grid derived from one MNI template, while
cortical time series are sampled on surfaces projected from the [Glasser2016]_ template.
If CIFTI outputs are requested (with the ``--cifti-outputs`` argument), the BOLD series are also
saved as ``dtseries.nii`` CIFTI2 files::

  sub-<subject_label>/
    asl/
      sub-<subject_label>_[specifiers]_asl.dtseries.nii

CIFTI output resolution can be specified as an optional parameter after ``--cifti-output``.
By default, '91k' outputs are produced and match up to the standard `HCP Pipelines`_ CIFTI
output (91282 grayordinates @ 2mm). However, '170k' outputs are also possible, and produce
higher resolution CIFTI output (170494 grayordinates @ 1.6mm).

**Extracted confounding time series**.
For each :abbr:`ASL (arterial spin labelling)` run processed with *ASLPrep*, an
accompanying *confounds* file will be generated.
Confounds_ are saved as a :abbr:`TSV (tab-separated value)` file::

  sub-<subject_label>/
    asl/
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


CBF qulaity control 
-------------------
*ASLPrep* produces a qc file for each asl run::

    sub-<subject_label>/
      asl/
      sub-<subject_label>_[specifiers]_cbfqc.csv

The following quality control (qc) measures were  estimated: framewise displacement 
and and relative root mean square (relRMS).Other qc meaure include dice and jaccard indices, 
cross-correlation and coverage that estimate the coregistration quality of ASL and T1w images 
and normalization quality of ASL to the template. Quality evaluation index (QEI) was 
also computed for CBF [Sudipto Dolui 2016]. The QEI isautomated for objective quality 
evaluation of CBF maps and measured the CBF quality based  on structural similarity,
spatial variability and the percentatge of voxels with negtaive CBF within Grey matter

Confounds
---------
The :abbr:`ASL (arterial spin labelling)` signal measured the cerebral blood flow (CBF) and noise .
ASL like BOLD fMRI, is subjected to artifacts induced by head motion (spin history effects). 
These artifacts introduce noise, which corrupts the measurement and cannot be corrected by simple realignment.
However, the calculation of CBF from ASL data involves subtracting subsequently acquired volumes. 
Thus, head motion over two frames contribute artifactual signal which exacerbates the confound.

*Confounds*  are variables (non-CBF fluctuations)  representing fluctuations with a potential 
for qulaity control or denoising.
Such non-CBF fluctuations may drive spurious results in fMRI data analysis,

The most well established confounding variables in neuroimaging are the six head-motion parameters
(three rotations and three translations) - the common output of the head-motion correction
(also known as *realignment*) of popular fMRI preprocessing software
such as SPM_ or FSL_.
Beyond the standard head-motion parameters, the alsprep pipeline generates a framewise displacement,`DVARS and std_dvars`.
Confounding variables calculated in *aslprep* are stored separately for each subject,
session and run in :abbr:`TSV (tab-separated value)` files - one column for each confound variable.

Confound regressors description
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Basic confounds**. The most commonly used confounding time series:

- Estimated head-motion parameters:
  ``trans_x``, ``trans_y``, ``trans_z``, ``rot_x``, ``rot_y``, ``rot_z`` - the 6 rigid-body motion
  parameters (3 translations and 3 rotation), estimated relative to a reference image;

**Outlier detection**.
These confounds can be used to detect potential outlier time points -
frames with sudden and large motion or intensity spikes.
- ``framewise_displacement`` - is a quantification of the estimated bulk-head motion calculated using
formula proposed by [Power2012]_;
- ``dvars`` - the derivative of RMS variance over voxels (or :abbr:`DVARS (derivative of
RMS variance over voxels)`) [Power2012]_;
- ``std_dvars`` - standardized :abbr:`DVARS (derivative of RMS variance over voxels)`;


Confounds and "carpet"-plot on the visual reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The visual reports provide several sections per task and run to aid designing
a denoising strategy for subsequent analysis.
Some of the estimated confounds are plotted with a "carpet" visualization of the
:abbr:`ASL (Areterial Spin labelling)` time series [Power2016]_.
An example of these plots follows:

    in progress
    


See implementation on :mod:`~aslprep.workflows.bold.confounds.init_bold_confs_wf`.

.. topic:: References


  .. [Friston1996] Friston KJ1, Williams S, Howard R, Frackowiak RS, Turner R,
     Movement‐Related effects in fMRI time‐series. Magnetic Resonance in Medicine. 1996.
     doi:`10.1002/mrm.191035031 <https://doi.org/10.1002/mrm.1910350312>`_

  .. [Power2016] Power JD, A simple but useful way to assess fMRI scan qualities. NeuroImage. 2016.
     doi:`10.1016/j.neuroimage.2016.08.009 <http://doi.org/10.1016/j.neuroimage.2016.08.009>`_
     
  .. more references
