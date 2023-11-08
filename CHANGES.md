
## 0.5.1

### 🐛 Bug Fixes
* Make M0 TR optional for BASIL by @tsalo in https://github.com/PennLINC/aslprep/pull/323
* Fix `BolusCutOffTechnique` typo by @tsalo in https://github.com/PennLINC/aslprep/pull/328
### Other Changes
* Remove unnecessary dependencies by @tsalo in https://github.com/PennLINC/aslprep/pull/325
* Document M0 use in reports by @tsalo in https://github.com/PennLINC/aslprep/pull/329


**Full Changelog**: https://github.com/PennLINC/aslprep/compare/0.5.0...0.5.1

## 0.5.0

### 🎉 Exciting New Features
* Add force-ge and force-no-ge parameters by @tsalo in https://github.com/PennLINC/aslprep/pull/313
* Write out QC metadata by @tsalo in https://github.com/PennLINC/aslprep/pull/317
### 🐛 Bug Fixes
* Fix nibabel import by @tsalo in https://github.com/PennLINC/aslprep/pull/301
* Fix the Docker image by @tsalo in https://github.com/PennLINC/aslprep/pull/312
* Remove M0 metadata from multi-PLD GE metadata dictionary by @tsalo in https://github.com/PennLINC/aslprep/pull/315
* Fix CombineMotionParameters to support single-volume volume types by @tsalo in https://github.com/PennLINC/aslprep/pull/316
* Patch regmotoasl by @tsalo in https://github.com/PennLINC/aslprep/pull/320
* Patch regmotoasl again by @tsalo in https://github.com/PennLINC/aslprep/pull/321
### Other Changes
* Replace niworkflows calls with dependency calls (second attempt) by @tsalo in https://github.com/PennLINC/aslprep/pull/299
* Adopt Nipreps configuration and maintenance docs by @tsalo in https://github.com/PennLINC/aslprep/pull/231
* Fix Docker build steps by @tsalo in https://github.com/PennLINC/aslprep/pull/309


**Full Changelog**: https://github.com/PennLINC/aslprep/compare/0.4.0...0.5.0

## 0.4.0

This release re-implements multi-delay PCASL and single-delay Q2TIPS PASL support.

Additionally, we have replaced certain BOLD-specific steps with ASL-compatible alternatives. One such case is slice timing correction. Rather than interpolate the ASL time series to perform slice timing correction (as in BOLD preprocessing), ASLPrep now shifts PLDs according to the slice timing. The other major change is to motion correction, which is now performed separately for each image type in the ASL time series (e.g., control, label, and M0 volumes).

Thanks to Sudipto Dolui, Jian Hu, Jan Petr, Manuel Taso, and Kay Jann for their help and feedback.

<!-- Release notes generated using configuration in .github/release.yml at main -->

### 🛠 Breaking Changes
* Use same TPM threshold for GE and non-GE data by @tsalo in https://github.com/PennLINC/aslprep/pull/263
* Remove slice-timing correction by @tsalo in https://github.com/PennLINC/aslprep/pull/269
### 🎉 Exciting New Features
* Shift PostLabelingDelay(s) by slice times by @tsalo in https://github.com/PennLINC/aslprep/pull/280
* Perform motion correction separately for each image type by @tsalo in https://github.com/PennLINC/aslprep/pull/275
* Support single-PLD Q2TIPS PASL and multi-PLD PCASL by @tsalo in https://github.com/PennLINC/aslprep/pull/268
### 🐛 Bug Fixes
* Support runs with SCORE/SCRUB disabled by @tsalo in https://github.com/PennLINC/aslprep/pull/279
* Reorganize atlases and add them to the package data by @tsalo in https://github.com/PennLINC/aslprep/pull/282
### Other Changes
* Refactor workflow connections by @tsalo in https://github.com/PennLINC/aslprep/pull/261
* Refactor more by @tsalo in https://github.com/PennLINC/aslprep/pull/264
* Remove unused workflows by @tsalo in https://github.com/PennLINC/aslprep/pull/271
* Remove test-made nifti files before storing artifacts by @tsalo in https://github.com/PennLINC/aslprep/pull/274
* Improve plots by @tsalo in https://github.com/PennLINC/aslprep/pull/286
* Remove multi-echo elements by @tsalo in https://github.com/PennLINC/aslprep/pull/294
* Document M0 scaling by @tsalo in https://github.com/PennLINC/aslprep/pull/292


**Full Changelog**: https://github.com/PennLINC/aslprep/compare/0.3.0...0.4.0

## 0.3.0

0.3.0 reflects renewed maintenance for ASLPrep.

The primary focuses of this release are:
(1) fixing easily-fixable bugs,
(2) disabling broken features with difficult-to-fix bugs,
(3) cleaning up the codebase,
(4) expanding test coverage,
and (5) ensuring that outputs are BIDS-compliant.

To that end, we have unfortunately had to temporarily drop support for multi-PostLabelingDelay data,
as well as PASL data with the Q2TIPS BolusCutOffTechnique.
We will work on fixing these features for the next release.

Additionally, this release includes a number of breaking changes.
We have renamed several of the outputs to ensure that they are BIDS-compliant.
These outputs may need to change again in the future,
but that will happen in the next minor release (0.4.0) at the earliest.
See the table below for a list of the changed filenames.

| Description | <0.3.0 | 0.3.0 |
|---|---|---|
| CBF time series | `_cbf.nii.gz` | `_desc-timeseries_cbf.nii.gz` |
| Mean CBF | `_mean_cbf.nii.gz` | `_cbf.nii.gz` |
| CBF time series after SCORE denoising | `_desc-score_cbf.nii.gz` | `_desc-scoreTimeseries_cbf.nii.gz` |
| Mean CBF after SCORE denoising | `_desc-score_mean_cbf.nii.gz` | `_desc-score_cbf.nii.gz` |
| Bolus arrival time/arterial transit time | `_desc-bat_cbf.nii.gz` | `_att.nii.gz` |

Additionally, we have changed the atlases that are used for the parcellated CBF files.
We have brought the atlases in line with those used by [XCP-D](https://xcp-d.readthedocs.io/en/latest/).
However, in the near future, we plan to update these atlases _again_, so that they are synchronized across
ASLPrep, XCP-D, and QSIPrep, so please be aware of that upcoming change.

### 🛠 Breaking Changes
* Rename BIDS-noncompliant derivatives by @tsalo in https://github.com/PennLINC/aslprep/pull/216
* Use correct bolus values in CBF calculation and temporarily disable Q2TIPS, multi-PLD, and no-BolusCutOff processing by @tsalo in https://github.com/PennLINC/aslprep/pull/235
* Disable SCORE/SCRUB for GE/short runs and get GE workflow working by @tsalo in https://github.com/PennLINC/aslprep/pull/248
* Parallelize parcellation and replace atlases by @tsalo in https://github.com/PennLINC/aslprep/pull/254
* Remove "mean" from mean CBF filenames and add "timeseries" to 4D CBF filenames by @tsalo in https://github.com/PennLINC/aslprep/pull/257

### 🎉 Exciting New Features
* Hardcode T1blood for some field strengths and use Zhang 2013 formula for others by @tsalo in https://github.com/PennLINC/aslprep/pull/243
* Estimate labeling efficiency based on ASL type and number of background suppression pulses by @tsalo in https://github.com/PennLINC/aslprep/pull/244

### 🐛 Bug Fixes
* Pin looseversion by @tsalo in https://github.com/PennLINC/aslprep/pull/211
* Correct ordering of probseg maps by @josephmje in https://github.com/PennLINC/aslprep/pull/192
* Convert inversion times to a list for compatility with BASIL interface by @tsalo in https://github.com/PennLINC/aslprep/pull/222
* Pin networkx version by @tsalo in https://github.com/PennLINC/aslprep/pull/227
* Use LabelingEfficiency as BASIL `--alpha` parameter by @tsalo in https://github.com/PennLINC/aslprep/pull/233
* Remove BIDS-noncompliant m0z and cbf files by @tsalo in https://github.com/PennLINC/aslprep/pull/236
* Replace RepetitionTime with RepetitionTimePreparation by @tsalo in https://github.com/PennLINC/aslprep/pull/245
* Index aslcontext rows with "cbf", not "CBF" by @tsalo in https://github.com/PennLINC/aslprep/pull/252

### Other Changes
* Use aslprep_build as base Docker image - attempt 2 by @tsalo in https://github.com/PennLINC/aslprep/pull/204
* Autoformat with black and isort by @tsalo in https://github.com/PennLINC/aslprep/pull/205
* Replace relative imports with absolute ones by @tsalo in https://github.com/PennLINC/aslprep/pull/206
* Remove unused functions and classes by @tsalo in https://github.com/PennLINC/aslprep/pull/207
* Rename PEP8-noncompliant classes by @tsalo in https://github.com/PennLINC/aslprep/pull/208
* Consistently use f-strings for string formatting by @tsalo in https://github.com/PennLINC/aslprep/pull/209
* Work on fixing docstrings and splitting modules by @tsalo in https://github.com/PennLINC/aslprep/pull/210
* Remove unused arguments throughout package by @tsalo in https://github.com/PennLINC/aslprep/pull/212
* Replace smriprep calls with dependency calls by @tsalo in https://github.com/PennLINC/aslprep/pull/217
* Replace pybids calls with dependency calls by @tsalo in https://github.com/PennLINC/aslprep/pull/213
* Use custom formatting for workflow connections by @tsalo in https://github.com/PennLINC/aslprep/pull/220
* Improve documentation by @tsalo in https://github.com/PennLINC/aslprep/pull/218
* Add step to run unit tests by @tsalo in https://github.com/PennLINC/aslprep/pull/221
* Collect associated files at beginning of workflow by @tsalo in https://github.com/PennLINC/aslprep/pull/246
* Refactor confounds and QC metrics by @tsalo in https://github.com/PennLINC/aslprep/pull/256

### New Contributors
* @tsalo made their first contribution in https://github.com/PennLINC/aslprep/pull/204
* @josephmje made their first contribution in https://github.com/PennLINC/aslprep/pull/192

**Full Changelog**: https://github.com/PennLINC/aslprep/compare/0.2.8...0.3.0

## 0.2.7

Bids validation

## 0.2.6

11/12/2020

ENH - Add  GE SCAN processing with deltam or cbf

ENH - supressed freesurfer processing, FSL FLIRT with BBR for registration

ENH - basil and scorescrub as options, cbf computation is only default
