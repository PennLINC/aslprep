# ASLPrep vs fMRIPrep: Processing Steps and Workflows

This document summarizes how ASLPrep’s processing steps and workflow structure compare to fMRIPrep’s. ASLPrep processes **ASL (perfusion)** data and depends on fMRIPrep and smriprep; it reuses anatomical and many functional patterns while adding ASL-specific steps (CBF quantification, parcellation, etc.). Keeping this in sync with fMRIPrep helps maintain compatibility and makes it easier to port fMRIPrep improvements.

---

## 1. High-level structure (aligned)

| Aspect | fMRIPrep | ASLPrep |
|--------|----------|---------|
| **Entry** | `init_fmriprep_wf()` → per-subject `init_single_subject_wf()` | `init_aslprep_wf()` → per-subject `init_single_subject_wf()` |
| **Anatomical** | `smriprep`: `init_anat_fit_wf`, template iterator, ds_anat_volumes, FreeSurfer (recon-all, surfaces, CIFTI if enabled) | **Same**: uses smriprep anatomical workflow and outputs |
| **Fieldmaps / SDC** | SDCflows: `find_estimators`, `init_fmap_preproc_wf`, PEPOLAR, ANAT/syn-SDC, precomputed | **Same**: SDCflows with ASL-specific PEPOLAR (e.g. `m0scan`, `sbref`); syn uses `auto_bold_nss=False` for ASL |
| **Per-run modality** | One workflow per **BOLD** run: `init_bold_wf(bold_series=...)` | One workflow per **ASL** run: `init_asl_wf(asl_file=...)` |
| **Data collection** | `collect_data(..., task=..., echo=...)` → `subject_data['bold']` | `collect_data(..., queries={..., 'asl': ...})` → `subject_data['asl']` |
| **Config / output dir** | `config.execution.fmriprep_dir` | `config.execution.aslprep_dir` |
| **Reports** | Subject summary, about, reports-only mode | Same pattern; ASLPrep uses `reports-spec.yml` and renames BOLD→ASL in boilerplate |

So at the top level, ASLPrep is “same layout as fMRIPrep, but with ASL instead of BOLD and extra CBF/parcellation steps.”

---

## 2. Anatomical pipeline (shared)

Both rely on **smriprep** for:

- T1w (and optional T2w/FLAIR) preprocessing, skull-stripping, normalization
- FreeSurfer `recon-all` (optional), surfaces, CIFTI (optional)
- Template iterator and standard-space anatomical outputs
- Surface derivatives (inflated, curv, fsLR resampling, grayordinates)

ASLPrep’s `workflows/base.py` and fMRIPrep’s `workflows/base.py` both call:

- `init_anat_fit_wf`, `init_template_iterator_wf`, `init_ds_anat_volumes_wf`
- FreeSurfer: `init_ds_fs_segs_wf`, `init_surface_derivatives_wf`, `init_ds_surfaces_wf`, `init_ds_surface_metrics_wf`
- CIFTI: `init_gifti_morphometrics_wf`, `init_hcp_morphometrics_wf`, `init_morph_grayords_wf`, `init_resample_surfaces_wf`, `init_ds_grayord_metrics_wf`, `init_ds_fsLR_surfaces_wf`

**ASLPrep-specific:** `workflows/segmentation.py` defines SynthStrip-related helpers (e.g. `init_dl_prep_wf`). Skull-stripping strategy is ultimately driven by smriprep/config (e.g. SynthStrip vs ANTs); this file is ASLPrep’s local hook for that stack.

---

## 3. Modality-specific workflow modules

### 3.1 fMRIPrep: BOLD (`fmriprep/workflows/bold/`)

| Module | Role |
|--------|------|
| **base** | Orchestrates BOLD pipeline: fit → native → volumetric resample (anat/std) → surface/CIFTI → confounds |
| **hmc** | Head motion correction |
| **stc** | Slice timing correction |
| **t2s** | T2* combination for multi-echo BOLD |
| **registration** | BOLD reference → T1w (used inside fit) |
| **reference** | Build BOLD reference image |
| **resampling** | Volumetric (`apply`), surface, fsLR, grayordinates |
| **apply** | `init_bold_volumetric_resample_wf` (resample to T1w/std) |
| **confounds** | Confounds extraction (e.g. FD, aCompCor), carpet plot |
| **fit** | `init_bold_fit_wf` (reference + HMC + STC + SDC + coreg), `init_bold_native_wf` (native space) |
| **outputs** | DataSinks for BOLD derivatives, timing parameters |

BOLD pipeline order (conceptually): **reference → HMC → STC (if used) → SDC → BOLD–anat registration → resample to T1w/std/surfaces → confounds**.

### 3.2 ASLPrep: ASL (`aslprep/workflows/asl/`)

| Module | Role |
|--------|------|
| **base** | Orchestrates ASL pipeline: fit → native → CBF → confounds → CBF reporting → volumetric/surface/CIFTI resample → carpet → parcellation |
| **reference** | ASL reference image |
| **hmc** | Head motion correction (skipped in “GE-style” when volumes ≤ 5) |
| **fit** | `init_asl_fit_wf` (reference + HMC + SDC + ASL–anat registration), `init_asl_native_wf` |
| **apply** | `init_bold_volumetric_resample_wf` (from fMRIPrep) for volumetric; `init_asl_cifti_resample_wf` for CIFTI |
| **resampling** | `init_asl_surf_wf` (surface resampling; ASLPrep-specific wrapper, no BOLD timing) |
| **confounds** | `init_asl_confounds_wf`, `init_cbf_confounds_wf`, carpet plot |
| **cbf** | CBF quantification (e.g. BASIL, SCORE/SCRUB), `init_cbf_wf`, `init_parcellate_cbf_wf` |
| **outputs** | DataSinks for ASL/CBF derivatives (`init_ds_asl_native_wf`, `init_ds_volumes_wf`, `init_ds_ciftis_wf`) |
| **plotting** | `init_cbf_reporting_wf` (CBF QC/reporting) |

ASL pipeline order (conceptually): **reference → HMC (unless GE-style) → SDC → ASL–anat registration → CBF quantification → confounds/CBF reporting → resample to T1w/std/surfaces/CIFTI → carpet → parcellate CBF**.

---

## 4. Processing steps: side-by-side

| Step | fMRIPrep (BOLD) | ASLPrep (ASL) |
|------|------------------|---------------|
| Reference image | BOLD reference from first/echo | ASL reference (or M0 for GE-style) |
| Motion correction | Yes (HMC) | Yes (HMC), or skipped for GE / low-volume runs |
| Slice timing | Yes (optional STC) | **No** (not implemented for ASL) |
| Multi-echo / T2* | Yes (t2s for multi-echo) | **No** (no multi-echo ASL) |
| Susceptibility distortion | SDCflows (fieldmaps / syn-SDC) | Same (SDCflows; ASL + m0scan/sbref in PEPOLAR) |
| Modality→T1w registration | BOLD ref → T1w (coreg) | ASL ref → T1w (coreg) |
| Resample to T1w / standard | `init_bold_volumetric_resample_wf` | **Same** (`init_bold_volumetric_resample_wf` from fMRIPrep) |
| Surface / CIFTI resampling | `init_bold_surf_wf`, fsLR, grayords | `init_asl_surf_wf`, `init_asl_cifti_resample_wf` (ASL-specific) |
| Confounds | BOLD confounds (FD, aCompCor, etc.) | ASL + CBF confounds |
| Downstream quantification | None (BOLD stays as signal) | **CBF**: BASIL, SCORE/SCRUB, parcellation, CBF reports |

So the main **conceptual** differences are: no STC, no multi-echo, and the addition of **CBF** and **parcellation** in ASLPrep. Volumetric resampling is shared; surface/CIFTI is adapted for ASL (e.g. no BOLD timing in ASLPrep’s surf resampling).

---

## 5. Where ASLPrep reuses fMRIPrep

- **Anatomical**: smriprep (same).
- **Fieldmaps / SDC**: SDCflows (same); ASLPrep passes ASL-specific options (e.g. `auto_bold_nss=False` in syn, m0scan/sbref for PEPOLAR).
- **Volumetric resampling**: `fmriprep.workflows.bold.apply.init_bold_volumetric_resample_wf` is called from ASLPrep’s `asl/base.py` for T1w and standard space.
- **Single-subject layout**: Same pattern (anat, then per-run modality workflows, template iterator, CIFTI if enabled).
- **BIDS/derivatives**: Same ideas; ASLPrep uses `datatype='perf'`, `suffix='asl'`, and its own derivative specs.

---

## 6. Where ASLPrep diverges or adds

- **No slice timing correction (STC)** for ASL.
- **No multi-echo / T2*** workflow.
- **CBF workflows**: `asl/cbf.py` (BASIL, SCORE/SCRUB, parcellation), CBF confounds and reporting (`plotting`, `confounds`).
- **ASL reference and fit**: `asl/reference.py`, `asl/fit.py` (ASL-specific reference and registration, no STC/t2s).
- **Surface/CIFTI**: ASLPrep uses `init_asl_surf_wf` and `init_asl_cifti_resample_wf` instead of fMRIPrep’s bold surf/grayords directly (no BOLD timing; ASL-specific wiring).
- **Config**: `aslprep_dir`, `asl2anat_init`, ASL-specific options (e.g. `scorescrub`, `basil`, `m0_scale`).
- **Data**: `subject_data['asl']` (and sbref, m0scan where relevant) instead of `subject_data['bold']`; no task/echo in the same way.
- **map_fieldmap_estimation**: ASLPrep passes a flat list of ASL files and ASL-relevant suffixes (`asl`, `m0scan`, `sbref`) for PEPOLAR/syn; fMRIPrep passes bold_series (list of lists for echoes) and `bold`/`sbref`.

---

## 7. Keeping in sync with fMRIPrep

When fMRIPrep changes, these are the main places to check in ASLPrep:

1. **Anatomical / smriprep**  
   Both use smriprep; keep smriprep (and niworkflows) version compatible and watch for API changes in `init_anat_fit_wf`, template iterator, and surface/CIFTI outputs.

2. **SDCflows / fieldmaps**  
   Same library; if fMRIPrep changes estimator logic or `init_fmap_preproc_wf` usage, mirror any relevant changes in ASLPrep’s `workflows/base.py` (and any ASL-specific syn/PEPOLAR handling).

3. **Volumetric resampling**  
   ASLPrep calls `init_bold_volumetric_resample_wf` from fMRIPrep. If that workflow’s inputs/outputs or name change, update `aslprep/workflows/asl/base.py` (and any DataSinks that consume its outputs).

4. **BOLD fit vs ASL fit**  
   fMRIPrep’s `init_bold_fit_wf` / `init_bold_native_wf` are the BOLD analogues of ASLPrep’s `init_asl_fit_wf` / `init_asl_native_wf`. If fMRIPrep adds new steps (e.g. another correction) or reorders them, consider whether the same step or order makes sense for ASL and update the ASL fit/native workflows accordingly.

5. **Config and execution**  
   If fMRIPrep adds new workflow or execution options, consider adding equivalent options in ASLPrep (e.g. `asl2anat_init` already mirrors `bold2anat_init`).

6. **Reports and boilerplate**  
   ASLPrep replaces “BOLD” with “ASL” in boilerplate and uses its own report spec; keep report structure and required inputs aligned so report generation keeps working after fMRIPrep/smriprep changes.

This file should be updated when either codebase gains or changes major processing steps so that alignment and divergence remain documented.
