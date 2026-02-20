# Pinned dependency versions: ASLPrep vs fMRIPrep

Comparison of pinned / minimum versions between ASLPrep and fMRIPrep. **Flagged** items are version differences that may be worth aligning or documenting.

---

## Python (pyproject.toml)

| Package | ASLPrep | fMRIPrep | Note |
|---------|---------|----------|------|
| **nireports** | >= 25.1.0 | >= 24.1.0 | **FLAG:** ASLPrep has a higher minimum (25.1 vs 24.1). |
| acres | >= 0.2.0 | >= 0.2.0 | Aligned. |
| looseversion | >= 1.3 | >= 1.3 | Aligned. |
| nibabel | >= 5.1.1 | >= 5.1.1 | Aligned. |
| nipype | >= 1.9.0 | >= 1.9.0 | Aligned. |
| nitransforms | >= 25.0.1 | >= 25.0.1 | Aligned. |
| niworkflows | >= 1.14.4 | >= 1.14.4 | Aligned. |
| numpy | ~= 2.2 | >= 2.0 (Pixi: 2.2.*) | Both use 2.2 in container. |
| packaging | >= 24 | >= 24 | Aligned. |
| pandas | >= 2.2 | >= 2.2 | Aligned. |
| psutil | >= 5.4 | >= 5.4 | Aligned. |
| pybids | >= 0.16 | >= 0.16 | Aligned. |
| requests | >= 2.27 | >= 2.27 | Aligned. |
| sdcflows | >= 2.15.0 | >= 2.15.0 | Aligned. |
| smriprep | >= 0.19.2 | >= 0.19.2 | Aligned. |
| templateflow | >= 24.2.2 | >= 24.2.2 | Aligned. |
| toml | >= 0.10 | >= 0.10 | Aligned. |
| transforms3d | >= 0.4.2 | >= 0.4.2 | Aligned. |

**ASLPrep-only (not in fMRIPrep):** fmriprep ~= 25.2.2, importlib_resources, indexed_gzip, networkx ~= 3.3, nilearn ~= 0.11.0, sentry-sdk.

**fMRIPrep-only (not in ASLPrep):** nitime, tedana, codecarbon, APScheduler.

---

## Container / Conda stack

- **ASLPrep:** `docker/environment.yml` (Conda) + Dockerfile.base (FreeSurfer, Workbench, AFNI, etc.).
- **fMRIPrep:** Pixi (`[tool.pixi.dependencies]` in pyproject.toml) + Dockerfile.base (FreeSurfer, MSM; no Workbench in base).

### Overlapping Conda/Pixi-style deps (aligned)

python 3.12, nodejs 20, mkl 2024.2.2, mkl-service 2.4.2, numpy 2.2, scipy 1.15, matplotlib 3.10, pandas 2.2, h5py 3.13, graphviz 12.2, pandoc 3.7, ants 2.6, fsl-bet2 2111.8, fsl-flirt 2111.4, fsl-fast4 2111.3, fsl-fugue 2201.5, fsl-mcflirt 2111.0, fsl-miscmaths 2412.4, fsl-topup 2203.5.

### Flagged differences (container / non-Python)

| Dependency | ASLPrep | fMRIPrep | Note |
|------------|---------|----------|------|
| **git-annex** | *=alldep* (conda) | >= 10.20250721 (container feature, pypi) | Different source and versioning. |

**Connectome Workbench:** Aligned — ASLPrep now installs `connectome-workbench-cli=2.0` from conda-forge in `docker/environment.yml` (same approach as fMRIPrep’s Pixi `connectome-workbench-cli 2.0.*`).

### ASLPrep-only (container)

FSL ASL stack: fsl-basil, fsl-fabber_core, fsl-fabber_models_asl, fsl-oxford_asl.
Deep learning: tensorflow 2.19.1, pytorch 2.8.0.
Other: templateflow, surfa (pip in env).

### fMRIPrep-only (Pixi)

nitime 0.11.*, scikit-image 0.25.*, scikit-learn 1.6.*, connectome-workbench-cli 2.0.*.

---

## Summary of flagged differences

1. **nireports (Python):** ASLPrep requires >= 25.1.0; fMRIPrep requires >= 24.1.0. ASLPrep is stricter.
2. **git-annex (container):** Different source (conda vs pypi) and version scheme; optional for templateflow/DataLad.

Connectome Workbench is now aligned: both use 2.0.x via conda (ASLPrep: `connectome-workbench-cli=2.0` in environment.yml; fMRIPrep: Pixi `connectome-workbench-cli 2.0.*`).
