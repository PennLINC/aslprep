# Packaging and maintenance: ASLPrep vs fMRIPrep

This document summarizes how ASLPrep’s packaging and maintenance compare to fMRIPrep after the alignment work (and with CircleCI as the sole Docker builder; no GitHub Actions docker workflow).

---

## Aligned with fMRIPrep

- **Two-image Docker layout**
  - **Base image** built from `Dockerfile.base` (runtime dependencies only, no Python stack).
  - **Main image** built from `Dockerfile` that uses Pixi to build the Python environment, then assembles the final image on top of the base.
- **Layout**
  - `docker/files/` (e.g. FreeSurfer exclude list), `scripts/fetch_templates.py`.
- **.dockerignore**
  - Present in both to keep build context small.
- **.maint**
  - Both have `update_authors.py`, `update_changes.sh`, `update_ignore_revs.sh`, CONTRIBUTORS/MAINTAINERS/FORMER/PIs (ASLPrep also has HOW_WE_WORK, ROADMAP, etc.).
- **Release flow**
  - Version in CITATION.cff, optional Zenodo/authors update, tag → build/push images, GitHub Release.

---

## Differences

| Area | fMRIPrep | ASLPrep |
|------|----------|--------|
| **Docker CI** | GitHub Actions (`.github/workflows/docker.yml`): build base if missing, then build and push main to **ghcr.io**. | **CircleCI** only: `build_and_deploy` builds base + main and pushes both to **Docker Hub**. No GHA Docker workflow. |
| **Python/env in image** | **Pixi** in the main repo (`pyproject.toml` + `pixi.lock`); base has no Python stack; main image built with Pixi in Dockerfile. | **Pixi** in the main repo (`pyproject.toml` + `pixi.lock`); base has no Python stack; main image built with Pixi in Dockerfile. |
| **Docker registry** | **ghcr.io** (e.g. `ghcr.io/nipreps/fmriprep`). | **Docker Hub** (`pennlinc/aslprep`, `pennlinc/aslprep_build`). |
| **Docker wrapper** | **fmriprep-docker** package under `wrapper/` for a CLI that runs `docker run` with the right mounts. | **None**; users run `docker run` (or equivalent) directly. |
| **Changelog for releases** | Changelog file in repo + release notes from it / labels. | **GitHub release notes only** (labels + “Generate release notes”); no changelog update step in release instructions. |
| **Base image tag** | Date-based (e.g. `ghcr.io/nipreps/fmriprep-base:20251006`), updated when the base changes. | **`pennlinc/aslprep_build:0.0.20`** (versioned; CI and Dockerfile use this; `latest` may be stale). Base is rebuilt and pushed on every CircleCI deploy. |

---

## Summary

ASLPrep now matches fMRIPrep’s **structure** (base + main Dockerfiles, `.maint` scripts, no wrapper) and uses **Pixi** for dependency management (same as fMRIPrep), but keeps **CircleCI** (not GHA) for Docker build/deploy, **Docker Hub** (not ghcr.io), **GitHub release notes as the sole changelog** for releases.
