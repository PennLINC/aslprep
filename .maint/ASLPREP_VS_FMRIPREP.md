# Packaging and maintenance: ASLPrep vs fMRIPrep

This document summarizes how ASLPrep’s packaging and maintenance compare to fMRIPrep after the alignment work.

---

## Aligned with fMRIPrep

- **Two-image Docker layout**
  - **Base image** built from `Dockerfile.base` (runtime dependencies only, no Python stack).
  - **Main image** built from `Dockerfile` that uses Pixi to build the Python environment, then assembles the final image on top of the base.
- **Pixi dependency management**
  - Both use `pyproject.toml` + `pixi.lock` for unified conda + PyPI dependency resolution.
  - Base image has no Python stack; the main Dockerfile builds the environment via Pixi.
- **Date-based base image tags**
  - Both use date-stamped tags (e.g. `fmriprep-base:20251006`, `aslprep-base:20260219`), updated only when the base changes. CI builds the base only if the tagged image is missing from the registry.
- **Layout**
  - `docker/files/` (e.g. FreeSurfer exclude list), `scripts/fetch_templates.py`.
- **.dockerignore**
  - Present in both to keep build context small.
- **.maint**
  - Both have `update_authors.py`, `update_ignore_revs.sh`, CONTRIBUTORS/MAINTAINERS/FORMER/PIs. fMRIPrep additionally has `update_changes.sh`. ASLPrep has extra docs (HOW_WE_WORK, ROADMAP, etc.).
- **Release flow**
  - Version in CITATION.cff, optional Zenodo/authors update, tag → build/push images, GitHub Release.

---

## Differences

| Area | fMRIPrep | ASLPrep |
|------|----------|--------|
| **Docker CI** | GitHub Actions (`.github/workflows/docker.yml`): build base if missing, then build and push main to **ghcr.io**. | **CircleCI** only: `build_and_deploy` builds base if missing, then builds and pushes main to **Docker Hub**. No GHA Docker workflow. |
| **Docker registry** | **ghcr.io** (e.g. `ghcr.io/nipreps/fmriprep`, `ghcr.io/nipreps/fmriprep-base`). | **Docker Hub** (`pennlinc/aslprep`, `pennlinc/aslprep-base`). |
| **Docker wrapper** | **fmriprep-docker** package under `wrapper/` for a CLI that runs `docker run` with the right mounts. | **None**; users run `docker run` (or equivalent) directly. |
| **Changelog for releases** | Changelog file in repo + release notes from it / labels (`update_changes.sh`). | **GitHub release notes only** (labels + “Generate release notes”); no changelog file or update step. |

---

## Summary

ASLPrep now matches fMRIPrep’s **structure** (base + main Dockerfiles, Pixi dependency management, date-based base image tags, `.maint` scripts, no wrapper), but keeps **CircleCI** (not GHA) for Docker build/deploy, **Docker Hub** (not ghcr.io), and **GitHub release notes as the sole changelog** for releases.
