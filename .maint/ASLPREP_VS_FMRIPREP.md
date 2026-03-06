# Packaging and maintenance: ASLPrep vs fMRIPrep

This document summarizes how ASLPrep’s packaging and maintenance compare to fMRIPrep.

---

## Aligned with fMRIPrep

- **Two-image Docker layout**
  - **Base image** built from `Dockerfile.base` (runtime dependencies only, no Python stack).
  - **Main image** built from `Dockerfile` that uses Pixi to build the Python environment, then assembles the final image on top of the base.
- **Pixi dependency management**
  - Both use `pyproject.toml` + `pixi.lock` for unified conda + PyPI dependency resolution.
  - Both use a multi-stage Docker build with Pixi environments and a template prefetch stage.
- **Date-based base image tags**
  - Both use date-stamped base tags, and CI only rebuilds the base image when the referenced tag is missing in the registry.
- **Maintenance scaffolding**
  - Both include `.maint/update_authors.py`, `.maint/update_ignore_revs.sh`, and maintainer tables (`CONTRIBUTORS.md`, `MAINTAINERS.md`, `PIs.md`, `FORMER.md`).
- **Release flow**
  - Both rely on version tags + CI image publication and GitHub Releases.

---

## Differences

| Area | fMRIPrep | ASLPrep |
|------|----------|--------|
| **Container CI platform** | GitHub Actions (`.github/workflows/docker.yml`) builds base/main images and publishes to GHCR on non-PR events. | CircleCI (`.circleci/config.yml`) builds `image_prep` + test image/cache flow and publishes via `deploy_docker`. |
| **Container registry** | **ghcr.io** (`ghcr.io/nipreps/fmriprep`, `ghcr.io/nipreps/fmriprep-base`). | **Docker Hub** (`pennlinc/aslprep`, `pennlinc/aslprep-base`). |
| **Lockfile automation trigger** | `pixi-lock.yml` runs lockfile update steps on every `pull_request_target` and attempts push-back unconditionally. | `pixi-lock.yml` now gates lockfile updates to latest-commit edits of `pyproject.toml` / `pixi.lock` and only pushes on same-repo PRs. |
| **Image rebuild trigger logic** | Buildx cache behavior in GHA Docker workflow; no CircleCI-style marker cache key. | `image_prep` rebuild trigger keyed to `Dockerfile` + `pixi.lock` checksum cache (`build-v3`) and `imageprep-success.marker`. |
| **Testing execution model in CI** | Primary test workflow (`.github/workflows/test.yml`) runs tox on Python matrices; Docker workflow is focused on image builds. | CircleCI runs integration/unit tests by invoking `pytest` inside the built `pennlinc/aslprep:test` image. |
| **Docker wrapper package** | Includes `wrapper/` (`fmriprep-docker`) for a helper CLI around `docker run`. | No analogous wrapper package in this repository. |
| **Changelog tooling** | Includes `.maint/update_changes.sh`. | No equivalent changelog-maintenance script in `.maint`. |

---

## Summary

ASLPrep and fMRIPrep are now closely aligned on Docker/Pixi structure and base-image maintenance strategy, but they still differ in CI stack (CircleCI vs GitHub Actions), registry target (Docker Hub vs GHCR), and lockfile/test orchestration details.
