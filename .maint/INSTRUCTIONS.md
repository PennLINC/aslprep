# Release instructions for ASLPrep maintainers

Follow these steps before creating a new ASLPrep release. Run all commands from the **repository root** unless noted otherwise.

---

## 1. Prepare the release branch

- Ensure the release branch (e.g. `main`) is up to date and all CI checks pass.
- Confirm that the version you plan to release does not already exist as a tag: `git tag -l`.

---

## 2. Update CITATION.cff

- In `CITATION.cff`, set:
  - **`version`** to the new release version (e.g. `0.7.6`).
  - **`date-released`** to the release date in `YYYY-MM-DD` format.

---

## 3. (Optional) Update authorship and Zenodo

- If you use Zenodo and have a `.zenodo.json` file, refresh it from `.maint` tables and git history:

  ```bash
  python .maint/update_authors.py zenodo
  ```

- Ensure `.maint/CONTRIBUTORS.md` (and optionally `MAINTAINERS.md`, `PIs.md`, `FORMER.md`) are up to date before running.

---

## 4. Base image and Pixi (if applicable)

- The main Docker image is built **FROM** ``pennlinc/aslprep-base:latest`` (see ``Dockerfile``).
- The **base image** is built from **Dockerfile.base** (minimal runtime: FreeSurfer, AFNI, Workbench, C3d, AtlasPack, MSM; no Python/Conda). The Python stack (Pixi env) is built in the main **Dockerfile** and copied into the final image.
- **Pixi** (``pyproject.toml`` ``[tool.pixi.*]`` + ``pixi.lock``) defines the container’s Python/Conda environment. After changing ``[tool.pixi.dependencies]`` or related sections, regenerate the lock file and commit it:
  - **On Linux:** ``pixi lock``
  - **From macOS/Windows (with Docker):** ``./.maint/pixi_lock_in_docker.sh``
  The Docker build requires ``pixi.lock`` to be present in the repo.
- **CircleCI** (job ``build_and_deploy``) builds both images on every run: first ``Dockerfile.base`` → ``pennlinc/aslprep-base:latest``, then **Dockerfile** (Pixi build + templates + base) → ``pennlinc/aslprep:latest``, then pushes both to Docker Hub (when ``DOCKERHUB_TOKEN`` is set).
- For local testing, run ``docker build -f Dockerfile.base -t pennlinc/aslprep-base:latest .`` then ``docker build -t pennlinc/aslprep:dev .``.
- To use a different base tag (e.g. a date tag), set ``ARG BASE_IMAGE=pennlinc/aslprep-base:<tag>`` in **Dockerfile** and ensure that image exists.

## 5. Commit and push the release preparation

- Stage and commit all release-related edits (e.g. `CITATION.cff`, base image tag in Dockerfile if changed).
- Use a clear message, e.g.: `Prepare release 0.7.6`.
- Push the release branch to the remote.

---

## 6. Create and push the version tag

- Create an annotated tag for the new version:

  ```bash
  git tag -a 0.7.6 -m "Release 0.7.6"
  ```

- Push the tag:

  ```bash
  git push origin 0.7.6
  ```

- Pushing the tag will trigger **CircleCI** ``build_and_deploy`` (on the configured branch): build base image (``pennlinc/aslprep-base:latest``), build main image (``pennlinc/aslprep``), and push both to Docker Hub when ``DOCKERHUB_TOKEN`` is set. Tags pushed: ``unstable`` (always), ``latest`` and ``<version>`` (when ``CIRCLE_TAG`` is set).

---

## 7. Create the GitHub Release

- On GitHub, open **Releases** → **Draft a new release**.
- Choose the tag you just pushed (e.g. `0.7.6`).
- The release title can be the version (e.g. `0.7.6`) or a short phrase.
- For the description you can:
  - Use **“Generate release notes”** so GitHub fills it from labels (aligned with `.github/release.yml`), or
  - Copy the new section from `CHANGES.md` for this version.
- Publish the release.

---

## 8. Post-release (optional)

- **Git blame ignores**: If you added new commits to `.git-blame-ignore-revs` (e.g. for bulk style/formatting), run:

  ```bash
  ./.maint/update_ignore_revs.sh
  ```

  Commit the updated `.git-blame-ignore-revs` if it changed.

- **PyPI**: If you publish the package to PyPI, trigger that workflow or upload the sdist/wheel from the tag (e.g. `python -m build` and `twine upload`), following your usual process.

---

## Quick checklist

- [ ] `CITATION.cff` version and date-released set
- [ ] (Optional) Zenodo / authorship updated
- [ ] (If needed) Base image rebuilt / Dockerfile `BASE_IMAGE` updated
- [ ] Changes committed and pushed
- [ ] Version tag created and pushed
- [ ] GitHub Release created and published
