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

## 4. Base image (if applicable)

- The main Docker image is built **FROM** ``pennlinc/aslprep_build:0.0.20`` (see ``Dockerfile``). CI and the Dockerfile use this versioned tag so they use the current base (the ``latest`` tag may be stale on the registry).
- The base image is built from **Dockerfile.base** (runtime + conda env; no aslprep package) and tagged as ``pennlinc/aslprep_build:latest`` and ``pennlinc/aslprep_build:0.0.20`` (or a newer version when released).
- **CircleCI** (job ``build_and_deploy``) builds both images on every run: first ``Dockerfile.base`` → ``pennlinc/aslprep_build:latest``, then **Dockerfile** → ``pennlinc/aslprep:latest``, then pushes to Docker Hub (when ``DOCKERHUB_TOKEN`` is set). The base is pushed only as ``latest``; tests and the Dockerfile use the pinned ``0.0.20`` tag so CI is stable.
- For local testing, run ``docker build -f Dockerfile.base -t pennlinc/aslprep_build:0.0.20 .`` then ``docker build -t pennlinc/aslprep:dev .``.
- When releasing a new base image, bump the tag (e.g. to ``0.0.21``) in **Dockerfile** (``ARG BASE_IMAGE``), **.circleci/config.yml** (``.dockersetup`` image), and this section; then build and push that tag to the registry so CI and local builds use it.

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

- Pushing the tag will trigger **CircleCI** ``build_and_deploy`` (on the configured branch): build base image (``pennlinc/aslprep_build:latest``), build main image (``pennlinc/aslprep``), and push to Docker Hub when ``DOCKERHUB_TOKEN`` is set. Base is pushed as ``latest`` only. Tags pushed: ``unstable`` (always), ``latest`` and ``<version>`` (when ``CIRCLE_TAG`` is set).

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
