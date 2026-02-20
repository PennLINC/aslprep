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

- The main Docker image is built **FROM** ``pennlinc/aslprep-base:<YYYYMMDD>`` (see ``Dockerfile``). The tag is a date stamp (e.g. ``20260219``) indicating when the base was last built.
- The base image is built from **Dockerfile.base** (runtime dependencies only; no Python stack) and pushed to Docker Hub with both the date tag and ``latest``.
- **CircleCI** (job ``build_and_deploy``) checks if the base image already exists in the registry. If it does, the base build is skipped. If not, it builds from ``Dockerfile.base`` and pushes it. This means the base is only rebuilt when you bump the date tag.
- To release a new base image:
  1. Update the date tag in **Dockerfile** (``ARG BASE_IMAGE=pennlinc/aslprep-base:YYYYMMDD``) and **.circleci/config.yml** (``.dockersetup`` image) to today's date.
  2. Commit and push. The next CI run will detect the missing image and build/push it.
- For local testing: ``docker build -f Dockerfile.base -t pennlinc/aslprep-base:$(date +%Y%m%d) .`` then ``docker build --target aslprep -t pennlinc/aslprep:dev .``.

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

- Pushing the tag will trigger **CircleCI** ``build_and_deploy``: build base image if missing, build main image (``pennlinc/aslprep``), and push to Docker Hub when ``DOCKERHUB_TOKEN`` is set. Tags pushed: ``unstable`` (always), ``latest`` and ``<version>`` (when ``CIRCLE_TAG`` is set).

---

## 7. Create the GitHub Release

- On GitHub, open **Releases** → **Draft a new release**.
- Choose the tag you just pushed (e.g. `0.7.6`).
- The release title can be the version (e.g. `0.7.6`) or a short phrase.
- For the description you can:
  - Use **"Generate release notes"** so GitHub fills it from labels (aligned with `.github/release.yml`), or
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
- [ ] (If needed) Base image date tag bumped in Dockerfile and .circleci/config.yml
- [ ] Changes committed and pushed
- [ ] Version tag created and pushed
- [ ] GitHub Release created and published
