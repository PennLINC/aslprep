# Maintenance instructions for ASLPrep

Run all commands from the **repository root** unless noted otherwise.

---

## Updating runtime dependencies in Dockerfile.base

`Dockerfile.base` contains non-Python runtime dependencies that are downloaded as
pre-built binaries: FreeSurfer, AFNI, MSM, AtlasPack, and system libraries.
These change rarely and live outside the pixi ecosystem.

1. **Edit `Dockerfile.base`** to change the relevant version, URL, or tag:
   - **FreeSurfer**: Update the version in the `curl` URL and the exclude file path.
   - **AFNI**: The current setup downloads the latest `linux_openmp_64.tgz` tarball
     (no pinned version). To pin, switch to a versioned URL.
   - **AtlasPack**: Update the tag (e.g. `FROM pennlinc/atlaspack:0.1.0 AS atlaspack`).
   - **Ubuntu base**: Update the tag (e.g. `ARG BASE_IMAGE=ubuntu:jammy-20250730`).
   - **System packages**: Add/remove `apt-get install` entries.

2. **Bump the date tag** in `Dockerfile` to today's date:

   ```dockerfile
   ARG BASE_IMAGE=pennlinc/aslprep-base:<YYYYMMDD>
   ```

3. **Commit and push.** The CI `build_and_deploy` job runs
   `docker manifest inspect` against the new tag. Since the image does not exist
   yet, it will build `Dockerfile.base`, push it with the date tag and `latest`,
   then build the main image on top of it.

4. **Verify** the CI `build_and_deploy` job succeeds and the new base image
   appears on Docker Hub at `pennlinc/aslprep-base:<YYYYMMDD>`.

For local testing:

```bash
docker build -f Dockerfile.base -t pennlinc/aslprep-base:$(date +%Y%m%d) .
docker build --target aslprep -t pennlinc/aslprep:dev .
```

---

## Updating individual Python dependencies in pyproject.toml

There are two sections that define Python-ecosystem dependencies:

| Section | Managed by | Examples |
|---------|-----------|----------|
| `[project.dependencies]` | PyPI (pip) | fmriprep, nipype, nibabel, sdcflows |
| `[tool.pixi.dependencies]` | conda (pixi) | python, numpy, FSL tools, ANTs, pytorch |

Both are resolved together by pixi when generating `pixi.lock`.

1. **Edit the version specifier** in `pyproject.toml`. For example:
   - PyPI: change `"niworkflows >= 1.14.4"` to `"niworkflows >= 1.15.0"`
     in `[project.dependencies]`.
   - conda: change `scipy = "1.15.*"` to `scipy = "1.16.*"`
     in `[tool.pixi.dependencies]`.

2. **Commit and open a pull request.**

3. **The `pixi-lock.yml` GitHub Action runs automatically** on the PR. It checks
   whether `pixi.lock` is stale and, if so, regenerates it and pushes the
   updated lockfile back to the PR branch.

4. **Review the lockfile update.** Because `pixi.lock` is marked as binary in
   `.gitattributes`, GitHub will not render a diff. You can inspect locally with
   `git diff HEAD~1 -- pixi.lock` if needed.

5. **Merge the PR** once CI passes.

> **Note:** `pixi lock` only runs on Linux (the project specifies
> `platforms = ["linux-64"]`). You cannot regenerate the lockfile locally on
> macOS. The CI workflow handles this for you.

---

## Updating all Python dependencies at once with pixi update

`pixi update` resolves every dependency to the latest version compatible with the
specifiers in `pyproject.toml` and rewrites `pixi.lock`. **This must be done on
a Linux machine** because the project only targets `linux-64`.

1. **SSH into a Linux machine** (or use a Linux container, GitHub Codespace, etc.).

2. **Clone and check out a branch:**

   ```bash
   git clone https://github.com/PennLINC/aslprep.git
   cd aslprep
   git checkout -b update-deps
   ```

3. **Install pixi** (if not already installed):

   ```bash
   curl -fsSL https://pixi.sh/install.sh | bash
   ```

4. **Run `pixi update`:**

   ```bash
   pixi update
   ```

   To update only a specific package:

   ```bash
   pixi update numpy
   ```

5. **Review, commit, and push:**

   ```bash
   git diff pixi.lock | head -200
   git add pixi.lock
   git commit -m "Update pixi lockfile"
   git push -u origin update-deps
   ```

6. **Open a PR.** The `pixi-lock.yml` action will confirm the lockfile is
   already up to date. CI tests will validate nothing is broken.

> **Tip:** `pixi update` can only update within the version ranges declared in
> `pyproject.toml`. To widen a constraint (e.g. `scipy = "1.15.*"` to
> `scipy = "1.16.*"`), edit `pyproject.toml` first, then run `pixi update`.

---

## Releasing a new version

### 1. Prepare the release branch

- Ensure the release branch (e.g. `main`) is up to date and all CI checks pass.
- Confirm that the version you plan to release does not already exist as a tag: `git tag -l`.

### 2. Update CITATION.cff

- In `CITATION.cff`, set:
  - **`version`** to the new release version (e.g. `0.7.6`).
  - **`date-released`** to the release date in `YYYY-MM-DD` format.

### 3. (Optional) Update authorship and Zenodo

- If you use Zenodo and have a `.zenodo.json` file, refresh it from `.maint` tables and git history:

  ```bash
  python .maint/update_authors.py zenodo
  ```

- Ensure `.maint/CONTRIBUTORS.md` (and optionally `MAINTAINERS.md`, `PIs.md`, `FORMER.md`) are up to date before running.

### 4. Base image (if applicable)

- The main Docker image is built **FROM** `pennlinc/aslprep-base:<YYYYMMDD>` (see `Dockerfile`). The tag is a date stamp (e.g. `20260219`) indicating when the base was last built.
- The base image is built from **Dockerfile.base** (runtime dependencies only; no Python stack) and pushed to Docker Hub with both the date tag and `latest`.
- **CircleCI** (job `build_and_deploy`) checks if the base image already exists in the registry. If it does, the base build is skipped. If not, it builds from `Dockerfile.base` and pushes it. This means the base is only rebuilt when you bump the date tag.
- To release a new base image:
  1. Update the date tag in **Dockerfile** (`ARG BASE_IMAGE=pennlinc/aslprep-base:YYYYMMDD`).
  2. Commit and push. The next CI run will detect the missing image and build/push it.

### 5. Commit and push the release preparation

- Stage and commit all release-related edits (e.g. `CITATION.cff`, base image tag in Dockerfile if changed).
- Use a clear message, e.g.: `Prepare release 0.7.6`.
- Push the release branch to the remote.

### 6. Create and push the version tag

- Create an annotated tag for the new version:

  ```bash
  git tag -a 0.7.6 -m "Release 0.7.6"
  ```

- Push the tag:

  ```bash
  git push origin 0.7.6
  ```

- Pushing the tag will trigger **CircleCI** `build_and_deploy`: build base image if missing, build main image (`pennlinc/aslprep`), and push to Docker Hub when `DOCKERHUB_TOKEN` is set. Tags pushed: `unstable` (always), `latest` and `<version>` (when `CIRCLE_TAG` is set).

### 7. Create the GitHub Release

- On GitHub, open **Releases** -> **Draft a new release**.
- Choose the tag you just pushed (e.g. `0.7.6`).
- The release title can be the version (e.g. `0.7.6`) or a short phrase.
- For the description you can:
  - Use **"Generate release notes"** so GitHub fills it from labels (aligned with `.github/release.yml`), or
  - Copy the new section from `CHANGES.md` for this version.
- Publish the release.

### 8. Post-release (optional)

- **Git blame ignores**: If you added new commits to `.git-blame-ignore-revs` (e.g. for bulk style/formatting), run:

  ```bash
  ./.maint/update_ignore_revs.sh
  ```

  Commit the updated `.git-blame-ignore-revs` if it changed.

- **PyPI**: If you publish the package to PyPI, trigger that workflow or upload the sdist/wheel from the tag (e.g. `python -m build` and `twine upload`), following your usual process.

### Release checklist

- [ ] `CITATION.cff` version and date-released set
- [ ] (Optional) Zenodo / authorship updated
- [ ] (If needed) Base image date tag bumped in Dockerfile
- [ ] Changes committed and pushed
- [ ] Version tag created and pushed
- [ ] GitHub Release created and published
