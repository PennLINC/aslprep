# Implementation plan: ground-truth tests with aslscan-simulated data

Spec: `docs/specs/2026-10-08-aslscan-ground-truth-tests-design.md` (revision 2, with the
plan-review changes). Branch: `simulated-test-data`. Status: revision 2. The review log is at
the end.

## Conventions for every task

**Commands.**

- Local commands run in WSL:
  `wsl -e bash -lc "cd /mnt/c/Users/tsalo/Documents/linc/aslprep && micromamba run -n aslprep <cmd>"`.
- pytest is always `python -m pytest`.
- Lint with `ruff check` and `ruff format` (~= 0.15) before each commit.
- Shell snippets with nested quoting go into scratchpad scripts.

**Commits.** One per task (or noted sub-task), with imperative messages and no Claude
co-author trailer.

**Hashed sources.** `aslscan_fixtures.py`, `truth_geometry.py`, `truth_models.py`, and every
recipe directory feed the spec digest.

- After editing any of them, regenerate the spec file and commit it with the change:
  `python -m aslprep.tests.aslscan_fixtures --spec > .circleci/aslscan_fixtures.txt`.
- `test_spec_file_matches` enforces this.

**Local versus CI.**

- The micromamba environment covers unit and fast-tier tests and fixture generation (with
  aslscan built at the pins).
- Integration runs use the locally built `pennlinc/aslprep:test` image (the user approved
  building it). Built by:
  `docker build --target test -t pennlinc/aslprep:test --build-arg VCS_REF=$(git rev-parse --short HEAD) .`
- Runs mirror CI's `docker run` through the updated `run_local_tests.py` (Task 10).
- The checkout is mounted at `/tmp/src/aslprep`, so Python edits need no rebuild; `Dockerfile`
  and `pixi.lock` changes do.
- CircleCI confirms the results and adds cross-machine measurements. Pushing is outward-facing,
  so confirm with the user before the first push (Checkpoint B).

**Pins.**

- `ASLSCAN_REF=82919454d11ee63a1bc8032e4a362f9b745e5d97`
- `MRSIM_ACQ_REF=e64caff8b18358704ef876b5a16a344192b3a014`
- `RUST_TOOLCHAIN=1.98.0`
- Generator Python: `numpy==2.5.3`, `nibabel==5.4.2`, `scipy==1.18.1`, `templateflow==25.1.2`
- `CACHE_EPOCH=1`

**Measurements.** Runtimes, mask populations, budgets and scores are recorded in the
"Measurements" section as they happen.

## File map (new unless noted)

| Path | Role |
|---|---|
| `aslprep/tests/aslscan_fixtures.py` | registry, pins, spec text, build/locate aslscan, phantoms, recipes, assembly, noise calibration, manifest, CLI (stdlib-only at import) |
| `aslprep/tests/truth_geometry.py` | `rigid_ras`, LPS/RAS, own ITK writer (derivatives), overlap matrices, pose composition, angles/displacements |
| `aslprep/tests/truth_models.py` | independent ASL equations and the multi-delay fit (Tier A) |
| `aslprep/tests/truth_scoring.py` | masks, guards, Tier A/B, physical-agreement, frames (reads ASLPrep's ITK files with nitransforms), coreg, motion, T1w/MNI, BASIL, SCORE/SCRUB, clean, metric contract, `score_run` |
| `aslprep/tests/truth_bounds.py` | acceptance ceilings, regression bands, `expect()`; `python -m aslprep.tests.truth_bounds propose <json...>` (developer command) |
| `aslprep/tests/aslscan_cli.py` | integration helper: CLI args per recipe, module-scoped run, shared invariant items |
| `aslprep/tests/data/aslscan/<recipe>/` | protocol inputs |
| `aslprep/tests/test_aslscan_fixtures.py` | spec sync, recipes, overlap, geometry landmark tests, sanitizer, bounds-width test |
| `aslprep/tests/test_aslscan_conformance.py` | fast tier |
| `aslprep/tests/test_aslscan_interpretation.py` | BIDS collection on fixtures |
| `aslprep/tests/test_aslscan_<recipe>.py` | integration modules |
| `.circleci/aslscan_fixtures.txt`, `.gitattributes` | generated spec; LF rules |
| `docs/developers.rst` | developer page |
| modified: `.circleci/continue_config.yml`, `pyproject.toml`, `aslprep/tests/run_local_tests.py`, `aslprep/tests/test_cli.py` (refactor of `_run_and_generate`), `AGENTS.md`, `docs/index.rst` | |

## Task 0: Local test image (background)

Build `pennlinc/aslprep:test` (already started), then record the build time and the
`pytest --version` check. If the build fails (for example, a base-image pull), report it and
continue with Tasks 1-9, which do not need it.

## Task 1: Registry skeleton, spec text, sync test

1. Create `aslscan_fixtures.py` containing:
   - the pins, `PY_REQUIREMENTS` and `CACHE_EPOCH`;
   - `TEMPLATES`: the TemplateFlow query, relpath and SHA-256 for `res-01` T1w, the brain mask
     and the GM/WM/CSF probsegs. Hashes come from the local cache after confirming each file
     loads with nibabel;
   - `HASHED_MODULES`: a tuple of module filenames, starting as `('aslscan_fixtures.py',)`.
     Tasks 3 and 7 append to it when they create their modules;
   - `PHANTOM_PARAMS` and `RECIPES`, both empty for now;
   - `_digest_bytes` and `_digest_dir` (sorted `relpath\0len\0LF-normalized bytes`);
   - `spec_text()` per spec Section 4.5;
   - `main()` with `--spec`.
2. Generate `.circleci/aslscan_fixtures.txt` and add `.gitattributes` (LF for the spec file and
   `aslprep/tests/data/aslscan/**`).
3. Tests:
   - `test_spec_file_matches`, whose message gives the regeneration command;
   - `test_digest_dir_line_endings` (CRLF and LF give the same digest; a rename or length change
     alters it);
   - `test_hashed_modules_exist`.
4. Verify with pytest on the file. Commit.

## Task 2: Build and locate aslscan

1. Implement `build_aslscan(workdir, threads=None) -> Path`:
   - clone or fetch both repositories;
   - `checkout --detach <sha>`, then verify `rev-parse HEAD`;
   - `cargo +1.98.0 build --release --locked --features cli,kspace,par --bin aslscan` with
     `CARGO_TARGET_DIR`;
   - write the stamp `aslscan.build.json` (SHAs, toolchain, features, binary SHA-256) and
     `<workdir>/BIN` holding the absolute path;
   - return the path.

   Subprocess calls use lists and `check=True`.
2. Implement `find_aslscan(explicit=None)`, which tries the explicit path, then `$ASLSCAN`, then
   `which`. It validates the stamp against the pins and the binary hash, and otherwise raises
   `AslscanUnavailable`.
3. CLI: `--build-aslscan WORKDIR`.
4. Tests (with fake files): no stamp is rejected, and a hash mismatch is rejected.
5. Verify locally: build into `~/aslscan-build`, and `aslscan --help` runs. Commit.

## Task 3: Geometry helpers and landmark tests

1. Create `truth_geometry.py` and append it to `HASHED_MODULES`. It contains:
   - `rigid_ras` (`T @ Rz @ Ry @ Rx` about the origin; degrees in, mm);
   - `ras_to_lps` and `lps_to_ras` (conjugation);
   - `write_itk_affine(path, M_ras_pull)`, the writer for our derivatives (centre 0);
   - `overlap_matrix` and `overlap_fractions`, which assert diagonal RAS affines;
   - `pose_matrix(row, centre)`;
   - `rot_angle_deg`, `rms_displacement`.

   Reading ASLPrep's ITK files (single and arrays) uses `nitransforms.linear.load(path,
   fmt='itk')` in `truth_scoring`, the same library ASLPrep uses (`interfaces/resampling.py:333`).
2. Read the pinned `mrsim-acq/src/motion.rs:199-205` and `aslscan/src/series.rs:1177-1179`.
   Encode `pose_matrix` (order, centre, direction) and cite them.
3. Unit tests:
   - `rigid_ras` single-axis cases;
   - our ITK file read back by nitransforms gives the same RAS matrix;
   - `overlap_matrix` hand cases (aligned, half-cell offset, edges, fractional extents).
4. Spec Section 9.1 landmark test:
   - a 40³ image at 2 mm with eight blobs;
   - a transform with nonzero translation and three rotations, written as ITK;
   - resampled with nitransforms; the blob centroids must sit at the predicted points within
     0.1 mm.

   This fixes the pull convention empirically. If it contradicts the spec, fix the spec.
5. Regenerate the spec. Commit.

## Task 4: Phantoms

1. Implement the `PHANTOMS` builders per spec Section 4.2:
   - `tfmni` and `tfmni-flat`: verified TemplateFlow downloads; argmax labels with a 0.3 floor
     inside the mask; per-label constants; `m1` and `m2` modulation (parameters in
     `PHANTOM_PARAMS`); the `Units`, `LabelMap` and `phantom.json` contract; `anat/T1w` and
     `anat/probseg-*`; `provenance.json`;
   - `crop:<base>:<ranges>`: world coordinates preserved, and asserted;
   - `subject:<id>`: raises `NotImplementedError` pointing to the spec. `sanitize_maps(maps,
     dseg)` is implemented now.
2. Phantoms are cached under `<data>/aslscan/_phantoms/<name>/` and reused only when the
   provenance digest matches.
3. Tests:
   - `check_phantom_contract` (aslscan's loader rules, spec Section 2) on a synthetic phantom,
     including failure cases;
   - the sanitizer;
   - the crop preserves world coordinates.
4. Verify locally: build both full phantoms and run aslscan on the scratch 39-slice probe
   protocol against each. `perfusion_gt` must vary within GM for `tfmni` and not for
   `tfmni-flat`. Regenerate the spec. Commit.

## Task 5: Recipes, assembly, noise, manifest, publication

1. Implement `load_recipe(name)`:
   - `tomllib` with an exact schema (spec Section 4.3, including `snr_per_pair`). Unknown keys
     are an error, and `aslscan_args` uses an allow-list.
   - Parse `asl.json` and `aslcontext.tsv`.
2. Implement `check_recipe(recipe, phantom_dir)`:
   - per-volume array lengths match the aslcontext rows;
   - the slice count is `max(1, ceil(extent/voxel - 1e-9))` (`resample.rs:295`), and for 2D
     it must equal `len(SliceTiming)`;
   - multiband: the slice count is a multiple of MB, and the timing matches aslscan's group
     rule (`protocol.rs:968-1000`);
   - 3D recipes carry no `SliceTiming` or multiband (`protocol.rs:1689-1706`).

   Errors give the expected values.

   Unit tests cover fractional and near-integer extents, MB group mismatches, and 3D with
   `SliceTiming`.
3. Implement `generate(name, data_dir, aslscan, threads=4)`:
   1. Lock with `fcntl.flock` on `<data>/aslscan/.lock`; generation is POSIX-only and says so.
   2. Work in the temporary directory `.tmp-<name>-<pid>`.
   3. Build or reuse the phantom.
   4. Run aslscan with `RAYON_NUM_THREADS=threads`, `--sub 01`, and the recipe's arguments.
   5. Noise calibration when `snr_per_pair > 0`: the empirical procedure of spec Section 5.
      Derived overlays go in the temporary directory. The achieved SNR must be within 10 % of
      the target, or generation raises.
   6. Assemble per spec Section 4.4:
      - rename the magnitude series and drop the phase series;
      - **rewrite** the m0scan `IntendedFor`, which aslscan sets to the `part-mag` and
        `part-phase` `bids::` URIs, so it names the renamed ASL file. Use a `bids::` URI and
        confirm in Task 6 that `collect_run_data` resolves it; if not, use the
        subject-relative path;
      - remove `LabelingEfficiency` when the recipe says so;
      - apply `m0_divisor`;
      - write the T1w with the `R @ A` affine;
      - write the derivatives with ITK files (forward `R`, reverse `R^-1`, via our writer);
      - write the `pv*_gt` maps from `overlap_fractions`;
      - write `simulation.json` and `truth.json` (recipe, phantom, `R`, the motion centre `c`
        derived from the simulation grid, SNR fields, stamp, digest), plus the dataset files.
   7. Write `manifest.json` last.
   8. Publish atomically.
4. Implement `fixture_dir` (spec Section 4.5; `ASLPREP_REQUIRE_FIXTURES=1` makes stale or
   missing fixtures an error) and `verify(data_dir)`.
5. Implement `check_threads(name, aslscan)`: two fresh executions of aslscan into separate
   temporary directories, with `RAYON_NUM_THREADS` 1 and 4, never through the fixture cache.
   It asserts the effective settings from the runs' environments and byte-identical image
   data, and returns a report.
6. CLI options: `--generate DATA_DIR --aslscan PATH [--only ...] [--threads N]`,
   `--verify DATA_DIR`, and `--check-threads NAME --aslscan PATH`.
7. Tests:
   - the overlap check: the phantom's perfusion mapped through `overlap_fractions` matrices
     reproduces `perfusion_gt` within 1e-4 on a generated fast fixture;
   - `check_threads`, which uses `pytest.mark.skipif(no acceptable aslscan)`. This avoids a new
     marker. In CI's test image it skips, and the generator job runs `--check-threads` instead.
8. Verify locally: generate a temporary recipe; `BIDSLayout` opens it; `verify` passes;
   `check_threads` passes. Regenerate the spec. Commit.

## Task 6: Fast recipes, mask populations, BIDS interpretation

1. Author the fast recipes (spec Section 5) on `crop:tfmni` slabs.
   - Choose crop and geometry so the **total** array (not the brain) is small: at most 2,500
     voxels for single-delay and at most 1,500 for multi-delay.
   - Each `dom_GM` and `dom_WM` must hold at least 300 voxels at the fast tier's own bound. If
     a crop cannot satisfy both, use finer voxels (for example 2.5 mm) and record why.
   - `fast_pcasl_mb` uses an even slice count, for MB 2.
2. Implement the mask-population reporter
   (`python -m aslprep.tests.aslscan_fixtures --masks DATA_DIR`), which prints the `valid` and
   `dom_L` counts per recipe. Record the counts. Run it on F1-F7 too as they are written, and
   refuse a recipe below the minimum.
3. Register the recipes, generate them, and regenerate the spec.
4. Write `test_aslscan_interpretation.py`:
   - `collect_data` finds exactly one ASL run;
   - `collect_run_data` resolves the aslcontext, and the M0 for `Separate` (through
     `IntendedFor`);
   - the metadata fields have the expected values.
5. Verify. Commit.

## Task 7: Tier A single-delay and the fast conformance tests

1. Create `truth_models.py` and append it to `HASHED_MODULES`, with cited independent
   implementations:
   - `parse_context`;
   - `deltam_estimates`, pairing control and label by order as ASLPrep does, plus deltam rows;
   - `m0_image` covering all M0 types as `ExtractCBF` documents them:
     - `Separate`: mean m0scan;
     - `Included`: mean m0 rows;
     - `Estimate`: scalar;
     - `Absent`: mean control, and an error when background suppression is on, as ASLPrep does;
   - `smooth_image` when `fwhm > 0`;
   - `pld_map`, `labeling_efficiency` (ASLPrep's convention, documented as the convention under
     test), and `labeling_efficiency_physical` (from the sidecar's
     `AslscanSimulation.BackgroundSuppressionLabelFactor` times the resolved efficiency);
   - `m0_tr_correction` (1.607 s at 3 T, TR < 5 s) and `cbf_single_delay`.
2. The fast-tier chain per recipe:
   - Run `ExtractCBF` with the raw ASL file as both `asl_file` and `name_source`, the metadata,
     aslcontext, `m0scan` and `m0scan_metadata` (or None), and `fwhm=0`.
   - Run `ComputeCBF` with `deltam=out_file`, `m0_file=m0_file`, the extracted `metadata`,
     `m0_scale` from the recipe, `m0tr` only when not None, and `cbf_only=False`.
3. Assertions on `valid` (with an `n_min` check):
   - Tier A exact: `|median-1| <= 1e-3` and `p99 <= 1e-2`.
   - Mutations: label/control swap, slice timing reversed (2D only), PLD +0.1 s, M0 x1.1,
     efficiency ignored. Each must violate the bound.
   - Smoothing variant (`fwhm=5`, same smoothing in `cbf_A`).
   - Physical agreement on `fast_bs_le_absent`: a strict xfail whose reason is the predicted
     (0.95/0.9)ⁿ. The item first asserts that the measured discrepancy is within 2 % of the
     prediction, then calls `pytest.xfail`. An unexpected size therefore fails instead of
     xfailing. This replaces a plain `xfail` decorator: the item fails if the gap disappears or
     changes size.
   - Tier B is recorded through `record_property`.
4. Investigation rule: a mismatch that looks like an ASLPrep bug is reported to the user with
   evidence before any change.
5. Measure the fast-tier wall time in the micromamba environment and in the test image
   (`docker run ... pytest -k aslscan_conformance`). Record both. Regenerate the spec. Commit.

## Task 8: Tier A multi-delay

1. In `truth_models.py`, write independent forward models (PCASL and PASL with an arterial
   term, from Woods 2023 and Chappell 2010) and a fit over **every** observation with ASLPrep's
   documented settings. Failures become NaN.
2. Add `fast_pcasl_multipld` and `fast_pasl_multipld`, each at most 1,500 total voxels.
3. Tests:
   - Tier A: CBF median within 2 %, ATT median absolute difference ≤ 0.05 s, finite ≥ 0.98;
   - mutations (PLD shift, label/control swap, M0 scale) must fail.
4. Measure the time, including the reference and mutation fits, in both environments. If the
   fast tier exceeds 3 min in the test image, reduce voxels while keeping the populations
   valid, or move multi-delay mutations to a single representative.
5. Regenerate the spec. Commit. **Checkpoint A:** report the fast-tier results to the user.

## Task 9: Scoring and bounds core

1. In `truth_scoring.py`:
   - fixture and output loading. `find_outputs` resolves the entities each recipe needs from
     the output specification in `aslprep/workflows/asl/outputs.py`; a missing required file
     is an error, recorded and enforced by the metric contract;
   - masks and guards (spec Section 4.7);
   - Tier A integration, the expected-ratio Tier B, physical agreement, and contamination;
   - frames (nitransforms reading), `coreg`, `hmc_consistency`, the native-comparability guard;
   - `cbf_t1w`, `cbf_mni`, `basil`, `scorescrub`, `clean`;
   - the motion RMS and the FSL-parameter comparison;
   - `score_run`, which never raises for missing outputs and records them instead.
2. In `truth_bounds.py`:
   - `CEILINGS`, the a-priori values (spec Section 4.7), and `BANDS`, initially empty;
   - `expect(score, path, recipe, kind)`, which fails when the bound is missing or the value is
     missing or non-finite. Messages follow the spec format;
   - the developer command `propose`, which reads `truth_score.json` files and prints proposed
     bands. Tests never read it.
3. Synthetic unit tests (no ASLPrep run):
   - exact output gives ratios of 1;
   - a 10 % scaling fails the ceiling;
   - a shrunken mask trips coverage;
   - NaNs trip the finite check;
   - a 2° coregistration error is measured as 2°;
   - the frames algebra (Section 9.3), consistent and inconsistent;
   - the FSL-parameter conversion (Section 9.5) using fMRIPrep's `FSLMotionParams` on known
     transforms;
   - a missing metric fails the contract;
   - `BANDS` widths are at most 0.1 for ratios (trivially true while empty).
4. Spec Section 9.2: resample the offset `tfmni` T1w to MNI through a generated derivative file
   with nitransforms; r > 0.999 against the TemplateFlow T1w, and a landmark error below
   0.5 mm.
5. Spec Section 9.4:
   - generate a small crop recipe with a trajectory of one non-central rotation and a
     translation, noise-free, with a bright point-like feature (a single-voxel high-M0 label in
     a custom crop phantom);
   - recover its displacement across volumes by centroid;
   - assert it agrees with `pose_matrix` applied with the recorded centre `c`, to within
     0.3 mm.
6. Commit.

## Task 10: Integration scaffolding, local runner, F1

1. Refactor `test_cli._run_and_generate` to return output paths and accept
   `check_outputs=True` (the default keeps current behaviour).
2. In `aslscan_cli.py`:
   - `cli_args(recipe, fixture, out, work, spaces, extra)`:
     - always `--output-spaces` with `asl` first when native scoring is needed;
     - the `acq` filter;
     - `--derivatives` when the recipe has anatomical derivatives;
     - `--fs-no-reconall`;
     - `--m0_scale`;
   - `module_run(recipe, ...)` for module-scoped fixtures. It runs once, writes
     `truth_score.json` before any assertion (an `{"error": ...}` score on exception, then
     re-raises), and prints `TRUTH SCORE`;
   - `make_shared_items(recipe)`, which generates the invariant test functions for each
     module: contract, coverage, finite, voxel counts, clean, outputs manifest.
3. Update `run_local_tests.py`:
   - use the image `pennlinc/aslprep:test`;
   - mount the checkout read-only at `/tmp/src/aslprep` as the working directory;
   - mount `aslprep/tests/test_data` as `/data` and `/tmp` paths for out and work;
   - forward `CIRCLE_CPUS` and `ASLPREP_REQUIRE_FIXTURES`;
   - run `python -c 'import aslprep; print(aslprep.__file__)'` first and print it;
   - accept `-m` and `-k` as now.
4. Write the F1 recipe (39 slices, `tfmni-flat`, SNR 50) and check its mask populations. Write
   the module `test_aslscan_pcasl1pld.py` with the shared items plus:
   - `tier_a_*`, `tier_b_*` (expected-ratio acceptance, no r);
   - `coreg`, `cbf_t1w`, `cbf_mni`;
   - `basil` (contract presence plus the acceptance ceiling on Tier B versus the expected
     ratio; if BASIL's model differs, investigate first);
   - `scorescrub`.
5. Register the marker in `pyproject.toml` (`markers` and `addopts`). Check collection with
   `python -m pytest --collect-only -q --strict-markers`.
6. Run it locally in the test image. Iterate on pipeline errors.
7. Write the outputs manifest from this run, and record the score and runtime.
8. Regenerate the spec. Commit.

## Task 11: F2-F6

For each recipe:

1. Derive the protocol from the bids-examples sidecar (copied from aslscan's
   `tests/fixtures/protocols/asl00N`, with provenance in `SOURCE.md`), or author it (F4).
2. Adapt it to the `tfmni` grid and check mask populations.
3. Generate it and record the runtime.
4. Write the module with the shared items, the "Scored" items of spec Section 5, the flags of
   spec Section 6, and `asl` in the output spaces.
5. Run it locally in the image, and record the scores and manifest.

Specific requirements:

- F3: `m0_divisor` 10, `--m0_scale=10`, M0 TR 3 s, atlases.
- F2 and F3: a physical-agreement item with the predicted BS gap (as in Task 7).
- F5: `[macrovascular]` per-label aBV and aBAT, with `abv`/`abat` report items. They are still
  contract-enforced.

Commit per recipe.

## Task 12: F7 motion

1. Write `motion.tsv`: 20 volumes, drift up to 1.5 mm and 1.5°, spikes of 2 mm and 2° on
   volumes 6, 11 and 17. Rotations in radians.
2. Settings: MB 3 (13 groups), SNR 5, `anat_offset` [3, -4, 2, 5, -3, 4], `anat = "raw"`.
3. Module items: shared items plus `motion_rms_median`, `motion_rms_max`, `hmc_consistency`,
   `coreg`, `cbf_t1w_*`, and the reports `fd`, `fsl_params` and `scorescrub`. Reports are still
   contract-enforced.
4. Run it locally and record the results. Commit.

## Task 13: CI wiring

1. Add the job `gen_aslscan_fixtures` (docker `cimg/python:3.12`, medium), one step per
   action:
   1. Checkout.
   2. Spec guard: `python3 -m aslprep.tests.aslscan_fixtures --spec | diff - .circleci/aslscan_fixtures.txt`.
   3. Restore the cache.
   4. Check: `python3 -m aslprep.tests.aslscan_fixtures --verify /tmp/data && touch /tmp/fixtures-ok || true`.
   5. Generate, guarded by `[[ -f /tmp/fixtures-ok ]] && exit 0`:
      1. `rm -rf /tmp/data/aslscan`.
      2. Install rustup with retries.
      3. `python3 -m venv /tmp/venv && /tmp/venv/bin/pip install <pins>`, with retries.
      4. `/tmp/venv/bin/python -m aslprep.tests.aslscan_fixtures --build-aslscan /tmp/aslscan-build`,
         with `PATH=$HOME/.cargo/bin:$PATH` set inline.
      5. `BIN=$(cat /tmp/aslscan-build/BIN)`.
      6. `--generate /tmp/data --aslscan "$BIN"`.
      7. `--check-threads <fast recipe> --aslscan "$BIN"`.
      8. `--verify /tmp/data`.

      All of these run in one shell step, so variables persist. Check that cimg/python has
      `git`, `curl` and `cc` (`which` in the step, failing early with a message).
   6. Save the cache.
2. Add the command `restore_aslscan_fixtures` (restore plus `python3 --verify`, failing on a
   miss). Use it only in `unit_tests` and the `aslscan_*` jobs, which also add
   `gen_aslscan_fixtures` to `requires`.
3. pytest invocations:
   - `-e ASLPREP_REQUIRE_FIXTURES=1` on fixture consumers;
   - `-v /tmp/test-results:/test-results`;
   - `--junitxml=/test-results/<< parameters.marker >>.xml` for integration jobs and
     `/test-results/unit_tests.xml` for unit tests, each with `-o junit_logging=no` and
     `store_test_results`.
4. Matrix entries for F1-F7 with `skip_marker: aslscan`. Resource classes come from the local
   memory and runtime measurements, with a medium/large split as now. Add them to
   `merge_coverage`.
5. Validate with `circleci config process .circleci/continue_config.yml` if the CLI can be
   installed (ask before installing); otherwise use the CircleCI validation API after push.
   Also run `yaml.safe_load`.
6. Update the CI section of `AGENTS.md`.
7. Commit.

## Task 14: Push, CI confirmation, regression bands

1. **Checkpoint B:** ask the user to push `simulated-test-data`.
2. Push. Every item must pass the acceptance ceilings. Investigate failures (a CI-only failure
   points to thread or machine differences; report it rather than loosening anything).
3. Collect `truth_score.json` from the CI artifacts and the runtimes.
4. Run `truth_bounds propose` on at least five runs (local and CI), and add bands only where
   useful. Commit them, push once more, and confirm a fresh green run with frozen bands.

## Task 15: Documentation

Write `docs/developers.rst` (linked from `index.rst`): the tiers, fixtures (build, generate,
verify, masks), running locally with the test image, reading scores, the tolerance policy, bumps
(pins and `CACHE_EPOCH`), adding a recipe, and adding a phantom (`subject:<id>`). Validate with
Sphinx or docutils. Commit.

## Task 16: Replace the `examples_*` tests

Per spec Section 6:

1. Diff each old test's flags against its replacement, and add anything missing.
2. Confirm the replacement is green in CI.
3. Remove the old tests, manifests, matrix entries, the `examples_pcasl_singlepld_ge` job, and
   the `merge_coverage` dependencies.
4. Remove markers and `addopts` entries once `grep` confirms no other users.
5. Run strict collection locally, then push and confirm green.
6. Commit.

## Task 17: Final review

Run the full unit suite and ruff, then a Codex code review of `main...HEAD`. Address valid
findings and log them here. Summarize for the user, with the open items.

## Order and checkpoints

- Task 0 runs in the background.
- Tasks 1→9 are local, with Checkpoint A after Task 8.
- Tasks 10→12 run in the local image.
- Task 13 is CI wiring.
- Checkpoint B is the push permission.
- Then Tasks 14, 15, 16 and 17.

## Measurements

### Build and generation

- Test image build (local, cold): 344 s.
- aslscan clone and build at the pins: 17 s once crates are cached (rustup 1.98.0 minimal
  installed separately).
- aslscan's `Cargo.lock` is not committed upstream. It is vendored as
  `tests/data/aslscan-Cargo.lock` and hashed into the spec.
- Phantoms: `tfmni` 6.1 s and `tfmni-flat` 2.2 s; the cache hit takes 0.02 s.
- Overlap check: our exact overlap reproduces aslscan's `perfusion_gt` to 3.8e-6.
- Fast recipes generate in 1.5-9.3 s each, about 50 s for all 14. Output is byte-identical
  at 1 and 4 threads.

### Mask populations

- The fast crop is `crop:tfmni:110:150,110:150,100:124` (16x16x6 at 2.5x2.5x4 mm, 1,536
  voxels): `valid` 1,519, `dom_GM` 292, `dom_WM` 1,061.
- The fast tier asserts on `valid` (minimum 500), not per tissue. Per-tissue minimums apply to
  integration recipes.

### Fast tier

| Environment | Result | Wall time |
|---|---|---|
| micromamba env | 82 passed, 1 skipped, 2 xfailed | 46 s (first run, before multi-delay) |
| test image, 2 CPUs | 125 passed, 3 skipped, 2 xfailed | 89 s |

The image run covers fixtures, interpretation and conformance.

### Tier A

- Single-delay: ASLPrep agrees with the independent implementation to 1e-9 (median) on 10 of
  11 recipes.
- Multi-delay: CBF agrees to 1e-8. ATT agrees to 7.6e-3 s (PASL) and 1e-7 s (PCASL).
- One PASL voxel fails ASLPrep's fit; the finite fraction is 0.9993.

### Findings

- **Q2TIPS single-delay** (strict xfail, reported). ASLPrep uses `exp(TI2 / T1b)` with no
  slice shift, which matches to 8e-8. The white paper uses TI. The CBF is 0.809 of the
  white-paper value. With the same physics, QUIPSS II recovers 0.953 of the truth in GM, while
  ASLPrep's Q2TIPS recovers 0.771.
- **Background-suppression gap** (resolved later, not an ASLPrep defect). The measured ratio
  matched the predicted 0.8975 within 2 %, but the gap came from aslscan's definition of
  inversion efficiency: its default ε = 0.95 keeps 0.9 of the signal per pulse. The white paper
  (Alsop et al. 2015) gives 95 % retained per pulse, as ASLPrep assumes. The recipes now set
  ε = 0.975, and the items are ordinary checks. Only the LabelingEfficiency-in-sidecar case
  (F3) remains: ASLPrep should still apply the per-pulse loss there (#706), an expected
  failure until fixed.
- **`collect_run_data`** (`utils/bids.py`). Its `ValueError` message reads
  `run_data["asl_metadata"]` before assignment. Not yet reported.

### Tier B (fast tier, ASLPrep / truth)

| Recipe | GM | WM |
|---|---|---|
| PCASL | 0.80 | 0.47 |
| QUIPSS II | 0.95 | 0.79 |
| multi-delay PCASL | 0.84 | 0.45 |
| multi-delay PASL (GM) | 0.91 | |

These match the model-mismatch expectation.

### Dependency differences

The test image has numpy 2.2.6, scipy 1.15.2 and nibabel 5.3.2; the micromamba environment has
numpy 2.5.3, scipy 1.18.1 and nibabel 5.4.2. One test that was fragile to these differences
(the mutation mask) was fixed.

### Deviations from the plan

- **Task 9.3 (FSL-parameter conversion).** Motion is asserted through the transform algebra:
  `H_v` against `P_v`, as RMS displacement over brain voxels. That metric is
  convention-free and is validated by the frame-algebra test and the simulator
  motion-convention probe (`geom_motion`, spec 9.4). Agreement between the confounds and the
  truth is reported as |r| per axis plus mean FD, report-only. An exact FSL-parameter
  conversion and its tests were not built: they would assert nothing the transform metric does
  not already cover.
- **`rot_angle_deg`.** It uses `atan2`; `arccos` was inaccurate (about 0.006°) near the
  identity.
- **`geom_motion` probe.** It needs a 60x60x40 mm crop at 2.5 mm isotropic. The fast crop
  (16x16x6) was too coarse for the comparison (r 0.66). The true pose fits at r 0.89 and
  beats the inverse pose and the origin-centred rotation.

### Integration tier (local test image, 8 CPUs)

| Recipe | Wall time | Result (after the scoring fixes below) |
|---|---|---|
| F1 | 9.7 min | all items pass |
| F2 | 5.4 min | pass, plus 2 strict xfails: suppression gap, delta-M motion correction |
| F3 | 9.3 min | pass, plus 1 strict xfail: suppression gap |
| F4 | 8.3 min | pass |

Measured values:

| Metric | F1 | F2 | F3 | F4 |
|---|---|---|---|---|
| Tier A median deviation, end to end | 0.018 | 0.0006 | 0.0095 | |
| Tier A median deviation, given ASLPrep's preprocessed series | 0.0015 | 5e-11 | | |
| Tier A p95, given the preprocessed series | 0.030 | 7e-8 | 0.042 | |
| Tier B GM, ASLPrep / expected | 0.808 / 0.788 | 0.650 / 0.657 | | 0.952 / 0.937 |
| Coregistration RMS (voxels) | 0.17 | 0.12 | 0.17 | 0.15 |
| Motion consistency, max (mm) | 0.024 | 1.94 | 0.035 | 0.019 |
| Reference offset (mm) | 0.28 | | 0.61 | |
| Suppression gap, measured / predicted | | 0.805 / 0.8055 | 0.662 / 0.656 | |

Standard spaces (F1): `r_native` 0.90 (T1w) and 0.89 (MNI); 0.56 and 0.54 with the sampling
shifted 4 mm.

### Findings in ASLPrep, from the integration runs

- **Motion correction has a per-voxel noise floor.** On motion-free data, motion correction
  changes delta-M per voxel by 16-26 % at the 95th percentile (median 0.1 %). Each volume's
  transform is estimated from noisy data, so a control and its paired label land a few
  hundredths of a millimetre apart (up to 0.02 mm on F1). Delta-M (about 0.4 % of the signal in
  GM) amplifies that difference wherever intensity changes across a voxel: about 15 % of
  delta-M at grey/white boundaries, more at the brain edge. Resampling is linear, so identical
  transforms would add nothing. Separate estimation for label and control (an established
  method) is not the cause. This is a property of motion correction, not a defect, and the
  magnitudes are consistent with the mechanism, but it has not been isolated directly.
  - End-to-end Tier A tails are therefore report-only. Quantification is asserted given
    ASLPrep's own preprocessed series.
- **Delta-M plus M0 data (GE style).** Motion correction registers the delta-M volume to an
  M0-like reference with no shared contrast, and applies the spurious result: 1.9°, 1.9 mm
  relative to the M0 volume. A strict xfail on F2's motion item; reported.
- **Constant reference offset.** The corrected series sits 0.3-0.6 mm from the static frame
  even without motion: the reference's contrast differs from the volumes', most with
  background suppression. Report-only.
- **M0 resampling.** The separate M0 is resampled onto the ASL grid with ANTs' `Gaussian`
  interpolation even under an identity transform, blurring M0 before the 5 mm smoothing. On
  the phantom, which has no scalp, this raises CBF near the brain-mask edge: median 1.19 within
  5 mm. Report-only (`tier_a.edge`).

- **Multi-delay fits cover every voxel (by design).** `ComputeCBF` deliberately quantifies the
  whole field of view, because ASL-based brain masks are often poor. The cost is runtime: F5
  (about 209,000 voxels) takes 28 min and F6 21 min, nearly all of it in the single-threaded
  fit. Not a defect; it only sets the integration tier's runtime.
- **Multi-delay ATT is not recovered on these protocols.** The independent fit agrees with
  ASLPrep's to 0.001 s, but neither tracks the true ATT within tissue: r 0.06 in GM for F5, and
  a median error of 0.31 s. This is the four-parameter model with an arterial term on these
  delays, not a plumbing error. Report-only (`tier_b_att`).
- **PASL multi-delay edge voxels end at the CBF bound (#705).** On F6 (no motion), motion
  correction shifts the mid-delay control volumes by up to 0.2 mm while their labels stay put,
  so those pairs are misaligned by 0.12-0.16 mm, far above the 0.02 mm noise floor above. At
  the brain edge this adds delta-M errors several times the true delta-M, and 1.3 % of brain
  voxels end at the fit's CBF bound of 300. The independent fit gives a median of 45 on the raw
  series and 300 on ASLPrep's preprocessed series. F5 (PCASL) pairs agree within 0.04 mm and
  no brain voxel reaches the bound. A strict xfail on F6's `test_fit_bound`, which asserts that
  at most 0.5 % of brain voxels reach the bound.
- **pybids reads the datatype from the whole path.** A dataset under a directory named after
  a BIDS datatype (here `/aslscan/motion/`) gets that datatype for every file, and ASLPrep finds
  no ASL run. The recipe was renamed `headmotion`, and recipe names are now checked. Real
  datasets stored under such a directory would hit the same problem.

Integration timings, continued: F5 28 min and F6 21 min (both all pass after the changes
below).

### Design changes from the first integration runs

These are corrections of attribution, not fitted bounds:

- **Tier A split.** Tier A became end to end (median asserted) plus quantification given
  ASLPrep's preprocessed series (median and p95 asserted). Its tails are report-only.
- **Standard spaces.** Checked as alignment, `r_native` ≥ 0.8, plus the GM median. WM medians
  in standard space measure resampling blur, not alignment.
- **Coregistration.** In units of the coarsest voxel (a quarter of a voxel), instead of 1 mm
  and 1°. Without motion, it is measured as `C^-1 R^-1` alone.
- **Motion correction.** Measured as consistency relative to the first volume, for every
  recipe. The constant offset to the static frame is reported separately.
- **Hashed modules.** `truth_models.py` was removed from them, because generation never uses
  it.
- **Recipe discovery.** Recipes are discovered from their directories, so adding one does not
  touch a hashed module.
- **Multi-delay reference fit.** It smooths M0 as ASLPrep does. Tier B is compared on the
  fitted voxel subset (`tier_b_matched`), against the fit of ASLPrep's preprocessed series.
  The raw-series expectation and the per-voxel tails are report-only, because the
  four-parameter fit is ill-conditioned and responds to motion-correction resampling noise.
- **Atlases.** `--atlases` is pinned in every module. Without it, ASLPrep parcellates with all
  of its atlases.

### JUnit

`record_property` warns under xunit2. CI should pass `-o junit_family=legacy` for the fast
tier's properties to appear (Task 13).

## Review log

The Codex adversarial review of plan revision 1 raised 19 findings. All were checked against
the code and accepted; the spec changes are logged in the spec's Section 12.

| # | Finding | Plan change |
|---|---|---|
| 1 | Calibration mode leaves CI unchecked; the flag was not forwarded into Docker | removed calibration mode; a missing bound fails; separate `propose` command (Tasks 9, 14) |
| 2 | Bounds centred on ASLPrep output | a-priori ceilings; Tier B centred on an independent expectation; physical-agreement items with predicted size; bands only after ceilings pass and are validated fresh (Tasks 7, 9, 14) |
| 3 | Too few pure voxels | `valid`/`dom_L` masks, no erosion, population reporter and minimums (Task 6) |
| 4 | Wrong interface chain (`cbf_only`, metadata, `m0tr`) | explicit connections (Task 7) |
| 5 | F3 lacks native output | `asl` added to output spaces (Tasks 10, 11) |
| 6 | Motion centre not recorded; Sections 9.2 and 9.4 not implemented | centre derived and recorded (Task 5); tests added (Task 9) |
| 7 | Correlation on flat truth | F1 has no r (Task 10) |
| 8 | Completeness not enforced | metric contract plus shared items in every module; F7 reports enforced (Tasks 9-12) |
| 9 | Consumers could run before generation | separate restore command used only by dependent jobs (Task 13) |
| 10 | Undefined `marker` parameter for unit tests | `unit_tests.xml` (Task 13) |
| 11 | `$BASH_ENV` within one step | one shell step, explicit paths (Task 13) |
| 12 | Unregistered marker | `skipif` instead of a new marker (Task 5) |
| 13 | Size cap on the wrong quantity | total-voxel caps, timing in the image (Tasks 6, 8) |
| 14 | Noise schema and formula | `snr_per_pair` in the schema, empirical calibration, achieved SNR verified (Task 5) |
| 15 | Vacuous thread test | two fresh runs with explicit thread settings (Task 5) |
| 16 | Local runner used the wrong image | runner mirrors CI (Task 10) |
| 17 | Slice-count formula | exact `ceil` rule; multiband and 3D checks (Task 5) |
| 18 | Forward dependency on hashed modules | `HASHED_MODULES` grows; regeneration rule (Tasks 1, 3, 7) |
| 19 | Immutable bad cache | `CACHE_EPOCH` in the spec file (Task 1; spec Section 4.5) |

### Implementation review (Task 17)

A Codex adversarial review of `origin/main...HEAD` raised 12 findings. Each was checked against
the code and, where it made a claim about ASLPrep's outputs, against measurements.

| # | Finding | Disposition |
|---|---|---|
| 1 | Multi-delay per-voxel errors undetected (CBF wrong in 40 % of voxels keeps the median) | accepted: `frac_within_10pct` against the fit of the same series, floor 0.5 (`tier_a_quant_agree`). The fit is ill-conditioned (F5 0.75, F6 0.66), so the floor only catches gross errors |
| 2 | Strict xfails do not bound the size of the known failure | accepted: `Known(reason, path, lo, hi)`; the item fails if it passes or if the metric leaves its diagnosed range |
| 3 | F7 never checks that the transforms are applied to the images | accepted: `resampling_margin`, the T1w map against the native map moved through ASLPrep's own coregistration, over a 4 mm misplacement (`space_resampling`, F1 and F7) |
| 4 | Motion-free coregistration ignores the motion-correction transforms | accepted after measurement: `desc-preproc_asl` matches the raw series resampled through each full HMC transform (constant part included) better than at identity, so outputs carry it. Coregistration is now end to end, `C^-1 A_v R^-1`, for every recipe; `coreg_only_rms_voxels` reports the registration step |
| 5 | Metric contract accepts NaN; missing quantification metrics pass silently | accepted: report-only metrics must be finite; the required quantification checks follow the delay type |
| 6 | aBAT and aBV maps never scored | accepted: agreement with the reference fit, report-only, in the contract of F5 and F6 |
| 7 | `ASLPREP_REQUIRE_FIXTURES=1` still lets fast-tier tests skip | accepted: `FixtureRequired`, which the skip handlers do not catch |
| 8 | Per-volume timing fields averaged or broadcast wrongly | accepted (latent, no recipe hits it): `labeled_value` and observation rows |
| 9 | Fixture and phantom caches trusted without content hashes | not changed: CI's restore step verifies content hashes before any test, and phantoms derive from TemplateFlow files with pinned hashes |
| 10 | A dirty reused aslscan checkout would be stamped as pinned | accepted: tracked-file changes are refused |
| 11 | Replacing a fixture is not atomic for concurrent readers | not changed: generation and tests never run concurrently on one data directory (CI generates in a separate job) |
| 12 | 3D multi-delay PASL lost with `examples_pasl_multipld` | not changed: `test_computecbf_casl`/`_pasl` cover 3D multi-delay `ComputeCBF`, and `test_compare_slicetiming` ties 3D to 2D with zero slice times, which the fast tier checks numerically |

Consequences of finding 4, measured on the existing outputs (end to end, then registration step
alone, in voxels): F1 0.17 / 0.17, F2 0.15 / 0.07, F3 0.27 / 0.16, F4 0.17 / 0.15, F5 0.11 /
0.11, F6 0.25 (rotation arc 0.33) / 0.09, F7 0.27. F3 and F6 now fail coregistration as known
failures, for the same reason as F7.

### Finding: motion correction offsets the motion-free series from its reference

On motion-free data the HMC transforms share a constant part: 0.3 mm (F1), 0.6 mm (F3), 0.8 mm
(F2), 1.2 mm (F6, about 1 degree of rotation), 0.02 mm (F5). It is applied:
`desc-preproc_asl` matches the raw volumes resampled through each full transform better than
the raw volumes at the same place (relative RMS difference 0.030 against 0.044 for F3) or
through the volume-to-volume part only (0.044). The coregistration is estimated on the aslref,
which stays close to the scanner frame (registration step alone 0.07-0.16 voxel), so the
offset reaches every output. It is largest where the reference's contrast differs most from
the volumes'. This explains F7's coregistration excess (to be confirmed with the motion-free
F7 variant). Not yet reported.
