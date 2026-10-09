# Ground-truth tests with aslscan-simulated data: design

Status: revision 2 (2026-10-08). This revision addresses the Codex adversarial review (Section 12)
and adds pluggable phantom sources so real-subject maps can be used later (Section 4.2).
Branch: `simulated-test-data`.

## 1. Goal

ASLPrep's integration tests only check that a run finishes and writes an expected list of
files. A run that produces CBF maps off by a factor of two, swaps label and control, or ignores
slice timing passes. This design adds a second kind of test: ASLPrep processes data simulated
by [aslscan](https://github.com/PennLINC/aslscan), and the outputs are scored against the
simulator's ground truth.

It follows the pattern of PennLINC/qsiprep#1173 (TRXScan fixtures) and covers the
digital-reference-object level of PennLINC/aslprep#679.

Goals:

1. Generate a catalogue of simulated BIDS datasets in CI, reproducibly and cheaply, from a
   pinned aslscan revision, pinned phantom inputs, and protocol files committed to this
   repository.
2. Score ASLPrep's quantitative outputs (CBF, ATT, motion, coregistration) against the truth,
   with tolerances given per metric and unit. Failures must be readable in the CircleCI Tests
   tab.
3. Run fast, preprocessing-free conformance checks of the CBF interfaces on every push, in the
   unit-test job.
4. Replace the five `examples_*` integration tests (bids-examples data that only check file
   names) once each of their tested branches has an explicit simulated replacement (Section 6).
5. Make the phantom a pluggable input, so a later change can simulate from a real subject's
   CBF, ATT, aBAT, aBV, T1, T2 and T2* maps without redesign.

Non-goals (this change):

- Scalar, equation-level reference tests against the OSIPI toolbox (#679, "unit-level").
- Public benchmark data (OSIPI challenge) and scheduled or release-only runs.
- Susceptibility distortion correction scoring.
- Building a real-subject phantom. The interface is defined here; the data and its hosting are a
  follow-up.
- Removing the real-data tests (`test_001`, `test_002`, `test_003_*`, `qtab`).
- Publishing aslscan binaries or wheels.

## 2. Facts this design relies on

These were measured with aslscan `82919454d11ee63a1bc8032e4a362f9b745e5d97` and mrsim-acq
`e64caff8b18358704ef876b5a16a344192b3a014`, both the current `main` of their public PennLINC
repositories.

**Build.**

- aslscan has no releases. The binary needs the `cli` feature; the default feature set is
  empty.
- The build command is
  `cargo build --release --locked --features cli,kspace,par --bin aslscan`.
- mrsim-acq is a path dependency at `../mrsim-acq`.
- A release build takes about 20 s locally once crates are fetched.

**Cost.** One full-brain simulation takes 5-6 s and 0.8 GB for 2D EPI PCASL (20 volumes, about
64x68x39), 9 s for 3D GRASE, and 23 s for 3D spiral. A dataset is about 20 MB.

**Outputs.**

- `sub-01/perf/sub-01_part-{mag,phase}_asl.nii.gz` with sidecars, plus `sub-01_aslcontext.tsv`.
- `sub-01_m0scan.nii.gz` when `M0Type` is `Separate`.
- `sub-01/perf/ground-truth/` holds `desc-{perfusion,att,T1map,T2map,M0map,dseg,deltam}_gt`,
  each with a JSON giving its `Units` and `Resampling`.
- Motion adds `desc-motion_gt.tsv` (translations in mm, rotations in radians, rotation
  `Rz Ry Rx` about the field-of-view centre; `mrsim-acq/src/motion.rs`) and
  `desc-motionEvents_gt.tsv`.
- `.bidsignore` hides `ground-truth/`.

**Geometry.**

- The ASL series, the M0 scan and the truth maps share one RAS grid, placed in the phantom's
  world coordinates.
- The slice count is set by the phantom's extent divided by the slice thickness. 193 mm at
  5 mm gives 39 slices, and `SliceTiming` must have exactly that many entries
  (`resample.rs:295`, `series.rs:597`).

**Truth maps.**

- `perfusion_gt` is a volume-weighted mean over the phantom voxels each acquisition voxel
  overlaps (exact cell overlap, `resample.rs:119`).
- `att_gt` is the same, over perfused phantom voxels only.
- `dseg_gt` is a majority vote.
- The acquired images also carry the acquisition's spatial response (oversampling, Fourier
  truncation), which the truth maps do not.
- `desc-deltam_gt` is populated on label and `deltam` rows and is zero on control rows. It is in
  simulator units before acquisition scaling: the image difference divided by it is about 87-94,
  varying with tissue. Scoring never compares image intensities with it directly.

**Phantom contract** (`src/phantom.rs`):

- Required: one NIfTI and JSON per map, with `Units` exactly as listed:
  - `perfusion` (`ml/100g/min`), `att` (`s`), `T1map` (`s`), `T2map` (`s`), `T2starmap` (`s`);
  - `M0map` (`arbitrary`);
  - `dseg` (`label indices`, plus a `LabelMap`).
- Optional: `abv` (`fraction`), `aatt` (`s`), `fieldmap` (`Hz`), and `phantom.json`
  (`LambdaBloodBrain`, `T1ArterialBlood`, `MagneticFieldStrength`).
- Background M0 must be 0. Foreground T1, T2, T2* and M0 must be positive, and T2* < T2,
  because T2' is derived from them.
- `class` T2 mode requires T2 and T2* to be constant within each label. Perfusion, ATT and T1
  may vary voxelwise.
- The 3D gradient-echo EPI readout refuses smoothly varying T1.

**Model mismatch.**

- aslscan uses the Buxton general kinetic model: bolus integration, a flow-dependent apparent
  tissue T1, and label in tissue relaxing with tissue T1.
- ASLPrep's single-delay equation (Alsop 2015) assumes blood T1 throughout.
- On noise-free 2D PCASL (PLD 1.8 s), the white-paper equation gives about 0.78x the true CBF in
  pure GM and about 0.33-0.45x in WM, the same on both phantoms tried.
- This is expected physics, not an ASLPrep bug, and it is why scoring has two tiers (Section 7).

**Background suppression.**

- ASLPrep's default efficiency, when the sidecar has no `LabelingEfficiency`, is
  `base x 0.95^n`, with n = `BackgroundSuppressionNumberPulses`, default 1 (`utils/asl.py`).
- aslscan's signed BS factor on the label is `(1 - 2 x 0.95)^n` at its default inversion
  efficiency: magnitude 0.9 per pulse, with the sign alternating (`longitudinal.rs:254`).
- These disagree physically, so Tier B scores and documents the gap (Section 8).

**ASLPrep input behaviour.**

- ASL runs are queried as `datatype=perf, suffix=asl`, so a `part-phase` series would become a
  second run.
- M0 scans are found through BIDS `InformedBy` associations, so the m0scan sidecar needs
  `IntendedFor`.
- `ComputeCBF` builds an all-ones mask and fits every voxel; it has no mask input.
- `ExtractCBF` smooths M0 with `fwhm` (default 5 mm, from `--smooth_kernel`) through
  `nibabel.processing.smooth_image`; `fwhm = 0` disables smoothing.
- The M0 TR correction below 5 s uses ASLPrep's fixed tissue T1 of 1.607 s at 3 T.
- The anatomical fast-track (`--derivatives`) needs, under `sub-01/anat/`:
  - `desc-preproc_T1w`, `desc-brain_mask`, `dseg` and `label-{GM,WM,CSF}_probseg`;
  - `.h5` or `.txt` `from-T1w_to-<template>` transforms and their reverses.

  Templates without transforms are registered by ASLPrep.
- Atlases are warped through TemplateFlow's `MNI152NLin2009cAsym` to `MNI152NLin6Asym`
  transform, not by a new registration.

## 3. Overview

```
.circleci/aslscan_fixtures.txt (generated: pins + every input digest; the cache key)
        |
gen_aslscan_fixtures job (cimg/python:3.12 + pinned rustup toolchain)
  1. spec-sync guard (always, before any cache decision)
  2. restore cache; accept only if every fixture's manifest matches the spec digest
  3. on miss: build aslscan at the pins; python -m aslprep.tests.aslscan_fixtures
       --generate /tmp/data --aslscan <path>
       a. phantoms (Section 4.2), verified against pinned template checksums
       b. per recipe: aslscan -> assemble BIDS dataset + truth -> manifest
  4. save_cache
        |
restore_build_and_data (every test job) restores the same cache
        |
unit_tests: fast conformance tier on small cropped fixtures + dataset interpretation
aslscan_* integration modules: one ASLPrep run per module (module-scoped fixture),
  score JSON written, then one pytest item per metric
        |
truth_score.json (artifact) + JUnit XML (Tests tab)
```

## 4. Components

### 4.1 Simulator pinning and build

The registry module (4.5) pins:

- `ASLSCAN_REF` and `MRSIM_ACQ_REF`: full SHAs;
- `RUST_TOOLCHAIN`, an exact stable version, for example `1.89.0`.

`build_aslscan(workdir)`:

1. Clones both repositories side by side.
2. Checks out the SHAs and verifies `git rev-parse HEAD` matches.
3. Runs `cargo +$RUST_TOOLCHAIN build --release --locked --features cli,kspace,par --bin aslscan`.
4. Returns the absolute binary path.

aslscan does not commit a `Cargo.lock`: its `.gitignore` lists it. ASLPrep therefore vendors
one, `aslprep/tests/data/aslscan-Cargo.lock`, taken from the build that produced the measured
facts. It is copied into the checkout before the `--locked` build, and its digest is part of
the spec file and the binary stamp.

The path is passed explicitly to `--generate` with `--aslscan PATH`. There is no reliance on
`PATH` persisting across CI steps.

`aslscan --version` reports only `0.0.0`, so a binary cannot prove its revision.
`build_aslscan` therefore writes `aslscan.build.json` next to the binary, recording both SHAs,
the toolchain, the features and the binary's SHA-256.

Locally, `find_aslscan()` resolves the binary from the `ASLSCAN` environment variable, then
from `PATH`. A binary found this way is accepted only if its stamp exists, matches the pins,
and matches the binary's hash. Otherwise generation stops with a message to run
`--build-aslscan`.

`test-hooks` is never enabled. `RAYON_NUM_THREADS` is fixed (4) during generation. The builder
verifies byte identity of one recipe at 1 and 4 threads in a fast-tier test.

### 4.2 Phantoms (pluggable)

A phantom is selected by name in each recipe (`recipe.toml: phantom = "<name>"`). `PHANTOMS` in
the registry maps names to builders. Every phantom produces the same directory contract:

```
_phantoms/<name>/
  perfusion.nii.gz/.json  att  T1map  T2map  T2starmap  M0map  dseg   (required, aslscan contract)
  abv.nii.gz/.json  aatt.nii.gz/.json                                 (optional)
  phantom.json                                                        (LambdaBloodBrain, T1ArterialBlood, MagneticFieldStrength)
  anat/T1w.nii.gz                                                     (anatomical image in the phantom's world)
  anat/probseg-{GM,WM,CSF}.nii.gz                                     (priors for derivatives; never used for scoring)
  provenance.json                                                     (source, input checksums, builder digest)
```

This change implements three builders.

**`tfmni` (default).** Built from TemplateFlow `MNI152NLin2009cAsym` `res-01`:

- `T1w`, `desc-brain_mask` and `label-{GM,WM,CSF}_probseg`.
- Each file's SHA-256 is pinned in the registry and verified after download. A mismatch is an
  error.
- Labels: the argmax of the probsegs inside the mask. A voxel whose maximum probability is
  below 0.3 is background (0). Labels are 1 grey_matter, 2 white_matter, 3 csf (the
  `LabelMap` names aslscan's overlays use).
- T1, T2, T2* and M0 are per-label constants taken from ASLDRO's `hrgt_icbm_2009a_nls_3t`
  (GM 1.33/0.080/0.066/74.62, WM 0.83/0.110/0.053/64.72, CSF 3.0/0.300/0.200/68.05).
  This is a new phantom that shares ASLDRO's constants, not ASLDRO's anatomy.
- Perfusion and ATT vary smoothly within tissue, so correlation metrics measure more than GM/WM
  discrimination:
  - `perfusion = base x (1 + 0.2 m1(x))`, base GM 60 and WM 20;
  - `att = base x (1 + 0.15 m2(x))`, base GM 0.8 and WM 1.2;
  - `m1` and `m2` are fixed low-frequency fields in world mm (products of cosines with periods
    of 90-140 mm, range [-1, 1], defined in code). CSF is perfusion 0, ATT 1000.

**`tfmni-flat`.** The same with `m1 = m2 = 0`. It is used where a constant truth simplifies
diagnosis (F1).

**`crop:<base>:<x0:x1,y0:y1,z0:z1>`.** A voxel-index crop of another phantom, used by the fast
tier.

The crop keeps the affine origin shifted correctly, so world coordinates are unchanged.

**Real-subject phantom (future, interface only).** A builder named `subject:<id>` downloads a
pinned tarball (URL and SHA-256 in the registry, as in `data_versions.txt`) containing
perfusion, ATT, aBAT (`aatt`), aBV (`abv`), T1, T2, T2*, M0 and dseg maps, the subject's T1w,
and optional anatomical derivatives.

The builder sanitizes the maps to the aslscan contract, and records every change it makes in
`provenance.json`:

- set M0 to 0 outside dseg > 0;
- clip foreground T1, T2, T2* and M0 to a positive floor;
- clip T2* to at most 0.95 x T2;
- set ATT to 1000 where perfusion is 0.

It requires `--t2-mode voxel` (a recipe-level aslscan argument). If the maps are coarser than
the T1w, it resamples them to the T1w grid. Such a phantom lacks sub-voxel partial-volume
detail, which the provenance records.

Everything downstream, including the T1w as anatomical image, `anat_offset`, scoring against
truth maps and macrovascular truth, works unchanged. When `abv` and `aatt` are present,
multi-delay recipes gain `abv` and `abat` metrics (Section 7).

Nothing in scoring may assume tissue constants. Expected values always come from the truth maps
under the same masks.

### 4.3 Recipes and protocol inputs

Each recipe directory, `aslprep/tests/data/aslscan/<name>/`, holds these committed files:

- `asl.json`: the BIDS sidecar given to aslscan.
- `aslcontext.tsv`.
- `overlay.toml`: always sets `seed`, `[acquisition] noise_variance` and
  `[acquisition] matrix`, plus any feature keys.
- `recipe.toml`: the settings ASLPrep's assembly step owns. Keys and defaults:

  ```toml
  phantom = "tfmni"            # Section 4.2
  acq = "pcasl1pld"            # BIDS acq- label for output names
  anat = "derivatives"         # or "raw"
  anat_offset = [0, 0, 0, 0, 0, 0]  # tx ty tz (mm) rx ry rz (deg), Section 4.4
  keep_labeling_efficiency = true
  m0_divisor = 1.0             # M0 image divided by this; test passes --m0_scale=<same>
  snr_per_pair = 0             # 0 = noise-free; else the target delta-M SNR in GM-dominant voxels (Section 5)
  t2_mode = "auto"             # aslscan --t2-mode (the only aslscan argument a recipe sets)
  ```

  `anat` may also be `"none"`, for fast-tier recipes on cropped phantoms.

- Optional `motion.tsv` (trajectory).
- `SOURCE.md` for protocols derived from bids-examples: the upstream file, its commit, and every
  change.

Every recipe commits its exact geometry (voxel sizes, `matrix`, `SliceTiming` with the slice
count implied by the phantom). The builder pre-computes the slice count from the phantom extent
and the slice thickness. It refuses a recipe whose `SliceTiming` length differs, before calling
aslscan, with a message giving the expected count.

### 4.4 Fixture assembly

Assembly runs after aslscan writes `<tmp>/raw/`. All generation happens in
`<data>/aslscan/.tmp-<name>-<pid>/`, under a lock file `<data>/aslscan/.lock` (fcntl on POSIX).
The finished directory is published by atomic rename to `<data>/aslscan/<name>/`. A stale
target is first renamed aside and then removed.

The layout:

```
<name>/
  dataset_description.json      (raw; GeneratedBy: aslscan + both SHAs; phantom name)
  README
  .bidsignore                   (**/ground-truth, derivatives/**)
  sub-01/anat/sub-01_T1w.nii.gz  (+ .json)
  sub-01/perf/sub-01_acq-<acq>_asl.nii.gz/.json        (the part-mag series, renamed)
  sub-01/perf/sub-01_acq-<acq>_aslcontext.tsv
  sub-01/perf/sub-01_acq-<acq>_m0scan.nii.gz/.json     (Separate only; IntendedFor: perf/sub-01_acq-<acq>_asl.nii.gz)
  sub-01/perf/ground-truth/
      sub-01_desc-*_gt.nii.gz/.json                    (aslscan's, unchanged)
      sub-01_desc-pv{GM,WM,CSF}_gt.nii.gz              (exact overlap fractions, below)
      simulation.json                                  (the AslscanSimulation block)
      truth.json                                       (Section 4.6 frames; offsets; recipe; digests)
  derivatives/anat/                                    (anat = "derivatives" only)
      dataset_description.json                         (DatasetType: derivative)
      sub-01/anat/sub-01_desc-preproc_T1w.nii.gz
      sub-01/anat/sub-01_desc-brain_mask.nii.gz
      sub-01/anat/sub-01_dseg.nii.gz                   (phantom labels)
      sub-01/anat/sub-01_label-{GM,WM,CSF}_probseg.nii.gz
      sub-01/anat/sub-01_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.txt
      sub-01/anat/sub-01_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.txt
  manifest.json                                        (written last; Section 4.5)
```

- **Phase series.** Dropped.
- **`AslscanSimulation` block.** It stays in the ASL sidecar, where ASLPrep ignores unknown
  keys.
- **`LabelingEfficiency`.**
  - `keep_labeling_efficiency = false` removes `LabelingEfficiency` from the ASL sidecar when
    `Resolved.LabelingEfficiency.Source == "Default"`. This exercises ASLPrep's
    background-suppression branch. It does not make the physics agree (Section 8).
  - Recipes with BS exist in both forms (F2 and F3).
- **`m0_divisor`.** It divides the M0 image (float32) to emulate vendors that store a scaled M0;
  the test passes the matching `--m0_scale`.
- **`desc-pv*_gt`: exact overlap fractions.**
  - The phantom and the acquisition grids are both axis-aligned RAS with diagonal affines. The
    builder asserts this and refuses otherwise.
  - Along each axis it builds a 1D matrix `W[a, p]` of the length of phantom cell `p` inside
    acquisition cell `a`, divided by the acquisition cell length.
  - The fraction of label L is the separable product `W_x ⊗ W_y ⊗ W_z` applied to
    `onehot(dseg == L)`.
  - This is the same volume-weighted overlap aslscan uses for `perfusion_gt`
    (`resample.rs:119`). A fast-tier test checks it reproduces `perfusion_gt` from the phantom
    perfusion map to 1e-4.
- **Anatomical offset.**
  - Notation: `A` is the T1w's phantom-world affine, and `R` is the 4x4 RAS rigid matrix from
    `anat_offset` (rotation `Rz Ry Rx` about the world origin, then translation).
  - The image data are untouched. The written T1w affine is `A' = R @ A`, with qform and sform
    both set. The same affine is applied to every derivative image.
  - Point mapping: a phantom-world point `p` appears in the anatomy's world at `R p`.
- **Derivative transforms.**
  - ITK image transforms used for resampling are pull mappings: a transform file
    `from-X_to-Y` maps points of Y's space to X's space.
  - With MNI space equal to phantom space (`tfmni`), `from-T1w_to-MNI152NLin2009cAsym` maps an
    MNI (phantom) point `p` to the T1w point `R p`. The file stores `R`, and its reverse stores
    `R^-1`, each converted RAS to LPS by conjugation with `diag(-1, -1, 1, 1)` (a basis change,
    not an inverse).
  - A real-subject phantom has no MNI transform. Its derivatives omit these files, and ASLPrep
    registers.
  - The convention is verified by landmark tests that resample through `nitransforms` (Section
    9), not by round trips.

### 4.5 Registry, spec file, and fixture manifest

`aslprep/tests/aslscan_fixtures.py` imports with the standard library only. Heavy imports are
local to functions.

`spec_text()` contents, with LF line endings:

```
# generated by python -m aslprep.tests.aslscan_fixtures --spec; do not edit
CACHE_EPOCH=<n>                             (bump to discard a bad cache entry; CircleCI caches are immutable)
ASLSCAN_REF=<sha>
MRSIM_ACQ_REF=<sha>
RUST_TOOLCHAIN=<ver>
PY_REQUIREMENTS=numpy==..;nibabel==..;scipy==..;templateflow==..
TEMPLATE <relpath> <sha256>                 (one per pinned TemplateFlow file)
BUILDER <sha256 of aslscan_fixtures.py>
BUILDER <sha256 of truth_scoring.py's geometry helpers module, if separate>
PHANTOM <name> <digest of builder parameters>
RECIPE <name> <digest>                      (one per recipe)
DIGEST <sha256 of all lines above>
```

A recipe digest is the SHA-256 over the sorted entries `relpath \0 length \0 LF-normalized
bytes` of its directory.

**Committed spec file.**

- `.circleci/aslscan_fixtures.txt` is `spec_text()`.
- `.gitattributes` marks it `text eol=lf`, and protocol files likewise.
- Its checksum is the CI cache key.

**Manifest.** Each fixture's `manifest.json` records:

- its fixture digest: the SHA-256 of the common spec lines (epoch, pins, lock file, Python
  requirements, templates, builder modules), its phantom's digest, and its recipe's digest.
  Adding or editing one recipe therefore leaves the other local fixtures valid. CI still
  regenerates everything when the spec file changes, because the spec file keys the cache;
- the spec `DIGEST` (provenance only);
- the recipe and phantom digests;
- the aslscan binary's SHA-256 and reported version;
- a SHA-256 of every file in the fixture.

**`fixture_dir(name, data_dir)`.**

- It returns the directory only if `manifest.json` exists and its fixture digest equals the
  current one.
- File hashes are checked by `--verify`, which CI runs once after restoring the cache, not on
  every open, for cost.
- If the manifest is stale or missing:
  - when `ASLPREP_REQUIRE_FIXTURES=1` (set in CI), it raises;
  - otherwise it generates on demand under the lock, if an acceptable aslscan is found;
  - otherwise the fast tier skips and integration tests fail.

**Command line.**

- `--spec`
- `--generate DATA_DIR --aslscan PATH [--only NAME...]`
- `--build-aslscan WORKDIR` (prints the path)
- `--verify DATA_DIR`

### 4.6 Frames and transforms used in scoring

All transforms are 4x4 RAS point mappings unless stated otherwise:

- `P_v`: the simulator's pose for volume `v`. It maps a static-phantom point to its moved
  position. It is built from `desc-motion_gt.tsv` row `v` as `T(c) Rz Ry Rx T(-c) + t`.
  - aslscan does not record the centre `c`. Motion is applied about the simulation grid's
    centre (`aslscan/src/series.rs:1177`, `mrsim-acq/src/motion.rs:199`), so the builder
    derives `c` from the simulation grid in `simulation.json` and records it in `truth.json`.
  - The exact composition comes from `motion.rs` and is pinned by the simulator landmark test
    (Section 9.4).
  - Without motion, `P_v = I`.
- `H_v`: ASLPrep's HMC transform for volume `v`, from
  `from-orig_to-aslref_mode-image_desc-hmc_xfm.txt`. As stored it is a pull mapping (an aslref
  point to a volume-`v` point); the scorer converts ITK LPS to RAS.
- `C`: ASLPrep's `from-aslref_to-T1w_mode-image_desc-coreg_xfm.txt`, a pull mapping (a T1w
  point to an aslref point).
- `R`: the anatomical offset (Section 4.4). The true mapping from a static-phantom point to its
  T1w point is `R`.

Derived quantities:

- The estimated mapping from a static-phantom point (seen in volume `v`) to the T1w is
  `E_v = C^-1 H_v^-1 P_v`. Here `P_v` moves the point to where volume `v` saw it, `H_v^-1`
  carries volume-`v` coordinates to the aslref, and `C^-1` carries aslref coordinates to the T1w.
- The error is `D_v = E_v R^-1`, which is the identity when everything is right. The `coreg`
  metric is the median over volumes of the rotation angle of `D_v` and the RMS displacement of
  `D_v` over the truth brain voxels' centres (mm).
- The aslref's pose relative to the static phantom is `A_ref = H_v^-1 P_v`, which should agree
  across `v`; the spread is reported as `hmc_consistency`.
- Native-space ("asl") outputs are compared with native truth only when
  `max_v |A_ref - I|` is below 0.1 mm or 0.1°. That holds for every no-motion recipe; the scorer
  asserts it rather than assuming it.
- The motion recipe is scored in T1w space through `R`, plus the motion metrics.

### 4.7 Scoring (`aslprep/tests/truth_scoring.py`)

The scorer uses only numpy, scipy, nibabel and, for tests only, nitransforms. It returns a
nested dict and writes it before any assertion.

**Masks.** Ratios are taken per voxel against that voxel's own truth, which already includes
partial volume. Tissue purity is therefore needed only to interpret results per tissue, not to
make the comparison valid. Strict purity with erosion leaves too few voxels at ASL resolution
(about 190 GM voxels on F1's grid after one erosion, and none after a smoothing erosion).

- `brain_truth = (pvGM + pvWM + pvCSF) >= 0.5` on the acquisition grid.
- `valid = brain_truth & (pvCSF < 0.5) & (denominator > floor)`. The floor is 5 ml/100g/min on
  `perfusion_gt` for Tier B and on `cbf_A` for Tier A. This excludes zero-perfusion CSF and
  near-zero denominators.
- `dom_L = valid & (pv_L >= 0.7)` for the tissue-dominant masks (GM, WM).
- No erosion. Smoothing is applied identically to `cbf_A` (Tier A), so it needs no mask
  margin.
- Each metric's mask is intersected with ASLPrep's output brain mask, resampled by nearest
  neighbour where grids differ.
- Every statistic asserts its mask is non-empty and at least `n_min` voxels, before computing.
- The populations of `valid` and `dom_L` are measured for each recipe when the recipe is
  written. A recipe is accepted only with at least 300 voxels in each `dom_L`.

**Coverage guards.** Each is its own metric, asserted:

- `mask_coverage`: the fraction of `brain_truth` inside ASLPrep's mask; ≥ 0.95.
- `finite_fraction`: the fraction of finite output values in the scoring mask; ≥ 0.999 for
  single-delay and ≥ 0.98 for multi-delay. ASLPrep's fit can fail per voxel.
- `n_voxels[L]`: ≥ 300 per `dom_L`, so a shrinking mask cannot pass.
- Nonzero denominators; `cbf_A` must be > 0 wherever a ratio is taken.

**Tier A, conformance.** These checks catch plumbing bugs: ordering, sign, timing, scaling and
M0 handling. Tier A is scored on raw inputs, in two variants.

- **Fast tier: exact conditional conformance.** The interfaces are run on the raw inputs:
  `ExtractCBF` with `fwhm=0`, then `ComputeCBF`.
  - `ComputeCBF` receives `ExtractCBF`'s `out_file`, `m0_file`, volume-selected `metadata`,
    and `m0tr` (only when not None), with `cbf_only=False`. `cbf_only=True` means the data are
    already CBF and bypasses quantification.
  - The scorer's independent implementation computes `cbf_A` from the same raw images. It
    covers:
    - aslcontext parsing;
    - Alsop 2015 PCASL and PASL QUIPSS, QUIPSS II and Q2TIPS single-delay equations;
    - the per-slice PLD shift from `SliceTiming` and `SliceEncodingDirection`;
    - labeling efficiency per ASLPrep's documented rule (sidecar value, else base x
      0.95^n_pulses);
    - the M0 TR < 5 s correction with T1 tissue 1.607 s at 3 T;
    - `m0_divisor` and `--m0_scale`.
  - The written equations cite the papers, not ASLPrep's code.
  - Bounds: `|median(cbf/cbf_A) - 1| <= 1e-3` and p99 of `|cbf/cbf_A - 1| <= 1e-2`.
  - Each recipe also gets **mutation tests**. The scorer's own `cbf_A` is recomputed with one
    deliberate error at a time: label and control swapped, slice timing reversed, PLD shifted
    by 0.1 s, M0 scaled by 1.1, labeling efficiency ignored. The ASLPrep output must fail the
    bound for each, which proves the check is sensitive.
  - Smoothing is tested separately: `fwhm=5` versus `cbf_A` computed with
    `nibabel.processing.smooth_image(m0, 5)`, with the same bounds.
- **Integration: end-to-end regression.** ASLPrep's native-space CBF is compared with `cbf_A`,
  computed from the raw inputs with ASLPrep's configured smoothing applied to M0 by the same
  nibabel call.
  - Preprocessing also registers and resamples M0 and averages volumes robustly, so this is a
    regression measure, not exact conformance.
  - Acceptance ceilings are fixed in advance, not fitted: 5 % on the median deviation and 15 %
    on p95.
- **Multi-delay Tier A.** There is no closed form. The scorer refits all observations (every
  label/control or deltam volume, not per-PLD means, so unequal repeats weight correctly). The
  fit uses its own implementation of the documented forward models:
  - PCASL: Buxton with an arterial component, as described in ASLPrep's documentation and
    Woods 2023; PASL likewise.
  - The same four parameters, initial values `(60, 1.2, 1, 0.02)`, bounds
    `((0, 0, 0, 0), (300, 5, 5, 0.1))`, and failed-fit-to-NaN handling as ASLPrep, stated in
    the scorer with citations.

  The fast tier runs it on cropped fixtures of at most 2,000 voxels. It is accepted when the
  median CBF ratio is within 2 % and ATT within 0.05 s. The same mutation tests apply.

  Its independence is limited to the forward equations. Tier B is what checks the physics.

**Tier B, accuracy.** These metrics compare against the physical truth, voxel by voxel. The
reference is always the truth map itself, never a tissue constant.

**Expected bias.** The scorer computes `expected_ratio[L] = median(cbf_A / perfusion_gt)` over
`dom_L` from the raw data, with no ASLPrep involvement. For multi-delay recipes it uses the
scorer's own fit. This is what an implementation of ASLPrep's documented model should show,
given the simulator's physics.

**Acceptance.** ASLPrep's `ratio_median[L]` must lie within the Tier A ceiling of
`expected_ratio[L]`. The bound is centred on an independent expectation, never on ASLPrep's
own output.

**Physical-agreement items.** These are a separate set of checks. ASLPrep's output is compared
with `cbf_phys`: `cbf_A` recomputed with the simulator's true labeling attenuation (the
sidecar's `AslscanSimulation.BackgroundSuppressionLabelFactor`) instead of ASLPrep's convention.
A known disagreement is a strict xfail with its a-priori predicted size (Section 8).

- `ratio_median[L]`: the median of `cbf / perfusion_gt` over `dom_L`.
- `ratio_p05`, `ratio_p95`: tails, to catch local failures.
- `r[L]`: Pearson r within each tissue, only on `tfmni` recipes, which vary perfusion within
  tissue. It is not computed for `tfmni-flat` (F1), where the only variation is partial-volume
  contamination.
- `gm_wm_contrast`: `median(cbf[pure_GM]) / median(cbf[pure_WM])` divided by the same ratio of
  the truth under the same masks.
- For ATT, outside CSF: `att_mae[L]` (s), `att_r[L]`.
- With macrovascular truth: `abv` and `abat` medians and correlations, report-only until
  measured.
- `contamination[L]`: the median of `perfusion_gt` over `dom_L` divided by the pure-tissue
  value from the phantom. It measures what the masks let through. It is report-only, used to
  interpret Tier B.

**Other metrics.**

- `coreg` and `hmc_consistency` (Section 4.6).
- `motion`. The primary metric compares relative transforms: the per-volume RMS displacement
  over brain voxels of `H_v^-1 P_v (H_0^-1 P_0)^-1` (mm), median and max.
- A secondary metric, report-only at first, compares ASLPrep's `trans_*`/`rot_*` confounds with
  the truth converted through the same FSL-parameter convention (fMRIPrep `FSLMotionParams`).
  Its conversion is validated by unit tests on known translations, combined rotations, a
  non-central origin and an oblique grid.
- `fd_mean` is compared with FD (Power) computed from the truth parameters converted the same
  way.
- `cbf_t1w`: `space-T1w` CBF sampled at acquisition-grid truth voxel centres mapped through `R`;
  Tier B metrics.
- `cbf_mni`: `space-MNI152NLin2009cAsym` CBF on `tfmni` (MNI equals phantom up to `R`); Tier B
  metrics. Skipped for real-subject phantoms.
- `basil`: Tier B on `desc-basil_cbf` (and the GM/WM variants where written). Report-only until
  measured.
- `scorescrub`: Tier B on `desc-score_cbf` and `desc-scrub_cbf`. On no-motion recipes they must
  be within 5 % of the mean CBF's Tier B ratio. On the motion recipe, report-only.
- `clean`: no `sub-*/log/*/crash-*.txt` file, and the report contains "No errors to report!".

**Metric contract.** Each recipe declares the metrics it must produce.

- Every integration module has a `test_metric_contract` item. It fails if any declared metric,
  including report-only ones, is missing or non-finite.
- Every module also has the shared invariants as items: completeness, coverage, finite values,
  voxel counts, the clean run, and the output manifest.

**Tolerances.** Bounds are of two kinds, in each metric's own unit (ratio, s, mm, deg,
unitless r):

1. **Acceptance ceilings.** Fixed a priori in the spec or code, from the method's documented
   expectations; never fitted to ASLPrep output:
   - Tier A: 5 % on the median and 15 % on p95 in integration; 0.1 % and 1 % in the fast tier.
   - Tier B: within the Tier A ceiling of the independent expected ratio.
   - Coregistration: ≤ 1.0° and ≤ 1.0 mm RMS (with `--sloppy`, as qsiprep).
   - HMC: median ≤ 0.5 mm and max ≤ 1.5 mm RMS displacement error.
   - ATT: ≤ 0.15 s median absolute error from the expected fit.
   - Coverage and finite-value thresholds as above.

   A ceiling may only be loosened in a commit that explains why the method cannot meet it.
2. **Regression bands.** Optional, separate `test_regression_*` items, added only after:
   - the acceptance ceilings pass;
   - at least five measured runs (local and CI) exist;
   - the frozen band then passes on a fresh run.

   A band is measured ± max(3 × spread, half the ceiling). Every band carries a comment with
   the measured values, the date and the aslscan SHA. A unit test asserts that every ratio
   band is at most 0.1 wide, so a 10 % scaling error always fails.

There is no calibration mode in the tests. A missing bound is always a failure. Proposed bands
come from a separate developer command that reads `truth_score.json` files.

### 4.8 Tests

**Fast tier** (`unit_tests`):

- `test_aslscan_fixtures.py`:
  - spec sync;
  - recipe validation (files, JSON, aslcontext versus array lengths, the slice-count rule,
    recipe keys);
  - the exact-overlap check against `perfusion_gt`;
  - geometry landmark tests (Section 9);
  - thread-count identity.
- `test_aslscan_conformance.py`: Tier A (exact and smoothing variants), the mutation tests, and
  Tier B reported but not asserted. It is parametrized over the fast recipes (Section 5),
  which are small cropped fixtures, because `ComputeCBF` has no mask input.
- `test_aslscan_interpretation.py`:
  - `collect_data` finds exactly one ASL run;
  - `collect_run_data` resolves the aslcontext and the M0 (through `IntendedFor`);
  - metadata `M0Type` handling.
- Budget: measured before merge; the target is ≤ 3 min added.

**Integration tier.** One module per recipe, `test_aslscan_<recipe>.py`, with
`pytestmark = pytest.mark.aslscan_<recipe>`.

- A module-scoped fixture runs ASLPrep once (CLI arguments from a shared helper), writes
  `truth_score.json`, and returns the score.
- Each metric is its own test function (`test_tier_a_gm`, `test_coreg`, ...), so one failure
  does not hide the others.
- Known gaps use `@pytest.mark.xfail(strict=True, reason=...)` on the single affected item.
- `test_outputs_manifest` keeps the expected-outputs check, because output names are part of
  the BIDS contract.
- `test_clean_run` covers the clean-run check.

### 4.9 CI (`.circleci/continue_config.yml`)

**`gen_aslscan_fixtures`.** Docker `cimg/python:3.12`, `resource_class: medium`. Steps:

1. Checkout.
2. Spec-sync guard: `python -m aslprep.tests.aslscan_fixtures --spec | diff - .circleci/aslscan_fixtures.txt`.
   It runs before any cache decision.
3. `restore_cache` on `aslscan-v1-{{ checksum ".circleci/aslscan_fixtures.txt" }}`.
4. If every recipe's manifest digest matches the spec, stop. Otherwise remove
   `/tmp/data/aslscan`.
5. Install the toolchain:
   `curl ... rustup-init -y --profile minimal --default-toolchain $RUST_TOOLCHAIN` (retried
   three times). Then `echo 'source $HOME/.cargo/env' >> $BASH_ENV`.
6. Create a venv at `/tmp/venv`, with `pip install` of the pinned requirements (retried), and
   `echo 'source /tmp/venv/bin/activate' >> $BASH_ENV`.
7. `--build-aslscan /tmp/aslscan-build`, capturing the path, then
   `--generate /tmp/data --aslscan <path>`. Network failures, including GitHub and TemplateFlow,
   fail the job after the retries; no partial cache is saved.
8. `--verify /tmp/data`.
9. `save_cache` (only reached on success).

Shell state does not carry between CircleCI steps, and `$BASH_ENV` additions take effect only
in later steps. Each step therefore calls `/tmp/venv/bin/python` and `$HOME/.cargo/bin/cargo`
by absolute path, and passes the binary path through a file (`/tmp/aslscan-build/BIN`).

**Other jobs.**

- A separate command, `restore_aslscan_fixtures`, restores the cache, runs `--verify` with
  `python3`, and fails on a miss.
- It is used only by jobs that require `gen_aslscan_fixtures`: `unit_tests` and every
  `aslscan_*` job. Legacy integration jobs neither restore nor verify fixtures.
- The fixture-consuming containers get `-e ASLPREP_REQUIRE_FIXTURES=1`.
- JUnit output goes to `--junitxml=/test-results/<< parameters.marker >>.xml` for integration
  jobs and to `/test-results/unit_tests.xml` for `unit_tests`, followed by `store_test_results`.
- The expanded configuration is checked with `circleci config process` (or the CircleCI
  config-validation API) before pushing.
- Integration jobs (`aslscan_*`) are added to the matrix. `[skip aslscan]` skips them all.
- `merge_coverage` requires the new jobs.
- `get_data` and `data_versions.txt` are unchanged.
- The `examples_*` jobs are removed only in the final task, after the equivalence map
  (Section 6) is satisfied in CI.

### 4.10 Local workflow

- Generation runs on the host:
  `python -m aslprep.tests.aslscan_fixtures --build-aslscan ~/aslscan-build`, then
  `--generate aslprep/tests/test_data --aslscan <path>`.
- `run_local_tests.py` is updated to mirror CI:
  - it uses the `pennlinc/aslprep:test` image (built locally with `--target test`);
  - it mounts the checkout read-only at `/tmp/src/aslprep` as the working directory;
  - it mounts `aslprep/tests/test_data` as `/data`;
  - it forwards `ASLPREP_REQUIRE_FIXTURES` and `CIRCLE_CPUS`.

  A check in the runner prints `aslprep.__file__` from inside the container, to confirm the
  mounted checkout is the code being tested.
- The container never builds aslscan.
- The developer page (`docs/developers.rst`, new, linked from the toctree) documents
  generation, scoring, the tolerance policy, and the bump procedure: change the pins,
  regenerate the spec, re-measure, and move bounds in the same PR.

## 5. Fixture catalogue

**Integration recipes** (full-brain `tfmni` unless noted; geometry committed per recipe). Every
recipe scored in native space adds `asl` to `--output-spaces`, because native derivatives are
written only when that space is requested (`workflows/asl/base.py:478`). A blank "Key flags"
cell means `--output-spaces asl` only.

| ID / marker | Protocol | Exercises | Anat | Key flags | Scored |
|---|---|---|---|---|---|
| F1 `aslscan_pcasl1pld` (`tfmni-flat`) | 2D EPI PCASL, LD 1.8, PLD 1.8, 3.5x3.5x5 mm, 39 slices ascending 30 ms, TR 5 s, 10 pairs control-first, separate M0 TR 8 s, no BS | single-delay PCASL, slice-time PLD, SCORE/SCRUB, BASIL, T1w and MNI outputs, atlases | derivatives (identity) | `--output-spaces asl T1w MNI152NLin2009cAsym`, `--scorescrub --basil`, atlases 4S156Parcels | A, B, coreg, cbf_t1w, cbf_mni, basil, scorescrub |
| F2 `aslscan_pcasl1pld_ge3d` | asl001-derived 3D spiral, `deltam` + included `m0scan`, BS, `LabelingEfficiency` removed, `m0_divisor` 96 | deltam path, included M0, BS branch, m0_scale | derivatives | `--m0_scale=96 --basil` | A, B (BS gap: strict xfail if measured), basil |
| F3 `aslscan_pcasl1pld_grase` | asl005-derived 3D GRASE, separate M0 with TR 3 s (short-TR correction), BS, `LabelingEfficiency` kept, `m0_divisor` 10 | 3D PCASL, M0 TR < 5 s correction, m0_scale, MNI output, atlases | derivatives | `--output-spaces MNI152NLin2009cAsym --basil --m0_scale=10`, atlases | A, B, cbf_mni |
| F4 `aslscan_pasl1pld` | 2D EPI PASL Q2TIPS (0.7, 1.6), PLD 1.8, separate M0, interleaved slice order | single-delay PASL Q2TIPS, interleaved timing | derivatives | | A, B |
| F5 `aslscan_pcasl_multipld` | asl004-derived 2D EPI, 6 PLDs, label-first, unequal repeats (8/8/8/6/6/6 pairs), separate M0, `m0_divisor` 10, arterial compartment (per-label aBV/aBAT) | multi-delay PCASL, label-first ordering, unequal repeats, abv/abat | derivatives | `--m0_scale=10 --scorescrub`, atlases | A (fit), B (CBF, ATT; abv/abat report) |
| F6 `aslscan_pasl_multipld` | asl003-derived 3D GRASE Q2TIPS, TIs 0.9-3.0 s, separate M0, `m0_divisor` 10 | multi-delay PASL | derivatives | `--m0_scale=10 --basil --scorescrub`, atlases | A (fit), B |
| F7 `aslscan_motion` | F1 protocol with multiband 3 (13 slice groups), committed trajectory (≤ 2 mm, ≤ 2°, spikes at three volumes), noise at GM ΔM SNR 5 per pair, T1w `anat_offset` [3, -4, 2, 5, -3, 4] | HMC, FD, coregistration under a known pose, smriprep from raw T1w, multiband | raw | `--output-spaces T1w --scorescrub` | motion, coreg, hmc_consistency, cbf_t1w, scorescrub (report) |

**Implemented deviations from the table** (recorded in each recipe's `SOURCE.md`):

- **F4** uses QUIPSS II, not Q2TIPS. Single-delay Q2TIPS has a known quantification
  discrepancy (ASLPrep uses TI2 in the decay term), which is pinned in the fast tier and goes
  to a separate fix.
- **F5** uses a total readout of 25 ms (asl004: 60 ms, which cannot precede its 14 ms echo
  time on a 68-line readout), slices 30 ms apart, and TR 4.5 s.
- **F6** is 2D EPI at 4x4x6 mm without background suppression. asl003's 3D GRASE geometry does
  not fit the full phantom with its `NumberShots`, 3D GRASE is covered by F3, and suppression
  is covered by the PCASL recipes.
- **F7** scores motion, coregistration, T1w-space alignment and coverage. Native-space Tier A
  and B are report-only there, because the motion-corrected series sits in the aslref frame.
- **Recipe names** drop the `aslscan_` prefix (for example `pcasl1pld`); the pytest markers keep
  it.

**Fast recipes** (`crop:tfmni` slabs of about 40x40x30 mm, small matrices; each about 1 s):

- `fast_pcasl_seq`, `fast_pcasl_rev` (descending slices, `SliceEncodingDirection k-`),
  `fast_pcasl_mb` (multiband 2, no motion);
- `fast_pasl_quipss2` (QUIPSS II), `fast_pasl_q2tips`;
- `fast_m0_included` (2D control/label + m0scan rows), `fast_m0_estimate` (`M0Estimate`),
  `fast_m0_absent` (control-derived);
- `fast_m0_shorttr` (separate M0 at TR 3 s);
- `fast_label_first`, `fast_deltam`;
- `fast_bs_le_absent` (BS 2 pulses, no `LabelingEfficiency`);
- `fast_pcasl_multipld`, `fast_pasl_multipld` (≤ 2,000 voxels).

**Noise.** It is specified as `snr_per_pair`, the target ΔM SNR in GM-dominant voxels per
ΔM estimate: one control-label difference, or one `deltam` row. The builder calibrates
aslscan's `noise_variance` (the variance per real and imaginary component of the reconstructed
complex image) empirically, so the same procedure works for paired and direct-ΔM acquisitions:

1. Run noise-free.
2. Run at a reference variance σ0².
3. Measure the standard deviation of the noisy-minus-noise-free ΔM estimates and the
   noise-free ΔM, in `dom_GM`.
4. Scale to σ² = σ0² × (ΔM / (SNR × s0))².
5. Run again and verify the achieved SNR is within 10 % of the target.

The target, σ² and the achieved SNR are recorded in `truth.json`.

- The fast recipes are noise-free.
- F1-F6 use SNR 50 per pair, so Tier A is not noise-limited and SCORE/SCRUB have finite
  variance.
- F7 uses SNR 5 per pair.

## 6. Replacement policy and equivalence map

| Old test | Branches it covered | Replacement |
|---|---|---|
| `examples_pasl_multipld` | multi-delay PASL Q2TIPS, BASIL, SCORE/SCRUB, `--m0_scale=10`, atlases, asl output | F6 |
| `examples_pcasl_multipld` | multi-delay PCASL, label-first, SCORE/SCRUB, `--m0_scale=10`, atlases | F5 |
| `examples_pcasl_singlepld_ge` | GE deltam + M0 included, BASIL, SCORE/SCRUB, `--m0_scale=96`, atlases | F2 (add `--scorescrub`, atlases) |
| `examples_pcasl_singlepld_philips` | 2D PCASL Philips, BASIL, SCORE/SCRUB, MNI output, atlases | F1 |
| `examples_pcasl_singlepld_siemens` | 3D GRASE PCASL, BASIL, `--m0_scale=10`, MNI output, atlases | F3 (`m0_divisor` 10, atlases) |

The final implementation task removes each old test only after its replacement passes in CI.
At that point the plan rechecks every flag in the old test's CLI against its replacement's.
Before deleting markers, it uses `grep` to confirm no other test uses them, and updates
`pyproject.toml` `addopts` and `markers`. The nonqtab download stays for unit tests and
`test_00*`.

## 7. Physics notes used by Tier B

Tier B ratios are expected to be below 1 for single-delay recipes, especially in WM, because of
the model mismatch in Section 2. The bounds encode the measured biases and are not claims of
agreement. A future aslscan option for a white-paper-consistent kinetic model (Section 11)
would allow sharp Tier B bounds for the single-delay equation.

## 8. Known gaps (expected strict xfails)

- **Background-suppression attenuation (resolved).** The white paper (Alsop et al. 2015) gives
  an inversion efficiency of "approximately 95%, so each inversion pulse reduces the ASL signal
  by approximately 5%". ASLPrep applies 0.95ⁿ. aslscan's `inversion_efficiency` ε is the
  fraction inverted, and scales the label difference by (1 − 2ε) per pulse. Its default
  ε = 0.95 keeps only 0.9 of the signal per pulse, so an earlier revision reported a spurious
  gap of (0.95/0.9)ⁿ.
  - Every recipe with suppression now sets ε = 0.975 (1 − 2ε = −0.95). The physical-agreement
    items of `fast_bs_le_absent` and F2 are ordinary checks.
- **LabelingEfficiency in the sidecar (open question).** When the sidecar has
  `LabelingEfficiency`, ASLPrep applies no suppression loss. BIDS defines that field as the
  labeling efficiency alone, so with suppression ASLPrep's CBF is 0.95ⁿ of the physically
  calibrated value (F3: 0.81 for four pulses).
  - It is an expected failure, asserted first at its predicted size, until it is decided how
    ASLPrep should read the field.
- **BS pulse count.** With `BackgroundSuppressionNumberPulses` absent, ASLPrep assumes 1. The
  recipes state the pulse count explicitly so that the items above are isolated.

Any other xfail requires investigation first. A suspected ASLPrep bug is reported with
evidence before any test or code change.

## 9. Geometry verification (fast tier)

1. **`rigid_ras` and ITK text.** Write a transform for `R` with a nonzero translation and a
   rotation about all three axes. Resample a landmark image (a 3D grid of small blobs at known
   world points) with `nitransforms` using that file. Assert each blob lands at the
   independently computed point, to within 0.1 mm.
2. **Derivative transforms.** Resample the `tfmni` phantom T1w (the offset version) to MNI
   through the generated `from-T1w_to-MNI` file. Compare with the TemplateFlow T1w: r > 0.999,
   and a landmark error below 0.5 mm.
3. **HMC and coregistration frames.** Build a synthetic `P_v`. Write ITK files in ASLPrep's
   conventions for a known `H_v` and `C`. Assert that `E_v R^-1 = I` when they are consistent
   and that the error is detected when they are not.
4. **Motion convention.** Run aslscan on a crop with a single known rotation about a
   non-central axis. Recover the pose from image landmarks and assert it matches
   `desc-motion_gt.tsv` under the pinned composition.
5. **FSL-parameter conversion.** Cover translations, combined rotations, a non-central origin
   and an oblique grid.

## 10. Risks and mitigations

- **Build cost on a cache miss.** Installing rustup, fetching crates and building take about
  2-4 min, and happen only when the spec changes.
- **Upstream changes.** Pinned SHAs are verified after checkout. A missing SHA fails the job;
  there are no fallbacks.
- **TemplateFlow drift or unavailability.** File checksums are pinned and downloads retried; a
  failure fails the job.
- **Thread nondeterminism in ASLPrep.** Coreg and motion bounds include margins measured over
  three CI runs. Generator determinism is verified at the fixed thread count.
- **Runtime.** It is remeasured on the 39-slice grid before resource classes are fixed. Estimates
  with derivatives and `--sloppy`: F1-F6 10-20 min each; F7 (smriprep) 30-45 min.
- **Scope.** Each recipe can be developed and landed independently. The plan sequences F1 and
  the fast tier first.

## 11. Follow-ups (out of scope)

- Real-subject phantoms (Section 4.2 interface): data hosting and consent, then add
  `subject:<id>` recipes.
- Equation-level OSIPI tests (#679, level 1).
- Susceptibility distortion correction scoring with a phantom B0 field and generated `fmap/`
  inputs.
- An upstream aslscan option for a white-paper-consistent kinetic model.
- An aslscan release or wheel, which would let CI skip the Rust build.

## 12. Review log

The Codex adversarial review of revision 1 raised 18 findings. All were verified against the
source and accepted:

| # | Finding | Resolution |
|---|---|---|
| 1 | Transform direction | pull/point conventions (4.4, 4.6) and landmark tests (9) |
| 2 | aslref frame | `P_v`, `H_v`, `C` algebra and guarded native comparison (4.6) |
| 3 | No mask in `ComputeCBF` | small cropped fast fixtures (5) |
| 4 | Smoothing and M0 T1 | `fwhm=0` exact tier, smoothing variant, T1 tissue 1.607 s (4.7) |
| 5 | Motion convention | relative transforms first, validated FSL conversion (4.7, 9) |
| 6 | PV construction | exact separable overlap, truth-derived references, contamination metric (4.4, 4.7) |
| 7 | Multi-delay Tier A | all observations, stated fit settings, mutation tests (4.7) |
| 8 | xfail scope | module-scoped run, one item per metric, `strict` xfail per item (4.8) |
| 9 | 39 slices | committed geometry and a pre-check (4.3) |
| 10 | BS factors | corrected (2, 8) |
| 11 | Silent coverage loss | coverage guards, metric manifest, quantiles, spatially varying truth, per-unit bounds (4.2, 4.7) |
| 12 | Replacement coverage | equivalence map, more variants, removal last (5, 6) |
| 13 | Cache and manifest | template checksums, per-fixture manifest, guard first, LF (4.5, 4.9) |
| 14 | Build command | features, `--bin`, explicit path (4.1) |
| 15 | Fixture enforcement | `ASLPREP_REQUIRE_FIXTURES`, lock and atomic publish, host-side generation (4.4, 4.5, 4.10) |
| 16 | Schemas | complete tree and contract (2, 4.2, 4.4) |
| 17 | Atlases | corrected (2) |
| 18 | Wording | corrected (2, 4.1, 5) |

Also added in revision 2 at the user's request: pluggable phantoms with a real-subject interface
(4.2).

The Codex adversarial review of plan revision 1 raised 19 findings. These changed the spec (all
accepted after checking the code):

| Plan finding | Spec change |
|---|---|
| Calibration mode and bounds fitted to output | a-priori acceptance ceilings, Tier B centred on an independent expectation, physical-agreement items, no calibration mode (4.7, 8) |
| Too few voxels | `valid` and `dom_L` masks, no erosion, populations measured per recipe (4.7) |
| `cbf_only` semantics | explicit connections (4.7) |
| Native outputs need the `asl` space | `asl` added to output spaces (5) |
| No recorded motion centre | derived and recorded (4.6) |
| r on flat truth | restricted to `tfmni` (4.7) |
| Metric completeness | metric contract and shared invariants (4.7) |
| CI dependencies, JUnit name, `$BASH_ENV`, immutable cache | `CACHE_EPOCH` and job wiring (4.5, 4.9) |
| Local runner | mirrors CI (4.10) |
| Noise schema | `snr_per_pair` (4.3) |

## 13. Acceptance criteria

1. On a clean machine, `--build-aslscan` then `--generate` produce every recipe, and `--verify`
   passes.
2. The committed spec matches `spec_text()`, and CI enforces it before using the cache.
3. The fast tier passes, includes the mutation tests, and adds at most 3 minutes, as measured.
4. Every `aslscan_*` item passes or strict-xfails with a measured reason. `truth_score.json` is
   in the artifacts, and the JUnit results are in the Tests tab.
5. The `examples_*` tests are removed only per Section 6. `test_00*` and `qtab` are untouched.
6. The developer documentation covers generation, scoring, tolerances, the bump procedure, and
   adding a phantom.
