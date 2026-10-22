# Evaluating joint multicontrast ASL-to-anatomical coregistration

Status: plan, not started. Branch `coreg-multicontrast-eval`, from `main` at 7f5d87f.

## Question

Can one rigid transform, estimated from several ASL/anatomical image pairs at once (for example
M0 to T1w by mutual information plus delta-M to grey-matter probability by cross-correlation),
align ASL data to the T1w more accurately, or fail less often, than ASLPrep's current
coregistration and the best single-contrast alternatives?

This is a development-only evaluation. Nothing in `aslprep/` changes unless the decision at the
end is positive, and integration is then its own PR.

The starting point is an agent-written proposal (saved outside the repository,
`aslprep-joint-multicontrast-coregistration.md`). This plan keeps its hypothesis, metric pairs
and ablation list, and changes the comparators, the simulation design and the real-data
outcomes for the reasons in "What changes from the proposal".

## What exists

### Production coregistration

- `init_asl_fit_wf` (`aslprep/workflows/asl/fit.py`) builds one coregistration reference
  (Stage 4): a single volume type chosen by `select_reference_volume_type` (M0 or separate M0 if
  prioritized, else CBF, delta-M, M0, control), bias-corrected and skull-stripped, unwarped if a
  fieldmap or gradient-nonlinearity correction applies.
- Stage 5 registers it with fMRIPrep's `init_bold_reg_wf`: `mri_coreg` initialization, then
  FLIRT-BBR (no FreeSurfer) or `bbregister`, 6 dof by default, with fMRIPrep's fallback to the
  initialization when BBR moves too far. BBR assumes T2-like contrast (grey brighter than white).
- Related issues: PennLINC/aslprep#186 (30-40 % of one site's runs misregistered; a mean-CBF
  reference fixed them), PennLINC/aslprep#343 (closed as a duplicate), and PennLINC/aslprep#672
  (an opt-in perfusion-derived reference, `--asl2anat-reference`, with an `auto` mode). This
  evaluation is the multicontrast extension of #672 and should be referenced there.

### The HMC probe already tested single perfusion references (D3/DA)

The HMC probe (`C:\Users\tsalo\Documents\linc\aslprep-hmc-probe`, branch `worktree-hmc-probe`,
plan `docs/plans/2026-10-13-hmc-probe.md`, data under `/home/tsalo/hmc-probe`) ran Task 5b and
Task 7 on exactly the neighbouring question:

- `scripts/hmc_probe/coreg_experiment.py` reproduces production's FSL path standalone
  (`mri_coreg`, then FLIRT-BBR with both BBR signs, fMRIPrep's fallback rule) and scores
  transforms against aslscan truth; `perfusion.py` builds the perfusion levels (P-mean,
  P-robust, P-winsor, P-masked) and the eligibility statistic.
- Simulation results: the conventional reference registers best on median (0.16-0.22 mm);
  perfusion images built from production's series reach about 0.5 mm. On the stack-of-spirals
  fixture every conventional reference fails (8.7-10.2 mm) from the BBR sign, and the perfusion
  image recovers it (0.23-0.35 mm). 3D GRASE's 1.06 mm error is also a sign failure (the opposite
  sign gives 0.45 mm).
- Real data: `/home/tsalo/hmc-probe/runs/baseline` holds full ASLPrep runs on the aslprep-zoo
  panel (26 datasets, about 78 ASL files, 14 GB); `coreg-5b-real/coreg_5b_real.tsv` has every
  conventional and perfusion source registered with both signs, with eligibility statistics.
- Lessons that carry over: split-half reproducibility rewarded a systematic artifact (it must
  not be read alone), and FLIRT between near-identical images prefers off-grid solutions.

### Truth-scored simulation

`aslprep/tests/truth_scoring.py` (`score_frames`) scores ASLPrep's `from-aslref_to-T1w`
transform against the aslscan truth (RMS displacement over brain voxels, in voxels of the
acquisition grid, and rotation as arc length at 70 mm). Only `headmotion` and `grasemotion`
set a nonzero `anat_offset` today. aslscan writes a T1w (and, in `derivatives` mode, tissue
probability maps), but no T2w.

## What changes from the proposal

1. **The comparator is production with its known defect fixed, and the best single perfusion
   reference.** Production's large failures in simulation are BBR sign failures. A joint method
   that beats sign-broken production proves little; it must beat (a) production with the BBR
   sign set per family as Task 5b found it, and (b) D3/DA from #672. The adoption bar is the
   better of those, as the probe set for DA.
2. **Engine and objective are separated.** ANTs versus FLIRT is a confound: every
   single-contrast arm is also run in ANTs with the same schedule as the joint arms, so the joint
   arms are compared first with their own single-metric ablations, then with production.
3. **Simulation cannot decide.** In aslscan, perfusion and the T1w are both generated from the
   same tissue labels, and `derivatives` mode hands ASLPrep the truth GM map. Delta-M to GM
   probability is then close to circular. Simulation uses `anat = 'raw'` (ASLPrep segments the
   T1w itself) and stress fixtures where perfusion departs from grey matter. It checks
   mechanics, capture range and failure modes; the decision rests on real data.
4. **Real-data outcomes do not reuse a registration objective.** Grey/white CBF contrast is what
   the delta-M term optimizes, so it is a diagnostic only. Primary outcomes are blinded
   adjudication and cross-run cycle consistency (below), which no arm optimizes.
5. **ANTs masks are per stage, not per metric.** `antsRegistration --masks` applies to every
   metric of a stage, so contrast-specific masking (vessel suppression for delta-M) is done by
   building the moving image (winsorizing as P-winsor does), not with masks.
6. **Many runs have no M0.** The static-tissue term uses the conventional reference
   (`desc-coreg_aslref`) whatever its source, and results are stratified by M0 type.
7. **T2w arms are real-data only.** A T2w synthesized from the phantom's labels would be
   perfectly aligned and trivially informative. Baseline runs have no T2w derivatives, so the
   probe registers each T2w to its preprocessed T1w itself (ANTs rigid, same session), and that
   step is checked visually before use.

## Arms

All arms consume the same baseline derivatives, so they share HMC, distortion correction and
anatomy. ASL-space inputs are on the coregistration reference's grid: `desc-coreg_aslref`, the
native preprocessed series `desc-preproc_asl` (already motion- and distortion-corrected; delta-M
and the probe's P-winsor are built from it), the native M0, and native CBF.

| ID | Engine | Source(s) -> target(s) | Role |
|---|---|---|---|
| P0 | FSL (production) | coreg reference -> T1w, BBR default sign | Production as run |
| P1 | FSL | coreg reference -> T1w, BBR sign per family (Task 5b) | Production, sign fixed |
| D3 | FSL | P-winsor -> T1w, BBR | Single perfusion reference (#672) |
| DA | FSL | P1 or D3, chosen by the probe's rule | Automatic single reference (#672) |
| A | ANTs | coreg reference -> T1w, MI | Static-only ablation |
| B | ANTs | delta-M (P-winsor) -> GM probability, CC | Perfusion-only ablation |
| D | ANTs | A + B in one stage | Primary joint candidate |
| C | ANTs | coreg reference -> T1w (MI) + -> T2w (MI) | Dual anatomical target (T2w subset) |
| E | ANTs | D + coreg reference -> T2w (MI) | Joint with T2w (T2w subset) |
| F | ANTs | D + control/label mean -> best anatomy (MI) | Residual static tissue (no or weak BS) |
| G | ANTs | A, then D initialized from A | Two-stage |
| H | ANTs | D with native CBF instead of delta-M | Model-based perfusion (multi-delay) |

DA is included only if the probe's DA rule was implemented and frozen; otherwise P1 and D3 stand
for it, and the bar is the better of the two.

Every ANTs arm starts from the same initialization (`mri_coreg`, as production does), so
differences come from the objective. Capture range is tested separately.

## Data

**Simulation (development only).**

- Still-head fixtures with a nonzero `anat_offset` and `anat = 'raw'` for each readout aslscan
  supports (2D EPI, 3D GRASE, stack-of-spirals, GE 3D), single- and multi-delay, with and without
  background suppression and separate M0. New recipes live in the probe's directory, not in
  `aslprep/tests/data/aslscan/`, so the fixture spec file is untouched.
- The existing motion fixtures (`headmotion`, `grasemotion`).
- Stress fixtures where perfusion departs from grey matter: focal macrovascular signal, an
  absent-perfusion territory, few pairs (2-4), long arrival with a short delay, and spiral
  ringing. The HMC probe planned these (its Task 4, "Stress cases"); reuse them if they exist,
  otherwise write them here and check each shows its effect before use.

**Real data.** The aslprep-zoo panel through the existing baseline runs (read in place; the zoo
is never modified). The split into development and held-out runs is by subject, stratified by
acquisition family, and fixed before any arm is run: the HMC probe's split if its panel manifest
records one, otherwise a new one written to the plan. T2w subset: `ds005752`, `ds007522`,
`ds004146`, `ds007116`, `ds008579`, `tarricone`. Multi-run subjects for cycle consistency:
`ds008020` (17 runs), `kashyap` (9), `ds007522` (8), `ds003823`, `ds000234`, `ds000235`.

## Outcomes

**Simulation.**

- Error against truth: `C^-1 R^-1` RMS displacement over the brain in acquisition voxels and the
  rotation arc, as `score_frames` computes it (motion-free fixtures: the registration's own share).
- Failure rate: error above 1 voxel.
- Capture range: from the truth perturbed by random rigid offsets at 2, 5, 10 and 15 mm/degrees
  (20 draws per level, fixed seeds), the fraction that converge within 0.5 voxel.

**Real data (primary).**

1. Blinded adjudication: for each development and held-out run where two arms differ by more
   than a tie distance (RMS over the T1w brain, in mm), a pairwise forced choice with labels
   randomized and hidden, showing the T1w white-matter boundary over the coreg reference and over
   delta-M. Ties allowed. The tie distance comes from repeat ratings of the same pairs.
2. Cycle consistency across runs of a session: for runs a and b, the inconsistency between
   `C_b^-1 C_a` and a direct same-contrast registration of aslref a to aslref b. It needs no
   anatomical target and is optimized by no arm.
3. Robustness: from each arm's own solution perturbed as in the capture-range test, the spread of
   the returned transforms.
4. Fallback and outlier counts: BBR fallbacks, transforms far from the `mri_coreg` start.

**Diagnostics only (optimized by some arm).** Local contrast across the white-matter boundary
for each source, the grey/white ratio of CBF in T1w space, each arm's cost at the other arms'
solutions, and agreement with production.

**Engineering.** Run time per arm, behaviour with missing M0 or T2w, and the inputs a production
version would need (delta-M before Stage 5 means moving part of the CBF workflow, as #343 noted).

## Decision rule (fixed before held-out data are scored)

Joint registration (the best of D, E, G, H, chosen on development data only) is worth integrating
as an opt-in method if, on held-out data:

1. adjudication prefers it over the better of P1 and D3 (or DA) by more than the bootstrap
   uncertainty, resampling subjects, with at least 20 informative disagreements (fewer means
   "insufficient evidence");
2. cycle-consistency inconsistency is not worse than that comparator's in any family with at
   least three runs, by more than its own bootstrap uncertainty;
3. its simulation failure rate is not above the comparator's, and its capture range is at least
   as wide;
4. its run time stays under 5 minutes per run in the ASLPrep container.

A default change needs the same result in every family with enough runs. A negative or
insufficient result is recorded on #672 with the numbers.

## Tasks

### Task 0: Setup

- Environment: `aslprep-pvc` for scoring (or `hmc-probe`, which already has the probe's code
  editable); ANTs, FSL and FreeSurfer through the ASLPrep container the baseline runs used,
  pinned by digest. Git stays on Windows (worktree gitdir).
- Working data under `/home/tsalo/coreg-mc`; delete per-arm work directories once scored.
- Inventory the baseline runs: per run, family, M0 type, background suppression, delays, T2w,
  run count per session, and the ASLPrep version that produced it.

### Task 1: Harness

Scripts in `scripts/coreg_mc/` (probe code, not package code):

- `inputs.py`: collect the ASL-space images from a baseline run, build delta-M and P-winsor
  (port `perfusion.py`'s builder and eligibility, citing it), smooth the GM probability map to
  the ASL resolution, and register the T2w to the T1w.
- `run_fsl.py`: P0, P1 and D3, ported from `coreg_experiment.py`.
- `run_ants.py`: one `antsRegistration` call per arm from a table of metric terms; outputs an ITK
  transform in the same convention as ASLPrep's `from-aslref_to-T1w`.
- `score.py`: truth error (reusing `truth_scoring`/`truth_geometry`), cycle consistency,
  perturbation tests, diagnostics.
- Self-tests: recover a known rigid transform applied to an aslref (per engine); check the
  transform direction and composition against `score_frames` on one fixture.

Checkpoint: Codex adversarial review of the harness and scorer.

### Task 2: Pilot and freeze ANTs settings

On three development fixtures and three development real runs: the multiresolution schedule,
MI bins and sampling, CC radius, whether the anatomical or the ASL image is the fixed image,
and the smoothing of the GM map. Frozen in this plan before Task 3. Weights are not chosen here.

### Task 3: Simulation screen (P0, P1, D3, A, B, D)

Truth error, failure rate and capture range on every simulation fixture. Gate to continue: D is
better than both A and B on at least one stress family without being worse than the better of
them elsewhere. If D only matches its better ablation everywhere, stop and report.

### Task 4: Real development screen and weights

P0, P1, D3, A, B and D on development runs. Sensitivity of D to the perfusion weight (0.2, 0.4,
0.6 of the total) and to the delta-M builder; the weight is chosen on development data by
adjudication and cycle consistency, then frozen.

### Task 5: Extensions (only if Task 3 and Task 4 support D)

C and E on the T2w subset; F on runs without or with weak background suppression; H on
multi-delay runs; G everywhere. Each judged against D on development data; at most two go on.

### Task 6: Held-out evaluation

The frozen arms on held-out runs, the blinded adjudication, and the decision rule.
Checkpoint: Codex adversarial review of the results before the decision is written.

### Task 7: Report

A walkthrough HTML report built by `scripts/coreg_mc/report.py` from the result files (every
number in the prose computed from data), published as an artifact and republished as results
change; a summary comment on #672.

## Risks

- **Old baselines.** The baseline runs predate recent changes (brain extraction choice #753,
  label-offset correction #755). All arms share them, so comparisons stay fair, but P0 is
  "production at that version". Re-running a few runs on current `main` checks that the
  coregistration reference has not changed materially.
- **Adjudication load.** Pairwise review across many arms is large; Task 3/4 gating keeps the
  held-out arm count to at most four.
- **Correlated terms.** M0 and the conventional reference are the same image when M0 is the
  reference; F and H are kept out of the primary candidate for that reason.
- **Distortion.** Rigid registration of distorted images; runs with fieldmaps are compared with
  and without unwarping only if Task 4 suggests distortion drives disagreements.

## Open questions

- Who adjudicates (one reviewer or two), and how many pairs is a reasonable load?
- Should the HMC probe's DA rule be finished first, so DA is the comparator, or is the better of
  P1 and D3 enough?
- Are the probe's stress fixtures written, or does this evaluation write them?
