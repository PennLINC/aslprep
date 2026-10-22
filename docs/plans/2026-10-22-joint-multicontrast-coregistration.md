# Plan: joint multicontrast coregistration probe

Branch: `coreg-multicontrast-eval`. Status: revision 2, modelled on the HMC probe
(`docs/plans/2026-10-13-hmc-probe.md` on branch `worktree-hmc-probe`); no Codex review yet.
Issues: #672 (scope), #186 (motivating failure), #343 (closed as a duplicate of the same idea).

## Goal

Decide whether ASLPrep should estimate the ASL-to-T1w transform from several image pairs at once
(one rigid transform, one optimization, several similarity terms), by testing it against the
current coregistration and the best single-contrast alternatives on simulated data with known
truth and on real data. "Better" means:

1. it lands closer to the true anatomical pose where truth exists;
2. it fails less often (fewer gross misregistrations, fewer fallbacks to the initialization);
3. it converges from a wider range of starting positions;
4. it does not depend on one contrast being good: it degrades gracefully when the static-tissue
   or the perfusion contrast is poor;
5. it is deterministic, fast enough, and needs no new heavy dependency.

The outcome may be "adopt behind a flag", "adopt for some acquisition families", or "do not
adopt". A negative result is recorded on #672 with its numbers.

## Current method

`init_asl_fit_wf` (`aslprep/workflows/asl/fit.py`):

- **Stage 4** builds one coregistration reference from the HMC reference: the volume type chosen
  by `select_reference_volume_type` (with `prioritize_m0`: an included M0, then a separate M0;
  otherwise CBF, delta-M, an included M0, a separate M0, control), bias-corrected and
  skull-stripped (`init_synthstrip_aslref_wf` or the legacy path, #753), and unwarped when a
  fieldmap or gradient-nonlinearity correction applies.
- **Stage 5** registers it with fMRIPrep's `init_bold_reg_wf`: `mri_coreg` initialization, then
  FLIRT-BBR without FreeSurfer or `bbregister` with it, 6 dof by default (`--asl2anat-dof`),
  initialized per `--asl2anat-init` (`t1w`, `t2w`, `header`), with fMRIPrep's fallback to the
  initialization when BBR moves too far (`--force bbr` / `no-bbr` override it).

BBR assumes T2-like contrast (grey brighter than white): `contrast_type='t2'` for `bbregister`,
FLIRT's default BBR slope. One image, one contrast assumption, one cost.

## What is already known

### Open issues

| Issue | What it reports |
|---|---|
| #186 | 30-40 % of one site's stack-of-spirals runs misregistered: ringing in single volumes defeats the reference; re-registering mean CBF with FLIRT, then BBR, fixed them |
| #343 | Closed as a duplicate: build the coregistration reference from CBF once it exists |
| #672 | An opt-in perfusion-derived reference (`--asl2anat-reference {aslref,perfusion}`), and an `auto` mode choosing per run |

### HMC probe results (2026-10-15 to 2026-10-17)

The HMC probe took #672 into its scope as policies D3 (a perfusion image as the coregistration
reference) and DA (choose per run). Its Task 5b and Task 7 are the closest prior evidence:

- `coreg_experiment.py` reproduces production's FSL path standalone (`mri_coreg`, then FLIRT-BBR
  with FSL's `bbr.sch`, both BBR signs, fMRIPrep's fallback rule as nipype's `_calc_norm_affine`
  against 15 mm) and matches ASLPrep's own registration.
- Simulation, conventional reference: best on median (0.16-0.22 mm). 3D GRASE: 1.06 mm, a BBR
  sign failure (the reference has white brighter than grey; BBR is rejected and fMRIPrep falls
  back to `mri_coreg`); the opposite sign gives 0.45 mm.
- Simulation, perfusion images (P-mean, P-robust, P-winsor, P-masked): about 0.5 mm when built
  from production's corrected series, about 0.25 mm from better-aligned series.
- Stack-of-spirals fixture: every conventional reference fails (8.7-10.2 mm, the BBR sign) and
  the perfusion image recovers it (0.23-0.35 mm).
- Real data: every conventional and perfusion source on the zoo panel was registered with both
  signs, with eligibility statistics (`/home/tsalo/hmc-probe/coreg-5b-real/coreg_5b_real.tsv`,
  750 rows). No blind review was run, so no real-data verdict exists.
- Lessons carried over: split-half reproducibility ranked a systematically broken result best
  (a systematic artifact reproduces in both halves), and FLIRT between near-identical images
  preferred spurious off-grid solutions. Neither is a safe decision metric on its own.

So the large failures seen so far are a contrast-sign problem in BBR, and a single perfusion
image fixes the spiral case. The question here is narrower than the attached proposal assumed:
does combining contrasts add anything once the sign is right and a perfusion reference exists?

### Literature

| Approach | Source | Note |
|---|---|---|
| Boundary-based registration | Greve and Fischl 2009 | Production's cost; sign-dependent |
| M0 or control to T1w, then CBF to a pseudo-CBF image built from GM and WM probabilities | ExploreASL (Mutsaerts et al. 2020) | Sequential, not joint; details from summaries, to be read in Task 1 |
| Calibration-image BBR, then perfusion-weighted refinement | HCP-ASL (Kirk et al. 2025) | Sequential; from the paper summary and code, to be read in Task 1 |
| Perfusion-weighted image registered with BBR | BASIL `oxford_asl` | Single image |
| Several metric terms, each with its own image pair, in one stage | ANTs `antsRegistration` | The joint optimizer used here; mask and weight behaviour confirmed in Task 1 |
| Multi-image, multi-metric registration | elastix | An alternative engine, not run unless ANTs fails Task 1 |

No published head-to-head comparison of joint against sequential or single-contrast ASL
coregistration was found. The attached proposal
(`aslprep-joint-multicontrast-coregistration.md`, agent-written, outside the repository) is the
source of the candidate pairs.

## Scope

In scope:

- which ASL-space images enter the cost (conventional reference, M0, delta-M or perfusion image,
  CBF, control/label mean) and which anatomical targets (T1w, T2w, GM probability, a pseudo-CBF
  image);
- the joint optimizer, its metrics, weights, masking and schedule;
- two-stage (static, then joint) variants;
- how a joint method compares with production, production with the sign fixed, D3 and DA.

Out of scope:

- head-motion correction (the HMC probe), the reference image construction itself except where
  a candidate needs a new input (delta-M before Stage 5);
- nonlinear or distortion-correcting registration;
- bbregister in simulation (the phantom has no FreeSurfer surfaces; HMC probe open decision 9);
- side-by-side registration reportlets (they belong to a production change).

### Proposals in the issues and the attached proposal

| Proposal | Source | Disposition |
|---|---|---|
| Perfusion-derived coregistration reference | #672, #343, #186 | R2 (the HMC probe's D3), a comparator |
| Choose between conventional and perfusion reference per run | #672 | RA (the HMC probe's DA), a comparator if its rule is frozen |
| Restrict the cost to informative tissue | #672 | A factor: masking in the joint arms |
| The BBR contrast configuration | #672 | R1: the sign per family, as the HMC probe measured it |
| M0 to T1w MI plus delta-M to GM probability CC | proposal | J1, the primary candidate |
| Add M0 to T2w | proposal | J2 (T2w subset, real data only) |
| M0 to T1w and to T2w | proposal | J3 |
| Add the control/label mean | proposal | J4 (no or weak background suppression) |
| Static initialization, then joint refinement | proposal | J5 |
| CBF instead of delta-M | proposal | J6 (multi-delay) |
| Synthetic tissue-contrast target | proposal | A factor of J1's perfusion target (pseudo-CBF), not a separate arm |
| T1w x T2w product target | proposal | Not run; the proposal itself argues against it |
| Uncertainty-weighted metrics | proposal | Not run; revisit only if J1 passes and fails on low-SNR runs |
| CBF-to-T1w BBR | proposal | Covered by R2 (perfusion image, BBR) |

## Baseline

The baseline is production at a pinned revision, called as it is. Every arm consumes the same
products of one baseline run per dataset (HMC, distortion correction, the coregistration
reference, the native series, anatomy), so arms differ only in the registration.

The HMC probe's baseline runs (`/home/tsalo/hmc-probe/runs/baseline`, the zoo panel, 14 GB) and
its fixture runs are reused if they pass Task 0's check: they predate #753 (brain extraction
choice) and #755 (label-offset correction), so a sample is rerun on the pinned revision and the
coregistration reference and production's transform compared. If they differ materially, the
panel is rerun. Success is defined against what production gets at the pinned revision.

## Candidates

| ID | Engine | Cost terms (moving -> fixed) | Role |
|---|---|---|---|
| R0 | FSL | coregistration reference -> T1w, BBR default sign, fallback | Production, as run |
| R1 | FSL | as R0, BBR sign per family | Production with the known defect fixed |
| R2 | FSL | perfusion image (best D3 level) -> T1w, BBR | Single perfusion reference (#672) |
| RA | FSL | R1 or R2, chosen per run by the HMC probe's DA rule | Automatic single reference (#672) |
| S1 | ANTs | coregistration reference -> T1w, MI | Static-only ablation |
| S2 | ANTs | perfusion image -> perfusion target, CC | Perfusion-only ablation |
| J1 | ANTs | S1 + S2 in one stage | Primary joint candidate |
| J2 | ANTs | J1 + coregistration reference -> T2w, MI | Joint with T2w |
| J3 | ANTs | coregistration reference -> T1w (MI) and -> T2w (MI) | Dual anatomical target |
| J4 | ANTs | J1 + control/label mean -> T1w, MI | Residual static tissue |
| J5 | ANTs | S1, then J1 from S1's result | Two-stage |
| J6 | ANTs | J1 with native CBF as the perfusion image | Model-based perfusion |

R-arms are comparators. S-arms isolate the engine: a J-arm is judged first against its own
S-arms (same engine, schedule and initialization), so that ANTs-versus-FLIRT differences are not
credited to the joint objective, and then against the better of R1, R2 and RA.

RA is a comparator only if the HMC probe's DA rule is implemented and frozen by Task 3; until
then the bar is the better of R1 and R2, and that is stated in the results.

Factors of the joint arms (each fixed in a list before results are seen):

| Factor | Levels |
|---|---|
| Perfusion image | P-winsor (default); P-robust; native CBF (J6 only) |
| Perfusion target | GM probability smoothed to the ASL resolution; pseudo-CBF (a weighted sum of GM and WM probabilities, smoothed, with ExploreASL's weights as read from its code in Task 1) |
| Perfusion weight | 0.2, 0.4, 0.6 of the total |
| Masking | none; the conventional brain mask on the moving side (applies to every term of a stage) |
| Fixed image | anatomical fixed (ANTs samples in the T1w domain); ASL fixed (the inverse, sampled on the ASL grid) |

Every ANTs arm starts from production's `mri_coreg` initialization, so differences come from the
objective. Capture range is measured separately (metric 3).

### Bounding the combinations

Crossing every factor for every J-arm is not feasible (J1 alone has 2 x 2 x 3 x 2 x 2 = 48
variants). The screening is staged, and the lists are fixed before results are seen:

1. **Engine pilot (Task 2):** S1 and S2 schedule, MI bins and sampling, CC radius, fixed-image
   direction, on three development fixtures and three development real runs; frozen before
   Task 3.
2. **Primary screen (Task 4):** R0-R2, S1, S2 and J1 over the perfusion-image, target, weight and
   masking levels on the development fixtures. At most two J1 variants go on.
3. **Extension screen (Task 5):** J2-J6 built on the J1 variants kept, each where it applies. At
   most two extensions go on, judged against J1.
4. **Real data (Task 6):** at most four arms plus R0, R1 and R2 (and RA).

If J1 does not pass gate 2 against its own S-arms, the extensions are not run, and the probe
reports a negative result.

## Simulation

### What aslscan provides

- Known rigid offset between the ASL acquisition and the T1w (`anat_offset`), and known per-volume
  motion.
- Background suppression with arbitrary pulse times; PCASL and PASL, single and multi-delay;
  2D EPI with multiband, 3D GRASE, stack-of-spirals.
- In `derivatives` mode, smriprep-style tissue probability maps and a brain mask with the T1w;
  in `raw` mode, only the T1w, so ASLPrep segments it.
- Truth: the anatomical offset `R` and the poses; `truth_scoring.score_frames` scores
  `from-aslref_to-T1w` against them.

### What it does not model

- **No T2w.** A T2w synthesized from the phantom labels would be perfectly aligned and trivially
  informative. J2 and J3 are real-data only.
- **Perfusion and anatomy share one set of tissue labels.** The T1w is synthesized from the
  labels that also define perfusion, and `derivatives` mode hands ASLPrep the truth GM map. S2's
  target is then nearly the generator of its moving image, so S2 and every J-arm have an
  advantage in simulation that real data do not give them. Simulation checks mechanics, capture
  range and failure modes; it cannot decide.
- **One anatomy.** Held-out simulations share it, so generalization across anatomy rests on real
  data.
- **No extracerebral tissue** in the ASL phantom, and a brain-only synthesized T1w: MI terms see
  an easier problem than on real data.
- No vessels outside `st_vessels`, and smooth perfusion modulation elsewhere.

Real data therefore carry the decision.

### Design

Reused from the HMC probe (recipes in `scripts/hmc_probe/recipes/`, phantoms in `fixtures.py`;
all with `anat_offset = [2.0, -3.0, 1.5, 3.0, -2.0, 2.5]`; `check_stress.py` already shows each
stress fixture's effect):

| Fixture | What it tests here |
|---|---|
| `rc_pcasl`, `rc_pcasl_flip`, `rc_pcasl_grid` | 2D PCASL; flipped orientation; separate M0 on another grid |
| `rc_grase` | Background-suppressed 3D GRASE with inverted static contrast (R0's sign failure) |
| `rc_pasl_multipld` | Multi-delay PASL (J6) |
| `st_spiral` | Spiral ringing (R0 fails, R2 recovers) |
| `st_vessels` | Macrovascular perfusion away from grey matter (S2 and J1 should suffer) |
| `st_noperf` | A territory with no perfusion (perfusion disagrees with GM) |
| `mo_pcaslshort2_drift`, `mo_pcaslshort4_drift` | Two and four pairs: a noisy perfusion image |
| `headmotion`, `grasemotion` (`aslprep/tests/data/aslscan/`) | Motion with a T1w offset |
| `mo_*_ho`, `mo_*_ho2` | Held-out trajectories and seeds |

New, in this probe's directory (not under `aslprep/tests/data/aslscan/`, so the CI fixture digest
is unchanged):

- **`raw` variants** of `rc_pcasl`, `rc_grase`, `st_spiral`, `st_vessels` and `st_noperf`: the same
  recipe with `anat = 'raw'`, so ASLPrep segments the T1w. The difference between a J-arm's score
  on the `derivatives` and `raw` forms bounds how much its simulation advantage is circularity.
- **A second anatomical offset** (larger rotation, about 8 degrees) on `rc_pcasl` and `rc_grase`,
  for the capture range at realistic misalignments.
- **Held-out still fixtures:** `rc_*` with new seeds and a different offset, generated before the
  development scores are seen. Thresholds are fixed before they are scored.

### Metrics

1. **Accuracy:** RMS displacement over the brain of `C^-1 R^-1` (the registration's share,
   motion-free fixtures) and of the end-to-end map through the HMC transforms (motion fixtures),
   in acquisition voxels, and the rotation as arc length at 70 mm, as `score_frames` computes them.
2. **Failure:** error above one acquisition voxel; and, for R-arms, whether BBR fell back.
3. **Capture range:** from the truth perturbed by random rigid offsets at 2, 5, 10 and 15 mm and
   degrees (20 draws per level, fixed seeds), the fraction that converge within 0.5 voxel.
4. **Contrast of each source**, measured from the truth (local contrast across the white-matter
   boundary and the grey/white ratio, as the HMC probe defines it): explains failures, not ranked.
5. **Cost:** run time and failures to run.

## Real data

No truth exists, so three designs are combined:

1. **Native:** blind review of alignment, and transform agreement, on data as acquired.
2. **Cross-run consistency:** in sessions with several runs, the anatomical transforms of runs a
   and b must agree with a direct registration of run a's reference to run b's (same contrast,
   same session): the inconsistency between `C_b^-1 C_a` and that direct transform. No arm
   optimizes it, and it needs no anatomical target.
3. **Semi-synthetic:** each run's T1w header moved by a known rigid offset (no resampling), so
   every arm should move its solution by exactly that offset. This tests robustness and capture
   range on real contrast, not accuracy.

Decision metrics, each with a direction fixed in advance:

- **Adjudicated alignment** on the reviewed subset: reviewers, blind to the arm, compare pairs of
  arms on the same run (forced choice, ties allowed) using the T1w white-matter boundary over the
  coregistration reference and over the perfusion image, in a common physical frame. Two
  reviewers; their disagreement is the uncertainty. A pair is shown only when the two transforms
  differ by more than a tie distance (RMS over the T1w brain, in mm), set in Task 3 from repeat
  ratings of the same pairs.
- **Cross-run inconsistency** (lower is better), on sessions with at least two comparable runs.
- **Semi-synthetic recovery:** error of the moved solution against the applied offset.
- **Gross failure:** an adjudicated failure, or a transform more than 5 mm (RMS over the brain)
  from every other arm's.

Diagnostics, with no assumed good direction:

- agreement with R0, and with each other arm;
- local contrast across the white-matter boundary of each source, through each arm's own
  transform (it is what BBR optimizes);
- the grey/white CBF ratio in T1w space (what the perfusion term optimizes);
- each arm's cost evaluated at the other arms' solutions;
- BBR fallback rates, by family.

Smoother maps, higher grey/white ratios and agreement with production alone do not count as
better accuracy. Results are stratified by acquisition family.

### The real-data panel

The HMC probe's panel (`/home/tsalo/hmc-probe/panel/manifest.tsv`, 78 runs from the
`aslprep-zoo` collection, read in place and never modified; adapted copies as `panel_config.py`
defines them). What matters for coregistration:

| Family | Datasets | Notes for this probe |
|---|---|---|
| 3D stack-of-spirals | `ds000234`, `ds000235`, `ds000236`, `ds000240`, `ds007116`, `ds008579`, `ds006131`, `shanmugan_hup6` | The #186 family; R0 expected to fail on ringing; `ds007116` and `ds008579` have T2w |
| 3D GRASE | `ds007908`, `ds006148`, `ds003823`, `kashyap`, `ds007522`, `jsr-aslprep-debug`, `lsa-aslprep-debug`, `tarricone`, `grotzinger_menopausal_asl` | Background-suppressed, inverted static contrast likely; `ds007522` and `tarricone` have T2w |
| 2D EPI and multiband | `ds008719`, `morys_bad_brainmask`, `ds004856`, `ds004146`, `ds004873`, `2d_mb_pcasl`, `ds000254` | Low resolution; `ds004146` has T2w |
| 3D EPI at 7 T | `ds008020` | 17 runs in two subjects: the main cross-run set |
| GE product output | `ds005752` | A single delta-M: the perfusion term has one volume and no eligibility statistic (reported separately, never pooled); has T2w |

- **T2w subset (J2, J3):** `ds005752`, `ds007522`, `ds004146`, `ds007116`, `ds008579`,
  `tarricone`. The baseline runs did not preprocess the T2w, so the probe registers each T2w to
  its preprocessed T1w (ANTs rigid, same session); every such registration is checked by eye once,
  before any arm is scored.
- **Cross-run sets:** `ds008020`, `kashyap`, `ds007522`, `ds003823`, `ds000234`, `ds000235`
  (within a session only; the list is confirmed in Task 0).
- **No M0:** `ds007908`, `ds004856` and others: the static term uses the conventional reference
  whatever its source, and results are stratified by M0 type.
- Families with one dataset (`ds008020`, `2d_mb_pcasl`, `ds000254`, `ds005752`) are reported but
  not judged.

#### Development and held-out data

Fixed in Task 0, before any arm is scored, following the HMC probe's axes; each held-out set is
labelled with what it is independent of, and a claim is made only along that axis:

| Axis | Held out | Independent of development data in |
|---|---|---|
| Dataset | `tarricone` (multi-delay GRASE), `ds008719` (2D EPI), and one spiral dataset chosen in Task 0 | Subject, site and protocol |
| Scanner | Datasets grouped by scanner serial number, each group on one side | Scanner |
| Run | Within multi-run subjects, half the sessions | Session; not subject |

The HMC probe never recorded its split, so this one is new; if the HMC probe later fixes one, the
two are reconciled before Task 6 and any overlap is reported.

#### The reviewed subset

Fixed in Task 0: at least one run per family; every spiral and GRASE run; named cases for each
condition (no M0, separate M0, inverted contrast, single delta-M, T2w available, fieldmap
available); then enlarged, during Task 6, with runs where the J-arm and the best comparator differ
by more than the tie distance, until at least 20 informative disagreements exist or the panel is
exhausted. An applicability table lists, for every run and every metric and gate, whether it
applies, is not available on that path, or is not tested; a metric that is not available is never
counted as a pass.

## Decision gates

Applied in order; an arm must pass each one. Each outcome (accuracy, failure, capture range,
real-data alignment) is judged separately; trade-offs are reported per metric, not folded into one
score.

1. **Implementation:** each engine recovers a known rigid transform applied to an aslref; ANTs and
   FLIRT outputs compose to ASLPrep's `from-aslref_to-T1w` convention and score as ASLPrep's own
   transform does on one fixture; repeated runs give identical transforms.
2. **Joint benefit (simulation, development):** a J-arm is better than both of its S-arms on at
   least one stress family (`st_*`, `rc_grase`, short series) and not worse than the better of them
   elsewhere by more than a margin, on the `raw` forms as well as the `derivatives` forms.
3. **Against the comparators (simulation):** failure rate not above the better of R1 and R2, capture
   range at least as wide, and median accuracy within a margin of it, on development and then
   held-out fixtures.
4. **Real data (held-out):** adjudication prefers the J-arm over the better of R1, R2 (and RA) by
   more than the bootstrap uncertainty, resampling subjects (sessions of a subject together), with
   at least 20 informative disagreements (fewer is "insufficient evidence"); cross-run
   inconsistency and semi-synthetic recovery not worse than that comparator's in any judged family
   by more than a margin; no more gross failures.
5. **Operation:** deterministic (fixed ANTs seed and thread count), under the run-time budget, and
   no new dependency (ANTs is already in the image).

Margins, the tie distance and the run-time budget are written down in Task 3, before any J-arm is
scored. Among arms that pass, the simplest is chosen. A default change additionally needs gate 4 in
every judged family.

## Tasks

### Task 0: Panel and fixture preflight

- Pin the baseline revision and the container image by digest.
- Rerun a sample of the panel (one run per family) and two fixtures on the pinned revision; compare
  the coregistration reference and R0's transform with the HMC probe's baseline. Reuse the baseline
  if they agree within the tie distance's starting guess (0.2 mm), otherwise rerun the panel.
- A manifest of every run: family, M0 type, background suppression, delays, reference type, T2w,
  fieldmap, runs per session, which arms and metrics apply.
- Fix the development and held-out lists, the cross-run sets and the reviewed subset.

### Task 1: Harness

Probe code in `scripts/coreg_mc/`, importing the HMC probe's modules where they already do the job
(`coreg_experiment.py` for R0-R2, `perfusion.py` for the perfusion images and eligibility),
ported rather than imported if the two branches must diverge:

- inputs: from one baseline run, the coregistration reference, the native preprocessed series, the
  native M0 and CBF, delta-M and the perfusion levels, the anatomical targets (smoothed GM
  probability, pseudo-CBF, T2w registered to the T1w);
- `run_ants.py`: one `antsRegistration` call per arm from a table of cost terms, writing an ITK
  transform in ASLPrep's convention;
- RA from the HMC probe's DA rule, if frozen;
- read the ExploreASL and HCP-ASL registration code and record what each does.

### Task 2: Scorer and engine pilot

- Metrics 1-5 in simulation, reusing `truth_scoring` and `truth_geometry`; cross-run inconsistency,
  semi-synthetic recovery and the diagnostics on real data.
- Validation: the identity, the truth, deliberately wrong transforms and a transform expressed in
  the inverse direction must score as expected (gate 1).
- The engine pilot of "Bounding the combinations", frozen in this plan.

Checkpoint: Codex adversarial review of the harness and scorer.

### Task 3: Fixtures and thresholds

- Generate the `raw` variants, the second-offset fixtures and the held-out still fixtures; check
  the `raw` segmentation against the truth once.
- Measure R0-R2 on every fixture; set the margins, the tie distance's starting value, and the
  run-time budget.

### Task 4: Primary screen

S1, S2 and the J1 variants on the development fixtures; gate 2; keep at most two J1 variants.

### Task 5: Extensions

J4, J5 and J6 in simulation where they apply (J2 and J3 wait for real data); gates 2 and 3 against
J1; at most two go on. Then gate 3 for the kept arms on the held-out fixtures.

### Task 6: Real data

- Write the real-data criteria (directions and margins) before any real-data J-arm score is seen.
- Run the kept arms and the comparators on the panel; J2 and J3 on the T2w subset.
- Repeat ratings for the tie distance; the blind review of the reviewed subset, enlarged with
  disagreements as described; cross-run and semi-synthetic metrics on the whole panel.
- Gate 4 on the held-out runs.

Checkpoint: Codex adversarial review of the results before the decision is written.

### Task 7: Decision and report

- Apply the gates; write the decision into this plan and as a comment on #672.
- A walkthrough HTML report (`scripts/coreg_mc/report.py`, every number computed from the result
  files), published as an artifact and republished as results change.
- If an arm is adopted: a separate PR adding it behind a flag (`--asl2anat-reference` or a new
  `--asl2anat-method`), with delta-M built before Stage 5, followed by an ordinary Codex code review.

### Review checkpoints

- Codex adversarial review of this plan before Task 0.
- Codex adversarial review of the harness and scorer after Task 2.
- Codex adversarial review of the results in Task 6, and of any implementation in Task 7, then an
  ordinary Codex code review as the final pass.

## Implementation rules

- Every arm reads the same baseline products; none reruns HMC, distortion correction or anatomy.
- Moving images are on the coregistration reference's grid and frame; a separate M0 is composed
  through its own transform, never assumed aligned.
- Build perfusion images from the native preprocessed series; resample measured ASL data once.
- Contrast-specific masking is built into the moving image (as P-winsor clips), because ANTs masks
  apply to every term of a stage.
- Record, per run and arm, the terms, weights, masks, initialization, the final cost per term, and
  any fallback.

## Conventions

- Local commands run in WSL in the `hmc-probe` micromamba environment (the HMC probe's code is
  editable there); pytest is `python -m pytest`.
- ANTs, FSL and FreeSurfer run in the pinned ASLPrep image.
- Git runs on Windows (the worktree's gitdir is a Windows path).
- Probe data live under `/home/tsalo/coreg-mc`; work directories are deleted once scored. The zoo
  is read in place and never modified.

## Open decisions

1. Who the two reviewers are, and the review load per reviewer.
2. Whether RA must be frozen in the HMC probe first, or the better of R1 and R2 is a sufficient bar.
3. Whether to run FreeSurfer on the phantom so that bbregister arms can be tested in simulation.
4. Whether to extend the panel (more spiral and GRASE subjects) if the disagreement count stays
   below 20.

## References

- Greve and Fischl 2009, boundary-based registration: https://doi.org/10.1016/j.neuroimage.2009.06.060
- Mutsaerts et al. 2020, ExploreASL: https://doi.org/10.1016/j.neuroimage.2020.117031
- Kirk et al. 2025, HCP-ASL: https://doi.org/10.1162/imag_a_00444; code:
  https://github.com/physimals/hcp-asl
- BASIL `oxford_asl`: https://asl-docs.readthedocs.io/en/latest/oxford_asl_userguide.html
- ANTs: https://github.com/ANTsX/ANTs
- elastix: https://elastix.dev/
- HMC probe plan: `docs/plans/2026-10-13-hmc-probe.md` (branch `worktree-hmc-probe`)

## Review log

**Revision 1 (2026-10-22).** First draft, from the attached proposal.

**Revision 2 (2026-10-22).** Rewritten on the HMC probe's model at the maintainer's request: goal
criteria, a proposal disposition table, named comparators (R-arms) and same-engine ablations
(S-arms), staged screening with fixed lists, the HMC probe's stress and reference-consumer fixtures
reused with `raw` variants as a circularity control, its real-data panel and blind-review design,
ordered decision gates, and review checkpoints. No Codex review yet.
