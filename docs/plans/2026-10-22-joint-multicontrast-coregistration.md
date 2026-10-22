# Plan: joint multicontrast coregistration probe

Branch: `coreg-multicontrast-eval`. Status: revision 3, modelled on the HMC probe
(`docs/plans/2026-10-13-hmc-probe.md` on branch `worktree-hmc-probe`), after one Codex
adversarial review (log at the end).
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
  otherwise CBF, delta-M, an included M0, a separate M0, control), N4 bias-corrected
  (`init_synthstrip_aslref_wf`, or the legacy path, #753), and unwarped when a fieldmap or
  gradient-nonlinearity correction applies. With SynthStrip the brain mask is a separate output:
  the reference written as `desc-coreg_aslref`, and passed to Stage 5, is not masked.
- **Stage 5** registers it with fMRIPrep's `init_bold_reg_wf`: `mri_coreg` initialization, then
  FLIRT-BBR without FreeSurfer or `bbregister` with it (whenever `run_reconall` is on), 6 dof by
  default (`--asl2anat-dof`), initialized per `--asl2anat-init` (`t1w`, `t2w`, `header`), with
  fMRIPrep's fallback to the initialization when BBR moves too far (`--force bbr` / `no-bbr`
  override it).

This probe covers the FSL path only: the HMC probe's baseline runs used `--fs-no-reconall`
(`run_panel.py`), and the phantom has no surfaces. Its conclusions do not extend to `bbregister`
runs; a bbregister comparator is open decision 3.

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
| Static registration, then perfusion-only refinement (sequential) | ExploreASL, HCP-ASL | Q1, an ablation every J-arm must beat |
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
| R0 | FSL | coregistration reference -> T1w, BBR default sign, fallback | Production (FSL path), as run |
| R1 | FSL | as R0, BBR sign per family by a rule frozen in Task 3 | Production with the known defect fixed |
| R2 | FSL | perfusion image -> T1w, BBR; one D3 level and eligibility rule frozen in Task 3 | Single perfusion reference (#672) |
| RA | FSL | R1 or R2, chosen per run by the HMC probe's DA rule | Automatic single reference (#672) |
| S1 | ANTs | coregistration reference -> T1w, MI | Static-only ablation |
| S2 | ANTs | perfusion image -> perfusion target, CC | Perfusion-only ablation |
| Q1 | ANTs | S1, then S2 initialized from S1 (sequential) | Sequential ablation, as ExploreASL and HCP-ASL work |
| J1 | ANTs | S1 + S2 in one stage | Primary joint candidate |
| J2 | ANTs | J1 + coregistration reference -> T2w, MI | Joint with T2w (exploratory) |
| J3 | ANTs | coregistration reference -> T1w (MI) and -> T2w (MI) | Dual anatomical target (exploratory) |
| J4 | ANTs | J1 + control/label mean -> T1w, MI | Residual static tissue |
| J5 | ANTs | S1, then J1 from S1's result | Two-stage joint |
| J6 | ANTs | J1 with native CBF as the perfusion image | Model-based perfusion |

R-arms are comparators. S-arms and Q1 isolate the engine and the sequencing: a J-arm is judged
first against S1, S2 and Q1 (same engine, schedules, images, masks, initialization and
resampling), so that neither ANTs-versus-FLIRT differences nor what a sequential strategy already
gives are credited to the joint objective, and then against the frozen comparator.

**The comparator is frozen in Task 3, on development data only**, before any J-arm is screened:
R1's sign rule, R2's D3 level and eligibility rule (the HMC probe eliminated none of its four
levels), and, per family, which of R1, R2 or RA is the bar (the better on development data).
RA is available only if the HMC probe's DA rule is implemented and frozen by then. R2 and RA
apply only where the perfusion image is eligible; each comparison is made on the runs where
both arms apply, and the excluded runs are counted.

J2 and J3 need a T2w, which no simulation has. They cannot pass gates 2 and 3, so they are
exploratory: reported on real data, never adopted under this protocol.

Factors of the joint arms (each fixed in a list before results are seen):

| Factor | Levels |
|---|---|
| Perfusion image | P-winsor (default); P-robust; native CBF (J6 only) |
| Perfusion target | GM probability smoothed to the ASL resolution; pseudo-CBF (a weighted sum of GM and WM probabilities, smoothed, with ExploreASL's weights as read from its code in Task 1) |
| Perfusion weight | 0.2, 0.4, 0.6 of the total |
| Masking | none; the conventional brain mask on the ASL side of the stage (see below) |
| Fixed image | anatomical fixed (ANTs samples in the T1w domain); ASL fixed (the inverse, sampled on the ASL grid) |

**Masks and directions.** `antsRegistration --masks [fixedMask,movingMask]` belongs to a stage and
applies to every metric term in it. Each ANTs arm's specification therefore names, per direction,
the fixed and moving image of every term and the fixed and moving mask of the stage: with the
anatomy fixed, the ASL brain mask is the moving mask and the T1w brain mask the fixed mask; with
ASL fixed, the two swap sides. A term that needs a different support from the others of its stage
is either given its own stage or uses an image built for it (masked to zero outside its support),
and that image's support is reported. Clipping (P-winsor caps values at the 99th percentile
inside the mask) is not masking.

Every ANTs arm starts, in the primary comparison, from production's `mri_coreg` initialization,
so differences come from the objective. Capture range is measured separately (metric 3).

### Bounding the combinations

Crossing every factor for every J-arm is not feasible (J1 alone has 2 x 2 x 3 x 2 x 2 = 48
variants). The screening is staged, and the lists are fixed before results are seen:

1. **Engine pilot (Task 2):** S1 and S2 schedule, MI bins and sampling, CC radius, fixed-image
   direction, on three development fixtures and three development real runs; frozen before
   Task 3.
2. **Primary screen (Task 4):** R0-R2, S1, S2 and J1 over the perfusion-image, target, weight and
   masking levels on the development fixtures. At most two J1 variants go on.
3. **Extension screen (Task 5):** J4-J6 built on the J1 variants kept, each where it applies. At
   most two extensions go on, judged against J1.
4. **Held-out simulation and real data (Tasks 5 and 6):** at most four J-arms, each with its own
   S1, S2 and Q1 ablations carried along, plus R0 and the frozen comparator; J2 and J3 as
   exploratory arms on the T2w subset.

If J1 does not pass gate 2 against S1, S2 and Q1, the extensions are not run, and the probe
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
  labels that also define perfusion (in both `derivatives` and `raw` modes), and `derivatives`
  mode also hands ASLPrep the truth GM map. S2's target is then nearly the generator of its
  moving image, so S2 and every J-arm have an advantage in simulation that real data do not give
  them. No fixture removes this; only real data (or an independently generated anatomy, HMC probe
  open decision 6) do. Simulation checks mechanics, capture range and failure modes; it cannot
  decide.
- **One anatomy.** Held-out simulations share it, so generalization across anatomy rests on real
  data.
- **No extracerebral tissue** in the ASL phantom, and a brain-only synthesized T1w: MI terms see
  an easier problem than on real data.
- No vessels outside `st_vessels`, and smooth perfusion modulation elsewhere.

Real data therefore carry the decision.

### Design

Reused from the HMC probe (recipes in `scripts/hmc_probe/recipes/`, phantoms in `fixtures.py`;
`check_stress.py` already shows each stress fixture's effect). Offsets differ by fixture, and the
`mo_*` recipes set none (the default is zero), so the offset of every fixture is listed in the
Task 0 manifest and fixtures without one are used only where the table says:

| Fixture | T1w offset | What it tests here |
|---|---|---|
| `rc_pcasl`, `rc_pcasl_flip`, `rc_pcasl_grid` | `[2, -3, 1.5, 3, -2, 2.5]` | 2D PCASL; flipped orientation; separate M0 on another grid |
| `rc_grase` | same | Background-suppressed 3D GRASE with inverted static contrast (R0's sign failure) |
| `rc_pasl_multipld` | same | Multi-delay PASL (J6) |
| `st_spiral` | same | Spiral ringing (R0 fails, R2 recovers) |
| `st_vessels` | same | Macrovascular perfusion away from grey matter (S2 and J1 should suffer) |
| `st_noperf` | same | A territory with no perfusion (perfusion disagrees with GM) |
| `headmotion` | `[3, -4, 2, 5, -3, 4]` | Motion with a T1w offset (end-to-end score) |
| `grasemotion` | `[-2, 3, 1.5, -3, 2, 4]` | Motion with a T1w offset, GRASE |
| `mo_pcaslshort2_drift`, `mo_pcaslshort4_drift` | none | Two and four pairs; regenerated here with the `rc_*` offset (new recipes) before use |
| `mo_*_ho`, `mo_*_ho2` | none | Held-out trajectories and seeds; only with an offset added, as above |

New, in this probe's directory (not under `aslprep/tests/data/aslscan/`, so the CI fixture digest
is unchanged):

- **`raw` variants** of `rc_pcasl`, `rc_grase`, `st_spiral`, `st_vessels` and `st_noperf`: the same
  recipe with `anat = 'raw'`, so ASLPrep segments the T1w. The difference between a J-arm's score
  on the two forms measures only how much it gains from being handed the truth tissue maps; the
  shared-label circularity remains in both. Each `raw` segmentation is scored against the truth
  (Dice per tissue, with a threshold set in Task 3) before the fixture is used.
- **A second anatomical offset** (larger rotation, about 8 degrees) on `rc_pcasl` and `rc_grase`,
  for the capture range at realistic misalignments.
- **Held-out still fixtures:** `rc_*` with new seeds and a different offset, generated before the
  development scores are seen. Thresholds are fixed before they are scored.

### Metrics

1. **Accuracy:** RMS displacement over the brain of `C^-1 R^-1` (the registration's share,
   motion-free fixtures) and of the end-to-end map through the HMC transforms (motion fixtures),
   in acquisition voxels, and the rotation as arc length at 70 mm, as `score_frames` computes them.
2. **Failure:** error above one acquisition voxel; and, for R-arms, whether BBR fell back.
3. **Capture range:** a separate protocol from the primary comparison. The initial pull transform
   (T1w point -> aslref point) is the truth `C* = R^-1` composed with a perturbation `P` applied
   on the T1w side, `C0 = C* P`, with `P` a random rigid transform about the brain centre at 2, 5,
   10 and 15 mm and degrees (20 draws per level, fixed seeds, the same draws for every arm). Every
   arm's own initialization is replaced by `C0`: for R-arms, FLIRT-BBR starts from `C0` with
   `mri_coreg` skipped; for ANTs arms, `--initial-moving-transform` is `C0` in the matching
   direction. The score is the fraction that converge within 0.5 voxel. The same draws applied to
   `mri_coreg`'s output, with `mri_coreg` kept, are reported separately as robustness of the
   production chain. A test checks that `C0 = C*` (no perturbation) scores zero for every arm's
   input convention.
4. **Contrast of each source**, measured from the truth (local contrast across the white-matter
   boundary and the grey/white ratio, as the HMC probe defines it): explains failures, not ranked.
5. **Cost:** run time and failures to run.

## Real data

No truth exists, so three designs are combined:

1. **Native:** blind review of alignment, and transform agreement, on data as acquired.
2. **Cross-run consistency:** in sessions with several runs, the anatomical transforms of runs a
   and b must agree with a direct registration of run a's reference (moving) to run b's (fixed),
   same contrast, same session. ASLPrep's stored `from-aslref_to-T1w` matrix `C` is a pull map
   (T1w point -> aslref point; `truth_scoring.py`), so the map from run b's aslref points to run
   a's implied by the two coregistrations is `C_a C_b^-1`, and the direct registration `D_ab` is
   the pull map in the same direction (fixed b point -> moving a point). The score is the RMS displacement of `D_ab^-1 C_a C_b^-1`
   over run b's brain. No arm optimizes it, and it needs no anatomical target. Tested in Task 2
   with the identity and a known relative transform on oblique headers.
3. **Semi-synthetic:** each run's T1w (and its masks and tissue maps) given a new header,
   `A' = M A`, for a known rigid world transform `M` (no resampling). A T1w point `x'` in the new
   world is `x = M^-1 x'` in the old, so an arm whose solution on the original is `C` should
   return `C M^-1` (pull convention). With `C'` its actual solution on the moved T1w, the score is
   the RMS displacement of `C'^-1 C M^-1` over the brain in the new world. This tests robustness and capture range on real contrast, not accuracy.

Decision metrics, each with a direction fixed in advance:

- **Adjudicated alignment** on the reviewed subset: reviewers, blind to the arm, compare pairs of
  arms on the same run (forced choice, ties allowed) using the T1w white-matter boundary over the
  coregistration reference and over the perfusion image, in a common physical frame. Two
  reviewers; their disagreement is the uncertainty. A pair is shown only when the two transforms
  differ by more than a tie distance (RMS over the T1w brain, in mm), set in Task 3 from repeat
  ratings of development pairs; pairs closer than that are counted as ties without review. The
  result is reported with its denominators: the runs in the set, the pairs that differ, the
  informative choices, and the win fraction among them.
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
| Dataset | `tarricone` (multi-delay GRASE), `ds008719` (2D EPI), and one spiral scanner group chosen in Task 0 (a spiral dataset goes with every dataset that shares its scanner) | Subject, site and protocol |
| Scanner | Datasets grouped by scanner serial number where the HMC probe's manifest records one (`ds006131`, `ds000240` and `ds008579` share 66044; `ds006148` and `ds007908` share 67080; `ds007116` and `shanmugan_hup6` share 167024), each group on one side. Datasets with no recorded serial (`ds008020`, `ds004873`, `ds004856`, `tarricone`, `ds000234`, `ds000235`, `ds000236`, `ds005752`, `ds008719`, `2d_mb_pcasl`, `ds000254`) are each their own group, labelled "dataset", not "scanner" | Scanner, only for the groups with a serial |
| Run | Within multi-run subjects, half the sessions | Session; not subject |

The HMC probe never recorded its split, so this one is new; if the HMC probe later fixes one, the
two are reconciled before Task 6 and any overlap is reported.

#### The reviewed subset

Fixed in Task 0, before any arm is run: at least one run per family; every spiral and GRASE run;
named cases for each condition (no M0, separate M0, inverted contrast, single delta-M, T2w
available, fieldmap available). Its held-out part is the **primary review set**, and gate 4 is
decided on it.

An **enrichment set** may be added in Task 6 if the primary set has fewer than 20 informative
disagreements: further held-out runs, drawn by a rule fixed in Task 3 that uses no review result
(the J-arm and the frozen comparator differ by more than the tie distance, in random order, until
20 informative disagreements or the panel is exhausted). It is reported separately, with its
selection rule and denominators, and it can confirm but not replace a primary-set result.
An applicability table lists, for every run and every metric and gate, whether it
applies, is not available on that path, or is not tested; a metric that is not available is never
counted as a pass.

## Decision gates

Applied in order; an adoptable arm (J1, J4, J5, J6) must pass each one. J2 and J3 are
exploratory and are not gated. Each outcome (accuracy, failure, capture range, real-data
alignment) is judged separately; trade-offs are reported per metric, not folded into one score.

1. **Implementation:** each engine recovers a known rigid transform applied to an aslref; ANTs and
   FLIRT outputs compose to ASLPrep's `from-aslref_to-T1w` convention and score as ASLPrep's own
   transform does on one fixture; the unperturbed capture-range start scores zero; repeated runs
   give identical transforms.
2. **Joint benefit (simulation):** the J-arm is better than each of S1, S2 and Q1 on at least one
   stress family (`st_*`, `rc_grase`, short series) and not worse than the best of them elsewhere
   by more than a margin, on the `raw` forms as well as the `derivatives` forms. Screened on the
   development fixtures, then **replicated on the held-out fixtures** with the same S1, S2 and Q1
   run alongside: the held-out result must show the same direction in the same stress families.
3. **Against the comparator (simulation):** failure rate not above the frozen comparator's, capture
   range at least as wide, and median accuracy within a margin of it, on development and then
   held-out fixtures.
4. **Real data (held-out, primary review set):** adjudication prefers the J-arm over the frozen
   comparator of each family by more than the bootstrap uncertainty, resampling subjects (sessions
   of a subject together), with at least 20 informative disagreements (fewer is "insufficient
   evidence"; an enrichment set is reported but cannot decide); the J-arm is not worse than Q1 by
   the same test; cross-run inconsistency and semi-synthetic recovery not worse than the
   comparator's in any judged family by more than a margin; no more gross failures.
5. **Operation:** deterministic (fixed ANTs seed and thread count), under the run-time budget, and
   no new dependency (ANTs is already in the image).

Margins, the tie distance, the comparator per family, the enrichment rule and the run-time budget
are written down in Task 3 from development data only, before any J-arm is scored and before any
held-out fixture or run is scored by any arm. Among arms that pass, the simplest is chosen (Q1, if
the J-arm does not beat it, is a result in its own right). A default change additionally needs
gate 4 in every judged family.

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

- Generate the `raw` variants (and score their segmentations), the offset versions of the short
  series, the second-offset fixtures and the held-out still fixtures. Held-out fixtures are
  generated now and not scored until Task 5.
- On development fixtures and development real runs only: measure R0, R1 (both signs) and the four
  D3 levels; freeze R1's sign rule, R2's level and eligibility rule, and the comparator per family;
  repeat ratings of development pairs for the tie distance; set the margins, the enrichment rule
  and the run-time budget.

### Task 4: Primary screen

S1, S2, Q1 and the J1 variants on the development fixtures; gate 2 (development part); keep at most
two J1 variants.

### Task 5: Extensions and held-out simulation

J4, J5 and J6 on the development fixtures where they apply; gates 2 and 3 against J1; at most two
go on. Then gates 2 and 3 on the held-out fixtures for the kept arms, with S1, S2, Q1, R0 and the
comparator run alongside.

### Task 6: Real data

- Write the real-data criteria (directions and margins) before any real-data J-arm score is seen.
- Run the kept arms, their S1, S2 and Q1, R0 and the comparator on the panel; J2 and J3 on the T2w
  subset, reported as exploratory.
- The blind review of the primary review set; the enrichment set only if the rule triggers;
  cross-run and semi-synthetic metrics on the whole panel.
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
- Clipping, image masking and registration masks are three different things, and each arm's
  specification says which it uses: P-winsor clips values; a term-specific support needs its own
  stage or an image masked to zero outside it, because ANTs masks apply to every term of a stage;
  stage masks follow the fixed/moving direction ("Masks and directions").
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
2. Whether RA must be frozen in the HMC probe first, or the better of R1 and R2 per family is a
   sufficient bar.
3. Whether to add `bbregister` comparators (FreeSurfer surfaces for the phantom, and reruns of the
   panel with FreeSurfer), so that the conclusions cover production's FreeSurfer path.
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
ordered decision gates, and review checkpoints.

**Codex adversarial review of revision 2 (2026-10-22).** Fourteen findings, each checked against
the source; all are addressed in revision 3, two in part (9, 13).

| # | Finding | Change |
|---|---|---|
| 1 | The coregistration reference is bias-corrected, not skull-stripped; SynthStrip's mask is a separate output | "Current method" corrected (`reference.py`) |
| 2 | R0 is FSL only, but production uses `bbregister` whenever FreeSurfer runs | Scope stated as the FSL path (the baselines used `--fs-no-reconall`); bbregister is open decision 3 |
| 3 | `C_b^-1 C_a` does not compose: the stored matrix is a pull map (T1w -> aslref) | Cross-run score is `D_ab^-1 C_a C_b^-1` with `D_ab` defined; semi-synthetic expectation written algebraically; both tested in Task 2 |
| 4 | The ASL mask cannot be the moving mask in both directions, and stage masks cannot differ per term | "Masks and directions" added: masks per direction, separate stages or masked images for term-specific support |
| 5 | P-winsor clips; it does not mask | Clipping, image masking and registration masks distinguished in the implementation rules |
| 6 | No sequential comparator, so a joint objective could take credit for what sequencing gives | Q1 (S1, then S2 from S1) added and carried through every gate |
| 7 | Capture range was not operationally defined | Protocol defined: `C0 = C* P`, the same draws for every arm, initialization replaced, production-chain robustness reported separately, zero-perturbation test |
| 8 | J2 and J3 cannot pass the simulation gates | Made exploratory: reported, never adopted under this protocol |
| 9 | Enlarging the review set by disagreements conditions the sample, and "best comparator" was undefined | Partly: the enrichment rule uses no review result, so it is kept, but the primary held-out review set is frozen in Task 0 and alone decides gate 4; the enrichment set is reported separately; the comparator is frozen per family on development data |
| 10 | The joint benefit was never replicated out of sample, and R0-R2 were measured on every fixture before margins were set | S1, S2 and Q1 carried to the held-out fixtures and runs, gate 2 replicated there; Task 3 measures on development data only |
| 11 | `raw` mode still synthesizes the T1w from the perfusion labels | Claim narrowed to leakage of the truth tissue maps; segmentations scored; circularity stated as not removable in simulation |
| 12 | R2 ("best D3 level") and the R1 sign rule were not frozen | Both frozen in Task 3 on development data, with eligibility and common-support comparisons |
| 13 | The scanner-serial split is not supported by the zoo catalog | Partly: the HMC probe's manifest records serials for some datasets (Codex could not read it), and the groups are now listed; datasets without one form their own "dataset" groups |
| 14 | Not every reused fixture has the stated offset; `mo_*` recipes have none | Offsets listed per fixture; short series regenerated with an offset |

Codex could not read the WSL-side probe data (manifest, baseline runs, `coreg_5b_real.tsv`), the
GitHub issues, or the attached proposal, and did not check the literature.
