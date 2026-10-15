# PSF-aware partial volume correction and weighted surface projection: design

Status: revision 4 (2026-10-15). This revision addresses the Codex adversarial reviews of
revisions 1, 2 and 3 (Section 11).
Branch: `worktree-pvc-spec`.
Issues: PennLINC/aslprep#681 (native PSF-aware PVC), PennLINC/aslprep#676 (PV- and
uncertainty-aware surface projection). Related: #671 and #726 (BASIL posterior SDs),
#682 (uncertainty as first-class output), #679 (conformance suite).

## 1. Goal

ASL voxels are large (typically 3-4 mm in-plane and 4-8 mm through-plane) and their effective
resolution is coarser still, so most cortical voxels mix gray matter (GM), white matter (WM)
and CSF. ASLPrep's only partial volume correction (PVC) today is BASIL spatial-VB PVC, which
runs whenever `--basil` is on. It receives GM and WM maps that were resampled to the ASL grid
with Gaussian interpolation and no model of the acquisition's point-spread function (PSF).
Surface projection ignores tissue composition and uncertainty altogether.

This design adds:

- **Part A (#681).**
  - A tissue partial-volume (PV) map workflow that models the acquisition's effective
    resolution.
  - An opt-in native PVC method for single-PLD data: local linear regression (LR) after
    Asllani et al. (2008), applied to the calibration-free signal q·ΔM rather than to
    quantified CBF.
  - Status, conditioning and noise diagnostics.
- **Part B (#676).**
  - Opt-in weighted volume-to-surface projection on the Connectome Workbench paths, using
    PV-based and precision-based weights.
  - Undilated per-vertex support outputs.

Both parts are opt-in. With no new flag, the workflow graph and the derivative list are
unchanged.

Non-goals (this change):

- Changing BASIL PVC, its inputs, or its coupling to `--basil`. A follow-up may feed BASIL the
  new PV maps (Section 9).
- LR PVC of multi-PLD data or of scanner-provided CBF (Section 5.1 explains why).
- A user-supplied PSF kernel image.
- Total-least-squares or trimmed regression, Bayesian PVC, and surface-based (Toblerone-style)
  PV estimation.
- Modelling the blur from head motion correction and SDC resampling (Section 4.1).
- PVC of the CBF time series, of ATT, or of aBV.
- A statistically correct vertex-level uncertainty (Section 6.5).
- Valid-observation (temporal) coverage weights (Section 6.6).
- Changing the FreeSurfer (`mri_vol2surf`) projection for fsnative/fsaverage outputs.

## 2. What the literature and other pipelines do

### 2.1 Methods

- **Linear regression (LR).** Asllani, Borogovac & Brown (2008, MRM 60:1362,
  doi:10.1002/mrm.21670) model each voxel's signal as a PV-weighted sum of tissue signals and
  solve for GM and WM by least squares over a local kernel.
  - The tissue values are assumed constant within the kernel, so the method smooths by
    construction.
  - Asllani et al. (2009, HBM 30:2927, doi:10.1002/hbm.20719) regress the mean control image
    (giving tissue M0) and the mean control−label image (giving tissue ΔM) separately, then
    quantify CBF per tissue.
  - Kernel sizes used in practice run from 3x3x1 to 11x11x1 voxels. 5x5x1 is the most common
    choice (Petr et al. 2018; ExploreASL; Zhao et al. 2017).
- **Spatial VB (BASIL).** Chappell et al. (2011, MRM 65:1173, doi:10.1002/mrm.22641) use the PV
  estimates and the different GM and WM kinetics with spatial priors.
  - Zhao et al. (2017, NeuroImage 162:384, doi:10.1016/j.neuroimage.2017.08.072) found that
    spatial VB keeps more detail.
  - LR was markedly less sensitive to noise and PV-map errors (RMSE 6.3 vs 23.4 at SNR 5).
- **Fixed GM:WM ratio.** Johnson et al. (2005, Radiology 234:851) and Du et al. (2006,
  Neurology 67:1215) divide by pGM + 0.4 pWM. This assumes a global ratio, which Chappell et
  al. (2021) criticise.
- **Modified LR.** Liang et al. (2013, MRM 69:531) use least trimmed squares to blur less.

### 2.2 Recommendations

- **Chappell et al. 2021** (NeuroImage 238:118236, doi:10.1016/j.neuroimage.2021.118236):
  - Analyse with and without PVC, and report both.
  - Prefer spatially regularized methods.
  - Bring high-resolution PV estimates to ASL resolution by summing tissue within each voxel,
    not by interpolation.
  - The PSF is usually ignored, which it should not be.
  - Use pGM > 0.7 as a pragmatic "pure GM" threshold.
  - Report the PV algorithm, the resampling, the registration and the thresholds, and share
    the PV maps.
- **Woods et al. 2024** (MRM 92:469, doi:10.1002/mrm.30091) encourage PVC as an *additional*
  analysis. They note that multi-timepoint data suit PVC methods that use the kinetics.
- **Suzuki et al. 2024** (ASL lexicon, MRM 91:1743, doi:10.1002/mrm.29815) define CBF_GM and
  CBF_WM, estimated either above a stated PV threshold or by explicit PVC.
- **Petr et al. 2018** (MAGMA 31:725, doi:10.1007/s10334-018-0691-y):
  - PV maps must match the *effective* resolution.
  - Measured effective FWHM: 3.0x3.2x7.4 mm for 2D EPI (3x3x7 mm acquired), and
    4.6x4.3x11.8 mm for one GE 3D spiral protocol (3.8x3.8x4 mm acquired).
  - Rigid registration of distorted EPI biased GM CBF by about 8%.
  - LR with a 5x5x1 kernel beat GM thresholding regardless of PV-map quality.
- **Vidorreta et al. 2013** (NeuroImage 66:662, doi:10.1016/j.neuroimage.2012.10.087) simulated a
  slice-direction PSF FWHM of 1.38 voxels for *single-shot* 3D GRASE.

### 2.3 Implementations

| Pipeline | PV maps to ASL space | PVC | Thresholds |
|---|---|---|---|
| oxford_asl / oxasl | FAST or user PVEs, PV < 0.1 zeroed; no PSF model | spatial VB with a forced spatial prior | reported means pGM >= 0.8, pWM >= 0.9 |
| ExploreASL | pGM/pWM smoothed to a sequence-specific effective resolution, then downsampled (`xASL_wrp_PreparePV`, `xASL_init_DefaultEffectiveResolution`) | native LR, off by default (`xASL_im_PVCkernel`); see below | voxel PV 0.01, kernel PV sum 0.1; outputs masked at pGM > 0.1 and pGM + pWM > 0.3 |
| HCP-ASL (Kirk et al. 2025, doi:10.1162/imag_a_00444) | surface/aseg-derived PVs at super-resolution | BASIL spatial VB | n/a |

ExploreASL's native LR uses a flat 5x5x1-voxel kernel or a Gaussian kernel of 10x10x4 mm FWHM,
solved with an SVD pseudo-inverse. Below rank 2 it falls back to a single tissue. A negative
coefficient is clipped to 0 and the other tissue is re-fitted.

ExploreASL's relative effective-resolution factors are:

- 2D EPI: [1, 1, 1].
- 3D GRASE: [1.1, 1.1, 1.38].
- 3D spiral: [1.13, 1.16, 2.53]. The code attributes these to Petr 2018, but they do not match
  that paper's published 4.6x4.3x11.8 mm.

### 2.4 Surface projection

- **fMRIPrep / ASLPrep (fsLR and CIFTI).**
  - Workbench `-volume-to-surface-mapping -ribbon-constrained` weights each voxel by its
    geometric overlap with the vertex's white-to-pial column.
  - `--project-goodvoxels` adds a *binary* `-volume-roi`.
  - fsnative/fsaverage outputs use `mri_vol2surf`, an unweighted average of 6 points across the
    ribbon.
- **HCP-ASL** (`scripts/VolumetoSurfaceASL.sh`) projects with `-volume-roi <1/variance>
  -weighted`, then writes `-bad-vertices-out`, dilates 10 mm, resamples and smooths 2 mm. It
  runs separately for PVC and non-PVC CBF and ATT.

## 3. Current ASLPrep behaviour this design builds on

All references are to this branch. Each item below was checked against the source during
review.

**CBF grid.** `mean_cbf` is computed on the native ASL reference grid ("aslref"), from
`asl_native_wf`'s output: the series with head motion correction and SDC applied, resampled
once (`aslprep/workflows/asl/base.py:429-490`).

**Quantification.**
- `ExtractCBF` smooths the M0 (or control-as-M0) image with `--smooth-kernel`
  (`aslprep/interfaces/cbf.py:163-202`). Only the smoothed image is output (`m0_file`, `:112`).
- `ComputeCBF` scales M0, corrects it for short TR, zeroes low M0 (`:507-522`), and divides ΔM
  by the scaled M0 voxel by voxel (`:737`).

**Tissue maps in ASL space.**
- `t1w_tpms` ([GM, WM, CSF]) are warped to aslref with `ApplyTransforms(interpolation='Gaussian')`
  through the inverse `aslref2anat_xfm` (`aslprep/workflows/asl/cbf.py:344-420`).
- They feed BASIL (`--pvgm/--pvwm`) and SCORE/SCRUB.

**BASIL PVC.**
- `desc-basilGM_cbf` and `desc-basilWM_cbf`, plus `...GMSD` and `...WMSD` (#726), come from
  `BASILCBF` (`aslprep/interfaces/cbf.py:1049-1128`).
- SD maps are zero outside BASIL's mask (`:1033-1044`).
- In `--sloppy` mode the PVC maps are copies of the uncorrected BASIL map.

**Derivative plumbing.** None of these is generic over `desc`:
- `BASE_INPUT_FIELDS` (`aslprep/workflows/asl/outputs.py:139-217`, indexed at `:737-760`).
- `cbf_3d_derivs` and `att_derivs` (`base.py:274-320`).
- The CBF confounds loop over `cbf_3d_derivs` (`base.py:605`) into explicitly enumerated fields
  (`aslprep/workflows/asl/confounds.py:559-566`).
- The reporting inputnode (`aslprep/workflows/asl/plotting.py:61-79`).
- Sink patterns (`aslprep/data/aslprep_bids_config.json:41-44`):
  - The perf NIfTI pattern allows only the suffixes
    `asl|aslref|abat|abv|att|cbf|coverage|mask` and has no `label` entity.
  - The GIFTI and CIFTI patterns allow only `asl|abat|abv|att|cbf`.

**Surface paths.**
- fsnative/fsaverage: `init_asl_surf_wf` (`aslprep/workflows/asl/resampling.py:40-300`).
- fsLR GIFTI: `init_asl_wb_surf_wf` (`resampling.py:303-562`), which calls fMRIPrep's
  `init_wb_vol_surf_wf`.
- CIFTI: `init_ds_ciftis_wf` (`aslprep/workflows/asl/outputs.py:1040-1233`), which calls
  `init_bold_fsLR_resampling_wf`.
- Both fMRIPrep 25.2.6 workflows build a node `volume_to_surface =
  VolumeToSurfaceMapping(method='ribbon-constrained')` (`fmriprep/workflows/bold/resampling.py:614,
  970`). The interface exposes `weighted` (requires `volume_roi`).
- The CIFTI path then dilates 10 mm, masks, and resamples with `ADAP_BARY_AREA`.
- Derivatives are warped to the anatomical grid with Lanczos interpolation (`resampling.py:482`,
  `outputs.py:1139`).

**Test truth.**
- aslscan fixtures provide `desc-pv{GM,WM,CSF}_gt` (exact box-overlap fractions on the
  acquisition grid), `desc-deltam_gt` and `desc-M0map_gt`.
- The flat phantom (`tfmni-flat`, recipe F1 `pcasl1pld`) has constant pure-tissue perfusion
  (GM 60, WM 20, CSF 0). GM and WM also differ in M0, T1, T2, T2* and ATT
  (`aslprep/tests/aslscan_fixtures.py:85-112`).
- The fixtures' anatomical probseg maps are binary 1 mm labels.
- Readouts include 2D EPI, 3D GRASE with four shots (`pcasl1pld_grase`, `NumberShots: 4`), and
  3D spiral (`pcasl1pld_ge3d`, 4x4x8 mm).

## 4. Part A: tissue PV maps on the ASL grid (`init_asl_pv_wf`)

A new workflow in `aslprep/workflows/asl/pv.py` builds the PV maps used by PVC and by PV-weighted
projection. It runs when `--pvc-method` is not `none` or `--surface-weighting` is not `none`. Every
weighting mode needs it, because the mean GM-fraction support map (§6.4) is written in all modes.

### 4.1 Signal model and what it does and does not model

Write *h* for the effective response, p_t for the high-resolution tissue fraction, and a_t for a
tissue's signal amplitude. For PVC the signal is y = q·ΔM (Section 5.1), and a_t = M0_t·f_t. The
signal is linear in the tissue amplitudes:

    y(x) = Σ_t (h ⊗ (p_t · a_t))(x) + ε(x),   t ∈ {GM, WM}.

CSF is taken to have zero perfusion, so it contributes nothing to y. If a_t is constant over the
support of *h* and of the regression kernel, the model becomes

    y(x) = Σ_t a_t · p̃_t(x),   where p̃_t = h ⊗ p_t.

The regressors are therefore p̃_t, sampled at aslref voxel centres. They are not the fractions
interpolated to voxel centres, which is what ASLPrep computes today.

**Scope of *h*.** *h* models the acquisition: integration over the voxel plus a stationary,
separable readout blur along the acquisition axes. It does **not** model the spatially varying
blur from motion-dependent and SDC resampling. That blur cannot be represented by one stationary
kernel. The omission is documented in the boilerplate and the sidecars. A user-supplied FWHM is
used as given; the help text says that it is the acquisition's effective resolution and that
preprocessing blur is not added on top. Forward-projecting tissue maps through the actual
resampling operators is a follow-up (Section 9).

### 4.2 Effective resolution

The option `--asl-effective-resolution` accepts:

- `auto` (default);
- `nominal`;
- three numbers: the **total** response FWHM in mm, given in the order (in-plane, in-plane,
  slice).

**Axes.** The slice axis is the one named in `SliceEncodingDirection`, or k if that field is
absent. The two in-plane factors in every built-in default are equal, so how the two in-plane
axes are assigned does not matter for `auto`. For user values, the first in-plane value applies
to the lower-numbered in-plane array axis.

**Voxel size.** "Voxel size" is the aslref voxel size. `AcquisitionVoxelSize` differs only when
the reconstruction interpolated. In that case it is recorded, and a warning says the defaults
assume no reconstruction interpolation.

**Kernel.**
- `nominal` is a box of one voxel. Its FWHM equals the voxel size, and it is what Chappell 2021
  recommends.
- A FWHM request *F_d* ≥ voxel size *v_d* defines a continuous kernel: the box of width *v_d*
  convolved with a Gaussian.
  - Its profile is closed form: k(x) ∝ erf((x + v/2)/(σ√2)) − erf((x − v/2)/(σ√2)).
  - σ_d is found by root-finding on that profile so that the continuous kernel's FWHM equals
    *F_d*. A request equal to *v_d* is the box alone.
  - FWHM is a property of this continuous kernel. The implemented discrete kernel is its exact
    integral over each subvoxel bin (from the erf antiderivative), so the discretization is the
    only remaining approximation. It is bounded by test 7 (§7.1), not by a separate discrete
    FWHM, which would depend on how one interpolates.
  - The kernel is truncated at ±(v/2 + 4σ) and renormalized. The slab halo (§4.3) is that
    half-width rounded up to whole subvoxels.
- Requests below *v_d* are rejected for user values. A response narrower than the voxel is not
  physical here.

**`auto` resolves as follows** (values are relative FWHM factors):

| Metadata | Factors (in-plane, in-plane, slice) | Source label |
|---|---|---|
| `MRAcquisitionType == '2D'` | nominal | `auto:2d` |
| `MRAcquisitionType == '3D'`, `PulseSequenceType` matches `/grase/i`, `NumberShots` absent or 1 | [1.1, 1.1, 1.38] (Vidorreta 2013) | `auto:grase-single-shot` |
| `MRAcquisitionType == '3D'`, `PulseSequenceType` matches `/spiral/i` | [1.17, 1.17, 2.95] (Petr 2018; in-plane = mean of 4.6/3.8 and 4.3/3.8) | `auto:spiral-petr2018` |
| anything else, including multi-shot GRASE | nominal, with a warning to set the value explicitly | `nominal-fallback` |

The spiral entry is an approximation from a single GE protocol and is labelled as such in the
warning, the sidecar and the docs. The GRASE factor is not applied to segmented GRASE because its
validity there is not established.

**Recorded provenance.** The resolved FWHM, the source label, the metadata fields used, and the
axis mapping go into the PV maps' sidecars and the boilerplate.

### 4.3 Construction

The source maps are sMRIPrep's `probseg` maps. They are treated as samples of a continuous
tissue-probability field at T1w voxel centres, reconstructed with linear interpolation. Outside
the T1w field of view (FOV) the field is 0.

1. **Supersampled grid.** Build a grid aligned with aslref with an odd integer subdivision *s_d*
   per axis. *s_d* is the smallest value whose subvoxel is ≤ min(0.5 × T1w voxel, 1 mm). Each
   aslref centre falls on a subvoxel centre.
2. **Warp.** Warp GM and WM onto it with the inverse `aslref2anat_xfm`, linear interpolation,
   and `cval=0`. Clip to [0, 1].
3. **Filter.** Apply the separable, bin-integrated kernel (§4.2) along each axis with
   zero-padding (`mode='constant'`).
4. **Sample.** Decimate at aslref centres.
5. **Chunk.**
   - Process in slabs of aslref slices along the slice axis, each with the §4.2 halo, so peak
     memory is bounded by the slab rather than the full supersampled volume.
   - `mem_gb` is declared from the slab size.
   - Slabs are stitched without overlap in the output.

The integration error is controlled by convergence, not asserted. Tests (§7.1.8) check:

- tripling *s_d* changes p̃ by less than 0.01 (absolute) on translated and oblique synthetic
  grids;
- slab and whole-volume results agree to 1e-6, including with the slice axis set to i or j.

p̃_CSF is not built. No regression uses it.

### 4.4 Derivatives (native ASL space only)

Space entities follow the existing native CBF derivatives.

    *_label-GM_desc-psf_probseg.nii.gz   (+ .json: FWHM, source label, axes, kernel, s_d)
    *_label-WM_desc-psf_probseg.nii.gz

The perf NIfTI sink pattern gains `[_label-{label}]` before `[_desc-{desc}]` and the suffix
`probseg`. These maps are grid-specific, so they are not resampled to other spaces.

## 5. Part A: local-regression PVC (`LocalRegressionPVC`)

### 5.1 Scope and estimand

**Single-PLD only.** LR PVC runs only for single-PLD PCASL/PASL data with label/control or ΔM
volumes. For multi-PLD data or scanner-provided CBF, `parse_args` and the workflow raise an
error that points to `--basil`.

- A kinetic fit of a mixed signal is not a linear mixture of the tissue fits, because GM and WM
  ATT differ.
- Scanner CBF has unknown calibration.

Tissue-specific kinetic fitting is what BASIL already does.

**What `ComputeCBF` does.** For single-PLD data it computes

    CBF_k = q_k · ΔM_k / m̄_k        (aslprep/interfaces/cbf.py:735-744)

where:

- q_k is `perfusion_factor`, which carries λ, labelling efficiency and exp(PLD_k/T1_blood), so it
  depends on the slice through the slice-timing PLD shift;
- m̄_k is `scaled_m0data`: M0 smoothed, short-TR-corrected, scaled, with `zero_low_m0` applied
  (`:512-522`);
- non-finite results are then set to 0.

CBF is not a linear mixture of tissue CBF over a neighbourhood, because m̄ varies across it. The
CBF/s inversion proposed in revision 2 is not well defined either, because the zeroed voxels
cannot be inverted.

**Regress the calibration-free signal y = q·ΔM.**

    y_k = q_k · ΔM_k = Σ_t p̃_t,k · d_t,k

What q removes, and what it does not:

- q_k is ASLPrep's *estimation* model. It uses the simplified single-compartment blood-T1
  equation, with each slice's own PLD.
- The physical tissue signal also depends on tissue T1 (T1′) and on ATT. These enter through the
  slice's acquisition time and do not cancel against q. Within a slice they are the same for
  every voxel of a tissue, so d_t is constant across a slice's in-plane neighbourhood. Across
  slices of a slice-timed acquisition it varies.
- Measured on the F1 fixture, single-slice fits are exact (relative residual about 3e-8). One
  coefficient pair across all slices leaves a 3.4% residual, and d_GM differs between slices 5
  and 30 by about 12% (revision-3 review, §11).

So:

- With `SliceTiming` that varies across slices, the kernel's slice extent must be 1. A larger
  value is rejected at workflow construction, with an explanation.
- 3D acquisitions, and 2D data without slice timing, may use through-plane kernels.
- Between slices the residual tissue-kinetic variation is model error. The sidecar states the
  assumption that d_t is constant within the kernel.
- No M0 enters the regression.

`ComputeCBF` gains three internal (non-derivative) outputs, taken from arrays it already holds:

| Output | Content |
|---|---|
| `mean_qdeltam` | mean over volumes of q_k·ΔM_k, before any sanitization |
| `m0_calib` | m̄, exactly the denominator used for CBF |
| `quant_valid` | mask where m̄ is finite and > 0 and q is finite |

**Calibrate with the centre voxel's M0.**

    CBF_t(c) = d_t / m̄_c.

- This is the same M0 and the same convention as the uncorrected CBF at *c*.
- What this gives, kept separate:
  - **Algebraically**, with noiseless data, a kernel of a single tissue t gives
    CBF_t(c) = y_c/m̄_c, which is the uncorrected CBF at *c*.
  - **With noise**, LR estimates a neighbourhood coefficient, so it agrees with uncorrected CBF
    only in expectation.
  - **Physiologically**, neither equals f_t unless m̄_c = M0_t. Because m̄ is smoothed
    (`--smooth-kernel`), that requires the calibration's own support, not just the regression
    kernel, to be a single tissue.
- In mixed voxels m̄_c is the mixed M0, as in any voxelwise-calibrated ASL map. ASLPrep runs
  oxford_asl with `-c <M0>` and `--cgain`. Whether that configuration uses voxelwise
  calibration is checked against oxford_asl's source during PR 1 and stated only for that
  configuration. A median LR-minus-BASIL difference is reported but is not taken as proof that
  the conventions match. The boilerplate states the convention.
- Tissue-specific M0 calibration (Asllani 2009) is a follow-up. A tissue M0 fit is not
  identifiable in single-tissue neighbourhoods, and dropping a collinear CSF column biases it,
  as the revision-2 review showed (§11).
- With `M0Type: Estimate`, m̄ is the supplied scalar after scaling, and the formula applies
  unchanged.
- Control-derived M0 without background suppression likewise applies unchanged.
- Background-suppressed data with `M0Type: Absent` is already rejected upstream.

**Which maps.**

- `mean_cbf`: y = `mean_qdeltam`.
- `mean_cbf_score`: SCORE's mean of retained volumes, with y = CBF_score·m̄ on `quant_valid`
  voxels. This is exact because q and m̄ are time-invariant. Test 6 asserts that SCORE's output is
  the retained-volume mean.
- **SCRUB is deferred.** Its estimator shrinks toward a global GM/WM prior fitted in CBF space
  (`aslprep/utils/cbf.py:371, 407, 469`), which is not a constant-coefficient PV mixture after any
  rescaling.

### 5.2 Estimator and numerical contract

For each voxel *c* in the brain mask, let K(c) be the box of `--pvc-kernel-size` voxels centred
on *c*. The size is given as (in-plane, in-plane, slice) with odd values; the default is
`5 5 1`. The valid set is

    V(c) = {k ∈ K(c) : in mask, quant_valid_k, y_k finite,
            p̃_GM,k and p̃_WM,k finite, p̃_GM,k + p̃_WM,k > 0}.

A genuine y_k = 0 is a valid observation. Invalid voxels are identified by `quant_valid`, not by
the value 0.

All valid voxels get uniform weight. With X the |V|×2 design,

    G = XᵀX = Σ_V p̃ p̃ᵀ,   b = Xᵀy,   β = G⁻¹ b,
    RSS = Σ_V y² − 2βᵀb + βᵀGβ.

Every sum is a box-filtered image, computed under this contract:

- Arrays are built with `np.where(valid, value, 0.0)`. Multiplying by a mask is not allowed,
  because 0·NaN is NaN.
- The same `valid` array enters every term: G, b, the count, the supports, and Σy².
- Accumulation is float64.
- Box filters use `scipy.ndimage.uniform_filter(..., mode='constant', cval=0)`, scaled by the
  kernel volume to give sums. There is no reflection at image edges.
- RSS is **computed directly**, not from the expansion identity.
  - Loop over the kernel offsets (25 for `5 5 1`), shifting y and p̃ with zero padding.
  - Accumulate valid · (y_{c+o} − p̃_{c+o}ᵀβ_c)² per offset.
  - This avoids cancellation and never discards a real positive residual.
  - It is still vectorized across voxels.
- The 2x2 solve is closed form, with det(G) ≤ 0 treated as singular.

- The solve is applied only where the status (§5.3) selects that fit. Singular matrices never
  reach a solver: the 2x2 and 1x1 solves are closed form and evaluated under the status mask.

There is no intercept and no non-negativity constraint. Clipping biases mean GM CBF upward, while
ASLPrep keeps negative CBF elsewhere and QCs it. Negative estimates are kept and counted.

### 5.3 Status, support and conditioning

The quantities used are:

- S_t = Σ_V p̃_t, the kernel's tissue volume in voxel units;
- the normalized Gram matrix Ĝ = D^{-1/2} G D^{-1/2}, with D = diag(G), and its condition number
  κ = κ₂(Ĝ) (equal to κ₂ of the column-normalized design, squared);
- the contamination ratio c_GM = Σ p̃_GM p̃_WM / Σ p̃_GM², which is the factor by which WM perfusion
  leaks into a GM-only fit (c_WM likewise).

Each voxel is assigned the first status, in this order, whose conditions hold. The guards come
first, so a fit is never attempted for an invalid centre:

| Order | Status | Fit | Conditions |
|---|---|---|---|
| 1 | 0 | none | centre outside the brain mask |
| 2 | 4 | none | centre not `quant_valid` (no calibration) |
| 3 | 1 | two-tissue | \|V\| ≥ 5, S_GM ≥ τ, S_WM ≥ τ, D > 0, κ ≤ κ_max |
| 4 | 2 | GM only | \|V\| ≥ 4, S_GM ≥ τ, c_GM ≤ c_max |
| 5 | 3 | WM only | \|V\| ≥ 4, S_WM ≥ τ, c_WM ≤ c_max |
| 6 | 4 | none | otherwise (not identifiable) |

**Defaults** (interface inputs, documented, not CLI flags in this change):

- τ = 0.5;
- κ_max = 100, equivalent to κ₂(X̂) ≤ 10;
- c_max = 0.05, so the leak is at most 5% of the other tissue's value.

The minimum counts leave at least three residual degrees of freedom. The singular case is the
D > 0 check: a zero column fails it.

**The 50:50 counterexample stays unidentifiable.** With p̃_GM = p̃_WM = 0.5 everywhere, κ is
infinite and c = 1, so the voxel gets status 4. A single-tissue fit would have reported 80 for
true GM 60 and WM 20.

**Output masking** follows oxford_asl and ExploreASL:

- CBF_GM is written where p̃_GM(c) ≥ 0.1 and the status is 1 or 2.
- CBF_WM is written where p̃_WM(c) ≥ 0.1 and the status is 1 or 3.
- Everything else is 0, matching how ASLPrep stores CBF outside the mask.

The status map distinguishes "not estimated" from a true 0.

### 5.4 Diagnostics and derivatives

Under i.i.d. noise of variance σ², Var(β̂) = σ² G⁻¹.

- **Noise amplification.** a_t = sqrt([G⁻¹]_tt), or 1/sqrt(Σ p̃_t²) for single-tissue fits.
- **Residual RMS.** r = sqrt(RSS / (|V| − n_params)).
- **Standard error.** se(CBF_t) = a_t · r / m̄_c.
  - It is **conditional on the calibration m̄_c**, which the uncorrected CBF also treats as
    known. No fitted denominator is involved.
  - `zero_low_m0` already removes near-zero m̄.
  - The sidecars say "conditional on M0".
- **Spatial correlation.** Neighbouring voxels are correlated by the PSF. The true covariance is
  G⁻¹XᵀΣXG⁻¹, which can be larger or smaller per coefficient. The sidecars call se an "i.i.d.
  approximation", and #682 tracks a better estimate.

Derivatives:

| File | Where it goes |
|---|---|
| `*_desc-lrGM_cbf.nii.gz`, `*_desc-lrWM_cbf.nii.gz` | native + resampled (normalized, §6.2) + parcellated with their own validity mask |
| `*_desc-lrGMSE_cbf.nii.gz`, `*_desc-lrWMSE_cbf.nii.gz` | **native only** |
| `*_desc-lr_pvstatus.nii.gz` (+ .json `Levels`) | native only |
| `*_desc-lr_pvcondition.nii.gz` | native only |
| `*_desc-lrGM_pvnoise.nii.gz`, `*_desc-lrWM_pvnoise.nii.gz` | native only |
| `*_desc-lr_pvrmse.nii.gz` | native only |

With `--scorescrub`, the SCORE map adds `desc-scoreLRGM`, `desc-scoreLRWM` and their SE, status and
diagnostics. SCRUB is deferred (§5.1).

**SE stays native.** Interpolating an SE image does not give the SE of the interpolated estimate:
two independent estimates with SE 2, averaged equally, have SE √2, not 2. LR estimates also
share neighbourhoods. The SE maps are therefore written in native space only, where they are the
i.i.d. SE of the estimate stored in that voxel. Their surface use is described in §6.2.

**Parcellation uses per-map validity.** Without it, failed estimates stored as 0 enter parcel
means: five valid 60s and five failures give 30. `init_parcellate_cbf_wf` therefore passes each
LR map's own validity mask to `NiftiParcellate`. That mask is the status-and-masking rule of §5.3
for that tissue, intersected with `asl_mask` and warped with the map. The coverage TSV then
reports estimation coverage. The existing maps keep `asl_mask`. A test mixes genuine zero
estimates with failed ones.

**Plumbing changes, all enumerated:**

- New `BASE_INPUT_FIELDS` entries (`outputs.py`).
- New `cbf_wf` outputnode fields (`cbf.py`).
- LR maps added to `cbf_3d_derivs`, with a normalized resampling branch (`base.py`).
- SE maps and diagnostics go in a new native-only list, written by a dedicated sink loop.
- Per-map validity masks reach `init_parcellate_cbf_wf` (`cbf.py:820-835`).
- The `Resampling` boilerplate (`base.py:332`) currently says every volumetric resampling uses
  Lanczos. It becomes conditional: it names the normalized linear scheme when LR maps or
  weighted projection are present.
- The LR maps are **excluded** from the CBF confounds loop at `base.py:605`. QEI and
  negative-CBF QC for LR maps are a follow-up, so `confounds.py`'s enumerated fields stay as
  they are.
- New reporting inputnode fields (`plotting.py`).
- New entries in `reports-spec.yml` and `reports-spec-asl.yml`.
- In `aslprep_bids_config.json`, the perf NIfTI pattern gains the ASLPrep-specific suffixes
  `pvstatus|pvcondition|pvnoise|pvrmse|probseg`. These are project-specific derivatives, not
  BIDS-standard suffixes, and the docs say so.

Every sidecar records the method, the kernel, the thresholds, the calibration convention, the PV
map files, the resolved FWHM and its source.

### 5.5 CLI and config

```text
--pvc-method {none,local-regression}              default none -> config.workflow.pvc_method
--asl-effective-resolution {auto,nominal | F F F} default auto -> config.workflow.asl_effective_resolution
--pvc-kernel-size N N N                            default 5 5 1 -> config.workflow.pvc_kernel_size
```

- All three go in `g_conf`, kebab-case.
- `parse_args` validates:
  - odd, positive kernel sizes;
  - a positive FWHM triple or a keyword;
  - `local-regression` is incompatible with multi-PLD data and `--force`d CBF-only processing.
    This check happens when the data are known, in workflow construction.
- `--pvc-method` is separate from `--basil`. BASIL PVC stays tied to `--basil`, and the help text
  says so. Adding `basil` as a value would change what `--basil` means, so it is a follow-up.
- `--pvc-tissues` from #681 is not added. A free CSF *perfusion* regressor adds a parameter whose
  expected value is zero, and it worsens conditioning. No M0 is regressed (§5.1), so CSF does
  not enter anywhere.
- LR PVC runs in `--sloppy` mode too, because it is cheap.

### 5.6 Reports

- `desc-lrGM_cbf.svg`: a `CBFSummaryPlot` like the one for `basilGM`.
- `desc-lr_pvc.svg` (the `pvc` figure suffix already exists): one figure with four panels:
  - p̃_GM next to the current Gaussian-interpolated pGM;
  - the status map;
  - a_GM;
  - uncorrected-vs-corrected histograms within p̃_GM ≥ 0.1.
- Text in `CBFSummary`: the fraction of in-mask voxels in each status, the resolved FWHM and its
  source, and the calibration convention.
- Boilerplate:
  - cites Asllani 2008/2009, Petr 2018 and Chappell 2021;
  - states that corrected values are an additional analysis;
  - states that motion and SDC resampling blur is not modelled;
  - states that the spiral and GRASE defaults are approximations.

## 6. Part B: weighted surface projection

### 6.1 Option

```text
--surface-weighting {none,pv,precision,pv-precision}   default none -> config.workflow.surface_weighting
```

- Applies only to the Workbench paths: fsLR GIFTI and CIFTI cortex.
- It is an error if no Workbench output is requested.
- CIFTI subcortical voxels come from the MNI6 volume and are not affected.
- fsnative/fsaverage `mri_vol2surf` projection is unchanged and documented as unweighted.

### 6.2 Weights

Weights are built on the same anatomical, ASL-resolution grid the Workbench node receives, using
the same `aslref2anat` transform as the data. Interpolation is **linear**, not Lanczos, so
weights and validity stay non-negative.

**Validity.** Each parent map carries a validity mask in aslref:

| Map | Valid where |
|---|---|
| internally quantified single-PLD CBF (`mean_cbf`, SCORE, SCRUB) | brain mask ∧ `quant_valid` ∧ finite |
| multi-PLD fitted CBF | brain mask ∧ finite ∧ not flagged in the fit's failure output |
| scanner-provided CBF | brain mask ∧ finite |
| LR maps | status ∈ {1, 2} for GM (1, 3 for WM) and inside the output mask |
| BASIL maps | BASIL's mask, i.e. SD > 0 and finite |

The mask is warped with linear interpolation and thresholded at 0.5. Invalid voxels get weight 0
in every mode.

**Validity-aware resampling of the data.** Weighting cannot undo contamination from resampling.
A target voxel that is 75% valid (value 60) and 25% invalid (stored 0) interpolates to 45 yet
passes the 0.5 threshold. So whenever `--surface-weighting` is not `none`, the map fed to
`volume_to_surface` is resampled by normalized interpolation:

    x' = warp(x · valid) / warp(valid)      where warp(valid) ≥ 0.5, else 0 with weight 0

- `warp` is linear interpolation with the same transform and reference.
- This replaces the Lanczos warp for the weighted surface input only. The unweighted default path
  is unchanged.
- The new LR maps use the same normalized scheme in their volumetric outputs to other spaces.
  These are new outputs, so no existing derivative changes.
- BASIL maps keep their existing volumetric resampling.

**Modes.**

- **`pv`.** w = p̃_GM (Section 4), clipped to [0, 1].
  - Ribbon overlap measures how much of a vertex's column lies in a voxel. p̃_GM measures how
    much of the voxel's signal is GM signal. These are different quantities, so multiplying them
    is not double counting. It does change the estimand: it down-weights voxels whose value is
    mostly WM or CSF signal.
  - For PVC maps, PV weighting is a heuristic, not a reliability weight.
- **`precision`.** w = 1 / max(v, v_floor). v is the native variance of *that map*
  (`basilSD²`, `basilGMSD²`, `basilWMSD²`, `lrGMSE²`, `lrWMSE²`; BASIL ATT uses
  `basilSD_att²`), resampled internally with the **same normalized scheme** as the data:
  v' = warp(v·valid)/warp(valid).
  - Without normalization, zero-filled invalid neighbours would lower v at boundaries. A valid
    variance of 4 next to an invalid 0 would become 3, so lost support would raise the precision.
  - Non-finite v marks the voxel invalid before warping.
  - Even normalized, v' is not the variance of the resampled estimate, which is diag(AΣAᵀ). The
    weight is labelled a *reliability field*, not an inverse variance, in the docs and sidecars.
  - The resampled v' is used only internally and is never written as an SE derivative.
  - v_floor is the 5th percentile of v over valid voxels with finite v > 0 inside the ribbon,
    and never below 1e-6 × the median of those voxels.
  - If no such voxel exists, the map falls back to `pv` if PV is available, or to unweighted
    otherwise.
  - Maps with no variance map use the same fallback.
  - The sidecar field `SurfaceWeighting` records the mode actually used.
- **`pv-precision`.** The product of the two, with the same fallback rules.

With `--project-goodvoxels`, the binary goodvoxels mask multiplies the weight.

### 6.3 Wiring

fMRIPrep's two workflows build an unweighted `volume_to_surface` node. ASLPrep sets
`wf.get_node('volume_to_surface').inputs.weighted = True` after constructing them, and connects
the per-map weight volume to `inputnode.volume_roi` **in place of** the goodvoxels connection,
not in addition to it. Goodvoxels has already been folded into the weight.

A unit test asserts that the node exists and that `volume_roi` has a single source in both paths.
A rename in fMRIPrep then fails loudly. Copying the workflows would cost more to maintain.

### 6.4 Per-vertex support outputs (fsLR GIFTI only)

Support maps are computed on the native surfaces with the plain ribbon-constrained mapping (no
ROI). For each vertex, write a_i for the ribbon overlap of voxel i.

**The denominator is the represented column.**

- The grid is the ASL-resolution anatomical reference from `GenerateSamplingReference`. It is
  cropped to the T1w brain-mask bounding box plus two voxels (`base.py:684`, `apply.py:167`,
  `niworkflows/interfaces/nibabel.py:627`). It does **not** span the whole T1w FOV.
- Voxels outside the ASL FOV but inside this grid have indicator 0 and count in the denominator.
  A column truncated by the ASL FOV therefore gets low coverage.
- Geometry outside the grid cannot count. Workbench builds each vertex's sampling polyhedron
  from the vertex and its neighbours on both surfaces (`niworkflows/interfaces/workbench.py:538`).
  That polyhedron lies inside the convex hull of those points, and the grid's bounding box is
  convex. So a vertex's column is fully represented exactly when every point of its 1-ring, on the
  white and the pial surface, lies inside the grid box.
- The check flags every vertex whose 1-ring has a point outside the grid. It is computed on the
  native surfaces, resampled to fsLR like the other support maps, thresholded at 0.5, and added
  to the insufficient-support mask. ASLPrep warns when any vertex is flagged.

Maps:

- **Coverage:** Σ a_i·valid_i / Σ a_i, from projecting the validity indicator.
- **Mean weight:** Σ a_i w_i / Σ a_i, from projecting the unnormalized weight.
- **Mean GM fraction:** Σ a_i p̃_GM,i / Σ a_i, from projecting p̃_GM. This is #676's "mean
  partial-volume fraction", and it is reported in every mode.
- **Effective weighted fraction:** (Σ a_i w_i)² / (Σ a_i · Σ a_i w_i²), from the mean weight and
  a projection of w².
  - It lies in [0, 1] and is defined as 0 when Σ a_i w_i = 0.
  - Low values mean a few dominant voxels.
  - A value of 1 means the weight is uniform over the represented column. That includes one voxel
    covering the whole column, which is full support at ASL resolution.
  - It is not an independent-sample count.

**Transport to fsLR.** The support maps are **not dilated**. They are resampled to fsLR with
`ADAP_BARY_AREA` and masked to the atlas cortex.

**Insufficient-support mask.** A vertex is masked if any of these holds:

- resampled coverage < 0.5;
- it is flagged by `-bad-vertices-out` (resampled and thresholded at 0.5);
- it is flagged by the out-of-grid check above.

The perfusion values there are left as
they are (dilated, as now) and flagged rather than set to NaN, so existing consumers do not
break.

Outputs (GIFTI only; CIFTI has no hemisphere entity and mixes in subcortical voxels):

    *_hemi-L_space-fsLR_den-32k_desc-<map>_projcoverage.func.gii
    *_hemi-L_space-fsLR_den-32k_desc-<map>_projweight.func.gii
    *_hemi-L_space-fsLR_den-32k_desc-<map>_projneff.func.gii
    *_hemi-L_space-fsLR_den-32k_desc-<map>_projpvgm.func.gii
    *_hemi-L_space-fsLR_den-32k_desc-<map>_mask.func.gii

The GIFTI sink pattern gains `projcoverage|projweight|projneff|projpvgm|mask`. When only CIFTI is
requested, these GIFTIs are still written alongside it.

### 6.5 Uncertainty caveat

A vertex value is a weighted mean of correlated voxels. Its variance needs the voxel covariance,
which neither BASIL nor LR provides, and Workbench normalizes the weights internally. The
existing BASIL SD surface outputs are labelled "projected voxel SD (weighted mean)", not a vertex
standard error, and use the same weights as their parent map. LR SE maps are not projected
(§5.4).

### 6.6 What #676 asked for that is deferred

- **Valid-observation coverage.** This is the fraction of volumes contributing per voxel. ASLPrep
  has only a **global** retained-volume fraction, from SCORE's outlier index
  (`aslprep/utils/cbf.py:344`), which is the same for every voxel. It could be recorded in the
  sidecars, but it carries no spatial information for weighting. Voxelwise valid-observation
  coverage does not exist yet, so it is deferred.
- **Distance-from-surface weights.** Workbench offers distinct distance weighting
  (`-gaussian`, `-thin-columns`). Choosing among them is a scope decision for a later change. It
  is not a duplicate of ribbon overlap.
- **Total contributing weight.** Workbench's ribbon mapping normalizes by Σ a_i, so Σ a_i w_i
  cannot be recovered from its outputs without per-vertex weight dumps. Mean weight is provided
  instead, and is not equivalent.

Covered: spatial coverage, the effective weighted fraction, mean weight, mean GM fraction, the
support mask, precision with a floor, and goodvoxels.

### 6.7 Reports

- One reportlet per weighted map: the weighted projection, coverage, mean weight, and the
  difference from the unweighted projection. The unweighted projection is computed in the
  working directory only, for the report.
- Text: the fraction of cortical vertices with insufficient support, and the weighting mode
  actually used.

## 7. Validation

### 7.1 Unit tests (no external data)

`aslprep/tests/test_utils_pvc.py` and `test_workflows_pvc.py`:

1. **Ideal recovery (exact-model test).** This tests the estimator under its own model. It is not
   evidence for the physical premise, which 7.2(a) measures.
   - Synthetic 1 mm phantom: smooth GM/WM boundaries, d_GM and d_WM constant per slice.
   - Forward-model y with a known kernel, decimate to 4 mm, and build p̃ with the same kernel.
   - Noiseless LR (`5 5 1`) must recover d_t to within 0.5% in status-1 voxels.
   - With uniform m̄ equal to M0_t, CBF_t equals f_t. That is the physiological case, with
     homogeneous calibration support.
   - **Algebraic case.** With non-uniform m̄, noiseless single-tissue kernels give
     CBF_t = y_c/m̄_c exactly.
   - **Noise case.** Over noise realisations, the mean of CBF_t matches the expected uncorrected
     CBF to within 2%.
   - **Slice timing.** With `SliceTiming` varying, a `5 5 3` kernel is rejected.
   - **3D without slice timing.** A `5 5 3` kernel recovers d_t to within 0.5%.
2. **PSF mismatch.** The same data with `nominal` p̃ gives a GM error larger than 2%. The test
   reports the size and sign for this phantom and claims no general direction.
3. **Identifiability.**
   - The 50:50 phantom gets status 4 and never reports 80.
   - Pure GM in the whole kernel gets status 2.
   - Empty support gets status 4.
   - An invalid centre (`quant_valid` false) surrounded by a well-conditioned, valid two-tissue
     neighbourhood gets status 4, and no division happens.
   - Status precedence is deterministic.
4. **Noise amplification.** Over 500 i.i.d. realisations added after blurring, the empirical SD of
   d̂_GM matches σ·a_GM to within 10%.
5. **Box-filter implementation.** It matches a brute-force per-voxel `lstsq` on a random volume
   that includes NaNs, genuine zeros, image corners, out-of-mask voxels and a 1e6 dynamic range.
   - β matches to a relative 1e-10.
   - RSS, r and SE match to a scale-aware tolerance (1e-10 × Σy² for RSS).
   - A near-exact fit with ±0.001 residuals around 60 keeps its positive RSS.
6. **Quantification plumbing.**
   - `mean_qdeltam` equals `mean_cbf`·`m0_calib` on `quant_valid` voxels, for PCASL, PASL, slice
     timing on and off, short-TR M0, and `M0Type` Separate, Estimate and control-derived.
   - Voxels zeroed by `zero_low_m0` are outside `quant_valid`, and a centre there gets status 4.
   - SCORE's mean equals the mean of its retained volumes.
   - `d_t/m̄_c` reproduces `ComputeCBF` for a single-tissue neighbourhood.
7. **Effective resolution.**
   - Each metadata case maps to the stated factors and source labels, including multi-shot
     GRASE → `nominal-fallback`.
   - The continuous box⊗Gaussian kernel's FWHM is within 0.1% of the request: 1.1, 1.21, 1.38
     and 2.95 voxels.
   - The bin-integrated discrete kernel, at the *s_d* used, differs from the continuous kernel
     (averaged over each bin) by less than 1e-3 in L1.
   - Truncation removes less than 1e-4 of the mass.
   - Values below the voxel size are rejected.
   - `SliceEncodingDirection` selects the slice axis for both the PSF and the LR kernel.
8. **PV-map construction.**
   - With translated and oblique synthetic grids, `nominal` p̃ matches an independent
     high-resolution reference integration to within 0.01.
   - Tripling *s_d* changes p̃ by less than 0.01.
   - Voxels at the FOV edge are zero-padded.
   - Slab and whole-volume results agree to 1e-6, with the slice axis set to i, j and k.
9. **Surface weights.**
   - The floor ignores zero and invalid variances.
   - An all-zero variance map falls back to `pv` (or unweighted).
   - Invalid voxels get weight 0.
   - `pv-precision` is the product.
   - Goodvoxels zeroes weights.
   - Normalized resampling of a constant 60 field with constant variance 4, next to invalid
     zeros, keeps both 60 and 4 wherever warp(valid) ≥ 0.5. Lanczos and plain linear
     interpolation do not.
   - Uncorrected CBF voxels outside `quant_valid` get weight 0.
   - `--surface-weighting precision` with `--pvc-method none` builds `init_asl_pv_wf` and writes
     `projpvgm` (graph test).
   - Parcellating an LR map that mixes genuine zero estimates with failed ones averages only the
     valid voxels and reports estimation coverage.
10. **Surface wiring.** In both paths and with goodvoxels on and off:
    - `volume_to_surface` exists and has `weighted=True`;
    - `volume_roi` has exactly one source;
    - BASIL ATT pairs with `basilSD_att`.
11. **Support maps.** On a synthetic ribbon with Workbench available (skipped otherwise):
    - a column truncated by the FOV gets coverage below 1;
    - a single dominant voxel gives a low effective fraction;
    - zero total weight gives an effective fraction of 0 and a masked vertex;
    - a surface vertex outside the reference grid is masked, with a warning;
    - a vertex inside the grid whose 1-ring neighbour is outside is also masked;
    - support maps are not dilated;
    - the medial wall is excluded.
12. **Default path unchanged.** Workflow construction without the new flags gives the same node
    set and connections as before (a graph comparison). No new sinks are added.
13. **Filenames.** Every new output and sidecar builds a valid path through ASLPrep's
    `DerivativesDataSink` config.
14. **CLI and config.** The options round-trip through ToML. Invalid kernels and resolutions are
    rejected. Multi-PLD and CBF-only data with `local-regression` raise an error.

### 7.2 aslscan integration (report-only first)

Each error source is scored separately, so a change in one is not read as a change in another.

**(a) Oracle estimator on noiseless truth.** This is the estimator alone, with no preprocessing
and no PV workflow. It runs in the unit-test job, where fixtures are available.

- **Signal.**
  - `desc-deltam_gt` is 4D and aligned with `aslcontext`. It is non-zero only in label rows.
  - Average **only** the rows whose context is `label` (or `deltam`), with the sign convention
    that the fixture's truth sidecar states.
  - Multiply by q, computed from the fixture's `asl.json` through the same helper as
    `perfusion_factor`, slice timing included.
- **Regressors.** `desc-pv{GM,WM}_gt`, the exact regressors for the truth geometry (box overlap,
  no readout PSF).
- **Target.** The flat phantom's tissue amplitudes, fitted **per slice**: d_t,s from a two-tissue
  least-squares fit over all voxels of slice *s* with pv_GM + pv_WM > 0.
  - Tissue T1′ and ATT make d_t vary between slices (§5.1), so a whole-volume fit is not the
    target.
  - The revision-3 review measured per-slice relative residuals of about 3e-8 on F1, which
    confirms the in-slice linear-mixture premise.
  - Each slice's residual is reported as a running check of that premise.
- **Assertion.** Local LR (`5 5 1`, `nominal`) matches d_t,s for its slice to within 1% in status-1
  voxels.
- **Report only.**
  - **Kinetic error:** the spread of d_t,s across slices.
  - **Calibration error:** d_t,s divided by the phantom's tissue M0 constants, compared with 60
    and 20. This measures the q model (blood-T1 simplification, ATT) and λ, not unmixing.

**(b) End to end, F1** (`pcasl1pld`, flat phantom, 2D EPI). Add `--pvc-method local-regression`.
The conditional pure-tissue truth is

    f_t^true(c) = box(p_t · f_t)(c) / box(p_t)(c)

with p_t and f_t from the 1 mm phantom and box the ASL voxel. It is undefined where
box(p_t) < 0.05. For the flat phantom it is exactly 60 and 20.

Reported metrics:

- **Physiological error.** `desc-lrGM` and `desc-lrWM` against f_t^true in status-1 voxels.
  This includes the mixed-M0 calibration convention.
- **PV-map error.** The written `label-GM_desc-psf_probseg` (`nominal` for 2D) against
  `desc-pvGM_gt`: median absolute difference. Registration and segmentation enter here.
- **Preprocessing error.** LR rerun on ASLPrep's `mean_qdeltam` with `desc-pv*_gt` as regressors,
  compared with (a).
- **Status share.** The fraction of status-1 voxels.
- **LR vs BASIL difference.** The median difference between `lrGM` and `basilGM`. This is
  descriptive, not evidence that the calibration conventions match (§5.1).

**(c) Readouts** (`pcasl1pld_ge3d`, `pcasl1pld_grase`).

- Report the GM ratio to f_t^true for `auto` and for `nominal`, together with the (a)-style oracle
  error and the local-constancy error. Local-constancy error is LR on the noiseless truth with
  exact regressors versus f_t^true, and it is non-zero for the cosine-modulated phantom.
- The auto-versus-nominal difference is then read *net* of those two terms. This measures whether
  the defaults suit aslscan's simulated readouts; it does not assume they do.
- The four-shot GRASE fixture resolves to `nominal-fallback`, so it is also run with an explicit
  `--asl-effective-resolution`.
- f_t^true for these non-flat phantoms needs `desc-perfusion{GM,WM}_gt`, computed as above in the
  fixture assembler. That changes the fixture spec, so `.circleci/aslscan_fixtures.txt` is
  regenerated.

**Ceilings.** None until the first CI measurements, the same policy as BASIL in the aslscan design.
They are then set in `truth_bounds.py`.

**Part B.** There is no surface truth in the fixtures, so Part B is covered by unit tests 9-11 and
a smoke run in one existing CIFTI-producing integration test. Surface truth is a follow-up.

## 8. Risks

| Risk | Mitigation |
|---|---|
| Effective-resolution defaults wrong for a protocol | Limited to single-shot GRASE and spiral, labelled as approximations; fallback to nominal otherwise; measured in 7.2(c); user override |
| Motion and SDC blur not modelled | Stated in boilerplate and sidecars; follow-up for forward-projected PV maps |
| PV-map errors from segmentation or registration (~8% bias, Petr 2018) | LR is the method least sensitive to PV-map error (Zhao 2017); the docs advise SDC; PV maps are shared |
| Noise amplification near boundaries | status precedence with contamination bound; condition and noise maps; output masks |
| Mixed-voxel M0 calibration differs from tissue-specific calibration | same convention as uncorrected CBF; stated in boilerplate; measured in 7.2(a)-(b); tissue M0 is a follow-up |
| Users treat LR maps as a replacement for uncorrected CBF | always additional; boilerplate and docs (Chappell 2021, Woods 2024) |
| fMRIPrep renames the `volume_to_surface` node | unit test 10 |
| Supersampling memory | slab processing; `mem_gb` declared |
| Precision weights dominated by a few voxels | validity-aware floor; effective-fraction map exposes it |

## 9. Follow-ups (out of scope)

- Feed BASIL the p̃ maps behind an option, and compare in the aslscan suite.
- Tissue-specific M0 calibration (Asllani 2009), with an identifiability analysis for
  single-tissue and CSF-collinear neighbourhoods.
- LR for SCRUB, through a signal-domain robust estimator.
- Add `basil` as a `--pvc-method` value, decoupled from `--basil`.
- Forward-project tissue maps through the actual HMC and SDC resampling operators.
- Gaussian-weighted kernels, LTS or TLS regression, and a user PSF kernel image.
- QEI and negative-CBF QC for LR maps.
- Surface-based PV (Toblerone) and surface ground truth in aslscan.
- A covariance-aware vertex standard error, and valid-observation coverage.
- A data-driven effective-resolution estimate (ExploreASL `x.ResolutionEstimation`).
- Possible upstreaming of weighted projection to fMRIPrep.
- Side note: the goodvoxels coefficient of variation is computed on the raw label/control
  series (`base.py:815-846`), so label/control alternation contributes to it. This is worth
  checking separately.

## 10. Implementation phases

1. **PR 1 (#681).** `init_asl_pv_wf`, the `ComputeCBF` internal outputs,
   `LocalRegressionPVC`, CLI and config, derivatives and plumbing, reports, unit tests 1-8 and
   12-14, and 7.2(a)-(b).
2. **PR 2 (#681).** The `desc-perfusion{GM,WM}_gt` fixtures and 7.2(c).
3. **PR 3 (#676).** Surface weighting, support outputs, reports, and unit tests 9-11.

## 11. Review log

The Codex adversarial review of revision 1 raised 17 findings (16 major, 1 minor). Each was checked
against the source and the mathematics. All were accepted, wholly or in part. Changes by finding:

1. **Box⊗Gaussian FWHM arithmetic was wrong** (verified: box variance is v²/12, not
   (v/2.3548)²). σ is now solved numerically for the actual FWHM (§4.2). The test tolerance is 1%
   on the actual FWHM (§7.1.7).
2. **"Whole-chain" PSF was overclaimed.** *h* is now scoped to the acquisition. Motion and SDC
   blur are stated as unmodelled (§4.1), and forward projection is a follow-up.
3. **`auto` was over-broad.** 2D now maps to nominal. GRASE factors apply only to single-shot
   acquisitions (the four-shot fixture was the counterexample). The spiral factors are labelled
   as a single-protocol approximation, with equal in-plane factors. Axis order is defined for both
   the PSF and the LR kernel (§4.2, §5.2).
4. **Linear interpolation is not exact box integration.** The TPM is now treated as a continuous
   field. Accuracy is enforced by a convergence test with translated and oblique grids, with zero
   padding and slab processing (§4.3, §7.1.8).
5. **M0 smoothing did not justify CBF-domain LR** (verified: `cbf.py:163-202`, `:737`). LR now
   regresses ΔM (via CBF/s) and unsmoothed M0 separately, then quantifies per tissue (§5.1). The
   claim of equivalence with BASIL is replaced by a measurement (§7.2b). *Superseded by
   revision-2 findings 1-3: LR now regresses y = q·ΔM and calibrates with m̄_c.*
6. **Multi-PLD is not a linear mixture.** LR is restricted to single-PLD data. Multi-PLD and
   scanner CBF raise an error that points to BASIL (§5.1).
7. **The single-tissue fallback was biased** (the 50:50 example gives 80). A fallback is now
   allowed only when the contamination ratio is ≤ 0.05, with deterministic precedence and a test
   (§5.3, §7.1.3).
8. **Degrees of freedom and κ were undefined.** Minimum counts now apply to every fit. κ is
   defined on D^{-1/2}GD^{-1/2}, and the singular cases are handled (§5.3).
9. **The box-filter numerical contract was missing.** The contract is now explicit: `np.where`,
   one validity array, `mode='constant'`, float64, and RSS clamping (§5.2, §7.1.5).
10. **Invalid estimates could get the largest precision** (verified: BASIL SD is zero outside its
    mask, `cbf.py:1033-1044`). Validity masks were added, and the floor ignores invalid and zero
    variances, with a fallback (§6.2).
11. **Resampled SD² is not the variance of the resampled map.** Variance is now resampled with
    linear interpolation and called a reliability field (§6.2).
12. **The plumbing is not generic over `desc`** (verified: `outputs.py:737-760`, `base.py:605`,
    `confounds.py:559-566`, `plotting.py:61-79`). Every change is enumerated. LR maps are excluded
    from the confounds loop, and `volume_roi` is replaced, not added (§5.4, §6.3).
13. **Filenames would not sink** (verified: config lines 41-44). Patterns and suffixes are
    specified, and support outputs are GIFTI-only (§4.4, §5.4, §6.4, §7.1.13).
14. **Coverage was not support.** The denominator is defined over the full T1w-FOV grid, and an
    effective weighted fraction was added. A per-voxel count is not possible without
    `-output-weights` per vertex, so the fraction is used instead (§6.4). *The full-FOV premise
    was wrong; see revision-2 finding 10.*
15. **Support transport and CIFTI were unspecified.** Support maps are undilated, resampled and
    masked; they are GIFTI-only and leave subcortical voxels unaffected (§6.4).
16. **Fixture truth was ambiguous.** Validation is now split into an estimator test on noiseless
    truth, an end-to-end test against a defined conditional truth, and readout measurement
    (§7.2).
17. **Wording, and items from #676.** Directional claims were qualified (§5.4, §7.1.2). The
    default-path check is now a graph comparison rather than a "byte-for-byte" claim. The
    deferred #676 items are listed (§6.6), and ATT/`basilSD_att` pairing was added (§6.2).

### Revision 2 review (13 findings: 9 major, 4 minor)

The Codex adversarial review of revision 2 focused on what revision 2 introduced. Each finding was
checked against the source and the mathematics, and all were accepted. Several of the revision-2
fixes above (5, 10, 14) are superseded as described here.

1. **Dropping CSF from the M0 fit biases calibration** (Codex's 25-voxel example: GM CBF about 27%
   low with zero residual). **Resolved by removing M0 regression entirely.** LR regresses
   y = q·ΔM, and calibration uses the centre voxel's m̄ (§5.1).
2. **The M0 solver failed for single-tissue neighbourhoods, and the batched solve and per-tissue
   status were inconsistent.** Removed along with M0 regression. Solves are now closed form under
   the status mask, and status is a single map again. Genuine zero ΔM is no longer excluded:
   validity comes from `quant_valid` (§5.2).
3. **CBF/s is not invertible where `ComputeCBF` sanitized results** (verified: `cbf.py:735-744`
   division and `nan_to_num`; `:512-522` short-TR correction and `zero_low_m0`). `ComputeCBF` now
   exports `mean_qdeltam`, `m0_calib` (the exact denominator) and `quant_valid` (§5.1).
4. **SCRUB/s is not a PV mixture** (verified: the prior term `mu * Globalprior` in CBF space,
   `utils/cbf.py:371, 407`). SCRUB LR is deferred. SCORE is kept, with a test that it is a
   retained-volume mean (§5.1, §7.1.6).
5. **Through-plane kernels broke under slice timing** (verified: slice-shifted PLDs at
   `cbf.py:524` enter q at `:739`). Regressing y = q·ΔM cancels q exactly, so any kernel extent is
   valid. Test 1 now uses slice timing and a 5x5x3 kernel. *Wrong; see revision-3 finding 1:
   tissue kinetics do not cancel, so slice-timed data are limited to one-slice kernels.*
6. **Constant `M0Estimate` does not give constant tissue coefficients.** Moot, because there is no
   M0 fit. m̄ is the scalar itself.
7. **The SE ignored fitted-M0 uncertainty.** Moot for the same reason. SE is now explicitly
   conditional on m̄_c, as the uncorrected CBF is, with no fitted denominator (§5.4).
8. **The 7.2(a) test could not run as written** (verified: no 1 mm ΔM file; `deltam_gt` is 4D and
   context-aligned). It now selects label rows by `aslcontext` and computes q from metadata. The
   oracle target comes from a whole-volume fit. The PV workflow is scored separately in 7.2(b).
9. **Zero-filled invalid estimates contaminate resampled valid ones.** Weighted projection and the
   new LR volumetric outputs now use normalized, validity-aware linear resampling, with a test
   (§6.2, §7.1.9).
10. **The reference grid does not span the T1w FOV** (verified: `GenerateSamplingReference` with
    `fov_mask=t1w_mask` crops to the bounding box plus 2 voxels). The text was corrected and a
    vertex-in-grid check added (§6.4).
11. **The discrete kernel FWHM was not guaranteed.** FWHM is now defined on the continuous kernel.
    The discrete kernel is its exact bin integral, with truncation and halo specified, and tests
    bound the L1 error and check slab equivalence on all slice axes (§4.2-4.3, §7.1.7-8).
12. **Conditional truth is not the LR estimand.** 7.2 now separates oracle, PV-map, preprocessing,
    local-constancy and calibration errors. Auto-versus-nominal is read net of them.
13. **Some surface diagnostics were missing or mislabelled.** A mean GM-fraction map was added.
    Zero support is defined. The effective-fraction interpretation is qualified. Total weight is
    explicitly deferred, with the reason. Distance weighting is now described as a scope decision
    (§6.4, §6.6).

### Revision 3 review (11 findings: 8 major, 3 minor)

The Codex adversarial review of revision 3 ran read-only numerical probes against the pinned F1
fixture and aslscan's kinetic model. Each finding was checked, and all were accepted.

1. **q does not cancel tissue kinetics across slices.** Evidence: `ComputeCBF` uses the blood-T1
   single-compartment q, while aslscan uses tissue T1′ and ATT per slice acquisition time. On F1,
   single-slice fits are exact (residual about 3e-8), but one whole-volume fit leaves 3.4%, and
   d_GM differs by about 12% between slices 5 and 30. With slice timing, LR kernels are now
   limited to one slice, and the oracle fits each slice separately. Test 1 is relabelled as an
   exact-model test (§5.1, §7.1.1, §7.2(a)).
2. **Parcellating LR maps counted failures as zeros** (verified: a shared `asl_mask` at
   `cbf.py:820-835`). Each LR map now gets its own validity mask, and coverage reports estimation
   coverage (§5.4).
3. **The reliability field was not normalized**, so lost support raised precision. Variance now
   uses the same normalized resampling, with a joint test (§6.2, §7.1.9).
4. **Uncorrected-CBF validity ignored `quant_valid`.** Validity rules are now given per CBF type
   (§6.2).
5. **Precision-only mode lacked p̃_GM for `projpvgm`.** `init_asl_pv_wf` now runs for any weighting
   mode, with a graph test (§4, §7.1.9).
6. **The vertex endpoint check was insufficient.** Workbench polyhedra use the 1-ring on both
   surfaces. The check is now exact via convexity, and its flag is carried into the mask (§6.4,
   §7.1.11).
7. **Resampled SE is not an SE.** LR SE maps are now native-only and are not projected. The surface
   reliability field is internal (§5.4, §6.2, §6.5).
8. **Status precedence could fit an invalid centre.** The guards (outside mask, invalid centre)
   now come first, with a test (§5.3, §7.1.3).
9. **Calibration claims were too broad.** Algebraic, noisy and physiological agreement are now
   separated. The oxford_asl statement is limited to ASLPrep's verified configuration, and the
   LR-vs-BASIL difference is descriptive (§5.1, §7.1.1, §7.2(b)).
10. **The RSS cutoff erased real residuals.** RSS is now computed directly over kernel offsets, and
    the tests cover positive near-exact residuals (§5.2, §7.1.5).
11. **Stale text.** The Lanczos boilerplate (`base.py:332`) becomes conditional. The temporal
    coverage text now notes SCORE's global retained-volume index (§5.4, §6.6).

## 12. Acceptance criteria

- **Defaults.** With no new flags, the workflow graph and the derivative list are unchanged (unit
  test 12), and existing integration tests pass.
- **`--pvc-method local-regression`** on single-PLD data:
  - writes the §4.4 and §5.4 derivatives through the sink config, with sidecars recording the
    method, the kernel, the thresholds, the FWHM and its source, and the calibration;
  - passes unit tests 1-8 and 12-14;
  - reports §7.2(a) and (b).
  - On multi-PLD or CBF-only data it fails with a clear message.
- **`--surface-weighting`**:
  - runs the Workbench nodes with `-weighted` and a single `volume_roi` source;
  - writes the support GIFTIs;
  - passes unit tests 9-11.
