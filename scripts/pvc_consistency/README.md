# Truth-free consistency tests of partial volume correction

Real data have no ground-truth CBF. These tests judge partial volume correction (PVC) by
consistency instead, on ASLPrep's own outputs:

- **Test-retest:** the same protocol in several sessions should give the same regional
  pattern.
- **Resolution:** the same subject scanned at two resolutions should give the same corrected
  gray matter CBF. Uncorrected CBF falls at the coarser resolution, because each voxel mixes
  in more white matter.
- **Contrast:** the GM/WM CBF ratio of each run cancels the run's calibration scale, so it
  can be compared across protocols.

`consistency.py` computes all three from an ASLPrep derivatives directory. It compares
ASLPrep's uncorrected CBF, BASIL's uncorrected CBF, ASLPrep's local-regression GM CBF
(`desc-lrGM`) and BASIL's partial-volume-corrected GM CBF (`desc-basilGM`).

## Running

The runs below used `pennlinc/aslprep:unstable` with this branch's package mounted over the
image's, and these options:

```
--basil --pvc-method local-regression --output-spaces asl T1w --fs-no-reconall
```

`--asl-effective-resolution` was left at its default, `nominal`.

| Dataset (aslprep-zoo) | Subject | Runs | BIDS filter |
|---|---|---|---|
| `kashyap` | `01` | sessions 20, 32 and 64, each with `acq-hires` (2 mm) and `acq-lores` (3.4 mm) 3D PCASL | `{"asl": {"run": 1}}` (PCASL only) |
| `ds007116` | `20721` | sessions 1-3, 3D stack-of-spirals PCASL; session 1 is the `acq-VARIANT...` acquisition | none |

Then, with any environment that has numpy, pandas and nibabel:

```
python consistency.py retest   <ds007116 derivatives> --subject 20721 --exclude VARIANT
python consistency.py retest   <kashyap derivatives> --subject 01 --group acq-lores
python consistency.py resolution <kashyap derivatives> --subject 01
python consistency.py contrast <kashyap derivatives> --subject 01
```

`--atlas` selects the parcellation (default `4S156Parcels`). Every method is scored on the
same parcels: those with a positive value, and at least 50% coverage, for every method in
every run compared.

## Results (Oct 2026)

### Test-retest

- **ICC** is ICC(3,1) of log parcel CBF, with parcels as targets and sessions as raters.
- **CV** is the median within-parcel coefficient of variation across sessions.
- **Scaled CV** is the same after dividing each session by its median.

**ds007116, sessions 2 and 3** (same protocol; 4S156Parcels, 109 parcels):

| Map | ICC | CV | Scaled CV |
|---|---|---|---|
| ASLPrep uncorrected | 0.86 | 4.4% | 3.1% |
| BASIL uncorrected | 0.85 | 5.4% | 2.7% |
| ASLPrep LR GM | **0.92** | **3.3%** | 3.3% |
| BASIL PVC GM | 0.86 | 4.6% | 4.0% |

At 4S456Parcels the order is the same: LR GM has ICC 0.90, against 0.84 uncorrected, 0.80
BASIL and 0.79 BASIL PVC.

With session 1 included, every map's raw CV rises to 21-25%, because that acquisition reads
about 35% lower overall. Its scaled CV is 7-9%, with LR GM lowest at 7.1%.

**kashyap, three sessions per resolution** (4S156Parcels):

| Map | ICC, 3.4 mm | Scaled CV, 3.4 mm | ICC, 2 mm | Scaled CV, 2 mm |
|---|---|---|---|---|
| ASLPrep uncorrected | 0.88 | 7.7% | 0.69 | 8.6% |
| BASIL uncorrected | 0.90 | 6.8% | 0.74 | 8.5% |
| ASLPrep LR GM | **0.96** | **5.3%** | **0.80** | 9.5% |
| BASIL PVC GM | 0.89 | 8.6% | 0.76 | 10.0% |

LR GM gives the most reproducible regional pattern in every comparison. Part of that is
smoothing: LR estimates one GM value per 5x5x1 neighbourhood, and smoother maps are more
reproducible whether or not they are more accurate (#739).

### Resolution (kashyap)

| GM/WM contrast | ses-20 2 mm | ses-20 3.4 mm | ses-32 2 mm | ses-32 3.4 mm | ses-64 2 mm | ses-64 3.4 mm |
|---|---|---|---|---|---|---|
| ASLPrep uncorrected | 1.29 | 2.77 | 1.29 | 2.59 | 1.74 | 2.71 |
| BASIL uncorrected | 1.23 | 2.66 | 1.26 | 2.54 | 1.60 | 2.60 |
| ASLPrep LR | 1.53 | 5.12 | 1.56 | 5.11 | 2.26 | 5.28 |
| BASIL PVC | 1.37 | 2.56 | 1.39 | 2.54 | 1.86 | 2.60 |

- **The 2 mm runs have lower GM/WM contrast than the 3.4 mm runs, even uncorrected.** Partial
  volume mixing alone predicts the opposite. So the 2 mm acquisition is effectively much
  blurrier (or noisier) than its voxel size, and nominal resolution does not describe it.
  Neither correction can make the resolutions agree while the PV maps assume 2 mm.
- **Parcel CBF lores/hires ratios:**
  - uncorrected: 0.97-1.17;
  - BASIL uncorrected: 0.99-1.18;
  - LR GM: 1.15-1.44;
  - BASIL PVC GM: 1.12-1.44.

  The two protocols also differ in labeling duration (1.6 against 1.8 s) and readout, so
  these absolute ratios mix calibration with partial volume effects. Both corrections raise
  the 3.4 mm runs relative to the 2 mm runs by about the same factor (1.15-1.27).
- **LR's contrast at 3.4 mm (about 5.1) is about twice BASIL's (2.5-2.6) and about twice
  ds007116's.** LR's WM estimate is low there. Without truth, this test can't say which is
  right.

### Contrast across sessions (ds007116)

The contrast is stable across sessions for every map: LR 1.94-2.09, BASIL PVC 1.36-1.44,
ASLPrep uncorrected 1.81-1.94. BASIL's corrected contrast is lower than its uncorrected one,
because its corrected WM is high.

## Interpretation and next steps

- On test-retest, LR is at least as reproducible as BASIL at the parcel level.
- The resolution test is inconclusive, because the 2 mm acquisition's effective resolution
  is much coarser than its voxel size. The natural next step is the effective-resolution
  sweep on kashyap's 2 mm runs: if a wider assumed FWHM brings their corrected contrast to
  the 3.4 mm value, that is a real-data estimate of the protocol's effective resolution
  (#745).
- These are single subjects. More subjects would be needed for between-subject ICCs.
