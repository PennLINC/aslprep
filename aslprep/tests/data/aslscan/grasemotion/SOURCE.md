# grasemotion

Written for ASLPrep (integration tier). F8. The F3 Siemens-style 3D GRASE PCASL protocol (four ky segments, each read after its own labeling, four suppression pulses) with 6 pairs, a separate M0 at TR 6 s, head motion between volumes and between the segments of a volume, delta-M SNR 10 per pair, and anatomical derivatives under a known rigid offset (-2, 3, 1.5 mm; -3, 2, 4 degrees).

Derived from the pcasl1pld_grase recipe, since replaced by its multi-delay version, pcasl_multipld_grase.

Changes:

- 6 pairs instead of 8; the M0 at TR 6 s and stored undivided
- `motion.tsv` is a timed trace (a `time` column, seconds from the first labeling), sampled at each
  segment's own time (`within_volume_sampling = "shot"`), so the segments of one volume can be
  read at different poses: k-space lines that disagree, which the reconstruction turns into
  ghosting and blur. Segment `s` of volume `v` is read at `3.8 + 4.95 (4 v + s)` s. The trace is
  piecewise linear, so it holds only its knots:
  - a drift from rest to (1.0, 0.6, -0.5 mm; 0.012, 0.008, -0.01 rad) over the 240 s run;
  - a step of (-0.8, 0.5, 0.7 mm; 0, 0.012, 0 rad) from 80 s to 81 s, between volumes 3 and 4;
  - a nod of (0, 1.2, 1.8 mm; 0.03, 0, 0 rad) that rises from 149 s to 150 s, between segments 1
    and 2 of volume 7, holds through segment 0 of volume 8, and falls before its segment 1.

`motion.tsv` was written by a script that sums the components above at the union of their knots.
