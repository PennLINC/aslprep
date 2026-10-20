# headmotion

Written for ASLPrep (integration tier). F7. The F1 protocol at multiband 3 (13 slice groups 80 ms apart) with head motion between and within volumes, delta-M SNR 5 per pair, the T1w header offset by a known rigid transform (3, -4, 2 mm; 5, -3, 4 degrees), and no anatomical derivatives, so ASLPrep preprocesses the anatomy itself.

Motion:

- `motion.tsv` is a timed trace (a `time` column, seconds from the first labeling), sampled at
  every slice group's own time (`within_volume_sampling = "shot"`). Slice group `s` of volume `v`
  is read at `5 v + 3.6 + 0.08 s` s. The trace is piecewise linear, so it holds only its knots:
  - a drift from rest to (1.2, -0.8, 0.4 mm; 0.017, -0.009, 0.012 rad) over the 100 s run;
  - a step of (0.8, 0, 0.6 mm; 0.012, 0, 0 rad) from 30.5 s to 31.5 s, between volumes 5 and 6;
  - a nod of (0, 1.5, 2.0 mm; 0.035, 0, 0 rad) that rises from 59.0 s to 59.3 s, during volume
    11's readout (slice groups 6 to 8 catch it rising), holds through volume 12, and falls between
    volumes 12 and 13;
  - a twitch of (1.5, 0, 0 mm; 0, 0, 0.03 rad) from 83.9 s to 84.3 s, inside volume 16's readout.
- `events.tsv` adds two slice-group events (the `desc-motionEvents_gt.tsv` format): volume 3,
  group 7 jumps (0, 0, 0.8 mm; 0.014, 0, 0 rad) with its slices at 70 % signal, and volume 14,
  group 4 jumps (-0.6, 0, 0 mm; 0, 0, -0.01 rad) at 80 %. A jump holds for the rest of its volume.
- `spin_history = true`: through-plane motion between slice excitations changes the tissue
  magnetization (aslscan's box-slab model).

`motion.tsv` was written by a script that sums the components above at the union of their knots.
