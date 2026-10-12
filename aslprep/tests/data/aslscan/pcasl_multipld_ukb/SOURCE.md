# pcasl_multipld_ukb

Written for ASLPrep (integration tier). F8. A UK Biobank-like acquisition: multi-delay PCASL with a 2-shot 3D GRASE readout, coarse voxels, heavy through-slice blurring, and one label/control pair per delay. Written to test whether the tissue-type pattern reported in https://neurostars.org/t/aslprep-csf-shows-elevated-cbf-in-basil-output-multi-pld-uk-biobank-data/36453 follows from the sequence.

Not derived from a BIDS example. The protocol values come from:

- UK Biobank Brain Imaging Documentation, sections 2.10 and 3.13 (https://biobank.ctsu.ox.ac.uk/crystal/crystal/docs/brain_mri.pdf);
- Okell et al., "Brain Perfusion Imaging of a Large Population: Arterial Spin Labelling MRI in UK Biobank", eLife reviewed preprint 108799 (https://elifesciences.org/reviewed-preprints/108799), Methods.

As in UK Biobank:

- PCASL, labeling duration 1.8 s, post-labeling delays 0.4, 0.8, 1.2, 1.6 and 2.0 s, one label/control pair per delay
- 3D GRASE, 2 shots, echo time 29 ms, refocusing flip angle 120 degrees
- voxels 3.4 x 3.4 x 4.5 mm on a 64 x 48 x 32 grid (field of view 218 x 163 x 144 mm; the phantom is cropped to it, which cuts 0.3 % of the brain at the bottom of the cerebellum)
- pre-saturation and two background-suppression inversions
- a separate M0 without labeling or suppression, repetition time 5.39 s
- the ASL series stored 10 times larger than the M0 (`m0_divisor`, run with `--m0_scale=10`)

Differences, from limits of aslscan at the pinned revision:

- UK Biobank uses 6/8 partial Fourier along the slice direction (12 echoes per shot). aslscan has no partial Fourier along kz, so each shot reads 16 echoes at the same 29 ms spacing. The four extra echoes carry under 1 % of the grey-matter signal, so the through-slice blurring is close.
- UK Biobank encodes phase left-right. aslscan encodes phase along the second image axis only, so phase runs anterior-posterior here, on the same matrix and field of view.
- UK Biobank times the two inversions to null T1 = 0.7 s and 1.4 s at 100 ms before excitation, which places pulses inside the labeling period. aslscan cannot keep the label's sign through such pulses, so both pulses follow labeling here, timed to minimise the static grey-matter, white-matter and CSF signal while keeping it positive.
- UK Biobank minimises the repetition time for each delay. Here it is 4.6 s for every volume; with pre-saturation the static signal does not depend on it.

Other choices: an arterial compartment (aBV 2 % in grey matter, 1 % in white matter), a per-pair grey-matter SNR of 8, and the anatomy given as derivatives.
