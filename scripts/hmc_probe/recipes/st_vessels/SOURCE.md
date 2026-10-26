# st_vessels

Stress case for the D3 perfusion image (plan: "Stress cases"): focal macrovascular signal.
`rc_pcasl` (2D PCASL, 1.8 s label, 1.8 s delay, separate posed M0, fieldmap) on the
`tfmni-b0-vessels` phantom, which adds arterial blood volume along the M1, proximal ACA,
pericallosal and P2 segments (peak 0.3, Gaussian across a 1.5-2.5 mm radius) with a 3.0 s
arterial transit time: a very slow-flow case, chosen so that label is still in the arteries when
every slice is read (3.6-4.74 s after labeling starts).
