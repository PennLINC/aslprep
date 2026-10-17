# st_noperf

Stress case for the D3 perfusion image (plan: "Stress cases"): a region of absent perfusion.
`rc_pcasl` (2D PCASL, separate posed M0, fieldmap) on the `tfmni-b0-noperf` phantom, whose
perfusion is zero inside a 30 mm sphere over the left middle cerebral artery territory
(MNI -40, -20, 20), as in an occlusion. Tissue contrast is unchanged, so the conventional
reference does not see the territory and the perfusion image does.
