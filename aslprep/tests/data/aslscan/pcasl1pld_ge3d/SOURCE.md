# pcasl1pld_ge3d

Written for ASLPrep (integration tier). F2. GE-style 3D spiral PCASL: one delta-M volume and an included M0, four background-suppression pulses, no LabelingEfficiency in the sidecar, M0 stored divided by 96.

Derived from bids-examples commit 7150efcf2c040465e9fcb6745f0b5aa4e84919c1, as copied into aslscan tests/fixtures/protocols/asl001/asl.json.

Changes:

- kept the acquisition fields aslscan reads; dropped descriptive fields (model, coil, software and labeling-location descriptions)
- the spiral interleaves, readout time and dwell time come from aslscan's asl001_p5 overlay (they are not in the sidecar)
- the M0 volume is divided by 96 after simulation (m0_divisor), as GE stores it scaled
