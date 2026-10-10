# pcasl1pld_grase

Written for ASLPrep (integration tier). F3. Siemens-style 3D GRASE PCASL with a separate M0 at TR 3 s (short-TR correction), four suppression pulses, LabelingEfficiency kept in the sidecar, M0 stored divided by 10.

Derived from bids-examples commit 7150efcf2c040465e9fcb6745f0b5aa4e84919c1, as copied into aslscan tests/fixtures/protocols/asl005/asl.json.

Changes:

- kept the acquisition fields aslscan reads; dropped descriptive fields
- segmentation and phase-encoding direction come from aslscan's asl005_p5 overlay
- the separate M0 is simulated at TR 3 s (asl005_p5: 6 s) to exercise ASLPrep's short-TR correction
- the M0 is divided by 10 after simulation (m0_divisor)
