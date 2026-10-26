# pcasl_multipld_grase

Written for ASLPrep (integration tier). F3. Siemens-style 3D GRASE multi-delay PCASL (LD 1.8 s; PLDs 0.5-4.0 s in 0.5 s steps, cycled twice), with a separate M0 at TR 3 s (short-TR correction), four suppression pulses, LabelingEfficiency kept in the sidecar, M0 stored divided by 10.

A 3D readout has no slice-timing offsets, so every voxel's earliest delay is 0.5 s, and the long delays sample the whole bolus. This is the protocol in which ASLPrep should locate ATT (issue #742).

Derived from bids-examples commit 7150efcf2c040465e9fcb6745f0b5aa4e84919c1, as copied into aslscan tests/fixtures/protocols/asl005/asl.json.

Changes:

- kept the acquisition fields aslscan reads; dropped descriptive fields
- segmentation and phase-encoding direction come from aslscan's asl005_p5 overlay
- eight delays (asl005: one, 2.0 s), cycled so that the volumes of a delay are not contiguous, with two pairs per delay
- TR 6.95 s (asl005: 4.95 s), the asl005 margin after the longest delay
- the four suppression pulses at 1.85-2.15 s (asl005: 2.29-3.705 s), before the earliest readout
- the separate M0 is simulated at TR 3 s (asl005_p5: 6 s) to exercise ASLPrep's short-TR correction
- the M0 is divided by 10 after simulation (m0_divisor)

Replaces pcasl1pld_grase, the single-delay version of this protocol.
