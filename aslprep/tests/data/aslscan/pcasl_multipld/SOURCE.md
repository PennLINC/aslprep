# pcasl_multipld

Written for ASLPrep (integration tier). F5. 2D EPI multi-delay PCASL (LD 1.4 s; PLDs 0.25-1.5 s), label first, unequal repeats per delay (8, 8, 8, 6, 6, 6 pairs), an arterial compartment, LabelingEfficiency 0.88 from the sidecar, M0 stored divided by 10.

Derived from bids-examples commit 7150efcf2c040465e9fcb6745f0b5aa4e84919c1, as copied into aslscan tests/fixtures/protocols/asl004/asl.json.

Changes:

- slices 30 ms apart (asl004: 45.2 ms) and TR 4.5 s (asl004: 4.05 s), so the 48 slices the phantom needs fit after the longest delay
- TotalReadoutTime 25 ms (asl004: 60 ms, which cannot precede its 14 ms echo time on a 68-line readout)
- unequal repeats per delay (asl004: 8 pairs at every delay)
- an arterial compartment (aBV 0.02/0.01, aBAT 0.7/1.0 s in GM/WM)
- the separate M0 is simulated at TR 8 s and divided by 10 after simulation
