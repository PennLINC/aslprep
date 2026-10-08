# pasl_multipld

Written for ASLPrep (integration tier). F6. 2D EPI multi-delay PASL, Q2TIPS (0.7, 1.6 s), TIs 0.9-3.0 s, label first, M0 stored divided by 10.

Derived from bids-examples commit 7150efcf2c040465e9fcb6745f0b5aa4e84919c1, as copied into aslscan tests/fixtures/protocols/asl003/asl.json.

Changes:

- TIs from 0.9 s (asl003 starts at 0.3 s, before its own 0.7 s bolus cut-off, which aslscan refuses)
- 2D EPI at 4 x 4 x 6 mm instead of 3D GRASE at 8 x 4 x 6 mm: the GRASE segmentation does not fit the full phantom's phase-encoding lines with NumberShots 2, and 3D GRASE is covered by pcasl1pld_grase
- no background suppression (asl003: two pulses, which for PASL need the bolus-position model); suppression is covered by the PCASL recipes
- TR 4 s (asl003: 3.5 s) for the 2D slice train
