# st_spiral

Stress case for the D3 perfusion image (plan: "Stress cases"): stack-of-spirals readout
artifacts. aslscan's `p8_spiral_kz` protocol (Vidorreta et al. 2017; Siemens SPIRAL_V20_HCP:
unbalanced pCASL, 1.8 s label and delay, 3.75 mm, 6/8 partial Fourier and acceleration 2 along kz,
an m0scan and 14 pairs, background suppression with presaturation) on the `tfmni-b0` phantom,
reconstructed by gridding with Nyquist-spaced turns, which leaves a background haze like that of
the real data (ds000235). Four receive coils in four rings along z, the minimum for the kz
acceleration (the scan used a 64-channel head coil). The real family's main problem, the
single-volume ringing that defeats the conventional reference (#503, #186), is only partly
reproduced: aslscan's haze does not follow the head as the real rings do.
