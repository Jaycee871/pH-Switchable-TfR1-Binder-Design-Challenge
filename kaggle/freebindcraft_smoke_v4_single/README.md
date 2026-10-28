# FreeBindCraft smoke v4 — single fixed-shape attempt

Prepared as a bounded fallback if the long v3 T4 run proves too slow.

Changes relative to v3:
- fixed binder length: **60 aa**
- official CLI hard cap: **--max-attempts 1**
- one trajectory maximum
- deterministic NumPy sampling seed: **20261029**
- same TfR1 apical target and conserved hotspots 343/344/347/350
- same Kaggle-specific NVIDIA driver-library restoration

This directory is intentionally **not connected to an automatic workflow yet**.
It should be launched only after the current v3 kernel finishes, so it does not
consume overlapping Kaggle GPU quota or overwrite the active kernel session.
