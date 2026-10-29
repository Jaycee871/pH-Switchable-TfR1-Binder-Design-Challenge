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


## v4 hard-hotspot update

The v3 investigation showed that hotspot syntax was valid, but the default ColabDesign
inter-chain contact objective used a 20 A cutoff. v4 therefore tests a stricter objective:

- `inter_contact_distance = 8.0 A`
- `weights_con_inter = 2.0`
- `inter_contact_number = 2`
- fixed 60-aa binder
- one attempt only

The run also independently recomputes the final 5 A target contact set. A FreeBindCraft
trajectory is not considered hotspot-engaging merely because its internal loss accepts it.
