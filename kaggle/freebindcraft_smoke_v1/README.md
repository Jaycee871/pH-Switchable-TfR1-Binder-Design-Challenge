# TfR1 FreeBindCraft Kaggle smoke

Purpose: determine whether FreeBindCraft can run on the available Kaggle T4 for the trimmed TfR1 apical target.

- Target: PDB 11QC chain A residues 189-383
- Hotspots: 343,344,347,350
- Binder length: 55-75 aa
- Relaxed filters
- 4 hallucination trajectories maximum
- 1 accepted design target
- PyRosetta disabled
- ipSAE ranking
- Smoke-stage iteration counts reduced

A successful environment run is useful even if zero designs are accepted. Production parameters must not inherit the reduced smoke-test iteration budget.
