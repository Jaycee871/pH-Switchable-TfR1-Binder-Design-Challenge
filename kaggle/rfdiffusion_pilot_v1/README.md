# TfR1 RFdiffusion pilot v1

This Kaggle GPU kernel generates a small first-pass library of de novo protein-binder backbones against the **human TfR1 apical domain (PDB 11QC, chain A residues 189–383)**.

The pilot intentionally uses a cropped generation target to keep RFdiffusion tractable on a Kaggle T4. Candidate sequences are **not** considered validated at this stage; they must subsequently pass ProteinMPNN sequence design, full-ectodomain human/mouse structural screening, pH-switch analysis, and novelty/developability filters.

Families:

- `cross_species_loop_5` — conserved/adjacent loop anchors
- `cross_species_edge_5` — broader conserved apical edge
- `functional_patch_mixed_5` — experimentally validated 11RN footprint, including species-variable positions
