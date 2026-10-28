# RFdiffusion smoke test 001 — QC decision

**Run:** GitHub Actions 40679577930  
**Artifact:** tfr1-rfdiffusion-smoke (artifact 12812486681)  
**Result:** infrastructure PASS; scientific gate REJECT

## Generated complex

- Binder chain: A
- Binder length: 82 aa
- Target chain: B
- Target length: 639 modeled residues
- RFdiffusion elapsed time: 120.1 s
- Requested hotspots: human 206, 209, 214, 344, 374

## Direct target contacts

Using a 5 Å heavy-atom cutoff, the generated binder contacted human-equivalent target residues:

`196, 209, 210, 211, 212`

This includes the conserved anchor **L209**, but also places the interface directly over **V210, Y211, and L212**, which are cross-species-risk positions in the Proteinbase structure-aware mapping.

## Decision

Do **not** spend ProteinMPNN/Boltz validation budget on this backbone.

The run proves that the hosted NVIDIA RFdiffusion path and repository secret are functioning. The next backbone should shift the hotspot constraint toward the conserved 343–374 apical subcluster to reduce dependence on the human-specific 210–212 region.
