# A29_03 — sequence-only Boltz2 cross-species QC

## Result

Candidate **TFR1_A29_MPNN_03** was evaluated against full human and mouse TfR1 ectodomain sequences using unconstrained / sequence-only Boltz2 complex prediction.

- Human overall confidence: **0.7512**
- Mouse overall confidence: **0.4431**

At first glance this looks like a substantial species gap. However, direct interface inspection shows that the predicted complexes do **not** preserve the RFdiffusion-designed apical pose.

### Human sequence-only pose

The 5 Å interface is dominated by human TfR1 residues in the approximately **612–657** region, with additional contacts elsewhere. This is far from the intended apical design neighborhood around the conserved 209 / 347 region.

### Mouse sequence-only pose

The predicted mouse complex similarly shifts to a different surface, with many contacts in the approximately **651–685** region.

## Interpretation

The human-vs-mouse confidence difference cannot yet be treated as a valid cross-reactivity estimate because the two sequence-only predictions are evaluating alternative docking poses rather than the intended RFdiffusion interface.

Therefore A29_03 is **not rejected on cross-species grounds from this test alone**.

## Next test

Use the RFdiffusion seed-20261029 complex geometry as an explicit Boltz2 structural template, with chain-specific template assignment, and then compare whether the intended apical interface can be retained for:

1. human TfR1 + A29_03
2. mouse TfR1 + A29_03

This separates **interface-geometry retention** from unconstrained global docking.
