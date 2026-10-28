# A29_03 — final cross-species structural QC

## Decision

**REJECT for downstream pH-switch optimization.**

A29_03 folds confidently as an isolated sequence in the predicted complexes, but neither sequence-only nor RFdiffusion-templated Boltz2 supports the intended TfR1 apical binding interface.

## Templated Boltz2 metrics

| Metric | Human TfR1 | Mouse TfR1 |
|---|---:|---:|
| Overall confidence | 0.7649 | 0.6829 |
| pTM | 0.8157 | 0.7247 |
| protein ipTM | 0.1437 | 0.1550 |
| complex pLDDT | 0.9202 | 0.8148 |
| interface-weighted pLDDT | 0.8827 | 0.7852 |
| complex ipDE (Å) | 18.06 | 15.42 |

The high per-complex pLDDT values indicate that the individual chains can be modeled confidently. The very low protein ipTM and large ipDE indicate that their **relative interface geometry is not trustworthy**.

## Contact-location check

The RFdiffusion design was intended to use the apical neighborhood around conserved positions including human L209 and G347.

After templated Boltz2 refinement:

- human contacts shifted primarily to approximately residues **526–654**;
- mouse contacts shifted primarily to approximately residues **523–657**;
- the intended 209 / 347 apical anchors were not retained in the final 5 Å interface.

Therefore the apparent human/mouse complex predictions do not validate the designed pose.

## Methodological implication

The first ProteinMPNN batch used `use_soluble_model=True`. All 16 designed sequences were highly enriched in charged E/K residues and contained no histidines. This family is useful as a negative/baseline observation, but its interface chemistry appears too weak to preserve the RFdiffusion docking pose.

## Next action

Keep the same RFdiffusion seed-20261029 backbone and redesign binder chain A with the **standard ProteinMPNN model** (`use_soluble_model=False`) while keeping TfR1 fixed.

This isolates the sequence-design variable without spending additional RFdiffusion budget.
