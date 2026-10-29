# Phase 4I — full human/mouse validation of hotspot-engaging FreeBindCraft trajectory

Date: 2026-10-29. Run: GitHub Actions `40816382080`.

Candidate: `TFR1_FBC343_l67_s385919` (67 aa), harvested from the cancelled
FreeBindCraft v3 smoke.

This trajectory was selected **not** because of FreeBindCraft's own score
(`i_pTM=0.19`), but because independent geometry analysis showed it was the only
one of eight harvested trajectories to make a real 5 Å contact with the requested
343-350 hotspot cluster.

Local FreeBindCraft-complex contacts after correcting target renumbering:

- preferred contacts: human **231, 350**
- hotspot-core contacts: **350**
- explicit species-risk contacts: none
- severe clashes <2 Å: 0

## Full-complex Boltz-2 result: FAIL

| Metric | Human | Mouse |
|---|---:|---:|
| protein ipTM | **0.393** | **0.196** |
| pTM | 0.657 | 0.602 |
| complex pLDDT | 0.804 | 0.620 |
| preferred contacts | **231, 350** | none |
| apical contact | true | **false** |
| per-species pose gate | fail | fail |

Human retained the intended local epitope remarkably well: 231 and 350 remained in
the 5 Å contact set after full-ectodomain prediction, with zero explicit risk contacts.
However, protein ipTM remained below the 0.50 gate.

Mouse collapsed completely away from the apical domain. Contacts moved to residues
~398, 444-448, 599-606 and 761-763. No preferred contact remained.

Cross-species homologous overlap: **0 residues**, Jaccard **0.0**.

## Interpretation

This is more informative than Phase 4H.

Phase 4H's high-FreeBindCraft-score trajectory (`l70_s574476`) had better nominal
confidence but never engaged the requested hotspot. Phase 4I proves that FreeBindCraft
*can* generate a trajectory that actually touches the desired 343-350 region, and that
the human full-complex model can preserve that local geometry. The remaining failure is
now specifically **cross-species pose stability / interface confidence**, not merely
hotspot targeting.

The v3 harvest is therefore exhausted:

- `l70_s574476`: high internal i_pTM, wrong local epitope, full-pose fail.
- `l67_s385919`: correct local hotspot contact, human epitope retained, mouse collapse.
- `l73_s257087`: low internal i_pTM and no hotspot-core contact.
- five LowConfidence trajectories: no hotspot-core contact.

No pH engineering should be performed on any v3 trajectory.

## Root-cause insight for v4

ColabDesign receives `target_hotspot_residues` correctly. Numeric hotspots such as
`343,344,347,350` are valid for a single target chain because `prep_pos()` defaults
to the first target chain when no chain prefix is supplied.

The weakness is the default hotspot-contact objective itself:

- inter-chain contact cutoff: **20 Å**
- contacts per hotspot position: **2**
- categorical (non-binary) contact loss

With a 20 Å cutoff, the optimisation can satisfy the hotspot-biased loss without
placing the final 5 Å interface on the hotspot. v4 should therefore tighten the contact
objective rather than merely change hotspot syntax.
