# TfR1 final-v1 freeze note

Date: 2026-10-31.

## Decision

The final-v1 portfolio contains **16 designs**, not 20.

Tracks 2/3 permit *at most* 20 designs; filling all 20 slots is not required.
The v0 portfolio contained four very close or lower-information variants that
would reduce within-submission diversity without adding a qualitatively new
mechanism or validation signal.

Removed from final-v1:

- **TF4_008_E47D** — near-duplicate of the v6/R42P sequence family.
- **TF4_011_F37** — near-duplicate of F36/F39; its cross-species overlap is distal.
- **TF4_018_R2_72** — redundant with the F39/R2 sparse-reconciliation cluster.
- **TF4_020_L67_XS07** — same l67 rescue family as XS05 with weaker distinguishing evidence.

## Why TF4_007_V5_T1 is retained despite one cysteine

TF4_007_V5_T1 has one unpaired cysteine, which is a real expression/developability
liability and is explicitly disclosed in metadata. It is retained because it is
sequence-diverse, originates from the F39 warm-start basin under mined hard-negative
design, and produces the largest V5 homologous overlap (7 residues; Jaccard 0.1489).
It contributes information not represented by a near-identical substitute.

## Final-v1 policy

No candidate will be added solely to reach 20 designs. A late replacement should
occur only if it provides qualitatively stronger evidence, such as:

- sequence-only full-target human+mouse binding;
- desired/apical engagement in both species;
- a credible pH-switch mechanism with preserved neutral-pH binding;
- or a substantial reduction in sequence liability without losing a unique design family.

The upload CSV remains minimal and uses only the required fields:
`name,sequence,molecule_class`.
