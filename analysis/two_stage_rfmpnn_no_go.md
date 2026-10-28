# Two-stage RFdiffusion → ProteinMPNN path: no-go for the TfR1 apical epitope

## Decision

The first two design families generated with a **two-stage backbone-first workflow**
(RFdiffusion followed by ProteinMPNN) are closed as a **no-go for this epitope**.

This is a method-level decision, not a claim that RFdiffusion or ProteinMPNN are
generally unsuitable for TfR1.

## Evidence

### Family 1 — 87-aa seed 20261029

RFdiffusion placed the backbone near the intended apical surface, including
human L209 and G347. After sequence design, however:

- ProteinMPNN repeatedly converged to strongly E/K-enriched helical sequences.
- Sequence-only Boltz2 predictions moved to unrelated TfR1 surfaces.
- RF-derived structural templates did not rescue the intended apical pose.
- A29_03 templated Boltz2 protein ipTM:
  - human: ~0.144
  - mouse: ~0.155
- Expected apical anchors were absent from the final 5 Å interfaces.

Position-specific ProteinMPNN bias improved local chemistry but still failed to
retain the designed epitope after independent complex prediction.

### Family 2 — compact 57-aa seed 20261102

This backbone passed a stronger backbone-level gate:

- length: 57 aa
- target contacts: 344 / 347 / 348 / 362
- conserved preferred anchors: 344 / 347
- radius of gyration: 14.88 Å
- maximum CA span: 51.02 Å

ProteinMPNN produced more localized interface chemistry at positions 1/2/5,
but independent human/mouse Boltz2 validation again moved away from the intended
344/347 footprint.

Representative results:

- C57_02 mouse protein ipTM ~0.425, but the predicted interface moved to a
  different TfR1 surface rather than preserving the human epitope.
- C57_05 mouse protein ipTM ~0.433 with the same failure mode.
- C57_12 human/mouse protein ipTM ~0.170 / ~0.150 and no retained 344/347 core.

A higher mouse ipTM therefore did **not** establish cross-reactivity: the human
and mouse predictions were not using the same epitope.

## Interpretation

Across two independent backbone families, the repeated failure occurs at the
**handoff between backbone generation and sequence design**. A backbone can be
geometrically placed against the desired epitope, yet the separately designed
sequence does not necessarily encode enough interaction energy to preserve that
pose during independent refolding/docking.

Further tuning of the same two-stage family would confound method failure with
candidate selection and consume additional validation budget.

## Pivot

The next campaign uses **RFD4-Proteina**, which jointly generates:

- binder backbone,
- amino-acid sequence,
- side-chain/all-atom structure,

while the TfR1 target is conditioned directly and the conserved apical epitope
is supplied as hotspot conditioning.

Initial RFD4 target:

- human TfR1 structural crop: residues 190–400 from 11QC chain A
- conserved hotspot subcluster: E343 / K344 / G347 / E350
- binder length: 80–100 aa
- design task: generation
- first request: one replicate only

The first RFD4 result is still an **unvalidated design candidate**. It must pass
independent human/mouse structural validation before any pH-switch histidine
scan is attempted.
