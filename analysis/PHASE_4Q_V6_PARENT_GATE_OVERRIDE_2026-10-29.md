# Phase 4Q — v6 tempered-binary parent, gate-override validation

Date: 2026-10-29 (run 23:28–23:35 Taipei).

## Why this run exists

`TFR1_FBC343_V6_l60_s289817` failed the strict local 5 A gate in Phase 4M **only** on
species-risk contacts (210, 211, 348, 371) and was therefore never sent to full
human/mouse Boltz-2. Phases 4N and 4O then spent Boltz-2 budget on *derivatives* of this
backbone while the parent itself — the only trajectory with zero severe clashes AND
genuine conserved hotspot-core engagement — stayed unmeasured.

This phase bypasses **only** the upstream local pre-filter. The downstream Boltz-2 gate is
byte-identical to Phase 4F/4G/4M. Recorded as a gate-override diagnostic in
`designs/family_D_freebindcraft/v6_gate_override_pose_screen/selection_decision.json`.

## Independent full-complex result

### Human TfR1

- protein iPTM: **0.690**
- pTM: 0.768
- confidence: 0.789
- complex pLDDT: 0.813
- preferred contacts: **231, 347, 350**
- `preferred_core_retained`: **True**
- species-risk contacts: **371 only**

### Mouse TfR1

- protein iPTM: **0.282**
- pTM: 0.645
- confidence: 0.556
- complex pLDDT: 0.624
- preferred contacts: **none**
- `preferred_core_retained`: **False**

### Cross-species

- homologous overlap (human numbering): **[]**
- overlap count: **0**
- Jaccard: **0.0**
- `cross_species_pose_pass`: **False**

## Findings

**1. This is the best human-side result the project has produced.**

Compared with the previous references:

| Candidate | Human iPTM | Mouse iPTM | core retained |
|---|---:|---:|---|
| `l70_s574476` (FreeBindCraft v3) | 0.432 | 0.677 | **False** |
| `MS_15` (4O dual-state) | 0.734 | 0.251 | — |
| **`V6_l60_s289817` (this run)** | **0.690** | **0.282** | **True** |

MS_15 scores marginally higher on human, but the v6 parent is the only candidate that
combines a high human iPTM with `preferred_core_retained = True` and essentially no
species-risk engagement.

**2. The local gate was partly over-pessimistic.** The trajectory-level gate rejected this
design for contacting 210, 211, 348 and 371. In the independent full-complex human
prediction, **only 371 survives**; 210, 211 and 348 are not contacted at all. The strict
local gate therefore discarded the project's best parent on the strength of three contacts
that do not reproduce under independent full-ectodomain prediction. Any remaining
FreeBindCraft budget should re-examine whether that pre-filter is costing real candidates.

**3. The mouse failure is worse than a low score — it is the wrong surface entirely.**
The mouse contacts are 401, 402, 442–450, 603–616, 668–675, 760, 761. The homologous
apical domain is ~191–385. The binder is not making a weak apical pose; it is docking to a
completely different region of the mouse ectodomain. This is consistent with Phase 4O's
conclusion that the bottleneck is backbone/docking geometry, and it is now demonstrated on
the parent backbone itself rather than only on its sequence derivatives.

## Status

Pose passes across the whole project remain **zero**. The v6 fixed-backbone family is
closed (Phase 4O). This run does not reopen it; it quantifies what was lost.
