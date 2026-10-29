# Phase 4P — Odin-Multi random-initialization smoke and seeded-v2 escalation

Date: 2026-10-30 (Asia/Taipei)

## Purpose

After FreeBindCraft v6 established a strong human apical pose but failed mouse
cross-reactivity, and fixed-backbone ProteinMPNN rescues also failed, Odin-Multi
was introduced to optimize one shared binder sequence against human and mouse
TfR1 contexts inside the same AF2/ColabDesign trajectory.

The first Odin campaign was intentionally a one-trajectory smoke test. Its role
was to prove the multi-context execution path and determine whether a random
initialization could enter a joint human/mouse binding basin.

## Smoke configuration

- binder length: 60 aa
- trajectories: 1
- human context: apical TfR1 residues 189-383
- mouse context: homologous apical TfR1 residues 191-385
- human hotspots: 343 / 344 / 347 / 350
- mouse hotspots: 345 / 346 / 349 / 352
- both contexts: on-target
- selection: `best_i_ptm`
- design stages saved by Odin:
  - 50 logits frames
  - 20 soft frames
  - 10 STE frames

The run completed successfully after packaging and environment fixes.

## Selected design

Design:

`human_tfr1_apical_t00000_l60_s190640`

Sequence:

`PPPPLPPAKSPGEAYHRERNLKREAELEAKVAGGEKLTPWEEKELERLKARRAEYMAEAA`

Odin termination:

`LowConfidence`

Selected frame:

- iteration: 78
- stage: STE
- worst-target iPTM: **0.143356**
- worst target: human TfR1

## Full trajectory analysis

The complete 80-frame `trajectory.pickle` was recovered from the 78 MB GitHub
Actions artifact and inspected directly.

### By stage

| Stage | Frames | max human iPTM | max mouse iPTM | max worst-target iPTM | mean human | mean mouse |
|---|---:|---:|---:|---:|---:|---:|
| logits | 50 | 0.169 | 0.215 | **0.157** | 0.118 | 0.127 |
| soft | 20 | 0.150 | 0.143 | 0.143 | 0.111 | 0.115 |
| STE | 10 | 0.143 | 0.153 | 0.143 | 0.128 | 0.135 |

The best joint frame in the entire run occurred very early:

- iteration 5, logits stage
- human iPTM: **0.157**
- mouse iPTM: **0.167**
- worst-target iPTM: **0.157**

The trajectory therefore never entered a strong joint-binding basin at any stage.
The problem is not a late-stage collapse from a previously successful solution.

### Final auxiliary losses

At the final frame:

**Human**

- iPTM: 0.143
- interface-contact loss (`i_con`): ~6.07
- normalized interface-PAE loss (`i_pae`): ~0.737
- pLDDT: ~0.837

**Mouse**

- iPTM: 0.151
- interface-contact loss: ~6.06
- normalized interface-PAE loss: ~0.733
- pLDDT: ~0.814

The two contexts are similarly poor. There is no evidence that one species'
gradient dominated the other. Instead, both contexts failed to establish the
desired interface while the individual-chain confidence remained relatively high.

## Independent Phase 4P validation

The selected Odin sequence was re-docked **from sequence only** against the full
human and mouse TfR1 ectodomains with Boltz-2. Odin's optimized pose was not used
as a template.

### Human

- protein iPTM: **0.255**
- interface: residues approximately 610-650
- preferred apical contacts: none
- apical gate: failed

### Mouse

- protein iPTM: **0.136**
- interface: predominantly residues 609-725
- preferred apical contacts: none
- apical gate: failed

### Cross-species

- homologous overlap count: 2
- Jaccard: 0.0488
- cross-species pose pass: **false**

## Interpretation

This smoke test proves that the Odin-Multi multi-context pipeline works, but the
single random initialization did not find the desired basin.

Three observations drive the next experiment:

1. **The target concept is not unsupported.** FreeBindCraft v6 independently
   achieved human protein iPTM 0.690 with conserved 347/350 contacts retained.

2. **The Odin smoke did not fail because of a human-vs-mouse trade-off.** Both
   contexts were poor throughout the trajectory.

3. **More iterations alone are not justified.** The trajectory did not show an
   improving joint-binding trend that merely needed additional optimization time.

The next experiment should therefore change initialization and the search objective,
not simply extend the same random trajectory.

## Escalation: seeded Odin-Multi v2

The v2 campaign is now active.

### Warm start

The shared 60-aa sequence is initialized from the independently validated
FreeBindCraft v6 parent:

`SPSVISTTTTTYTYHGKQVTRIRTTINPGEKIEPAPGERQMRGRDWETGYIIVHTYKVNS`

The Odin/ColabDesign auxiliary-target forward path uses the same learnable
`params["seq"]`, so the warm-start is genuinely shared across human and mouse
contexts.

### Conserved-core objective

Only the contacts that reproduced in the independent human full-complex prediction
are promoted to core hotspots:

- human: **347 / 350**
- mouse homologues: **349 / 352**

### Loss rebalance

- human gradient weight: 1.0
- mouse gradient weight: 1.25
- interface-contact weight: 1.5
- interface-PAE weight: 0.15
- iPTM weight: 0.15
- ProteinMPNN autoregressive CE: 0.5

The mouse context receives a modestly larger gradient because human already has a
validated binding basin while mouse remains the dominant failure mode.

### Search budget

- 3 independent trajectories
- 60 logits iterations
- 30 soft iterations
- 15 STE iterations
- selection by `best_i_ptm`

## Downstream gate

Phase 4R is connected to the seeded-v2 workflow. It validates all three selected
Odin trajectories, ranked by Odin worst-target iPTM, using independent sequence-only
full-ectodomain Boltz-2 predictions.

A neutral-pH parent must satisfy:

- human protein iPTM >= 0.50;
- mouse protein iPTM >= 0.50;
- preferred/apical contacts in both species;
- homologous human/mouse contact overlap >= 3;
- homologous epitope Jaccard >= 0.15.

No pH-switch claim will be advanced before this neutral-pH cross-species gate passes.
