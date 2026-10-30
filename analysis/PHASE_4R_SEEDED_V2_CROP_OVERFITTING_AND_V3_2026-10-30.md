# Phase 4R — Seeded Odin-Multi reveals crop overfitting; anti-escape v3 escalation

Date: 2026-10-30

## Seeded V2 objective

Seeded V2 warm-started Odin-Multi from the independently validated FreeBindCraft v6
sequence and jointly optimized one 60-aa binder against two cropped apical contexts:

- human TfR1: conserved core hotspots 347 / 350
- mouse TfR1: homologues 349 / 352

The mouse target received a modestly larger gradient (1.25 vs 1.0). Three independent
trajectories were run and selected by `best_i_ptm`.

## Odin-Multi result

All three trajectories converged strongly inside the cropped-apical design problem:

| Design | Odin worst-target iPTM |
|---|---:|
| s767446 | **0.766** |
| s482896 | **0.747** |
| s734966 | **0.738** |

This is a major optimization improvement over the previous random-initialization smoke,
whose best worst-target iPTM across the complete trajectory was only ~0.157.

Therefore the v6 warm-start, conserved-core hotspot reduction and mouse-weight increase
successfully solved the **local two-context optimization problem**.

## Independent Phase 4R full-ectodomain validation

Odin poses were discarded. Each selected sequence was re-docked from sequence only against
the complete human and mouse TfR1 ectodomains with Boltz-2.

### s734966

- human protein iPTM: **0.330**
- mouse protein iPTM: **0.136**
- human interface: distal ~619-650
- mouse: diffuse non-homologous contacts
- homologous overlap: **0**
- pass: false

### s482896

- human protein iPTM: **0.189**
- mouse protein iPTM: **0.132**
- human interface: distal ~615-658
- mouse interface includes distal ~529-706
- homologous overlap: **0**
- pass: false

### s767446

- human protein iPTM: **0.130**
- mouse protein iPTM: **0.142**
- human interface includes 633-640 / 760
- mouse interface includes N-terminal ~153-186 and distal ~518-573
- homologous overlap: **0**
- pass: false

No sequence passed the strict neutral-pH cross-species gate.

## Interpretation: context / crop overfitting

The important observation is the discrepancy between the design context and the independent
full-target context:

- cropped-apical Odin worst-target iPTM: **0.738-0.766**
- full-ectodomain human iPTM: **0.130-0.330**
- full-ectodomain mouse iPTM: **0.132-0.142**
- homologous overlap: **0 for all three**

Thus Seeded V2 did not fail because human and mouse gradients were incompatible. It learned
a sequence that performs well in the restricted cropped-apical world but is not epitope-specific
when the rest of TfR1 is available.

The recurrent full-target escape basins are now empirically identifiable:

- N-terminal mouse surface: approximately 153-188
- distal human surface: approximately 615-658
- distal mouse surfaces: approximately 518-573 and 637-706

The next design objective should therefore expose these competing surfaces during optimization
instead of discovering them only after generation.

## Anti-Escape V3

V3 keeps the validated v6 warm-start and conserved apical targets but introduces explicit
negative design through Odin-Multi's native off-target role.

### On-targets

- human apical: 347 / 350, gradient weight 1.0
- mouse apical: 349 / 352, gradient weight 1.25

### Escape off-targets

Each species receives one off-target context containing two spatially preserved fragments as
separate target chains:

**Human**

- chain A: 122-188 (11QC has no coordinate for residue 121)
- chain B: 500-760

**Mouse**

- chain A: 123-188
- chain B: 500-763

These regions cover both recurrent N-terminal and distal full-target escape modes.

### Negative loss

Per off-target context:

- gradient weight: **0.20**
- interface-contact weight: **-1.0**
- interface PAE weight: **-0.10**
- iPTM weight: **-0.05**
- contact distance: **8 A**
- ProteinMPNN autoregressive CE: **-0.5**
- standard Odin specificity clipping retained

The two desired target contexts still dominate the merged gradient.

### Diagnostic budget

Because the off-target contexts are substantially larger than the 195-aa apical targets,
the first V3 run is deliberately bounded to **one trajectory**. Odin supports extension of
the same immutable run later if this diagnostic improves full-target transfer.

Selection uses `best_clipped_i_pae_ratio`, which balances the weakest target interaction
against the strongest off-target interaction.

## Downstream gate

Phase 4S is connected automatically. It discards the Odin pose and performs independent
sequence-only full human and full mouse Boltz-2 prediction.

The candidate still must satisfy:

- human protein iPTM >= 0.50
- mouse protein iPTM >= 0.50
- apical/preferred contacts retained in both species
- homologous contact overlap >= 3
- homologous-epitope Jaccard >= 0.15

No pH-switch engineering is promoted until this neutral-pH gate is passed.
