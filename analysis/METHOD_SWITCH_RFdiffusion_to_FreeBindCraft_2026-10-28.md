# Method switch: RFdiffusion exploration → FreeBindCraft pose-optimised design

**Date:** 2026-10-28  
**Challenge:** Anthropic × Adaptyv 2026 — TfR1 / CD71  
**Design objective:** human TfR1 ON at pH 7.4, mouse TfR1 ON at pH 7.4, human TfR1 OFF/reduced at pH 6.0.

## Executive decision

RFdiffusion is retained as an exploration/history lane but is no longer the primary
backbone generator for this target.

The project switches to a pose-optimising binder generator (FreeBindCraft first)
because repeated RFdiffusion campaigns show that:

1. compact backbones can be generated;
2. conserved apical contacts can be generated;
3. zero-risk contacts can be generated;
4. but these properties are not reliably achieved **together**, and sequence-designed
   RFdiffusion scaffolds do not preserve the intended human/mouse full-ectodomain pose.

This is a method-selection decision based on accumulated negative results, not an
assumption that RFdiffusion is generally unsuitable for protein-binder design.

---

## Evidence from the RFdiffusion lane

### Old compact seed / C57 lane

The seed-20261102 scaffold repeatedly collapsed away from the intended apical epitope
after full-complex validation.

The explicit V2H rescue experiment did not change that conclusion:

| Candidate | Human protein ipTM | Mouse protein ipTM | Preferred apical core retained |
|---|---:|---:|---|
| V2H_03 | 0.183 | 0.169 | No |
| V2H_05 | 0.143 | 0.147 | No |

The predicted interfaces relocated to the distal ~507–657 surface.

**Decision:** close the seed-20261102 / C57 primary lane.

### Kaggle RFdiffusion pilot

Six backbones were generated successfully. The original workflow was marked failed
only because its final Git push was non-fast-forward; the compute itself succeeded.

After correcting target-residue renumbering in the scientific gate:

- strict pass: **0 / 6**
- near-miss `cross_species_edge_5/design_0`:
  - preferred 343/347
  - zero explicit risk contacts
  - Rg 18.812 Å
  - CA span 64.808 Å
- near-miss `functional_patch_mixed_5/design_0`:
  - preferred 343/347
  - risk 210
  - Rg 17.427 Å
  - CA span 68.065 Å

The pilot therefore demonstrated useful epitope sampling but insufficient compactness.

### NVIDIA apical-only Phase 2D

Target was restricted to 11QC chain A residues 189–383.

One of six backbones passed the local structural gate:

- family: conserved_343
- seed: 20261124
- binder length: 86 aa
- preferred contacts: 209, 347
- risk contact: 210
- Rg: 17.05 Å
- CA span: 54.22 Å
- severe clashes: 0

This was the strongest RFdiffusion scaffold and was carried into sequence design.

### Pose-first ProteinMPNN + full human/mouse validation

For the Phase-2D scaffold, only one of eight initial ProteinMPNN sequences passed the
composition prefilter.

Candidate `TFR1_NB01_conserved_343_s20261124_08`:

- human protein ipTM: **0.202**
- human contacts retained: 207, 209, 214, 344, 374
- human risk contacts: 208, 210, 211, 212
- mouse protein ipTM: **0.173**
- mouse preferred apical contacts: none
- mouse pose: collapsed to distal ~614–728 surface

This was an important partial result: the redesigned scaffold could remain on the
human apical surface, but not on the homologous mouse surface.

### Cross-species tolerant sequence rescue

Binder positions 74–76 lie across the human L209 / V210 region, where V210 corresponds
to a chemically different mouse D213 environment. A narrow ProteinMPNN rescue biased
those positions toward neutral/tolerant residues (A/G/Q/N/S/T).

Two sequences passed the composition gate:

- `TFR1_XS343_14`: 74–76 fully tolerant; charged fraction 0.523
- `TFR1_XS343_12`: two of three tolerant; charged fraction 0.465

Full human/mouse validation still failed:

| Candidate | Human protein ipTM | Mouse protein ipTM | Homologous apical pose |
|---|---:|---:|---|
| XS343_14 | 0.120 | 0.212 | Fail |
| XS343_12 | 0.130 | 0.123 | Fail |

**Decision:** close the sequence-rescue lane. The failure is not explained by a single
74–76 chemical mismatch.

### Phase 2E: zero-risk 343-cluster generation

Six 50–70 aa backbones were sampled with hotspots 343/344/347/350.

Strict pass: **0 / 6**.

The result separated the design requirements:

- seeds 20261231 / 20261232:
  - preferred 343/347
  - zero risk contacts
  - too extended (Rg 25.28 / 24.63 Å)
- seed 20261234:
  - preferred 347/350
  - zero risk
  - too extended (Rg 27.56 Å)
- seed 20261236:
  - 55 aa
  - zero risk
  - Rg 12.04 Å
  - span 43.37 Å
  - but only one preferred anchor (347)

Thus dual-anchor chemistry and compactness were both possible, but not together.

### Phase 2F: compact zero-risk refinement

Binder length narrowed to 45–60 aa and hotspot pairs/triplets were sampled.

Strict pass: **0 / 6**.

Best compact example:

- seed 20261312
- 55 aa
- Rg 12.906 Å
- span 40.508 Å
- zero risk
- only preferred anchor 347

The shorter binder did not solve hotspot control.

### Phase 2G: local 330–383 epitope window

The target was reduced further to remove the previously observed escape surfaces near
210 and 291.

An initial analysis bug incorrectly inferred target/binder chain roles when both chains
were 54 aa. Re-analysis identified the gly-rich generated chain as binder and the
native-sequence chain as target.

Corrected strict pass: **0 / 6**.

Several backbones were compact (Rg ~11–13 Å, span ~35–50 Å), but the interfaces drifted
toward residues ~356–367 rather than the requested 343–350 cluster.

**Interpretation:** shrinking the target improved compactness but did not make RFdiffusion
hotspot control reliable enough for this particular design problem.

---

## Why FreeBindCraft is the next primary lane

FreeBindCraft optimises binder-target geometry inside an AlphaFold2 hallucination loop,
then performs ProteinMPNN refinement and AF2 re-validation. This changes the role of
structure prediction from a downstream filter to part of the design objective.

Initial smoke configuration:

- target: 11QC chain A residues 189–383
- target hotspots: 343, 344, 347, 350
- binder range: 55–75 aa (future strict smoke should use a fixed 60 aa shape)
- PyRosetta: disabled
- final accepted-design target: 1
- rank metric: ipSAE
- relaxed filters for environment/feasibility smoke only
- full human/mouse Boltz-2 validation remains mandatory

### Kaggle smoke v1

Environment:
- 2 × Tesla T4
- 15,360 MiB per GPU
- NVIDIA driver 580.178.04
- host CUDA 13.0

The FreeBindCraft Python/CUDA stack and ~5.3 GB AlphaFold parameter archive installed,
but the installer ended at its own `nvidia-smi` self-check.

This was an infrastructure failure, **not** a model or VRAM failure.

### Kaggle smoke v2

The installer was run with its documented `--skip-gpu-check` option and an explicit
JAX device probe was added.

Installation completed successfully with:

- jax 0.11.2
- jaxlib 0.11.2
- jax-cuda12-plugin 0.11.2
- flax 0.12.10
- colabdesign 1.1.3
- OpenMM 8.6.1

The JAX probe then failed at `cuInit(0)` with CUDA error 303.

Inspection of FreeBindCraft's conda activation script showed that it deliberately
unsets `LD_LIBRARY_PATH`. On Kaggle, that also hides the host-mounted NVIDIA driver
libraries.

This was again an environment-library visibility failure, **not** an OOM result.

### Kaggle smoke v3

The runner now restores only NVIDIA **driver** library directories while retaining
FreeBindCraft's pinned pip CUDA 12.9 user-space runtime.

The current v3 run has already exceeded both v1 and v2 failure times, consistent with
having progressed beyond the earlier installer/JAX blockers. Final confirmation must
come from the completed Kaggle log.

---

## Non-negotiable downstream gate

A FreeBindCraft-accepted structure will not be treated as a competition candidate until
it passes:

1. full human TfR1 Boltz-2 prediction;
2. human protein ipTM >= 0.50;
3. conserved apical contact retained;
4. full mouse TfR1 Boltz-2 prediction;
5. mouse protein ipTM >= 0.50;
6. homologous human/mouse epitope overlap;
7. only then: pH-switch engineering;
8. explicit pH 6.0 negative-state evaluation.

This ordering is intentional: **pose first, pH chemistry second**.
