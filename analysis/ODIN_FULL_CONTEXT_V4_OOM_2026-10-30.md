# Odin full-context v4 — GPU memory failure analysis

Date: 2026-10-30.

## Intended experiment

V4 was the first Odin-Multi campaign to expose the optimizer directly to the
complete human and mouse TfR1 ectodomains:

- human: residues 122-760, hotspots 347 / 350;
- mouse: residues 123-763, homologous hotspots 349 / 352;
- shared 60-aa v6 warm-start;
- human gradient weight 1.0;
- mouse gradient weight 1.25;
- AF2 multimer design family;
- one bounded trajectory.

The scientific question was whether eliminating context truncation would prevent
the distal-surface escape seen in the cropped-context Odin campaigns.

## Result: no scientific trajectory was produced

The Kaggle kernel completed environment setup successfully:

- Odin-Multi installed;
- pinned Odin commit verified;
- full human and mouse target structures validated;
- JAX detected the Kaggle T4 GPU;
- v6 warm-start patch applied;
- human/mouse loss configuration loaded.

The failure occurred immediately after entering:

`Stage 1: Test Logits`

Before a usable design frame was saved, JAX/XLA reported:

`RESOURCE_EXHAUSTED: Out of memory while trying to allocate 36422684968 bytes`

Requested allocation:

**36,422,684,968 bytes (~33.9 GiB / 36.4 GB decimal)**

The available Kaggle T4 has approximately 16 GB device memory.

The design process exited with code 1 and no `trajectory.pickle`,
`selection.csv`, or selected sequence was produced.

## Interpretation

This is an **engineering / hardware feasibility failure**, not evidence that
full-ectodomain joint optimization is scientifically ineffective.

The experiment never entered a valid AF2 optimization trajectory, so it provides
no information about:

- human vs mouse gradient compatibility;
- full-context hotspot localization;
- cross-species iPTM;
- sequence quality;
- selector behavior.

Phase 4U and Phase 4V were therefore correctly skipped.

Phase 5A / 5B workflow runs that appear green after this failure are gate-only
no-op workflows: no cross-species parent or pH-switch candidate was generated.

## Next diagnostic: V4B

The next experiment retains the complete ectodomain contexts and all scientific
target definitions but changes the AF2 design family:

- `use_multimer_design: true -> false`;
- monomer/PTM binder model family;
- `sample_models: false`;
- one trajectory;
- only 3 logits / 2 soft / 1 STE iterations.

V4B is strictly a memory-feasibility probe.

### Decision rule

If V4B also OOMs before producing a trajectory:

**stop attempting full-ectodomain gradient backpropagation on Kaggle T4.**

The next architecture should use cropped/apical gradient generation coupled to
full-ectodomain inference or rescoring in the loop, rather than continuing to
shrink an intrinsically over-memory gradient model.

If V4B succeeds:

use its measured runtime and memory feasibility to define a bounded V4C
full-context monomer/PTM campaign, followed by independent full-complex Boltz-2.
