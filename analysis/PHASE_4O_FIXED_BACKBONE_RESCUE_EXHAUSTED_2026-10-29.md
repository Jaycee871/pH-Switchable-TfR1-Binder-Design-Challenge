# Phase 4O — fixed-backbone cross-species rescue exhausted

Date: 2026-10-29.

## Phase 3J / 4N: two-position cleanup

The FreeBindCraft v6 parent was the first relaxed trajectory with genuine conserved
hotspot engagement and zero severe clashes. Its remaining human species-risk contacts
were concentrated around binder B46/B47, so a strict ProteinMPNN cleanup fixed all other
58 binder positions and redesigned only B46/B47.

Three valid variants were tested:

- SS
- SN
- ST

All three retained strong human apical binding. The best human protein iPTM values were
~0.70-0.71.

However, all three mouse predictions left the homologous apical surface:

- SS: mouse protein iPTM 0.467
- SN: mouse protein iPTM 0.492
- ST: mouse protein iPTM 0.331

All three had zero homologous human/mouse epitope overlap.

**Conclusion:** B46/B47 chemistry is not the root cause of mouse pose collapse.

## Phase 3K / 4O: tied human+mouse dual-state ProteinMPNN

To test sequence-level multistate design more directly, two structural states were
combined in one ProteinMPNN input:

- A+B = human apical state + binder
- C+D = mouse homologous apical state + binder

The two complexes were spatially separated. Every binder position was tied:

`B_i == D_i, i = 1..60`

Only interface positions were designable; non-interface positions remained fixed to the
relaxed v6 parent.

An initial prefilter incorrectly rejected all 30 sequences because the parent itself
contains a fixed >4-residue homopolymer run. The filter was corrected to reject only a
candidate-introduced worsening. No ProteinMPNN rerun was required.

Top four sequences were validated independently with full human and mouse TfR1 Boltz-2.

## Independent results

### MS_08

- human protein iPTM: **0.583**
- mouse protein iPTM: **0.162**
- homologous overlap: **0**

### MS_38

- human protein iPTM: **0.598**
- mouse protein iPTM: **0.162**
- homologous overlap: **0**

### MS_15

- human protein iPTM: **0.734**
- human hotspot contacts: **347, 350**
- mouse protein iPTM: **0.251**
- homologous overlap: **0**

### MS_16

- human protein iPTM: **0.617**
- mouse protein iPTM: **0.192**
- homologous overlap: **0**

No candidate passed the cross-species pose gate.

## Interpretation

The failure is now demonstrated at three sequence-design levels:

1. post-hoc single-state ProteinMPNN;
2. narrow two-position liability cleanup;
3. tied human/mouse dual-state ProteinMPNN.

In every case, the v6 backbone remains highly designable for human TfR1 but does not
support a stable homologous mouse apical pose under independent full-complex prediction.

The bottleneck is therefore **backbone / docking geometry**, not a small set of binder
side chains.

## Decision

**Close the v6 fixed-backbone sequence-rescue family.**

Do not spend additional ProteinMPNN or Boltz-2 budget on incremental variants of this
backbone.

The next generator must optimize backbone/pose and sequence under human and mouse target
contexts simultaneously.

## Escalation: Odin-Multi

Odin-Multi is selected because its native cross-reactivity mode evaluates the same binder
sequence independently against multiple on-target structures inside every design
trajectory, then merges the target gradients to update the shared sequence/pose.

The first TfR1 smoke uses:

- human TfR1 apical context: hotspots 343/344/347/350;
- mouse TfR1 apical context: hotspots 345/346/349/352;
- both contexts as on-targets;
- fixed 60-aa binder;
- one joint trajectory;
- 8 A contact objective;
- `best_i_ptm` selection (maximize the weakest on-target iPTM);
- independent full human/mouse Boltz-2 as the downstream gate.

No pH-switch engineering is justified until this neutral-pH cross-species parent exists.
