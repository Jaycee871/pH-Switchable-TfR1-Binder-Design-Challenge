# Phase 4K — targeted cross-species rescue of FreeBindCraft l67

Date: 2026-10-29.

Parent: `TFR1_FBC343_l67_s385919`.

Phase 4I established that the parent has a real human apical contact at 231/350,
with R10 providing the cleanest conserved E350/G351 anchor, but the mouse
full-complex prediction collapses away from the apical domain.

A narrow ProteinMPNN rescue therefore redesigned the main sequence liabilities
(W7/Y8/Y11; secondary H14/K18) while preserving R10.

Sixteen sequences passed the composition/anchor prefilter. The two best
low-MPNN-score candidates were validated with the unchanged full human/mouse gate.

## Result: FAIL

### TFR1_FBC_L67_XS_07

- human protein ipTM: **0.364**
- human preferred contacts: none
- mouse protein ipTM: **0.183**
- mouse apical contact: false
- homologous overlap: **0**
- cross-species pose pass: false

The rescue lost the parent's useful human 231/350 epitope.

### TFR1_FBC_L67_XS_05

- human protein ipTM: **0.397**
- human preferred contacts: **231, 350**
- human risk contacts: none
- mouse protein ipTM: **0.203**
- mouse apical contact: false
- homologous overlap: **0**
- cross-species pose pass: false

This candidate preserves the parent's human apical geometry better than XS_07,
but does not improve cross-species stability. Mouse again relocates to distal
surfaces rather than the homologous apical patch.

## Interpretation

The l67 mouse-collapse phenotype is not rescued by neutralising the obvious
non-conserved peripheral contacts while holding R10 fixed.

That negative result matters: the cross-species failure is likely encoded in
a broader backbone/interface geometry or cooperative interaction pattern, not
a small set of W7/Y8/Y11/H14/K18 side-chain liabilities.

## Decision

**Close the l67 sequence-rescue lane.**

Do not spend additional Boltz-2 budget on incremental ProteinMPNN variants of
this parent.

The active primary route is now the FreeBindCraft v5 binary hard-hotspot
experiment, which changes the generator objective itself rather than trying to
repair an already unstable cross-species parent.

No pH engineering is justified until a parent first passes the full neutral-pH
human/mouse pose gate.
