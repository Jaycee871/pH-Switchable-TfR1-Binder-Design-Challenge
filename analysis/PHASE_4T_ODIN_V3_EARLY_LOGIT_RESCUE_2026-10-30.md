# Phase 4T — Odin v3 early-logit rescue

Date: 2026-10-30.

## Question

The anti-escape v3 selector chose STE frame 73, even though the complete trajectory
contained substantially better joint target/off-target metrics during the logits stage.

To test whether the failure was merely a selector/discretization artifact, the four
best early logits frames were converted to discrete sequences by argmax over the
saved 20-AA logits using ColabDesign/AlphaFold's internal residue order:

`ARNDCQEGHILKMFPSTWYV`

Frames tested:

- 39: worst-target Odin iPTM 0.764, clipped specificity ratio 1.900
- 37: worst-target Odin iPTM 0.761, clipped specificity ratio 1.806
- 31: worst-target Odin iPTM 0.715, clipped specificity ratio 1.557
- 36: worst-target Odin iPTM 0.713, clipped specificity ratio 1.625

All four were independently re-docked from sequence only against the full human
and mouse TfR1 ectodomains with Boltz-2.

## Result: no cross-species parent

### F39

- human protein iPTM: **0.206**
- human preferred contacts: none
- human escapes predominantly to 619-760
- mouse protein iPTM: **0.223**
- mouse preferred contacts: **345, 346, 349, 352**
- mouse nevertheless makes a very broad multi-surface interface
- homologous overlap: 1, Jaccard 0.015

This is the most informative rescue: the mouse sequence can still touch all four desired
apical residues, but the full-complex interaction is globally weak and non-specific while
the human state escapes distally.

### F37

- human protein iPTM: **0.330**
- mouse protein iPTM: **0.135**
- homologous overlap count: **9**
- Jaccard: **0.141**

The apparent cross-species overlap is on the wrong distal surface around 627-654, not
the desired apical epitope.

### F31

- human protein iPTM: 0.150
- mouse protein iPTM: 0.138
- homologous overlap: 1
- no desired cross-species pose

### F36

- human protein iPTM: 0.311
- mouse protein iPTM: 0.178
- homologous overlap: 0
- no desired cross-species pose

## Interpretation

The failure of anti-escape v3 is **not** explained by the official selector choosing
a suboptimal late STE frame.

The trajectory contains a high-scoring cropped-context basin, but discretizing its best
early frames does not transfer that basin to the full ectodomain.

This establishes a stronger form of context overfitting:

1. cropped positive contexts can report high joint human/mouse iPTM;
2. cropped decoy contexts can report strong apparent off-target rejection;
3. even the best high-scoring early-logit sequences still escape or become non-specific
   in full-ectodomain prediction.

Therefore adding more cropped escape decoys is not justified.

## Decision

Close the cropped-positive + cropped-decoy Odin family.

The next diagnostic must expose Odin directly to the **full ectodomain structural
context during optimization**, so the distal surfaces cannot disappear from the model's
world.

Because full ectodomain contexts are substantially larger, the next campaign should be
bounded to one warm-started trajectory and a reduced iteration budget sized to Kaggle's
time limit. Independent sequence-only full human/mouse Boltz-2 remains the downstream gate.
