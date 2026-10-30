# Full-context inference-loop audit through IL4 and escalation to mined hard negatives

Date: 2026-10-31 (Asia/Taipei)

## Why the interpretation changed

The project originally treated the FreeBindCraft v6 human full-complex prediction
(protein iPTM 0.690 with 347/350 contacts) as evidence that the v6 sequence itself
encoded a strong human apical binder.

The later **IL0 sequence-only control** showed that interpretation was too strong.

When the identical v6 sequence was re-docked without a structural template:

- human protein iPTM: **0.318**
- mouse protein iPTM: **0.380**
- human preferred/apical contacts: none
- mouse preferred/apical contacts: none
- both predictions moved to distal surfaces

Therefore the earlier 0.690 human result contains substantial **pose/template
conditioning**. It is useful as a proposal state, but it is not proof that the
sequence independently selects the apical epitope.

All subsequent inference-loop rounds are interpreted against the sequence-only
full-ectodomain control.

## IL1 — broad full-context tied ProteinMPNN

Proposal states:

- human: Phase 4Q v6 full-human pose;
- mouse: F39 full-mouse pose that contacts 345/346/349/352.

The union of the two interface sets (29 binder positions) was redesigned with a
single tied sequence.

72 unique sequences were generated. Top four changed 12-15 residues.

Independent sequence-only full-ectodomain validation:

- best human iPTM: **0.354**
- best mouse iPTM: **0.237**
- no candidate retained the human preferred apical core
- homologous overlap remained negligible

Interpretation: redesigning the whole inferred interface destroyed too much of
the parent sequence and did not encode epitope specificity.

## IL2 — sparse seven-position reconciliation

The design space was reduced to the seven v6/F39 sequence disagreements:

10, 14, 42, 46, 47, 51, 52

All other 53 positions were fixed.

Top results:

- R2_99: human 0.225 / mouse 0.262
- R2_58: human 0.181 / mouse 0.169
- R2_43: human 0.292 / mouse 0.142
- R2_72: human **0.473** / mouse 0.130

R2_58 recovered several mouse apical contacts including 345/346/349, but the
human prediction still escaped. No sequence passed.

Interpretation: sparse reconciliation preserved more parent behavior than IL1,
but no combination among the seven differences simultaneously encoded the human
and mouse apical epitope.

## IL3 — one-position causal flips

Single substitutions were tested to isolate which v6/F39 differences drive the
full-context behavior.

Notable results:

- **R42P**: human 0.384 / mouse 0.360, overlap count 3
- **E47D**: mouse 0.441
- **W46A**: mouse 0.390
- **T10E**: human 0.271 / mouse 0.346

However, the improved scores were still generated on non-apical surfaces.
R42P is the most balanced score-wise but has essentially no desired apical
engagement.

Interpretation: some single positions modulate global docking confidence, but
none is sufficient to encode the target epitope.

## IL4 — human-anchor encoding test

IL4 asked a narrower question: can ProteinMPNN encode the known human proposal
pose into sequence before asking for mouse cross-reactivity?

A larger human-interface design was generated and human sequence-only predictions
were screened first.

Best result:

- IL4_246 human protein iPTM: **0.574**

But:

- preferred contacts: none
- apical fraction: **0**
- interface moved to the distal 600-series surface

No IL4 sequence passed the human anchor gate, so the downstream mouse step was
correctly skipped.

Interpretation: a high protein iPTM alone is not enough. ProteinMPNN can create a
confident binder for TfR1 while changing the epitope completely.

## Core conclusion through IL4

The dominant failure mode is now **epitope specificity**, not raw complex confidence.

Across IL0-IL4, sequence-only full-ectodomain predictions repeatedly find
alternative TfR1 surfaces, especially:

### Human recurrent escape basins

- approximately 606-678
- secondary N-terminal basin around 145-187

### Mouse recurrent escape basins

- approximately 516-582
- approximately 629-680

These locations recur across multiple independent sequences and mutation sets.

This explains why:

- cropped-apical Odin can reach high internal iPTM;
- fixed-backbone ProteinMPNN can score well on proposal structures;
- template-conditioned human predictions can look strong;

while the same sequence re-docks elsewhere when the complete receptor is exposed.

## Why Anti-Escape V3 was not sufficient

V3 already used negative design, but its off-target contexts were very broad:

- human 122-188 + 500-760
- mouse 123-188 + 500-763

and each species-level escape context had gradient weight only 0.20.

The negative signal was therefore distributed across hundreds of residues instead
of being concentrated on the actual recurrent escape patches.

## Escalation: mined hard-negative Odin V5

V5 is the first generator round whose negative contexts are mined directly from
the sequence-only full-complex errors above.

### Positive contexts

- human apical 189-383; hotspots 347/350; gradient 1.0
- mouse apical homolog; hotspots 349/352; gradient 1.25

### Localized hard negatives

- human 606-678; gradient 0.55
- human 145-187; gradient 0.25
- mouse 516-582; gradient 0.45
- mouse 629-680; gradient 0.45

Each escape patch is a separate Odin off-target context.

### Dual warm-start experiment

Two deterministic trajectories test two known basins under the same objective:

1. v6 sequence — human-side proposal parent
2. F39 sequence — mouse full-complex sequence with observed 345/346/349/352 contacts

This isolates whether the same mined-negative objective behaves differently when
initialized from the human-side versus mouse-side basin.

### Validation

Odin internal scores are not accepted as candidate evidence.

Both V5 selections must be independently re-docked **sequence-only** against the
complete human and mouse ectodomains in Phase 4W.

A real parent still requires:

- human protein iPTM >= 0.50
- mouse protein iPTM >= 0.50
- preferred/apical engagement in both species
- homologous epitope overlap >= 3
- Jaccard >= 0.15

No pH-switch work is promoted unless this neutral-pH gate passes.
