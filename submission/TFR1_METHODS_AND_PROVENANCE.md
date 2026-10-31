# TfR1 Challenge 04 — Methods, ranking, and provenance

## Submission snapshot

This package contains a ranked 20-design portfolio for the 2026 Anthropic × Adaptyv Proteinbase TfR1 challenge.

All submitted sequences are de novo computational protein designs. No experimentally validated TfR1 binder sequence was used as a starting sequence. Known TfR1 binders and viral complexes were used only to define and calibrate structural epitope constraints. Iterative refinement and warm-starting used only de novo sequences generated inside this challenge project.

Repository: https://github.com/Jaycee871/pH-Switchable-TfR1-Binder-Design-Challenge

## Design objective

The intended molecular phenotype is:

1. bind human TfR1 at pH 7.4;
2. bind mouse TfR1 at pH 7.4;
3. reduce or lose human TfR1 binding at pH 6.0.

The preferred design region was the TfR1 apical domain. Human core residues 347/350 and homologous mouse residues 349/352 were used repeatedly as conserved structural anchors. Species-variable positions around the known apical patch were explicitly tracked.

## Methods used

The portfolio deliberately spans several independent computational families instead of ranking by a single model score.

### 1. Backbone generation and hotspot design

- RFdiffusion was used in early backbone-generation experiments.
- FreeBindCraft / AF2-based hallucination was used to generate human TfR1 apical binders.
- Binary hotspot losses were tested to localize the interface to the 343–350 region.
- ProteinMPNN was used for sequence design and constrained sequence rescue.

### 2. Human–mouse multistate design

Several increasingly strict approaches were tested:

- post-hoc mouse validation of human-designed binders;
- sparse ProteinMPNN sequence rescue;
- a tied-state human+mouse ProteinMPNN objective;
- Odin-Multi joint optimization with human and mouse as simultaneous on-target contexts.

### 3. Full-context validation

A major finding was that cropped apical-domain design scores did not reliably transfer to the full ectodomain.

Candidates were therefore re-docked against full human and mouse ectodomains with Boltz-2. Where possible, sequence-only full-ectodomain prediction was used so that the design pose could not simply be reproduced by template conditioning.

This uncovered a recurrent failure mode: apparently strong apical-domain designs often escaped to distal TfR1 surfaces in the complete receptor.

### 4. Inference-in-the-loop negative design

Full-ectodomain prediction failures were treated as new training/design information.

Recurrent escape regions included approximately:

- human 606–678;
- human 145–187;
- mouse 516–582;
- mouse 629–680.

The final Odin V5 campaign used these empirically mined surfaces as explicit localized off-target contexts while retaining human and mouse apical positive contexts. Two trajectories were initialized from different de novo basins: a human-side FreeBindCraft sequence and a mouse-apical Odin sequence.

### 5. pH mechanism

The challenge is explicitly pH dependent. Histidine-containing chemistries were tracked because protonation near endosomal pH can provide a plausible conditional-binding mechanism.

However, this submission does **not** claim that an in-silico pH switch has been proven. A planned PROPKA / pKa-ranking stage was intentionally gated behind a strict neutral-pH cross-species parent criterion, and no candidate cleared that internal gate. The submitted portfolio therefore preserves sequence diversity and histidine-content diversity for experimental testing rather than reporting an unsupported pH-selectivity prediction.

## Internal validation policy

A deliberately strict internal cross-species pose gate required:

- human protein ipTM >= 0.50;
- mouse protein ipTM >= 0.50;
- desired/apical engagement in both species;
- homologous contact overlap >= 3;
- Jaccard >= 0.15.

No final candidate satisfied all of these requirements simultaneously.

This negative result is retained in the provenance rather than hidden. The ranked portfolio is therefore an uncertainty-aware experimental shortlist, not a claim that the three challenge objectives are computationally solved.

## Why these 20 designs were selected

The final ranking is an ensemble judgement over:

- full-ectodomain human and mouse confidence;
- evidence of apical or preferred contacts;
- homologous cross-species contact overlap;
- sequence-only versus pose-conditioned evidence;
- generator / backbone diversity;
- negative-design provenance;
- histidine content and sequence liabilities;
- avoidance of submitting only near-identical variants.

Examples of complementary evidence:

- **TF4_001_L70** has the strongest full-target mouse ipTM in the portfolio (0.677) and 10 homologous contact overlaps, although the predicted epitope is non-apical.
- **TF4_002_XS08** retains a strong human apical pose (human ipTM 0.703; contacts 231/347/350) and reaches mouse ipTM 0.492, but the mouse pose is non-apical.
- **TF4_003_F39** directly contacts mouse 345/346/349/352 in a sequence-only full-target prediction, providing the strongest qualitative mouse-apical evidence.
- **TF4_004_V5_T0** and **TF4_007_V5_T1** represent the final mined-hard-negative method.
- **TF4_005_R42P** is a sparse causal variant with relatively balanced sequence-only human/mouse scores.
- **TF4_013_L67**, **TF4_014_L67_XS05**, and **TF4_020_L67_XS07** provide a structurally distinct 67-aa family rather than only variants of the dominant 60-aa scaffold.

## Important interpretation notes

1. Protein ipTM values are model outputs, not experimental affinities.
2. Some early FreeBindCraft validation is pose-conditioned; this is marked in the metadata file.
3. Sequence-only redocking often produced lower scores and different epitopes, which is why those results were given substantial weight.
4. Functional/apical binding was preferred during ranking, but non-apical candidates were not automatically discarded because the challenge assays binding experimentally and the official challenge does not require the apical epitope.
5. All experimental conclusions remain to be determined by the Proteinbase wet-lab assay.

## Files

- `tfr1_submission_upload_v0.csv` — minimal upload-safe CSV with only required columns.
- `tfr1_submission_metadata_v0.csv` — ranking evidence and provenance by design.
- `TFR1_METHODS_AND_PROVENANCE.md` — this methods note.
- `TFR1_SUBMISSION_CHECKLIST.md` — final review steps before upload.

