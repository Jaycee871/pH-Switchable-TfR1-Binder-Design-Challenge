# TfR1 Challenge 04 — Official methodology statement

Author: Pack Kwan Low (Jaycee), GitHub `Jaycee871`.
Date: 2026-10-31 (Asia/Taipei).
Status: This is the official methodology statement for the project, as written by the author
and submitted alongside the Challenge 04 TfR1 entry. The text in the section below is
reproduced verbatim and has not been edited, reworded, or summarised.

---

## Methodology statement (verbatim, as written by the author)

We used an iterative, structure-guided design pipeline to generate de novo TfR1 binders with three objectives: binding to human TfR1 at physiological pH, cross-reactivity with mouse TfR1, and a plausible route toward reduced binding under acidic endosomal conditions.

Our workflow combined RFdiffusion, FreeBindCraft, ProteinMPNN, Odin-Multi, and independent Boltz-2 structure prediction. Human and mouse TfR1 were first aligned structurally to identify a conserved apical binding region, with human residues 347/350 and the homologous mouse residues 349/352 used as recurring structural anchors. Designs were generated using both single-target and multistate objectives, including tied human–mouse sequence optimization.

A key part of our methodology was independent full-ectodomain validation. We found that high-confidence designs generated against cropped apical domains could sometimes re-dock to unrelated distal surfaces when evaluated against the complete receptor. Rather than ignoring these failures, we used them as negative-design information. Recurrent full-context escape regions were identified from sequence-only Boltz-2 predictions and fed back into later Odin-Multi rounds as explicit off-target contexts. The final generation stage therefore combined positive human/mouse apical objectives with empirically mined negative surfaces.

The submitted set is intentionally a diverse experimental portfolio rather than a collection of near-identical top-scoring variants. Designs were selected using multiple signals, including human and mouse full-ectodomain confidence, conservation of desired contacts, homologous cross-species contact overlap, sequence diversity, method diversity, and known sequence liabilities. We did not rank designs using a single confidence metric.

We also explicitly distinguished pose-conditioned predictions from sequence-only re-docking. Sequence-only full-receptor predictions were given greater weight when assessing whether a sequence genuinely encoded its intended binding preference.

For the pH-switch objective, histidine-containing chemistries were retained as a mechanistic design feature because histidine protonation can change substantially between physiological and endosomal pH. However, we do not claim that pH-selective release has been computationally proven. The submitted sequences should therefore be viewed as uncertainty-aware candidates for experimental testing rather than as computationally validated pH switches.

All submitted sequences were generated de novo within this project. Known TfR1 binders and structural complexes were used only to define structural constraints and epitope hypotheses; experimentally validated binder sequences were not used as starting sequences.


---

## Maintainer's fact-check appendix (NOT part of the author's statement)

> **This appendix is an annotation added by the repository maintainer on 2026-10-31
> (Asia/Taipei). It is not part of, and does not modify, the author's statement above.
> It records which claims in the statement could be verified against files in this
> repository, which could not, and where the repository adds material context.**

### Verified against the repository

- **FreeBindCraft, ProteinMPNN, Odin-Multi and Boltz-2 were all used.** Confirmed by
  `submission/tfr1_submission_metadata_final_v1.csv` (`method_family` column), the
  `odin_multi/` and `freebindcraft/` trees, and the Boltz-2 pose screens under
  `predictions/`. Of the 16 submitted designs: 5 are FreeBindCraft-derived
  (TF4_001, TF4_002, TF4_006, TF4_013, TF4_014), 6 are Odin-derived
  (TF4_003, TF4_004, TF4_007, TF4_009, TF4_015, TF4_016), 4 are inference-loop
  (TF4_005, TF4_010, TF4_017, TF4_019), and 1 is tied human+mouse ProteinMPNN (TF4_012).
- **Human 347/350 and mouse 349/352 as recurring anchors.** Confirmed. Stated in
  `submission/TFR1_METHODS_AND_PROVENANCE.md`, used as the explicit hotspot set in the
  Odin-Multi seeded v2 campaign described in `README.md`, and visible in the
  `human_preferred_contacts` column (`231;347;350` for TF4_002, TF4_006; `347;350` for
  TF4_012) and `mouse_preferred_contacts` (`345;346;349;352` for TF4_003).
- **Tied human-mouse sequence optimization was actually run.** Confirmed. TF4_012_MS15
  carries `method_family = Tied-human-mouse_ProteinMPNN`, and `README.md` documents it as
  Phase 3K / 4O, with all 60 binder positions tied across explicit human and mouse states.
- **Mined negative surfaces were fed back into Odin-Multi.** Confirmed, and specifically.
  `kaggle/odin_multi_hard_negative_v5/` embeds four escape decoys named
  `human_escape_606_678`, `human_escape_145_187`, `mouse_escape_516_582` and
  `mouse_escape_629_680`, passed to the optimizer as explicit `--context` off-targets with
  per-decoy gradient weights. These match the recurrent escape regions recorded in
  `TFR1_METHODS_AND_PROVENANCE.md` exactly. The earlier broad anti-escape V3 campaign is
  present at `kaggle/odin_multi_anti_escape_v3/`. TF4_004 and TF4_007 come from the V5
  mined-hard-negative stage; TF4_015 comes from the V3 broad anti-escape stage.
- **Each of the 16 submitted designs traces to a stated lineage.** Confirmed. Every row of
  `tfr1_submission_metadata_final_v1.csv` carries a `source_id` and a `method_family`, and
  the four designs dropped between v0 (20) and final-v1 (16) are itemised with reasons in
  `submission/TFR1_FINAL_V1_FREEZE_NOTE.md`.
- **Pose-conditioned vs sequence-only evidence was tracked separately.** Confirmed. The
  `evidence_scope` column distinguishes `pose-conditioned full-ectodomain Boltz-2` from
  `sequence-only full-ectodomain Boltz-2`, and TF4_006 additionally records a contradicting
  sequence-only control (`h=0.318, m=0.380`) against its pose-conditioned human iPTM of 0.690.
- **No known binder used as a starting sequence.** Consistent with
  `TFR1_METHODS_AND_PROVENANCE.md`. Note this is corroboration from another document in the
  same repository, not independent verification; nothing in the repository contradicts it.

### Partly verified, or needing qualification

- **RFdiffusion.** RFdiffusion was used in the project - `structures/rfdiffusion/` contains
  several campaigns (`nvidia_apical_only`, `nvidia_local_343_window`, `seeded_selection`,
  `nvidia_zero_risk_343`, `nvidia_zero_risk_compact_refine`, `smoke_001`) with selection
  summaries. **However, none of the 16 submitted designs has an RFdiffusion `method_family`.**
  The statement's "our workflow combined RFdiffusion..." is accurate about the project as a
  whole, but no submitted sequence traces to an RFdiffusion backbone.
- **"Histidine-containing chemistries were retained."** True for most of the portfolio but
  not all of it. Per `tfr1_submission_liability_audit_v1.csv`, **TF4_014_L67_XS05 and
  TF4_017_IL4_246 contain zero histidines**, and TF4_001, TF4_012, TF4_015 and TF4_016
  contain only one. Histidine content was tracked and diversified rather than uniformly
  retained.
- **"Structurally aligned human and mouse TfR1."** The resulting residue mapping and
  conserved apical epitope masks are present and checked off in `README.md`, but the
  alignment run itself is not reproducible from a single committed artefact.

### Material context the statement does not mention

- **No design passed the project's own Stage 4 cross-species pose gate.** The internal gate
  required human protein iPTM >= 0.50 **and** mouse protein iPTM >= 0.50 **and** retained
  apical/preferred contacts in both species **and** homologous contact overlap >= 3 **and**
  homologous Jaccard >= 0.15. **All 25 `pose_pass_ids.json` files under `predictions/` are
  empty** (`{"pose_pass_ids": []}`), covering every FreeBindCraft, Odin-Multi, inference-loop
  and NVIDIA campaign in the repository. Not one design in this project cleared the gate.
- **The closest near-miss.** `TF4_007_V5_T1` has the highest cross-species homologous Jaccard
  anywhere in the portfolio at **0.1489**, below the 0.15 threshold, with homologous overlap
  7 (which does clear the >= 3 sub-criterion), human protein iPTM **0.2946** and mouse protein
  iPTM **0.1374** - both far below the 0.50 requirement, and its shared epitope is distal
  rather than apical. `TF4_001_L70` has the largest homologous overlap (10) but a Jaccard of
  only 0.1205. No design clears the gate on any combination of criteria.
- **Highest scores in the portfolio are single-species and partly pose-conditioned.** The
  highest human iPTM is 0.7341 (TF4_012_MS15) with mouse at 0.2509 and homologous overlap 0;
  the highest mouse iPTM is 0.6769 (TF4_001_L70) on a non-apical epitope. The portfolio
  contains no design that is simultaneously confident in both species.
- **The author's statement is honest about not claiming a proven pH switch.** That caution is
  supported: `TFR1_METHODS_AND_PROVENANCE.md` records that the planned PROPKA / pKa-ranking
  stage was gated behind a neutral-pH cross-species parent criterion that no candidate
  cleared, so the pH-switch stage was never reached. No pH 6.0 negative-state evaluation was
  implemented at all.

### Could not be verified either way

- Whether the structural alignment specifically "identified" the conserved apical region, as
  opposed to the region being taken from the challenge's recommended apical epitope. The
  README lists the apical domain as the competition's own recommended functional epitope.
- The provenance of inputs outside this repository (model weights, versions, exact
  hyperparameters for every run). Campaign scripts are committed, but run environments are not
  fully pinned; `requirements.txt` contains only `biopython==1.85`.
