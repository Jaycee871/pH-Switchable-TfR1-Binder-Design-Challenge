# pH-Switchable TfR1 Binder Design Challenge

**Cross-species, pH-switchable TfR1 binder design for the Anthropic × Adaptyv Protein Design Competition 2026.**

> **Design objective:** Human TfR1 pH 7.4 **ON** · Mouse TfR1 pH 7.4 **ON** · Human TfR1 pH 6.0 **OFF**

## Challenge

This repository documents our computational workflow for **Challenge 04: TfR1 (CD71)** in the Anthropic × Adaptyv Protein Design Competition 2026.

- **Target:** Transferrin receptor 1 (TfR1 / CD71)
- **Human target construct:** soluble ectodomain homodimer, residues 121–760
- **Structure reference:** PDB **11QC**, chains A–B
- **Recommended functional epitope:** apical domain
- **Submission window:** October 19 – November 1, 2026
- **Challenge page:** https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1

## Experimental objectives

The challenge evaluates the same design under multiple conditions:

1. **Human TfR1 binding at pH 7.4**
   - measurable binding is required;
   - stronger affinity is preferred.

2. **Mouse TfR1 cross-reactivity at pH 7.4**
   - designs should bind both human and mouse TfR1;
   - among human binders with KD < 1 μM, the mouse-to-human KD ratio should be as close to 1 as possible.

3. **pH-selective release**
   - bind human TfR1 at pH 7.4;
   - ideally show **no detectable binding at pH 6.0**.

### Ranking priority

**pH-selective binding → mouse cross-reactivity → human affinity**

This makes the project a genuine **multi-state design problem**, rather than a conventional affinity-only binder design task.

## Design strategy

We will explore complementary de novo design families rather than optimize a single scaffold family:

### Family A — Conserved-apical binders
Target structurally conserved regions of the human/mouse apical-domain surface to prioritize cross-species recognition.

### Family B — Multi-state binders
Optimize one sequence against three design states:

- Human TfR1 at pH 7.4 — positive state
- Mouse TfR1 at pH 7.4 — positive state
- Human TfR1 at pH 6.0 — negative/release state

### Family C — Explicit pH-switch designs
Explore interface chemistries, including histidine-mediated interaction networks, intended to destabilize the complex as the environment becomes acidic.

## Planned computational workflow

```text
Target preparation
      ↓
Human ↔ Mouse sequence/structure mapping
      ↓
Conserved epitope analysis
      ↓
De novo backbone / binder generation
      ↓
Sequence design
      ↓
Human pH 7.4 complex prediction
      ↓
Mouse pH 7.4 cross-reactivity prediction
      ↓
pH-switch / interface analysis
      ↓
Interface-confidence and developability filters
      ↓
Structural clustering + diversity selection
      ↓
Final ranked submission set
```

Candidate tools may include **RFdiffusion**, **ProteinMPNN**, **Boltz/BoltzGen**, structure-based interface analysis, and complementary physics- or sequence-based filters. Exact models, versions, parameters, and provenance will be recorded as the workflow is executed.

## Repository structure

```text
.
├── README.md
├── data/
│   ├── human_TfR1/
│   └── mouse_TfR1/
├── structures/
├── notebooks/
├── scripts/
├── designs/
│   ├── family_A_conserved/
│   ├── family_B_multistate/
│   └── family_C_pH_switch/
├── predictions/
├── analysis/
├── metadata/
└── submissions/
```

## Submission constraints

According to the competition rules:

- **de novo designs only**
- **10–250 amino acids**
- accepted formats include single-chain proteins, nanobodies, scFvs, and Fabs
- submission CSV requires at minimum:
  - `name`
  - `sequence`
  - `molecule_class`

Additional computational metrics, predicted structures, methodology, and provenance will be retained to support design selection and reproducibility.

## Current findings (2026-10-30)

The project has moved beyond the 2026-10-28 fixed-backbone stage. The central result is now
clearer: **human apical binding is achievable, but no design has yet demonstrated a homologous
mouse apical pose under independent full-ectodomain prediction.**

### 1. FreeBindCraft v6 established a real human-side signal

The tempered-binary v6 parent `TFR1_FBC343_V6_l60_s289817` is the strongest mechanistically
credible human parent produced so far. In the independent Phase 4Q full-complex gate-override
diagnostic it achieved:

- human protein iPTM **0.690**;
- preferred contacts **231, 347, 350**;
- conserved core **347/350 retained**;
- only one reproduced species-risk contact (**371**).

This also showed that the earlier local geometry gate was partly over-conservative: three
trajectory-level risk contacts (210/211/348) did not reproduce in the independent human
full-complex prediction.

However, the same v6 parent failed mouse cross-reactivity decisively:

- mouse protein iPTM **0.282**;
- predicted mouse interface moved outside the homologous apical region;
- homologous human/mouse contact overlap **0**.

See `analysis/PHASE_4Q_V6_PARENT_GATE_OVERRIDE_2026-10-29.md`.

### 2. Fixed-backbone sequence rescue is exhausted

Two sequence-level rescue strategies were tested on the v6 backbone.

**Phase 3J / 4N — strict B46/B47 cleanup**

Only B46/B47 were redesigned while the remaining 58 residues were fixed. Human binding
remained strong (protein iPTM approximately 0.70), and the best mouse score reached **0.492**,
but the mouse prediction still left the homologous apical epitope and homologous overlap
remained zero.

**Phase 3K / 4O — tied human+mouse dual-state ProteinMPNN**

The same sequence was constrained against explicit human and mouse structural states with
all 60 binder positions tied across states. The best human candidate reached protein iPTM
**0.734**, while mouse remained only **0.16-0.25** with homologous overlap **0**.

Therefore the remaining bottleneck is not a small set of side chains. The v6 backbone/docking
geometry itself is not a stable cross-species solution. Additional incremental ProteinMPNN
variants of this backbone are closed.

See `analysis/PHASE_4O_FIXED_BACKBONE_RESCUE_EXHAUSTED_2026-10-29.md`.

### 3. Odin-Multi smoke proved the true multi-context pipeline, but not a candidate

A true multi-context Odin-Multi design was then run with human and mouse TfR1 as simultaneous
on-targets sharing one binder sequence. The complete GPU run succeeded, but its single random
trajectory terminated **LowConfidence**.

Odin selection:

- selected worst-target iPTM: **0.143**;
- the entire trajectory never entered a strong joint-binding basin;
- maximum worst-target iPTM across saved frames was only approximately **0.157**.

Independent sequence-only full-ectodomain Boltz-2 validation then failed both species:

- human protein iPTM **0.255**, non-apical;
- mouse protein iPTM **0.136**, non-apical;
- homologous overlap count **2**, Jaccard **0.0488**;
- `cross_species_pose_pass = false`.

This is interpreted as a failed initialization / search trajectory, not evidence that
multi-context optimization is intrinsically unsuitable.

### 4. Current active campaign: Odin-Multi seeded v2

The active v2 campaign uses the independently validated FreeBindCraft v6 sequence as the
shared Odin-Multi warm-start instead of random initialization.

Controlled changes:

- human conserved core hotspots: **347 / 350**;
- mouse homologues: **349 / 352**;
- human gradient weight **1.0**;
- mouse gradient weight **1.25**;
- interface-contact weight **1.5**;
- interface-PAE weight **0.15**;
- iPTM weight **0.15**;
- ProteinMPNN autoregressive prior reduced to **0.5**;
- **3 independent trajectories**;
- 60 logits / 30 soft / 15 STE iterations.

The campaign selects each trajectory by `best_i_ptm`, maximizing its weakest on-target
interface pTM. Phase 4R is already connected and will independently re-dock **all three**
selected sequences against the full human and full mouse ectodomains using sequence-only
Boltz-2. Odin's own optimized pose is not accepted as validation evidence.

### Current scientific gate

A candidate is not considered a neutral-pH cross-species parent unless independent validation
shows, for both human and mouse:

- protein iPTM >= 0.50;
- retained apical / preferred contacts;
- a homologous human-mouse epitope overlap (not merely good scores on different surfaces).

**pH-switch engineering remains downstream of this gate.** Histidine / protonation work will
not be promoted to a final candidate until a neutral-pH cross-species parent exists.

## Status

- [x] Repository initialized
- [x] Retrieve and validate human TfR1 target structure (PDB 11QC, 11RN)
- [x] Prepare mouse TfR1 target (UniProt Q62351, residues 123-763)
- [x] Build human-mouse residue mapping (Proteinbase structure-aware correspondence)
- [x] Define conserved apical epitope masks (Tier A / Tier B / RISK)
- [x] Generate first design library (60+ sequences across the A29, BIAS and C57 families)
- [x] Introduce histidine at the interface (Phase 6: 24 designs, 18 passing, 4-10 His each)
      - *caveat: the +4.0 bias saturated rather than selected, and the designs are unvalidated*
- [ ] Run multi-state structural evaluation - *partial: human pH 7.4 done for all 16 final-v1
  designs and mouse pH 7.4 for 15 of them; the pH 6.0 negative state is not implemented*
- [x] Rank and cluster candidates - *ranked and clustered on multiple signals (human/mouse
  protein iPTM, preferred-contact retention, homologous overlap and Jaccard, sequence liability
  and redundancy) and a 16-design final-v1 portfolio was selected and submitted; no design
  passed the Stage 4 cross-species pose gate*
- [x] Prepare final Proteinbase submission package - *final-v1 is 16 designs; see
  `submission/TFR1_FINAL_V1_FREEZE_NOTE.md`*
- [x] Submit to Proteinbase Challenge 04 - *submitted 2026-10-31; "Submission successful"
  observed on screen at 20:39 Asia/Taipei. No Proteinbase submission ID was captured. See
  `submission/TFR1_SUBMISSION_RECORD.md`*

---

**Anthropic × Adaptyv Protein Design Competition 2026 — Challenge 04: TfR1**
