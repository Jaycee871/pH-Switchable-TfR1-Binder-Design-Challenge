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

## Status

- [x] Repository initialized
- [ ] Retrieve and validate human TfR1 target structure
- [ ] Prepare mouse TfR1 target
- [ ] Build human–mouse residue mapping
- [ ] Define conserved apical epitope masks
- [ ] Generate first design library
- [ ] Run multi-state structural evaluation
- [ ] Rank and cluster candidates
- [ ] Prepare final Proteinbase submission package

---

**Anthropic × Adaptyv Protein Design Competition 2026 — Challenge 04: TfR1**
