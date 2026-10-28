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

## Current findings (2026-10-28)

Full detail in [`analysis/PHASE_STATUS_2026-10-28.md`](analysis/PHASE_STATUS_2026-10-28.md).
All computational results have now been recovered from GitHub Actions artifacts and committed;
see [`analysis/ARTIFACT_PROVENANCE.md`](analysis/ARTIFACT_PROVENANCE.md).

**There is currently no submittable candidate.** Two conclusions drive the remaining work:

1. **Pose collapse is universal.** For every structurally validated design, the predicted human
   interface lands at approximately residues 507-760 rather than the intended apical patch
   (202-215 / 343-374). Interface pTM is 0.14-0.30 and the intended apical anchors are not
   retained in any case, even when the RFdiffusion geometry is supplied to Boltz-2 as an
   explicit mmCIF template. High overall confidence and complex pLDDT in these runs reflect
   confident folding of the individual chains, not a correct docking pose.

2. **No pH-switch mechanism exists yet.** Not one of the 60+ designed sequences contains a
   histidine. The Phase 3D interface bias left histidine *available* but never *rewarded* it, so
   ProteinMPNN never selected it. Since histidine protonation is the proposed switch mechanism,
   the designs currently have no physical basis for pH selectivity.

A third point deserves emphasis because the raw numbers are misleading: `TFR1_C57_02` and
`TFR1_C57_05` show the project's highest mouse interface pTM (0.425 / 0.433), but their predicted
mouse interface (residues 148-200) is a **different surface** from their own predicted human
interface (507-651). That is not cross-reactivity; it is two unrelated docking artifacts.

### Still missing

- a validated binding pose for any candidate;
- histidine-containing designs;
- a pH 6.0 negative-state evaluation (requires new code — nothing in the pipeline evaluates pH 6.0);
- structural diversity (all designs trace back to just two backbones);
- developability filters (pI, net charge, aggregation, low-complexity);
- the submission CSV.

The only backbone to have passed its own gate is the compact apical **seed 20261102** (57 aa,
preferred contacts 344 and 347, Rg 14.88 Å), which is the planned starting point for the next
design round.

> **Deadline: 2026-11-01.** The submission window closes on November 1, 2026.

## Status

- [x] Repository initialized
- [x] Retrieve and validate human TfR1 target structure (PDB 11QC, 11RN)
- [x] Prepare mouse TfR1 target (UniProt Q62351, residues 123-763)
- [x] Build human-mouse residue mapping (Proteinbase structure-aware correspondence)
- [x] Define conserved apical epitope masks (Tier A / Tier B / RISK)
- [x] Generate first design library (60+ sequences across the A29, BIAS and C57 families)
- [ ] Run multi-state structural evaluation - *partial: human and mouse pH 7.4 done for 6
  candidates; the pH 6.0 negative state is not implemented*
- [ ] Rank and cluster candidates - *blocked: no candidate currently passes the interface pose
  gate, and no diversity clustering has been run*
- [ ] Prepare final Proteinbase submission package - *not started; `submissions/` is empty*

---

**Anthropic × Adaptyv Protein Design Competition 2026 — Challenge 04: TfR1**
