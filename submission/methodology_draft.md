# Submission methodology draft — TfR1 conditional binder

## Design objective

Generate de novo protein minibinders that satisfy the three experimental objectives of the
2026 Anthropic × Adaptyv TfR1 challenge:

1. bind human TfR1 at pH 7.4;
2. bind mouse TfR1 at pH 7.4 with a comparable epitope/pose;
3. release human TfR1 under acidic endosomal conditions.

The primary design region is the TfR1 apical domain. The workflow deliberately avoids
treating any single in-silico metric as a sufficient selection criterion.

## Stage 1 — structure-aware cross-species epitope definition

Human TfR1 PDB 11QC and human/mouse ectodomain sequences were aligned and mapped to
UniProt numbering. The apical surface was partitioned into:

- conserved preferred anchors, including human 343/344/347/350;
- supporting conserved contacts;
- species-risk residues, including human 348 and 371;
- the highly variable Tyr211-centred patch.

The preferred 343–350 cluster was selected because its anchor residues remain chemically
conserved across human and mouse, while the immediately adjacent N348/K350-equivalent
position is treated as a risk rather than an anchor.

## Stage 2 — generator exploration and negative-evidence gating

RFdiffusion was used as an exploratory backbone generator across progressively restricted
target geometries:

- full / apical target;
- conserved-cluster hotspot campaigns;
- compact 45–60 aa campaigns;
- local 330–383 epitope window.

Each generated backbone was independently analysed for:

- binder length and compactness;
- preferred target contacts;
- species-risk contacts;
- severe clashes;
- interface localisation.

Repeated campaigns demonstrated that compactness, zero-risk contacts and dual-anchor
engagement could each be achieved, but were not reliably obtained together. ProteinMPNN
sequence rescue of the strongest scaffold also failed to preserve the homologous mouse
apical pose during full-complex validation.

Rather than increasing RFdiffusion seed count after this repeated failure mode, the
workflow switched generator class.

## Stage 3 — pose-optimised de novo generation

FreeBindCraft is used as the next primary generator because target–binder geometry is
optimised within an AlphaFold2 hallucination loop rather than being checked only after
backbone generation.

Initial configuration:

- target: human TfR1 11QC chain A, residues 189–383;
- target hotspots: 343, 344, 347, 350;
- binder class: de novo protein minibinder;
- no starting binder sequence or known binder scaffold;
- ProteinMPNN refinement;
- AF2 revalidation;
- open-source relaxation/scoring path;
- ipSAE used as a generation-stage ranking signal.

FreeBindCraft acceptance is explicitly treated as a **generation-stage filter**, not as
evidence of cross-species binding or pH switching.

## Stage 4 — strict full human/mouse pose gate

Top generated sequences are predicted independently against the full human and mouse TfR1
ectodomains with Boltz-2.

A sequence proceeds only when:

- human protein ipTM >= 0.50;
- mouse protein ipTM >= 0.50;
- the human complex retains a preferred apical contact;
- the mouse complex retains the homologous apical surface;
- human/mouse contact sets show explicit homologous epitope overlap.

This gate was introduced because earlier sequence-designed RFdiffusion candidates could
show plausible human apical binding yet relocate to a distal surface on mouse TfR1.

## Stage 5 — mechanistic pH-switch engineering

pH-switch engineering begins **only after** a parent sequence passes the cross-species
pose gate.

Interface residues are ranked for histidine substitution using the predicted human
complex. The scan favours positions where protonated histidine could destabilise the
interface at acidic pH, especially contacts near conserved basic TfR1 residues, while
penalising placements next to acidic target residues that could instead stabilise a
protonated histidine.

The search includes:

- single-His variants for mechanistic attribution;
- constrained double-His variants;
- constrained triple-His clusters, motivated by prior experimental TfR studies in which
  multiple interface histidines produced stronger pH-dependent dissociation than isolated
  substitutions.

These are hypotheses only until explicit protonation / pKa / pH-dependent energetic
evaluation is completed.

## Stage 6 — final ranking

Final submission ranking will integrate, rather than collapse into a single score:

- human neutral-pH pose confidence;
- mouse neutral-pH pose confidence;
- homologous epitope retention;
- interface quality / ipSAE;
- sequence liabilities and structural compactness;
- explicit pH-switch evidence;
- sequence and structural diversity;
- provenance and generator diversity.

The workflow is intentionally failure-aware: negative validation results are used to
change design strategy rather than being hidden or overridden by a favourable metric.

## Provenance

All submitted proteins are de novo designs. Known TfR binders and published pH-sensitive
TfR studies are used only to define biological constraints, interpret mechanisms and
calibrate evaluation logic; no known binder sequence is used as a starting scaffold.
