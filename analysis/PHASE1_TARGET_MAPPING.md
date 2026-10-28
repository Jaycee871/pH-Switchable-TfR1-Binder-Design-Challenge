# Phase 1 — TfR1 target preparation and cross-species mapping

## Ground-truth target constructs

- **Human TfR1:** P02786-1 residues **121–760**, 640 aa, soluble ectodomain homodimer.
- **Mouse TfR1:** Q62351-1 residues **123–763**, 641 aa, soluble ectodomain homodimer.
- **Human structure:** PDB **11QC**, chains A–B.
- **Recommended design region:** apical domain.
- **Reference apical interface:** the 2026 VP1u–TfR1 structure (PDB **11RN**) contacts human TfR1 primarily at residues **208–212**, plus **I202, N348, and K371**.

A BLOSUM62 global alignment of the two assay ectodomains gives **499 identical residues across 639 paired positions = 78.09% identity**, consistent with the competition's stated ~78% ectodomain identity.

## Critical mapping caveat

A naive global alignment is **not sufficient** for the Tyr211 patch because a short insertion/deletion makes several equally plausible sequence alignments possible.

Proteinbase provides explicit challenge-level correspondences for the species-variable hotspot. These are treated as authoritative overrides:

| Human | Mouse | Interpretation |
|---|---|---|
| D204 | Q207 | substitution |
| R208 | N211 | substitution |
| V210 | D213 | substitution |
| Y211 | P214 | substitution |
| L212 | — | no counterpart |
| N348 | K350 | substitution |
| K371 | Q373 | substitution |

The reproducible mapping script keeps the raw global-alignment mapping and applies these challenge-specific overrides before downstream use.

## Design consequence

The exact VP1u footprint is **not a good cross-species anchor by itself**: most of its highlighted contact chemistry sits in the species-variable patch. The first design campaign should therefore use a broader apical surface and split generation into two complementary masks:

1. **Conserved-apical mask** — prioritise conserved residues adjacent to the known functional patch and penalise dependence on the explicitly variable positions.
2. **Switch-engineering mask** — permit peripheral contact with the functional patch while introducing pH-responsive chemistry primarily on the binder side.

This preserves the biological advantage of the apical domain while reducing the risk of generating a human-only binder.

## Allele note

PDB 11QC and the challenge human construct correspond to P02786 with **Gly142**. Older TfR1 structures can contain the alternative **Ser142** allele and should not be mixed into this workflow without correction.
