# Phase 1 — TfR1 target and epitope definition

## Target state

The competition assay uses soluble TfR1 ectodomain homodimers:

- Human TfR1: UniProt P02786-1, residues 121–760 (640 aa)
- Mouse TfR1: UniProt Q62351, residues 123–763 (641 aa)
- Human structural reference: PDB 11QC
- Apical-domain reference complex: PDB 11RN

## Design problem

The binder must be treated as a three-state design problem:

| State | Desired outcome |
|---|---|
| Human TfR1, pH 7.4 | ON |
| Mouse TfR1, pH 7.4 | ON |
| Human TfR1, pH 6.0 | OFF / strongly weakened |

The competition ranks successful designs by pH selectivity first, mouse cross-reactivity second, and human affinity third.

## Apical-domain cross-species risk

Proteinbase identifies the functional apical patch around human Tyr211 (human residues 202–215, 348, 371). Six of nine key positions differ between human and mouse.

| Human | Mouse | Consequence |
|---|---|---|
| D204 | Q207 | Avoid human-only charge anchoring |
| R208 | N211 | Avoid obligatory Arg salt bridge |
| V210 | D213 | Large hydrophobic→acidic chemistry change |
| Y211 | P214 | Central hotspot is not side-chain conserved |
| L212 | — | No direct mouse counterpart |
| N348 | K350 | Polar→positive change |
| K371 | Q373 | Positive→polar change |

### Important numbering note

A conventional global sequence alignment is locally ambiguous around the indel near human 210–212. For design masks and interpretation, use the **Proteinbase structure-aware correspondence above**, not an automatically inferred residue-number offset.

## Initial design-mask policy

### Cross-species positive mask
Prefer contacts that:
- involve conserved apical-domain backbone geometry;
- distribute binding energy over multiple conserved residues;
- tolerate alternative side-chain chemistry at the seven risk positions above.

### Human-specific penalty mask
Down-rank designs whose predicted interface depends strongly on:
- D204 charge,
- R208 charge,
- V210 hydrophobic packing,
- Y211 aromatic packing,
- L212 geometry,
- N348 side-chain identity,
- K371 charge.

These residues are not forbidden contacts. They are **risk positions**: a design can touch them, but should not make one of them its only energetic anchor.

### pH-switch design layer
The pH-switch objective should be introduced on the **binder side** where possible, so cross-species target recognition is not sacrificed merely to create acid sensitivity. Histidine-centered networks are one candidate mechanism, but will be evaluated as one family rather than assumed to work.

## Next phase

1. Retrieve 11QC and 11RN coordinates.
2. Identify solvent-exposed conserved residues surrounding the apical patch.
3. Build explicit positive and penalty residue masks.
4. Generate separate conserved-apical and pH-switch design families.
5. Score every sequence against the three states before ranking.
