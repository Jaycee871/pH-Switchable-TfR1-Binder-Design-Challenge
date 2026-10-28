# Structure-derived apical design mask

This mask was derived from the Phase 1 target assets generated in GitHub Actions.

## Reference-interface check

Using PDB 11RN, TfR1 residues within 5 Å of VP1u were measured independently for both symmetry-related interfaces.

Both interfaces returned the same ten human TfR1 contact residues:

`202, 203, 204, 208, 209, 210, 211, 212, 348, 371`

This independently recovers the functional patch described on the Proteinbase challenge page.

## Surface-neighborhood analysis

The structural neighborhood was expanded around those reference-interface residues in PDB 11QC. Residue solvent accessibility was estimated on the homodimer, then combined with human–mouse sequence conservation.

The resulting design mask uses three classes:

- **Tier A — preferred anchor:** conserved human/mouse identity and strong surface accessibility near the apical patch.
- **Tier B — supporting contact:** conserved and structurally nearby, but less exposed.
- **RISK — cross-species-sensitive:** official Proteinbase hotspot positions whose mouse correspondence changes identity or geometry.

The mask is intentionally not a hard prohibition list. A generated binder may contact a RISK residue, but candidates that concentrate most of their predicted interface energy on one or more RISK residues should be down-ranked.

## Preferred anchor pool

Primary cross-species anchors:

`N206, G207, L209, E214, K231, E343, K344, G347, E350, K374`

Supporting conserved contacts:

`P216, V228, A340, G351, N372`

## Strong cross-species risk positions

`D204, R208, V210, Y211, L212, N348, K371`

For the local indel around human 210–212, the Proteinbase structure-aware mapping is authoritative.

## Design implication

The first generation campaign should not ask a model merely to "bind around Y211". Instead, it should encourage an interface spanning several **conserved anchor residues** while treating the human-specific hotspot chemistry as optional/peripheral.

This is the basis for Family A (conserved-apical binders) and the positive-state component of Family B (multi-state binders).
