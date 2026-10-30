# Odin-Multi anti-escape v3

Seeded V2 solved the **cropped-apical** human/mouse objective (worst-target iPTM
0.738-0.766) but all three selected sequences escaped to non-apical surfaces when
re-docked against full ectodomains.

V3 therefore introduces explicit negative design against the recurrent escape basins.

## On-target contexts

- Human apical TfR1: hotspots 347 / 350, gradient 1.0
- Mouse apical TfR1: homologues 349 / 352, gradient 1.25

## Off-target escape contexts

Each off-target PDB contains two spatially preserved target fragments as separate chains:

- Human: A=122-188, B=500-760
- Mouse: A=123-188, B=500-763

These cover the N-terminal and distal surfaces repeatedly selected by independent
full-target Boltz-2 predictions.

Each off-target context contributes a modest repulsive gradient (0.20), negative
interface-contact / iPTM terms and specificity clipping. The desired apical targets
still dominate the total gradient.

## Search

- FreeBindCraft v6 60-aa sequence warm-start
- 1 diagnostic trajectory
- 60 logits / 30 soft / 15 STE / 5 greedy
- 8 A target and off-target contact distance
- selection: best_clipped_i_pae_ratio

Phase 4S independently re-docks the selected sequence from sequence only against
the full human and mouse ectodomains. Odin's own predicted pose is not validation.

If the anti-escape objective improves full-target transfer, the same Odin run can be extended to additional trajectory indices without changing the immutable contexts or restarting completed work.
