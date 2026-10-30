# Odin-Multi full-context v4

## Why this campaign exists

Seeded V2 solved the cropped human/mouse apical objective but failed independent
full-ectodomain transfer. Anti-Escape V3 then added cropped off-target decoys;
the decoys looked strongly rejected inside Odin, yet independent Boltz-2 still
placed the binder on distal full-target surfaces. Phase 4T also ruled out the
possibility that V3 merely selected the wrong late trajectory frame.

Therefore the remaining hypothesis is **context truncation itself**.

V4 removes that abstraction. Both Odin target contexts are the complete
ectodomains:

- human 11QC chain A: residues 122-760, hotspots 347/350
- mouse full predicted chain A: residues 123-763, hotspots 349/352

No cropped positive context and no cropped decoy context is used.

## Runtime-bounded diagnostic

Full contexts are much larger than the V3 four-context crop system. V3 required
~245 min of design time for 75 saved frames. A naive 75-frame full-context run
is expected to exceed Kaggle's runtime budget.

V4 is therefore intentionally bounded:

- v6 60-aa sequence warm-start
- one trajectory
- 20 logits / 10 soft / 5 STE / 0 greedy iterations
- one design recycle
- human gradient 1.0
- mouse gradient 1.25
- interface contact weight 2.0
- iPTM weight 0.20
- interface PAE weight 0.20
- selection by best_i_ptm
- Odin pinned to commit `e4567e487d43a4de36524f491e51ad1ea9730d1e`

This is a diagnostic of whether **seeing the whole ectodomain during design**
changes the basin. It is not yet a production campaign.

The selected sequence must still pass independent sequence-only full human/full
mouse Boltz-2 validation (Phase 4U).
