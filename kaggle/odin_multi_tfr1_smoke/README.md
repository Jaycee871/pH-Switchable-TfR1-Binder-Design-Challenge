# Odin-Multi TfR1 cross-reactivity smoke

True multi-context design after the FreeBindCraft v6 backbone and all subsequent
single-state / fixed-backbone sequence rescues failed to preserve the homologous
mouse apical pose.

## Contexts

- Human TfR1 apical domain, residues 189-383
  - hotspots: 343, 344, 347, 350
- Mouse TfR1 homologous apical domain, residues 191-385
  - hotspots: 345, 346, 349, 352

Both contexts use Odin-Multi's **target** role, so their gradients jointly update
one shared binder sequence.

## Smoke budget

- binder: fixed 60 aa
- 1 trajectory
- 30 soft / 20 temporary / 10 hard / 5 greedy iterations
- contact cutoff: 8 A
- inter-contact weight: 1.25
- one required contact per hotspot
- selection: `best_i_ptm` (maximize the weakest on-target iPTM)

This run tests whether true multi-context optimisation can solve the cross-species
pose failure that single-target FreeBindCraft and post-hoc ProteinMPNN could not.

Odin-Multi selection is not the final gate. Selected sequences must still pass
the repository's independent full human/mouse Boltz-2 validation before any
pH-switch engineering.
