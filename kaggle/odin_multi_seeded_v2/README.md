# Odin-Multi seeded v2 — TfR1 cross-reactivity

This campaign follows the failed one-trajectory random-initialization smoke.

## Why v2 changes the initialization

The smoke trajectory never entered a joint binding basin: its maximum
worst-target iPTM was only ~0.157 and the final human/mouse interface losses
remained high. In contrast, the FreeBindCraft v6 parent has an independently
validated human full-complex pose with protein iPTM 0.690 and retained
conserved contacts 347/350.

v2 therefore warm-starts Odin-Multi from the **v6 60-aa sequence** while still
allowing the shared sequence to be optimized jointly against human and mouse.

## Objective

- Human on-target: hotspots **347 / 350**, gradient weight 1.0
- Mouse on-target: homologous hotspots **349 / 352**, gradient weight 1.25
- Shared 60-aa sequence
- 3 independent trajectories
- Interface contact cutoff 8 A
- interface contact weight 1.5
- interface PAE weight 0.15
- iPTM weight 0.15
- ProteinMPNN autoregressive prior reduced to 0.5
- 60 logits / 30 soft / 15 STE iterations
- selection: `best_i_ptm`

All selected trajectories must be independently re-docked against the **full
human and full mouse ectodomains** with sequence-only Boltz-2. Odin's own
optimized pose is not accepted as validation evidence.
