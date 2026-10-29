# FreeBindCraft / ColabDesign hotspot-objective audit

Date: 2026-10-29.

This audit was performed after the FreeBindCraft v3 harvest showed a discrepancy between
the requested TfR1 hotspot and the final 5 A trajectory interface.

## 1. Hotspot syntax was valid

FreeBindCraft passes `target_hotspot_residues` directly to ColabDesign:

`af_model.prep_inputs(... hotspot=target_hotspot_residues ...)`

ColabDesign `prep_pos()` resolves a numeric hotspot without a chain prefix against the
first target chain. Therefore, for target chain A:

`343,344,347,350`

is valid and resolves to the original PDB residue numbers 343, 344, 347 and 350.
Adding explicit `A` prefixes would be clearer but would not fix the v3 behavior by itself.

## 2. The default objective is soft and spatially permissive

FreeBindCraft rewrites ColabDesign's inter-chain contact options to:

- `num = inter_contact_number` (default 2)
- `cutoff = inter_contact_distance` (default 20 A)
- `binary = False`

When a hotspot mask is present, ColabDesign computes `i_con` only between hotspot target
positions and binder positions. `get_con_loss()` takes the best `num` binder contacts per
hotspot position, then averages over the hotspot target positions (`num_pos` defaults to
infinity).

Thus the default objective is *not* equivalent to requiring a physical 5 A interface on
each hotspot residue. A trajectory can obtain a favourable hotspot-biased loss while the
final heavy-atom 5 A interface is elsewhere on the target.

This explains why v3 trajectories could show `Trajectory_WrongHotspot = 0` while independent
geometry showed no 343/344/347/350 contact.

## 3. The requested hotspot patch is geometrically feasible

11QC chain A C-alpha distances:

| Pair | CA distance |
|---|---:|
| 343-344 | 3.82 A |
| 343-347 | 6.07 A |
| 344-347 | 5.02 A |
| 347-350 | 8.35 A |
| 343-350 | 12.55 A |
| 344-350 | 12.67 A |

The four residues therefore form one compact surface neighborhood that a 60-aa binder can
plausibly cover. Failure to engage the patch should not be attributed to an impossible
geometric target.

## 4. v4 change

The v4 single-attempt experiment keeps the same valid hotspot specification and changes the
objective instead:

- fixed binder length: 60 aa
- one attempt only
- hotspot residues: 343/344/347/350
- `inter_contact_distance = 8.0 A`
- `weights_con_inter = 2.0`
- `inter_contact_number = 2`
- independent final heavy-atom 5 A hotspot gate

The independent gate requires:

- >=1 contact with 343/344/347/350
- >=2 preferred apical contacts
- <=1 explicit species-risk contact
- <=10 severe atom clashes (<2 A)

Only a trajectory passing this independent geometry gate becomes eligible for full
human/mouse Boltz-2 validation.

## 5. If v4 still fails

The next escalation should change the loss form, not simply add random seeds. ColabDesign
supports a binary contact loss, but FreeBindCraft currently hardcodes `binary=False` for
`i_con`. A future experiment could patch that behavior so the objective directly maximizes
probability of distances below a tight cutoff.

That experiment should be justified only if v4 shows that an 8 A categorical contact
objective remains insufficient.
