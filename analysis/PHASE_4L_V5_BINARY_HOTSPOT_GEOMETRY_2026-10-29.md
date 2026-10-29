# Phase 4L — FreeBindCraft v5 binary-hotspot geometry analysis

Date: 2026-10-29.

## Why v5 was run

FreeBindCraft v4 tightened the hotspot-contact cutoff from the default 20 A to 8 A,
but retained ColabDesign's categorical contact loss. Its single controlled trajectory
did not survive FreeBindCraft's clash gate.

v5 changed exactly one design-objective feature:

`i_con binary=False -> binary=True`

The target, hotspot set (343/344/347/350), 60-aa binder length, single-attempt budget,
8 A cutoff, contact weight 2.0 and two contacts per hotspot were retained.

## Workflow caveat

The first GitHub monitor used the wrong Kaggle slug
(`tfr1-freebindcraft-smoke-v5-binary`). The actual notebook created by Kaggle was:

`packkwanlow/tfr1-freebindcraft-v5-binary-hotspot`

This caused the GitHub monitor to fail even though the Kaggle notebook itself completed.
The raw output was subsequently re-harvested using the correct slug. No model rerun was
required.

## FreeBindCraft result

The v5 notebook completed normally (`bindcraft_exit_code = 0`) and generated one
trajectory:

`TFR1_FBC343_V5_l60_s289817`

FreeBindCraft classified it as `Clashing`, so no trajectory row or MPNN design was
carried forward.

The failure table shows:

- Trajectory_Clashes = 1
- Trajectory_WrongHotspot = 0
- other trajectory filters = 0

Thus the design was not rejected for missing the requested hotspot.

## Independent heavy-atom interface analysis

The v5 clashing trajectory makes genuine 5 A contacts with:

- preferred residues: **209, 214, 343, 344, 347**
- hotspot core: **343, 344, 347**
- risk residues: **210, 211, 212, 348, 371**

Minimum heavy-atom distances to the four requested hotspots:

| Human TfR1 residue | Minimum distance |
|---|---:|
| 343 | 2.05 A |
| 344 | 2.05 A |
| 347 | 4.06 A |
| 350 | 9.22 A |

This is a large improvement in hotspot localisation compared with v4.

## v4 vs v5 controlled comparison

Both runs used a 60-aa binder and the same deterministic CLI seed.

### v4 categorical 8 A objective

- 5 A hotspot contacts: **350 only**
- preferred contacts: 350
- explicit risk contacts: none
- severe heavy-atom pairs <2 A: **150**
- inter-chain CA clashes <2.5 A: **6**
- binder self-CA clashes <2.5 A: **0**

Minimum hotspot distances:

- 343: 15.81 A
- 344: 15.98 A
- 347: 12.96 A
- 350: 4.83 A

### v5 binary 8 A objective

- 5 A hotspot contacts: **343, 344, 347**
- preferred contacts: 209, 214, 343, 344, 347
- explicit risk contacts: 210, 211, 212, 348, 371
- severe heavy-atom pairs <2 A: **52**
- inter-chain CA clashes <2.5 A: **0**
- binder self-CA clashes <2.5 A: **2**

Binder self-clashes:

- B15 GLY ↔ B53 ALA: 2.336 A
- B15 GLY ↔ B54 GLU: 1.777 A

The hotspot-facing segment is around binder residues 35-39. The self-clashing segment
(B15 vs B53/B54) is spatially distinct from the hotspot-contacting segment.

Direct binder CA displacement between v4 and v5 is large (RMSD ~41 A), confirming that
binary contact loss materially changed the generated pose rather than merely tweaking
a local side-chain arrangement.

## Interpretation

**The binary hotspot objective worked.**

v5 should not be interpreted as another hotspot-targeting failure. It moved the physical
interface onto three of the four requested conserved hotspot residues and eliminated the
inter-chain CA clashes observed in v4.

The remaining failure is primarily binder backbone self-crossing, plus excessive
peripheral engagement of species-risk surfaces.

One likely contributor is that v5 still used the deliberately shortened smoke-test
optimisation budget:

- 25 soft
- 15 temporary
- 3 hard
- 5 greedy iterations

compared with the normal FreeBindCraft 4-stage budget:

- 75 soft
- 45 temporary
- 5 hard
- 15 greedy iterations

v5 also used a stronger-than-default compactness term (Rg weight 0.40) and a strong
binary contact term (weight 2.0, two contacts per hotspot).

## Decision: v6 tempered binary

Do **not** abandon binary hotspot targeting.

v6 is a controlled optimisation-quality follow-up:

- retain binary i_con
- fixed 60 aa
- same hotspot set
- same 8 A cutoff
- one controlled attempt
- hard-target preset / predict_initial_guess
- restore full 75/45/5/15 optimisation
- inter-contact weight: 2.0 -> **1.25**
- contacts per hotspot: 2 -> **1**
- Rg weight: 0.40 -> **0.30**
- iPTM loss weight: 0.10 -> **0.05**

The goal is to preserve the successful 343/344/347 localisation while allowing enough
geometric freedom to eliminate binder self-crossing and reduce species-risk peripheral
contacts.

Full human/mouse validation remains downstream of the independent 5 A hotspot/risk/clash
gate. No pH engineering is justified before that gate passes.
