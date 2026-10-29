# Phase 4M — FreeBindCraft v6 tempered-binary result

Date: 2026-10-29.

## Controlled objective

v6 was designed to test whether the successful v5 binary hotspot objective could be
retained while eliminating the v5 binder self-crossing.

Changes relative to v5:

- binary inter-chain contact loss retained;
- full 4-stage optimisation restored (75 / 45 / 5 / 15);
- hard-target preset / initial guess enabled;
- contact weight reduced from 2.0 to 1.25;
- required contacts per hotspot reduced from 2 to 1;
- radius-of-gyration weight reduced from 0.40 to 0.30;
- iPTM loss returned to 0.05;
- fixed 60-aa binder and one controlled attempt retained.

## Result

FreeBindCraft produced one **Relaxed** trajectory:

`TFR1_FBC343_V6_l60_s289817`

Generation-stage metrics:

- pLDDT: **0.75**
- pTM: **0.76**
- i_pTM: **0.35**
- i_pAE: **0.50**
- unrelaxed clashes: 30
- relaxed clashes: **0**
- shape complementarity: **0.61**
- PackStat: **0.65**
- interface residues: **18**
- trajectory time: approximately **60 minutes**

The v5 binder self-crossing failure was therefore resolved.

## Independent 5 A geometry

Target contacts:

`210, 211, 231, 243-250, 274, 346-352, 354, 357, 364-366, 368-369, 371`

Preferred contacts:

- **231**
- **347**
- **350**

Requested hotspot-core contacts:

- **347**
- **350**

Species-risk contacts:

- **210**
- **211**
- **348**
- **371**

Severe heavy-atom clashes below 2 A:

- **0**

The independent strict local gate therefore failed only because the trajectory makes
too many explicit species-risk contacts. It was not sent to full human/mouse Boltz-2.

## Residue-level liability mapping

The main risk liability is highly concentrated:

### Binder B46 = W

B46 contacts:

- human 210
- human 211
- human 347
- human 348
- human 371

Thus one large aromatic side chain spans both the desired conserved hotspot and several
species-variable surfaces.

### Binder B47 = E

B47 contacts:

- human 371

By contrast, hotspot 350 has multiple independent supports:

- B23
- B42
- B43
- B44
- B52
- B54

and hotspot 347 is also supported by B44/B45.

## Interpretation

v6 demonstrates a second independent improvement over the earlier generator campaigns:

1. binary contact loss can place the interface on the desired 343-350 region;
2. tempered optimisation can remove the v5 backbone self-crossing;
3. the remaining cross-species liability is no longer distributed throughout the binder;
   it is concentrated primarily in B46/B47.

This supports a narrow sequence-level cleanup rather than another generator restart.

## Decision: strict two-position cleanup

Phase 3J fixes every binder position except B46/B47 using ProteinMPNN
`fixed_positions_jsonl`.

Parent:

- B46 = W
- B47 = E

Preferred cleanup chemistry:

- B46: small/neutral `A/G/N/Q/S/T/V`
- B47: small/neutral `A/G/N/Q/S/T`

B44/B45 and all other hotspot-supporting positions remain fixed.

No pH engineering is performed. Any cleanup variant must pass the unchanged full
human/mouse neutral-pH pose gate before histidine switching is considered.
