# Cross-species liability map — FreeBindCraft l67 trajectory

Candidate: `TFR1_FBC343_l67_s385919`

This note explains why the human full-complex pose can retain the intended apical
epitope while the mouse prediction collapses.

## Conserved anchors

- Human K231 -> mouse K233
- Human E350 -> mouse E352
- Human G351 -> mouse G353

Binder residue **R10** is the cleanest conserved anchor:
- R10 -> human E350: 3.30 A
- R10 -> human G351: 3.77 A

Binder residue **Y21** also contacts only conserved human E350 in the local geometry.

## Main cross-species liabilities

The human trajectory additionally relies on several chemically non-conserved peripheral
contacts:

- human 205 K -> mouse 208 S
- 247 Y -> 249 S
- 248 T -> 250 Y
- 249 P -> 251 S
- 352 D -> 354 S
- 365 M -> 367 L
- 366 V -> 368 E
- 369 E -> 371 Q

The binder positions most exposed to these mismatches are:

- **W7**: contacts human 248, 352, 365, 366 in addition to conserved 350/351.
- **Y8**: strongly dependent on the 247/365 neighborhood.
- **Y11**: contacts conserved K231 but also human-specific 247/248/249.
- H14 and K18 contact human 369 E -> mouse 371 Q.

By contrast, **R10** is dominated by the conserved 350/351 pair.

## Rescue implication

If a sequence-level rescue is attempted after the hard-hotspot v4 experiment:

1. keep R10 fixed;
2. preserve the E350/G351 anchor geometry;
3. redesign W7/Y8/Y11 toward chemically tolerant residues rather than bulky aromatics;
4. treat H14/K18 as secondary liabilities because the 369 E -> Q substitution changes
   local electrostatics;
5. revalidate human and mouse full ectodomains before any pH engineering.

This is a fallback path only. The v4 hard-hotspot generator test has priority.
