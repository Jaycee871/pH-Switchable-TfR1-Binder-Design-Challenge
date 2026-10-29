# FreeBindCraft v6 — tempered binary hotspot

Controlled follow-up to v5.

## What v5 proved

With the same 60-aa design and seed family, the binary hotspot loss moved the
actual interface onto human TfR1 343/344/347. The v5 clashing trajectory had:

- hotspot-core contacts: 343, 344, 347
- no inter-chain CA clash below 2.5 A
- two binder self-CA clashes (B15 with B53/B54)
- strong contact with species-risk residues around 210-212, 348 and 371

Thus binary hotspot targeting worked, but the reduced smoke optimisation produced
an over-compressed / self-crossing binder geometry.

## v6 controlled changes

Retained:
- binary i_con = true
- target 11QC chain A residues 189-383
- hotspots 343/344/347/350
- fixed 60-aa binder
- 8 A contact cutoff
- one attempt
- same deterministic CLI seed

Changed:
- hard-target preset (predict_initial_guess=true)
- full 4-stage optimisation: 75 / 45 / 5 / 15 iterations
- inter-contact weight: 2.0 -> **1.25**
- contacts per hotspot: 2 -> **1**
- radius-of-gyration weight: 0.40 -> **0.30**
- iPTM loss weight: 0.10 -> **0.05**
- 8 ProteinMPNN sequences / max 2 downstream sequences

The goal is not weaker hotspot targeting. It is to retain the successful v5
343-350 localisation while giving the binder enough optimisation freedom to
avoid self-crossing and reduce species-risk peripheral contacts.

Independent 5 A hotspot/risk/clash gating remains mandatory before full human/mouse
Boltz-2 validation.
