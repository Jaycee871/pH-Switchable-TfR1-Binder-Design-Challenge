# Phase 6 — histidine-enriched design, results and assessment (2026-10-28)

**Run:** GitHub Actions 40691292024 — success, 43 s
**Artifact:** tfr1-his-enriched-seed20261102-proteinmpnn (committed under
`designs/family_C_pH_switch/his_enriched_seed20261102/`)
**Backbone:** RFdiffusion compact apical seed 20261102, the only backbone in this project to
have passed its own scientific gate.

## Outcome in one line

The histidine bias worked, and **overshot**. Every design now carries histidines at the
interface, but the bias saturated those positions rather than selecting among them, and the
dominant feature produced is an N-terminal poly-histidine stretch that resembles a purification
tag more than a designed pH-dependent network.

## Filter results

24 sequences designed (12 per target × 2 sampling temperatures). **18 of 24 pass.**

| Metric | Value |
|---|---|
| Designs with zero histidine | **0** (previously 60+ of 60+) |
| Histidine count | min 4, max 10, mean 8.1 |
| Designs clearing His ≥ 3 | **24 of 24** |
| Designs passing the full filter | **18 of 24** |
| Rejections | 6, all solely for `homopolymer_run_6` |

Note carefully: the His ≥ 3 criterion was met by **every** sequence. The filter did no
discriminating work on histidine content at all; the only rejections came from the composition
sanity bound. That is itself evidence the bias was too strong.

## Where the histidines landed

Positional histidine frequency across all 24 designs:

| Position | Tier | Distance to target | H frequency |
|---:|---|---:|---:|
| 1 | core | 3.0 Å | **24/24** |
| 2 | core | 3.7 Å | 7/24 |
| 3 | rim | 5.7 Å | 18/24 |
| 4 | rim | 5.7 Å | 21/24 |
| 5 | core | 4.2 Å | **24/24** |
| 6 | rim | 5.8 Å | 18/24 |
| 40 | outer | 7.7 Å | 1/24 |
| 41 | core | 4.9 Å | **22/24** |
| 44 | rim | 5.8 Å | 21/24 |
| 45 | core | 4.8 Å | 20/24 |
| 46 | outer | 7.2 Å | 1/24 |
| 55 | outer | 7.1 Å | 18/24 |

Mean interface (core + rim) histidines per design: **7.4 of a possible 9**. Off-site histidines
at genuinely non-interface positions: **zero** in every design.

So the histidines are emphatically **not** scattered surface residues. They are exactly where
they were directed. But that is the problem, not the achievement: at positions 1, 5, 41, 4 and
44 the residue is effectively fixed across all 24 designs. ProteinMPNN was not choosing
histidine on structural merit; a +4.0 logit bias (~55× odds) simply overwhelmed its preference.
Only position 2 (7/24) shows the model exercising real discretion. **Effective interface
diversity is close to 1.**

### Two distinct clusters, with very different credibility

**Cluster 1 — N-terminal residues 1-6. Largely an artifact.**
The core tier (1, 2, 5) and the rim tier (3, 4, 6) together tile residues 1-6 contiguously in
sequence. Biasing both therefore produced a continuous histidine stretch. Six designs are
literally `HHHHHH`; most of the rest are `HVHHHH` or `HVEHHH`. An N-terminal His6 is a standard
affinity tag, and is a poor candidate for a cooperative switch: it sits at the most flexible
terminus of the binder, and a run of adjacent identical residues cannot form the specific
His–acidic pairings that give a pH switch its cooperativity. This is a consequence of how the
interface tiers were defined and was not anticipated when the bias was set.

**Cluster 2 — residues 41/44/45, the recurring `HLEHH` motif. Genuinely promising.**
This cluster sits in a structured part of the binder rather than at a terminus, and positions 41
and 45 are true ≤5 Å core contacts (4.9 Å and 4.8 Å). Histidines here are adjacent to glutamates
in the local motif, which is the arrangement that can actually produce pH-dependent behaviour:
neutral His hydrogen-bonding to Glu at pH 7.4, protonated His experiencing altered electrostatics
at pH 6.0. If any part of this design round carries a real mechanism, it is this cluster.

Position 55 (18/24) is 7.1 Å from the target — peripheral, and unlikely to contribute.

## Composition and charge sanity

| Metric | Value |
|---|---|
| Length | 57 aa, all designs (within the 10-250 aa submission limit) |
| Histidine fraction | **14.3%** (natural proteins average ≈2.3%) |
| Net charge, pH 7.4 | min −12, max +4, mean −4.5 |
| Net charge, pH 6.0, all His protonated | min −3, max +12, mean +3.6 |
| Mean charge swing on protonation | **+8.1** |
| Free cysteines | none (C omitted at design time) |

A histidine fraction of 14.3% is roughly six times the natural average. That is a developability
concern in its own right — unusual composition affects expression, solubility and potential
immunogenicity — and it is a liability the Phase 9 filter will have to weigh honestly rather
than wave through.

The +8.1 mean charge swing is large and is in the direction one would want for acid-triggered
release. But the sign of the effect on binding depends on the electrostatics of the TfR1 surface
being contacted, which has not been evaluated. A charge swing is a *prerequisite* for a switch,
not evidence of one.

## Diversity among survivors

Pairwise sequence identity across the 18 survivors: min 43.9%, **mean 63.5%**, max 80.7%. All 18
are unique sequences; 34 of 153 pairs exceed 70% identity and 3 exceed 80%.

This looks healthier than it is. The variation is concentrated in the middle of the scaffold
(roughly residues 7-40), which faces away from the target. The interface positions that actually
matter are near-invariant. For the purpose of sampling distinct binding chemistries, this set
should be treated as **one design, not eighteen**.

## Honest read: is this a pH switch?

**No — not yet, and the histidine count alone does not demonstrate one.**

What Phase 6 established is narrow but real: the pipeline can now place histidines at chosen
interface positions, and the project has gone from zero histidines in 60+ designs to 24 designs
carrying 4-10 each. That removes an absolute blocker. Three reasons it is not yet a mechanism:

1. **Saturation is not selection.** Positions 1 and 5 are histidine in 100% of designs. The model
   was not asked which positions *should* be histidine; it was told. Nothing here constitutes
   structural evidence that histidine belongs at these positions.
2. **The dominant motif is tag-like.** The N-terminal run is the single most prominent feature of
   this batch and is the least plausible switch element in it. Six designs were rejected for it
   outright, and the surviving designs are milder versions of the same artifact.
3. **The upstream pose problem is untouched.** Every validated design in this project still docks
   at human residues ~507-760 rather than the apical patch, with interface pTM 0.14-0.30. These
   24 sequences have not been structurally validated at all. Histidines on a binder that does not
   engage the intended epitope give a pH-insensitive non-binder, not a pH-selective one.

A genuine demonstration requires showing that the complex is predicted to be stable at pH 7.4 and
destabilised at pH 6.0 — that is Phase 7, and it has not been run.

## Recommendation before spending Phase 7 budget

Phase 6 should be re-run once with a corrected bias before its output is carried forward:

- **Lower the bias** from +4.0/+3.0 to roughly +1.5/+2.0. Phase 3D showed +1.25 is too weak and
  Phase 6 shows +4.0 is too strong; the useful window lies between, and a bias in that range
  would let ProteinMPNN choose *which* interface positions genuinely tolerate histidine.
- **Break up the contiguous N-terminal block.** Bias a non-adjacent subset — for example
  positions 2, 5, 41 and 45 — so the result is a distributed network rather than a poly-His run.
- **Forbid histidine runs explicitly** rather than catching them after the fact in the filter.
- Prefer designs pairing interface histidines with nearby acidic residues, which is what makes a
  switch cooperative.

This costs one ProteinMPNN call. Given that Phase 7 costs roughly two Boltz-2 calls per candidate
and the deadline is **2026-11-01**, re-running Phase 6 cheaply is better value than scoring a
batch whose main feature is a purification tag.

## Rate limiting

No rate limiting occurred. The Phase 4E backoff (catch 429, honour `Retry-After`, linear backoff)
is implemented in `scripts/nvidia_proteinmpnn_his_enriched.py` but its retry path was never
entered. One API call total.
