# Project status — 2026-10-28

Honest assessment of the state of Challenge 04 (TfR1 / CD71) work, written after recovering all
GitHub Actions run artifacts into the repository. Submission deadline: **2026-11-01**.

Every number below is taken directly from committed artifact JSON. Provenance for each file is
in `analysis/ARTIFACT_PROVENANCE.md`.

## Headline

**There is currently no submittable candidate.** Every candidate that has been structurally
validated binds the wrong surface of TfR1, and no design produced so far contains a single
histidine, so no pH-switch mechanism exists.

## Candidate inventory

Human target numbering follows UniProt P02786. The intended apical epitope is approximately
residues 202-215 and 343-374. `iptm` is the protein interface pTM score reported by Boltz-2;
it measures confidence in the *relative geometry of the two chains*. "Core retained" is whether
the final 5 Å interface actually contacted the intended apical anchor residues.

| Candidate | Len | Source run | Human conf / iptm | Mouse conf / iptm | Human contacts | Mouse contacts | Core retained | His | State |
|---|---:|---|---|---|---|---|---|---:|---|
| TFR1_A29_MPNN_03 (seq-only) | 87 | 40683143306 | 0.751 / n/a | 0.443 / n/a | ~612-657 | ~651-685 | No | 0 | superseded by templated test |
| TFR1_A29_MPNN_03 (templated) | 87 | 40683978064 | 0.765 / 0.144 | 0.683 / 0.155 | ~526-654 | ~523-657 | No | 0 | **rejected** (`A29_03_final_qc.md`) |
| TFR1_A29_STD_01..16 | 87 | 40684384435 | not validated | not validated | — | — | — | 0 (all 16) | designed only |
| TFR1_A29_BIAS_02 | 87 | 40684975855 | 0.787 / 0.243 | 0.719 / 0.304 | 525-760 | shifted | No | 0 | fails pose gate |
| TFR1_A29_BIAS_04 | 87 | 40684975855 | 0.764 / 0.142 | 0.686 / 0.207 | 525-760 | shifted | No | 0 | fails pose gate |
| TFR1_A29_BIAS_01,03,05-12 | 87 | 40684712205 | not validated | not validated | — | — | — | 0 (all) | designed only |
| TFR1_C57_02 | 57 | 40686634117 | 0.766 / 0.149 | 0.741 / 0.425 | 525-651, 759 | **152-200, 384** | No | 0 | fails pose gate |
| TFR1_C57_05 | 57 | 40686634117 | 0.767 / 0.136 | 0.749 / 0.433 | 507-651 | **148-200, 384** | No | 0 | fails pose gate |
| TFR1_C57_12 | 57 | 40687192094 | 0.778 / 0.170 | 0.686 / 0.150 | 507-651 | 526-657 | No | 0 | fails pose gate |
| TFR1_C57_01,03,04,06-11 | 57 | 40686426665 | not validated | not validated | — | — | — | 0 (all) | designed only |
| TFR1_A002_MPNN_01..17 | 87/89 | 40682591391 | — | — | — | — | — | 0 (all) | early batch; `_01` is a poly-glycine artifact row |

### Backbone history

- **Smoke 001** (82 aa, run 40679577930): rejected. Interface sat directly on cross-species risk
  positions V210 / Y211 / L212 (`analysis/rfdiffusion_smoke_001_qc.md`).
- **Seeds 20261028 (89 aa) and 20261029 (87 aa)** (run 40682591391): **both failed** the strict
  gate (`requires_tier_A_contact: true`, `max_risk_contacts: 0`); `selection_summary.json`
  records `"selected": null` and `strict_cross_species_pass: false` for both. Seed 20261029 was
  nevertheless carried forward as the parent of the A29 / STD / BIAS families despite retaining a
  risk contact at residue 210. This is a known deviation from the stated gate.
- **Compact apical** (run 40685364707): seed 20261031 (64 aa) failed on compactness
  (Rg 19.57 Å > 17 Å limit); seed 20261101 (57 aa) failed with no preferred contacts;
  **seed 20261102 (57 aa) passed** (`scientific_pass: true`, preferred contacts 344 and 347,
  Rg 14.88 Å, CA span 51.02 Å). This is the only backbone that has passed its gate on its own
  terms, and it is the parent of the C57 family.

## Finding 1 — pose collapse is universal

For every validated design, the Boltz-2 human interface lands in approximately residues
**507-760**, not the apical patch at 202-215 / 343-374. `protein_iptm` is between **0.14 and
0.30** in all cases, and `preferred_core_retained` / `core_retained` is **false** in all cases.

This holds even when the RFdiffusion complex geometry is supplied to Boltz-2 as an explicit
mmCIF structural template with chain-specific assignment (Phase 4B and later). The designed
docking pose is therefore not surviving structure prediction.

Note the specific trap in these numbers: overall confidence (0.68-0.79) and complex pLDDT
(0.80-0.93) look healthy, because the individual chains fold confidently. Those values say
nothing about whether the chains are docked correctly. The low `iptm` and large `complex_ipde`
(13-18 Å) are the metrics that matter here, and they indicate the interface geometry is not
trustworthy for any candidate.

## Finding 2 — the C57_02 / C57_05 mouse numbers are an artifact, not cross-reactivity

`TFR1_C57_02` and `TFR1_C57_05` show the highest mouse `protein_iptm` in the project (0.425 and
0.433). **This must not be reported as cross-species binding.**

Their predicted mouse interface is at mouse residues **148-200 and 384**, while their own
predicted human interface is at human residues **507-651**. These are two different, unrelated
surfaces of the receptor. The predictions are not describing one epitope recognised in both
species; they are two independent docking artifacts that happen to score moderately. Genuine
cross-reactivity requires the *same* epitope in both species, which no candidate currently shows.

## Finding 3 — zero histidines across every design

All 60+ designed sequences contain **no histidine at all** (`histidine_count: 0` in every row of
every `candidates.csv`; `interface_histidines: 0` throughout the Phase 3D batch).

Phase 3D applied a position-specific bias matrix at the interface, and its code comment states
the intent to "keep polar/H-bond options and histidine available". Keeping histidine *available*
is not the same as *rewarding* it: the bias values for H were left at 0.0, so ProteinMPNN had no
incentive to select it and never did. Histidine protonation near pH 6.0 is the primary proposed
mechanism for the pH switch, so the designs presently have no physical basis for pH selectivity.

## Finding 4 — the pH 6.0 negative state does not exist

The challenge ranks **pH-selective binding first**, ahead of mouse cross-reactivity and human
affinity. No workflow in this repository evaluates pH 6.0.

`scripts/prepare_targets.py` declares `"negative": ["human_pH6.0"]` as metadata, and nothing
consumes it. A search across `scripts/` and `.github/workflows/` finds histidine only in
sequence-counting code and in three-letter-residue lookup tables. There is no protonation model,
no pH-6.0 complex prediction, and no scoring of the difference between the pH 7.4 and pH 6.0
states. This is missing code, not merely an unrun workflow.

## Finding 5 — Phase 4D failure (HTTP 429) and the Phase 4E fix

Run 40686634117 failed with:

```
RuntimeError: Boltz2 HTTP 429: {"status":429,"title":"Too Many Requests"}
  scripts/nvidia_boltz2_validate_compact_c57.py:240 -> post_json() -> :127
```

`nvidia_boltz2_validate_compact_c57.py` issues 6 sequential Boltz-2 calls (3 candidates × 2
species) separated only by a fixed `time.sleep(15)`, with **no handling for HTTP 429**. Calls
1-4 (C57_02 and C57_05, human and mouse) succeeded; call 5, C57_12 against human TfR1, hit the
NVIDIA NIM rate limit and the unhandled exception terminated the job.

Phase 4E (`scripts/nvidia_boltz2_retry_c57_12.py`) resolved this correctly: it catches
`exc.code == 429`, honours the `Retry-After` response header when present, and otherwise backs
off linearly by `20 * attempt` seconds. The retry run completed in 51 s and its backoff path was
never actually entered, indicating the limit is a short-window burst limit rather than an
exhausted daily quota.

**Operational rule for future phases:** any new script making NVIDIA API calls must reuse the
Phase 4E backoff pattern, and batches should be kept to roughly 6-8 calls per run with at least
15 s spacing.

### Other infrastructure failures on record

- Run 40684259747 (Phase 3C, first attempt): malformed ProteinMPNN call; fixed in `d75d6d7`.
- Run 40683720323 (Phase 4B, first attempt): fixed by emitting mmCIF with explicit chain IDs.
- Runs 40681292701 and 40681815723 (Phase 3): failed **by design** with
  `SystemExit("No seeded RFdiffusion backbone passed the strict gate")`. These are correct gate
  rejections, not defects.
- Runs 40688159092 and 40688311957 (Phase 5, RFD4-Proteina): HTTP 422 Unprocessable Entity, then
  HTTP 500 Internal Server Error after a `condition_spec` wire-format fix. This integration is
  unproven and is not currently on the critical path to a submission.

## What is missing before a submission can be assembled

1. A validated binding pose for any candidate — all current candidates fail the pose check.
2. Histidine-containing designs, without which pH selectivity cannot be achieved.
3. A pH 6.0 negative-state evaluation (requires new code).
4. Diversity. All 87-mers derive from one backbone (seed 20261029) and all 57-mers from one
   backbone (seed 20261102); within each family the sequences are near-identical E/K-rich
   helical bundles. Effective structural diversity is 2, not 60.
5. Developability filtering. No pI, net-charge, aggregation, cysteine, or low-complexity screen
   has been applied. By inspection the designs are extreme polyampholytes — `TFR1_A29_STD_11`
   has 41 acidic and 18 basic residues in 87 aa; `TFR1_C57_02` has 22 acidic and 20 basic in
   57 aa.
6. A submission CSV. `submissions/` is empty; the challenge requires at minimum `name`,
   `sequence`, `molecule_class`. All current candidates (57-87 aa) are within the 10-250 aa limit.

## Planned next phases

1. **Phase 6** — redesign on the validated seed-20261102 compact backbone with a *positive*
   histidine bias at interface positions; filter to `histidine_count >= 3`.
2. **Phase 7** — pH 6.0 negative-state scoring (new code), using a protonated-mimic state to
   estimate pH selectivity. This is an approximation and must be documented as such.
3. **Phase 8** — enforce the pose gate that is currently only *reported* (reject
   `core_retained: false` or `iptm < 0.5`), then cluster survivors for diversity.
4. **Phase 9** — developability filters and submission CSV assembly (no API calls).

If Phases 6-8 yield nothing that retains the designed pose, a submission should still be made
from the best-characterised candidates with their metrics reported honestly.


## Update — Phase 6B / 6C / 7A

Subsequent work sharpened the conclusion above.

- **Phase 6B:** lowering histidine bias to +2.0/+1.3 at a non-contiguous set of interface positions produced **0 histidines in 24/24 sequences**. The scaffold therefore shows a steep histidine-selection threshold rather than a smooth dose response.
- **Phase 6C:** +4.0 bias at binder position 2 alone still produced **V2 in 24/24 sequences**. Position-2 histidine seen occasionally in the broad Phase-6 batch was context-dependent, not a stable local ProteinMPNN preference.
- **Phase 7A:** two explicitly engineered V2H mutants were validated against full human and mouse TfR1. All four predictions again collapsed to the distal ~507-657 surface, with protein ipTM **0.143-0.183**, no preferred apical-core contacts, and no apical pose pass.

**Decision:** the compact seed-20261102 / C57 lane is closed as the primary route. New work is pose-first: independent apical-domain RFdiffusion backbones must pass a strict conserved-contact/risk/compactness gate, then survive full human/mouse apical pose validation before any pH-switch chemistry is engineered.

See `analysis/PHASE7A_K344_V2H_2026-10-28.md`.
