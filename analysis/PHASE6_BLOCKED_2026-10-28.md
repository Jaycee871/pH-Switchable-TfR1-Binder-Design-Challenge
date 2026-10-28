# Phase 6 — authentication blocker (RESOLVED) (2026-10-28)

> **SUPERSEDED 2026-10-28.** This blocker is resolved. The `workflow` OAuth scope was
> granted, `.github/workflows/phase6-his-enriched-proteinmpnn.yml` was pushed in commit
> `13136a1`, and Phase 6 ran successfully as GitHub Actions run **40691292024**.
> Results and assessment: `analysis/PHASE6_HIS_ENRICHED_2026-10-28.md` and
> `designs/family_C_pH_switch/his_enriched_seed20261102/`.
> The document below is retained only as a record of the authentication problem and how it
> was fixed. **Phase 6 is not blocked.**

Phase 6 (histidine-enriched ProteinMPNN) is **written, validated offline, and committed, but has
not been executed.** No design sequences exist yet and no assessment of them can be made.

## Blocker

The companion workflow `.github/workflows/phase6-his-enriched-proteinmpnn.yml` could not be
pushed:

```
! [remote rejected] main -> main (refusing to allow an OAuth App to create or
  update workflow .github/workflows/phase6-his-enriched-proteinmpnn.yml
  without `workflow` scope)
```

The authenticated `gh` token holds scopes `gist, read:org, repo` — it lacks `workflow`, which
GitHub requires for any push that creates or modifies a file under `.github/workflows/`. The
separately configured GitHub connector cannot reach this private repository either (HTTP 404),
so it is not an alternative route.

Because the NVIDIA API key is a repository secret available only inside Actions, Phase 6 cannot
be run locally instead. **The workflow file is the only missing piece.**

## Resolution

The repository owner needs to grant the scope once, then the workflow can be added and
dispatched:

```
gh auth refresh -h github.com -s workflow
```

This is an interactive browser/device-code flow and must be performed by the account owner.

## What is already in place

- `scripts/nvidia_proteinmpnn_his_enriched.py` — committed in `24eb700`.
- The workflow YAML — written and YAML-validated, but present only on the working machine at
  `.github/workflows/phase6-his-enriched-proteinmpnn.yml`; it is untracked and unpushed.

### Offline validation already performed

- Bias matrix builds correctly: 57 rows × 21-letter alphabet, H at +4.0 (core), +3.0 (rim),
  +1.5 (outer), with D/E/K/R penalties at core and rim and all non-interface rows left at zero.
- Interface tiers were recomputed directly from `selected_backbone.pdb` and reproduce the
  Phase 3E positions exactly: core ≤5 Å = [1, 2, 5, 41, 45], rim 5-6 Å = [3, 4, 6, 44].
- The His ≥ 3 filter and the sanity bounds behave correctly on constructed test sequences,
  including rejection of a poly-glycine control.

One filter threshold was corrected during that offline testing. An initial
`max_single_aa_fraction` of 0.40 rejected every realistic sequence, because this backbone family
is natively 39-46% glutamate; it would have failed designs for a reason unrelated to histidine
enrichment. The hard threshold is now 0.50, and polyampholyte character is recorded as a
non-fatal note, which is properly the subject of the Phase 9 developability filter.

## Status

Phase 6 remains **not run**. No histidine-enriched sequences exist, so the project still has no
mechanism for the top-ranked pH-selectivity criterion, exactly as described in
`analysis/PHASE_STATUS_2026-10-28.md`. The deadline is 2026-11-01.
