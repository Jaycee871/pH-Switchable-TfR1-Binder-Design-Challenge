# Artifact provenance

All computational results in this repository were produced by GitHub Actions workflow runs in
this repository. The pipeline historically wrote its outputs **only** to Actions run artifacts
and did not commit them. Those artifacts are uploaded with `retention-days: 14`, so they
**expire 14 days after their run** (the earliest, from 2026-10-28, expire on 2026-11-11).

This commit recovers that evidence into the repository so the work is reproducible and
auditable after artifact expiry. Files were downloaded with `gh run download` on 2026-10-28 and
copied unmodified; no values were edited, recomputed, or regenerated.

## Run index

| Run ID | Workflow | Conclusion | Artifact | Committed to |
|---|---|---|---|---|
| 40679036292 | Phase 1 Target Preparation | success | tfr1-phase1-target-assets | `metadata/target_manifest.json`, `structures/11QC.pdb`, `structures/11RN.pdb` |
| 40679577930 | Phase 2 RFdiffusion Smoke Test | success | tfr1-rfdiffusion-smoke | `structures/rfdiffusion/smoke_001/` |
| 40682591391 | Phase 3B ProteinMPNN on Seed 20261029 | success | tfr1-seed20261029-proteinmpnn | `designs/family_A_conserved/proteinmpnn_seeded_seed20261029/`, `structures/rfdiffusion/seeded_selection/` |
| 40683143306 | Phase 4 Boltz2 Cross-Species Smoke | success | tfr1-A29-03-boltz2-human-mouse | `predictions/A29_03_seq_only/` |
| 40683978064 | Phase 4B Templated Boltz2 Cross-Species | success | tfr1-A29-03-boltz2-templated-human-mouse | `predictions/A29_03_templated/` |
| 40684384435 | Phase 3C Standard ProteinMPNN | success | tfr1-seed20261029-proteinmpnn-standard | `designs/family_A_conserved/proteinmpnn_standard_seed20261029/` |
| 40684712205 | Phase 3D Interface-Biased ProteinMPNN | success | tfr1-seed20261029-proteinmpnn-interface-biased | `designs/family_A_conserved/proteinmpnn_interface_biased_seed20261029/` |
| 40684975855 | Phase 4C Validate Interface-Biased Candidates | success | tfr1-interface-biased-boltz2-validation | `predictions/interface_biased_validation/` |
| 40685364707 | Phase 2B Compact Apical RFdiffusion | success | tfr1-compact-apical-rfdiffusion | `structures/rfdiffusion/compact_apical/` |
| 40686426665 | Phase 3E Compact Backbone ProteinMPNN | success | tfr1-compact-seed20261102-proteinmpnn | `designs/family_A_conserved/compact_seed20261102_proteinmpnn/` |
| 40686634117 | Phase 4D Validate Compact C57 Candidates | **failure (HTTP 429)** | tfr1-compact-C57-boltz2-validation | `predictions/compact_C57_validation/` |
| 40687192094 | Phase 4E Retry C57_12 Boltz2 | success | tfr1-C57-12-boltz2-retry | `predictions/C57_12_retry/` |

### Note on run 40686634117 (Phase 4D)

This run **failed**, but its artifact is still partially valid and is committed here. The job
uploaded results before exiting because `actions/upload-artifact` ran unconditionally with
`if-no-files-found: warn`. The artifact therefore contains **complete** human and mouse
predictions for `TFR1_C57_02` and `TFR1_C57_05`, and **only the input template** for
`TFR1_C57_12` — the C57_12 human prediction call is where the HTTP 429 occurred. The missing
C57_12 predictions were produced separately by Phase 4E (run 40687192094) and are committed
under `predictions/C57_12_retry/`.

## Secret handling

Every committed file was scanned for `nvapi-`, `Bearer`, `api_key`, `authorization`, and
`NVIDIA_API` before staging. **No matches were found and nothing had to be excluded.**

The `request_metadata.json` files are committed deliberately: they record endpoint URLs and
model parameters (contigs, hotspots, sampling temperatures, seeds, `use_soluble_model`, omitted
amino acids) but contain no credentials and no raw HTTP request bodies. They are the primary
record of how each model call was parameterised. The `NVIDIA_API_KEY` is supplied to workflows
as a repository secret and is masked in logs; it does not appear in any artifact.

## Not committed

- Raw Actions logs (large, and reproducible via `gh run view <id> --log` until they rotate).
- `designs/family_B_multistate/` and `designs/family_C_pH_switch/` remain empty because **no
  multi-state or explicit pH-switch designs have been produced yet**. See
  `analysis/PHASE_STATUS_2026-10-28.md`.
