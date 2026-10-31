# TfR1 Challenge 04 — submission record

Date of record: 2026-10-31 (Asia/Taipei).

> **Amendment — 2026-10-31 (Asia/Taipei).** After this record was first written, the two stale
> `protein` values it describes below were updated to `single_chain` at the author's
> instruction: `submission/TFR1_SUBMISSION_CHECKLIST.md` and
> `submission/tfr1_submission_metadata_final_v1.csv` were both changed on 2026-10-31 and merged
> as squash commit `e723912` (PR #2). The original text of this record is preserved unchanged
> below, with bracketed amendment notes marking the two passages that no longer describe the
> current state of the repository. What the repository looked like at submission time is
> unaffected: `submission/tfr1_submission_as_submitted_2026-10-31.csv` remains authoritative for
> `molecule_class`.

## What was submitted

The Challenge 04 TfR1 entry was submitted to the **Anthropic × Adaptyv Protein Design
Competition, round 4**, on Proteinbase:

<https://proteinbase.com/competitions/anthropic-adaptyv-2026>

The submitted file contained **16 designs** plus a header row, with the columns
`name,sequence,molecule_class`, and `molecule_class = single_chain` on every row.

The exact bytes that were uploaded are preserved in this repository as
`submission/tfr1_submission_as_submitted_2026-10-31.csv`.

## Evidence for the submission — and its strength

The only evidence that the submission completed is an **on-screen "Submission successful"
message, observed at 2026-10-31 20:39 Asia/Taipei**, reported by the author via a screenshot
shown in conversation.

This is a **user-reported screen confirmation**. It is deliberately *not* described as a
receipt. Specifically:

- it is **not** an API response or callback from Proteinbase;
- it is **not** an emailed confirmation;
- it was **not** captured from the Proteinbase account by any automated process;
- the screenshot itself is not committed to this repository.

**No submission ID, confirmation number, entry number, reviewer, queue position, or
Proteinbase-side timestamp was captured.** The 20:39 Asia/Taipei time above is the time the
confirmation was observed by the author, not a timestamp issued by Proteinbase. Nothing in
this record should be read as Proteinbase-side metadata, because none was collected.

## The 16 design names submitted

1. TF4_001_L70
2. TF4_002_XS08
3. TF4_003_F39
4. TF4_004_V5_T0
5. TF4_005_R42P
6. TF4_006_V6
7. TF4_007_V5_T1
8. TF4_009_F36
9. TF4_010_R2_58
10. TF4_012_MS15
11. TF4_013_L67
12. TF4_014_L67_XS05
13. TF4_015_ANTI_V3
14. TF4_016_SEED734966
15. TF4_017_IL4_246
16. TF4_019_R2_99

The non-contiguous numbering is expected: TF4_008, TF4_011, TF4_018 and TF4_020 were dropped
when the v0 20-design portfolio was pruned to final-v1, for the reasons itemised in
`submission/TFR1_FINAL_V1_FREEZE_NOTE.md`.

## Relationship to `tfr1_submission_upload_final_v1.csv`

**The submitted content is identical to the frozen file.** All 16 names, all 16 sequences and
all 16 `molecule_class` values match exactly, in the same row order. Verified by byte
comparison.

The **only** difference between the file the author uploaded and
`submission/tfr1_submission_upload_final_v1.csv` is **line endings**: the uploaded file used
CRLF (`\r\n`) on all 17 lines; the frozen file in the repository uses LF. Both end with a
trailing newline and neither carries a byte-order mark. There is no content difference of any
kind. `submission/tfr1_submission_as_submitted_2026-10-31.csv` is stored with LF endings and
is therefore byte-identical to the frozen file; the CRLF detail is recorded here rather than
in the file itself.

### The `single_chain` vs `protein` discrepancy

This discrepancy is real, but it is **not** between the submitted file and the frozen CSV. It
is between the submitted file and other documents in the repository:

- `submission/tfr1_submission_upload_final_v1.csv` — says `single_chain`. **Matches what was
  submitted.** It was changed from `protein` to `single_chain` in commit `9ae2a9b`
  ("Fix TfR1 molecule_class to Proteinbase single_chain"), authored 2026-10-31 20:24
  Asia/Taipei, roughly fifteen minutes before the submission confirmation was observed.
- `submission/TFR1_SUBMISSION_CHECKLIST.md` — still contains the line
  "`molecule_class=protein` for all designs in this package". **This no longer describes what
  was submitted.** *[Amended 2026-10-31: no longer true of the current repository — this line
  now reads `molecule_class=single_chain`, changed in `e723912` (PR #2). It was accurate as of
  the date of record.]*
- `submission/tfr1_submission_metadata_final_v1.csv` — still carries `molecule_class = protein`
  in all 16 rows. **This no longer describes what was submitted.** *[Amended 2026-10-31: no
  longer true of the current repository — all 16 rows now carry `single_chain`, changed in
  `e723912` (PR #2). It was accurate as of the date of record.]*
- `submission/tfr1_submission_upload_v0.csv` and `tfr1_submission_metadata_v0.csv` — the
  superseded 20-design package; not consulted here.

These stale `protein` values are **left in place deliberately**. The point of this record is an
accurate account of what went out, not a retroactive tidy-up of the frozen package or of the
history. Anyone reading the checklist or the metadata CSV should treat
`tfr1_submission_as_submitted_2026-10-31.csv` as authoritative for `molecule_class`.

*[Amended 2026-10-31 (Asia/Taipei): this paragraph described the position as of the date of
record and is superseded. The author subsequently instructed that the checklist and the
final-v1 metadata CSV be brought into line with what was submitted, and both were updated from
`protein` to `single_chain` the same day, merged as squash commit `e723912` (PR #2). The two
files named above therefore no longer carry stale `protein` values. The frozen upload CSV and
`tfr1_submission_as_submitted_2026-10-31.csv` were not altered, and the latter remains
authoritative for `molecule_class`.]*

## Internal validation status at time of submission

Stated plainly, because it is easy to lose in the surrounding detail:

**No design in this project passed the Stage 4 / Boltz-2 cross-species pose gate.** The gate
required human protein iPTM >= 0.50, mouse protein iPTM >= 0.50, retained apical/preferred
contacts in both species, homologous contact overlap >= 3, and homologous Jaccard >= 0.15.
Every one of the 25 `pose_pass_ids.json` files under `predictions/` is empty
(`{"pose_pass_ids": []}`).

The closest result anywhere in the portfolio is `TF4_007_V5_T1`, with homologous Jaccard
**0.1489** (below the 0.15 threshold), homologous overlap 7, human protein iPTM **0.2946** and
mouse protein iPTM **0.1374** — both an order of consideration below the 0.50 requirement — and
on a distal rather than apical epitope. No pH 6.0 negative-state evaluation was implemented,
so no pH-selectivity claim is made or supported.

The submitted portfolio is an uncertainty-aware experimental shortlist. It is not a set of
computationally validated cross-species or pH-switchable binders.

## Related files

- `submission/tfr1_submission_as_submitted_2026-10-31.csv` — the exact content uploaded.
- `submission/tfr1_submission_upload_final_v1.csv` — the frozen final-v1 CSV (identical content).
- `submission/TFR1_METHODS_STATEMENT.md` — the author's official methodology statement, plus a
  separately-headed maintainer fact-check appendix.
- `submission/TFR1_FINAL_V1_FREEZE_NOTE.md` — why the portfolio is 16 designs, not 20.
- `submission/TFR1_METHODS_AND_PROVENANCE.md` — the longer technical methods and provenance note.
- `submission/TFR1_SUBMISSION_CHECKLIST.md` — pre-upload checklist (the `protein` line noted here at the date of record was updated to `single_chain` on 2026-10-31; see the amendment above).
