# TfR1 submission checklist

## Official hard requirements

- [ ] Challenge 04 TfR1 is still open before uploading.
- [ ] Use at most 20 designs for Track 2 / Track 3.
- [ ] CSV contains `name`, `sequence`, `molecule_class`.
- [ ] `molecule_class=protein` for all designs in this package.
- [ ] Every sequence is unique.
- [ ] Every sequence is 10–250 aa.
- [ ] No C-terminal Twin-Strep tag is included in any submitted sequence.
- [ ] Designs are de novo; no known TfR1 binder sequence was used as the starting sequence.

## Local package checks already satisfied in v0

- [x] 20 unique sequences.
- [x] All sequences use the 20 canonical amino-acid alphabet.
- [x] Length range: 60–70 aa.
- [x] Required CSV has exactly the three minimum upload columns.
- [x] Extended metadata is separated from the upload-safe CSV.
- [x] Candidate order is intentional: top row = highest internal preference.

## Manual review before final upload

- [ ] Review the top 20 order one final time after any remaining computation finishes.
- [ ] Decide whether any late candidate clearly displaces a bottom-5 design.
- [ ] Inspect TF4_004_V5_T1: it contains one cysteine; retain only if that liability is acceptable.
- [ ] Confirm the Proteinbase track/account used for this challenge.
- [ ] Upload `tfr1_submission_upload_v0.csv`.
- [ ] Add the GitHub repository as the method/provenance link.
- [ ] Provide `TFR1_METHODS_AND_PROVENANCE.md` or paste its condensed method summary in the submission form.
- [ ] Do not claim computational proof of pH selectivity.
- [ ] Do not claim that all submitted designs target the apical epitope; metadata distinguishes predicted epitopes.
- [ ] Final human review before pressing Submit.

## Freeze policy

After the final submission CSV is uploaded, do not replace it for a marginal model-score change. Replace a design only if a new candidate adds qualitatively stronger evidence (for example, sequence-only full-target human+mouse binding at the desired homologous epitope) or fixes a clear sequence liability.
