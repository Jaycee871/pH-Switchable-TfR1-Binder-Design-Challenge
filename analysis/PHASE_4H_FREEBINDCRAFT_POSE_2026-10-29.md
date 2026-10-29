# Phase 4H — Boltz-2 cross-species pose validation of the FreeBindCraft trajectory design

Date: 2026-10-29. Run: GitHub Actions `40814596579`. Deadline: **2026-11-01 23:59 AoE**.

Candidate: `TFR1_FBC343_l70_s574476`, the best of the three trajectories harvested from the
cancelled Kaggle FreeBindCraft smoke run (`freebindcraft/kaggle_smoke_v1_output/`, commit
`cebb5ae`).

Sequence (70 aa, **1 histidine**, position 27):

```
MEVQTLAVAEMMEILKNPKTWWGYGRHPSGTNTGRWEIVFFTEQEDPRQDWEITVVEVEEEQWREFMAQT
```

Recovered independently from two sources that agree exactly: the `trajectory_stats.csv` row and
a CA-trace extraction of chain B of
`designs/Trajectory/Relaxed/TFR1_FBC343_l70_s574476.pdb`.

## Result: FAIL

`predictions/freebindcraft_trajectory_pose_screen/pose_pass_ids.json` is `{"pose_pass_ids": []}`.

The gate is the **unchanged** Phase 4F/4G gate. `scripts/nvidia_boltz2_validate_freebindcraft_pose.py`
is a copy of `scripts/nvidia_boltz2_validate_rescue_pose.py` with only the input CSV path, the
output directory, and a hydrogen-stripping line in `make_template` changed (the OpenMM-relaxed
FreeBindCraft PDB carries 2038 hydrogens; earlier RFdiffusion inputs were heavy-atom only, so
stripping them makes the Boltz-2 template construction identical to earlier phases). No threshold
was altered: `protein_iptm >= 0.50`, `preferred_core_retained`, apical contact,
`overlap >= 3 and jaccard >= 0.15`.

| Metric | Human | Mouse | Gate |
|---|---:|---:|---|
| protein ipTM | **0.432** | **0.677** | >= 0.50 |
| pTM | 0.714 | 0.444 | — |
| complex pLDDT | 0.814 | 0.612 | — |
| complex ipDE | 12.40 | 12.64 | — |
| Preferred apical contacts | **none** | **none** | >= 1 required |
| `preferred_core_retained` | **false** | **false** | must be true |
| Species-risk contacts | none | none | — |
| `apical_contact` (189-383) | true | true | required |
| Per-species pose gate | **fail** | **fail** | — |

Human 5 Å target contacts (UniProt P02786): 154, 161, 193, 195, 197, 215, 216, 217, 290, 291,
292, 294, 295, 296, **565, 569, 570, 572, 574, 575, 578, 579**.

Mouse 5 Å target contacts (human numbering): 152-158, 198, 219, 261-266, 285-302, 386, 411-423,
458-462, **489-539**, 571-573 — **73 residues**.

Cross-species epitope: homologous overlap 10 residues (154, 217, 290-296, 569, 570),
**Jaccard 0.1205 < 0.15**, `same_epitope_pass` false. `cross_species_pose_pass` false.

## What this actually means

1. **The apical anchor is still not engaged.** Zero preferred contacts in either species. The
   intended anchors 206/207/209/214/231/343/344/347/350/374 are untouched. Human does clip the
   far edge of the 202-215 band (193, 195, 197, 215, 216, 217) but never the core residues, and
   touches nothing in 343-374 at all.
2. **Pose collapse is not cured, only diluted.** Human still puts 8 of 22 contacts in the
   565-579 region — the same protease-like/helical-domain surface that every NVIDIA-lane
   candidate collapsed onto.
3. **The mouse ipTM of 0.677 is not a pass in disguise.** It clears the numeric threshold, but
   pTM is 0.444 and complex pLDDT is 0.612, and the binder smears across 73 target residues
   spanning 152-573 — essentially the whole ectodomain. That is a non-specific wrap with
   confident relative placement, not a defined epitope. Taking it as encouraging would be
   self-deception.
4. **Boltz-2 does not reproduce FreeBindCraft's own score.** FreeBindCraft reported i_pTM 0.69
   for this trajectory; Boltz-2 gives 0.432 for human. That gap is expected and is the point of
   the gate: the 0.69 is an AF2 hallucination loop scoring the objective it was optimising, so
   it is not independent evidence. This design also never cleared FreeBindCraft's *own* final
   filters — `final_design_stats.csv` and `mpnn_design_stats.csv` are header-only and 12 MPNN
   redesigns were rejected (11 on i_pTM, 12 on i_pAE).
5. **The hotspot specification did not take effect.** `target_settings.json` requested
   `target_hotspot_residues: "343,344,347,350"`, yet the relaxed trajectory PDB itself — before
   Boltz-2 is involved at all — contacts only 217, 218, 259-264 and 283-300 of the target, with
   zero contacts at 343/344/347/350 and zero in 202-215. `Trajectory_WrongHotspot` in
   `failure_csv.csv` is 0, i.e. BindCraft never flagged this, which means the hotspot restraint
   was not being enforced on the trajectory. Any future FreeBindCraft run must be verified on
   this point first.

## Comparison with the previous best

| | `TFR1_NB01_conserved_343_s20261124_08` | `TFR1_FBC343_l70_s574476` |
|---|---|---|
| Human ipTM | 0.202 | **0.432** |
| Mouse ipTM | 0.173 | **0.677** |
| Human `core_retained` | **true** (207, 209, 214, 344, 374) | **false** |
| Mouse `core_retained` | false | false |
| Epitope overlap / Jaccard | 0 / 0.00 | **10 / 0.121** |
| Risk contacts | 208, 210, 211, 212 | none |
| Pose gate | fail | fail |

Honest reading: this candidate is **better on confidence and worse on specificity**. It roughly
doubles the best human ipTM ever recorded in this project and is the first design to produce any
cross-species epitope overlap at all, but it is also the first to lose `preferred_core_retained`
on the human side, which NB01_08 had. On the criterion that actually decides the challenge —
does it dock at the apical patch — it is a regression. It is not a near miss.

## Standing

Eight candidates have now been through the cross-species pose gate. `pose_pass_ids` is empty in
all three screens (`nvidia_apical_pose_screen`, `nvidia_apical_rescue_pose_screen`,
`freebindcraft_trajectory_pose_screen`). No pH 6.0 state has been modelled anywhere in the
repository. There is still no submittable candidate, and histidine engineering on this design
would be engineering a pH switch into an interface that does not exist.
