#!/usr/bin/env python3
"""Sequence-level developability gate for pose-valid TfR1 binders.

This stage is intentionally downstream of the structural pose gate. It does not
attempt to rescue candidates that bind the wrong epitope.

The script is generic: provide a candidate CSV and a JSON file containing
pose-pass candidate IDs. It computes transparent sequence-level metrics only.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from Bio.SeqUtils.ProtParam import ProteinAnalysis


HYDRO = set("AILMVFWY")


def longest_run(seq):
    if not seq:
        return 0
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i-1] else 1
        best = max(best, run)
    return best


def max_hydrophobic_window(seq, window=9):
    if len(seq) <= window:
        return sum(x in HYDRO for x in seq), len(seq)
    best = 0
    for i in range(len(seq)-window+1):
        best = max(best, sum(x in HYDRO for x in seq[i:i+window]))
    return best, window


def load_pose_ids(path):
    obj = json.loads(path.read_text(encoding="utf-8"))
    for key in ("pose_pass_ids", "cross_species_pose_pass_ids", "passing_ids"):
        if key in obj:
            return set(obj[key])
    if isinstance(obj, list):
        return set(obj)
    raise ValueError(f"Could not find pose-pass ID list in {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--pose-ids", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    pose_ids = load_pose_ids(args.pose_ids)
    with args.candidates.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    selected = [r for r in rows if r.get("candidate_id") in pose_ids]
    args.out.mkdir(parents=True, exist_ok=True)

    results = []
    for row in selected:
        seq = row.get("sequence") or row.get("binder_sequence") or ""
        seq = seq.strip().upper()
        if not seq:
            continue

        pa = ProteinAnalysis(seq)
        acidic = seq.count("D") + seq.count("E")
        basic = seq.count("K") + seq.count("R")
        charged_fraction = (acidic + basic) / len(seq)
        max_aa_fraction = max(seq.count(x) for x in set(seq)) / len(seq)
        hyd_n, hyd_w = max_hydrophobic_window(seq)

        hard = []
        warnings = []
        if not 10 <= len(seq) <= 250:
            hard.append("submission_length")
        if longest_run(seq) > 5:
            hard.append("homopolymer_run_gt_5")
        if max_aa_fraction > 0.40:
            hard.append("single_aa_fraction_gt_0.40")
        if charged_fraction > 0.70:
            hard.append("charged_fraction_gt_0.70")
        if hyd_w >= 7 and hyd_n / hyd_w >= 0.89:
            hard.append("extreme_local_hydrophobicity")

        pI = pa.isoelectric_point()
        gravy = pa.gravy()
        instability = pa.instability_index()
        q74 = pa.charge_at_pH(7.4)
        q60 = pa.charge_at_pH(6.0)

        if pI < 4.5 or pI > 10.0:
            warnings.append("extreme_pI")
        if gravy > 0.5:
            warnings.append("high_GRAVY")
        if instability > 60:
            warnings.append("high_instability_index")
        if abs(q74) > 15:
            warnings.append("large_net_charge_pH7.4")
        if seq.count("C"):
            warnings.append("contains_cysteine")

        results.append({
            **row,
            "dev_length": len(seq),
            "dev_molecular_weight_Da": round(pa.molecular_weight(), 2),
            "dev_pI": round(pI, 3),
            "dev_GRAVY": round(gravy, 3),
            "dev_instability_index": round(instability, 2),
            "dev_charge_pH7.4": round(q74, 3),
            "dev_charge_pH6.0": round(q60, 3),
            "dev_charge_shift_7.4_to_6.0": round(q60-q74, 3),
            "dev_histidine_count": seq.count("H"),
            "dev_charged_fraction": round(charged_fraction, 3),
            "dev_max_single_aa_fraction": round(max_aa_fraction, 3),
            "dev_longest_homopolymer": longest_run(seq),
            "dev_max_hydrophobic_window": f"{hyd_n}/{hyd_w}",
            "developability_hard_flags": ";".join(hard),
            "developability_warnings": ";".join(warnings),
            "developability_pass": not hard,
        })

    out_csv = args.out / "developability.csv"
    if results:
        with out_csv.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=results[0].keys())
            w.writeheader()
            w.writerows(results)
    else:
        out_csv.write_text("candidate_id,developability_pass\n", encoding="utf-8")

    summary = {
        "pose_pass_input_count": len(pose_ids),
        "matched_candidate_count": len(selected),
        "evaluated_count": len(results),
        "developability_pass_count": sum(bool(r["developability_pass"]) for r in results),
        "pass_ids": [r["candidate_id"] for r in results if r["developability_pass"]],
        "note": (
            "These are sequence-level screening heuristics, not experimental "
            "developability measurements. pH-dependent charge is descriptive only."
        ),
    }
    (args.out / "developability_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
