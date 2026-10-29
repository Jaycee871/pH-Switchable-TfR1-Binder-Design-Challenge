#!/usr/bin/env python3
"""Aggregate independent Boltz-2 validation outputs for Odin-Multi seeded v2."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "predictions/odin_multi_seeded_v2_pose_screen"

def main():
    results=[]
    for p in sorted(OUT.glob("rank_*/batch_summary.json")):
        d=json.loads(p.read_text())
        d["validation_dir"]=str(p.parent.relative_to(ROOT))
        results.append(d)

    if not results:
        raise SystemExit("No rank_*/batch_summary.json files found")

    results.sort(
        key=lambda d: min(
            float(d["species"]["human"]["protein_iptm"] or 0),
            float(d["species"]["mouse"]["protein_iptm"] or 0),
        ),
        reverse=True,
    )

    (OUT/"batch_summary.json").write_text(
        json.dumps({"candidates":results},indent=2)+"\n"
    )
    winners=[d["candidate_id"] for d in results if d.get("cross_species_pose_pass")]
    (OUT/"pose_pass_ids.json").write_text(
        json.dumps({"pose_pass_ids":winners},indent=2)+"\n"
    )

    report={
        "ranked":[{
            "id":d["candidate_id"],
            "odin_worst_target_i_ptm":d.get("odin_selection",{}).get("worst_target_i_ptm"),
            "human_ipTM":d["species"]["human"]["protein_iptm"],
            "mouse_ipTM":d["species"]["mouse"]["protein_iptm"],
            "overlap":d["cross_species_epitope"]["overlap_count"],
            "jaccard":d["cross_species_epitope"]["jaccard"],
            "pass":d["cross_species_pose_pass"],
        } for d in results],
        "pose_pass_ids":winners,
    }
    print(json.dumps(report,indent=2))

if __name__=="__main__":
    main()
