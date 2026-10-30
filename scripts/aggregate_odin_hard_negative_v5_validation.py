#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"predictions/odin_hard_negative_v5_pose_screen"

def main():
    results=[]
    for p in sorted(OUT.glob("rank_*/batch_summary.json")):
        d=json.loads(p.read_text())
        d["validation_dir"]=str(p.parent.relative_to(ROOT))
        results.append(d)
    if not results: raise SystemExit("No V5 rank validation outputs")
    results.sort(
        key=lambda d: (
            bool(d.get("cross_species_pose_pass")),
            min(float(d["species"]["human"]["protein_iptm"] or 0),float(d["species"]["mouse"]["protein_iptm"] or 0)),
            d["cross_species_epitope"]["jaccard"],
        ),
        reverse=True,
    )
    winners=[d["candidate_id"] for d in results if d.get("cross_species_pose_pass")]
    summary={
        "candidates":results,
        "ranked":[{
            "id":d["candidate_id"],
            "odin_worst_target_i_ptm":d.get("odin_selection",{}).get("worst_target_i_ptm"),
            "human_ipTM":d["species"]["human"]["protein_iptm"],
            "mouse_ipTM":d["species"]["mouse"]["protein_iptm"],
            "human_preferred":d["species"]["human"]["contacts"]["preferred_contacts"],
            "mouse_preferred":d["species"]["mouse"]["contacts"]["preferred_contacts"],
            "overlap_count":d["cross_species_epitope"]["overlap_count"],
            "jaccard":d["cross_species_epitope"]["jaccard"],
            "pass":d["cross_species_pose_pass"],
        } for d in results],
        "pose_pass_ids":winners,
    }
    (OUT/"batch_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (OUT/"pose_pass_ids.json").write_text(json.dumps({"pose_pass_ids":winners},indent=2)+"\n")
    print(json.dumps(summary["ranked"],indent=2))

if __name__=="__main__":
    main()
