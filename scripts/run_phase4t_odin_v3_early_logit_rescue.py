#!/usr/bin/env python3
"""Phase 4T: independently validate four rescued Odin v3 early-logit sequences."""
from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SELECTION=ROOT/"odin_multi/v3_early_logit_rescue/selection.csv"
OUT=ROOT/"predictions/odin_v3_early_logit_pose_screen"
VALIDATOR=ROOT/"scripts/nvidia_boltz2_validate_odin_multi_selected.py"

def main():
    with SELECTION.open(newline="") as fh:
        rows=[r for r in csv.DictReader(fh) if r.get("selection_status")=="selected"]
    rows.sort(key=lambda r: float(r.get("worst_target_i_ptm") or "-inf"), reverse=True)
    if not rows:
        raise SystemExit("No early-logit rescue candidates")

    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True,exist_ok=True)

    for rank,row in enumerate(rows):
        env=os.environ.copy()
        env["ODIN_SELECTION_PATH"]=str(SELECTION)
        env["ODIN_SELECTION_RANK"]=str(rank)
        env["ODIN_VALIDATION_OUT"]=str(OUT/f"rank_{rank}")
        print(f"=== {row['design_id']} (rank {rank}) ===",flush=True)
        subprocess.run([sys.executable,str(VALIDATOR)],check=True,env=env)

    results=[]
    for p in sorted(OUT.glob("rank_*/batch_summary.json")):
        d=json.loads(p.read_text())
        d["validation_dir"]=str(p.parent.relative_to(ROOT))
        results.append(d)

    results.sort(
        key=lambda d: min(
            float(d["species"]["human"]["protein_iptm"] or 0),
            float(d["species"]["mouse"]["protein_iptm"] or 0),
        ),
        reverse=True,
    )
    winners=[d["candidate_id"] for d in results if d.get("cross_species_pose_pass")]

    (OUT/"batch_summary.json").write_text(
        json.dumps({"candidates":results},indent=2)+"\n"
    )
    (OUT/"pose_pass_ids.json").write_text(
        json.dumps({"pose_pass_ids":winners},indent=2)+"\n"
    )

    print(json.dumps({
        "ranked":[{
            "id":d["candidate_id"],
            "odin_worst_target_i_ptm":d.get("odin_selection",{}).get("worst_target_i_ptm"),
            "human_ipTM":d["species"]["human"]["protein_iptm"],
            "mouse_ipTM":d["species"]["mouse"]["protein_iptm"],
            "human_preferred":d["species"]["human"]["contacts"]["preferred_contacts"],
            "mouse_preferred":d["species"]["mouse"]["contacts"]["preferred_contacts"],
            "overlap":d["cross_species_epitope"]["overlap_count"],
            "jaccard":d["cross_species_epitope"]["jaccard"],
            "pass":d["cross_species_pose_pass"],
        } for d in results],
        "pose_pass_ids":winners,
    },indent=2))

if __name__=="__main__":
    main()
