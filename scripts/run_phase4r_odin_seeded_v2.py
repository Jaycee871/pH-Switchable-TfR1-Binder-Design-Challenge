#!/usr/bin/env python3
"""Run independent Boltz-2 validation for all selected Odin-Multi seeded-v2 designs."""
from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SELECTION=ROOT/"odin_multi/kaggle_seeded_v2/selection.csv"
OUT=ROOT/"predictions/odin_multi_seeded_v2_pose_screen"
VALIDATOR=ROOT/"scripts/nvidia_boltz2_validate_odin_multi_selected.py"
AGGREGATOR=ROOT/"scripts/aggregate_odin_seeded_v2_validation.py"

def main():
    if not SELECTION.exists():
        raise SystemExit(f"Missing selection: {SELECTION}")

    with SELECTION.open(newline="") as fh:
        rows=[r for r in csv.DictReader(fh) if r.get("selection_status")=="selected"]
    rows.sort(key=lambda r: float(r.get("worst_target_i_ptm") or "-inf"), reverse=True)
    n=min(3,len(rows))
    if n<1:
        raise SystemExit("No selected seeded-v2 Odin designs")

    if OUT.exists():
        import shutil
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True,exist_ok=True)

    print(f"Validating {n} seeded-v2 designs in descending Odin worst-target iPTM order.")
    for rank in range(n):
        env=os.environ.copy()
        env["ODIN_SELECTION_PATH"]=str(SELECTION)
        env["ODIN_SELECTION_RANK"]=str(rank)
        env["ODIN_VALIDATION_OUT"]=str(OUT/f"rank_{rank}")
        print(f"=== independent Boltz-2 validation rank {rank} ===",flush=True)
        subprocess.run([sys.executable,str(VALIDATOR)],check=True,env=env)

    subprocess.run([sys.executable,str(AGGREGATOR)],check=True)

if __name__=="__main__":
    main()
