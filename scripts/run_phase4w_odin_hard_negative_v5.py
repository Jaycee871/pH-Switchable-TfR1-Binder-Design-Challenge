#!/usr/bin/env python3
from __future__ import annotations
import csv, os, subprocess, sys, shutil
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SELECTION=ROOT/"odin_multi/kaggle_hard_negative_v5/selection.csv"
OUT=ROOT/"predictions/odin_hard_negative_v5_pose_screen"
VALIDATOR=ROOT/"scripts/nvidia_boltz2_validate_odin_multi_selected.py"
AGG=ROOT/"scripts/aggregate_odin_hard_negative_v5_validation.py"

def main():
    if not SELECTION.exists():
        raise SystemExit(f"Missing V5 selection: {SELECTION}")
    with SELECTION.open(newline="") as fh:
        rows=[r for r in csv.DictReader(fh) if r.get("selection_status")=="selected"]
    rows.sort(key=lambda r: float(r.get("worst_target_i_ptm") or "-inf"), reverse=True)
    n=min(2,len(rows))
    if n<1:
        raise SystemExit("No selected V5 designs")
    if OUT.exists(): shutil.rmtree(OUT)
    OUT.mkdir(parents=True,exist_ok=True)
    for rank in range(n):
        env=os.environ.copy()
        env["ODIN_SELECTION_PATH"]=str(SELECTION)
        env["ODIN_SELECTION_RANK"]=str(rank)
        env["ODIN_VALIDATION_OUT"]=str(OUT/f"rank_{rank}")
        subprocess.run([sys.executable,str(VALIDATOR)],check=True,env=env)
    subprocess.run([sys.executable,str(AGG)],check=True)

if __name__=="__main__":
    main()
