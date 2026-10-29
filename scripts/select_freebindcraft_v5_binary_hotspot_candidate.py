#!/usr/bin/env python3
"""Select a v5 FreeBindCraft trajectory for full human/mouse validation.

Consumes the synced v5 Kaggle output. It only emits a candidate if the trajectory
already passes the independently recomputed local 5 A hotspot gate.

This keeps Boltz-2 budget downstream of a real geometric requirement rather than
FreeBindCraft's internal acceptance alone.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "freebindcraft/kaggle_v5_binary_results"
SUMMARY = RUN / "smoke_summary.json"
DESIGNS = RUN / "designs"
OUT = ROOT / "designs/family_D_freebindcraft/v5_binary_hotspot_pose_screen"
OUT.mkdir(parents=True, exist_ok=True)


def load_trajectory_rows():
    path = DESIGNS / "trajectory_stats.csv"
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as fh:
        return {r["Design"]: r for r in csv.DictReader(fh)}


def main():
    if not SUMMARY.exists():
        raise SystemExit(f"Missing v5 summary: {SUMMARY}")

    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    gates = summary.get("trajectory_local_gates", [])
    passing = [g for g in gates if g.get("true_hotspot_gate_pass")]

    decision = {
        "v5_trajectory_count": len(gates),
        "true_hotspot_gate_pass_count": len(passing),
        "selected_candidate": None,
        "reason": None,
    }

    out_csv = OUT / "pose_screen_candidates.csv"
    if not passing:
        out_csv.write_text(
            "candidate_id,backbone_id,source_backbone,binder_chain,target_chain,"
            "binder_length,sequence,mpnn_score,histidine_count,source_stage,"
            "fbc_i_ptm,fbc_ptm,fbc_plddt,fbc_i_pae,prefilter_pass\n",
            encoding="utf-8",
        )
        decision["reason"] = "No v5 trajectory passed the independent 5A hotspot gate."
        (OUT / "selection_decision.json").write_text(
            json.dumps(decision, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(decision, indent=2))
        return

    rows = load_trajectory_rows()
    # Single-attempt v5 normally has at most one trajectory, but sort defensively:
    # hotspot-core count, preferred count, fewer risk contacts.
    passing.sort(
        key=lambda g: (
            len(g.get("hotspot_core_contacts", [])),
            len(g.get("preferred_contacts", [])),
            -len(g.get("risk_contacts", [])),
        ),
        reverse=True,
    )
    gate = passing[0]
    pdb_rel = gate["pdb"]
    design = Path(pdb_rel).stem
    row = rows.get(design)
    if row is None:
        raise RuntimeError(f"No trajectory_stats row for {design}")

    seq = row["Sequence"].strip()
    hist = seq.count("H")
    source = RUN / pdb_rel

    fields = [
        "candidate_id","backbone_id","source_backbone","binder_chain","target_chain",
        "binder_length","sequence","mpnn_score","histidine_count","source_stage",
        "fbc_i_ptm","fbc_ptm","fbc_plddt","fbc_i_pae",
        "local_preferred_contacts","local_hotspot_core_contacts",
        "local_risk_contacts","prefilter_pass",
    ]
    candidate = {
        "candidate_id": design,
        "backbone_id": design.replace("TFR1_",""),
        "source_backbone": str(source.relative_to(ROOT)),
        "binder_chain": "B",
        "target_chain": "A",
        "binder_length": len(seq),
        "sequence": seq,
        "mpnn_score": "",
        "histidine_count": hist,
        "source_stage": "freebindcraft_v5_binary_hotspot_trajectory",
        "fbc_i_ptm": row.get("i_pTM",""),
        "fbc_ptm": row.get("pTM",""),
        "fbc_plddt": row.get("pLDDT",""),
        "fbc_i_pae": row.get("i_pAE",""),
        "local_preferred_contacts": ";".join(map(str, gate.get("preferred_contacts",[]))),
        "local_hotspot_core_contacts": ";".join(map(str, gate.get("hotspot_core_contacts",[]))),
        "local_risk_contacts": ";".join(map(str, gate.get("risk_contacts",[]))),
        "prefilter_pass": "True",
    }

    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerow(candidate)

    decision["selected_candidate"] = candidate["candidate_id"]
    decision["reason"] = "Passed independent 5A hotspot gate; eligible for full human/mouse validation."
    decision["gate"] = gate
    (OUT / "selection_decision.json").write_text(
        json.dumps(decision, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(decision, indent=2))


if __name__ == "__main__":
    main()
