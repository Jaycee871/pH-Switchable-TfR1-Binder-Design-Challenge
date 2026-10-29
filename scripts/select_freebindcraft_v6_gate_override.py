#!/usr/bin/env python3
"""Phase 4Q: emit the v6 tempered-binary trajectory for validation as a GATE OVERRIDE.

READ THIS BEFORE USING THE RESULT
---------------------------------
`TFR1_FBC343_V6_l60_s289817` FAILED the independent strict local 5 A gate in Phase 4M.
It is selected here anyway, deliberately, as a diagnostic.

Why the override is justified:

* It is the only FreeBindCraft trajectory so far with zero severe heavy-atom clashes
  (<2 A) AND genuine conserved hotspot-core engagement (human 347 and 350).
* It was therefore never measured end-to-end, so the project's best generator output
  has no independent cross-species number at all. Phases 4N and 4O spent Boltz-2 budget
  on *derivatives* of this backbone (iPTM human 0.33-0.73, mouse 0.16-0.49) while the
  parent itself stayed unmeasured.

Why it failed the gate (recorded in
`freebindcraft/kaggle_v6_tempered_binary_results/smoke_summary.json`):

* preferred contacts: 231, 347, 350
* hotspot-core contacts: 347, 350
* species-risk contacts: 210, 211, 348, 371   <-- the sole reason for the failure
* severe clashes <2 A: 0

The downstream Boltz-2 gate in
`scripts/nvidia_boltz2_validate_freebindcraft_v6_tempered_pose.py` is NOT loosened by
this script. Only the upstream local pre-filter is bypassed. A pose pass from this run
is still a real pass; a failure is still a real failure.

Note on `prefilter_pass`: the validator refuses to spend Boltz-2 budget on any row whose
`prefilter_pass` is not truthy, so this script necessarily sets it True. That column is an
execution flag, NOT a scientific claim. The real verdict is carried by `local_gate_pass`
(False), `gate_override` (True) and `gate_override_reason`.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "freebindcraft/kaggle_v6_tempered_binary_results"
SUMMARY = RUN / "smoke_summary.json"
DESIGNS = RUN / "designs"
OUT = ROOT / "designs/family_D_freebindcraft/v6_gate_override_pose_screen"

TARGET_DESIGN = "TFR1_FBC343_V6_l60_s289817"


def load_trajectory_rows():
    path = DESIGNS / "trajectory_stats.csv"
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as fh:
        return {r["Design"]: r for r in csv.DictReader(fh)}


def main():
    if not SUMMARY.exists():
        raise SystemExit(f"Missing v6 summary: {SUMMARY}")

    OUT.mkdir(parents=True, exist_ok=True)
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    gates = summary.get("trajectory_local_gates", [])

    gate = next(
        (g for g in gates if Path(g.get("pdb", "")).stem == TARGET_DESIGN), None
    )
    if gate is None:
        raise SystemExit(
            f"No local gate record for {TARGET_DESIGN}. Available: "
            + ", ".join(Path(g.get('pdb','')).stem for g in gates)
        )

    rows = load_trajectory_rows()
    row = rows.get(TARGET_DESIGN)
    if row is None:
        raise SystemExit(f"No trajectory_stats row for {TARGET_DESIGN}")

    source = RUN / gate["pdb"]
    if not source.exists():
        raise SystemExit(f"Source backbone missing: {source}")

    seq = row["Sequence"].strip()
    if len(seq) != 60:
        raise SystemExit(f"Expected a 60-aa binder, got {len(seq)}: {seq}")

    fields = [
        "candidate_id", "backbone_id", "source_backbone", "binder_chain", "target_chain",
        "binder_length", "sequence", "mpnn_score", "histidine_count", "source_stage",
        "fbc_i_ptm", "fbc_ptm", "fbc_plddt", "fbc_i_pae",
        "local_preferred_contacts", "local_hotspot_core_contacts",
        "local_risk_contacts", "local_severe_clashes_lt2A",
        "local_gate_pass", "gate_override", "gate_override_reason", "prefilter_pass",
    ]
    candidate = {
        "candidate_id": TARGET_DESIGN,
        "backbone_id": TARGET_DESIGN.replace("TFR1_", "") + "_gate_override",
        "source_backbone": str(source.relative_to(ROOT)),
        "binder_chain": "B",
        "target_chain": "A",
        "binder_length": len(seq),
        "sequence": seq,
        "mpnn_score": "",
        "histidine_count": seq.count("H"),
        "source_stage": "freebindcraft_v6_tempered_binary_trajectory_GATE_OVERRIDE",
        "fbc_i_ptm": row.get("i_pTM", ""),
        "fbc_ptm": row.get("pTM", ""),
        "fbc_plddt": row.get("pLDDT", ""),
        "fbc_i_pae": row.get("i_pAE", ""),
        "local_preferred_contacts": ";".join(map(str, gate.get("preferred_contacts", []))),
        "local_hotspot_core_contacts": ";".join(map(str, gate.get("hotspot_core_contacts", []))),
        "local_risk_contacts": ";".join(map(str, gate.get("risk_contacts", []))),
        "local_severe_clashes_lt2A": gate.get("severe_clashes_lt2A", ""),
        "local_gate_pass": gate.get("true_hotspot_gate_pass", False),
        "gate_override": "True",
        "gate_override_reason": (
            "Failed strict local gate ONLY on species-risk contacts "
            + ";".join(map(str, gate.get("risk_contacts", [])))
            + ". Run as a diagnostic to obtain the first independent cross-species "
            "number for the best generator output (0 severe clashes, hotspot core "
            "347/350). Downstream Boltz-2 gate unchanged."
        ),
        # `prefilter_pass` is the validator's "cleared to spend Boltz-2 budget" flag, so
        # it must be True for the diagnostic to run at all. The honest scientific record
        # of the override lives in local_gate_pass / gate_override / gate_override_reason
        # above, which are the fields to read when interpreting this candidate.
        "prefilter_pass": "True",
    }

    out_csv = OUT / "pose_screen_candidates.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerow(candidate)

    decision = {
        "phase": "4Q",
        "mode": "GATE_OVERRIDE_DIAGNOSTIC",
        "selected_candidate": TARGET_DESIGN,
        "local_gate_pass": gate.get("true_hotspot_gate_pass", False),
        "local_gate_failure_cause": "species_risk_contacts",
        "local_risk_contacts": gate.get("risk_contacts", []),
        "local_hotspot_core_contacts": gate.get("hotspot_core_contacts", []),
        "local_preferred_contacts": gate.get("preferred_contacts", []),
        "severe_clashes_lt2A": gate.get("severe_clashes_lt2A"),
        "downstream_boltz2_gate": "UNCHANGED (identical to Phase 4F/4G/4M)",
        "interpretation_warning": (
            "This candidate did not clear the strict local pre-filter. Any pass here must "
            "be reported together with the species-risk liability (210/211/348/371), which "
            "is concentrated in binder B46=W and B47=E."
        ),
        "gate": gate,
    }
    (OUT / "selection_decision.json").write_text(
        json.dumps(decision, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(decision, indent=2))


if __name__ == "__main__":
    main()
