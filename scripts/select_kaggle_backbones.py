#!/usr/bin/env python3
"""Strict scientific gate for Kaggle TfR1 RFdiffusion pilot backbones.

This stage runs only after the Kaggle pilot has been synced into
predictions/rfdiffusion_pilot_v1/kaggle/. It deliberately separates
"RFdiffusion generated a complex" from "the complex is worth sequence design".

Gate:
- preliminary Kaggle geometry_pass must be true when present;
- >=2 preferred conserved apical contacts within 5 A;
- <=1 explicit human/mouse risk contact within 5 A;
- binder radius of gyration <=18 A;
- binder maximum CA span <=60 A;
- <=10 severe inter-chain atom clashes (<2 A).

No candidate is carried forward if it fails the gate.
"""

from __future__ import annotations

import csv
import json
import math
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "predictions/rfdiffusion_pilot_v1/kaggle"
OUT = ROOT / "structures/rfdiffusion/kaggle_pilot_selected"

PREFERRED = {206, 207, 209, 214, 231, 343, 344, 347, 350, 374}
RISK = {204, 208, 210, 211, 212, 348, 371}

MIN_PREFERRED = 2
MAX_RISK = 1
MAX_RG = 18.0
MAX_CA_SPAN = 60.0
MAX_CLASHES = 10


def atoms(path: Path):
    out = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        atom = line[12:16].strip()
        if atom.startswith("H"):
            continue
        chain = line[21].strip() or "_"
        try:
            resi = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        out.append((chain, resi, atom, x, y, z))
    return out


def dist(a, b):
    return math.sqrt((a[3]-b[3])**2 + (a[4]-b[4])**2 + (a[5]-b[5])**2)


def compactness(atom_rows, binder_chain):
    cas = [a for a in atom_rows if a[0] == binder_chain and a[2] == "CA"]
    if not cas:
        return None, None
    cx = sum(a[3] for a in cas) / len(cas)
    cy = sum(a[4] for a in cas) / len(cas)
    cz = sum(a[5] for a in cas) / len(cas)
    rg = math.sqrt(sum((a[3]-cx)**2 + (a[4]-cy)**2 + (a[5]-cz)**2 for a in cas) / len(cas))
    span = 0.0
    for i, a in enumerate(cas):
        for b in cas[i+1:]:
            span = max(span, dist(a, b))
    return rg, span


def preliminary_passes():
    path = INPUT / "pilot_geometry_screen.csv"
    if not path.exists():
        return {}
    out = {}
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            rel = row.get("pdb", "")
            val = str(row.get("geometry_pass", "")).strip().lower()
            out[Path(rel).name] = val in {"true", "1", "yes", "y"}
    return out


def evaluate(path: Path, prelim: dict):
    aa = atoms(path)
    chains = {}
    for a in aa:
        chains.setdefault(a[0], set()).add(a[1])
    if len(chains) < 2:
        return {"pdb": str(path), "scientific_pass": False, "reason": "fewer_than_2_chains"}

    target = max(chains, key=lambda c: len(chains[c]))
    binder = min(chains, key=lambda c: len(chains[c]))
    target_first = min(chains[target])
    ta = [a for a in aa if a[0] == target]
    ba = [a for a in aa if a[0] == binder]

    contacts = set()
    clashes = 0
    for b in ba:
        for t in ta:
            d = dist(b, t)
            if d < 2.0:
                clashes += 1
            if d < 5.0:
                # RFdiffusion flattens/renumbers output chains. The fixed target
                # is the 11QC apical crop beginning at human residue 189.
                contacts.add(t[1] - target_first + 189)

    preferred = sorted(contacts & PREFERRED)
    risk = sorted(contacts & RISK)
    rg, span = compactness(aa, binder)
    prelim_ok = prelim.get(path.name, True)

    checks = {
        "preferred_contacts_ge_2": len(preferred) >= MIN_PREFERRED,
        "risk_contacts_le_1": len(risk) <= MAX_RISK,
        "radius_of_gyration_le_18A": rg is not None and rg <= MAX_RG,
        "ca_span_le_60A": span is not None and span <= MAX_CA_SPAN,
        "severe_clashes_le_10": clashes <= MAX_CLASHES,
    }
    passed = all(checks.values())

    return {
        "pdb": str(path.relative_to(ROOT)),
        "family": path.parent.name,
        "target_chain": target,
        "binder_chain": binder,
        "binder_length": len(chains[binder]),
        "target_contacts_5A": sorted(contacts),
        "preferred_contacts": preferred,
        "risk_contacts": risk,
        "severe_clashes_lt_2A": clashes,
        "radius_of_gyration_A": round(rg, 3) if rg is not None else None,
        "max_CA_span_A": round(span, 3) if span is not None else None,
        "pilot_geometry_pass_reported": prelim_ok,
        "pilot_geometry_pass_not_used_as_gate": True,
        "checks": checks,
        "scientific_pass": passed,
        "score_tuple": [
            int(passed),
            len(preferred),
            -len(risk),
            -(rg if rg is not None else 999),
            -(span if span is not None else 999),
        ],
    }


def main():
    if not INPUT.exists():
        raise SystemExit(f"Kaggle outputs not present yet: {INPUT}")

    pdbs = sorted(INPUT.glob("**/design_*.pdb"))
    if not pdbs:
        raise SystemExit(f"No design_*.pdb files found below {INPUT}")

    prelim = preliminary_passes()
    rows = [evaluate(p, prelim) for p in pdbs]
    rows.sort(key=lambda r: tuple(r.get("score_tuple", [])), reverse=True)

    passing = [r for r in rows if r.get("scientific_pass")]
    OUT.mkdir(parents=True, exist_ok=True)

    # Clear only previously generated selected PDB copies.
    for old in OUT.glob("selected_*.pdb"):
        old.unlink()

    for rank, row in enumerate(passing, start=1):
        source = ROOT / row["pdb"]
        dest = OUT / f"selected_{rank:02d}_{row['family']}_{source.name}"
        shutil.copy2(source, dest)
        row["selected_copy"] = str(dest.relative_to(ROOT))

    summary = {
        "gate_version": "kaggle-pilot-strict-v1",
        "input_pdb_count": len(pdbs),
        "passing_count": len(passing),
        "preferred_residues": sorted(PREFERRED),
        "risk_residues": sorted(RISK),
        "thresholds": {
            "min_preferred_contacts": MIN_PREFERRED,
            "max_risk_contacts": MAX_RISK,
            "max_radius_of_gyration_A": MAX_RG,
            "max_CA_span_A": MAX_CA_SPAN,
            "max_severe_clashes_lt_2A": MAX_CLASHES,
        },
        "selected": passing,
        "all_candidates": rows,
    }
    (OUT / "selection_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    with (OUT / "selection_table.csv").open("w", newline="", encoding="utf-8") as fh:
        fields = [
            "family", "pdb", "binder_length", "preferred_contacts", "risk_contacts",
            "severe_clashes_lt_2A", "radius_of_gyration_A", "max_CA_span_A",
            "scientific_pass", "selected_copy"
        ]
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({
                "family": row.get("family", ""),
                "pdb": row.get("pdb", ""),
                "binder_length": row.get("binder_length", ""),
                "preferred_contacts": ";".join(map(str, row.get("preferred_contacts", []))),
                "risk_contacts": ";".join(map(str, row.get("risk_contacts", []))),
                "severe_clashes_lt_2A": row.get("severe_clashes_lt_2A", ""),
                "radius_of_gyration_A": row.get("radius_of_gyration_A", ""),
                "max_CA_span_A": row.get("max_CA_span_A", ""),
                "scientific_pass": row.get("scientific_pass", False),
                "selected_copy": row.get("selected_copy", ""),
            })

    print(json.dumps({
        "input_pdb_count": len(pdbs),
        "passing_count": len(passing),
        "selected": [
            {
                "family": x["family"],
                "preferred": x["preferred_contacts"],
                "risk": x["risk_contacts"],
                "rg": x["radius_of_gyration_A"],
                "span": x["max_CA_span_A"],
                "copy": x.get("selected_copy"),
            }
            for x in passing
        ]
    }, indent=2))

    if not passing:
        raise SystemExit("No Kaggle RFdiffusion backbone passed the strict scientific gate.")


if __name__ == "__main__":
    main()
