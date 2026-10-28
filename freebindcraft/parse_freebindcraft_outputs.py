#!/usr/bin/env python3
"""Parse and independently gate FreeBindCraft accepted TfR1 designs.

Usage:
  python freebindcraft/parse_freebindcraft_outputs.py \
      --design-dir path/to/freebindcraft/designs \
      --out-dir freebindcraft/parsed_candidates

This is a post-processing gate only. It does not treat FreeBindCraft acceptance
as evidence of human/mouse cross-reactivity or pH switching.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

PREFERRED = {206, 207, 209, 214, 231, 343, 344, 347, 350, 374}
RISK = {204, 208, 210, 211, 212, 348, 371}
HOTSPOT_CORE = {343, 344, 347, 350}


def atom_rows(path: Path):
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM") or len(line) < 54:
            continue
        atom = line[12:16].strip()
        if atom.startswith("H"):
            continue
        chain = line[21].strip() or "_"
        try:
            resid = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        rows.append((chain, resid, atom, x, y, z))
    return rows


def dist(a, b):
    return math.sqrt((a[3]-b[3])**2 + (a[4]-b[4])**2 + (a[5]-b[5])**2)


def contact_metrics(path: Path):
    atoms = atom_rows(path)
    chains = {}
    for a in atoms:
        chains.setdefault(a[0], set()).add(a[1])
    if len(chains) < 2:
        return {
            "pdb": str(path),
            "contact_error": "fewer_than_two_chains",
        }

    # FreeBindCraft target settings use TfR1 chain A. Accepted complexes normally
    # keep the target as A and the designed binder as B. Fall back to the longest
    # chain as target only if A is absent.
    target = "A" if "A" in chains else max(chains, key=lambda c: len(chains[c]))
    non_target = [c for c in chains if c != target]
    binder = "B" if "B" in non_target else min(non_target, key=lambda c: len(chains[c]))

    ta = [a for a in atoms if a[0] == target]
    ba = [a for a in atoms if a[0] == binder]
    contacts = set()
    binder_contacts = set()
    clashes = 0

    for b in ba:
        hit = False
        for t in ta:
            d = dist(b, t)
            if d < 2.0:
                clashes += 1
            if d <= 5.0:
                contacts.add(t[1])
                hit = True
        if hit:
            binder_contacts.add(b[1])

    preferred = sorted(contacts & PREFERRED)
    risk = sorted(contacts & RISK)
    hotspot_core = sorted(contacts & HOTSPOT_CORE)

    return {
        "target_chain": target,
        "binder_chain": binder,
        "target_length": len(chains[target]),
        "binder_length_from_pdb": len(chains[binder]),
        "target_contacts_5A": sorted(contacts),
        "binder_contact_positions_5A": sorted(binder_contacts),
        "preferred_contacts": preferred,
        "risk_contacts": risk,
        "hotspot_core_contacts": hotspot_core,
        "severe_atom_clashes_lt_2A": clashes,
        "hotspot_core_ge_1": len(hotspot_core) >= 1,
        "preferred_ge_2": len(preferred) >= 2,
        "risk_le_1": len(risk) <= 1,
    }


def find_accepted_pdb(design_dir: Path, design: str):
    accepted = design_dir / "Accepted"
    candidates = [
        accepted / f"{design}.pdb",
        *accepted.glob(f"**/{design}.pdb") if accepted.exists() else [],
    ]
    for p in candidates:
        if p.exists():
            return p
    hits = list(design_dir.glob(f"**/{design}.pdb"))
    return hits[0] if hits else None


def fnum(row, key):
    val = (row.get(key) or "").strip()
    try:
        return float(val)
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    final_csv = args.design_dir / "final_design_stats.csv"
    if not final_csv.exists():
        raise SystemExit(f"Missing {final_csv}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    with final_csv.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    parsed = []
    for row in rows:
        design = row.get("Design", "").strip()
        if not design:
            continue
        pdb = find_accepted_pdb(args.design_dir, design)
        cm = contact_metrics(pdb) if pdb else {"contact_error": "accepted_pdb_not_found"}

        seq = (row.get("Sequence") or "").strip()
        iptm = fnum(row, "Average_i_pTM")
        ipsae = fnum(row, "Average_ipSAE")
        plddt = fnum(row, "Average_pLDDT")
        hotspot_rmsd = fnum(row, "Average_Hotspot_RMSD")

        independent_pass = bool(
            pdb
            and cm.get("hotspot_core_ge_1")
            and cm.get("preferred_ge_2")
            and cm.get("risk_le_1")
            and cm.get("severe_atom_clashes_lt_2A", 999) <= 10
        )

        parsed.append({
            "rank": row.get("Rank", ""),
            "design": design,
            "sequence": seq,
            "length": len(seq) if seq else row.get("Length", ""),
            "freebindcraft_average_i_pTM": iptm,
            "freebindcraft_average_ipSAE": ipsae,
            "freebindcraft_average_pLDDT": plddt,
            "freebindcraft_hotspot_RMSD_A": hotspot_rmsd,
            "accepted_pdb": str(pdb) if pdb else "",
            **cm,
            "independent_apical_gate_pass": independent_pass,
            "full_human_mouse_validation_required": True,
        })

    (args.out_dir / "parsed_candidates.json").write_text(
        json.dumps(parsed, indent=2) + "\n", encoding="utf-8"
    )

    flat = []
    for r in parsed:
        x = dict(r)
        for key in (
            "target_contacts_5A", "binder_contact_positions_5A",
            "preferred_contacts", "risk_contacts", "hotspot_core_contacts",
        ):
            if key in x and isinstance(x[key], list):
                x[key] = ";".join(map(str, x[key]))
        flat.append(x)

    if flat:
        fields = []
        for r in flat:
            for k in r:
                if k not in fields:
                    fields.append(k)
        with (args.out_dir / "parsed_candidates.csv").open(
            "w", newline="", encoding="utf-8"
        ) as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(flat)

    summary = {
        "freebindcraft_accepted_count": len(rows),
        "parsed_count": len(parsed),
        "independent_apical_gate_pass_count": sum(
            bool(x["independent_apical_gate_pass"]) for x in parsed
        ),
        "independent_apical_gate_pass_ids": [
            x["design"] for x in parsed if x["independent_apical_gate_pass"]
        ],
        "note": (
            "Passing this local apical gate is necessary but not sufficient. "
            "Full human/mouse ectodomain validation and pH-state analysis remain required."
        ),
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
