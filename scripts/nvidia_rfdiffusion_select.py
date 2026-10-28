#!/usr/bin/env python3
"""Seeded RFdiffusion selection for the conserved-apical TfR1 family.

Runs seeds sequentially and stops at the first backbone satisfying the strict
cross-species contact gate:
- at least one Tier-A conserved contact
- zero contacts to the Proteinbase cross-species-risk set at 5 Å

Only two seeds are attempted initially to control API usage.
"""

from __future__ import annotations

import json
import math
import os
import shutil
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "structures/11QC.pdb"
OUTDIR = ROOT / "designs/family_A_conserved/seeded_selection"
URL = "https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate"

HOTSPOTS = ["A343", "A344", "A347", "A350", "A374"]
CONTIGS = "A122-760/0 70-90"
SEEDS = [20261028, 20261029]
TIER_A = {206, 207, 209, 214, 231, 343, 344, 347, 350, 374}
RISK = {204, 208, 210, 211, 212, 348, 371}


def extract_chain_a(pdb_text: str) -> str:
    out = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            if len(line) >= 22 and line[21] == "A":
                out.append(line)
        elif line.startswith(("CRYST1", "REMARK", "TITLE ")):
            out.append(line)
    out.extend(["TER", "END"])
    return "\n".join(out) + "\n"


def post_json(payload: dict, api_key: str) -> dict:
    req = Request(
        URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {api_key}",
        },
        method="POST",
    )
    try:
        with urlopen(req, timeout=900) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"RFdiffusion HTTP {exc.code}: {detail}") from exc


def parse_atoms(pdb_text: str):
    atoms = {}
    residues = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM  ") or len(line) < 54:
            continue
        chain = line[21].strip() or "_"
        try:
            resid = int(line[22:26])
            xyz = tuple(float(line[i:j]) for i, j in ((30, 38), (38, 46), (46, 54)))
        except ValueError:
            continue
        atoms.setdefault(chain, []).append((resid, xyz))
        residues.setdefault(chain, set()).add(resid)
    return atoms, residues


def distance(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def analyze_contacts(pdb_text: str, cutoff: float = 5.0) -> dict:
    atoms, residues = parse_atoms(pdb_text)
    if len(residues) < 2:
        raise RuntimeError("Expected binder + target chains")

    ordered = sorted(residues, key=lambda c: len(residues[c]))
    binder_chain = ordered[0]
    target_chain = ordered[-1]
    target_first = min(residues[target_chain])
    binder_xyz = [xyz for _, xyz in atoms[binder_chain]]

    min_dist_by_original = {}
    for out_resid, xyz in atoms[target_chain]:
        original = out_resid - target_first + 122
        best = min(distance(xyz, bxyz) for bxyz in binder_xyz)
        if original not in min_dist_by_original or best < min_dist_by_original[original]:
            min_dist_by_original[original] = best

    contacts = sorted(r for r, d in min_dist_by_original.items() if d <= cutoff)
    tier_a = sorted(set(contacts) & TIER_A)
    risk = sorted(set(contacts) & RISK)
    strict_pass = bool(tier_a) and len(risk) == 0

    return {
        "binder_chain": binder_chain,
        "binder_length": len(residues[binder_chain]),
        "target_chain": target_chain,
        "target_length": len(residues[target_chain]),
        "contact_cutoff_A": cutoff,
        "human_equivalent_contacts": contacts,
        "tier_A_contacts": tier_a,
        "risk_contacts": risk,
        "strict_cross_species_pass": strict_pass,
    }


def main() -> None:
    api_key = os.getenv("NVIDIA_API_KEY") or os.getenv("NGC_API_KEY")
    if not api_key:
        raise RuntimeError("Missing NVIDIA_API_KEY / NGC_API_KEY")
    if not INPUT.exists():
        raise RuntimeError("Missing structures/11QC.pdb")

    target_pdb = extract_chain_a(INPUT.read_text(encoding="utf-8"))
    OUTDIR.mkdir(parents=True, exist_ok=True)

    attempts = []
    selected = None

    for seed in SEEDS:
        payload = {
            "input_pdb": target_pdb,
            "contigs": CONTIGS,
            "hotspot_res": HOTSPOTS,
            "diffusion_steps": 50,
            "random_seed": seed,
        }
        result = post_json(payload, api_key)
        output_pdb = result.get("output_pdb")
        if not output_pdb:
            raise RuntimeError(
                f"No output_pdb for seed {seed}; keys={sorted(result.keys())}"
            )

        pdb_path = OUTDIR / f"rfdiffusion_seed_{seed}.pdb"
        pdb_path.write_text(output_pdb, encoding="utf-8")
        contact = analyze_contacts(output_pdb)
        record = {
            "seed": seed,
            "elapsed_ms": result.get("elapsed_ms"),
            "pdb": str(pdb_path.relative_to(ROOT)),
            "contact_analysis": contact,
        }
        attempts.append(record)
        print(json.dumps(record, indent=2))

        if contact["strict_cross_species_pass"]:
            selected = record
            shutil.copyfile(pdb_path, OUTDIR / "selected_backbone.pdb")
            break

    summary = {
        "hotspot_res": HOTSPOTS,
        "contigs": CONTIGS,
        "seeds_requested": SEEDS,
        "attempts_run": len(attempts),
        "attempts": attempts,
        "selected": selected,
        "strict_gate": {
            "requires_tier_A_contact": True,
            "max_risk_contacts": 0,
        },
    }
    (OUTDIR / "selection_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    if selected is None:
        raise SystemExit("No seeded RFdiffusion backbone passed the strict gate")

    print("Selected RFdiffusion backbone")
    print(json.dumps(selected, indent=2))


if __name__ == "__main__":
    main()
