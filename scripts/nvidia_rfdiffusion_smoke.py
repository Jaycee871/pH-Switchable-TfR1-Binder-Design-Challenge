#!/usr/bin/env python3
"""Hosted RFdiffusion smoke test for Challenge 04 TfR1.

Attempt 002 shifts the design pressure away from the human-specific
210-212 hotspot and toward a conserved apical subcluster.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "structures/11QC.pdb"
OUTDIR = ROOT / "designs/family_A_conserved/smoke"
URL = "https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate"

# Conserved/exposed apical subcluster selected from the Phase-1 mask.
HOTSPOTS = ["A343", "A344", "A347", "A350", "A374"]
CONTIGS = "A122-760/0 70-90"
OUTPUT_NAME = "rfdiffusion_backbone_002.pdb"

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


def post_json(url: str, payload: dict, api_key: str) -> dict:
    body = json.dumps(payload).encode("utf-8")
    req = Request(
        url,
        data=body,
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
        raise RuntimeError("Expected binder + target chains in RFdiffusion output")

    ordered = sorted(residues, key=lambda c: len(residues[c]))
    binder_chain = ordered[0]
    target_chain = ordered[-1]
    target_first = min(residues[target_chain])

    binder_xyz = [xyz for _, xyz in atoms[binder_chain]]
    contacts = []
    min_dist_by_original = {}

    for out_resid, xyz in atoms[target_chain]:
        original = out_resid - target_first + 122
        best = min(distance(xyz, bxyz) for bxyz in binder_xyz)
        if original not in min_dist_by_original or best < min_dist_by_original[original]:
            min_dist_by_original[original] = best

    for original, best in sorted(min_dist_by_original.items()):
        if best <= cutoff:
            contacts.append(original)

    tier_a_contacts = sorted(set(contacts) & TIER_A)
    risk_contacts = sorted(set(contacts) & RISK)

    return {
        "binder_chain": binder_chain,
        "binder_length": len(residues[binder_chain]),
        "target_chain": target_chain,
        "target_length": len(residues[target_chain]),
        "target_first_output_residue": target_first,
        "contact_cutoff_A": cutoff,
        "human_equivalent_contacts": contacts,
        "tier_A_contacts": tier_a_contacts,
        "risk_contacts": risk_contacts,
        "qc": {
            "has_tier_A_contact": bool(tier_a_contacts),
            "risk_contact_count": len(risk_contacts),
            "advance_to_sequence_design": bool(tier_a_contacts) and len(risk_contacts) <= 1,
        },
    }


def summarize_pdb(pdb_text: str) -> dict:
    _, residues = parse_atoms(pdb_text)
    return {
        chain: {
            "residue_count": len(ids),
            "first_residue": min(ids) if ids else None,
            "last_residue": max(ids) if ids else None,
        }
        for chain, ids in residues.items()
    }


def main() -> None:
    api_key = os.getenv("NVIDIA_API_KEY") or os.getenv("NGC_API_KEY")
    if not api_key:
        raise RuntimeError("Missing NVIDIA_API_KEY / NGC_API_KEY")

    if not INPUT.exists():
        raise RuntimeError(
            "structures/11QC.pdb not found. Run scripts/prepare_targets.py first."
        )

    target_pdb = extract_chain_a(INPUT.read_text(encoding="utf-8"))
    payload = {
        "input_pdb": target_pdb,
        "contigs": CONTIGS,
        "hotspot_res": HOTSPOTS,
        "diffusion_steps": 50,
    }

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "request_metadata_002.json").write_text(
        json.dumps(
            {
                "endpoint": URL,
                "contigs": CONTIGS,
                "hotspot_res": HOTSPOTS,
                "diffusion_steps": 50,
                "target": "11QC chain A",
                "family": "A_conserved_apical",
                "purpose": "cross-species-focused smoke test 002",
                "design_rationale": "shift away from human 210-212 risk patch",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    result = post_json(URL, payload, api_key)
    output_pdb = result.get("output_pdb")
    if not output_pdb:
        raise RuntimeError(
            f"No output_pdb returned. Response keys: {sorted(result.keys())}"
        )

    out_pdb = OUTDIR / OUTPUT_NAME
    out_pdb.write_text(output_pdb, encoding="utf-8")

    contact_analysis = analyze_contacts(output_pdb)
    summary = {
        "elapsed_ms": result.get("elapsed_ms"),
        "output_pdb": str(out_pdb.relative_to(ROOT)),
        "output_structure_summary": summarize_pdb(output_pdb),
        "hotspot_res": HOTSPOTS,
        "contigs": CONTIGS,
        "contact_analysis": contact_analysis,
    }
    (OUTDIR / "result_summary_002.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("RFdiffusion smoke test 002 complete")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
