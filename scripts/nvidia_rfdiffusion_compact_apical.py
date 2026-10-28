#!/usr/bin/env python3
"""Generate a compact de novo TfR1 apical binder backbone with RFdiffusion.

This second-generation campaign addresses the failure mode of the first
87-residue backbone family:
- shorter binder (55-70 aa)
- hotspot focus restricted to the conserved 343-350 apical subcluster
- backbone-level interface and compactness gates before sequence design
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
OUTDIR = ROOT / "designs/family_A_conserved/compact_apical_backbones"
URL = "https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate"

HOTSPOTS = ["A343", "A344", "A347", "A350"]
CONTIGS = "A122-760/0 55-70"
SEEDS = [20261031, 20261101, 20261102]

PREFERRED = {343, 344, 347, 350}
SUBCLUSTER = set(range(343, 351))
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


def parse_structure(pdb_text: str):
    atoms = {}
    residues = {}
    ca = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM  ") or len(line) < 54:
            continue
        chain = line[21].strip() or "_"
        atom_name = line[12:16].strip()
        try:
            resid = int(line[22:26])
            xyz = tuple(float(line[i:j]) for i, j in ((30, 38), (38, 46), (46, 54)))
        except ValueError:
            continue
        atoms.setdefault(chain, []).append((resid, atom_name, xyz))
        residues.setdefault(chain, set()).add(resid)
        if atom_name == "CA":
            ca.setdefault(chain, {})[resid] = xyz
    return atoms, residues, ca


def dist(a, b):
    return math.sqrt(sum((x-y)**2 for x, y in zip(a, b)))


def compactness(ca_xyz):
    pts = list(ca_xyz)
    n = len(pts)
    centroid = tuple(sum(p[i] for p in pts)/n for i in range(3))
    rg = math.sqrt(sum(dist(p, centroid)**2 for p in pts)/n)
    diameter = 0.0
    for i in range(n):
        for j in range(i+1, n):
            diameter = max(diameter, dist(pts[i], pts[j]))
    end_to_end = dist(pts[0], pts[-1]) if n >= 2 else 0.0
    return {
        "radius_of_gyration_A": round(rg, 3),
        "max_CA_span_A": round(diameter, 3),
        "end_to_end_A": round(end_to_end, 3),
    }


def analyze(pdb_text: str, cutoff=5.0):
    atoms, residues, ca = parse_structure(pdb_text)
    if len(residues) < 2:
        raise RuntimeError("Expected binder + target chains")

    ordered = sorted(residues, key=lambda c: len(residues[c]))
    binder_chain = ordered[0]
    target_chain = ordered[-1]
    target_first = min(residues[target_chain])

    binder_xyz = [xyz for _,_,xyz in atoms[binder_chain]]
    contacts = set()
    for out_resid, _, txyz in atoms[target_chain]:
        if any(dist(txyz, bxyz) <= cutoff for bxyz in binder_xyz):
            # Target input begins at human UniProt 122.
            contacts.add(out_resid - target_first + 122)

    contacts = sorted(contacts)
    subcluster = sorted(set(contacts) & SUBCLUSTER)
    preferred = sorted(set(contacts) & PREFERRED)
    risk = sorted(set(contacts) & RISK)
    geom = compactness([ca[binder_chain][r] for r in sorted(ca[binder_chain])])

    scientific_pass = (
        55 <= len(residues[binder_chain]) <= 70
        and len(subcluster) >= 2
        and len(preferred) >= 1
        and len(risk) <= 1
        and geom["radius_of_gyration_A"] <= 17.0
        and geom["max_CA_span_A"] <= 55.0
    )

    return {
        "binder_chain": binder_chain,
        "binder_length": len(residues[binder_chain]),
        "target_chain": target_chain,
        "target_length": len(residues[target_chain]),
        "contact_cutoff_A": cutoff,
        "human_target_contacts": contacts,
        "subcluster_343_350_contacts": subcluster,
        "preferred_contacts": preferred,
        "risk_contacts": risk,
        "compactness": geom,
        "scientific_pass": scientific_pass,
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
            raise RuntimeError(f"No output_pdb for seed {seed}")

        path = OUTDIR / f"compact_apical_seed_{seed}.pdb"
        path.write_text(output_pdb, encoding="utf-8")
        qc = analyze(output_pdb)
        record = {
            "seed": seed,
            "elapsed_ms": result.get("elapsed_ms"),
            "pdb": str(path.relative_to(ROOT)),
            "qc": qc,
        }
        attempts.append(record)
        print(json.dumps(record, indent=2))

        if qc["scientific_pass"]:
            selected = record
            shutil.copyfile(path, OUTDIR / "selected_backbone.pdb")
            break

    summary = {
        "hotspots": HOTSPOTS,
        "contigs": CONTIGS,
        "seeds": SEEDS,
        "gate": {
            "binder_length": "55-70",
            "min_subcluster_contacts": 2,
            "min_preferred_contacts": 1,
            "max_risk_contacts": 1,
            "max_radius_of_gyration_A": 17.0,
            "max_CA_span_A": 55.0,
        },
        "attempts": attempts,
        "selected": selected,
    }
    (OUTDIR / "selection_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    if selected is None:
        raise SystemExit("No compact apical backbone passed the scientific gate")

    print("Selected compact apical backbone")
    print(json.dumps(selected, indent=2))


if __name__ == "__main__":
    main()
