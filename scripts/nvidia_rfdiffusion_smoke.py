#!/usr/bin/env python3
"""Single hosted RFdiffusion smoke test for Challenge 04 TfR1.

This intentionally generates only ONE Family-A backbone to validate the
target/hotspot configuration before a larger campaign is launched.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "structures/11QC.pdb"
OUTDIR = ROOT / "designs/family_A_conserved/smoke"
URL = "https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate"

# Conserved / cross-species-friendlier anchor set around the apical patch.
HOTSPOTS = ["A206", "A209", "A214", "A344", "A374"]
CONTIGS = "A122-760/0 70-90"


def extract_chain_a(pdb_text: str) -> str:
    out = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            if len(line) >= 22 and line[21] == "A":
                # Keep protein atoms and any target-associated HETATM on chain A.
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


def summarize_pdb(pdb_text: str) -> dict:
    residues = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM  ") or len(line) < 27:
            continue
        chain = line[21].strip() or "_"
        try:
            resid = int(line[22:26])
        except ValueError:
            continue
        residues.setdefault(chain, set()).add(resid)
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
    (OUTDIR / "request_metadata.json").write_text(
        json.dumps(
            {
                "endpoint": URL,
                "contigs": CONTIGS,
                "hotspot_res": HOTSPOTS,
                "diffusion_steps": 50,
                "target": "11QC chain A",
                "family": "A_conserved_apical",
                "purpose": "single-backbone smoke test",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    result = post_json(URL, payload, api_key)
    output_pdb = result.get("output_pdb")
    if not output_pdb:
        safe_keys = sorted(result.keys())
        raise RuntimeError(f"No output_pdb returned. Response keys: {safe_keys}")

    out_pdb = OUTDIR / "rfdiffusion_backbone_001.pdb"
    out_pdb.write_text(output_pdb, encoding="utf-8")

    summary = {
        "elapsed_ms": result.get("elapsed_ms"),
        "output_pdb": str(out_pdb.relative_to(ROOT)),
        "output_structure_summary": summarize_pdb(output_pdb),
        "hotspot_res": HOTSPOTS,
        "contigs": CONTIGS,
    }
    (OUTDIR / "result_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("RFdiffusion smoke test complete")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
