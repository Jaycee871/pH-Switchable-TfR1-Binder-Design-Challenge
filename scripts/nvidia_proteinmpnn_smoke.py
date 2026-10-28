#!/usr/bin/env python3
"""ProteinMPNN sequence-design stage for an RFdiffusion TfR1 binder backbone.

Designs binder chain A only; target chain B is kept fixed.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
URL = "https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict"
DEFAULT_BACKBONE = ROOT / "designs/family_A_conserved/seeded_selection/selected_backbone.pdb"
OUTDIR = ROOT / "designs/family_A_conserved/proteinmpnn_seeded"


def post_json(url: str, payload: dict, api_key: str) -> dict:
    req = Request(
        url,
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
        raise RuntimeError(f"ProteinMPNN HTTP {exc.code}: {detail}") from exc


def parse_fasta(text: str):
    records = []
    header = None
    seq = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header = line[1:]
            seq = []
        else:
            seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def main() -> None:
    api_key = os.getenv("NVIDIA_API_KEY") or os.getenv("NGC_API_KEY")
    if not api_key:
        raise RuntimeError("Missing NVIDIA_API_KEY / NGC_API_KEY")

    backbone = Path(os.getenv("BACKBONE_PDB", str(DEFAULT_BACKBONE)))
    if not backbone.exists():
        raise RuntimeError(f"Backbone not found: {backbone}")

    payload = {
        "input_pdb": backbone.read_text(encoding="utf-8"),
        "input_pdb_chains": ["A"],
        "num_seq_per_target": 8,
        "sampling_temp": [0.1, 0.2],
        "use_soluble_model": True,
        "ca_only": False,
        "omit_AAs": ["C"],
        "random_seed": 20261028,
    }

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "request_metadata.json").write_text(
        json.dumps(
            {
                "endpoint": URL,
                "backbone": str(backbone),
                "designed_chain": "A",
                "fixed_target_chain": "B",
                "num_seq_per_target": 8,
                "sampling_temp": [0.1, 0.2],
                "use_soluble_model": True,
                "omit_AAs": ["C"],
                "random_seed": 20261028,
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    result = post_json(URL, payload, api_key)
    mfasta = result.get("mfasta", "")
    if not mfasta.strip():
        raise RuntimeError(f"No mfasta returned. Response keys: {sorted(result.keys())}")

    (OUTDIR / "designed_sequences.fasta").write_text(mfasta, encoding="utf-8")

    records = parse_fasta(mfasta)
    designed = [
        (h, s)
        for h, s in records
        if not h.lower().startswith("input")
        and "native" not in h.lower()
        and not h.lower().startswith("wt")
    ]
    scores = result.get("scores", [])

    rows = []
    for i, (header, sequence) in enumerate(designed, start=1):
        # ProteinMPNN may separate chains with '/' in FASTA output.
        parts = [p for p in sequence.split("/") if p]
        binder = parts[0] if parts else sequence
        score_match = re.search(r"(?:^|, )score=([0-9.]+)", header)
        global_match = re.search(r"global_score=([0-9.]+)", header)
        temp_match = re.search(r"T=([0-9.]+)", header)
        sample_match = re.search(r"sample=(\\d+)", header)
        rows.append(
            {
                "candidate_id": f"TFR1_A29_MPNN_{i:02d}",
                "binder_sequence": binder,
                "binder_length": len(binder),
                "temperature": temp_match.group(1) if temp_match else "",
                "sample": sample_match.group(1) if sample_match else "",
                "mpnn_score": score_match.group(1) if score_match else "",
                "global_score": global_match.group(1) if global_match else "",
                "header": header,
            }
        )

    if not rows:
        raise RuntimeError("No designed sequences parsed from ProteinMPNN output")

    with (OUTDIR / "candidates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "backbone": str(backbone.relative_to(ROOT)) if backbone.is_relative_to(ROOT) else str(backbone),
        "designed_chain": "A",
        "candidate_count": len(rows),
        "candidate_lengths": sorted(set(r["binder_length"] for r in rows)),
        "scores_returned": len(scores),\n        "input_rows_excluded": len(records) - len(designed),
        "output_fasta": str((OUTDIR / "designed_sequences.fasta").relative_to(ROOT)),
        "output_csv": str((OUTDIR / "candidates.csv").relative_to(ROOT)),
    }
    (OUTDIR / "result_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print("ProteinMPNN sequence design complete")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
