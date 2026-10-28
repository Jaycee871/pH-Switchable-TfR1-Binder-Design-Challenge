#!/usr/bin/env python3
"""ProteinMPNN sequence design for the compact TfR1 apical backbone.

Input is the selected 57-aa RFdiffusion seed 20261102 backbone.
Only binder chain A is redesigned; TfR1 chain B remains fixed.
"""

from __future__ import annotations

import csv
import json
import os
import re
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
URL = "https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict"
DEFAULT_BACKBONE = ROOT / "recovered/selected_backbone.pdb"
OUTDIR = ROOT / "designs/family_A_conserved/compact_seed20261102_proteinmpnn"

# Binder-side geometry from the selected RFdiffusion backbone.
CORE_INTERFACE = [1, 2, 5]
SECONDARY_CONTACT = [41, 45]
SHELL = [3, 4, 6, 44]


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
        "num_seq_per_target": 6,
        "sampling_temp": [0.10, 0.20],
        "use_soluble_model": False,
        "ca_only": False,
        "omit_AAs": ["C"],
        "random_seed": 20261102,
    }

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "request_metadata.json").write_text(
        json.dumps(
            {
                "backbone": str(backbone),
                "designed_chain": "A",
                "fixed_target_chain": "B",
                "binder_length": 57,
                "num_seq_per_target": 6,
                "sampling_temp": [0.10, 0.20],
                "use_soluble_model": False,
                "omit_AAs": ["C"],
                "random_seed": 20261102,
                "core_interface_positions": CORE_INTERFACE,
                "secondary_contact_positions": SECONDARY_CONTACT,
                "shell_positions": SHELL,
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    result = post_json(payload, api_key)
    mfasta = result.get("mfasta", "")
    if not mfasta.strip():
        raise RuntimeError(f"No mfasta returned; keys={sorted(result.keys())}")
    (OUTDIR / "designed_sequences.fasta").write_text(mfasta, encoding="utf-8")

    records = parse_fasta(mfasta)
    designed = [
        (h, s) for h, s in records
        if not h.lower().startswith("input")
        and "native" not in h.lower()
        and not h.lower().startswith("wt")
    ]

    rows = []
    for i, (header, sequence) in enumerate(designed, start=1):
        binder = [p for p in sequence.split("/") if p][0]
        if len(binder) != 57:
            raise RuntimeError(f"Unexpected binder length {len(binder)}")
        score_match = re.search(r"(?:^|, )score=([0-9.]+)", header)
        global_match = re.search(r"global_score=([0-9.]+)", header)
        temp_match = re.search(r"T=([0-9.]+)", header)
        sample_match = re.search(r"sample=(\\d+)", header)

        core_sig = "".join(binder[p - 1] for p in CORE_INTERFACE)
        contact_sig = "".join(
            binder[p - 1] for p in CORE_INTERFACE + SECONDARY_CONTACT
        )

        rows.append({
            "candidate_id": f"TFR1_C57_{i:02d}",
            "binder_sequence": binder,
            "binder_length": len(binder),
            "temperature": temp_match.group(1) if temp_match else "",
            "sample": sample_match.group(1) if sample_match else "",
            "mpnn_score": score_match.group(1) if score_match else "",
            "global_score": global_match.group(1) if global_match else "",
            "histidine_count": binder.count("H"),
            "aromatic_count": sum(binder.count(x) for x in "FWY"),
            "hydrophobic_count": sum(binder.count(x) for x in "AILMVFWY"),
            "acidic_count": binder.count("D") + binder.count("E"),
            "basic_count": binder.count("K") + binder.count("R"),
            "core_interface_signature": core_sig,
            "contact_signature": contact_sig,
            "header": header,
        })

    if not rows:
        raise RuntimeError("No designed sequences parsed")

    with (OUTDIR / "candidates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "candidate_count": len(rows),
        "candidate_lengths": sorted(set(r["binder_length"] for r in rows)),
        "input_rows_excluded": len(records) - len(designed),
        "core_interface_positions": CORE_INTERFACE,
        "secondary_contact_positions": SECONDARY_CONTACT,
        "shell_positions": SHELL,
        "histidine_counts": [r["histidine_count"] for r in rows],
        "aromatic_counts": [r["aromatic_count"] for r in rows],
        "hydrophobic_counts": [r["hydrophobic_count"] for r in rows],
    }
    (OUTDIR / "result_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print("Compact-backbone ProteinMPNN design complete")
    print(json.dumps(summary, indent=2))
    for row in rows:
        print(
            row["candidate_id"],
            "score="+str(row["mpnn_score"]),
            "H="+str(row["histidine_count"]),
            "arom="+str(row["aromatic_count"]),
            "hydro="+str(row["hydrophobic_count"]),
            "core="+row["core_interface_signature"],
            "contacts="+row["contact_signature"],
            row["binder_sequence"],
        )


if __name__ == "__main__":
    main()
