#!/usr/bin/env python3
"""Position-biased ProteinMPNN exploration on the retained TfR1 backbone.

Only binder-side residues that geometrically contact the RFdiffusion apical
interface are biased. The remainder of the 87-aa scaffold is left to
ProteinMPNN.

ProteinMPNN alphabet order for bias_by_res_jsonl:
ACDEFGHIKLMNPQRSTVWYX
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
DEFAULT_BACKBONE = ROOT / "recovered/seeded_selection/rfdiffusion_seed_20261029.pdb"
OUTDIR = ROOT / "designs/family_A_conserved/proteinmpnn_interface_biased_seed29"

ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
LENGTH = 87

# RFdiffusion geometry, 5-6 Å binder-side interface positions.
CORE_G347 = {14}
CORE_209_210 = {67, 70, 71}
SHELL = {15, 17, 18, 72, 74}


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


def blank_row():
    return [0.0] * len(ALPHABET)


def add_bias(row, residues, value):
    for aa in residues:
        row[ALPHABET.index(aa)] += value


def build_bias_matrix():
    matrix = [blank_row() for _ in range(LENGTH)]

    # Position 14 faces the conserved G347 neighborhood. Encourage packing
    # residues but keep polar/H-bond options and histidine available.
    for pos in CORE_G347:
        row = matrix[pos - 1]
        add_bias(row, "FWYLIVM", 1.20)
        add_bias(row, "HQNST", 0.60)
        add_bias(row, "DEKR", -0.45)

    # Positions facing human L209/V210 and the mouse-shifted chemistry should
    # favor adaptable polar/aromatic chemistry over a hard hydrophobic anchor.
    for pos in CORE_209_210:
        row = matrix[pos - 1]
        add_bias(row, "HQNSTY", 1.25)
        add_bias(row, "FWLIVM", 0.25)
        add_bias(row, "EK", -0.55)

    # Nearby shell: mild encouragement only, preserving model freedom.
    for pos in SHELL:
        row = matrix[pos - 1]
        add_bias(row, "HQNSTYFWLIVM", 0.45)
        add_bias(row, "EK", -0.20)

    return matrix


def main() -> None:
    api_key = os.getenv("NVIDIA_API_KEY") or os.getenv("NGC_API_KEY")
    if not api_key:
        raise RuntimeError("Missing NVIDIA_API_KEY / NGC_API_KEY")

    backbone = Path(os.getenv("BACKBONE_PDB", str(DEFAULT_BACKBONE)))
    if not backbone.exists():
        raise RuntimeError(f"Backbone not found: {backbone}")

    matrix = build_bias_matrix()
    bias_jsonl = json.dumps({"input": {"A": matrix}})

    payload = {
        "input_pdb": backbone.read_text(encoding="utf-8"),
        "input_pdb_chains": ["A"],
        "num_seq_per_target": 6,
        "sampling_temp": [0.15, 0.25],
        "use_soluble_model": False,
        "ca_only": False,
        "omit_AAs": ["C"],
        "bias_by_res_jsonl": bias_jsonl,
        "random_seed": 20261030,
    }

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "bias_spec.json").write_text(
        json.dumps(
            {
                "alphabet": ALPHABET,
                "core_G347": sorted(CORE_G347),
                "core_209_210": sorted(CORE_209_210),
                "shell": sorted(SHELL),
                "matrix": matrix,
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
        score_match = re.search(r"(?:^|, )score=([0-9.]+)", header)
        global_match = re.search(r"global_score=([0-9.]+)", header)
        temp_match = re.search(r"T=([0-9.]+)", header)
        sample_match = re.search(r"sample=(\\d+)", header)
        interface = "".join(binder[p-1] for p in sorted(CORE_G347 | CORE_209_210 | SHELL))

        rows.append({
            "candidate_id": f"TFR1_A29_BIAS_{i:02d}",
            "binder_sequence": binder,
            "binder_length": len(binder),
            "temperature": temp_match.group(1) if temp_match else "",
            "sample": sample_match.group(1) if sample_match else "",
            "mpnn_score": score_match.group(1) if score_match else "",
            "global_score": global_match.group(1) if global_match else "",
            "histidine_count": binder.count("H"),
            "interface_signature": interface,
            "interface_histidines": sum(binder[p-1] == "H" for p in CORE_G347 | CORE_209_210 | SHELL),
            "acidic_count": binder.count("D") + binder.count("E"),
            "basic_count": binder.count("K") + binder.count("R"),
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
        "core_G347": sorted(CORE_G347),
        "core_209_210": sorted(CORE_209_210),
        "shell": sorted(SHELL),
        "histidine_counts": [r["histidine_count"] for r in rows],
        "interface_histidine_counts": [r["interface_histidines"] for r in rows],
    }
    (OUTDIR / "result_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print("Interface-biased ProteinMPNN design complete")
    print(json.dumps(summary, indent=2))
    for row in rows:
        print(
            row["candidate_id"],
            row["mpnn_score"],
            row["global_score"],
            row["histidine_count"],
            row["interface_histidines"],
            row["interface_signature"],
            row["binder_sequence"],
        )


if __name__ == "__main__":
    main()
