#!/usr/bin/env python3
"""Histidine-enriched ProteinMPNN design on the compact TfR1 apical backbone.

Phase 6. Input is the gated 57-aa RFdiffusion seed 20261102 backbone, which is
the only backbone in this project to have passed its own scientific gate. Only
binder chain A is redesigned; TfR1 chain B remains fixed.

Motivation
----------
The competition ranks pH-selective release first. Histidine protonation near
pH 6.0 is the proposed switch mechanism, but all 60+ sequences designed so far
contain zero histidines. Phase 3D applied a bias of only +0.60 to +1.25 to H at
interface positions and still returned H=0 for every one of its 12 designs, so
that magnitude is empirically insufficient to overcome this backbone's strong
preference for E/K.

This phase therefore applies a substantially stronger, explicitly positive H
bias at the binder-side interface positions, together with a mild penalty on
the charged residues that currently saturate those positions.

Interface definition
--------------------
Binder-side positions are derived from the backbone geometry itself (minimum
heavy-atom distance from chain A residue to chain B):

  CORE  (<= 5 A): 1, 2, 5, 41, 45
  RIM   (5-6 A) : 3, 4, 6, 44
  OUTER (6-8 A) : 37, 38, 40, 42, 46, 48, 49, 52, 55, 56

CORE and RIM reproduce exactly the positions used in Phase 3E.

ProteinMPNN alphabet order for bias_by_res_jsonl:
ACDEFGHIKLMNPQRSTVWYX
"""

from __future__ import annotations

import csv
import json
import os
import re
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
URL = "https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict"
DEFAULT_BACKBONE = ROOT / "structures/rfdiffusion/compact_apical/selected_backbone.pdb"
OUTDIR = ROOT / "designs/family_C_pH_switch/his_enriched_seed20261102"

ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
LENGTH = 57

CORE = [1, 2, 5, 41, 45]
RIM = [3, 4, 6, 44]
OUTER = [37, 38, 40, 42, 46, 48, 49, 52, 55, 56]

# Bias magnitudes. ProteinMPNN adds these to the per-position logits, so a
# bias of b multiplies the odds of that residue by roughly exp(b).
#
# Phase 3D demonstrated that +1.25 is not enough (H=0 in 12/12 designs).
# +4.0 is ~e^4 = 55x odds, chosen as the smallest round value comfortably
# above the magnitude already shown to fail. It is a bias and not a hard
# constraint: ProteinMPNN can still decline H where it is structurally
# untenable, which is the desired behaviour.
H_BIAS_CORE = 4.0
H_BIAS_RIM = 3.0
H_BIAS_OUTER = 1.5
CHARGE_PENALTY_CORE = -1.0
CHARGE_PENALTY_RIM = -0.6

MIN_HISTIDINES = 3
# Hard rejection only for genuinely degenerate composition. This backbone
# family is natively E/K-rich (existing C57 designs run 39-46% E), so a
# threshold near 0.40 would reject every sequence for a reason unrelated to
# histidine enrichment. Polyampholyte character is recorded as a non-fatal
# note here and is properly the subject of the Phase 9 developability filter.
MAX_SINGLE_AA_FRACTION = 0.50
MAX_HOMOPOLYMER_RUN = 5
SUBMISSION_MIN_AA = 10
SUBMISSION_MAX_AA = 250


def post_json(payload: dict, api_key: str, max_attempts: int = 4) -> dict:
    """POST with the Phase 4E rate-limit handling: honour Retry-After, back off."""
    for attempt in range(1, max_attempts + 1):
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
            with urlopen(req, timeout=1200) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")
            if exc.code == 429 and attempt < max_attempts:
                retry_after = exc.headers.get("Retry-After")
                sleep_s = (
                    int(retry_after)
                    if retry_after and retry_after.isdigit()
                    else 20 * attempt
                )
                print(f"ProteinMPNN 429 on attempt {attempt}; retrying after {sleep_s}s")
                time.sleep(sleep_s)
                continue
            raise RuntimeError(f"ProteinMPNN HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("ProteinMPNN retry loop exhausted")


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

    # Core contacts: strongly encourage histidine, and relieve the E/K/D/R
    # saturation that currently occupies these positions.
    for pos in CORE:
        row = matrix[pos - 1]
        add_bias(row, "H", H_BIAS_CORE)
        add_bias(row, "QNSTY", 0.30)
        add_bias(row, "DEKR", CHARGE_PENALTY_CORE)

    # Rim: encourage histidine somewhat less, so the network can extend
    # outward from the core without the whole rim collapsing to H.
    for pos in RIM:
        row = matrix[pos - 1]
        add_bias(row, "H", H_BIAS_RIM)
        add_bias(row, "QNSTY", 0.20)
        add_bias(row, "DEKR", CHARGE_PENALTY_RIM)

    # Outer shell: mild only. A supporting histidine here is useful, but these
    # positions are 6-8 A from the target and should not be forced.
    for pos in OUTER:
        row = matrix[pos - 1]
        add_bias(row, "H", H_BIAS_OUTER)

    return matrix


def longest_run(seq: str) -> int:
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best


def composition_notes(binder: str) -> list[str]:
    """Non-fatal observations, reported but not used to reject."""
    notes = []
    acidic = binder.count("D") + binder.count("E")
    basic = binder.count("K") + binder.count("R")
    if (acidic + basic) / len(binder) > 0.55:
        notes.append(f"high_charge_fraction_{round((acidic + basic) / len(binder), 2)}")
    counts = {aa: binder.count(aa) for aa in set(binder)}
    top_aa = max(counts, key=counts.get)
    if counts[top_aa] / len(binder) > 0.35:
        notes.append(f"enriched_{top_aa}_{round(counts[top_aa] / len(binder), 2)}")
    return notes


def composition_flags(binder: str) -> list[str]:
    """Fatal flags: genuinely degenerate or non-submittable sequences."""
    flags = []
    if not (SUBMISSION_MIN_AA <= len(binder) <= SUBMISSION_MAX_AA):
        flags.append(f"length_outside_{SUBMISSION_MIN_AA}_{SUBMISSION_MAX_AA}")
    counts = {aa: binder.count(aa) for aa in set(binder)}
    top_aa = max(counts, key=counts.get)
    if counts[top_aa] / len(binder) > MAX_SINGLE_AA_FRACTION:
        flags.append(f"dominant_residue_{top_aa}_{counts[top_aa]}")
    run = longest_run(binder)
    if run > MAX_HOMOPOLYMER_RUN:
        flags.append(f"homopolymer_run_{run}")
    if "C" in binder:
        flags.append("free_cysteine")
    return flags


def main() -> None:
    api_key = os.getenv("NVIDIA_API_KEY") or os.getenv("NGC_API_KEY")
    if not api_key:
        raise RuntimeError("Missing NVIDIA_API_KEY / NGC_API_KEY")

    backbone = Path(os.getenv("BACKBONE_PDB", str(DEFAULT_BACKBONE)))
    if not backbone.exists():
        raise RuntimeError(f"Backbone not found: {backbone}")

    matrix = build_bias_matrix()
    bias_jsonl = json.dumps({"input": {"A": matrix}})

    num_seq = 12
    temps = [0.10, 0.20]
    seed = 20261106

    payload = {
        "input_pdb": backbone.read_text(encoding="utf-8"),
        "input_pdb_chains": ["A"],
        "num_seq_per_target": num_seq,
        "sampling_temp": temps,
        "use_soluble_model": False,
        "ca_only": False,
        "omit_AAs": ["C"],
        "bias_by_res_jsonl": bias_jsonl,
        "random_seed": seed,
    }

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "request_metadata.json").write_text(
        json.dumps(
            {
                "phase": "6",
                "endpoint": URL,
                "backbone": str(backbone),
                "backbone_source": "structures/rfdiffusion/compact_apical/selected_backbone.pdb (RFdiffusion seed 20261102, run 40685364707, scientific_pass=true)",
                "designed_chain": "A",
                "fixed_target_chain": "B",
                "binder_length": LENGTH,
                "num_seq_per_target": num_seq,
                "sampling_temp": temps,
                "use_soluble_model": False,
                "omit_AAs": ["C"],
                "random_seed": seed,
                "core_positions": CORE,
                "rim_positions": RIM,
                "outer_positions": OUTER,
                "histidine_bias": {
                    "core": H_BIAS_CORE,
                    "rim": H_BIAS_RIM,
                    "outer": H_BIAS_OUTER,
                    "charge_penalty_core": CHARGE_PENALTY_CORE,
                    "charge_penalty_rim": CHARGE_PENALTY_RIM,
                    "rationale": (
                        "Phase 3D used +0.60 to +1.25 for H and produced zero "
                        "histidines in 12/12 designs. Bias adds to logits, so "
                        "+4.0 is ~exp(4)=55x odds; chosen as the smallest round "
                        "magnitude well above the value already shown to fail, "
                        "while remaining a bias rather than a hard constraint."
                    ),
                },
                "filter": {
                    "min_histidines": MIN_HISTIDINES,
                    "max_single_aa_fraction": MAX_SINGLE_AA_FRACTION,
                    "max_homopolymer_run": MAX_HOMOPOLYMER_RUN,
                    "submission_length_range": [SUBMISSION_MIN_AA, SUBMISSION_MAX_AA],
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (OUTDIR / "bias_spec.json").write_text(
        json.dumps(
            {
                "alphabet": ALPHABET,
                "length": LENGTH,
                "core": CORE,
                "rim": RIM,
                "outer": OUTER,
                "histidine_bias_core": H_BIAS_CORE,
                "histidine_bias_rim": H_BIAS_RIM,
                "histidine_bias_outer": H_BIAS_OUTER,
                "matrix": matrix,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    result = post_json(payload, api_key)
    mfasta = result.get("mfasta", "")
    if not mfasta.strip():
        raise RuntimeError(f"No mfasta returned; keys={sorted(result.keys())}")
    (OUTDIR / "designed_sequences.fasta").write_text(mfasta, encoding="utf-8")

    records = parse_fasta(mfasta)
    designed = [
        (h, s)
        for h, s in records
        if not h.lower().startswith("input")
        and "native" not in h.lower()
        and not h.lower().startswith("wt")
    ]

    interface_positions = sorted(CORE + RIM)
    rows = []
    for i, (header, sequence) in enumerate(designed, start=1):
        binder = [p for p in sequence.split("/") if p][0]
        if len(binder) != LENGTH:
            raise RuntimeError(f"Unexpected binder length {len(binder)}")
        score_match = re.search(r"(?:^|, )score=([0-9.]+)", header)
        global_match = re.search(r"global_score=([0-9.]+)", header)
        temp_match = re.search(r"T=([0-9.]+)", header)
        sample_match = re.search(r"sample=(\d+)", header)

        his_positions = [p for p in range(1, LENGTH + 1) if binder[p - 1] == "H"]
        core_his = [p for p in CORE if binder[p - 1] == "H"]
        rim_his = [p for p in RIM if binder[p - 1] == "H"]
        outer_his = [p for p in OUTER if binder[p - 1] == "H"]
        off_his = [p for p in his_positions if p not in CORE + RIM + OUTER]

        acidic = binder.count("D") + binder.count("E")
        basic = binder.count("K") + binder.count("R")
        flags = composition_flags(binder)
        notes = composition_notes(binder)
        passes = binder.count("H") >= MIN_HISTIDINES and not flags

        rows.append(
            {
                "candidate_id": f"TFR1_H57_{i:02d}",
                "binder_sequence": binder,
                "binder_length": len(binder),
                "temperature": temp_match.group(1) if temp_match else "",
                "sample": sample_match.group(1) if sample_match else "",
                "mpnn_score": score_match.group(1) if score_match else "",
                "global_score": global_match.group(1) if global_match else "",
                "histidine_count": binder.count("H"),
                "histidine_positions": ";".join(str(p) for p in his_positions),
                "core_histidines": ";".join(str(p) for p in core_his),
                "rim_histidines": ";".join(str(p) for p in rim_his),
                "outer_histidines": ";".join(str(p) for p in outer_his),
                "offsite_histidines": ";".join(str(p) for p in off_his),
                "interface_histidine_count": len(core_his) + len(rim_his),
                "interface_signature": "".join(
                    binder[p - 1] for p in interface_positions
                ),
                "acidic_count": acidic,
                "basic_count": basic,
                "net_charge_pH7": basic - acidic,
                "net_charge_pH6_all_his_protonated": basic - acidic + binder.count("H"),
                "aromatic_count": sum(binder.count(x) for x in "FWY"),
                "hydrophobic_count": sum(binder.count(x) for x in "AILMVFWY"),
                "composition_flags": ";".join(flags),
                "composition_notes": ";".join(notes),
                "passes_filter": passes,
                "header": header,
            }
        )

    if not rows:
        raise RuntimeError("No designed sequences parsed")

    with (OUTDIR / "candidates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    survivors = [r for r in rows if r["passes_filter"]]
    with (OUTDIR / "candidates_passing.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(survivors)

    his_counts = [r["histidine_count"] for r in rows]
    summary = {
        "phase": "6",
        "candidate_count": len(rows),
        "candidate_lengths": sorted(set(r["binder_length"] for r in rows)),
        "input_rows_excluded": len(records) - len(designed),
        "core_positions": CORE,
        "rim_positions": RIM,
        "outer_positions": OUTER,
        "histidine_bias": {
            "core": H_BIAS_CORE,
            "rim": H_BIAS_RIM,
            "outer": H_BIAS_OUTER,
        },
        "histidine_counts": his_counts,
        "histidine_count_min": min(his_counts),
        "histidine_count_max": max(his_counts),
        "histidine_count_mean": round(sum(his_counts) / len(his_counts), 2),
        "designs_with_zero_histidine": sum(1 for c in his_counts if c == 0),
        "interface_histidine_counts": [r["interface_histidine_count"] for r in rows],
        "min_histidines_required": MIN_HISTIDINES,
        "passing_count": len(survivors),
        "passing_ids": [r["candidate_id"] for r in survivors],
        "rejected": [
            {
                "candidate_id": r["candidate_id"],
                "histidine_count": r["histidine_count"],
                "composition_flags": r["composition_flags"],
            }
            for r in rows
            if not r["passes_filter"]
        ],
        "core_position_occupancy": {
            str(p): "".join(r["binder_sequence"][p - 1] for r in rows) for p in CORE
        },
        "rim_position_occupancy": {
            str(p): "".join(r["binder_sequence"][p - 1] for r in rows) for p in RIM
        },
    }
    (OUTDIR / "result_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print("His-enriched ProteinMPNN design complete")
    print(json.dumps(summary, indent=2))
    for row in rows:
        print(
            row["candidate_id"],
            "score=" + str(row["mpnn_score"]),
            "H=" + str(row["histidine_count"]),
            "ifaceH=" + str(row["interface_histidine_count"]),
            "Hpos=" + (row["histidine_positions"] or "-"),
            "pass=" + str(row["passes_filter"]),
            row["binder_sequence"],
        )

    if not survivors:
        print(
            "WARNING: no design met the His>=%d filter; see result_summary.json"
            % MIN_HISTIDINES
        )


if __name__ == "__main__":
    main()
