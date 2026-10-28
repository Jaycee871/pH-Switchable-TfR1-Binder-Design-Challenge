#!/usr/bin/env python3
"""Pose-first ProteinMPNN design for Kaggle-selected TfR1 apical backbones.

Only RFdiffusion backbones that passed scripts/select_kaggle_backbones.py are used.
This phase intentionally does NOT force histidine. The immediate question is whether
sequence design can preserve a full-ectodomain human/mouse apical binding pose.
pH-switch engineering is deferred until a pose-valid scaffold exists.
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
SEL_DIR = ROOT / "structures/rfdiffusion/kaggle_pilot_selected"
SUMMARY = SEL_DIR / "selection_summary.json"
OUT = ROOT / "designs/family_A_conserved/kaggle_pose_first"
URL = "https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict"

MAX_BACKBONES = 3
NUM_SEQ_PER_TARGET = 4
TEMPS = [0.10, 0.20]


def post_json(payload: dict, api_key: str, max_attempts: int = 4) -> dict:
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
                sleep_s = int(retry_after) if retry_after and retry_after.isdigit() else 20 * attempt
                print(f"ProteinMPNN 429 attempt {attempt}; sleeping {sleep_s}s", flush=True)
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


def chain_info(pdb: Path):
    order = []
    residues = {}
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        ch = line[21].strip() or "_"
        try:
            r = int(line[22:26])
        except ValueError:
            continue
        if ch not in residues:
            residues[ch] = set()
            order.append(ch)
        residues[ch].add(r)
    if len(residues) < 2:
        raise RuntimeError(f"Expected target+binder chains in {pdb}: {residues}")
    binder = min(residues, key=lambda x: len(residues[x]))
    target = max(residues, key=lambda x: len(residues[x]))
    return order, {k: len(v) for k, v in residues.items()}, binder, target


def binder_from_record(sequence: str, chain_order: list[str], binder_chain: str) -> str:
    parts = [x for x in sequence.split("/") if x]
    if len(parts) == len(chain_order) and binder_chain in chain_order:
        return parts[chain_order.index(binder_chain)]
    if len(parts) == 1:
        return parts[0]
    # Fallback: the binder is the shortest chain.
    return min(parts, key=len)


def header_value(header: str, pattern: str):
    m = re.search(pattern, header)
    return m.group(1).strip() if m else ""


def longest_run(seq: str):
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i-1] else 1
        best = max(best, run)
    return best


def main():
    api_key = (
        os.getenv("NVIDIA_API_KEY")
        or os.getenv("NVIDIA_API_KEY1")
        or os.getenv("NVIDIA_API_TOKEN")
        or os.getenv("NGC_API_KEY")
    )
    if not api_key:
        raise RuntimeError("Missing NVIDIA API credential")
    if not SUMMARY.exists():
        raise SystemExit(f"Selection summary not present: {SUMMARY}")

    selection = json.loads(SUMMARY.read_text(encoding="utf-8"))
    selected = selection.get("selected", [])[:MAX_BACKBONES]
    if not selected:
        raise SystemExit("No strict-gate Kaggle backbone is available for ProteinMPNN.")

    OUT.mkdir(parents=True, exist_ok=True)
    all_rows = []
    backbone_summaries = []

    for bidx, item in enumerate(selected, start=1):
        rel = item["selected_copy"]
        pdb = ROOT / rel
        if not pdb.exists():
            raise RuntimeError(f"Selected backbone copy missing: {pdb}")

        chain_order, counts, binder_chain, target_chain = chain_info(pdb)
        expected_len = counts[binder_chain]
        design_id = f"KB{bidx:02d}_{item['family']}"
        bdir = OUT / design_id
        bdir.mkdir(parents=True, exist_ok=True)

        payload = {
            "input_pdb": pdb.read_text(encoding="utf-8"),
            "input_pdb_chains": [binder_chain],
            "num_seq_per_target": NUM_SEQ_PER_TARGET,
            "sampling_temp": TEMPS,
            "use_soluble_model": False,
            "ca_only": False,
            "omit_AAs": ["C"],
            "random_seed": 20261120 + bidx,
        }
        (bdir / "request_metadata.json").write_text(
            json.dumps({
                "source_backbone": rel,
                "source_family": item["family"],
                "chain_order": chain_order,
                "chain_lengths": counts,
                "designed_binder_chain": binder_chain,
                "fixed_target_chain": target_chain,
                "expected_binder_length": expected_len,
                "preferred_contacts": item.get("preferred_contacts", []),
                "risk_contacts": item.get("risk_contacts", []),
                "num_seq_per_target": NUM_SEQ_PER_TARGET,
                "sampling_temp": TEMPS,
                "design_policy": "pose-first; no forced histidine",
            }, indent=2) + "\n",
            encoding="utf-8",
        )

        result = post_json(payload, api_key)
        mfasta = result.get("mfasta", "")
        if not mfasta.strip():
            raise RuntimeError(f"No mfasta returned for {design_id}; keys={sorted(result)}")
        (bdir / "designed_sequences.fasta").write_text(mfasta, encoding="utf-8")

        records = [
            (h, s) for h, s in parse_fasta(mfasta)
            if not h.lower().startswith("input")
            and "native" not in h.lower()
            and not h.lower().startswith("wt")
        ]

        local = []
        for i, (header, sequence) in enumerate(records, start=1):
            binder = binder_from_record(sequence, chain_order, binder_chain)
            if len(binder) != expected_len:
                raise RuntimeError(
                    f"{design_id} candidate {i}: expected {expected_len} aa, got {len(binder)}"
                )
            acidic = binder.count("D") + binder.count("E")
            basic = binder.count("K") + binder.count("R")
            charge_fraction = (acidic + basic) / len(binder)
            top_count = max(binder.count(aa) for aa in set(binder))
            flags = []
            if charge_fraction > 0.65:
                flags.append("charge_fraction_gt_0.65")
            if top_count / len(binder) > 0.40:
                flags.append("single_residue_fraction_gt_0.40")
            if longest_run(binder) > 4:
                flags.append("homopolymer_run_gt_4")
            if "C" in binder:
                flags.append("free_cysteine")

            row = {
                "candidate_id": f"TFR1_{design_id}_{i:02d}",
                "backbone_id": design_id,
                "source_family": item["family"],
                "source_backbone": rel,
                "binder_chain": binder_chain,
                "binder_length": len(binder),
                "sequence": binder,
                "temperature": header_value(header, r"(?:^|,\s*)T=([^,]+)"),
                "sample": header_value(header, r"(?:^|,\s*)sample=([^,]+)"),
                "mpnn_score": header_value(header, r"(?:^|,\s*)score=([^,]+)"),
                "global_score": header_value(header, r"(?:^|,\s*)global_score=([^,]+)"),
                "histidine_count": binder.count("H"),
                "acidic_count": acidic,
                "basic_count": basic,
                "charge_fraction": round(charge_fraction, 4),
                "aromatic_count": sum(binder.count(x) for x in "FWY"),
                "hydrophobic_count": sum(binder.count(x) for x in "AILMVFWY"),
                "prefilter_flags": ";".join(flags),
                "prefilter_pass": not flags,
                "header": header,
            }
            local.append(row)
            all_rows.append(row)

        with (bdir / "candidates.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=local[0].keys())
            w.writeheader()
            w.writerows(local)

        backbone_summaries.append({
            "backbone_id": design_id,
            "source": rel,
            "candidate_count": len(local),
            "prefilter_pass_count": sum(r["prefilter_pass"] for r in local),
            "histidine_counts": [r["histidine_count"] for r in local],
        })
        time.sleep(10)

    # Select at most two pose-screen candidates per backbone by low MPNN score,
    # after coarse composition filtering. This keeps downstream Boltz-2 calls small.
    chosen = []
    for b in backbone_summaries:
        rows = [r for r in all_rows if r["backbone_id"] == b["backbone_id"] and r["prefilter_pass"]]
        rows.sort(key=lambda r: float(r["mpnn_score"]) if r["mpnn_score"] else 999.0)
        chosen.extend(rows[:2])

    fields = all_rows[0].keys()
    with (OUT / "all_candidates.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(all_rows)
    with (OUT / "pose_screen_candidates.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(chosen)

    (OUT / "result_summary.json").write_text(
        json.dumps({
            "policy": "pose-first",
            "backbones_used": len(backbone_summaries),
            "candidate_count": len(all_rows),
            "pose_screen_candidate_count": len(chosen),
            "backbones": backbone_summaries,
            "pose_screen_ids": [r["candidate_id"] for r in chosen],
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    print((OUT / "result_summary.json").read_text())


if __name__ == "__main__":
    main()
