#!/usr/bin/env python3
"""Pose-first ProteinMPNN design for NVIDIA apical-only RFdiffusion backbones.

Consumes structures/rfdiffusion/nvidia_apical_only/selection_summary.json.
Only strict-gate selected backbones are designed. Histidine is not forced here;
the purpose is to identify sequence/backbone pairs that preserve the intended
apical pose in subsequent full human/mouse validation.
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
SEL = ROOT / "structures/rfdiffusion/nvidia_apical_only/selection_summary.json"
OUT = ROOT / "designs/family_A_conserved/nvidia_apical_pose_first"
URL = "https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict"

MAX_BACKBONES = 4
NUM_SEQ_PER_TARGET = 4
TEMPS = [0.10, 0.20]


def post_json(payload, api_key, max_attempts=4):
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
                ra = exc.headers.get("Retry-After")
                sleep_s = int(ra) if ra and ra.isdigit() else 20 * attempt
                print(f"ProteinMPNN 429 attempt {attempt}; sleep {sleep_s}s", flush=True)
                time.sleep(sleep_s)
                continue
            raise RuntimeError(f"ProteinMPNN HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("ProteinMPNN retries exhausted")


def parse_fasta(text):
    out, header, seq = [], None, []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                out.append((header, "".join(seq)))
            header, seq = line[1:], []
        else:
            seq.append(line)
    if header is not None:
        out.append((header, "".join(seq)))
    return out


def chain_info(path):
    order, residues = [], {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        ch = line[21].strip() or "_"
        try:
            r = int(line[22:26])
        except ValueError:
            continue
        if ch not in residues:
            order.append(ch)
            residues[ch] = set()
        residues[ch].add(r)
    if len(residues) < 2:
        raise RuntimeError(f"Expected target+binder chains in {path}")
    binder = min(residues, key=lambda x: len(residues[x]))
    target = max(residues, key=lambda x: len(residues[x]))
    return order, {k: len(v) for k, v in residues.items()}, binder, target


def extract_binder(sequence, chain_order, binder_chain):
    parts = [x for x in sequence.split("/") if x]
    if len(parts) == len(chain_order) and binder_chain in chain_order:
        return parts[chain_order.index(binder_chain)]
    if len(parts) == 1:
        return parts[0]
    return min(parts, key=len)


def hv(header, regex):
    m = re.search(regex, header)
    return m.group(1).strip() if m else ""


def longest_run(seq):
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
    if not SEL.exists():
        raise SystemExit("NVIDIA apical-only selection summary not present")

    selection = json.loads(SEL.read_text(encoding="utf-8"))
    selected = selection.get("selected", [])[:MAX_BACKBONES]
    if not selected:
        raise SystemExit("No NVIDIA apical-only backbone passed the strict gate")

    OUT.mkdir(parents=True, exist_ok=True)
    rows_all = []
    backbone_meta = []

    for idx, item in enumerate(selected, start=1):
        pdb = ROOT / item["selected_copy"]
        order, lengths, binder_chain, target_chain = chain_info(pdb)
        expected = lengths[binder_chain]
        bid = f"NB{idx:02d}_{item['family']}_s{item['seed']}"
        bdir = OUT / bid
        bdir.mkdir(parents=True, exist_ok=True)

        payload = {
            "input_pdb": pdb.read_text(encoding="utf-8"),
            "input_pdb_chains": [binder_chain],
            "num_seq_per_target": NUM_SEQ_PER_TARGET,
            "sampling_temp": TEMPS,
            "use_soluble_model": False,
            "ca_only": False,
            "omit_AAs": ["C"],
            "random_seed": 20261200 + idx,
        }
        result = post_json(payload, api_key)
        mfasta = result.get("mfasta", "")
        if not mfasta.strip():
            raise RuntimeError(f"No mfasta returned for {bid}")
        (bdir / "designed_sequences.fasta").write_text(mfasta, encoding="utf-8")

        designed = [
            (h, s) for h, s in parse_fasta(mfasta)
            if not h.lower().startswith("input")
            and "native" not in h.lower()
            and not h.lower().startswith("wt")
        ]
        local = []
        for j, (header, sequence) in enumerate(designed, start=1):
            binder = extract_binder(sequence, order, binder_chain)
            if len(binder) != expected:
                raise RuntimeError(f"{bid}: expected {expected} aa, got {len(binder)}")
            acidic = binder.count("D") + binder.count("E")
            basic = binder.count("K") + binder.count("R")
            charge_fraction = (acidic + basic) / len(binder)
            top_fraction = max(binder.count(x) for x in set(binder)) / len(binder)
            flags = []
            if charge_fraction > 0.65:
                flags.append("charge_fraction_gt_0.65")
            if top_fraction > 0.40:
                flags.append("single_residue_fraction_gt_0.40")
            if longest_run(binder) > 4:
                flags.append("homopolymer_run_gt_4")
            if "C" in binder:
                flags.append("free_cysteine")

            row = {
                "candidate_id": f"TFR1_{bid}_{j:02d}",
                "backbone_id": bid,
                "source_family": item["family"],
                "source_seed": item["seed"],
                "source_backbone": item["selected_copy"],
                "binder_chain": binder_chain,
                "target_chain": target_chain,
                "binder_length": len(binder),
                "sequence": binder,
                "temperature": hv(header, r"(?:^|,\s*)T=([^,]+)"),
                "sample": hv(header, r"(?:^|,\s*)sample=([^,]+)"),
                "mpnn_score": hv(header, r"(?:^|,\s*)score=([^,]+)"),
                "global_score": hv(header, r"(?:^|,\s*)global_score=([^,]+)"),
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
            rows_all.append(row)
            local.append(row)

        with (bdir / "candidates.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=local[0].keys())
            w.writeheader()
            w.writerows(local)

        (bdir / "request_metadata.json").write_text(
            json.dumps({
                "source_backbone": item["selected_copy"],
                "source_family": item["family"],
                "source_seed": item["seed"],
                "source_qc": item["qc"],
                "chain_order": order,
                "chain_lengths": lengths,
                "designed_binder_chain": binder_chain,
                "fixed_target_chain": target_chain,
                "policy": "pose-first; no forced histidine",
            }, indent=2) + "\n",
            encoding="utf-8",
        )
        backbone_meta.append({
            "backbone_id": bid,
            "source": item["selected_copy"],
            "candidate_count": len(local),
            "prefilter_pass_count": sum(r["prefilter_pass"] for r in local),
        })
        time.sleep(10)

    chosen = []
    for meta in backbone_meta:
        rr = [x for x in rows_all if x["backbone_id"] == meta["backbone_id"] and x["prefilter_pass"]]
        rr.sort(key=lambda x: float(x["mpnn_score"]) if x["mpnn_score"] else 999.0)
        chosen.extend(rr[:2])

    fields = rows_all[0].keys()
    for filename, rows in [
        ("all_candidates.csv", rows_all),
        ("pose_screen_candidates.csv", chosen),
    ]:
        with (OUT / filename).open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

    (OUT / "result_summary.json").write_text(
        json.dumps({
            "policy": "pose-first",
            "backbones_used": len(backbone_meta),
            "candidate_count": len(rows_all),
            "pose_screen_candidate_count": len(chosen),
            "pose_screen_ids": [x["candidate_id"] for x in chosen],
            "backbones": backbone_meta,
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    print((OUT / "result_summary.json").read_text(), flush=True)


if __name__ == "__main__":
    main()
