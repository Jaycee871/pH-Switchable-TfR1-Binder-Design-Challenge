#!/usr/bin/env python3
"""Enumerate mechanistically motivated histidine-switch variants after pose validation.

This script is deliberately downstream of the strict human/mouse pose gate.
It does NOT score pH-dependent binding and must not be used to claim a pH switch.

Rationale
---------
Histidine is partially protonated near endosomal pH. At an interface, a His
substitution near a conserved basic target residue (Lys/Arg) can become less
favourable at acidic pH through electrostatic repulsion, while remaining mostly
neutral around pH 7.4. Conversely, His placed next to acidic Asp/Glu can become
more favourable when protonated, which is the wrong direction for this challenge.

The script therefore:
- scans only binder residues that actually contact human TfR1 in a full-complex
  pose that already passed the cross-species gate;
- favours positions contacting conserved/basic TfR1 residues such as K344/K374;
- penalises positions contacting acidic E/D residues;
- excludes Gly/Pro/Cys and already-His positions from mutation;
- emits singles and conservative double-His combinations for later explicit
  protonation / pKa / energy evaluation.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from Bio.PDB import MMCIFParser

ROOT = Path(__file__).resolve().parents[1]
POSE_ROOT = ROOT / "predictions/freebindcraft_full_pose_screen"
PASS_IDS = POSE_ROOT / "pose_pass_ids.json"
BATCH = POSE_ROOT / "batch_summary.json"
OUT = ROOT / "designs/pH_switch/freebindcraft_histidine_scan"

HUMAN_START = 121
# Conserved, exposed/basic positions are attractive for acidic-pH repulsion.
BASIC_FAVOUR = {207, 214, 344, 374}
# Acidic target contacts are penalised because protonated His may strengthen them.
ACIDIC_AVOID = {343, 350}
# Species-risk positions should not become the sole switch anchor.
SPECIES_RISK = {204, 208, 210, 211, 212, 348, 371}
EXCLUDED_SOURCE_AA = {"G", "P", "C", "H"}

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}


def distance(a, b):
    return math.sqrt(sum((float(x)-float(y))**2 for x,y in zip(a.coord,b.coord)))


def load_model(path: Path):
    return next(MMCIFParser(QUIET=True).get_structure(path.stem, str(path)).get_models())


def binder_sequence_from_structure(chain):
    seq = []
    positions = []
    for r in chain:
        if r.id[0] != " ":
            continue
        aa = AA3_TO_1.get(r.resname)
        if aa:
            seq.append(aa)
            positions.append(r.id[1])
    return "".join(seq), positions


def scan_candidate(cid: str, sequence: str):
    cif = POSE_ROOT / cid / "human/structure_1.cif"
    if not cif.exists():
        raise FileNotFoundError(cif)

    model = load_model(cif)
    target = model["A"]
    binder = model["B"]

    struct_seq, struct_pos = binder_sequence_from_structure(binder)
    if len(struct_seq) != len(sequence):
        raise RuntimeError(
            f"{cid}: binder length mismatch structure={len(struct_seq)} sequence={len(sequence)}"
        )

    binder_index = {resid:i+1 for i,resid in enumerate(struct_pos)}
    rows = []

    target_residues = [r for r in target if r.id[0] == " "]
    binder_residues = [r for r in binder if r.id[0] == " "]

    for br in binder_residues:
        bpos = binder_index[br.id[1]]
        src = sequence[bpos-1]
        if src in EXCLUDED_SOURCE_AA:
            continue

        contacted = set()
        min_d = {}
        for tr in target_residues:
            best = 999.0
            for ba in br:
                if ba.element == "H":
                    continue
                for ta in tr:
                    if ta.element == "H":
                        continue
                    d = distance(ba, ta)
                    if d < best:
                        best = d
            if best <= 5.0:
                uniprot = tr.id[1] + HUMAN_START - 1
                contacted.add(uniprot)
                min_d[uniprot] = round(best, 3)

        if not contacted:
            continue

        basic = sorted(contacted & BASIC_FAVOUR)
        acidic = sorted(contacted & ACIDIC_AVOID)
        risk = sorted(contacted & SPECIES_RISK)

        # A transparent heuristic only for prioritising mutations to evaluate.
        score = 0.0
        score += 4.0 * len(basic)
        score -= 3.0 * len(acidic)
        score -= 2.0 * len(risk)
        score += min(2.0, 0.15 * len(contacted))

        rows.append({
            "candidate_id": cid,
            "binder_position": bpos,
            "source_aa": src,
            "mutation": f"{src}{bpos}H",
            "contacted_target_residues": sorted(contacted),
            "basic_favour_contacts": basic,
            "acidic_avoid_contacts": acidic,
            "species_risk_contacts": risk,
            "min_distances_A": min_d,
            "heuristic_priority_score": round(score, 3),
        })

    rows.sort(
        key=lambda r: (
            r["heuristic_priority_score"],
            len(r["basic_favour_contacts"]),
            -len(r["species_risk_contacts"]),
        ),
        reverse=True,
    )
    return rows


def mutate(seq: str, positions):
    chars = list(seq)
    for pos in positions:
        chars[pos-1] = "H"
    return "".join(chars)


def main():
    if not PASS_IDS.exists() or not BATCH.exists():
        raise SystemExit(
            "No pose-pass output yet. Run FreeBindCraft generation and full human/mouse validation first."
        )

    passed = json.loads(PASS_IDS.read_text())["pose_pass_ids"]
    batch = json.loads(BATCH.read_text())
    OUT.mkdir(parents=True, exist_ok=True)

    all_ranked = []
    variants = []

    for cid in passed:
        sequence = batch[cid]["binder_sequence"]
        ranked = scan_candidate(cid, sequence)
        all_ranked.extend(ranked)

        # Singles: top 8 interface positions.
        for row in ranked[:8]:
            p = int(row["binder_position"])
            variants.append({
                "parent_id": cid,
                "variant_id": f"{cid}_H{p}",
                "positions": [p],
                "mutations": [row["mutation"]],
                "sequence": mutate(sequence, [p]),
                "priority_score": row["heuristic_priority_score"],
                "mechanism": "single_His_interface_scan",
            })

        # Doubles: only among the top 6, avoid adjacent positions and require at
        # least one position with a favourable basic-target contact.
        top = ranked[:6]
        for i in range(len(top)):
            for j in range(i+1, len(top)):
                a, b = top[i], top[j]
                pa, pb = int(a["binder_position"]), int(b["binder_position"])
                if abs(pa-pb) <= 2:
                    continue
                if not (a["basic_favour_contacts"] or b["basic_favour_contacts"]):
                    continue
                positions = sorted([pa,pb])
                variants.append({
                    "parent_id": cid,
                    "variant_id": f"{cid}_H{positions[0]}_H{positions[1]}",
                    "positions": positions,
                    "mutations": [a["mutation"],b["mutation"]],
                    "sequence": mutate(sequence, positions),
                    "priority_score": round(
                        a["heuristic_priority_score"] + b["heuristic_priority_score"], 3
                    ),
                    "mechanism": "double_His_interface_scan",
                })

    (OUT / "interface_histidine_positions.json").write_text(
        json.dumps(all_ranked, indent=2) + "\n"
    )
    variants.sort(key=lambda x:x["priority_score"], reverse=True)
    (OUT / "histidine_variants.json").write_text(
        json.dumps(variants, indent=2) + "\n"
    )

    if variants:
        with (OUT / "histidine_variants.csv").open("w", newline="") as fh:
            fields = [
                "parent_id","variant_id","positions","mutations","sequence",
                "priority_score","mechanism",
            ]
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            for x in variants:
                y = dict(x)
                y["positions"] = ";".join(map(str,y["positions"]))
                y["mutations"] = ";".join(y["mutations"])
                w.writerow(y)

    summary = {
        "pose_pass_parent_count": len(passed),
        "ranked_interface_positions": len(all_ranked),
        "generated_variants": len(variants),
        "warning": (
            "These variants are hypotheses only. Explicit pKa/protonation and "
            "pH-dependent energetic evaluation are required before ranking for submission."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
