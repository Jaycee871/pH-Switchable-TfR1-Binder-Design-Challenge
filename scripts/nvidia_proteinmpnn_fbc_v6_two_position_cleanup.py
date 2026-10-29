#!/usr/bin/env python3
"""Phase 3J: two-position cleanup of the FreeBindCraft v6 tempered-binary binder.

v6 is the first FreeBindCraft parent that is both relaxed and genuinely engages
the conserved hotspot core (347/350) with zero severe clashes. Its remaining
species-risk contacts are geometrically concentrated at binder residues B46/B47:

- B46=W contacts human 210/211 + 347/348 + 371
- B47=E contacts human 371

The rest of the useful hotspot interface is supported independently by B23,
B42-B45, B52 and B54, especially B44/B45 around 347/350.

This script therefore freezes EVERY binder residue except 46 and 47, and asks
ProteinMPNN only for a narrow chemistry cleanup at those two positions.

It is cross-species pose cleanup, not pH engineering.
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
BACKBONE = ROOT / (
    "freebindcraft/kaggle_v6_tempered_binary_results/designs/Trajectory/Relaxed/"
    "TFR1_FBC343_V6_l60_s289817.pdb"
)
OUT = ROOT / "designs/family_B_multistate/fbc_v6_two_position_cleanup"
URL = "https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict"

ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
LENGTH = 60
DESIGN_POSITIONS = [46, 47]
FIXED_POSITIONS = [p for p in range(1, LENGTH + 1) if p not in DESIGN_POSITIONS]

# Parent residues: B46=W, B47=E.
# W46 is broad/aromatic and bridges several species-risk surfaces.
# E47 touches only K371 (mouse Q373). Prefer compact/neutral chemistry.
P46_PREFERRED = set("AGNQSTV")
P47_PREFERRED = set("AGNQST")


def post(payload, key, max_attempts=4):
    for attempt in range(1, max_attempts + 1):
        req = Request(
            URL,
            data=json.dumps(payload).encode(),
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {key}",
            },
            method="POST",
        )
        try:
            with urlopen(req, timeout=1200) as r:
                return json.loads(r.read().decode())
        except HTTPError as exc:
            detail = exc.read().decode(errors="replace")
            if exc.code == 429 and attempt < max_attempts:
                ra = exc.headers.get("Retry-After")
                sleep = int(ra) if ra and ra.isdigit() else 20 * attempt
                print(f"429 attempt {attempt}; sleep {sleep}s", flush=True)
                time.sleep(sleep)
                continue
            raise RuntimeError(f"ProteinMPNN HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("ProteinMPNN retries exhausted")


def strip_hydrogens(text):
    out = []
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom = line[12:16].strip()
        elem = line[76:78].strip() if len(line) >= 78 else ""
        if atom.startswith("H") or elem == "H":
            continue
        out.append(line)
    return "\n".join(out) + "\nTER\nEND\n"


def blank():
    return [0.0] * len(ALPHABET)


def add(row, aas, value):
    for aa in aas:
        row[ALPHABET.index(aa)] += value


def bias_matrix():
    m = [blank() for _ in range(LENGTH)]

    # B46: shrink/remove broad aromatic footprint while retaining interface-capable
    # neutral chemistry. Strongly penalise aromatics and charged residues.
    add(m[45], "AGNQSTV", 2.0)
    add(m[45], "FWY", -2.0)
    add(m[45], "DEKR", -1.0)
    add(m[45], "HILM", -0.5)

    # B47: remove E371/K371 dependence. Neutral/small residues preferred.
    add(m[46], "AGNQST", 2.0)
    add(m[46], "DEKR", -1.5)
    add(m[46], "FWYH", -0.8)

    return m


def parse_fasta(text):
    out, h, s = [], None, []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if h is not None:
                out.append((h, "".join(s)))
            h, s = line[1:], []
        else:
            s.append(line)
    if h is not None:
        out.append((h, "".join(s)))
    return out


def hv(header, pat):
    m = re.search(pat, header)
    return m.group(1).strip() if m else ""


def binder_segment(sequence):
    parts = [x for x in sequence.split("/") if x]
    matches = [x for x in parts if len(x) == LENGTH]
    if len(matches) == 1:
        return matches[0]
    if len(parts) == 1 and len(parts[0]) == LENGTH:
        return parts[0]
    raise RuntimeError(
        f"Could not identify {LENGTH}-aa binder segment: {[len(x) for x in parts]}"
    )


def parent_sequence_from_pdb():
    aa3 = {
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q",
        "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
        "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    }
    seen = {}
    for line in BACKBONE.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM") or line[21].strip() != "B":
            continue
        r = int(line[22:26])
        if r not in seen:
            seen[r] = aa3[line[17:20].strip()]
    seq = "".join(seen[k] for k in sorted(seen))
    if len(seq) != LENGTH:
        raise RuntimeError(f"Expected {LENGTH}-aa parent, got {len(seq)}")
    return seq


def main():
    key = (
        os.getenv("NVIDIA_API_KEY")
        or os.getenv("NVIDIA_API_KEY1")
        or os.getenv("NVIDIA_API_TOKEN")
        or os.getenv("NGC_API_KEY")
    )
    if not key:
        raise RuntimeError("Missing NVIDIA API credential")
    if not BACKBONE.exists():
        raise RuntimeError(f"Missing v6 backbone: {BACKBONE}")

    parent = parent_sequence_from_pdb()
    if parent[45] != "W" or parent[46] != "E":
        raise RuntimeError(
            f"Unexpected parent B46/B47: {parent[45]}{parent[46]} (expected WE)"
        )

    bias = bias_matrix()
    fixed = json.dumps({"input": {"B": FIXED_POSITIONS}})
    payload = {
        "input_pdb": strip_hydrogens(BACKBONE.read_text(errors="replace")),
        "input_pdb_chains": ["B"],
        "num_seq_per_target": 12,
        "sampling_temp": [0.10, 0.20, 0.30],
        "use_soluble_model": False,
        "ca_only": False,
        "omit_AAs": ["C"],
        "fixed_positions_jsonl": fixed,
        "bias_by_res_jsonl": json.dumps({"input": {"B": bias}}),
        "random_seed": 20261101,
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "design_spec.json").write_text(
        json.dumps(
            {
                "parent_id": "TFR1_FBC343_V6_l60_s289817",
                "parent_sequence": parent,
                "design_positions": DESIGN_POSITIONS,
                "fixed_positions": FIXED_POSITIONS,
                "parent_46_47": parent[45:47],
                "preferred_46": sorted(P46_PREFERRED),
                "preferred_47": sorted(P47_PREFERRED),
                "structural_rationale": {
                    "B46_W": "contacts risk 210/211/348/371 while also contacting 347",
                    "B47_E": "contacts risk 371 only",
                    "preserved_hotspot_support": "B23/B42-B45/B52/B54 support 347/350 independently",
                },
            },
            indent=2,
        )
        + "\n"
    )

    result = post(payload, key)
    mfasta = result.get("mfasta", "")
    if not mfasta.strip():
        raise RuntimeError("No mfasta returned")
    (OUT / "designed_sequences.fasta").write_text(mfasta)

    records = [
        (h, s)
        for h, s in parse_fasta(mfasta)
        if not h.lower().startswith("input")
        and "native" not in h.lower()
        and not h.lower().startswith("wt")
    ]

    rows = []
    seen_seq = set()
    for i, (header, sequence) in enumerate(records, 1):
        binder = binder_segment(sequence)
        if binder in seen_seq:
            continue
        seen_seq.add(binder)

        changed = [p for p in range(1, LENGTH + 1) if binder[p - 1] != parent[p - 1]]
        only_46_47 = all(p in DESIGN_POSITIONS for p in changed)
        p46, p47 = binder[45], binder[46]

        flags = []
        if not only_46_47:
            flags.append("mutation_outside_46_47")
        if p46 not in P46_PREFERRED:
            flags.append("B46_not_small_neutral")
        if p47 not in P47_PREFERRED:
            flags.append("B47_not_small_neutral")
        if p46 == "H" or p47 == "H":
            flags.append("premature_histidine")
        if binder[43] != parent[43] or binder[44] != parent[44]:
            flags.append("B44_B45_anchor_changed")

        rows.append(
            {
                "candidate_id": f"TFR1_FBC_V6_XS_{i:02d}",
                "sequence": binder,
                "mpnn_score": hv(header, r"(?:^|,\s*)score=([^,]+)"),
                "global_score": hv(header, r"(?:^|,\s*)global_score=([^,]+)"),
                "temperature": hv(header, r"(?:^|,\s*)T=([^,]+)"),
                "B46": p46,
                "B47": p47,
                "changed_positions": ";".join(map(str, changed)),
                "only_46_47_changed": only_46_47,
                "B44_anchor": binder[43],
                "B45_anchor": binder[44],
                "histidine_count": binder.count("H"),
                "prefilter_flags": ";".join(flags),
                "prefilter_pass": not flags,
                "header": header,
            }
        )

    if not rows:
        raise RuntimeError("No designed sequences parsed")

    fields = list(rows[0].keys())
    with (OUT / "all_candidates.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    passing = [x for x in rows if x["prefilter_pass"]]
    passing.sort(
        key=lambda x: (
            float(x["mpnn_score"]) if x["mpnn_score"] else 999,
            x["B46"],
            x["B47"],
        )
    )
    chosen = passing[:4]

    with (OUT / "pose_screen_candidates.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(chosen)

    (OUT / "result_summary.json").write_text(
        json.dumps(
            {
                "candidate_count": len(rows),
                "unique_candidate_count": len(seen_seq),
                "prefilter_pass_count": len(passing),
                "pose_screen_ids": [x["candidate_id"] for x in chosen],
                "pose_screen_46_47": {
                    x["candidate_id"]: x["B46"] + x["B47"] for x in chosen
                },
                "note": (
                    "Only B46/B47 are allowed to change. Full human/mouse pose "
                    "validation is mandatory; no pH claim is made."
                ),
            },
            indent=2,
        )
        + "\n"
    )
    print((OUT / "result_summary.json").read_text())


if __name__ == "__main__":
    main()
