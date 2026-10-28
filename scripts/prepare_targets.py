#!/usr/bin/env python3
"""Prepare and validate Challenge 04 TfR1 target assets.

No third-party Python packages are required.

Outputs:
- structures/11QC.pdb        human TfR1 ectodomain reference
- structures/11RN.pdb        human TfR1 + VP1u apical-domain reference
- data/human_TfR1/...fasta   assay construct sequence
- data/mouse_TfR1/...fasta   assay construct sequence
- metadata/target_manifest.json

The official Proteinbase structure-aware hotspot correspondence should be used
for apical-domain cross-species design. Sequence alignment alone is ambiguous
around the local indel near human residues 210-212.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]

HUMAN_ACCESSION = "P02786"
MOUSE_ACCESSION = "Q62351"
HUMAN_START, HUMAN_END = 121, 760
MOUSE_START, MOUSE_END = 123, 763

URLS = {
    "human_fasta": f"https://rest.uniprot.org/uniprotkb/{HUMAN_ACCESSION}.fasta",
    "mouse_fasta": f"https://rest.uniprot.org/uniprotkb/{MOUSE_ACCESSION}.fasta",
    "11QC": "https://files.rcsb.org/download/11QC.pdb",
    "11RN": "https://files.rcsb.org/download/11RN.pdb",
}

OFFICIAL_HOTSPOT_MAP = [
    {"human": "D204", "mouse": "Q207"},
    {"human": "R208", "mouse": "N211"},
    {"human": "V210", "mouse": "D213"},
    {"human": "Y211", "mouse": "P214"},
    {"human": "L212", "mouse": None},
    {"human": "N348", "mouse": "K350"},
    {"human": "K371", "mouse": "Q373"},
]


def fetch_text(url: str) -> str:
    req = Request(url, headers={"User-Agent": "TfR1-Challenge-TargetPrep/1.0"})
    with urlopen(req, timeout=60) as response:
        return response.read().decode("utf-8")


def parse_fasta(text: str) -> str:
    return "".join(
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.startswith(">")
    )


def slice_uniprot(sequence: str, start: int, end: int) -> str:
    # UniProt numbering is 1-based and inclusive.
    return sequence[start - 1 : end]


def sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def fasta(header: str, seq: str, width: int = 80) -> str:
    rows = [seq[i : i + width] for i in range(0, len(seq), width)]
    return f">{header}\n" + "\n".join(rows) + "\n"


def main() -> None:
    human_full = parse_fasta(fetch_text(URLS["human_fasta"]))
    mouse_full = parse_fasta(fetch_text(URLS["mouse_fasta"]))

    human = slice_uniprot(human_full, HUMAN_START, HUMAN_END)
    mouse = slice_uniprot(mouse_full, MOUSE_START, MOUSE_END)

    assert len(human) == 640, f"Unexpected human construct length: {len(human)}"
    assert len(mouse) == 641, f"Unexpected mouse construct length: {len(mouse)}"

    human_path = ROOT / "data/human_TfR1/human_TfR1_P02786_121-760.fasta"
    mouse_path = ROOT / "data/mouse_TfR1/mouse_TfR1_Q62351_123-763.fasta"

    write(
        human_path,
        fasta(
            "human_TfR1|P02786-1|ectodomain|121-760|Proteinbase_Challenge04",
            human,
        ),
    )
    write(
        mouse_path,
        fasta(
            "mouse_TfR1|Q62351|ectodomain|123-763|Proteinbase_Challenge04",
            mouse,
        ),
    )

    pdb_status = {}
    for pdb_id in ("11QC", "11RN"):
        try:
            pdb_text = fetch_text(URLS[pdb_id])
            write(ROOT / f"structures/{pdb_id}.pdb", pdb_text)
            pdb_status[pdb_id] = {
                "downloaded": True,
                "sha256": sha256(pdb_text),
                "bytes": len(pdb_text.encode("utf-8")),
            }
        except Exception as exc:
            pdb_status[pdb_id] = {
                "downloaded": False,
                "error": str(exc),
            }

    manifest = {
        "challenge": "Challenge 04: TfR1",
        "human": {
            "accession": "P02786-1",
            "residues": [HUMAN_START, HUMAN_END],
            "length": len(human),
            "sha256": sha256(human),
            "structure": "11QC",
        },
        "mouse": {
            "accession": "Q62351",
            "residues": [MOUSE_START, MOUSE_END],
            "length": len(mouse),
            "sha256": sha256(mouse),
        },
        "reference_complex": "11RN",
        "official_hotspot_correspondence": OFFICIAL_HOTSPOT_MAP,
        "design_states": {
            "positive": ["human_pH7.4", "mouse_pH7.4"],
            "negative": ["human_pH6.0"],
        },
        "pdb_downloads": pdb_status,
        "sources": URLS,
    }

    write(
        ROOT / "metadata/target_manifest.json",
        json.dumps(manifest, indent=2) + "\n",
    )

    print("Target preparation complete")
    print(f"  human construct: {len(human)} aa")
    print(f"  mouse construct: {len(mouse)} aa")
    print("  hotspot mapping: Proteinbase structure-aware correspondence")
    for pdb_id, status in pdb_status.items():
        print(f"  {pdb_id}: {'OK' if status['downloaded'] else 'download failed'}")


if __name__ == "__main__":
    main()
