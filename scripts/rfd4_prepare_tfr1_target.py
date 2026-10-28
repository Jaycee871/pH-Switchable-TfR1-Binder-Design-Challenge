#!/usr/bin/env python3
"""Prepare a chemistry-valid mmCIF crop of the human TfR1 apical domain.

Source: PDB 11QC chain A.
Crop: human UniProt/PDB residue numbers 190-400.

The crop retains original residue numbering so RFD4 condition selections can
refer directly to E343, K344, G347 and E350.
"""

from __future__ import annotations

import json
from pathlib import Path

from Bio.PDB import MMCIFIO, PDBParser, Select

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "structures/11QC.pdb"
OUTDIR = ROOT / "structures/rfd4"
OUT = OUTDIR / "human_TfR1_11QC_A_190-400.cif"
META = OUTDIR / "human_TfR1_11QC_A_190-400.json"


class CropSelect(Select):
    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        return (
            residue.id[0] == " "
            and 190 <= residue.id[1] <= 400
        )


def main() -> None:
    if not SOURCE.exists():
        raise RuntimeError(
            "structures/11QC.pdb is missing; run scripts/prepare_targets.py first"
        )

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("11QC", str(SOURCE))

    selected = []
    for model in structure:
        for chain in model:
            if chain.id != "A":
                continue
            for residue in chain:
                if residue.id[0] == " " and 190 <= residue.id[1] <= 400:
                    selected.append(
                        {
                            "resid": residue.id[1],
                            "resname": residue.resname,
                            "atoms": sorted(a.name for a in residue),
                        }
                    )
        break

    if not selected:
        raise RuntimeError("No residues selected from 11QC chain A")

    present = {r["resid"] for r in selected}
    for required in (343, 344, 347, 350):
        if required not in present:
            raise RuntimeError(f"Required hotspot residue {required} is absent")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(OUT), CropSelect())

    hotspot_details = {
        r["resid"]: {
            "resname": r["resname"],
            "atoms": r["atoms"],
        }
        for r in selected
        if r["resid"] in {343, 344, 347, 350}
    }

    metadata = {
        "source": "PDB 11QC chain A",
        "residue_range": [190, 400],
        "residue_count": len(selected),
        "first_residue": min(present),
        "last_residue": max(present),
        "hotspots": hotspot_details,
        "output": str(OUT.relative_to(ROOT)),
    }
    META.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
