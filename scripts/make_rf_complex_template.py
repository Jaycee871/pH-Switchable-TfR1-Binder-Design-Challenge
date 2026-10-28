#!/usr/bin/env python3
"""Convert an RFdiffusion TfR1 complex into a sequence-compatible Boltz2 template.

RFdiffusion output layout for seed 20261029:
- chain A: 87-residue poly-G binder
- chain B: human TfR1 residues 122-760, renumbered 88-726

Boltz2 query layout:
- chain A: target (human 121-760 or mouse 123-763)
- chain B: designed 87-residue binder

The binder's backbone coordinates are retained while residue names are changed
from GLY to the ProteinMPNN A29_03 sequence. Only atoms present in the RF
backbone are retained; Boltz2 rebuilds the final all-atom structure.
"""

from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

BINDER_SEQUENCE = (
    "ALDEEVEKLLKEVQALEKEVAALKKKKKELEEKAKKEKNEEAEKIKEEAKKLGKEKEEKEEKL"
    "KELEKKLSEKAKKYGDKLEEEKKK"
)

AA3 = {
    "A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","E":"GLU","Q":"GLN",
    "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
    "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL",
}

DEFAULT_INPUT = ROOT / "recovered/seeded_selection/rfdiffusion_seed_20261029.pdb"
DEFAULT_OUTPUT = ROOT / "predictions/templates/A29_03_rf_seed20261029_template.pdb"


def rewrite_atom_line(line: str, new_chain: str, new_resid: int, new_resname: str | None = None) -> str:
    chars = list(line.rstrip("\n"))
    if len(chars) < 80:
        chars.extend(" " * (80 - len(chars)))
    if new_resname is not None:
        res = f"{new_resname:>3}"
        chars[17:20] = list(res)
    chars[21] = new_chain
    resid = f"{new_resid:4d}"
    chars[22:26] = list(resid)
    return "".join(chars).rstrip() + "\n"


def main() -> None:
    src = Path(os.getenv("RF_PDB", str(DEFAULT_INPUT)))
    dst = Path(os.getenv("TEMPLATE_PDB", str(DEFAULT_OUTPUT)))
    if not src.exists():
        raise RuntimeError(f"RFdiffusion PDB not found: {src}")
    if len(BINDER_SEQUENCE) != 87:
        raise RuntimeError("Expected 87-aa A29_03 binder")

    target_lines = []
    binder_lines = []

    for line in src.read_text(encoding="utf-8").splitlines(True):
        if not line.startswith("ATOM"):
            continue
        chain = line[21]
        try:
            old_resid = int(line[22:26])
        except ValueError:
            continue

        if chain == "A":
            # RF binder -> Boltz binder chain B, residues 1..87.
            if not (1 <= old_resid <= len(BINDER_SEQUENCE)):
                continue
            aa = BINDER_SEQUENCE[old_resid - 1]
            binder_lines.append(
                rewrite_atom_line(line, "B", old_resid, AA3[aa])
            )
        elif chain == "B":
            # RF target starts at human residue 122. In the full human 121-760
            # query this maps to query positions 2..640.
            new_resid = old_resid - 86
            target_lines.append(
                rewrite_atom_line(line, "A", new_resid, None)
            )

    if not target_lines or not binder_lines:
        raise RuntimeError("Failed to recover target/binder chains from RF PDB")

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(
        "".join(target_lines) + "TER\n" + "".join(binder_lines) + "TER\nEND\n",
        encoding="utf-8",
    )

    print(f"Wrote template: {dst}")
    print(f"  target atoms: {len(target_lines)}")
    print(f"  binder atoms: {len(binder_lines)}")
    print("  target template query positions: 2-640")
    print("  binder template query positions: 1-87")


if __name__ == "__main__":
    main()
