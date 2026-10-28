#!/usr/bin/env python3
"""Phase 7A: templated Boltz-2 validation of engineered K344-facing V2H mutants.\n\nTwo low-ProteinMPNN-score parents are mutated at binder position 2 (V2H).\nThe mutation is an explicit engineering hypothesis, not a ProteinMPNN-selected residue.\nEach mutant is evaluated against full human and mouse TfR1 with the gated RFdiffusion geometry as template.\n"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from Bio.PDB import MMCIFIO, MMCIFParser, NeighborSearch, PDBParser

ROOT = Path(__file__).resolve().parents[1]
URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
RF_PDB = Path(os.getenv("RF_PDB", str(ROOT / "structures/rfdiffusion/compact_apical/selected_backbone.pdb")))
OUTDIR = ROOT / "predictions/k344_v2h_validation"

CANDIDATES = {
    "TFR1_K344_V2H_03": "AHEENEKELEKKKKKEIKEAIKKGKEKIEELEKEYKKEKESLEELKKKMEEEEKLMK",
    "TFR1_K344_V2H_05": "AHEENRKELEKKKKKEVEKAKKEGEEKIEEIEKEYKKEIESLEKLENRMKEEEKAQK",
}

TARGETS = {
    "human": {
        "fasta": ROOT / "data/human_TfR1/human_TfR1_P02786_121-760.fasta",
        "offset": 120,
        "preferred_core": {343, 344, 347},
        "risk_core": {348},
    },
    "mouse": {
        "fasta": ROOT / "data/mouse_TfR1/mouse_TfR1_Q62351_123-763.fasta",
        "offset": 122,
        "preferred_core": {345, 346, 349},
        "risk_core": {350},
    },
}

AA3 = {
    "A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","E":"GLU","Q":"GLN",
    "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
    "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL",
}


def parse_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    )


def rewrite_atom_line(line: str, new_chain: str, new_resid: int, new_resname=None) -> str:
    chars = list(line.rstrip("\n"))
    if len(chars) < 80:
        chars.extend(" " * (80 - len(chars)))
    if new_resname is not None:
        chars[17:20] = list(f"{new_resname:>3}")
    chars[21] = new_chain
    chars[22:26] = list(f"{new_resid:4d}")
    return "".join(chars).rstrip() + "\n"


def make_template(sequence: str, prefix: Path) -> Path:
    if len(sequence) != 57:
        raise RuntimeError(f"Expected 57-aa binder, got {len(sequence)}")

    target_lines = []
    binder_lines = []

    for line in RF_PDB.read_text(encoding="utf-8").splitlines(True):
        if not line.startswith("ATOM"):
            continue
        chain = line[21]
        old_resid = int(line[22:26])

        if chain == "A":
            aa = sequence[old_resid - 1]
            binder_lines.append(rewrite_atom_line(line, "B", old_resid, AA3[aa]))
        elif chain == "B":
            # RF target is human 122-760. Full human query is 121-760,
            # therefore RF target residue 122 maps to query position 2.
            binder_len = 57
            target_first_output = binder_len + 1
            query_pos = old_resid - target_first_output + 2
            target_lines.append(rewrite_atom_line(line, "A", query_pos))

    if not target_lines or not binder_lines:
        raise RuntimeError("Failed to build complex template")

    pdb_path = prefix.with_suffix(".pdb")
    cif_path = prefix.with_suffix(".cif")
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        "".join(target_lines) + "TER\n" + "".join(binder_lines) + "TER\nEND\n",
        encoding="utf-8",
    )

    structure = PDBParser(QUIET=True).get_structure(prefix.stem, str(pdb_path))
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(cif_path))
    return cif_path


def post_json(payload: dict, api_key: str, max_attempts: int = 5) -> dict:
    """POST with rate-limit-aware backoff inherited from Phase 4E."""
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
                sleep_s = int(retry_after) if retry_after and retry_after.isdigit() else 25 * attempt
                print(f"Boltz2 429 attempt {attempt}; retrying after {sleep_s}s", flush=True)
                time.sleep(sleep_s)
                continue
            if exc.code in {500, 502, 503, 504} and attempt < max_attempts:
                sleep_s = 20 * attempt
                print(f"Boltz2 HTTP {exc.code} attempt {attempt}; retrying after {sleep_s}s", flush=True)
                time.sleep(sleep_s)
                continue
            raise RuntimeError(f"Boltz2 HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("Boltz2 retry loop exhausted")


def contact_map(
    cif_path: Path,
    offset: int,
    preferred_core: set[int],
    risk_core: set[int],
    cutoff: float = 5.0,
):
    structure = MMCIFParser(QUIET=True).get_structure("pred", str(cif_path))
    model = next(structure.get_models())
    target = model["A"]
    binder = model["B"]

    target_atoms = [
        atom
        for residue in target
        if residue.id[0] == " "
        for atom in residue
        if atom.element != "H"
    ]
    search = NeighborSearch(target_atoms)

    target_contacts = set()
    binder_contacts = set()
    for residue in binder:
        if residue.id[0] != " ":
            continue
        for atom in residue:
            if atom.element == "H":
                continue
            close = search.search(atom.coord, cutoff, level="A")
            if close:
                binder_contacts.add(residue.id[1])
            for target_atom in close:
                target_contacts.add(target_atom.get_parent().id[1] + offset)

    target_contacts = sorted(target_contacts)
    preferred_hits = sorted(set(target_contacts) & preferred_core)
    risk_hits = sorted(set(target_contacts) & risk_core)

    return {
        "cutoff_A": cutoff,
        "target_contacts_uniprot": target_contacts,
        "binder_contact_positions": sorted(binder_contacts),
        "preferred_core": sorted(preferred_core),
        "preferred_core_contacts": preferred_hits,
        "risk_core": sorted(risk_core),
        "risk_core_contacts": risk_hits,
        "preferred_core_retained": bool(preferred_hits),
    }


def main() -> None:
    api_key = (os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_API_KEY1") or os.getenv("NVIDIA_API_TOKEN") or os.getenv("NGC_API_KEY"))
    if not api_key:
        raise RuntimeError("Missing NVIDIA API credential")
    if not RF_PDB.exists():
        raise RuntimeError(f"Missing RF backbone: {RF_PDB}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    batch = {}

    for candidate_id, binder in CANDIDATES.items():
        cdir = OUTDIR / candidate_id
        template_path = make_template(
            binder, cdir / f"{candidate_id}_rf_seed20261102_template"
        )
        template_text = template_path.read_text(encoding="utf-8")

        cresult = {
            "binder_sequence": binder,
            "binder_length": len(binder),
            "core_signature_positions_1_2_5": binder[0] + binder[1] + binder[4],
            "species": {},
        }

        for species, spec in TARGETS.items():
            target = parse_fasta(spec["fasta"])
            payload = {
                "polymers": [
                    {
                        "id": "A",
                        "molecule_type": "protein",
                        "sequence": target,
                        "structural_templates": [{
                            "structure": template_text,
                            "format": "cif",
                            "name": f"{candidate_id}_target",
                            "chain_id": "A",
                        }],
                    },
                    {
                        "id": "B",
                        "molecule_type": "protein",
                        "sequence": binder,
                        "structural_templates": [{
                            "structure": template_text,
                            "format": "cif",
                            "name": f"{candidate_id}_binder",
                            "chain_id": "B",
                        }],
                    },
                ],
                "recycling_steps": 3,
                "sampling_steps": 50,
                "diffusion_samples": 1,
                "step_scale": 1.638,
                "output_format": "mmcif",
                "write_full_pae": False,
            }

            result = post_json(payload, api_key)
            structures = result.get("structures", [])
            if not structures:
                raise RuntimeError(f"No Boltz2 structure for {candidate_id}/{species}")

            spdir = cdir / species
            spdir.mkdir(parents=True, exist_ok=True)
            cif_path = spdir / "structure_1.cif"
            cif_path.write_text(structures[0]["structure"], encoding="utf-8")

            contacts = contact_map(
                cif_path,
                spec["offset"],
                spec["preferred_core"],
                spec["risk_core"],
            )
            protein_iptm = result.get("protein_iptm_scores", [None])[0]
            apical_contact = any(
                (202 <= r <= 215) or (343 <= r <= 374)
                for r in contacts["target_contacts_uniprot"]
            )
            pose_pass = bool(
                protein_iptm is not None
                and float(protein_iptm) >= 0.50
                and contacts["preferred_core_retained"]
                and apical_contact
            )
            summary = {
                "confidence": result.get("confidence_scores", [None])[0],
                "ptm": result.get("ptm_scores", [None])[0],
                "protein_iptm": protein_iptm,
                "complex_plddt": result.get("complex_plddt_scores", [None])[0],
                "complex_iplddt": result.get("complex_iplddt_scores", [None])[0],
                "complex_pde": result.get("complex_pde_scores", [None])[0],
                "complex_ipde": result.get("complex_ipde_scores", [None])[0],
                "pair_chains_iptm": result.get("pair_chains_iptm_scores", []),
                "contacts": contacts,
                "apical_contact": apical_contact,
                "pose_gate": {
                    "min_protein_iptm": 0.50,
                    "requires_preferred_core_retained": True,
                    "pass": pose_pass,
                },
            }
            (spdir / "result_summary.json").write_text(
                json.dumps(summary, indent=2) + "\n", encoding="utf-8"
            )
            cresult["species"][species] = summary
            print(candidate_id, species, json.dumps(summary, indent=2), flush=True)
            time.sleep(20)

        h = cresult["species"].get("human", {})
        m = cresult["species"].get("mouse", {})
        cresult["cross_species_pose_pass"] = bool(
            h.get("pose_gate", {}).get("pass")
            and m.get("pose_gate", {}).get("pass")
        )
        batch[candidate_id] = cresult

    (OUTDIR / "batch_summary.json").write_text(
        json.dumps(batch, indent=2) + "\n", encoding="utf-8"
    )
    print("Phase 7A K344-facing V2H Boltz2 validation complete")


if __name__ == "__main__":
    main()
