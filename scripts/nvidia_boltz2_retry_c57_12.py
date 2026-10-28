#!/usr/bin/env python3
"""Rate-limit-safe Boltz2 validation for compact candidate C57_12 only."""

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
RF_PDB = Path(os.getenv("RF_PDB", str(ROOT / "recovered/selected_backbone.pdb")))
OUTDIR = ROOT / "predictions/C57_12_retry"

CANDIDATE_ID = "TFR1_C57_12"
BINDER = "GVEEAEKELEKRKEEELKKAKEEGEEKIKEIEEKYKKEKESLEELRERMEEQKKRLE"

TARGETS = {
    "human": {
        "fasta": ROOT / "data/human_TfR1/human_TfR1_P02786_121-760.fasta",
        "offset": 120,
        "preferred_core": {344, 347},
        "risk_core": {348},
    },
    "mouse": {
        "fasta": ROOT / "data/mouse_TfR1/mouse_TfR1_Q62351_123-763.fasta",
        "offset": 122,
        "preferred_core": {346, 349},
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
        line.strip() for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    )


def rewrite(line: str, chain: str, resid: int, resname=None) -> str:
    chars=list(line.rstrip("\n"))
    if len(chars)<80:
        chars.extend(" "*(80-len(chars)))
    if resname is not None:
        chars[17:20]=list(f"{resname:>3}")
    chars[21]=chain
    chars[22:26]=list(f"{resid:4d}")
    return "".join(chars).rstrip()+"\n"


def build_template() -> Path:
    target_lines=[]
    binder_lines=[]
    target_first_output=58  # 57-aa binder occupies output residues 1..57

    for line in RF_PDB.read_text(encoding="utf-8").splitlines(True):
        if not line.startswith("ATOM"):
            continue
        chain=line[21]
        old=int(line[22:26])
        if chain=="A":
            binder_lines.append(rewrite(line,"B",old,AA3[BINDER[old-1]]))
        elif chain=="B":
            query_pos=old-target_first_output+2
            target_lines.append(rewrite(line,"A",query_pos))

    pdb=OUTDIR/f"{CANDIDATE_ID}_template.pdb"
    cif=OUTDIR/f"{CANDIDATE_ID}_template.cif"
    OUTDIR.mkdir(parents=True, exist_ok=True)
    pdb.write_text("".join(target_lines)+"TER\n"+"".join(binder_lines)+"TER\nEND\n",encoding="utf-8")
    structure=PDBParser(QUIET=True).get_structure("tmpl",str(pdb))
    io=MMCIFIO()
    io.set_structure(structure)
    io.save(str(cif))
    return cif


def post_json(payload: dict, api_key: str, max_attempts: int = 4) -> dict:
    for attempt in range(1,max_attempts+1):
        req=Request(
            URL,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type":"application/json","Authorization":f"Bearer {api_key}"},
            method="POST",
        )
        try:
            with urlopen(req,timeout=1200) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except HTTPError as exc:
            detail=exc.read().decode("utf-8",errors="replace")
            if exc.code==429 and attempt<max_attempts:
                retry_after=exc.headers.get("Retry-After")
                sleep_s=int(retry_after) if retry_after and retry_after.isdigit() else 20*attempt
                print(f"Boltz2 429 on attempt {attempt}; retrying after {sleep_s}s")
                time.sleep(sleep_s)
                continue
            raise RuntimeError(f"Boltz2 HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("Boltz2 retry loop exhausted")


def contacts(cif_path: Path, offset: int, preferred: set[int], risk: set[int], cutoff=5.0):
    structure=MMCIFParser(QUIET=True).get_structure("pred",str(cif_path))
    model=next(structure.get_models())
    target=model["A"]
    binder=model["B"]
    target_atoms=[
        atom for residue in target if residue.id[0]==" "
        for atom in residue if atom.element!="H"
    ]
    search=NeighborSearch(target_atoms)
    tset=set()
    bset=set()
    for residue in binder:
        if residue.id[0]!=" ":
            continue
        for atom in residue:
            if atom.element=="H":
                continue
            near=search.search(atom.coord,cutoff,level="A")
            if near:
                bset.add(residue.id[1])
            for a in near:
                tset.add(a.get_parent().id[1]+offset)
    t=sorted(tset)
    return {
        "target_contacts_uniprot":t,
        "binder_contact_positions":sorted(bset),
        "preferred_core_contacts":sorted(set(t)&preferred),
        "risk_core_contacts":sorted(set(t)&risk),
        "preferred_core_retained":bool(set(t)&preferred),
    }


def main():
    api_key=os.getenv("NVIDIA_API_KEY") or os.getenv("NGC_API_KEY")
    if not api_key:
        raise RuntimeError("Missing NVIDIA_API_KEY / NGC_API_KEY")
    if not RF_PDB.exists():
        raise RuntimeError(f"Missing RF backbone: {RF_PDB}")

    template_path=build_template()
    template=template_path.read_text(encoding="utf-8")
    batch={"candidate_id":CANDIDATE_ID,"sequence":BINDER,"species":{}}

    for idx,(species,spec) in enumerate(TARGETS.items()):
        if idx:
            time.sleep(15)
        target=parse_fasta(spec["fasta"])
        payload={
            "polymers":[
                {
                    "id":"A","molecule_type":"protein","sequence":target,
                    "structural_templates":[{
                        "structure":template,"format":"cif",
                        "name":"C57_12_rf_seed20261102","chain_id":"A"
                    }],
                },
                {
                    "id":"B","molecule_type":"protein","sequence":BINDER,
                    "structural_templates":[{
                        "structure":template,"format":"cif",
                        "name":"C57_12_rf_seed20261102","chain_id":"B"
                    }],
                },
            ],
            "recycling_steps":3,
            "sampling_steps":50,
            "diffusion_samples":1,
            "step_scale":1.638,
            "output_format":"mmcif",
            "write_full_pae":False,
        }

        result=post_json(payload,api_key)
        structure_text=result["structures"][0]["structure"]
        spdir=OUTDIR/species
        spdir.mkdir(parents=True,exist_ok=True)
        cif=spdir/"structure_1.cif"
        cif.write_text(structure_text,encoding="utf-8")

        summary={
            "confidence":result.get("confidence_scores",[None])[0],
            "ptm":result.get("ptm_scores",[None])[0],
            "protein_iptm":result.get("protein_iptm_scores",[None])[0],
            "complex_plddt":result.get("complex_plddt_scores",[None])[0],
            "complex_iplddt":result.get("complex_iplddt_scores",[None])[0],
            "complex_ipde":result.get("complex_ipde_scores",[None])[0],
            "contacts":contacts(
                cif,spec["offset"],spec["preferred_core"],spec["risk_core"]
            ),
        }
        (spdir/"result_summary.json").write_text(
            json.dumps(summary,indent=2)+"\n",encoding="utf-8"
        )
        batch["species"][species]=summary
        print(species,json.dumps(summary,indent=2))

    (OUTDIR/"summary.json").write_text(json.dumps(batch,indent=2)+"\n",encoding="utf-8")
    print("C57_12 retry validation complete")


if __name__=="__main__":
    main()
