#!/usr/bin/env python3
"""Phase 4P: independent sequence-only full human/mouse validation of Odin-Multi selection.

Odin-Multi is allowed to optimize both backbone/pose and sequence jointly across
human and mouse apical contexts. This validator deliberately discards Odin's
optimized structural pose and asks Boltz-2 to predict each full ectodomain +
selected binder complex from sequence only.

Strict gate:
- human protein ipTM >= 0.50
- mouse protein ipTM >= 0.50
- each species retains a preferred apical contact
- each species contacts the apical domain
- >=3 homologous human/mouse contact positions and Jaccard >=0.15

No pH claim is made here.
"""

from __future__ import annotations

import csv
import json
import os
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, NeighborSearch

ROOT = Path(__file__).resolve().parents[1]
SELECTION = Path(os.getenv("ODIN_SELECTION_PATH", str(ROOT / "odin_multi/kaggle_crossreactivity_smoke/selection.csv")))
HUMAN_FASTA = ROOT / "data/human_TfR1/human_TfR1_P02786_121-760.fasta"
MOUSE_FASTA = ROOT / "data/mouse_TfR1/mouse_TfR1_Q62351_123-763.fasta"
MASK_CSV = ROOT / "analysis/apical_design_mask.csv"
OUT = Path(os.getenv("ODIN_VALIDATION_OUT", str(ROOT / "predictions/odin_multi_crossreactivity_pose_screen")))
URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"

HUMAN_START = 121
MOUSE_START = 123

PREFERRED_H = {206,207,209,214,231,343,344,347,350,374}
RISK_H = {204,208,210,211,212,348,371}

OFFICIAL_OVERRIDES = {
    204:207,
    208:211,
    210:213,
    211:214,
    212:None,
    348:350,
    371:373,
}


def fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    )


def human_mouse_map():
    hs, ms = fasta(HUMAN_FASTA), fasta(MOUSE_FASTA)
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    aln = a.align(hs, ms)[0]
    c = aln.coordinates
    mapping = {}
    for k in range(c.shape[1]-1):
        h0,h1 = int(c[0,k]),int(c[0,k+1])
        m0,m1 = int(c[1,k]),int(c[1,k+1])
        dh,dm = h1-h0,m1-m0
        if dh and dm:
            for j in range(min(dh,dm)):
                mapping[HUMAN_START+h0+j] = MOUSE_START+m0+j
        elif dh:
            for j in range(dh):
                mapping[HUMAN_START+h0+j] = None
    mapping.update(OFFICIAL_OVERRIDES)

    if MASK_CSV.exists():
        with MASK_CSV.open(newline="",encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                try:
                    hp=int(row["human_residue"])
                except (ValueError,TypeError):
                    continue
                mp=(row.get("mouse_residue") or "").strip()
                mapping[hp]=int(mp) if mp else None
    return mapping


def read_selection():
    if not SELECTION.exists():
        raise SystemExit(f"Odin selection not present: {SELECTION}")
    with SELECTION.open(newline="",encoding="utf-8") as fh:
        rows=list(csv.DictReader(fh))
    selected=[r for r in rows if (r.get("selection_status") or "").strip()=="selected"]
    if not selected:
        raise SystemExit("Odin selection.csv contains no selected design")
    selected.sort(
        key=lambda r: float(r.get("worst_target_i_ptm") or "-inf"),
        reverse=True,
    )
    rank=int(os.getenv("ODIN_SELECTION_RANK","0"))
    if rank < 0 or rank >= len(selected):
        raise SystemExit(f"ODIN_SELECTION_RANK={rank} outside selected row count {len(selected)}")
    row=selected[rank]
    seq=(row.get("sequence") or "").strip()
    if not seq or len(seq)!=60:
        raise RuntimeError(f"Unexpected Odin binder length: {len(seq)}")
    return row,seq


def post_json(payload,key,max_attempts=5):
    for attempt in range(1,max_attempts+1):
        req=Request(
            URL,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type":"application/json","Authorization":f"Bearer {key}"},
            method="POST",
        )
        try:
            with urlopen(req,timeout=1200) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except HTTPError as exc:
            detail=exc.read().decode("utf-8",errors="replace")
            if exc.code==429 and attempt<max_attempts:
                ra=exc.headers.get("Retry-After")
                delay=int(ra) if ra and ra.isdigit() else 25*attempt
                time.sleep(delay)
                continue
            transient=exc.code in {500,502,503,504} or (
                exc.code==422 and any(
                    x in detail for x in ("InfrastructureError","AcceleratorError","CUDA error")
                )
            )
            if transient and attempt<max_attempts:
                time.sleep(25*attempt)
                continue
            raise RuntimeError(f"Boltz2 HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("Boltz2 retries exhausted")


def contacts(cif_path: Path, offset: int, preferred: set[int], risk: set[int]):
    st=MMCIFParser(QUIET=True).get_structure("pred",str(cif_path))
    model=next(st.get_models())
    target,binder=model["A"],model["B"]
    target_atoms=[
        atom for residue in target if residue.id[0]==" "
        for atom in residue if atom.element!="H"
    ]
    search=NeighborSearch(target_atoms)
    tc,bc=set(),set()
    for residue in binder:
        if residue.id[0]!=" ":
            continue
        for atom in residue:
            if atom.element=="H":
                continue
            close=search.search(atom.coord,5.0,level="A")
            if close:
                bc.add(residue.id[1])
            for x in close:
                tc.add(x.get_parent().id[1]+offset)
    tc=sorted(tc)
    return {
        "target_contacts_uniprot":tc,
        "binder_contact_positions":sorted(bc),
        "preferred_contacts":sorted(set(tc)&preferred),
        "risk_contacts":sorted(set(tc)&risk),
        "preferred_core_retained":bool(set(tc)&preferred),
    }


def main():
    key=(os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_API_KEY1")
         or os.getenv("NVIDIA_API_TOKEN") or os.getenv("NGC_API_KEY"))
    if not key:
        raise RuntimeError("Missing NVIDIA API credential")

    row,binder=read_selection()
    hm=human_mouse_map()
    mh={m:h for h,m in hm.items() if m is not None}
    preferred_m={hm[h] for h in PREFERRED_H if hm.get(h) is not None}
    risk_m={hm[h] for h in RISK_H if hm.get(h) is not None}

    targets={
        "human":{
            "sequence":fasta(HUMAN_FASTA),
            "offset":HUMAN_START-1,
            "preferred":PREFERRED_H,
            "risk":RISK_H,
        },
        "mouse":{
            "sequence":fasta(MOUSE_FASTA),
            "offset":MOUSE_START-1,
            "preferred":preferred_m,
            "risk":risk_m,
        },
    }

    OUT.mkdir(parents=True,exist_ok=True)
    result_all={
        "candidate_id":row.get("design_id") or "odin_selected",
        "binder_sequence":binder,
        "odin_selection":{
            "design_index":row.get("design_index"),
            "iteration":row.get("iteration"),
            "stage":row.get("stage"),
            "score":row.get("score"),
            "worst_target_i_ptm":row.get("worst_target_i_ptm"),
            "worst_target_context":row.get("worst_target_context"),
        },
        "validation_mode":"sequence_only_full_ectodomain_boltz2",
        "species":{},
    }

    for species,spec in targets.items():
        payload={
            "polymers":[
                {"id":"A","molecule_type":"protein","sequence":spec["sequence"]},
                {"id":"B","molecule_type":"protein","sequence":binder},
            ],
            "recycling_steps":3,
            "sampling_steps":50,
            "diffusion_samples":1,
            "step_scale":1.638,
            "output_format":"mmcif",
            "write_full_pae":False,
        }
        result=post_json(payload,key)
        structs=result.get("structures",[])
        if not structs:
            raise RuntimeError(f"No Boltz2 structure for {species}")
        spdir=OUT/species
        spdir.mkdir(parents=True,exist_ok=True)
        cif=spdir/"structure_1.cif"
        cif.write_text(structs[0]["structure"],encoding="utf-8")

        cm=contacts(cif,spec["offset"],spec["preferred"],spec["risk"])
        iptm=result.get("protein_iptm_scores",[None])[0]
        apical=any(189<=r<=383 if species=="human" else 191<=r<=385
                   for r in cm["target_contacts_uniprot"])
        pose_pass=bool(
            iptm is not None and float(iptm)>=0.50
            and cm["preferred_core_retained"] and apical
        )
        summary={
            "confidence":result.get("confidence_scores",[None])[0],
            "ptm":result.get("ptm_scores",[None])[0],
            "protein_iptm":iptm,
            "complex_plddt":result.get("complex_plddt_scores",[None])[0],
            "complex_iplddt":result.get("complex_iplddt_scores",[None])[0],
            "complex_pde":result.get("complex_pde_scores",[None])[0],
            "complex_ipde":result.get("complex_ipde_scores",[None])[0],
            "contacts":cm,
            "apical_contact":apical,
            "pose_gate":{
                "min_protein_iptm":0.50,
                "requires_preferred_contact":True,
                "pass":pose_pass,
            },
        }
        (spdir/"result_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
        result_all["species"][species]=summary
        time.sleep(20)

    human=result_all["species"]["human"]
    mouse=result_all["species"]["mouse"]
    hc=set(human["contacts"]["target_contacts_uniprot"])
    mc_as_h={mh[x] for x in mouse["contacts"]["target_contacts_uniprot"] if x in mh}
    overlap=sorted(hc&mc_as_h)
    union=hc|mc_as_h
    jaccard=len(overlap)/len(union) if union else 0.0
    same=len(overlap)>=3 and jaccard>=0.15
    result_all["cross_species_epitope"]={
        "human_mouse_homologous_overlap_human_numbering":overlap,
        "overlap_count":len(overlap),
        "jaccard":round(jaccard,4),
        "same_epitope_pass":same,
    }
    result_all["cross_species_pose_pass"]=bool(
        human["pose_gate"]["pass"] and mouse["pose_gate"]["pass"] and same
    )

    (OUT/"batch_summary.json").write_text(json.dumps(result_all,indent=2)+"\n")
    winners=[result_all["candidate_id"]] if result_all["cross_species_pose_pass"] else []
    (OUT/"pose_pass_ids.json").write_text(
        json.dumps({"pose_pass_ids":winners},indent=2)+"\n"
    )
    print(json.dumps({
        "candidate":result_all["candidate_id"],
        "odin_worst_target_i_ptm":result_all["odin_selection"]["worst_target_i_ptm"],
        "human_ipTM":human["protein_iptm"],
        "mouse_ipTM":mouse["protein_iptm"],
        "overlap":overlap,
        "jaccard":round(jaccard,4),
        "cross_species_pose_pass":result_all["cross_species_pose_pass"],
    },indent=2))


if __name__=="__main__":
    main()
