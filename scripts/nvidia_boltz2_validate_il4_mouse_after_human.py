#!/usr/bin/env python3
"""IL4-M: mouse and homologous-epitope validation for human-pass IL4 sequences.

Only sequences that already passed IL4-H are evaluated here. Mouse prediction
is template-free against the complete mouse TfR1 ectodomain. The final
cross-species gate requires:
- mouse protein iPTM >= 0.50;
- preferred/apical mouse contact;
- >=3 homologous human/mouse contacts;
- homologous-contact Jaccard >= 0.15.
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

ROOT=Path(__file__).resolve().parents[1]
HUMAN_BATCH=ROOT/"predictions/inference_loop_human_anchor_il4_human/batch_summary.json"
HUMAN_PASS=ROOT/"predictions/inference_loop_human_anchor_il4_human/human_pose_pass_ids.json"
MOUSE_FASTA=ROOT/"data/mouse_TfR1/mouse_TfR1_Q62351_123-763.fasta"
HUMAN_FASTA=ROOT/"data/human_TfR1/human_TfR1_P02786_121-760.fasta"
MASK_CSV=ROOT/"analysis/apical_design_mask.csv"
OUT=ROOT/"predictions/inference_loop_human_anchor_il4_cross_species"
URL="https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"

HUMAN_START=121
MOUSE_START=123
PREFERRED_H={206,207,209,214,231,343,344,347,350,374}
RISK_H={204,208,210,211,212,348,371}
OFFICIAL_OVERRIDES={204:207,208:211,210:213,211:214,212:None,348:350,371:373}


def fasta(path):
    return "".join(
        x.strip() for x in path.read_text().splitlines()
        if x.strip() and not x.startswith(">")
    )


def mapping():
    hs,ms=fasta(HUMAN_FASTA),fasta(MOUSE_FASTA)
    a=Align.PairwiseAligner();a.mode="global"
    a.substitution_matrix=substitution_matrices.load("BLOSUM62")
    a.open_gap_score=-10;a.extend_gap_score=-0.5
    aln=a.align(hs,ms)[0];c=aln.coordinates
    out={}
    for k in range(c.shape[1]-1):
        h0,h1=int(c[0,k]),int(c[0,k+1])
        m0,m1=int(c[1,k]),int(c[1,k+1])
        dh,dm=h1-h0,m1-m0
        if dh and dm:
            for j in range(min(dh,dm)):
                out[HUMAN_START+h0+j]=MOUSE_START+m0+j
        elif dh:
            for j in range(dh):
                out[HUMAN_START+h0+j]=None
    out.update(OFFICIAL_OVERRIDES)
    if MASK_CSV.exists():
        with MASK_CSV.open(newline="") as fh:
            for row in csv.DictReader(fh):
                try: h=int(row["human_residue"])
                except Exception: continue
                m=(row.get("mouse_residue") or "").strip()
                out[h]=int(m) if m else None
    return out


def post(payload,key,max_attempts=5):
    for attempt in range(1,max_attempts+1):
        req=Request(
            URL,data=json.dumps(payload).encode(),
            headers={"Content-Type":"application/json","Authorization":f"Bearer {key}"},
            method="POST",
        )
        try:
            with urlopen(req,timeout=1200) as r:
                return json.loads(r.read().decode())
        except HTTPError as exc:
            detail=exc.read().decode(errors="replace")
            if exc.code==429 and attempt<max_attempts:
                ra=exc.headers.get("Retry-After")
                time.sleep(int(ra) if ra and ra.isdigit() else 25*attempt)
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


def contacts(cif,preferred,risk):
    st=MMCIFParser(QUIET=True).get_structure("p",str(cif))
    model=next(st.get_models())
    target,binder=model["A"],model["B"]
    atoms=[a for r in target if r.id[0]==" " for a in r if a.element!="H"]
    ns=NeighborSearch(atoms)
    tc=set();bc=set()
    for r in binder:
        if r.id[0]!=" ": continue
        for a in r:
            if a.element=="H": continue
            close=ns.search(a.coord,5.0,level="A")
            if close: bc.add(r.id[1])
            for x in close:
                tc.add(x.get_parent().id[1]+MOUSE_START-1)
    tc=sorted(tc)
    apical=[x for x in tc if 191<=x<=385]
    preferred_contacts=sorted(set(tc)&preferred)
    return {
        "target_contacts_uniprot":tc,
        "binder_contact_positions":sorted(bc),
        "preferred_contacts":preferred_contacts,
        "risk_contacts":sorted(set(tc)&risk),
        "apical_contacts":apical,
        "apical_fraction":round(len(apical)/len(tc),4) if tc else 0.0,
        "preferred_core_retained":bool(preferred_contacts),
    }


def main():
    key=(os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_API_KEY1")
         or os.getenv("NVIDIA_API_TOKEN") or os.getenv("NGC_API_KEY"))
    if not key:
        raise RuntimeError("Missing NVIDIA API credential")
    if not HUMAN_BATCH.exists() or not HUMAN_PASS.exists():
        raise SystemExit("IL4-H outputs missing")

    human_data=json.loads(HUMAN_BATCH.read_text())
    pass_ids=json.loads(HUMAN_PASS.read_text()).get("human_pose_pass_ids",[])
    if not pass_ids:
        print(json.dumps({
            "status":"skipped",
            "reason":"No sequence-intrinsic human apical parent passed IL4-H",
        },indent=2))
        return

    by_id={x["candidate_id"]:x for x in human_data.get("candidates",[])}
    hm=mapping();mh={m:h for h,m in hm.items() if m is not None}
    preferred_m={hm[h] for h in PREFERRED_H if hm.get(h) is not None}
    risk_m={hm[h] for h in RISK_H if hm.get(h) is not None}
    mouse_seq=fasta(MOUSE_FASTA)

    OUT.mkdir(parents=True,exist_ok=True)
    results=[]

    for cid in pass_ids:
        hrow=by_id.get(cid)
        if not hrow:
            raise RuntimeError(f"Human pass ID missing from batch: {cid}")
        seq=hrow["binder_sequence"]
        payload={
            "polymers":[
                {"id":"A","molecule_type":"protein","sequence":mouse_seq},
                {"id":"B","molecule_type":"protein","sequence":seq},
            ],
            "recycling_steps":3,
            "sampling_steps":50,
            "diffusion_samples":1,
            "step_scale":1.638,
            "output_format":"mmcif",
            "write_full_pae":False,
        }
        data=post(payload,key)
        structs=data.get("structures",[])
        if not structs:
            raise RuntimeError(f"No mouse structure for {cid}")
        cdir=OUT/cid
        cdir.mkdir(parents=True,exist_ok=True)
        cif=cdir/"mouse_structure_1.cif"
        cif.write_text(structs[0]["structure"])

        cm=contacts(cif,preferred_m,risk_m)
        iptm=data.get("protein_iptm_scores",[None])[0]
        mouse_pose=bool(
            iptm is not None
            and float(iptm)>=0.50
            and cm["preferred_core_retained"]
            and bool(cm["apical_contacts"])
        )

        hc=set(hrow["human"]["contacts"]["target_contacts_uniprot"])
        mc_as_h={mh[x] for x in cm["target_contacts_uniprot"] if x in mh}
        overlap=sorted(hc&mc_as_h)
        union=hc|mc_as_h
        jac=len(overlap)/len(union) if union else 0.0
        same=len(overlap)>=3 and jac>=0.15
        cross=bool(mouse_pose and same)

        result={
            "candidate_id":cid,
            "binder_sequence":seq,
            "human":hrow["human"],
            "mouse":{
                "protein_iptm":iptm,
                "confidence":data.get("confidence_scores",[None])[0],
                "ptm":data.get("ptm_scores",[None])[0],
                "complex_plddt":data.get("complex_plddt_scores",[None])[0],
                "complex_iplddt":data.get("complex_iplddt_scores",[None])[0],
                "contacts":cm,
                "pose_gate_pass":mouse_pose,
            },
            "cross_species_epitope":{
                "overlap_human_numbering":overlap,
                "overlap_count":len(overlap),
                "jaccard":round(jac,4),
                "same_epitope_pass":same,
            },
            "cross_species_pose_pass":cross,
        }
        (cdir/"result_summary.json").write_text(json.dumps(result,indent=2)+"\n")
        results.append(result)
        time.sleep(15)

    results.sort(
        key=lambda x:(
            bool(x["cross_species_pose_pass"]),
            min(
                float(x["human"]["protein_iptm"] or 0),
                float(x["mouse"]["protein_iptm"] or 0),
            ),
            x["cross_species_epitope"]["jaccard"],
        ),
        reverse=True,
    )
    winners=[x["candidate_id"] for x in results if x["cross_species_pose_pass"]]
    summary={
        "candidates":results,
        "cross_species_pose_pass_ids":winners,
        "ranked":[{
            "candidate_id":x["candidate_id"],
            "human_ipTM":x["human"]["protein_iptm"],
            "mouse_ipTM":x["mouse"]["protein_iptm"],
            "mouse_preferred":x["mouse"]["contacts"]["preferred_contacts"],
            "mouse_apical_fraction":x["mouse"]["contacts"]["apical_fraction"],
            "overlap_count":x["cross_species_epitope"]["overlap_count"],
            "jaccard":x["cross_species_epitope"]["jaccard"],
            "pass":x["cross_species_pose_pass"],
        } for x in results],
    }
    (OUT/"batch_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (OUT/"cross_species_pose_pass_ids.json").write_text(
        json.dumps({"cross_species_pose_pass_ids":winners},indent=2)+"\n"
    )
    print(json.dumps(summary["ranked"],indent=2))


if __name__=="__main__":
    main()
