#!/usr/bin/env python3
"""Inference Loop R2: template-free full human/mouse Boltz-2 validation.

Each R1 ProteinMPNN sequence is predicted from sequence only against the complete
human and mouse TfR1 ectodomains. The proposal-state structures used by ProteinMPNN
are NOT supplied as templates.

Strict parent gate:
- human protein iPTM >= 0.50;
- mouse protein iPTM >= 0.50;
- preferred/apical contact retained in both species;
- >=3 homologous contact positions and Jaccard >=0.15.

A continuous search_score is also reported so that a failed R2 can seed the next iteration using
the best full-context partial solution rather than restarting blindly.
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
CANDIDATES=ROOT/"designs/inference_loop/full_context_r2/pose_screen_candidates.csv"
HUMAN_FASTA=ROOT/"data/human_TfR1/human_TfR1_P02786_121-760.fasta"
MOUSE_FASTA=ROOT/"data/mouse_TfR1/mouse_TfR1_Q62351_123-763.fasta"
MASK_CSV=ROOT/"analysis/apical_design_mask.csv"
OUT=ROOT/"predictions/inference_loop_full_context_r2"
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
        h0,h1=int(c[0,k]),int(c[0,k+1]);m0,m1=int(c[1,k]),int(c[1,k+1])
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
            method="POST"
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
                exc.code==422 and any(x in detail for x in ("InfrastructureError","AcceleratorError","CUDA error"))
            )
            if transient and attempt<max_attempts:
                time.sleep(25*attempt);continue
            raise RuntimeError(f"Boltz2 HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("Boltz2 retries exhausted")


def contacts(cif,offset,preferred,risk,apical_bounds):
    st=MMCIFParser(QUIET=True).get_structure("p",str(cif))
    model=next(st.get_models());target,binder=model["A"],model["B"]
    atoms=[a for r in target if r.id[0]==" " for a in r if a.element!="H"]
    ns=NeighborSearch(atoms);tc=set();bc=set()
    for r in binder:
        if r.id[0]!=" ": continue
        for a in r:
            if a.element=="H": continue
            close=ns.search(a.coord,5.0,level="A")
            if close: bc.add(r.id[1])
            for x in close: tc.add(x.get_parent().id[1]+offset)
    tc=sorted(tc)
    apical=[x for x in tc if apical_bounds[0]<=x<=apical_bounds[1]]
    return {
        "target_contacts_uniprot":tc,
        "binder_contact_positions":sorted(bc),
        "preferred_contacts":sorted(set(tc)&preferred),
        "risk_contacts":sorted(set(tc)&risk),
        "apical_contacts":apical,
        "apical_fraction":round(len(apical)/len(tc),4) if tc else 0.0,
        "preferred_core_retained":bool(set(tc)&preferred),
    }


def main():
    key=(os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_API_KEY1")
         or os.getenv("NVIDIA_API_TOKEN") or os.getenv("NGC_API_KEY"))
    if not key: raise RuntimeError("Missing NVIDIA API credential")
    if not CANDIDATES.exists(): raise SystemExit(f"Missing {CANDIDATES}")

    with CANDIDATES.open(newline="") as fh:
        rows=list(csv.DictReader(fh))
    rows=[r for r in rows if str(r.get("prefilter_pass","")).lower() in {"true","1","yes"}][:4]
    if not rows: raise SystemExit("No R2 candidates")

    hm=mapping();mh={m:h for h,m in hm.items() if m is not None}
    preferred_m={hm[h] for h in PREFERRED_H if hm.get(h) is not None}
    risk_m={hm[h] for h in RISK_H if hm.get(h) is not None}
    targets={
        "human":dict(sequence=fasta(HUMAN_FASTA),offset=HUMAN_START-1,preferred=PREFERRED_H,risk=RISK_H,bounds=(189,383)),
        "mouse":dict(sequence=fasta(MOUSE_FASTA),offset=MOUSE_START-1,preferred=preferred_m,risk=risk_m,bounds=(191,385)),
    }

    OUT.mkdir(parents=True,exist_ok=True)
    results=[]
    for row in rows:
        cid=row["candidate_id"];seq=row["sequence"];cres={"candidate_id":cid,"binder_sequence":seq,"source":row,"species":{}}
        for species,spec in targets.items():
            payload={
                "polymers":[
                    {"id":"A","molecule_type":"protein","sequence":spec["sequence"]},
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
            if not structs: raise RuntimeError(f"No structure for {cid}/{species}")
            sdir=OUT/cid/species;sdir.mkdir(parents=True,exist_ok=True)
            cif=sdir/"structure_1.cif";cif.write_text(structs[0]["structure"])
            cm=contacts(cif,spec["offset"],spec["preferred"],spec["risk"],spec["bounds"])
            iptm=data.get("protein_iptm_scores",[None])[0]
            pose=bool(iptm is not None and float(iptm)>=0.50 and cm["preferred_core_retained"] and cm["apical_contacts"])
            summ={
                "protein_iptm":iptm,
                "confidence":data.get("confidence_scores",[None])[0],
                "ptm":data.get("ptm_scores",[None])[0],
                "complex_plddt":data.get("complex_plddt_scores",[None])[0],
                "complex_iplddt":data.get("complex_iplddt_scores",[None])[0],
                "contacts":cm,
                "pose_gate_pass":pose,
            }
            (sdir/"result_summary.json").write_text(json.dumps(summ,indent=2)+"\n")
            cres["species"][species]=summ
            time.sleep(15)

        h=cres["species"]["human"];m=cres["species"]["mouse"]
        hc=set(h["contacts"]["target_contacts_uniprot"])
        mc_as_h={mh[x] for x in m["contacts"]["target_contacts_uniprot"] if x in mh}
        overlap=sorted(hc&mc_as_h);union=hc|mc_as_h
        jac=len(overlap)/len(union) if union else 0.0
        same=len(overlap)>=3 and jac>=0.15
        fullpass=bool(h["pose_gate_pass"] and m["pose_gate_pass"] and same)

        hi=float(h["protein_iptm"] or 0);mi=float(m["protein_iptm"] or 0)
        hp=1.0 if h["contacts"]["preferred_core_retained"] else 0.0
        mp=1.0 if m["contacts"]["preferred_core_retained"] else 0.0
        haf=h["contacts"]["apical_fraction"];maf=m["contacts"]["apical_fraction"]
        # continuous objective for next-loop selection; exact pass still uses strict gate.
        search_score=(
            2.0*min(hi,mi)
            +0.35*(hp+mp)
            +0.30*(haf+maf)
            +0.75*jac
            +0.15*min(len(overlap),5)
        )
        cres["cross_species_epitope"]={
            "overlap_human_numbering":overlap,
            "overlap_count":len(overlap),
            "jaccard":round(jac,4),
            "same_epitope_pass":same,
        }
        cres["cross_species_pose_pass"]=fullpass
        cres["search_score"]=round(search_score,4)
        results.append(cres)

    results.sort(key=lambda x:(bool(x["cross_species_pose_pass"]),x["search_score"]),reverse=True)
    winners=[x["candidate_id"] for x in results if x["cross_species_pose_pass"]]
    summary={
        "candidates":results,
        "ranked":[
            {
                "candidate_id":x["candidate_id"],
                "human_ipTM":x["species"]["human"]["protein_iptm"],
                "mouse_ipTM":x["species"]["mouse"]["protein_iptm"],
                "human_preferred":x["species"]["human"]["contacts"]["preferred_contacts"],
                "mouse_preferred":x["species"]["mouse"]["contacts"]["preferred_contacts"],
                "human_apical_fraction":x["species"]["human"]["contacts"]["apical_fraction"],
                "mouse_apical_fraction":x["species"]["mouse"]["contacts"]["apical_fraction"],
                "overlap_count":x["cross_species_epitope"]["overlap_count"],
                "jaccard":x["cross_species_epitope"]["jaccard"],
                "search_score":x["search_score"],
                "pass":x["cross_species_pose_pass"],
            } for x in results
        ],
        "pose_pass_ids":winners,
        "best_partial_id":results[0]["candidate_id"] if results else None,
    }
    (OUT/"batch_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (OUT/"pose_pass_ids.json").write_text(json.dumps({"pose_pass_ids":winners},indent=2)+"\n")
    print(json.dumps(summary["ranked"],indent=2))


if __name__=="__main__":
    main()
