#!/usr/bin/env python3
"""Phase 5B: neutral-pH structural retention screen for histidine variants.

Input:
  designs/pH_switch/odin_full_context_v4_histidine_scan/histidine_variants.csv

The top mechanistically ranked variants are independently predicted against the
full human and mouse TfR1 ectodomains with Boltz-2. This stage does NOT model pH.
It only asks whether a proposed His switch preserves the neutral-pH cross-species
parent geometry strongly enough to remain worth evaluating.

Strict neutral gate:
- human protein iPTM >= 0.50
- mouse protein iPTM >= 0.50
- preferred/apical contact retained in both species
- homologous contact overlap >=3 and Jaccard >=0.15
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
VARIANTS=ROOT/"designs/pH_switch/odin_full_context_v4_histidine_scan/histidine_variants.csv"
HUMAN_FASTA=ROOT/"data/human_TfR1/human_TfR1_P02786_121-760.fasta"
MOUSE_FASTA=ROOT/"data/mouse_TfR1/mouse_TfR1_Q62351_123-763.fasta"
MASK_CSV=ROOT/"analysis/apical_design_mask.csv"
OUT=ROOT/"predictions/odin_v4_histidine_neutral_screen"
URL="https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"

HUMAN_START=121
MOUSE_START=123
MAX_VARIANTS=8

PREFERRED_H={206,207,209,214,231,343,344,347,350,374}
RISK_H={204,208,210,211,212,348,371}

OFFICIAL_OVERRIDES={
    204:207,
    208:211,
    210:213,
    211:214,
    212:None,
    348:350,
    371:373,
}


def fasta(path:Path)->str:
    return "".join(
        x.strip() for x in path.read_text().splitlines()
        if x.strip() and not x.startswith(">")
    )


def hm_map():
    hs,ms=fasta(HUMAN_FASTA),fasta(MOUSE_FASTA)
    a=Align.PairwiseAligner()
    a.mode="global"
    a.substitution_matrix=substitution_matrices.load("BLOSUM62")
    a.open_gap_score=-10
    a.extend_gap_score=-0.5
    aln=a.align(hs,ms)[0]
    c=aln.coordinates
    mapping={}
    for k in range(c.shape[1]-1):
        h0,h1=int(c[0,k]),int(c[0,k+1])
        m0,m1=int(c[1,k]),int(c[1,k+1])
        dh,dm=h1-h0,m1-m0
        if dh and dm:
            for j in range(min(dh,dm)):
                mapping[HUMAN_START+h0+j]=MOUSE_START+m0+j
        elif dh:
            for j in range(dh):
                mapping[HUMAN_START+h0+j]=None
    mapping.update(OFFICIAL_OVERRIDES)
    if MASK_CSV.exists():
        with MASK_CSV.open(newline="") as fh:
            for row in csv.DictReader(fh):
                try:
                    h=int(row["human_residue"])
                except Exception:
                    continue
                m=(row.get("mouse_residue") or "").strip()
                mapping[h]=int(m) if m else None
    return mapping


def post(payload,key,max_attempts=5):
    for attempt in range(1,max_attempts+1):
        req=Request(
            URL,
            data=json.dumps(payload).encode(),
            headers={"Content-Type":"application/json","Authorization":f"Bearer {key}"},
            method="POST",
        )
        try:
            with urlopen(req,timeout=1200) as resp:
                return json.loads(resp.read().decode())
        except HTTPError as exc:
            detail=exc.read().decode(errors="replace")
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


def contact_metrics(cif_path:Path,offset:int,preferred:set[int],risk:set[int]):
    st=MMCIFParser(QUIET=True).get_structure("pred",str(cif_path))
    model=next(st.get_models())
    target,binder=model["A"],model["B"]
    ta=[
        atom for r in target if r.id[0]==" "
        for atom in r if atom.element!="H"
    ]
    search=NeighborSearch(ta)
    tc,bc=set(),set()
    for r in binder:
        if r.id[0]!=" ":
            continue
        for atom in r:
            if atom.element=="H":
                continue
            close=search.search(atom.coord,5.0,level="A")
            if close:
                bc.add(r.id[1])
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


def parse_variants():
    if not VARIANTS.exists():
        raise SystemExit(f"Missing Phase 5A variants: {VARIANTS}")
    with VARIANTS.open(newline="") as fh:
        rows=list(csv.DictReader(fh))
    rows.sort(key=lambda r:float(r.get("priority_score") or "-999"),reverse=True)
    return rows[:MAX_VARIANTS]


def main():
    key=(os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_API_KEY1")
         or os.getenv("NVIDIA_API_TOKEN") or os.getenv("NGC_API_KEY"))
    if not key:
        raise RuntimeError("Missing NVIDIA API credential")

    variants=parse_variants()
    if not variants:
        raise SystemExit("Phase 5A produced no histidine variants")

    mapping=hm_map()
    reverse={m:h for h,m in mapping.items() if m is not None}
    preferred_m={mapping[h] for h in PREFERRED_H if mapping.get(h) is not None}
    risk_m={mapping[h] for h in RISK_H if mapping.get(h) is not None}

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
    results=[]

    for idx,row in enumerate(variants):
        vid=row["variant_id"]
        binder=row["sequence"].strip()
        vdir=OUT/vid
        result={
            "variant_id":vid,
            "parent_id":row.get("parent_id"),
            "sequence":binder,
            "positions":[int(x) for x in (row.get("positions") or "").split(";") if x],
            "mutations":[x for x in (row.get("mutations") or "").split(";") if x],
            "priority_score":float(row.get("priority_score") or 0),
            "mechanism":row.get("mechanism"),
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
            data=post(payload,key)
            structs=data.get("structures",[])
            if not structs:
                raise RuntimeError(f"No Boltz2 structure: {vid} {species}")
            sdir=vdir/species
            sdir.mkdir(parents=True,exist_ok=True)
            cif=sdir/"structure_1.cif"
            cif.write_text(structs[0]["structure"])

            cm=contact_metrics(
                cif,spec["offset"],spec["preferred"],spec["risk"]
            )
            iptm=data.get("protein_iptm_scores",[None])[0]
            apical_range=(189,383) if species=="human" else (191,385)
            apical=any(apical_range[0]<=r<=apical_range[1] for r in cm["target_contacts_uniprot"])
            pose_pass=bool(
                iptm is not None and float(iptm)>=0.50
                and cm["preferred_core_retained"] and apical
            )
            summary={
                "protein_iptm":iptm,
                "confidence":data.get("confidence_scores",[None])[0],
                "ptm":data.get("ptm_scores",[None])[0],
                "complex_plddt":data.get("complex_plddt_scores",[None])[0],
                "complex_iplddt":data.get("complex_iplddt_scores",[None])[0],
                "contacts":cm,
                "apical_contact":apical,
                "pose_gate_pass":pose_pass,
            }
            (sdir/"result_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
            result["species"][species]=summary
            time.sleep(15)

        hc=set(result["species"]["human"]["contacts"]["target_contacts_uniprot"])
        mc_as_h={
            reverse[x]
            for x in result["species"]["mouse"]["contacts"]["target_contacts_uniprot"]
            if x in reverse
        }
        overlap=sorted(hc&mc_as_h)
        union=hc|mc_as_h
        jaccard=len(overlap)/len(union) if union else 0.0
        same=len(overlap)>=3 and jaccard>=0.15
        cross=bool(
            result["species"]["human"]["pose_gate_pass"]
            and result["species"]["mouse"]["pose_gate_pass"]
            and same
        )
        result["cross_species_epitope"]={
            "overlap_human_numbering":overlap,
            "overlap_count":len(overlap),
            "jaccard":round(jaccard,4),
            "same_epitope_pass":same,
        }
        result["neutral_cross_species_pose_pass"]=cross
        results.append(result)

    results.sort(
        key=lambda r:(
            bool(r["neutral_cross_species_pose_pass"]),
            min(
                float(r["species"]["human"]["protein_iptm"] or 0),
                float(r["species"]["mouse"]["protein_iptm"] or 0),
            ),
            r["cross_species_epitope"]["jaccard"],
            r["priority_score"],
        ),
        reverse=True,
    )

    winners=[r["variant_id"] for r in results if r["neutral_cross_species_pose_pass"]]
    summary={
        "screened_variant_count":len(results),
        "neutral_pose_pass_ids":winners,
        "ranked":[{
            "variant_id":r["variant_id"],
            "human_ipTM":r["species"]["human"]["protein_iptm"],
            "mouse_ipTM":r["species"]["mouse"]["protein_iptm"],
            "overlap_count":r["cross_species_epitope"]["overlap_count"],
            "jaccard":r["cross_species_epitope"]["jaccard"],
            "neutral_pass":r["neutral_cross_species_pose_pass"],
            "priority_score":r["priority_score"],
        } for r in results],
        "warning":"This is a neutral-pH structural retention screen only; it does not establish pH-selective release.",
    }
    (OUT/"batch_summary.json").write_text(json.dumps({"candidates":results},indent=2)+"\n")
    (OUT/"neutral_pose_pass_ids.json").write_text(
        json.dumps({"neutral_pose_pass_ids":winners},indent=2)+"\n"
    )
    (OUT/"screen_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
