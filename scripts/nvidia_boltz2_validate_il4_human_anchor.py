#!/usr/bin/env python3
"""IL4-H: template-free full-human Boltz-2 screen for human-anchor sequences.

A candidate passes only if the sequence itself, without any structural template,
recovers a full-human TfR1 pose with:
- protein iPTM >= 0.50;
- at least one preferred apical contact;
- at least one contact inside the 189-383 apical domain.

This is the first true sequence-intrinsic human parent gate.
"""

from __future__ import annotations

import csv
import json
import os
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from Bio.PDB import MMCIFParser, NeighborSearch

ROOT=Path(__file__).resolve().parents[1]
CANDIDATES=ROOT/"designs/inference_loop/human_anchor_il4/pose_screen_candidates.csv"
HUMAN_FASTA=ROOT/"data/human_TfR1/human_TfR1_P02786_121-760.fasta"
OUT=ROOT/"predictions/inference_loop_human_anchor_il4_human"
URL="https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"

HUMAN_START=121
PREFERRED={206,207,209,214,231,343,344,347,350,374}
RISK={204,208,210,211,212,348,371}


def fasta(path):
    return "".join(
        x.strip() for x in path.read_text().splitlines()
        if x.strip() and not x.startswith(">")
    )


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


def contacts(cif):
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
                tc.add(x.get_parent().id[1]+HUMAN_START-1)
    tc=sorted(tc)
    apical=[x for x in tc if 189<=x<=383]
    preferred=sorted(set(tc)&PREFERRED)
    return {
        "target_contacts_uniprot":tc,
        "binder_contact_positions":sorted(bc),
        "preferred_contacts":preferred,
        "risk_contacts":sorted(set(tc)&RISK),
        "apical_contacts":apical,
        "apical_fraction":round(len(apical)/len(tc),4) if tc else 0.0,
        "preferred_core_retained":bool(preferred),
    }


def main():
    key=(os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_API_KEY1")
         or os.getenv("NVIDIA_API_TOKEN") or os.getenv("NGC_API_KEY"))
    if not key:
        raise RuntimeError("Missing NVIDIA API credential")
    if not CANDIDATES.exists():
        raise SystemExit(f"Missing {CANDIDATES}")

    with CANDIDATES.open(newline="") as fh:
        rows=list(csv.DictReader(fh))
    rows=[
        r for r in rows
        if str(r.get("prefilter_pass","")).lower() in {"true","1","yes"}
    ][:12]
    if not rows:
        raise SystemExit("No IL4 human candidates")

    target_seq=fasta(HUMAN_FASTA)
    OUT.mkdir(parents=True,exist_ok=True)
    results=[]

    for row in rows:
        cid=row["candidate_id"];seq=row["sequence"]
        payload={
            "polymers":[
                {"id":"A","molecule_type":"protein","sequence":target_seq},
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
            raise RuntimeError(f"No structure for {cid}")
        cdir=OUT/cid
        cdir.mkdir(parents=True,exist_ok=True)
        cif=cdir/"structure_1.cif"
        cif.write_text(structs[0]["structure"])

        cm=contacts(cif)
        iptm=data.get("protein_iptm_scores",[None])[0]
        pose=bool(
            iptm is not None
            and float(iptm)>=0.50
            and cm["preferred_core_retained"]
            and bool(cm["apical_contacts"])
        )
        score=(
            2.0*float(iptm or 0)
            +0.80*(1.0 if cm["preferred_core_retained"] else 0.0)
            +0.50*cm["apical_fraction"]
            +0.10*min(len(cm["preferred_contacts"]),3)
        )
        result={
            "candidate_id":cid,
            "binder_sequence":seq,
            "source":row,
            "human":{
                "protein_iptm":iptm,
                "confidence":data.get("confidence_scores",[None])[0],
                "ptm":data.get("ptm_scores",[None])[0],
                "complex_plddt":data.get("complex_plddt_scores",[None])[0],
                "complex_iplddt":data.get("complex_iplddt_scores",[None])[0],
                "contacts":cm,
                "pose_gate_pass":pose,
            },
            "human_anchor_score":round(score,4),
        }
        (cdir/"result_summary.json").write_text(json.dumps(result,indent=2)+"\n")
        results.append(result)
        time.sleep(15)

    results.sort(
        key=lambda x:(
            bool(x["human"]["pose_gate_pass"]),
            x["human_anchor_score"],
        ),
        reverse=True,
    )
    winners=[x["candidate_id"] for x in results if x["human"]["pose_gate_pass"]]
    summary={
        "candidates":results,
        "ranked":[{
            "candidate_id":x["candidate_id"],
            "human_ipTM":x["human"]["protein_iptm"],
            "preferred_contacts":x["human"]["contacts"]["preferred_contacts"],
            "apical_fraction":x["human"]["contacts"]["apical_fraction"],
            "human_anchor_score":x["human_anchor_score"],
            "pass":x["human"]["pose_gate_pass"],
        } for x in results],
        "human_pose_pass_ids":winners,
        "best_partial_id":results[0]["candidate_id"] if results else None,
    }
    (OUT/"batch_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (OUT/"human_pose_pass_ids.json").write_text(
        json.dumps({"human_pose_pass_ids":winners},indent=2)+"\n"
    )
    print(json.dumps(summary["ranked"],indent=2))


if __name__=="__main__":
    main()
