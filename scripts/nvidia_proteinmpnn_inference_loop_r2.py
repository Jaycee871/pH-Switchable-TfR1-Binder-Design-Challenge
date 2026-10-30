#!/usr/bin/env python3
"""Inference Loop R2: tied ProteinMPNN over two full-context inferred states.

State H (human):
  Independently predicted full-human TfR1 + v6 binder.
  Protein iPTM 0.690, preferred contacts 231/347/350.

State M (mouse):
  Independently predicted full-mouse TfR1 + Odin V3 early-logit F39 binder.
  Protein iPTM 0.223, but it physically contacts mouse 345/346/349/352.

The two full complexes are placed far apart in one ProteinMPNN input. Binder
chains B and D are tied position-by-position, forcing one shared 60-aa sequence.
Only the UNION of the two observed full-context binder interfaces is designable.

This is full-context inference-in-the-loop sequence design, not pH engineering.
No new histidines or cysteines are allowed at designable positions.
"""

from __future__ import annotations

import csv
import json
import os
import re
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from Bio.PDB import MMCIFParser, PDBIO

ROOT=Path(__file__).resolve().parents[1]
HUMAN_CIF=ROOT/(
    "predictions/freebindcraft_v6_gate_override_pose_screen/"
    "TFR1_FBC343_V6_l60_s289817/human/structure_1.cif"
)
MOUSE_CIF=ROOT/(
    "predictions/odin_v3_early_logit_pose_screen/rank_0/"
    "mouse/structure_1.cif"
)
OUT=ROOT/"designs/inference_loop/full_context_r2"
URL="https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict"

LENGTH=60
TRANSLATE_X=260.0

PARENT="SPSVISTTTTTYTYHGKQVTRIRTTINPGEKIEPAPGERQMRGRDWETGYIIVHTYKVNS"
F39="SPSVISTTTETYTPHGKQVTRIRTTINPGEKIEPAPGERQMPGRDADTGYTVVHTYKVNS"

# R1 showed that redesigning the full 29-position interface union destroyed the
# validated human basin. R2 therefore permits changes ONLY at the seven positions
# where the good-human v6 parent and apical-touching mouse F39 sequence differ.
RECONCILIATION_POSITIONS={10,14,42,46,47,51,52}
HUMAN_INTERFACE={8,12,14,15,19,21,23,39,41,42,43,44,46,47,52,54,56,58}
MOUSE_INTERFACE={10,12,14,15,17,19,21,23,28,29,31,33,34,38,39,40,41,42,43,44,45,46,47,51,52,54,56,58}
DESIGN_POSITIONS=sorted(RECONCILIATION_POSITIONS)
FIXED_POSITIONS=[i for i in range(1,LENGTH+1) if i not in DESIGN_POSITIONS]

ALPHABET="ACDEFGHIKLMNPQRSTVWYX"
AA3={
    "A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","E":"GLU","Q":"GLN",
    "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
    "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL",
}


def post(payload,key,max_attempts=5):
    for attempt in range(1,max_attempts+1):
        req=Request(
            URL,
            data=json.dumps(payload).encode(),
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
                delay=int(ra) if ra and ra.isdigit() else 20*attempt
                time.sleep(delay)
                continue
            if exc.code in {500,502,503,504} and attempt<max_attempts:
                time.sleep(20*attempt)
                continue
            raise RuntimeError(f"ProteinMPNN HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("ProteinMPNN retries exhausted")


def cif_to_pdb_lines(cif:Path):
    st=MMCIFParser(QUIET=True).get_structure(cif.stem,str(cif))
    tmp=OUT/(cif.stem+"_converted.pdb")
    io=PDBIO()
    io.set_structure(st)
    io.save(str(tmp))
    return [
        x for x in tmp.read_text(errors="replace").splitlines()
        if x.startswith("ATOM") and not x[12:16].strip().startswith("H")
    ]


def rewrite(line,chain,dx=0.0,resname=None):
    chars=list(line.rstrip("\n"))
    if len(chars)<80:
        chars.extend(" "*(80-len(chars)))
    chars[21]=chain
    if dx:
        x=float(line[30:38])+dx
        chars[30:38]=list(f"{x:8.3f}")
    if resname is not None:
        chars[17:20]=list(f"{resname:>3}")
    return "".join(chars).rstrip()


def binder_residue_order(lines,chain="B"):
    order=[];seen=set()
    for line in lines:
        if line[21].strip()!=chain:
            continue
        rid=int(line[22:26])
        if rid not in seen:
            seen.add(rid);order.append(rid)
    return order


def build_dual_state():
    OUT.mkdir(parents=True,exist_ok=True)
    h=cif_to_pdb_lines(HUMAN_CIF)
    m=cif_to_pdb_lines(MOUSE_CIF)

    h_binder=binder_residue_order(h,"B")
    m_binder=binder_residue_order(m,"B")
    if len(h_binder)!=LENGTH or len(m_binder)!=LENGTH:
        raise RuntimeError(
            f"Expected 60-aa binders; human={len(h_binder)} mouse={len(m_binder)}"
        )
    h_index={r:i+1 for i,r in enumerate(h_binder)}
    m_index={r:i+1 for i,r in enumerate(m_binder)}

    out=[]
    for line in h:
        ch=line[21].strip()
        if ch=="A":
            out.append(rewrite(line,"A"))
        elif ch=="B":
            pos=h_index[int(line[22:26])]
            out.append(rewrite(line,"B",resname=AA3[PARENT[pos-1]]))
    out.append("TER")

    for line in m:
        ch=line[21].strip()
        if ch=="A":
            out.append(rewrite(line,"C",TRANSLATE_X))
        elif ch=="B":
            pos=m_index[int(line[22:26])]
            # Coordinates come from the F39 full-context inference state, but the
            # design chain identity is normalized to the human-v6 parent before
            # tied redesign. ProteinMPNN then optimizes one sequence over both.
            out.append(rewrite(line,"D",TRANSLATE_X,resname=AA3[PARENT[pos-1]]))
    out.extend(["TER","END"])
    text="\n".join(out)+"\n"
    (OUT/"dual_full_context_states.pdb").write_text(text)
    return text


def parse_fasta(text):
    out=[];h=None;s=[]
    for raw in text.splitlines():
        line=raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if h is not None:
                out.append((h,"".join(s)))
            h=line[1:];s=[]
        else:
            s.append(line)
    if h is not None:
        out.append((h,"".join(s)))
    return out


def header_value(header,pat):
    m=re.search(pat,header)
    return m.group(1).strip() if m else ""


def binder_segment(sequence):
    parts=[x for x in sequence.split("/") if x]
    sixty=[x for x in parts if len(x)==LENGTH]
    if len(sixty)>=2:
        if len(set(sixty))!=1:
            raise RuntimeError(f"Tied B/D sequences differ: {sixty}")
        return sixty[0]
    if len(sixty)==1:
        return sixty[0]
    raise RuntimeError(f"Could not identify tied binder in {[len(x) for x in parts]}")


def longest_run(seq):
    best=run=1
    for i in range(1,len(seq)):
        run=run+1 if seq[i]==seq[i-1] else 1
        best=max(best,run)
    return best


def bias_matrix():
    """Strongly favor the two empirically observed identities at each reconciliation site."""
    mat=[[0.0]*len(ALPHABET) for _ in range(LENGTH)]
    for pos in DESIGN_POSITIONS:
        pa=PARENT[pos-1]
        fa=F39[pos-1]
        # Preserve the validated human basin unless the mouse-derived identity is
        # structurally favored in both full-context states.
        mat[pos-1][ALPHABET.index(pa)]+=2.0
        if fa!=pa:
            mat[pos-1][ALPHABET.index(fa)]+=1.25
        # Mildly discourage unrelated chemistry at these already-informative sites.
        for aa in ALPHABET:
            if aa not in {pa,fa,"X"}:
                mat[pos-1][ALPHABET.index(aa)]-=0.35
    return mat


def main():
    key=(os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_API_KEY1")
         or os.getenv("NVIDIA_API_TOKEN") or os.getenv("NGC_API_KEY"))
    if not key:
        raise RuntimeError("Missing NVIDIA API credential")
    if not HUMAN_CIF.exists() or not MOUSE_CIF.exists():
        raise RuntimeError(
            f"Missing inference state: human={HUMAN_CIF.exists()} mouse={MOUSE_CIF.exists()}"
        )

    dual=build_dual_state()

    fixed=json.dumps({
        "input":{"B":FIXED_POSITIONS,"D":FIXED_POSITIONS}
    })
    tied=json.dumps({
        "input":[{"B":[i],"D":[i]} for i in range(1,LENGTH+1)]
    })
    bias=bias_matrix()
    bias_json=json.dumps({"input":{"B":bias,"D":bias}})

    payload={
        "input_pdb":dual,
        "input_pdb_chains":["B","D"],
        "num_seq_per_target":40,
        "sampling_temp":[0.08,0.15,0.22],
        "use_soluble_model":False,
        "ca_only":False,
        "omit_AAs":["C","H"],
        "fixed_positions_jsonl":fixed,
        "tied_positions_jsonl":tied,
        "bias_by_res_jsonl":bias_json,
        "random_seed":20261109,
    }

    (OUT/"design_spec.json").write_text(json.dumps({
        "method":"sparse full-context inference-loop reconciliation",
        "human_state":"Phase 4Q full-human v6 prediction, ipTM 0.690",
        "mouse_state":"Phase 4T F39 full-mouse prediction, ipTM 0.223; contacts 345/346/349/352",
        "parent_sequence":PARENT,
        "mouse_source_sequence":F39,
        "design_positions":DESIGN_POSITIONS,
        "fixed_positions":FIXED_POSITIONS,
        "human_interface_positions":sorted(HUMAN_INTERFACE),
        "mouse_interface_positions":sorted(MOUSE_INTERFACE),
        "tied_positions":"B_i == D_i for i=1..60",
        "new_histidine_allowed":False,
        "new_cysteine_allowed":False,
        "reconciliation_positions":sorted(RECONCILIATION_POSITIONS),
        "note":"R2 allows only the seven v6/F39 sequence differences to change; 53/60 positions are fixed.",
    },indent=2)+"\n")

    result=post(payload,key)
    mfasta=result.get("mfasta","")
    if not mfasta.strip():
        raise RuntimeError("ProteinMPNN returned no mfasta")
    (OUT/"designed_sequences.fasta").write_text(mfasta)

    records=[
        (h,s) for h,s in parse_fasta(mfasta)
        if not h.lower().startswith("input")
        and "native" not in h.lower()
        and not h.lower().startswith("wt")
    ]

    rows=[];seen=set()
    parent_run=longest_run(PARENT)
    for i,(header,seqout) in enumerate(records,1):
        seq=binder_segment(seqout)
        if seq in seen:
            continue
        seen.add(seq)
        changed=[p for p in range(1,LENGTH+1) if seq[p-1]!=PARENT[p-1]]
        outside=[p for p in changed if p not in DESIGN_POSITIONS]
        new_h=[p for p in changed if seq[p-1]=="H" and PARENT[p-1]!="H"]
        new_c=[p for p in changed if seq[p-1]=="C" and PARENT[p-1]!="C"]
        flags=[]
        if outside: flags.append("mutation_outside_full_context_interface")
        if new_h: flags.append("new_histidine_before_pH_stage")
        if new_c: flags.append("new_cysteine")
        if longest_run(seq)>parent_run: flags.append("homopolymer_worse_than_parent")
        charge=(seq.count("D")+seq.count("E")+seq.count("K")+seq.count("R"))/LENGTH
        if charge>0.50: flags.append("charge_fraction_gt_0.50")

        mpnn=header_value(header,r"(?:^|,\s*)score=([^,]+)")
        gscore=header_value(header,r"(?:^|,\s*)global_score=([^,]+)")
        temp=header_value(header,r"(?:^|,\s*)T=([^,]+)")
        rows.append({
            "candidate_id":f"TFR1_ILOOP_R2_{i:02d}",
            "sequence":seq,
            "mpnn_score":mpnn,
            "global_score":gscore,
            "temperature":temp,
            "changed_positions":";".join(map(str,changed)),
            "changed_count":len(changed),
            "distance_from_f39":sum(a!=b for a,b in zip(seq,F39)),
            "charge_fraction":round(charge,4),
            "prefilter_flags":";".join(flags),
            "prefilter_pass":not flags,
            "header":header,
        })

    if not rows:
        raise RuntimeError("No designed sequences parsed")

    fields=list(rows[0].keys())
    with (OUT/"all_candidates.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields);w.writeheader();w.writerows(rows)

    passing=[x for x in rows if x["prefilter_pass"]]
    # Preserve human state by mildly penalizing excessive mutation count.
    def rank(x):
        score=float(x["mpnn_score"]) if x["mpnn_score"] else 999.0
        return (score+0.06*x["changed_count"],x["changed_count"])
    passing.sort(key=rank)

    # Diversity-aware top 4: do not select exact near-clones when alternatives exist.
    chosen=[]
    for row in passing:
        seq=row["sequence"]
        if chosen and min(sum(a!=b for a,b in zip(seq,x["sequence"])) for x in chosen)<2:
            continue
        chosen.append(row)
        if len(chosen)>=4:
            break
    if len(chosen)<4:
        for row in passing:
            if row not in chosen:
                chosen.append(row)
            if len(chosen)>=4:
                break

    with (OUT/"pose_screen_candidates.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields);w.writeheader();w.writerows(chosen)

    summary={
        "candidate_count":len(rows),
        "unique_candidate_count":len(seen),
        "prefilter_pass_count":len(passing),
        "pose_screen_ids":[x["candidate_id"] for x in chosen],
        "pose_screen_mutation_counts":{x["candidate_id"]:x["changed_count"] for x in chosen},
        "method":"sparse seven-position full-context reconciliation",
    }
    (OUT/"result_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
