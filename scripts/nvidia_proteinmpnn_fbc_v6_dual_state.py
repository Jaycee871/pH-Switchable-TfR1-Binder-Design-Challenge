#!/usr/bin/env python3
"""Phase 3K: true human+mouse dual-state ProteinMPNN design on the v6 backbone.

Builds two structural states in one PDB:
  A+B = human TfR1 apical state + binder
  C+D = mouse homologous apical state + binder copy

The mouse state is translated far away so the two complexes do not interact.
ProteinMPNN designs chains B and D simultaneously, with every B_i tied to D_i,
forcing one identical binder sequence to be compatible with both species.

Only interface positions are designable. Non-interface binder positions are fixed
to the v6 parent to preserve the relaxed fold.

This is cross-species sequence design, not pH engineering.
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

ROOT = Path(__file__).resolve().parents[1]
HUMAN = ROOT / (
    "freebindcraft/kaggle_v6_tempered_binary_results/designs/Trajectory/Relaxed/"
    "TFR1_FBC343_V6_l60_s289817.pdb"
)
MOUSE = ROOT / (
    "predictions/freebindcraft_v6_two_position_cleanup_pose_screen/"
    "TFR1_FBC_V6_XS_08/TFR1_FBC_V6_XS_08_mouse_apical_template.pdb"
)
OUT = ROOT / "designs/family_B_multistate/fbc_v6_dual_state_mpnn"
URL = "https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict"

LENGTH = 60
TRANSLATE_X = 140.0

# Human v6 independent 5 A interface positions. We add a small shell around
# the hotspot-facing segment to let the joint objective adapt chemistry while
# preserving the parent fold.
HUMAN_INTERFACE = {
    12,14,15,17,19,21,23,39,40,41,42,43,44,45,46,47,52,54,56,58
}
DESIGN_POSITIONS = sorted(HUMAN_INTERFACE)
FIXED_POSITIONS = [p for p in range(1, LENGTH+1) if p not in HUMAN_INTERFACE]

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}


def post(payload,key,max_attempts=4):
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
                sleep=int(ra) if ra and ra.isdigit() else 20*attempt
                print(f"429 attempt {attempt}; sleep {sleep}s",flush=True)
                time.sleep(sleep)
                continue
            raise RuntimeError(f"ProteinMPNN HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("ProteinMPNN retries exhausted")


def strip_hydrogen_atom_lines(text):
    out=[]
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom=line[12:16].strip()
        elem=line[76:78].strip() if len(line)>=78 else ""
        if atom.startswith("H") or elem=="H":
            continue
        out.append(line)
    return out


def rewrite_chain_xyz(line, chain, dx=0.0):
    chars=list(line)
    if len(chars)<80:
        chars.extend(" "*(80-len(chars)))
    chars[21]=chain
    if dx:
        x=float(line[30:38])+dx
        chars[30:38]=list(f"{x:8.3f}")
    return "".join(chars).rstrip()


def build_dual_state_pdb():
    human_lines=strip_hydrogen_atom_lines(HUMAN.read_text(errors="replace"))
    mouse_lines=strip_hydrogen_atom_lines(MOUSE.read_text(errors="replace"))

    out=[]
    for line in human_lines:
        ch=line[21].strip()
        if ch=="A":
            out.append(rewrite_chain_xyz(line,"A"))
        elif ch=="B":
            out.append(rewrite_chain_xyz(line,"B"))

    out.append("TER")

    for line in mouse_lines:
        ch=line[21].strip()
        if ch=="A":
            out.append(rewrite_chain_xyz(line,"C",TRANSLATE_X))
        elif ch=="B":
            out.append(rewrite_chain_xyz(line,"D",TRANSLATE_X))

    out.extend(["TER","END"])
    return "\n".join(out)+"\n"


def chain_sequence(pdb_text,chain):
    seen={}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM") or line[21].strip()!=chain:
            continue
        r=int(line[22:26])
        if r not in seen:
            seen[r]=AA3_TO_1[line[17:20].strip()]
    return "".join(seen[k] for k in sorted(seen))


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


def hv(header,pat):
    m=re.search(pat,header)
    return m.group(1).strip() if m else ""


def longest_run(seq):
    best=run=1
    for i in range(1,len(seq)):
        run=run+1 if seq[i]==seq[i-1] else 1
        best=max(best,run)
    return best


def designed_binder_from_output(sequence):
    parts=[x for x in sequence.split("/") if x]
    sixty=[x for x in parts if len(x)==LENGTH]
    if len(sixty)>=2:
        # B and D are tied and must be identical.
        if len(set(sixty))!=1:
            raise RuntimeError(f"Tied B/D chains differ: {sixty}")
        return sixty[0]
    if len(sixty)==1:
        return sixty[0]
    if len(parts)==1 and len(parts[0])==LENGTH:
        return parts[0]
    raise RuntimeError(f"Cannot identify tied 60-aa binder: {[len(x) for x in parts]}")


def main():
    key=(os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_API_KEY1")
         or os.getenv("NVIDIA_API_TOKEN") or os.getenv("NGC_API_KEY"))
    if not key:
        raise RuntimeError("Missing NVIDIA API credential")
    if not HUMAN.exists() or not MOUSE.exists():
        raise RuntimeError(f"Missing state template: human={HUMAN.exists()} mouse={MOUSE.exists()}")

    dual=build_dual_state_pdb()
    parent=chain_sequence(dual,"B")
    mouse_binder=chain_sequence(dual,"D")
    if len(parent)!=LENGTH or len(mouse_binder)!=LENGTH:
        raise RuntimeError(f"Unexpected binder lengths: B={len(parent)} D={len(mouse_binder)}")

    # Force both copies to start from the same parent sequence before design.
    # The D sequence in the structural template came from the SN cleanup variant;
    # only coordinates matter. Rewrite D residue names to parent identity so the
    # two tied states have the same starting sequence.
    parent3 = {
        "A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","E":"GLU","Q":"GLN",
        "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
        "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL",
    }
    d_residues=[]
    seen=set()
    for line in dual.splitlines():
        if line.startswith("ATOM") and line[21].strip()=="D":
            r=int(line[22:26])
            if r not in seen:
                d_residues.append(r);seen.add(r)
    d_index={r:i+1 for i,r in enumerate(d_residues)}
    rewritten=[]
    for line in dual.splitlines():
        if line.startswith("ATOM") and line[21].strip()=="D":
            r=int(line[22:26]); pos=d_index[r]
            chars=list(line)
            chars[17:20]=list(f"{parent3[parent[pos-1]]:>3}")
            line="".join(chars)
        rewritten.append(line)
    dual="\n".join(rewritten)+"\n"

    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"dual_state_input.pdb").write_text(dual)

    fixed=json.dumps({
        "input":{
            "B":FIXED_POSITIONS,
            "D":FIXED_POSITIONS,
        }
    })
    tied=json.dumps({
        "input":[
            {"B":[i],"D":[i]} for i in range(1,LENGTH+1)
        ]
    })

    payload={
        "input_pdb":dual,
        "input_pdb_chains":["B","D"],
        "num_seq_per_target":16,
        "sampling_temp":[0.10,0.20,0.30],
        "use_soluble_model":False,
        "ca_only":False,
        "omit_AAs":["C"],
        "fixed_positions_jsonl":fixed,
        "tied_positions_jsonl":tied,
        "random_seed":20261102,
    }

    (OUT/"design_spec.json").write_text(json.dumps({
        "human_state":"A+B",
        "mouse_state":"C+D",
        "translation_x_A":TRANSLATE_X,
        "binder_length":LENGTH,
        "design_positions":DESIGN_POSITIONS,
        "fixed_positions":FIXED_POSITIONS,
        "tied_positions":"B_i == D_i for i=1..60",
        "parent_sequence":parent,
        "note":"Two structural states jointly constrain one identical binder sequence."
    },indent=2)+"\n")

    result=post(payload,key)
    mfasta=result.get("mfasta","")
    if not mfasta.strip():
        raise RuntimeError("No mfasta returned")
    (OUT/"designed_sequences.fasta").write_text(mfasta)

    records=[
        (h,s) for h,s in parse_fasta(mfasta)
        if not h.lower().startswith("input")
        and "native" not in h.lower()
        and not h.lower().startswith("wt")
    ]

    rows=[];seen_seq=set()
    for i,(header,sequence) in enumerate(records,1):
        binder=designed_binder_from_output(sequence)
        if binder in seen_seq:
            continue
        seen_seq.add(binder)

        changed=[p for p in range(1,LENGTH+1) if binder[p-1]!=parent[p-1]]
        off_interface=[p for p in changed if p not in HUMAN_INTERFACE]

        acidic=binder.count("D")+binder.count("E")
        basic=binder.count("K")+binder.count("R")
        charge_fraction=(acidic+basic)/LENGTH
        flags=[]
        if off_interface:
            flags.append("changed_fixed_noninterface_position")
        if charge_fraction>0.55:
            flags.append("charge_fraction_gt_0.55")
        if longest_run(binder)>4:
            flags.append("homopolymer_run_gt_4")
        if binder.count("C"):
            flags.append("contains_cysteine")

        rows.append({
            "candidate_id":f"TFR1_FBC_V6_MS_{i:02d}",
            "sequence":binder,
            "mpnn_score":hv(header,r"(?:^|,\s*)score=([^,]+)"),
            "global_score":hv(header,r"(?:^|,\s*)global_score=([^,]+)"),
            "temperature":hv(header,r"(?:^|,\s*)T=([^,]+)"),
            "changed_positions":";".join(map(str,changed)),
            "changed_count":len(changed),
            "charge_fraction":round(charge_fraction,4),
            "histidine_count":binder.count("H"),
            "prefilter_flags":";".join(flags),
            "prefilter_pass":not flags,
            "header":header,
        })

    if not rows:
        raise RuntimeError("No dual-state sequences parsed")

    fields=list(rows[0].keys())
    with (OUT/"all_candidates.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields);w.writeheader();w.writerows(rows)

    passing=[x for x in rows if x["prefilter_pass"]]
    passing.sort(key=lambda x:(
        float(x["mpnn_score"]) if x["mpnn_score"] else 999,
        x["changed_count"],
    ))
    chosen=passing[:4]

    with (OUT/"pose_screen_candidates.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields);w.writeheader();w.writerows(chosen)

    (OUT/"result_summary.json").write_text(json.dumps({
        "candidate_count":len(rows),
        "unique_candidate_count":len(seen_seq),
        "prefilter_pass_count":len(passing),
        "pose_screen_ids":[x["candidate_id"] for x in chosen],
        "note":"True human+mouse tied-sequence multistate ProteinMPNN. Full independent Boltz-2 validation required."
    },indent=2)+"\n")
    print((OUT/"result_summary.json").read_text())


if __name__=="__main__":
    main()
