#!/usr/bin/env python3
"""IL4: encode the Phase 4Q human apical pose into sequence with ProteinMPNN.

The Phase 4Q v6 human complex is structurally plausible and template-guided
(Boltz-2 protein iPTM ~0.690, preferred contacts 231/347/350), but the exact same
v6 sequence fails to recover that pose when predicted sequence-only
(full-human protein iPTM ~0.318, distal interface).

Therefore IL4 treats the Phase 4Q human complex as a structural design state,
not as independent binding evidence. ProteinMPNN redesigns only the actual
18-residue binder interface observed in that state while the other 42 binder
positions remain fixed.

No new His or Cys is allowed. The next gate is template-free full-human
Boltz-2; only sequences that recover the human apical pose without a structural
template are allowed to proceed to mouse testing.
"""

from __future__ import annotations

import csv
import json
import os
import re
import time
from pathlib import Path
from http.client import RemoteDisconnected
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from Bio.PDB import MMCIFParser, PDBIO

ROOT=Path(__file__).resolve().parents[1]
HUMAN_CIF=ROOT/(
    "predictions/freebindcraft_v6_gate_override_pose_screen/"
    "TFR1_FBC343_V6_l60_s289817/human/structure_1.cif"
)
OUT=ROOT/"designs/inference_loop/human_anchor_il4"
URL="https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict"

PARENT="SPSVISTTTTTYTYHGKQVTRIRTTINPGEKIEPAPGERQMRGRDWETGYIIVHTYKVNS"
LENGTH=60
DESIGN_POSITIONS={8,12,14,15,19,21,23,39,41,42,43,44,46,47,52,54,56,58}
FIXED_POSITIONS=[i for i in range(1,LENGTH+1) if i not in DESIGN_POSITIONS]
ALPHABET="ACDEFGHIKLMNPQRSTVWYX"


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
        except (RemoteDisconnected, URLError, TimeoutError, ConnectionResetError) as exc:
            if attempt<max_attempts:
                delay=20*attempt
                print(f"ProteinMPNN transport retry {attempt}: {exc}; sleep {delay}s",flush=True)
                time.sleep(delay)
                continue
            raise RuntimeError(f"ProteinMPNN transport failure: {exc}") from exc
    raise RuntimeError("ProteinMPNN retries exhausted")


def cif_to_pdb_text(cif:Path)->str:
    OUT.mkdir(parents=True,exist_ok=True)
    st=MMCIFParser(QUIET=True).get_structure(cif.stem,str(cif))
    tmp=OUT/"human_anchor_full_complex.pdb"
    io=PDBIO();io.set_structure(st);io.save(str(tmp))
    lines=[]
    for line in tmp.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        atom=line[12:16].strip()
        elem=line[76:78].strip() if len(line)>=78 else ""
        if atom.startswith("H") or elem=="H":
            continue
        lines.append(line)
    return "\n".join(lines)+"\nTER\nEND\n"


def parse_fasta(text):
    out=[];h=None;s=[]
    for raw in text.splitlines():
        line=raw.strip()
        if not line: continue
        if line.startswith(">"):
            if h is not None: out.append((h,"".join(s)))
            h=line[1:];s=[]
        else:
            s.append(line)
    if h is not None: out.append((h,"".join(s)))
    return out


def header_value(header,pat):
    m=re.search(pat,header)
    return m.group(1).strip() if m else ""


def binder_segment(seq):
    parts=[x for x in seq.split("/") if x]
    matches=[x for x in parts if len(x)==LENGTH]
    if len(matches)==1:
        return matches[0]
    if len(parts)==1 and len(parts[0])==LENGTH:
        return parts[0]
    raise RuntimeError(f"Could not identify 60-aa binder in {[len(x) for x in parts]}")


def longest_run(seq):
    if not seq: return 0
    best=run=1
    for i in range(1,len(seq)):
        run=run+1 if seq[i]==seq[i-1] else 1
        best=max(best,run)
    return best


def bias_matrix():
    # Mild parent regularization only. The purpose is to let ProteinMPNN encode
    # the observed interface, not to force recovery of the old v6 chemistry.
    mat=[[0.0]*len(ALPHABET) for _ in range(LENGTH)]
    for pos in DESIGN_POSITIONS:
        aa=PARENT[pos-1]
        mat[pos-1][ALPHABET.index(aa)]+=0.25
    return mat


def main():
    key=(os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_API_KEY1")
         or os.getenv("NVIDIA_API_TOKEN") or os.getenv("NGC_API_KEY"))
    if not key:
        raise RuntimeError("Missing NVIDIA API credential")
    if not HUMAN_CIF.exists():
        raise FileNotFoundError(HUMAN_CIF)

    pdb_text=cif_to_pdb_text(HUMAN_CIF)

    fixed=json.dumps({"input":{"B":FIXED_POSITIONS}})
    payload={
        "input_pdb":pdb_text,
        "input_pdb_chains":["B"],
        "num_seq_per_target":40,
        "sampling_temp":[0.08,0.12,0.18,0.25],
        "use_soluble_model":False,
        "ca_only":False,
        "omit_AAs":["C","H"],
        "fixed_positions_jsonl":fixed,
        "bias_by_res_jsonl":json.dumps({"input":{"B":bias_matrix()}}),
        "random_seed":20261110,
    }

    spec={
        "method":"human-anchor interface encoding from Phase 4Q full-human pose",
        "structural_state":str(HUMAN_CIF.relative_to(ROOT)),
        "parent_sequence":PARENT,
        "design_positions":sorted(DESIGN_POSITIONS),
        "fixed_positions":FIXED_POSITIONS,
        "template_guided_state_metrics":{
            "human_protein_iptm":0.6898922324180603,
            "preferred_contacts":[231,347,350],
        },
        "sequence_only_parent_control":{
            "human_protein_iptm":0.3181271255016327,
            "preferred_contacts":[],
            "apical_fraction":0.0,
        },
        "no_new_histidine":True,
        "no_new_cysteine":True,
        "proteinmpnn_batches":[{"seed":20261110,"n":40},{"seed":20261111,"n":40}],
        "next_gate":"template-free full-human Boltz-2",
    }
    (OUT/"design_spec.json").write_text(json.dumps(spec,indent=2)+"\n")

    mfasta_parts=[]
    for batch_index,seed in enumerate((20261110,20261111),start=1):
        payload["random_seed"]=seed
        print(f"ProteinMPNN IL4 batch {batch_index}/2, seed={seed}, n=40",flush=True)
        result=post(payload,key)
        part=result.get("mfasta","")
        if not part.strip():
            raise RuntimeError(f"ProteinMPNN batch {batch_index} returned no mfasta")
        (OUT/f"designed_sequences_batch_{batch_index}.fasta").write_text(part)
        mfasta_parts.append(part.rstrip())

    mfasta="\n".join(mfasta_parts)+"\n"
    (OUT/"designed_sequences.fasta").write_text(mfasta)

    records=[
        (h,s) for h,s in parse_fasta(mfasta)
        if not h.lower().startswith("input")
        and "native" not in h.lower()
        and not h.lower().startswith("wt")
    ]

    rows=[];seen=set();parent_run=longest_run(PARENT)
    for i,(header,rawseq) in enumerate(records,1):
        seq=binder_segment(rawseq)
        if seq in seen: continue
        seen.add(seq)
        changed=[p for p in range(1,LENGTH+1) if seq[p-1]!=PARENT[p-1]]
        outside=[p for p in changed if p not in DESIGN_POSITIONS]
        flags=[]
        if outside: flags.append("mutation_outside_human_interface")
        if any(seq[p-1]=="H" and PARENT[p-1]!="H" for p in changed):
            flags.append("new_histidine_before_pH_stage")
        if any(seq[p-1]=="C" and PARENT[p-1]!="C" for p in changed):
            flags.append("new_cysteine")
        if longest_run(seq)>parent_run:
            flags.append("homopolymer_worse_than_parent")
        charge=(seq.count("D")+seq.count("E")+seq.count("K")+seq.count("R"))/LENGTH
        if charge>0.50:
            flags.append("charge_fraction_gt_0.50")

        rows.append({
            "candidate_id":f"TFR1_ILOOP_IL4_{i:03d}",
            "sequence":seq,
            "mpnn_score":header_value(header,r"(?:^|,\s*)score=([^,]+)"),
            "global_score":header_value(header,r"(?:^|,\s*)global_score=([^,]+)"),
            "temperature":header_value(header,r"(?:^|,\s*)T=([^,]+)"),
            "changed_positions":";".join(map(str,changed)),
            "changed_count":len(changed),
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
    def score(x):
        mp=float(x["mpnn_score"]) if x["mpnn_score"] else 999.0
        # very mild mutation-count regularization to avoid gratuitous interface rewrite
        return mp+0.015*int(x["changed_count"])
    passing.sort(key=score)

    chosen=[]
    for row in passing:
        seq=row["sequence"]
        if chosen:
            mind=min(sum(a!=b for a,b in zip(seq,x["sequence"])) for x in chosen)
            if mind<3:
                continue
        chosen.append(row)
        if len(chosen)>=12:
            break
    if len(chosen)<12:
        for row in passing:
            if row not in chosen:
                chosen.append(row)
            if len(chosen)>=12:
                break

    with (OUT/"pose_screen_candidates.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields);w.writeheader();w.writerows(chosen)

    summary={
        "candidate_count":len(rows),
        "unique_candidate_count":len(seen),
        "prefilter_pass_count":len(passing),
        "human_screen_count":len(chosen),
        "human_screen_ids":[x["candidate_id"] for x in chosen],
        "mutation_counts":{x["candidate_id"]:x["changed_count"] for x in chosen},
    }
    (OUT/"result_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
