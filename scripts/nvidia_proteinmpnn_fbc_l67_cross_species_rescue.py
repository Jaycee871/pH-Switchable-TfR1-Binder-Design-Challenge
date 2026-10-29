#!/usr/bin/env python3
"""Fallback cross-species sequence rescue for FreeBindCraft l67 trajectory.

Do NOT run unless the v4 hard-hotspot generator test fails to yield a superior
cross-species parent.

The l67 trajectory has a real human apical contact at E350/G351, but mouse
full-complex prediction collapses. Interface mapping identifies binder W7/Y8/Y11
as the main residues depending on non-conserved human peripheral chemistry.
R10 is the cleanest conserved E350/G351 anchor and must be retained.

This script redesigns binder chain B only with a narrow ProteinMPNN bias.
It is pose rescue, not pH engineering.
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

ROOT=Path(__file__).resolve().parents[1]
BACKBONE=ROOT/"freebindcraft/kaggle_smoke_v1_output/tfr1_freebindcraft_smoke_v1/designs/Trajectory/Relaxed/TFR1_FBC343_l67_s385919.pdb"
OUT=ROOT/"designs/family_B_multistate/fbc_l67_cross_species_rescue"
URL="https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict"
ALPHABET="ACDEFGHIKLMNPQRSTVWYX"
LENGTH=67

PRIMARY=[7,8,11]      # W7/Y8/Y11: strongest variable-periphery liabilities
SECONDARY=[14,18]     # H14/K18: contact human E369 -> mouse Q371
ANCHOR=10             # R10 -> conserved human E350/G351


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


def strip_hydrogens(text):
    out=[]
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom=line[12:16].strip()
        elem=line[76:78].strip() if len(line)>=78 else ""
        if atom.startswith("H") or elem=="H":
            continue
        out.append(line)
    return "\n".join(out)+"\nTER\nEND\n"


def blank():
    return [0.0]*len(ALPHABET)


def add(row,aas,value):
    for aa in aas:
        row[ALPHABET.index(aa)] += value


def matrix():
    m=[blank() for _ in range(LENGTH)]

    # Broadly tolerant chemistry for variable human/mouse peripheral contacts.
    for p in PRIMARY:
        add(m[p-1],"AGQNST",1.7)
        add(m[p-1],"FWY",-1.2)
        add(m[p-1],"DEKRILMV",-0.35)

    # Secondary electrostatic liabilities: softer redesign pressure.
    for p in SECONDARY:
        add(m[p-1],"AGQNST",0.8)
        add(m[p-1],"DEKR",-0.45)

    # Preserve the experimentally useful conserved anchor geometry.
    add(m[ANCHOR-1],"R",3.0)
    add(m[ANCHOR-1],"K",0.5)
    return m


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


def binder_segment(sequence):
    parts=[x for x in sequence.split("/") if x]
    matches=[x for x in parts if len(x)==LENGTH]
    if len(matches)==1:
        return matches[0]
    if len(parts)==1 and len(parts[0])==LENGTH:
        return parts[0]
    raise RuntimeError(f"Could not identify 67-aa binder segment: {[len(x) for x in parts]}")


def main():
    key=(os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_API_KEY1")
         or os.getenv("NVIDIA_API_TOKEN") or os.getenv("NGC_API_KEY"))
    if not key:
        raise RuntimeError("Missing NVIDIA API credential")
    if not BACKBONE.exists():
        raise RuntimeError(f"Missing l67 backbone: {BACKBONE}")

    bias=matrix()
    payload={
        "input_pdb":strip_hydrogens(BACKBONE.read_text(errors="replace")),
        "input_pdb_chains":["B"],
        "num_seq_per_target":8,
        "sampling_temp":[0.10,0.20],
        "use_soluble_model":False,
        "ca_only":False,
        "omit_AAs":["C"],
        "bias_by_res_jsonl":json.dumps({"input":{"B":bias}}),
        "random_seed":20261031,
    }

    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"bias_spec.json").write_text(json.dumps({
        "primary_liability_positions":PRIMARY,
        "secondary_liability_positions":SECONDARY,
        "conserved_anchor_position":ANCHOR,
        "required_anchor_residue":"R",
        "preferred_tolerant_residues":"AGQNST",
        "rationale":"Reduce dependence on non-conserved human peripheral chemistry while preserving R10->E350/G351."
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

    rows=[]
    for i,(header,sequence) in enumerate(records,1):
        binder=binder_segment(sequence)
        acidic=binder.count("D")+binder.count("E")
        basic=binder.count("K")+binder.count("R")
        charge_fraction=(acidic+basic)/len(binder)
        primary_tolerant=sum(binder[p-1] in set("AGQNST") for p in PRIMARY)
        anchor_retained=binder[ANCHOR-1]=="R"
        flags=[]
        if charge_fraction>0.65:
            flags.append("charge_fraction_gt_0.65")
        if max(binder.count(x) for x in set(binder))/len(binder)>0.40:
            flags.append("single_aa_fraction_gt_0.40")
        if longest_run(binder)>4:
            flags.append("homopolymer_run_gt_4")
        if not anchor_retained:
            flags.append("R10_anchor_not_retained")
        if primary_tolerant<2:
            flags.append("fewer_than_2_primary_tolerant_positions")

        rows.append({
            "candidate_id":f"TFR1_FBC_L67_XS_{i:02d}",
            "sequence":binder,
            "mpnn_score":hv(header,r"(?:^|,\s*)score=([^,]+)"),
            "global_score":hv(header,r"(?:^|,\s*)global_score=([^,]+)"),
            "temperature":hv(header,r"(?:^|,\s*)T=([^,]+)"),
            "R10_anchor_retained":anchor_retained,
            "primary_7_8_11":"".join(binder[p-1] for p in PRIMARY),
            "secondary_14_18":"".join(binder[p-1] for p in SECONDARY),
            "primary_tolerant_count":primary_tolerant,
            "charge_fraction":round(charge_fraction,4),
            "histidine_count":binder.count("H"),
            "prefilter_flags":";".join(flags),
            "prefilter_pass":not flags,
            "header":header,
        })

    fields=rows[0].keys()
    with (OUT/"all_candidates.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields);w.writeheader();w.writerows(rows)

    passing=[x for x in rows if x["prefilter_pass"]]
    passing.sort(key=lambda x:(
        -x["primary_tolerant_count"],
        float(x["mpnn_score"]) if x["mpnn_score"] else 999
    ))
    chosen=passing[:2]
    with (OUT/"pose_screen_candidates.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields);w.writeheader();w.writerows(chosen)

    (OUT/"result_summary.json").write_text(json.dumps({
        "candidate_count":len(rows),
        "prefilter_pass_count":len(passing),
        "pose_screen_ids":[x["candidate_id"] for x in chosen],
        "note":"Fallback only; full human/mouse pose validation required before pH engineering."
    },indent=2)+"\n")
    print((OUT/"result_summary.json").read_text())


if __name__=="__main__":
    main()
