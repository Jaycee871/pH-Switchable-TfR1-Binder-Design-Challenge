#!/usr/bin/env python3
"""Phase 3H: cross-species-tolerant ProteinMPNN rescue on seed 20261124.

The Phase 4F parent retains the human apical surface but mouse collapses.
Backbone geometry shows binder positions 74-76 contact both conserved human L209
and species-risk V210 (mouse D213). This redesign biases that narrow strip toward
neutral, chemically tolerant residues rather than charged or bulky hydrophobic ones.

This is a pose-rescue experiment, not pH engineering.
"""

from __future__ import annotations
import csv,json,os,re,time
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request,urlopen

ROOT=Path(__file__).resolve().parents[1]
BACKBONE=ROOT/"structures/rfdiffusion/nvidia_apical_only/selected/selected_01_conserved_343_seed_20261124.pdb"
OUT=ROOT/"designs/family_B_multistate/nvidia_apical_rescue"
URL="https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict"
ALPHABET="ACDEFGHIKLMNPQRSTVWYX"
LENGTH=86

# Binder positions from direct geometry on the RFdiffusion complex.
PRIMARY_TOLERANCE=[74,75,76]  # bridge L209 and V210/D213
SUPPORT=[73,77]               # same local strip, slightly more peripheral
G347_CONTACT=[71]             # conserved anchor

def post(payload,key,max_attempts=4):
    for attempt in range(1,max_attempts+1):
        req=Request(URL,data=json.dumps(payload).encode(),headers={"Content-Type":"application/json","Authorization":f"Bearer {key}"},method="POST")
        try:
            with urlopen(req,timeout=1200) as r:return json.loads(r.read().decode())
        except HTTPError as exc:
            detail=exc.read().decode(errors="replace")
            if exc.code==429 and attempt<max_attempts:
                ra=exc.headers.get("Retry-After");sleep=int(ra) if ra and ra.isdigit() else 20*attempt
                print(f"429 attempt {attempt}; sleep {sleep}s",flush=True);time.sleep(sleep);continue
            raise RuntimeError(f"ProteinMPNN HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("ProteinMPNN retries exhausted")

def parse_fasta(text):
    out=[];h=None;s=[]
    for raw in text.splitlines():
        line=raw.strip()
        if not line:continue
        if line.startswith(">"):
            if h is not None:out.append((h,"".join(s)))
            h=line[1:];s=[]
        else:s.append(line)
    if h is not None:out.append((h,"".join(s)))
    return out

def blank():return [0.0]*len(ALPHABET)
def add(row,aav,value):
    for aa in aav:row[ALPHABET.index(aa)]+=value

def matrix():
    m=[blank() for _ in range(LENGTH)]
    tolerant="AGQNST"
    disfavor="DEKRFWYILVM"
    for p in PRIMARY_TOLERANCE:
        add(m[p-1],tolerant,1.6)
        add(m[p-1],disfavor,-0.8)
    for p in SUPPORT:
        add(m[p-1],tolerant,0.9)
        add(m[p-1],disfavor,-0.35)
    for p in G347_CONTACT:
        add(m[p-1],"AGSTNQ",1.0)
        add(m[p-1],"FWY",-0.5)
    return m

def hv(header,pat):
    x=re.search(pat,header);return x.group(1).strip() if x else ""

def longest(seq):
    best=run=1
    for i in range(1,len(seq)):
        run=run+1 if seq[i]==seq[i-1] else 1;best=max(best,run)
    return best

def main():
    key=os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_API_KEY1") or os.getenv("NVIDIA_API_TOKEN") or os.getenv("NGC_API_KEY")
    if not key:raise RuntimeError("Missing NVIDIA API credential")
    if not BACKBONE.exists():raise RuntimeError("Missing seed 20261124 backbone")
    bias=matrix()
    payload={
      "input_pdb":BACKBONE.read_text(),
      "input_pdb_chains":["A"],
      "num_seq_per_target":8,
      "sampling_temp":[0.10,0.20],
      "use_soluble_model":False,
      "ca_only":False,
      "omit_AAs":["C"],
      "bias_by_res_jsonl":json.dumps({"input":{"A":bias}}),
      "random_seed":20261301,
    }
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"bias_spec.json").write_text(json.dumps({
      "alphabet":ALPHABET,"primary_tolerance":PRIMARY_TOLERANCE,"support":SUPPORT,
      "g347_contact":G347_CONTACT,"preferred_residues":"AGQNST",
      "rationale":"neutral tolerance across conserved L209 and human V210 / mouse D213"
    },indent=2)+"\n")
    result=post(payload,key)
    mfasta=result.get("mfasta","")
    if not mfasta.strip():raise RuntimeError("No mfasta returned")
    (OUT/"designed_sequences.fasta").write_text(mfasta)

    records=[(h,s) for h,s in parse_fasta(mfasta) if not h.lower().startswith("input") and "native" not in h.lower() and not h.lower().startswith("wt")]
    rows=[]
    for i,(header,seq) in enumerate(records,1):
        binder=[x for x in seq.split("/") if x][0]
        if len(binder)!=LENGTH:raise RuntimeError(f"Expected 86 aa, got {len(binder)}")
        acidic=binder.count("D")+binder.count("E");basic=binder.count("K")+binder.count("R")
        cf=(acidic+basic)/len(binder)
        flags=[]
        if cf>0.65:flags.append("charge_fraction_gt_0.65")
        if max(binder.count(x) for x in set(binder))/len(binder)>0.40:flags.append("single_aa_fraction_gt_0.40")
        if longest(binder)>4:flags.append("homopolymer_run_gt_4")
        strip="".join(binder[p-1] for p in [71,73,74,75,76,77])
        tolerant_count=sum(binder[p-1] in set("AGQNST") for p in PRIMARY_TOLERANCE)
        rows.append({
          "candidate_id":f"TFR1_XS343_{i:02d}","backbone_id":"NB01_conserved_343_s20261124",
          "source_backbone":str(BACKBONE.relative_to(ROOT)),"binder_chain":"A","target_chain":"B",
          "binder_length":len(binder),"sequence":binder,
          "temperature":hv(header,r"(?:^|,\s*)T=([^,]+)"),"sample":hv(header,r"(?:^|,\s*)sample=([^,]+)"),
          "mpnn_score":hv(header,r"(?:^|,\s*)score=([^,]+)"),"global_score":hv(header,r"(?:^|,\s*)global_score=([^,]+)"),
          "histidine_count":binder.count("H"),"acidic_count":acidic,"basic_count":basic,"charge_fraction":round(cf,4),
          "interface_strip_71_73_77":strip,"primary_tolerant_count_74_76":tolerant_count,
          "prefilter_flags":";".join(flags),"prefilter_pass":not flags and tolerant_count>=2,"header":header
        })

    fields=rows[0].keys()
    with (OUT/"all_candidates.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields);w.writeheader();w.writerows(rows)
    passing=[r for r in rows if r["prefilter_pass"]]
    passing.sort(key=lambda r:(-r["primary_tolerant_count_74_76"],float(r["mpnn_score"]) if r["mpnn_score"] else 999))
    chosen=passing[:2]
    with (OUT/"pose_screen_candidates.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields);w.writeheader();w.writerows(chosen)
    (OUT/"result_summary.json").write_text(json.dumps({
      "candidate_count":len(rows),"prefilter_pass_count":len(passing),
      "pose_screen_ids":[x["candidate_id"] for x in chosen],
      "primary_tolerance_positions":PRIMARY_TOLERANCE
    },indent=2)+"\n")
    print((OUT/"result_summary.json").read_text(),flush=True)

if __name__=="__main__":main()
