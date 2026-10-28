#!/usr/bin/env python3
"""Phase 2F: compact zero-risk TfR1 343-350 RFdiffusion refinement.

Phase 2E established that zero-risk dual-anchor backbones are feasible, but its
dual-anchor successes were too extended. The compact seed (20261236) was 55 aa
with Rg 12.0 A / CA span 43.4 A but retained only G347.

This refinement narrows the binder to 45-60 aa and samples complementary hotspot
sets inside the conserved 343-350 apical cluster. The objective is explicitly:
compact + >=2 conserved anchors + zero known cross-species-risk contacts.
"""

from __future__ import annotations
import json, math, os, shutil, time
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/"structures/11QC.pdb"
OUT=ROOT/"structures/rfdiffusion/nvidia_zero_risk_compact_refine"
URL="https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate"

DOMAIN_START,DOMAIN_END=189,383
CONTIGS="A189-383/0 45-60"

FAMILIES={
  "pair_343_347":["A343","A347"],
  "pair_344_350":["A344","A350"],
  "triplet_343_347_350":["A343","A347","A350"],
}
SEEDS={
  "pair_343_347":[20261311,20261312],
  "pair_344_350":[20261313,20261314],
  "triplet_343_347_350":[20261315,20261316],
}

PREFERRED={343,344,347,350,374}
RISK={204,208,210,211,212,348,371}

def crop(text):
    out=[];rs=set()
    for line in text.splitlines():
        if not line.startswith("ATOM  ") or len(line)<54 or line[21]!="A": continue
        try:r=int(line[22:26])
        except ValueError:continue
        if DOMAIN_START<=r<=DOMAIN_END:
            out.append(line);rs.add(r)
    needed={343,344,347,350}
    if not needed<=rs: raise RuntimeError(f"Crop missing {sorted(needed-rs)}")
    return "\n".join(out)+"\nTER\nEND\n"

def post(payload,key,max_attempts=5):
    for attempt in range(1,max_attempts+1):
        req=Request(URL,data=json.dumps(payload).encode(),headers={"Content-Type":"application/json","Authorization":f"Bearer {key}"},method="POST")
        try:
            with urlopen(req,timeout=1200) as r:return json.loads(r.read().decode())
        except HTTPError as exc:
            detail=exc.read().decode(errors="replace")
            if exc.code==429 and attempt<max_attempts:
                ra=exc.headers.get("Retry-After")
                sleep=int(ra) if ra and ra.isdigit() else 25*attempt
                print(f"429 attempt {attempt}; sleep {sleep}s",flush=True);time.sleep(sleep);continue
            transient=exc.code in {500,502,503,504} or (
                exc.code==422 and any(x in detail for x in ("InfrastructureError","AcceleratorError","CUDA error"))
            )
            if transient and attempt<max_attempts:
                sleep=25*attempt
                print(f"Transient {exc.code}; sleep {sleep}s",flush=True);time.sleep(sleep);continue
            raise RuntimeError(f"RFdiffusion HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("RFdiffusion retries exhausted")

def parse(text):
    atoms={};res={};ca={}
    for line in text.splitlines():
        if not line.startswith("ATOM  ") or len(line)<54: continue
        atom=line[12:16].strip()
        if atom.startswith("H"): continue
        ch=line[21].strip() or "_"
        try:r=int(line[22:26]);xyz=(float(line[30:38]),float(line[38:46]),float(line[46:54]))
        except ValueError:continue
        atoms.setdefault(ch,[]).append((r,atom,xyz));res.setdefault(ch,set()).add(r)
        if atom=="CA":ca.setdefault(ch,{})[r]=xyz
    return atoms,res,ca

def dist(a,b): return math.sqrt(sum((x-y)**2 for x,y in zip(a,b)))

def compact(points):
    n=len(points)
    ctr=tuple(sum(p[i] for p in points)/n for i in range(3))
    rg=math.sqrt(sum(dist(p,ctr)**2 for p in points)/n)
    span=max((dist(a,b) for i,a in enumerate(points) for b in points[i+1:]),default=0.0)
    return round(rg,3),round(span,3)

def analyze(text):
    atoms,res,ca=parse(text)
    if len(res)<2: raise RuntimeError("Expected target+binder")
    binder=min(res,key=lambda c:len(res[c]));target=max(res,key=lambda c:len(res[c]))
    tf=min(res[target]);contacts=set();clashes=0
    for _,_,bxyz in atoms[binder]:
        for rr,_,txyz in atoms[target]:
            d=dist(bxyz,txyz)
            if d<2:clashes+=1
            if d<=5:contacts.add(rr-tf+DOMAIN_START)
    pref=sorted(contacts&PREFERRED);risk=sorted(contacts&RISK)
    rg,span=compact([ca[binder][r] for r in sorted(ca[binder])])
    checks={
      "binder_45_60":45<=len(res[binder])<=60,
      "preferred_ge_2":len(pref)>=2,
      "risk_eq_0":len(risk)==0,
      "rg_le_15_5":rg<=15.5,
      "span_le_50":span<=50.0,
      "clashes_le_10":clashes<=10,
    }
    return {
      "binder_chain":binder,"target_chain":target,
      "binder_length":len(res[binder]),"target_length":len(res[target]),
      "contacts":sorted(contacts),"preferred_contacts":pref,"risk_contacts":risk,
      "radius_of_gyration_A":rg,"max_CA_span_A":span,
      "severe_clashes_lt_2A":clashes,"checks":checks,
      "scientific_pass":all(checks.values())
    }

def main():
    key=os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_API_KEY1") or os.getenv("NVIDIA_API_TOKEN") or os.getenv("NGC_API_KEY")
    if not key: raise RuntimeError("Missing NVIDIA API credential")
    target=crop(INPUT.read_text())
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"target_11QC_apical_189_383.pdb").write_text(target)

    attempts=[]
    for fam,hotspots in FAMILIES.items():
        fdir=OUT/fam;fdir.mkdir(parents=True,exist_ok=True)
        for seed in SEEDS[fam]:
            result=post({
              "input_pdb":target,"contigs":CONTIGS,"hotspot_res":hotspots,
              "diffusion_steps":50,"random_seed":seed
            },key)
            pdb=result.get("output_pdb")
            if not pdb: raise RuntimeError(f"No output_pdb for {fam}/{seed}")
            path=fdir/f"{fam}_seed_{seed}.pdb";path.write_text(pdb)
            qc=analyze(pdb)
            rec={"family":fam,"hotspots":hotspots,"seed":seed,"elapsed_ms":result.get("elapsed_ms"),
                 "pdb":str(path.relative_to(ROOT)),"qc":qc}
            attempts.append(rec);print(json.dumps(rec,indent=2),flush=True);time.sleep(15)

    passing=[x for x in attempts if x["qc"]["scientific_pass"]]
    passing.sort(key=lambda x:(
      len(x["qc"]["preferred_contacts"]),
      -x["qc"]["radius_of_gyration_A"],
      -x["qc"]["max_CA_span_A"]
    ),reverse=True)
    sel=OUT/"selected"
    if sel.exists():shutil.rmtree(sel)
    sel.mkdir(parents=True)
    for i,x in enumerate(passing[:4],1):
        dst=sel/f"selected_{i:02d}_{x['family']}_seed_{x['seed']}.pdb"
        shutil.copy2(ROOT/x["pdb"],dst);x["selected_copy"]=str(dst.relative_to(ROOT))

    summary={
      "campaign":"zero_risk_compact_refine_v1","contigs":CONTIGS,
      "families":FAMILIES,"seeds":SEEDS,
      "gate":{"preferred_ge_2":True,"risk_eq_0":True,"rg_le_15_5":True,"span_le_50":True},
      "attempt_count":len(attempts),"passing_count":len(passing),
      "selected_count":min(4,len(passing)),"selected":passing[:4],"attempts":attempts
    }
    (OUT/"selection_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps({"passing_count":len(passing),"selected":passing[:4]},indent=2),flush=True)
    if not passing: raise SystemExit("No compact zero-risk dual-anchor backbone passed")

if __name__=="__main__":main()
