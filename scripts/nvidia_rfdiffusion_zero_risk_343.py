#!/usr/bin/env python3
"""Phase 2E: zero-risk compact RFdiffusion campaign on TfR1 conserved 343-350 cluster.

Rationale:
Phase 2D seed 20261124 preserved a human apical pose after sequence design, but its
backbone still contacted human risk residue V210 and the mouse complex collapsed.
This campaign intentionally shortens the binder and concentrates generation on the
conserved 343/344/347/350 subcluster.

Scientific gate is stricter than Phase 2D:
- 50-70 aa binder
- >=2 preferred contacts
- ZERO explicit cross-species risk contacts
- Rg <=16.5 A
- CA span <=52 A
- <=10 severe clashes
"""

from __future__ import annotations

import json
import math
import os
import shutil
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/"structures/11QC.pdb"
OUT=ROOT/"structures/rfdiffusion/nvidia_zero_risk_343"
URL="https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate"

DOMAIN_START,DOMAIN_END=189,383
CONTIGS="A189-383/0 50-70"
HOTSPOTS=["A343","A344","A347","A350"]
SEEDS=[20261231,20261232,20261233,20261234,20261235,20261236]

PREFERRED={343,344,347,350,374}
RISK={204,208,210,211,212,348,371}

def crop(text):
    out=[]; residues=set()
    for line in text.splitlines():
        if not line.startswith("ATOM  ") or len(line)<54 or line[21]!="A": continue
        try:r=int(line[22:26])
        except ValueError:continue
        if DOMAIN_START<=r<=DOMAIN_END:
            out.append(line);residues.add(r)
    missing={343,344,347,350}-residues
    if missing: raise RuntimeError(f"Crop missing hotspots {sorted(missing)}")
    return "\n".join(out)+"\nTER\nEND\n"

def post(payload,key,max_attempts=5):
    for attempt in range(1,max_attempts+1):
        req=Request(URL,data=json.dumps(payload).encode(),headers={"Content-Type":"application/json","Authorization":f"Bearer {key}"},method="POST")
        try:
            with urlopen(req,timeout=1200) as r:return json.loads(r.read().decode())
        except HTTPError as exc:
            detail=exc.read().decode(errors="replace")
            if exc.code==429 and attempt<max_attempts:
                ra=exc.headers.get("Retry-After"); sleep=int(ra) if ra and ra.isdigit() else 25*attempt
                print(f"429 attempt {attempt}; sleep {sleep}s",flush=True);time.sleep(sleep);continue
            transient=exc.code in {500,502,503,504} or (exc.code==422 and any(x in detail for x in ("InfrastructureError","AcceleratorError","CUDA error")))
            if transient and attempt<max_attempts:
                sleep=25*attempt;print(f"Transient {exc.code}; sleep {sleep}s",flush=True);time.sleep(sleep);continue
            raise RuntimeError(f"RFdiffusion HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("RFdiffusion retries exhausted")

def parse(text):
    atoms={};res={};ca={}
    for line in text.splitlines():
        if not line.startswith("ATOM  ") or len(line)<54:continue
        atom=line[12:16].strip()
        if atom.startswith("H"):continue
        ch=line[21].strip() or "_"
        try:
            r=int(line[22:26]);xyz=(float(line[30:38]),float(line[38:46]),float(line[46:54]))
        except ValueError:continue
        atoms.setdefault(ch,[]).append((r,atom,xyz));res.setdefault(ch,set()).add(r)
        if atom=="CA":ca.setdefault(ch,{})[r]=xyz
    return atoms,res,ca

def dist(a,b):return math.sqrt(sum((x-y)**2 for x,y in zip(a,b)))

def geom(points):
    n=len(points); ctr=tuple(sum(p[i] for p in points)/n for i in range(3))
    rg=math.sqrt(sum(dist(p,ctr)**2 for p in points)/n)
    span=max((dist(a,b) for i,a in enumerate(points) for b in points[i+1:]),default=0.0)
    return round(rg,3),round(span,3)

def analyze(text):
    atoms,res,ca=parse(text)
    if len(res)<2:raise RuntimeError("Expected target+binder")
    binder=min(res,key=lambda c:len(res[c])); target=max(res,key=lambda c:len(res[c]))
    target_first=min(res[target]); contacts=set();clashes=0
    for _,_,bxyz in atoms[binder]:
        for out_r,_,txyz in atoms[target]:
            d=dist(bxyz,txyz)
            if d<2:clashes+=1
            if d<=5:contacts.add(out_r-target_first+DOMAIN_START)
    pref=sorted(contacts&PREFERRED);risk=sorted(contacts&RISK)
    rg,span=geom([ca[binder][r] for r in sorted(ca[binder])])
    checks={
      "binder_50_70":50<=len(res[binder])<=70,
      "preferred_ge_2":len(pref)>=2,
      "risk_eq_0":len(risk)==0,
      "rg_le_16_5":rg<=16.5,
      "span_le_52":span<=52,
      "clashes_le_10":clashes<=10,
    }
    return {
      "binder_chain":binder,"target_chain":target,"binder_length":len(res[binder]),"target_length":len(res[target]),
      "contacts":sorted(contacts),"preferred_contacts":pref,"risk_contacts":risk,
      "radius_of_gyration_A":rg,"max_CA_span_A":span,"severe_clashes_lt_2A":clashes,
      "checks":checks,"scientific_pass":all(checks.values())
    }

def main():
    key=os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_API_KEY1") or os.getenv("NVIDIA_API_TOKEN") or os.getenv("NGC_API_KEY")
    if not key:raise RuntimeError("Missing NVIDIA API credential")
    target=crop(INPUT.read_text())
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"target_11QC_apical_189_383.pdb").write_text(target)

    attempts=[]
    for seed in SEEDS:
        result=post({"input_pdb":target,"contigs":CONTIGS,"hotspot_res":HOTSPOTS,"diffusion_steps":50,"random_seed":seed},key)
        pdb=result.get("output_pdb")
        if not pdb:raise RuntimeError(f"No output_pdb seed {seed}")
        path=OUT/f"zero_risk_343_seed_{seed}.pdb";path.write_text(pdb)
        qc=analyze(pdb)
        rec={"seed":seed,"elapsed_ms":result.get("elapsed_ms"),"pdb":str(path.relative_to(ROOT)),"qc":qc}
        attempts.append(rec);print(json.dumps(rec,indent=2),flush=True);time.sleep(15)

    passing=[x for x in attempts if x["qc"]["scientific_pass"]]
    passing.sort(key=lambda x:(len(x["qc"]["preferred_contacts"]),-x["qc"]["radius_of_gyration_A"],-x["qc"]["max_CA_span_A"]),reverse=True)
    sel=OUT/"selected"
    if sel.exists():shutil.rmtree(sel)
    sel.mkdir(parents=True)
    for i,x in enumerate(passing[:3],1):
        dst=sel/f"selected_{i:02d}_seed_{x['seed']}.pdb"
        shutil.copy2(ROOT/x["pdb"],dst);x["selected_copy"]=str(dst.relative_to(ROOT))

    summary={"campaign":"zero_risk_343_v1","contigs":CONTIGS,"hotspots":HOTSPOTS,"seeds":SEEDS,
             "attempt_count":len(attempts),"passing_count":len(passing),"selected_count":min(3,len(passing)),
             "selected":passing[:3],"attempts":attempts}
    (OUT/"selection_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps({"passing_count":len(passing),"selected":passing[:3]},indent=2),flush=True)
    if not passing:raise SystemExit("No zero-risk compact backbone passed")

if __name__=="__main__":main()
