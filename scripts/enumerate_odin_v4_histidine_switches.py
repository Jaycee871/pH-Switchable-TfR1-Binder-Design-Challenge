#!/usr/bin/env python3
"""Phase 5A: enumerate mechanistically motivated histidine-switch variants.

This stage is HARD-GATED on an independently validated Odin full-context v4
cross-species parent. It does not run if Phase 4U produced no pose-pass ID.

The goal is candidate generation, not a pH claim. Histidine substitutions are
ranked from the actual human full-complex interface:

- contacts to Lys/Arg target residues are favoured because His protonation near
  pH 6 can introduce positive-positive electrostatic repulsion;
- contacts to Asp/Glu are penalised because protonated His can instead be
  stabilised by a negative target residue;
- known human/mouse species-risk positions are penalised;
- Gly/Pro/Cys/His parent residues are not mutated;
- single, constrained double, and constrained triple-His hypotheses are emitted.

Explicit neutral-pH revalidation and pKa/protonation analysis remain downstream.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from Bio.PDB import MMCIFParser

ROOT=Path(__file__).resolve().parents[1]
POSE_ROOT=ROOT/"predictions/odin_full_context_v4_pose_screen"
PASS_IDS=POSE_ROOT/"pose_pass_ids.json"
BATCH=POSE_ROOT/"batch_summary.json"
HUMAN_CIF=POSE_ROOT/"human/structure_1.cif"
OUT=ROOT/"designs/pH_switch/odin_full_context_v4_histidine_scan"

HUMAN_START=121
PREFERRED_H={206,207,209,214,231,343,344,347,350,374}
SPECIES_RISK={204,208,210,211,212,348,371}
BASIC_RESN={"LYS","ARG"}
ACIDIC_RESN={"ASP","GLU"}
EXCLUDED_SOURCE_AA={"G","P","C","H"}

AA3_TO_1={
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}


def distance(a,b):
    return math.sqrt(sum((float(x)-float(y))**2 for x,y in zip(a.coord,b.coord)))


def load_model(path:Path):
    return next(MMCIFParser(QUIET=True).get_structure(path.stem,str(path)).get_models())


def binder_sequence(chain):
    seq=[]
    positions=[]
    for r in chain:
        if r.id[0]!=" ":
            continue
        aa=AA3_TO_1.get(r.resname)
        if aa:
            seq.append(aa)
            positions.append(r.id[1])
    return "".join(seq),positions


def mutate(seq,positions):
    x=list(seq)
    for p in positions:
        x[p-1]="H"
    return "".join(x)


def main():
    if not PASS_IDS.exists() or not BATCH.exists() or not HUMAN_CIF.exists():
        raise SystemExit("Phase 4U pose-pass files are not present")

    pass_ids=json.loads(PASS_IDS.read_text()).get("pose_pass_ids",[])
    if not pass_ids:
        print(json.dumps({
            "status":"skipped",
            "reason":"Phase 4U produced no cross-species pose-pass parent",
        },indent=2))
        return

    batch=json.loads(BATCH.read_text())
    candidate_id=batch.get("candidate_id")
    parent_seq=batch.get("binder_sequence")
    if candidate_id not in pass_ids:
        raise SystemExit(
            f"Phase 4U batch candidate {candidate_id!r} is not in pose_pass_ids {pass_ids}"
        )
    if not parent_seq:
        raise SystemExit("Phase 4U batch_summary has no binder_sequence")

    model=load_model(HUMAN_CIF)
    target=model["A"]
    binder=model["B"]

    struct_seq,struct_pos=binder_sequence(binder)
    if len(struct_seq)!=len(parent_seq):
        raise RuntimeError(
            f"Binder length mismatch: structure={len(struct_seq)} summary={len(parent_seq)}"
        )

    binder_index={rid:i+1 for i,rid in enumerate(struct_pos)}
    target_res=[r for r in target if r.id[0]==" "]
    binder_res=[r for r in binder if r.id[0]==" "]

    ranked=[]
    for br in binder_res:
        bpos=binder_index[br.id[1]]
        src=parent_seq[bpos-1]
        if src in EXCLUDED_SOURCE_AA:
            continue

        contacts=[]
        for tr in target_res:
            best=999.0
            for ba in br:
                if ba.element=="H":
                    continue
                for ta in tr:
                    if ta.element=="H":
                        continue
                    best=min(best,distance(ba,ta))
            if best<=5.0:
                up=tr.id[1]+HUMAN_START-1
                contacts.append({
                    "target_residue":up,
                    "target_resname":tr.resname,
                    "distance_A":round(best,3),
                })

        if not contacts:
            continue

        contacted={x["target_residue"] for x in contacts}
        basic=sorted(
            x["target_residue"] for x in contacts if x["target_resname"] in BASIC_RESN
        )
        acidic=sorted(
            x["target_residue"] for x in contacts if x["target_resname"] in ACIDIC_RESN
        )
        risk=sorted(contacted & SPECIES_RISK)
        preferred=sorted(contacted & PREFERRED_H)

        score=0.0
        score+=4.0*len(basic)
        score-=3.0*len(acidic)
        score-=2.5*len(risk)
        score+=0.5*len(preferred)
        score+=min(1.5,0.10*len(contacted))

        # His is aromatic; F/Y/W -> H can sometimes preserve interface packing
        # better than replacing a small aliphatic residue. This is only a mild
        # ranking bonus, never evidence of a pH switch.
        if src in {"F","Y","W"}:
            score+=0.75

        ranked.append({
            "candidate_id":candidate_id,
            "binder_position":bpos,
            "source_aa":src,
            "mutation":f"{src}{bpos}H",
            "contacted_target_residues":sorted(contacted),
            "preferred_contacts":preferred,
            "basic_contacts":basic,
            "acidic_contacts":acidic,
            "species_risk_contacts":risk,
            "contact_detail":contacts,
            "heuristic_priority_score":round(score,3),
        })

    ranked.sort(
        key=lambda x:(
            x["heuristic_priority_score"],
            len(x["basic_contacts"]),
            -len(x["acidic_contacts"]),
            -len(x["species_risk_contacts"]),
        ),
        reverse=True,
    )

    variants=[]

    # Single substitutions allow mechanistic attribution.
    for row in ranked[:10]:
        p=row["binder_position"]
        variants.append({
            "parent_id":candidate_id,
            "variant_id":f"{candidate_id}_H{p}",
            "positions":[p],
            "mutations":[row["mutation"]],
            "sequence":mutate(parent_seq,[p]),
            "priority_score":row["heuristic_priority_score"],
            "mechanism":"single_His_interface_scan",
        })

    # Double variants require at least one basic-target contact and exclude
    # combinations dominated by acidic/risk contacts.
    top=ranked[:7]
    for i in range(len(top)):
        for j in range(i+1,len(top)):
            pair=[top[i],top[j]]
            if not any(x["basic_contacts"] for x in pair):
                continue
            if sum(bool(x["acidic_contacts"]) for x in pair)>1:
                continue
            if sum(bool(x["species_risk_contacts"]) for x in pair)>1:
                continue
            pos=sorted(x["binder_position"] for x in pair)
            if abs(pos[0]-pos[1])<=1:
                continue
            variants.append({
                "parent_id":candidate_id,
                "variant_id":f"{candidate_id}_H{pos[0]}_H{pos[1]}",
                "positions":pos,
                "mutations":[x["mutation"] for x in pair],
                "sequence":mutate(parent_seq,pos),
                "priority_score":round(sum(x["heuristic_priority_score"] for x in pair),3),
                "mechanism":"double_His_interface_scan",
            })

    # Triple-His clusters are deliberately small and conservative.
    tri=ranked[:6]
    for i in range(len(tri)):
        for j in range(i+1,len(tri)):
            for k in range(j+1,len(tri)):
                group=[tri[i],tri[j],tri[k]]
                if not any(x["basic_contacts"] for x in group):
                    continue
                if sum(bool(x["acidic_contacts"]) for x in group)>1:
                    continue
                if sum(bool(x["species_risk_contacts"]) for x in group)>1:
                    continue
                pos=sorted(x["binder_position"] for x in group)
                variants.append({
                    "parent_id":candidate_id,
                    "variant_id":f"{candidate_id}_H{pos[0]}_H{pos[1]}_H{pos[2]}",
                    "positions":pos,
                    "mutations":[x["mutation"] for x in group],
                    "sequence":mutate(parent_seq,pos),
                    "priority_score":round(sum(x["heuristic_priority_score"] for x in group),3),
                    "mechanism":"triple_His_interface_cluster_scan",
                })

    # Deduplicate by sequence and retain highest-scoring mechanism.
    best={}
    for v in variants:
        if v["sequence"] not in best or v["priority_score"]>best[v["sequence"]]["priority_score"]:
            best[v["sequence"]]=v
    variants=sorted(best.values(),key=lambda x:x["priority_score"],reverse=True)

    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"interface_histidine_positions.json").write_text(
        json.dumps(ranked,indent=2)+"\n"
    )
    (OUT/"histidine_variants.json").write_text(
        json.dumps(variants,indent=2)+"\n"
    )

    fields=[
        "parent_id","variant_id","positions","mutations","sequence",
        "priority_score","mechanism",
    ]
    with (OUT/"histidine_variants.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields)
        w.writeheader()
        for row in variants:
            x=dict(row)
            x["positions"]=";".join(map(str,x["positions"]))
            x["mutations"]=";".join(x["mutations"])
            w.writerow(x)

    summary={
        "status":"generated",
        "parent_id":candidate_id,
        "parent_sequence":parent_seq,
        "ranked_interface_positions":len(ranked),
        "generated_unique_variants":len(variants),
        "top_variant_ids":[x["variant_id"] for x in variants[:12]],
        "warning":(
            "These are mechanistic pH-switch hypotheses only. They require neutral-pH "
            "human/mouse pose revalidation and explicit pKa/protonation analysis."
        ),
    }
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
