#!/usr/bin/env python3
"""Phase 5C: PROPKA-based mechanistic ranking of neutral-pose His variants.

This stage does NOT estimate binding affinity at pH 6.0. It ranks engineered
histidines by whether their predicted pKa makes protonation change strongly
between pH 7.4 and pH 6.0 while the variant has already passed the independent
neutral-pH human/mouse structural gate.

For each engineered His:
    f_H+(pH) = 1 / (1 + 10**(pH - pKa))

A useful switch candidate should have a large:
    delta_protonation = f_H+(6.0) - f_H+(7.4)

PROPKA is used only as a mechanistic proxy. Experimental binding remains decisive.
"""

from __future__ import annotations

import csv
import json
import math
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO

ROOT=Path(__file__).resolve().parents[1]
SCREEN=ROOT/"predictions/odin_v4_histidine_neutral_screen"
BATCH=SCREEN/"batch_summary.json"
PASS_IDS=SCREEN/"neutral_pose_pass_ids.json"
OUT=ROOT/"analysis/odin_v4_histidine_propka"

PKA_RE=re.compile(
    r"^\s*HIS\s+(-?\d+)\s+([A-Za-z0-9])\s+(-?\d+(?:\.\d+)?)\b"
)


def fraction_protonated(pka:float,ph:float)->float:
    return 1.0/(1.0+10.0**(ph-pka))


def cif_to_pdb(cif:Path,pdb:Path):
    st=MMCIFParser(QUIET=True).get_structure(cif.stem,str(cif))
    io=PDBIO()
    io.set_structure(st)
    io.save(str(pdb))


def run_propka(pdb:Path,work:Path)->Path:
    cmd=["propka3",str(pdb),"--quiet"]
    proc=subprocess.run(cmd,cwd=work,text=True,capture_output=True)
    (work/"propka_stdout.txt").write_text(proc.stdout or "")
    (work/"propka_stderr.txt").write_text(proc.stderr or "")
    if proc.returncode!=0:
        raise RuntimeError(
            f"PROPKA failed for {pdb.name}: return={proc.returncode}\n{proc.stderr}"
        )
    pka=work/(pdb.stem+".pka")
    if not pka.exists():
        hits=list(work.glob("*.pka"))
        if len(hits)==1:
            pka=hits[0]
        else:
            raise RuntimeError(f"No unique .pka output in {work}")
    return pka


def parse_histidine_pkas(path:Path):
    out={}
    for line in path.read_text(errors="replace").splitlines():
        m=PKA_RE.match(line)
        if not m:
            continue
        resid=int(m.group(1))
        chain=m.group(2)
        pka=float(m.group(3))
        out[(chain,resid)]=pka
    return out


def main():
    if not PASS_IDS.exists() or not BATCH.exists():
        raise SystemExit("Phase 5B neutral-screen outputs are missing")

    pass_ids=json.loads(PASS_IDS.read_text()).get("neutral_pose_pass_ids",[])
    if not pass_ids:
        print(json.dumps({
            "status":"skipped",
            "reason":"No neutral-pose His variants passed Phase 5B",
        },indent=2))
        return

    data=json.loads(BATCH.read_text())
    candidates=data.get("candidates",[])
    by_id={x["variant_id"]:x for x in candidates}

    OUT.mkdir(parents=True,exist_ok=True)
    results=[]

    for vid in pass_ids:
        row=by_id.get(vid)
        if not row:
            raise RuntimeError(f"Neutral-pass ID {vid} absent from Phase 5B batch")
        cif=SCREEN/vid/"human/structure_1.cif"
        if not cif.exists():
            raise FileNotFoundError(cif)

        positions=[int(x) for x in row.get("positions",[])]
        vout=OUT/vid
        vout.mkdir(parents=True,exist_ok=True)

        pdb=vout/"human_complex.pdb"
        cif_to_pdb(cif,pdb)

        pka_file=run_propka(pdb,vout)
        shutil.copy2(pka_file,vout/"human_complex.pka.txt")
        pkas=parse_histidine_pkas(pka_file)

        engineered=[]
        missing=[]
        for pos in positions:
            key=("B",pos)
            if key not in pkas:
                missing.append(pos)
                continue
            pka=pkas[key]
            f6=fraction_protonated(pka,6.0)
            f74=fraction_protonated(pka,7.4)
            delta=f6-f74
            engineered.append({
                "binder_position":pos,
                "predicted_pKa":round(pka,3),
                "fraction_protonated_pH6":round(f6,4),
                "fraction_protonated_pH7_4":round(f74,4),
                "delta_protonation":round(delta,4),
                "switch_window_5_5_to_7_0":bool(5.5<=pka<=7.0),
            })

        if not engineered:
            mechanistic=0.0
            mean_delta=0.0
            in_window=0
        else:
            mean_delta=sum(x["delta_protonation"] for x in engineered)/len(engineered)
            in_window=sum(x["switch_window_5_5_to_7_0"] for x in engineered)
            # Prefer strong protonation change and multiple engineered His in the
            # switchable pKa window, while retaining high neutral-pH cross-species confidence.
            neutral_floor=min(
                float(row["species"]["human"]["protein_iptm"] or 0),
                float(row["species"]["mouse"]["protein_iptm"] or 0),
            )
            mechanistic=2.0*mean_delta + 0.35*in_window + 0.5*neutral_floor

        result={
            "variant_id":vid,
            "sequence":row["sequence"],
            "mutations":row.get("mutations",[]),
            "engineered_positions":positions,
            "engineered_histidines":engineered,
            "missing_engineered_pKa_positions":missing,
            "neutral_human_ipTM":row["species"]["human"]["protein_iptm"],
            "neutral_mouse_ipTM":row["species"]["mouse"]["protein_iptm"],
            "homologous_jaccard":row["cross_species_epitope"]["jaccard"],
            "mean_engineered_delta_protonation":round(mean_delta,4),
            "engineered_His_in_switch_window":in_window,
            "mechanistic_rank_score":round(mechanistic,4),
        }
        (vout/"propka_summary.json").write_text(json.dumps(result,indent=2)+"\n")
        results.append(result)

    results.sort(
        key=lambda x:(
            x["mechanistic_rank_score"],
            x["engineered_His_in_switch_window"],
            x["mean_engineered_delta_protonation"],
            min(float(x["neutral_human_ipTM"]),float(x["neutral_mouse_ipTM"])),
        ),
        reverse=True,
    )

    recommended=results[:min(5,len(results))]
    summary={
        "neutral_pose_pass_count":len(results),
        "recommended_variant_ids":[x["variant_id"] for x in recommended],
        "ranked":results,
        "interpretation":(
            "Ranking combines neutral-pH cross-species structural retention with "
            "PROPKA-predicted engineered-His protonation changes between pH 7.4 and 6.0. "
            "It is not a pH-dependent binding free-energy calculation."
        ),
    }
    (OUT/"batch_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (OUT/"recommended_for_pH_testing.json").write_text(
        json.dumps({"variant_ids":[x["variant_id"] for x in recommended]},indent=2)+"\n"
    )

    with (OUT/"ranking.csv").open("w",newline="") as fh:
        fields=[
            "variant_id","neutral_human_ipTM","neutral_mouse_ipTM",
            "homologous_jaccard","engineered_His_in_switch_window",
            "mean_engineered_delta_protonation","mechanistic_rank_score",
        ]
        w=csv.DictWriter(fh,fieldnames=fields)
        w.writeheader()
        for x in results:
            w.writerow({k:x[k] for k in fields})

    print(json.dumps({
        "neutral_pose_pass_count":len(results),
        "recommended_variant_ids":[x["variant_id"] for x in recommended],
        "top_scores":[{
            "id":x["variant_id"],
            "score":x["mechanistic_rank_score"],
            "delta":x["mean_engineered_delta_protonation"],
            "window_His":x["engineered_His_in_switch_window"],
        } for x in recommended],
        "warning":"PROPKA ranking is mechanistic prioritization, not proof of pH-selective binding.",
    },indent=2))


if __name__=="__main__":
    main()
