#!/usr/bin/env python3
"""Audit the complete Odin-Multi full-context v4 trajectory.

This script is diagnostic only. It does not select or validate a competition
candidate. It reads Odin's saved trajectory.pickle and reports whether the
full-ectodomain optimisation ever entered a joint human/mouse binding basin.

Usage:
    python scripts/audit_odin_full_context_v4_trajectory.py \
        --root /tmp/odin_v4_output \
        --out analysis/odin_full_context_v4_trajectory_audit
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pickle
from pathlib import Path


def as_float(x):
    try:
        if hasattr(x, "item"):
            return float(x.item())
        return float(x)
    except Exception:
        return None


def find_one(root: Path, name: str) -> Path:
    hits=list(root.rglob(name))
    if not hits:
        raise SystemExit(f"Could not find {name} under {root}")
    if len(hits)>1:
        # Prefer a file beneath 01_designs/t00000.
        preferred=[p for p in hits if "01_designs" in p.parts and "t00000" in p.parts]
        if preferred:
            return preferred[0]
    return hits[0]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",required=True)
    ap.add_argument("--out",required=True)
    args=ap.parse_args()

    root=Path(args.root)
    out=Path(args.out)
    out.parent.mkdir(parents=True,exist_ok=True)

    traj_path=find_one(root,"trajectory.pickle")
    with traj_path.open("rb") as fh:
        traj=pickle.load(fh)

    contexts=[]
    for name,v in traj.items():
        if name.startswith("__") or not isinstance(v,dict):
            continue
        if "i_ptm" in v and "stage" in v and "iteration" in v:
            contexts.append(name)

    if len(contexts)<2:
        raise SystemExit(f"Expected >=2 Odin contexts; found {contexts}")

    # Prefer explicit human/mouse names where available.
    human=next((x for x in contexts if "human" in x.lower()),contexts[0])
    mouse=next((x for x in contexts if "mouse" in x.lower()),contexts[1])

    n=min(
        len(traj[human]["i_ptm"]),
        len(traj[mouse]["i_ptm"]),
        len(traj[human]["stage"]),
        len(traj[human]["iteration"]),
    )
    rows=[]
    for i in range(n):
        hi=as_float(traj[human]["i_ptm"][i])
        mi=as_float(traj[mouse]["i_ptm"][i])
        if hi is None or mi is None:
            continue
        stage=str(traj[human]["stage"][i])
        iteration=traj[human]["iteration"][i]
        try:
            iteration=int(iteration)
        except Exception:
            iteration=str(iteration)
        rows.append({
            "frame":i,
            "iteration":iteration,
            "stage":stage,
            "human_i_ptm":hi,
            "mouse_i_ptm":mi,
            "worst_target_i_ptm":min(hi,mi),
            "mean_target_i_ptm":(hi+mi)/2.0,
        })

    if not rows:
        raise SystemExit("No usable trajectory frames")

    # stage summaries
    stages={}
    for row in rows:
        stages.setdefault(row["stage"],[]).append(row)

    stage_summary={}
    for stage,rr in stages.items():
        best=max(rr,key=lambda x:x["worst_target_i_ptm"])
        stage_summary[stage]={
            "frame_count":len(rr),
            "max_human_i_ptm":max(x["human_i_ptm"] for x in rr),
            "max_mouse_i_ptm":max(x["mouse_i_ptm"] for x in rr),
            "max_worst_target_i_ptm":best["worst_target_i_ptm"],
            "best_joint_frame":best["frame"],
            "best_joint_iteration":best["iteration"],
            "mean_human_i_ptm":sum(x["human_i_ptm"] for x in rr)/len(rr),
            "mean_mouse_i_ptm":sum(x["mouse_i_ptm"] for x in rr)/len(rr),
        }

    ranked=sorted(rows,key=lambda x:x["worst_target_i_ptm"],reverse=True)
    best=ranked[0]
    last=rows[-1]

    summary={
        "trajectory_path":str(traj_path),
        "human_context":human,
        "mouse_context":mouse,
        "frame_count":len(rows),
        "best_joint_frame":best,
        "final_frame":last,
        "stage_summary":stage_summary,
        "top_joint_frames":ranked[:10],
        "diagnostic":{
            "ever_worst_target_ge_0_50":any(x["worst_target_i_ptm"]>=0.50 for x in rows),
            "ever_worst_target_ge_0_60":any(x["worst_target_i_ptm"]>=0.60 for x in rows),
            "ever_both_targets_ge_0_50":any(
                x["human_i_ptm"]>=0.50 and x["mouse_i_ptm"]>=0.50 for x in rows
            ),
            "best_vs_final_delta":best["worst_target_i_ptm"]-last["worst_target_i_ptm"],
        },
    }

    # Optional final auxiliary loss snapshot.
    try:
        aux_path=find_one(root,"auxiliary.pickle")
        with aux_path.open("rb") as fh:
            aux=pickle.load(fh)
        final=aux.get("final_aux_all",{}) if isinstance(aux,dict) else {}
        aux_out={}
        for ctx in (human,mouse):
            v=final.get(ctx,{}) if isinstance(final,dict) else {}
            if not isinstance(v,dict):
                continue
            ctx_out={}
            for key in ("i_ptm","i_pae","i_con","plddt","ptm","pae"):
                if key in v:
                    val=as_float(v[key])
                    if val is not None and math.isfinite(val):
                        ctx_out[key]=val
            aux_out[ctx]=ctx_out
        if aux_out:
            summary["final_auxiliary"]=aux_out
    except Exception as exc:
        summary["auxiliary_note"]=f"Could not parse auxiliary.pickle: {exc}"

    out_json=out.with_suffix(".json")
    out_csv=out.with_suffix(".csv")
    out_md=out.with_suffix(".md")

    out_json.write_text(json.dumps(summary,indent=2)+"\n")

    with out_csv.open("w",newline="") as fh:
        fields=[
            "frame","iteration","stage","human_i_ptm","mouse_i_ptm",
            "worst_target_i_ptm","mean_target_i_ptm",
        ]
        w=csv.DictWriter(fh,fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    md=[
        "# Odin-Multi full-context v4 trajectory audit",
        "",
        f"- Frames: **{len(rows)}**",
        f"- Human context: `{human}`",
        f"- Mouse context: `{mouse}`",
        f"- Best joint frame: **{best['frame']}** ({best['stage']}, iteration {best['iteration']})",
        f"- Best human iPTM: **{best['human_i_ptm']:.3f}**",
        f"- Best mouse iPTM: **{best['mouse_i_ptm']:.3f}**",
        f"- Best worst-target iPTM: **{best['worst_target_i_ptm']:.3f}**",
        f"- Final worst-target iPTM: **{last['worst_target_i_ptm']:.3f}**",
        "",
        "## Diagnostic",
        "",
        f"- Ever worst-target iPTM >= 0.50: **{summary['diagnostic']['ever_worst_target_ge_0_50']}**",
        f"- Ever worst-target iPTM >= 0.60: **{summary['diagnostic']['ever_worst_target_ge_0_60']}**",
        f"- Ever both targets >= 0.50: **{summary['diagnostic']['ever_both_targets_ge_0_50']}**",
        f"- Best-to-final worst-target drop: **{summary['diagnostic']['best_vs_final_delta']:.3f}**",
        "",
        "## Interpretation rule",
        "",
        "- If the trajectory never reaches joint iPTM 0.50, full-context optimisation failed to find a shared basin; extending the same selector is unlikely to rescue it.",
        "- If an early/intermediate frame reaches strong joint iPTM but the selected/final frame is much worse, inspect that frame before changing the generator.",
        "- Independent Phase 4U full-ectodomain Boltz-2 remains the candidate gate.",
        "",
    ]
    out_md.write_text("\n".join(md))
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
