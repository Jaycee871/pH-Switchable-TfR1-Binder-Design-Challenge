#!/usr/bin/env python3
"""Kaggle T4 smoke test for FreeBindCraft on the TfR1 apical domain.

Goal: test whether a small TfR1 target (11QC chain A, residues 189-383)
plus a 55-75 aa miniprotein can run FreeBindCraft on Kaggle T4.

This is intentionally NOT a production design campaign.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

WORK=Path("/kaggle/working")
TMP=Path("/kaggle/temp/tfr1_freebindcraft")
MINICONDA=TMP/"miniconda"
FBC=TMP/"FreeBindCraft"
OUT=WORK/"tfr1_freebindcraft_smoke_v1"
OUT.mkdir(parents=True,exist_ok=True)
TMP.mkdir(parents=True,exist_ok=True)

def run(cmd,cwd=None,env=None,check=True):
    print("+"," ".join(map(str,cmd)),flush=True)
    return subprocess.run(list(map(str,cmd)),cwd=cwd,env=env,check=check,text=True)

def setup_conda():
    conda=MINICONDA/"bin/conda"
    if conda.exists():
        return conda
    installer=TMP/"miniconda.sh"
    print("Downloading Miniconda",flush=True)
    urllib.request.urlretrieve(
        "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh",
        installer,
    )
    run(["bash",str(installer),"-b","-p",str(MINICONDA)])
    return conda

def setup_freebindcraft(conda):
    env=os.environ.copy()
    env["PATH"]=f"{MINICONDA}/bin:"+env.get("PATH","")
    if not FBC.exists():
        run([
            "git","clone","--depth","1","--branch","master",
            "https://github.com/cytokineking/FreeBindCraft.git",str(FBC)
        ],env=env)
    marker=FBC/"params/params_model_1_multimer_v3.npz"
    if not marker.exists():
        run([
            "bash","install_bindcraft.sh",
            "--pkg_manager","conda",
            "--env-name","BindCraft",
            "--no-pyrosetta",
            "--skip-gpu-check"
        ],cwd=FBC,env=env)
    else:
        print("FreeBindCraft environment/weights already present.",flush=True)
    return env

def prepare_target():
    full=TMP/"11QC.pdb"
    if not full.exists():
        urllib.request.urlretrieve("https://files.rcsb.org/download/11QC.pdb",full)
    kept=[];res=set()
    for line in full.read_text().splitlines():
        if not line.startswith("ATOM  ") or len(line)<54 or line[21]!="A":
            continue
        try:r=int(line[22:26])
        except ValueError:continue
        if 189<=r<=383:
            kept.append(line);res.add(r)
    if len(res)<180:
        raise RuntimeError(f"Unexpected apical crop: {len(res)} residues")
    target=WORK/"tfr1_apical_189_383_chainA.pdb"
    target.write_text("\n".join(kept)+"\nTER\nEND\n")
    return target

def prepare_settings(target):
    design_dir=OUT/"designs"
    settings={
        "design_path":str(design_dir),
        "binder_name":"TFR1_FBC343",
        "starting_pdb":str(target),
        "chains":"A",
        "target_hotspot_residues":"343,344,347,350",
        "lengths":[55,75],
        "number_of_final_designs":1,
    }
    settings_path=OUT/"target_settings.json"
    settings_path.write_text(json.dumps(settings,indent=2)+"\n")

    advanced=json.loads((FBC/"settings_advanced/default_4stage_multimer.json").read_text())
    # Smoke-test budget only. Production settings must be restored after this
    # environment/geometry test succeeds.
    advanced.update({
        "soft_iterations":25,
        "temporary_iterations":15,
        "hard_iterations":3,
        "greedy_iterations":5,
        "save_design_animations":False,
        "save_design_trajectory_plots":False,
        "zip_animations":False,
        "zip_plots":False,
        "num_seqs":4,
        "max_mpnn_sequences":1,
        "max_trajectories":4,
        "enable_rejection_check":False,
        "weights_iptm":0.10,
        "weights_rg":0.40,
    })
    advanced_path=OUT/"advanced_smoke.json"
    advanced_path.write_text(json.dumps(advanced,indent=2)+"\n")
    return settings_path,advanced_path

def summarize():
    design_dir=OUT/"designs"
    files=[str(p.relative_to(OUT)) for p in sorted(design_dir.glob("**/*")) if p.is_file()] if design_dir.exists() else []
    stats={}
    for name in ["trajectory_stats.csv","mpnn_design_stats.csv","final_design_stats.csv","failure_csv.csv"]:
        p=design_dir/name
        if p.exists():
            lines=p.read_text(errors="replace").splitlines()
            stats[name]={"rows":max(0,len(lines)-1),"bytes":p.stat().st_size}
            shutil.copy2(p,OUT/name)
    accepted=list(design_dir.glob("Accepted/**/*.pdb")) if design_dir.exists() else []
    summary={
        "target":"11QC chain A residues 189-383",
        "hotspots":[343,344,347,350],
        "binder_lengths":[55,75],
        "smoke_max_trajectories":4,
        "accepted_pdb_count":len(accepted),
        "stats":stats,
        "output_file_count":len(files),
        "environment_test_completed":True,
        "production_ready_if_no_oom":True,
    }
    (OUT/"smoke_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)

if __name__=="__main__":
    t0=time.time()
    run(["nvidia-smi"],check=False)
    print("Host Python",sys.version,flush=True)
    conda=setup_conda()
    env=setup_freebindcraft(conda)
    # The installer GPU probe relies on nvidia-smi visibility inside the conda
    # environment, which can fail on Kaggle even when the host GPU is healthy.
    # Verify the actual JAX backend directly before spending hallucination time.
    gpu_probe=[
        str(conda),"run","--no-capture-output","-n","BindCraft",
        "python","-c",
        "import jax; print('JAX devices:', jax.devices()); "
        "assert any(d.platform == 'gpu' for d in jax.devices()), 'JAX GPU backend unavailable'"
    ]
    run(gpu_probe,cwd=FBC,env=env)

    target=prepare_target()
    settings,advanced=prepare_settings(target)

    cmd=[
        str(conda),"run","--no-capture-output","-n","BindCraft",
        "python","-u",str(FBC/"bindcraft.py"),
        "--settings",str(settings),
        "--filters",str(FBC/"settings_filters/relaxed_filters.json"),
        "--advanced",str(advanced),
        "--no-pyrosetta",
        "--no-plots",
        "--no-animations",
        "--verbose",
        "--rank-by","ipSAE",
    ]
    proc=run(cmd,cwd=FBC,env=env,check=False)
    (OUT/"bindcraft_exit_code.txt").write_text(str(proc.returncode)+"\n")
    summarize()
    shutil.make_archive(str(WORK/"tfr1_freebindcraft_smoke_v1"),"zip",root_dir=OUT)
    print(f"elapsed_min={(time.time()-t0)/60:.1f}",flush=True)
    if proc.returncode!=0:
        raise SystemExit(proc.returncode)
