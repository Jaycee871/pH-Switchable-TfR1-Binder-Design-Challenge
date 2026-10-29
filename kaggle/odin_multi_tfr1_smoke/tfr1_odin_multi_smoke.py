#!/usr/bin/env python3
"""Kaggle T4 Odin-Multi cross-reactivity smoke for human+mouse TfR1.

This is the first true multi-context backbone+sequence design stage in the project.
The same binder sequence is optimized simultaneously against:
  - human TfR1 apical domain, hotspots 343/344/347/350
  - mouse TfR1 homologous apical domain, hotspots 345/346/349/352

One 60-aa trajectory is generated and selected by best_i_ptm (maximize the weakest
on-target iPTM). Independent full-human/full-mouse Boltz-2 validation remains
outside Odin-Multi.
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
TMP=Path("/kaggle/temp/tfr1_odin_multi")
MINICONDA=TMP/"miniconda"
ODIN=TMP/"odin_multi"
OUT=WORK/"tfr1_odin_multi_crossreactivity_smoke"
HERE=Path(__file__).resolve().parent
OUT.mkdir(parents=True,exist_ok=True)
TMP.mkdir(parents=True,exist_ok=True)

def run(cmd,cwd=None,env=None,check=True,capture=False):
    print("+"," ".join(map(str,cmd)),flush=True)
    return subprocess.run(
        list(map(str,cmd)),cwd=cwd,env=env,check=check,text=True,
        capture_output=capture,
    )

def setup_conda():
    conda=MINICONDA/"bin/conda"
    if conda.exists():
        return conda
    installer=TMP/"miniconda.sh"
    urllib.request.urlretrieve(
        "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh",
        installer,
    )
    run(["bash",str(installer),"-b","-p",str(MINICONDA)])
    return conda

def locate_driver_dirs():
    candidates=[]
    for raw in os.environ.get("LD_LIBRARY_PATH","").split(":"):
        if not raw:
            continue
        p=Path(raw)
        if p.exists() and any((p/n).exists() for n in ("libcuda.so.1","libcuda.so")):
            candidates.append(str(p))
    probe=subprocess.run(
        ["bash","-lc",
         "ldconfig -p 2>/dev/null | "
         "grep -E 'libcuda\\.so\\.1|libnvidia-ptxjitcompiler\\.so\\.1|libnvidia-nvvm\\.so' || true"],
        capture_output=True,text=True,check=False,
    )
    print("Driver probe:\n"+probe.stdout,flush=True)
    for line in probe.stdout.splitlines():
        if "=>" in line:
            p=Path(line.split("=>",1)[1].strip()).parent
            if p.exists():
                candidates.append(str(p))
    for fallback in ("/usr/local/nvidia/lib64","/usr/local/nvidia/lib","/usr/lib/x86_64-linux-gnu"):
        p=Path(fallback)
        if p.exists() and any((p/n).exists() for n in ("libcuda.so.1","libcuda.so")):
            candidates.append(str(p))
    out=[]
    for p in candidates:
        if p not in out:
            out.append(p)
    if not out:
        raise RuntimeError("Could not locate Kaggle NVIDIA driver libraries")
    return out,probe.stdout

def setup_odin(conda):
    env=os.environ.copy()
    env["PATH"]=f"{MINICONDA}/bin:"+env.get("PATH","")
    if not ODIN.exists():
        run([
            "git","clone","--recurse-submodules","--depth","1",
            "https://github.com/DigBioLab/odin_multi.git",str(ODIN)
        ],env=env)
    marker=MINICONDA/"envs/Odin-Multi/bin/python"
    if not marker.exists():
        run([
            "bash","install_odin_multi.sh",
            "--cuda","12.4",
            "--pkg-manager","conda",
        ],cwd=ODIN,env=env)
    else:
        print("Odin-Multi conda environment already present.",flush=True)
    return env

def prepare_inputs():
    human=OUT/"human_tfr1_apical_189_383.pdb"
    mouse=OUT/"mouse_tfr1_apical_homolog.pdb"
    shutil.copy2(HERE/"human_tfr1_apical_189_383.pdb",human)
    shutil.copy2(HERE/"mouse_tfr1_apical_homolog.pdb",mouse)

    human_json={
        "binder_name":"human_tfr1_apical",
        "starting_pdb":str(human),
        "chains":"A",
        "target_hotspot_residues":"A343,A344,A347,A350",
    }
    mouse_json={
        "binder_name":"mouse_tfr1_apical",
        "starting_pdb":str(mouse),
        "chains":"A",
        "target_hotspot_residues":"A345,A346,A349,A352",
    }

    target_loss={
        "role":"target",
        "gradient_weight":1.0,
        "rm_template_seq_design":False,
        "rm_template_sc_design":False,
        "weights_plddt":0.1,
        "weights_pae_intra":0.4,
        "weights_pae_inter":0.1,
        "weights_con_intra":1.0,
        "weights_con_inter":1.25,
        "intra_contact_distance":14.0,
        "inter_contact_distance":8.0,
        "intra_contact_number":2,
        "inter_contact_number":1,
        "intra_contact_num_pos":1000,
        "inter_contact_num_pos":1000,
        "weights_helicity":-0.3,
        "use_i_ptm_loss":True,
        "weights_iptm":0.05,
        "use_ptm_energy_craft_loss":False,
        "weights_ptm_energy_craft":0.0,
        "weights_rg":0.3,
        "use_rg_loss":True,
        "weights_termini_loss":0.1,
        "use_termini_distance_loss":False,
        "mpnn_backprop_num_samples":4,
        "weights_mpnn_structure_confidence":0.0,
        "weights_mpnn_sequence_kl":0.0,
        "weights_mpnn_autoregressive_ce":1.0,
        "mpnn_interface_only":False,
    }

    general={
        "lengths":[60,60],
        "use_multimer_design":True,
        "af_params_dir":"",
        "num_recycles_design":1,
        "sample_models":True,
        "ensemble_discrete_stages":False,
        "design_algorithm":"3stage_ste",
        "soft_iterations":30,
        "temporary_iterations":20,
        "hard_iterations":10,
        "greedy_iterations":5,
        "greedy_percentage":5,
        "use_early_stopping":False,
        "optimise_beta":True,
        "optimise_beta_extra_soft":0,
        "optimise_beta_extra_temp":0,
        "optimise_beta_recycles_design":3,
        "mpnn_weights":"original",
        "mpnn_use_ligand_context":False,
        "random_helicity":False,
        "dssp_path":"",
        "save_design_animations":False,
        "save_design_trajectory_plots":False,
        "remove_unrelaxed_trajectory":False,
        "omit_AAs":"C",
    }

    paths={}
    for name,obj in [
        ("human_target.json",human_json),
        ("mouse_target.json",mouse_json),
        ("target_loss.json",target_loss),
        ("general_smoke.json",general),
    ]:
        p=OUT/name
        p.write_text(json.dumps(obj,indent=2)+"\n")
        paths[name]=p
    return paths

def main():
    t0=time.time()
    run(["nvidia-smi"],check=False)
    print("Host Python",sys.version,flush=True)

    conda=setup_conda()
    env=setup_odin(conda)

    driver_dirs,probe=locate_driver_dirs()
    env["LD_LIBRARY_PATH"]=":".join(driver_dirs + [x for x in env.get("LD_LIBRARY_PATH","").split(":") if x])
    env["CUDA_VISIBLE_DEVICES"]="0"
    (OUT/"driver_probe.txt").write_text(
        "selected_driver_dirs="+":".join(driver_dirs)+"\n"+probe
    )

    # Verify actual JAX GPU visibility rather than relying on nvidia-smi inside conda.
    run([
        str(conda),"run","--no-capture-output","-n","Odin-Multi",
        "python","-c",
        "import jax; print('JAX devices:',jax.devices()); "
        "assert any(d.platform=='gpu' for d in jax.devices()), 'No JAX GPU'"
    ],cwd=ODIN,env=env)

    paths=prepare_inputs()
    rundir=OUT/"run"

    design_cmd=[
        str(conda),"run","--no-capture-output","-n","Odin-Multi",
        "python","-u",str(ODIN/"odin_multi.py"),"design",
        "--run-dir",str(rundir),
        "--context",str(paths["human_target.json"]),str(paths["target_loss.json"]),
        "--context",str(paths["mouse_target.json"]),str(paths["target_loss.json"]),
        "--advanced",str(paths["general_smoke.json"]),
        "--base-seed","20261103",
        "--num-designs","1",
    ]
    proc=run(design_cmd,cwd=ODIN,env=env,check=False)
    (OUT/"design_exit_code.txt").write_text(str(proc.returncode)+"\n")
    if proc.returncode!=0:
        raise SystemExit(proc.returncode)

    select_cmd=[
        str(conda),"run","--no-capture-output","-n","Odin-Multi",
        "python","-u",str(ODIN/"odin_multi.py"),"select",
        "--run-dir",str(rundir),
        "--method","best_i_ptm",
    ]
    sel=run(select_cmd,cwd=ODIN,env=env,check=False)
    (OUT/"select_exit_code.txt").write_text(str(sel.returncode)+"\n")
    if sel.returncode!=0:
        raise SystemExit(sel.returncode)

    selection=rundir/"02_selections/best_i_ptm/selection.csv"
    fasta=rundir/"02_selections/best_i_ptm/sequences.fasta"
    summary={
        "contexts":[
            {"name":"human_tfr1_apical","hotspots":[343,344,347,350]},
            {"name":"mouse_tfr1_apical","hotspots":[345,346,349,352]},
        ],
        "binder_length":60,
        "trajectory_count":1,
        "selection_method":"best_i_ptm",
        "design_iterations":{"soft":30,"temporary":20,"hard":10,"greedy":5},
        "selection_csv_exists":selection.exists(),
        "selection_fasta_exists":fasta.exists(),
        "elapsed_min":round((time.time()-t0)/60,1),
    }
    if selection.exists():
        shutil.copy2(selection,OUT/"selection.csv")
    if fasta.exists():
        shutil.copy2(fasta,OUT/"selected_sequences.fasta")
    (OUT/"smoke_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)

    shutil.make_archive(
        str(WORK/"tfr1_odin_multi_crossreactivity_smoke"),
        "zip",
        root_dir=OUT,
    )

if __name__=="__main__":
    main()
