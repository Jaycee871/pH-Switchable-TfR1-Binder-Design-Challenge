#!/usr/bin/env python3
"""Odin-Multi full-ectodomain v4b monomer memory probe for cross-reactive TfR1 binding.

Unlike all previous Odin campaigns, both positive contexts contain the COMPLETE
TfR1 ectodomain. Distal escape surfaces therefore remain present during every
design step rather than being represented by cropped positive/negative contexts.

Human: 11QC chain A, residues 122-760, hotspots 347/350.
Mouse: full predicted chain A, residues 123-763, homologous hotspots 349/352.

The shared 60-aa binder is warm-started from FreeBindCraft v6. One deliberately
bounded trajectory is run to test whether full-context optimisation can retain
the desired homologous apical basin within Kaggle T4 runtime.
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
TMP=Path("/kaggle/temp/tfr1_odin_full_context_v4b_probe")
MINICONDA=TMP/"miniforge"
ODIN=TMP/"odin_multi"
OUT=WORK/"tfr1_odin_full_context_v4b_probe"
OUT.mkdir(parents=True,exist_ok=True)
TMP.mkdir(parents=True,exist_ok=True)

ODIN_COMMIT="e4567e487d43a4de36524f491e51ad1ea9730d1e"
V6_WARM_START_SEQUENCE="SPSVISTTTTTYTYHGKQVTRIRTTINPGEKIEPAPGERQMRGRDWETGYIIVHTYKVNS"

TARGET_URLS={
    "human_tfr1_full_ectodomain.pdb":
      "https://raw.githubusercontent.com/Jaycee871/pH-Switchable-TfR1-Binder-Design-Challenge/main/structures/targets/human_tfr1_full_ectodomain_11QC_chainA.pdb",
    "mouse_tfr1_full_ectodomain.pdb":
      "https://raw.githubusercontent.com/Jaycee871/pH-Switchable-TfR1-Binder-Design-Challenge/main/structures/targets/mouse_tfr1_full_ectodomain_predicted_chainA.pdb",
}

EXPECTED={
    "human_tfr1_full_ectodomain.pdb":(122,760,639),
    "mouse_tfr1_full_ectodomain.pdb":(123,763,641),
}

def run(cmd,cwd=None,env=None,check=True,capture=False):
    print("+"," ".join(map(str,cmd)),flush=True)
    return subprocess.run(
        list(map(str,cmd)),cwd=cwd,env=env,check=check,text=True,
        capture_output=capture,
    )

def ca_residues(path: Path):
    nums=[]
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("ATOM") and len(line)>=54 and line[12:16].strip()=="CA":
            try:
                nums.append(int(line[22:26]))
            except ValueError:
                pass
    return nums

def materialize_targets():
    print("=== Stage 0: downloading and validating full ectodomain targets ===",flush=True)
    for name,url in TARGET_URLS.items():
        dest=OUT/name
        urllib.request.urlretrieve(url,dest)
        nums=ca_residues(dest)
        lo,hi,n=EXPECTED[name]
        if not nums:
            raise SystemExit(f"FATAL: no CA atoms in {name}")
        if min(nums)!=lo or max(nums)!=hi or len(nums)!=n:
            raise SystemExit(
                f"FATAL: {name} CA range/count {min(nums)}-{max(nums)} n={len(nums)}; "
                f"expected {lo}-{hi} n={n}"
            )
        print(f"  OK {name}: CA {min(nums)}-{max(nums)} n={len(nums)}",flush=True)

    h=set(ca_residues(OUT/"human_tfr1_full_ectodomain.pdb"))
    m=set(ca_residues(OUT/"mouse_tfr1_full_ectodomain.pdb"))
    for x in (347,350):
        if x not in h:
            raise SystemExit(f"FATAL: human hotspot {x} absent")
    for x in (349,352):
        if x not in m:
            raise SystemExit(f"FATAL: mouse hotspot {x} absent")

def setup_conda():
    conda=MINICONDA/"bin/conda"
    if conda.exists():
        return conda
    installer=TMP/"miniforge.sh"
    urllib.request.urlretrieve(
        "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh",
        installer,
    )
    run(["bash",str(installer),"-b","-p",str(MINICONDA)])
    return conda

def locate_driver_dirs():
    candidates=[]
    for raw in os.environ.get("LD_LIBRARY_PATH","").split(":"):
        if raw:
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
    env.pop("PYTHONPATH",None)
    env.pop("PYTHONHOME",None)

    if not ODIN.exists():
        run([
            "git","clone","--recurse-submodules",
            "https://github.com/DigBioLab/odin_multi.git",str(ODIN)
        ],env=env)
        run(["git","checkout",ODIN_COMMIT],cwd=ODIN,env=env)
        run(["git","submodule","update","--init","--recursive"],cwd=ODIN,env=env)

    got=run(["git","rev-parse","HEAD"],cwd=ODIN,env=env,capture=True).stdout.strip()
    if got!=ODIN_COMMIT:
        raise RuntimeError(f"Odin commit mismatch: {got} != {ODIN_COMMIT}")
    (OUT/"odin_commit.txt").write_text(got+"\n")

    marker=MINICONDA/"envs/Odin-Multi/bin/python"
    if not marker.exists():
        run([
            "bash","install_odin_multi.sh","--cuda","12.4","--pkg-manager","conda"
        ],cwd=ODIN,env=env)
    return env

def patch_warm_start():
    path=ODIN/"functions/colabdesign_utils.py"
    text=path.read_text()
    old='''    if aux_params_list:
        af_model.prep_offtarget_inputs(aux_params_list)

    af_model._args["gradient_weight"] = float(main_advanced["gradient_weight"])
'''
    new=f'''    if aux_params_list:
        af_model.prep_offtarget_inputs(aux_params_list)

    warm_start_seq = "{V6_WARM_START_SEQUENCE}"
    if len(warm_start_seq) != length:
        raise ValueError(f"Warm-start length {{len(warm_start_seq)}} != requested binder length {{length}}")
    print(f"[WARM START] full-context monomer-probe shared binder initialized from v6 ({{length}} aa)")
    af_model.restart(seed=seed, seq=warm_start_seq)

    af_model._args["gradient_weight"] = float(main_advanced["gradient_weight"])
'''
    if old not in text:
        raise RuntimeError("Verified Odin warm-start insertion point not found")
    path.write_text(text.replace(old,new,1))
    (OUT/"odin_warm_start_patch.txt").write_text(
        "Shared full-context warm-start sequence: "+V6_WARM_START_SEQUENCE+"\n"
    )

def prepare_inputs():
    human=OUT/"human_tfr1_full_ectodomain.pdb"
    mouse=OUT/"mouse_tfr1_full_ectodomain.pdb"

    human_json={
        "binder_name":"human_tfr1_full",
        "starting_pdb":str(human),
        "chains":"A",
        "target_hotspot_residues":"A347,A350",
    }
    mouse_json={
        "binder_name":"mouse_tfr1_full",
        "starting_pdb":str(mouse),
        "chains":"A",
        "target_hotspot_residues":"A349,A352",
    }

    human_loss={
        "role":"target",
        "gradient_weight":1.0,
        "rm_template_seq_design":False,
        "rm_template_sc_design":False,
        "weights_plddt":0.10,
        "weights_pae_intra":0.25,
        "weights_pae_inter":0.20,
        "weights_con_intra":0.90,
        "weights_con_inter":2.00,
        "intra_contact_distance":14.0,
        "inter_contact_distance":8.0,
        "intra_contact_number":2,
        "inter_contact_number":1,
        "intra_contact_num_pos":1000,
        "inter_contact_num_pos":1000,
        "weights_helicity":-0.3,
        "use_i_ptm_loss":True,
        "weights_iptm":0.20,
        "use_ptm_energy_craft_loss":False,
        "weights_ptm_energy_craft":0.0,
        "weights_rg":0.25,
        "use_rg_loss":True,
        "weights_termini_loss":0.10,
        "use_termini_distance_loss":False,
        "mpnn_backprop_num_samples":4,
        "weights_mpnn_structure_confidence":0.0,
        "weights_mpnn_sequence_kl":0.0,
        "weights_mpnn_autoregressive_ce":0.35,
        "mpnn_interface_only":False,
    }
    mouse_loss=dict(human_loss)
    mouse_loss["gradient_weight"]=1.25

    general={
        "lengths":[60,60],
        "use_multimer_design":False,
        "af_params_dir":"",
        "num_recycles_design":1,
        "sample_models":False,
        "ensemble_discrete_stages":False,
        "design_algorithm":"3stage_ste",
        "soft_iterations":3,
        "temporary_iterations":2,
        "hard_iterations":1,
        "greedy_iterations":0,
        "greedy_percentage":5,
        "use_early_stopping":False,
        "optimise_beta":True,
        "optimise_beta_extra_soft":0,
        "optimise_beta_extra_temp":0,
        "optimise_beta_recycles_design":1,
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
        ("human_loss.json",human_loss),
        ("mouse_loss.json",mouse_loss),
        ("general_full_context_v4b_probe.json",general),
    ]:
        p=OUT/name
        p.write_text(json.dumps(obj,indent=2)+"\n")
        paths[name]=p
    return paths

def main():
    t0=time.time()
    run(["nvidia-smi"],check=False)
    materialize_targets()

    conda=setup_conda()
    env=setup_odin(conda)
    patch_warm_start()

    driver_dirs,probe=locate_driver_dirs()
    env["LD_LIBRARY_PATH"]=":".join(driver_dirs+[x for x in env.get("LD_LIBRARY_PATH","").split(":") if x])
    env["CUDA_VISIBLE_DEVICES"]="0"
    (OUT/"driver_probe.txt").write_text(
        "selected_driver_dirs="+":".join(driver_dirs)+"\n"+probe
    )

    run([
        str(conda),"run","--no-capture-output","-n","Odin-Multi",
        "python","-c",
        "import jax; print('JAX devices:',jax.devices()); "
        "assert any(d.platform=='gpu' for d in jax.devices()), 'No JAX GPU'"
    ],cwd=ODIN,env=env)

    paths=prepare_inputs()
    rundir=OUT/"run"

    cmd=[
        str(conda),"run","--no-capture-output","-n","Odin-Multi",
        "python","-u",str(ODIN/"odin_multi.py"),"design",
        "--run-dir",str(rundir),
        "--context",str(paths["human_target.json"]),str(paths["human_loss.json"]),
        "--context",str(paths["mouse_target.json"]),str(paths["mouse_loss.json"]),
        "--advanced",str(paths["general_full_context_v4b_probe.json"]),
        "--base-seed","20261107",
        "--num-designs","1",
    ]
    proc=run(cmd,cwd=ODIN,env=env,check=False)
    (OUT/"design_exit_code.txt").write_text(str(proc.returncode)+"\n")
    if proc.returncode!=0:
        raise SystemExit(proc.returncode)

    sel=run([
        str(conda),"run","--no-capture-output","-n","Odin-Multi",
        "python","-u",str(ODIN/"odin_multi.py"),"select",
        "--run-dir",str(rundir),"--method","best_i_ptm"
    ],cwd=ODIN,env=env,check=False)
    (OUT/"select_exit_code.txt").write_text(str(sel.returncode)+"\n")
    if sel.returncode!=0:
        raise SystemExit(sel.returncode)

    selection=rundir/"02_selections/best_i_ptm/selection.csv"
    fasta=rundir/"02_selections/best_i_ptm/sequences.fasta"
    if selection.exists():
        shutil.copy2(selection,OUT/"selection.csv")
    if fasta.exists():
        shutil.copy2(fasta,OUT/"selected_sequences.fasta")

    summary={
        "contexts":[
            {"name":"human_tfr1_full","range":[122,760],"hotspots":[347,350],"gradient_weight":1.0},
            {"name":"mouse_tfr1_full","range":[123,763],"hotspots":[349,352],"gradient_weight":1.25},
        ],
        "binder_length":60,
        "trajectory_count":1,
        "use_multimer_design":False,
        "probe_purpose":"test whether full-ectodomain gradients fit on Kaggle T4 with monomer/PTM binder model",
        "selection_method":"best_i_ptm",
        "warm_start":"FreeBindCraft v6 sequence",
        "design_iterations":{"logits":3,"soft":2,"ste":1,"greedy":0},
        "odin_commit":ODIN_COMMIT,
        "selection_csv_exists":selection.exists(),
        "selection_fasta_exists":fasta.exists(),
        "elapsed_min":round((time.time()-t0)/60,1),
    }
    (OUT/"smoke_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)

if __name__=="__main__":
    main()
