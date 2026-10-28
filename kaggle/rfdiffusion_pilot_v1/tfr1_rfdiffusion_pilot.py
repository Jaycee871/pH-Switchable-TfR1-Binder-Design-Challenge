#!/usr/bin/env python3
"""Kaggle GPU runner for TfR1 RFdiffusion pilot v1.

Generation uses the human TfR1 apical domain (11QC chain A, residues 189-383)
rather than the full 1,280-aa ectodomain dimer. Full human/mouse ectodomains
are reserved for downstream cross-species complex screening.

This stage generates backbones only. It does not establish affinity,
cross-reactivity, or pH selectivity.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

WORK = Path("/kaggle/working")
TMP = Path("/kaggle/temp/tfr1_rf_runtime")
TMP.mkdir(parents=True, exist_ok=True)
VENV = TMP / "rfd_venv"
VPY = VENV / "bin" / "python"
RFD = TMP / "RFdiffusion"
TARGET_FULL = TMP / "11QC.pdb"
TARGET = TMP / "11QC_chainA_apical_189_383.pdb"
OUT = WORK / "tfr1_rfdiffusion_pilot_v1"
OUT.mkdir(parents=True, exist_ok=True)

FAMILIES = {
    "cross_species_loop_5": ["A202", "A206", "A209", "A213", "A214"],
    "cross_species_edge_5": ["A202", "A209", "A214", "A219", "A346"],
    "functional_patch_mixed_5": ["A202", "A208", "A209", "A348", "A371"],
}

CONTIG = "A189-383/0 70-100"


def run(cmd, cwd=None, env=None, check=True):
    print("+", " ".join(map(str, cmd)), flush=True)
    return subprocess.run(
        list(map(str, cmd)), cwd=cwd, env=env, check=check, text=True
    )


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_verified(url: str, path: Path, expected_sha256: str) -> None:
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.unlink(missing_ok=True)
    req = urllib.request.Request(
        url, headers={"User-Agent": "tfr1-rfdiffusion-kaggle/1.0"}
    )
    with urllib.request.urlopen(req, timeout=300) as response, tmp.open("wb") as out:
        while True:
            chunk = response.read(8 * 1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    got = sha256_file(tmp)
    if got.lower() != expected_sha256.lower():
        tmp.unlink(missing_ok=True)
        raise RuntimeError(
            f"Checkpoint SHA256 mismatch: expected {expected_sha256}, got {got}"
        )
    tmp.replace(path)


def download_checkpoint_with_fallback(path: Path) -> None:
    expected = "76e4e260aefee3b582bd76b77ab95d2592e64f00c51bf344968ab9239f3250bc"
    urls = [
        (
            "https://huggingface.co/GlandVergil/RFdiffusion/resolve/"
            "3cdaa7d9e22dbdf085abbd16f17b4dc31995ce4d/"
            "models/Complex_base_ckpt.pt"
        ),
        (
            "https://huggingface.co/datasets/LevinHarness/rfpeptides-assets/"
            "resolve/main/models/rfpeptides/Complex_base_ckpt.pt"
        ),
    ]

    if path.exists() and sha256_file(path).lower() == expected.lower():
        print("RFdiffusion checkpoint already present and verified.")
        return
    path.unlink(missing_ok=True)

    errors = []
    for url in urls:
        try:
            print("Downloading RFdiffusion checkpoint from:", url)
            download_verified(url, path, expected)
            print("Checkpoint SHA256 verified:", expected)
            return
        except Exception as exc:
            path.unlink(missing_ok=True)
            errors.append(f"{url}: {exc!r}")
            print("Checkpoint source failed:", repr(exc))
    raise RuntimeError("All checkpoint mirrors failed:\n" + "\n".join(errors))


def setup_environment():
    print("=== Environment ===")
    run(["nvidia-smi"], check=False)
    print("Host Python:", sys.version)

    uv_bin = shutil.which("uv")
    if uv_bin is None:
        run([sys.executable, "-m", "pip", "install", "--quiet", "uv"])
        uv_bin = shutil.which("uv")
    if uv_bin is None:
        raise RuntimeError("uv executable not found")

    run([uv_bin, "python", "install", "3.11"])
    if not VENV.exists():
        run([uv_bin, "venv", str(VENV), "--python", "3.11"])
    if not VPY.exists():
        raise RuntimeError(f"Python 3.11 virtualenv was not created: {VPY}")

    run([
        uv_bin, "pip", "install", "--python", str(VPY),
        "torch==2.4.1",
        "--index-url", "https://download.pytorch.org/whl/cu124",
    ])
    run([
        uv_bin, "pip", "install", "--python", str(VPY),
        "jedi", "omegaconf", "hydra-core", "icecream", "pyrsistent",
        "pynvml", "decorator", "opt_einsum_fx", "biopython",
        "numpy<2", "scipy", "pandas", "networkx", "requests", "psutil",
        "tqdm", "torchdata==0.8.0",
    ])
    run([
        uv_bin, "pip", "install", "--python", str(VPY),
        "git+https://github.com/NVIDIA/dllogger#egg=dllogger",
    ])
    run([
        uv_bin, "pip", "install", "--python", str(VPY),
        "--no-deps", "dgl",
        "--find-links", "https://data.dgl.ai/wheels/torch-2.4/cu124/repo.html",
    ])
    run([
        uv_bin, "pip", "install", "--python", str(VPY),
        "--no-deps", "e3nn==0.5.5", "opt_einsum_fx",
    ])

    if not RFD.exists():
        run([
            "git", "clone", "--depth", "1",
            "https://github.com/sokrypton/RFdiffusion.git", str(RFD)
        ])
    run([
        uv_bin, "pip", "install", "--python", str(VPY), "."
    ], cwd=RFD / "env" / "SE3Transformer")

    models = RFD / "models"
    models.mkdir(exist_ok=True)
    download_checkpoint_with_fallback(models / "Complex_base_ckpt.pt")

    env = os.environ.copy()
    env["DGLBACKEND"] = "pytorch"
    run([
        str(VPY), "-c",
        "import torch, dgl, e3nn; "
        "print('torch', torch.__version__, 'cuda', torch.cuda.is_available()); "
        "print('gpu', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'none'); "
        "print('dgl', dgl.__version__, 'e3nn', e3nn.__version__)"
    ], env=env)


def download_and_crop_target():
    if not TARGET_FULL.exists():
        urllib.request.urlretrieve(
            "https://files.rcsb.org/download/11QC.pdb", TARGET_FULL
        )

    kept = []
    residues = set()
    for line in TARGET_FULL.read_text().splitlines():
        if not line.startswith("ATOM  "):
            continue
        chain = line[21].strip()
        if chain != "A":
            continue
        try:
            resi = int(line[22:26])
        except ValueError:
            continue
        if 189 <= resi <= 383:
            kept.append(line)
            residues.add(resi)

    required_hotspots = {
        int(token[1:]) for values in FAMILIES.values() for token in values
    }
    missing = sorted(required_hotspots - residues)
    if missing:
        raise RuntimeError(f"11QC crop is missing required hotspot residues: {missing}")
    if len(residues) < 180:
        raise RuntimeError(
            f"11QC apical crop unexpectedly contains only {len(residues)} residues"
        )

    TARGET.write_text("\n".join(kept) + "\nTER\nEND\n")
    print(
        "Target crop:", TARGET,
        "residues=", len(residues),
        "range=", min(residues), max(residues),
        "bytes=", TARGET.stat().st_size,
    )


def run_diffusion():
    env = os.environ.copy()
    env["DGLBACKEND"] = "pytorch"

    for family, hotspots in FAMILIES.items():
        family_dir = OUT / family
        family_dir.mkdir(parents=True, exist_ok=True)
        prefix = family_dir / "design"

        cmd = [
            str(VPY),
            str(RFD / "run_inference.py"),
            f"inference.input_pdb={TARGET}",
            f"inference.output_prefix={prefix}",
            "inference.num_designs=2",
            "inference.write_trajectory=False",
            "denoiser.noise_scale_ca=0.5",
            "denoiser.noise_scale_frame=0.5",
            f"contigmap.contigs=[{CONTIG}]",
            f"ppi.hotspot_res=[{','.join(hotspots)}]",
        ]
        print("\n=== RFdiffusion:", family, "===")
        run(cmd, cwd=RFD, env=env)


def parse_atoms(pdb_path: Path):
    atoms = []
    with pdb_path.open() as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom = line[12:16].strip()
            if atom.startswith("H"):
                continue
            chain = line[21].strip() or "_"
            try:
                resi = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            atoms.append((chain, resi, atom, x, y, z))
    return atoms


def distance(a, b):
    return math.sqrt(
        (a[3]-b[3])**2 + (a[4]-b[4])**2 + (a[5]-b[5])**2
    )


def screen_design(pdb_path: Path, family: str):
    atoms = parse_atoms(pdb_path)
    residue_sets = {}
    for a in atoms:
        residue_sets.setdefault(a[0], set()).add(a[1])
    if len(residue_sets) < 2:
        return {
            "family": family,
            "pdb": str(pdb_path.relative_to(WORK)),
            "geometry_pass": False,
            "screen_note": f"Expected target+binder chains, found {residue_sets}",
        }

    target_chain = max(residue_sets, key=lambda c: len(residue_sets[c]))
    binder_chain = min(residue_sets, key=lambda c: len(residue_sets[c]))
    target_atoms = [a for a in atoms if a[0] == target_chain]
    binder_atoms = [a for a in atoms if a[0] == binder_chain]

    clashes = 0
    target_iface = set()
    binder_iface = set()
    for ba in binder_atoms:
        for ta in target_atoms:
            d = distance(ba, ta)
            if d < 2.0:
                clashes += 1
            if d < 5.0:
                target_iface.add(ta[1])
                binder_iface.add(ba[1])

    hit = 0
    hotspot_distances = []
    missing = []
    for token in FAMILIES[family]:
        source_resi = int(token[1:])
        hs_atoms = [a for a in target_atoms if a[1] == source_resi]
        if not hs_atoms:
            missing.append(token)
            hotspot_distances.append(f"{token}:missing")
            continue
        mind = min(distance(h, b) for h in hs_atoms for b in binder_atoms)
        hotspot_distances.append(f"{token}:{mind:.2f}")
        if mind <= 10.0:
            hit += 1

    passed = (
        hit >= 3
        and len(target_iface) >= 5
        and len(binder_iface) >= 5
        and clashes <= 10
    )
    note = ""
    if missing:
        note = "Missing output hotspot residue IDs: " + ",".join(missing)

    return {
        "family": family,
        "pdb": str(pdb_path.relative_to(WORK)),
        "target_chain": target_chain,
        "binder_chain": binder_chain,
        "target_residue_count": len(residue_sets[target_chain]),
        "binder_residue_count": len(residue_sets[binder_chain]),
        "severe_clashes_lt_2A": clashes,
        "interface_target_residues_lt_5A": len(target_iface),
        "interface_binder_residues_lt_5A": len(binder_iface),
        "hotspots_within_10A": hit,
        "hotspots_total": len(FAMILIES[family]),
        "hotspot_distances_A": ";".join(hotspot_distances),
        "geometry_pass": passed,
        "screen_note": note,
    }


def screen_outputs():
    rows = []
    for family in FAMILIES:
        for pdb in sorted((OUT / family).glob("design_*.pdb")):
            rows.append(screen_design(pdb, family))

    csv_path = OUT / "pilot_geometry_screen.csv"
    if rows:
        fields = sorted({k for row in rows for k in row})
        with csv_path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    summary = {
        "generated_pdbs": len(rows),
        "geometry_pass": sum(bool(r.get("geometry_pass")) for r in rows),
        "families": FAMILIES,
        "generation_target": "11QC chain A residues 189-383",
        "full_target_screening_required": True,
        "screening_is_preliminary": True,
        "next_stage": "ProteinMPNN, then full human/mouse TfR1 Boltz-2 screening",
    }
    (OUT / "pilot_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


def package_outputs():
    archive = shutil.make_archive(
        str(WORK / "tfr1_rfdiffusion_pilot_v1"), "zip", root_dir=OUT
    )
    print("Archive:", archive)


if __name__ == "__main__":
    t0 = time.time()
    setup_environment()
    download_and_crop_target()
    run_diffusion()
    screen_outputs()
    package_outputs()
    print(f"Total elapsed: {(time.time()-t0)/60:.1f} min")
