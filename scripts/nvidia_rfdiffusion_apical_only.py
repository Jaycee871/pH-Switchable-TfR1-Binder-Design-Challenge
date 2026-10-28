#!/usr/bin/env python3
"""Phase 2D: NVIDIA RFdiffusion on an apical-domain-only TfR1 target.

The previous compact campaign supplied nearly the full TfR1 ectodomain, and every
sequence-tested scaffold later collapsed to distal surfaces in full-complex prediction.
This campaign removes those alternative surfaces during backbone generation by cropping
11QC chain A to residues 189-383 (the apical domain).

Three hotspot families are sampled to create structural diversity:
- conserved_loop: 206/209/214/231
- conserved_343: 343/344/347/350
- bridge: 209/214/343/347

All outputs are retained and independently gated. The script does not stop after the
first pass because diversity is an explicit design objective.
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

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "structures/11QC.pdb"
OUT = ROOT / "structures/rfdiffusion/nvidia_apical_only"
URL = "https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate"

DOMAIN_START = 189
DOMAIN_END = 383
CONTIGS = "A189-383/0 60-90"

FAMILIES = {
    "conserved_loop": ["A206", "A209", "A214", "A231"],
    "conserved_343": ["A343", "A344", "A347", "A350"],
    "bridge": ["A209", "A214", "A343", "A347"],
}

SEEDS = {
    "conserved_loop": [20261121, 20261122],
    "conserved_343": [20261123, 20261124],
    "bridge": [20261125, 20261126],
}

PREFERRED = {206, 207, 209, 214, 231, 343, 344, 347, 350, 374}
RISK = {204, 208, 210, 211, 212, 348, 371}


def crop_chain_a_apical(text: str) -> str:
    kept = []
    residues = set()
    for line in text.splitlines():
        if not line.startswith("ATOM  "):
            continue
        if len(line) < 54 or line[21] != "A":
            continue
        try:
            r = int(line[22:26])
        except ValueError:
            continue
        if DOMAIN_START <= r <= DOMAIN_END:
            kept.append(line)
            residues.add(r)
    if not residues:
        raise RuntimeError("Failed to crop 11QC chain A apical domain")
    required = {206, 209, 214, 231, 343, 344, 347, 350}
    missing = sorted(required - residues)
    if missing:
        raise RuntimeError(f"Apical crop missing required hotspot residues: {missing}")
    return "\n".join(kept) + "\nTER\nEND\n"


def post_json(payload: dict, api_key: str, max_attempts: int = 5) -> dict:
    for attempt in range(1, max_attempts + 1):
        req = Request(
            URL,
            data=json.dumps(payload).encode("utf-8"),
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {api_key}",
            },
            method="POST",
        )
        try:
            with urlopen(req, timeout=1200) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")
            if exc.code == 429 and attempt < max_attempts:
                retry_after = exc.headers.get("Retry-After")
                sleep_s = int(retry_after) if retry_after and retry_after.isdigit() else 25 * attempt
                print(f"RFdiffusion 429 attempt {attempt}; sleeping {sleep_s}s", flush=True)
                time.sleep(sleep_s)
                continue
            transient = (
                exc.code in {500, 502, 503, 504}
                or (
                    exc.code == 422
                    and (
                        "InfrastructureError" in detail
                        or "AcceleratorError" in detail
                        or "CUDA error" in detail
                    )
                )
            )
            if transient and attempt < max_attempts:
                sleep_s = 25 * attempt
                print(
                    f"RFdiffusion transient HTTP {exc.code} attempt {attempt}; "
                    f"sleeping {sleep_s}s",
                    flush=True,
                )
                time.sleep(sleep_s)
                continue
            raise RuntimeError(f"RFdiffusion HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("RFdiffusion retry loop exhausted")


def parse_structure(pdb_text: str):
    atoms = {}
    residues = {}
    ca = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM  ") or len(line) < 54:
            continue
        chain = line[21].strip() or "_"
        atom = line[12:16].strip()
        if atom.startswith("H"):
            continue
        try:
            resid = int(line[22:26])
            xyz = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
        except ValueError:
            continue
        atoms.setdefault(chain, []).append((resid, atom, xyz))
        residues.setdefault(chain, set()).add(resid)
        if atom == "CA":
            ca.setdefault(chain, {})[resid] = xyz
    return atoms, residues, ca


def dist(a, b):
    return math.sqrt(sum((x-y)**2 for x, y in zip(a, b)))


def compactness(points):
    n = len(points)
    if not n:
        return {"radius_of_gyration_A": None, "max_CA_span_A": None}
    ctr = tuple(sum(p[i] for p in points)/n for i in range(3))
    rg = math.sqrt(sum(dist(p, ctr)**2 for p in points)/n)
    span = 0.0
    for i, a in enumerate(points):
        for b in points[i+1:]:
            span = max(span, dist(a, b))
    return {
        "radius_of_gyration_A": round(rg, 3),
        "max_CA_span_A": round(span, 3),
    }


def analyze(pdb_text: str):
    atoms, residues, ca = parse_structure(pdb_text)
    if len(residues) < 2:
        raise RuntimeError(f"Expected target+binder chains, got {residues}")

    binder = min(residues, key=lambda c: len(residues[c]))
    target = max(residues, key=lambda c: len(residues[c]))
    target_first = min(residues[target])

    ba = [x for x in atoms[binder]]
    ta = [x for x in atoms[target]]
    contacts = set()
    clashes = 0
    for _, _, bxyz in ba:
        for out_res, _, txyz in ta:
            d = dist(bxyz, txyz)
            if d < 2.0:
                clashes += 1
            if d <= 5.0:
                contacts.add(out_res - target_first + DOMAIN_START)

    preferred = sorted(contacts & PREFERRED)
    risk = sorted(contacts & RISK)
    geom = compactness([ca[binder][r] for r in sorted(ca.get(binder, {}))])

    checks = {
        "binder_length_60_90": 60 <= len(residues[binder]) <= 90,
        "preferred_contacts_ge_2": len(preferred) >= 2,
        "risk_contacts_le_1": len(risk) <= 1,
        "rg_le_18A": geom["radius_of_gyration_A"] is not None and geom["radius_of_gyration_A"] <= 18.0,
        "ca_span_le_60A": geom["max_CA_span_A"] is not None and geom["max_CA_span_A"] <= 60.0,
        "clashes_le_10": clashes <= 10,
    }

    return {
        "binder_chain": binder,
        "target_chain": target,
        "binder_length": len(residues[binder]),
        "target_length": len(residues[target]),
        "human_target_contacts_5A": sorted(contacts),
        "preferred_contacts": preferred,
        "risk_contacts": risk,
        "severe_clashes_lt_2A": clashes,
        "compactness": geom,
        "checks": checks,
        "scientific_pass": all(checks.values()),
    }


def main():
    api_key = (
        os.getenv("NVIDIA_API_KEY")
        or os.getenv("NVIDIA_API_KEY1")
        or os.getenv("NVIDIA_API_TOKEN")
        or os.getenv("NGC_API_KEY")
    )
    if not api_key:
        raise RuntimeError("Missing NVIDIA API credential")
    if not INPUT.exists():
        raise RuntimeError("Missing structures/11QC.pdb")

    target = crop_chain_a_apical(INPUT.read_text(encoding="utf-8"))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "target_11QC_chainA_apical_189_383.pdb").write_text(target, encoding="utf-8")

    attempts = []
    for family, hotspots in FAMILIES.items():
        fdir = OUT / family
        fdir.mkdir(parents=True, exist_ok=True)
        for seed in SEEDS[family]:
            payload = {
                "input_pdb": target,
                "contigs": CONTIGS,
                "hotspot_res": hotspots,
                "diffusion_steps": 50,
                "random_seed": seed,
            }
            result = post_json(payload, api_key)
            output = result.get("output_pdb")
            if not output:
                raise RuntimeError(f"No output_pdb for {family} seed {seed}")
            path = fdir / f"{family}_seed_{seed}.pdb"
            path.write_text(output, encoding="utf-8")
            qc = analyze(output)
            attempts.append({
                "family": family,
                "hotspots": hotspots,
                "seed": seed,
                "elapsed_ms": result.get("elapsed_ms"),
                "pdb": str(path.relative_to(ROOT)),
                "qc": qc,
            })
            print(json.dumps(attempts[-1], indent=2), flush=True)
            time.sleep(15)

    passing = [x for x in attempts if x["qc"]["scientific_pass"]]
    passing.sort(
        key=lambda x: (
            len(x["qc"]["preferred_contacts"]),
            -len(x["qc"]["risk_contacts"]),
            -x["qc"]["compactness"]["radius_of_gyration_A"],
        ),
        reverse=True,
    )

    selected_dir = OUT / "selected"
    if selected_dir.exists():
        shutil.rmtree(selected_dir)
    selected_dir.mkdir(parents=True, exist_ok=True)
    for i, x in enumerate(passing[:4], start=1):
        src = ROOT / x["pdb"]
        dst = selected_dir / f"selected_{i:02d}_{x['family']}_seed_{x['seed']}.pdb"
        shutil.copy2(src, dst)
        x["selected_copy"] = str(dst.relative_to(ROOT))

    summary = {
        "campaign": "nvidia_apical_only_v1",
        "generation_target": "11QC chain A residues 189-383 only",
        "contigs": CONTIGS,
        "families": FAMILIES,
        "seeds": SEEDS,
        "preferred_residues": sorted(PREFERRED),
        "risk_residues": sorted(RISK),
        "attempt_count": len(attempts),
        "passing_count": len(passing),
        "selected_count": min(4, len(passing)),
        "selected": passing[:4],
        "attempts": attempts,
    }
    (OUT / "selection_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print(json.dumps({
        "attempt_count": len(attempts),
        "passing_count": len(passing),
        "selected": [
            {
                "family": x["family"],
                "seed": x["seed"],
                "preferred": x["qc"]["preferred_contacts"],
                "risk": x["qc"]["risk_contacts"],
                "rg": x["qc"]["compactness"]["radius_of_gyration_A"],
                "span": x["qc"]["compactness"]["max_CA_span_A"],
                "copy": x.get("selected_copy"),
            }
            for x in passing[:4]
        ],
    }, indent=2), flush=True)

    if not passing:
        raise SystemExit("No NVIDIA apical-only RFdiffusion backbone passed the strict gate.")


if __name__ == "__main__":
    main()
