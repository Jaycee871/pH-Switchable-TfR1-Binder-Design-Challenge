#!/usr/bin/env python3
"""Phase 2G: local-epitope RFdiffusion around the conserved TfR1 343-350 cluster.

Why this exists
---------------
Phase 2D-2F showed that even an apical-domain-only target (189-383) leaves enough
surface for RFdiffusion to escape to residues ~210 or ~291. Compact backbones are
now achievable, but dual-anchor engagement is inconsistent.

This phase crops the fixed target further to human TfR1 residues 330-383. The
343/344/347/350 cluster retains its native 11QC geometry and flanking context,
while the major alternative surfaces observed in prior campaigns are removed.

Gate:
- binder 45-60 aa
- >=2 preferred conserved contacts
- zero risk contact to N348/K350-equivalent risk position (human 348) or K371/Q373
- <=1 contact in artificial crop-edge zones 330-334 / 379-383
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

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "structures/11QC.pdb"
OUT = ROOT / "structures/rfdiffusion/nvidia_local_343_window"
URL = "https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate"

WINDOW_START, WINDOW_END = 330, 383
CONTIGS = "A330-383/0 45-60"

FAMILIES = {
    "core4": ["A343", "A344", "A347", "A350"],
    "core3": ["A343", "A347", "A350"],
    "inner_pair": ["A344", "A347"],
}
SEEDS = {
    "core4": [20261401, 20261402],
    "core3": [20261403, 20261404],
    "inner_pair": [20261405, 20261406],
}

PREFERRED = {343, 344, 347, 350, 374}
RISK = {348, 371}
EDGE = set(range(330, 335)) | set(range(379, 384))


def crop(text: str) -> str:
    kept, residues = [], set()
    for line in text.splitlines():
        if not line.startswith("ATOM  ") or len(line) < 54 or line[21] != "A":
            continue
        try:
            r = int(line[22:26])
        except ValueError:
            continue
        if WINDOW_START <= r <= WINDOW_END:
            kept.append(line)
            residues.add(r)
    need = {343, 344, 347, 348, 350, 371, 374}
    missing = need - residues
    if missing:
        raise RuntimeError(f"Local 11QC window missing residues {sorted(missing)}")
    return "\n".join(kept) + "\nTER\nEND\n"


def post(payload: dict, key: str, max_attempts: int = 5) -> dict:
    for attempt in range(1, max_attempts + 1):
        req = Request(
            URL,
            data=json.dumps(payload).encode("utf-8"),
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {key}",
            },
            method="POST",
        )
        try:
            with urlopen(req, timeout=1200) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")
            if exc.code == 429 and attempt < max_attempts:
                ra = exc.headers.get("Retry-After")
                sleep_s = int(ra) if ra and ra.isdigit() else 25 * attempt
                print(f"429 attempt {attempt}; sleep {sleep_s}s", flush=True)
                time.sleep(sleep_s)
                continue
            transient = exc.code in {500, 502, 503, 504} or (
                exc.code == 422
                and any(x in detail for x in ("InfrastructureError", "AcceleratorError", "CUDA error"))
            )
            if transient and attempt < max_attempts:
                sleep_s = 25 * attempt
                print(f"Transient {exc.code}; sleep {sleep_s}s", flush=True)
                time.sleep(sleep_s)
                continue
            raise RuntimeError(f"RFdiffusion HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("RFdiffusion retries exhausted")


def parse(text: str):
    atoms, residues, ca = {}, {}, {}
    for line in text.splitlines():
        if not line.startswith("ATOM  ") or len(line) < 54:
            continue
        atom = line[12:16].strip()
        if atom.startswith("H"):
            continue
        chain = line[21].strip() or "_"
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
    ctr = tuple(sum(p[i] for p in points)/n for i in range(3))
    rg = math.sqrt(sum(dist(p, ctr)**2 for p in points)/n)
    span = max(
        (dist(a, b) for i, a in enumerate(points) for b in points[i+1:]),
        default=0.0,
    )
    return round(rg, 3), round(span, 3)


def analyze(text: str):
    atoms, residues, ca = parse(text)
    if len(residues) < 2:
        raise RuntimeError("Expected target+binder chains")
    binder = min(residues, key=lambda c: len(residues[c]))
    target = max(residues, key=lambda c: len(residues[c]))
    target_first = min(residues[target])

    contacts, clashes = set(), 0
    for _, _, bxyz in atoms[binder]:
        for out_res, _, txyz in atoms[target]:
            d = dist(bxyz, txyz)
            if d < 2.0:
                clashes += 1
            if d <= 5.0:
                contacts.add(out_res - target_first + WINDOW_START)

    preferred = sorted(contacts & PREFERRED)
    risk = sorted(contacts & RISK)
    edge = sorted(contacts & EDGE)
    rg, span = compactness([ca[binder][r] for r in sorted(ca[binder])])

    checks = {
        "binder_45_60": 45 <= len(residues[binder]) <= 60,
        "preferred_ge_2": len(preferred) >= 2,
        "risk_eq_0": len(risk) == 0,
        "edge_contacts_le_1": len(edge) <= 1,
        "rg_le_16_5": rg <= 16.5,
        "span_le_52": span <= 52.0,
        "clashes_le_10": clashes <= 10,
    }
    return {
        "binder_chain": binder,
        "target_chain": target,
        "binder_length": len(residues[binder]),
        "target_length": len(residues[target]),
        "contacts": sorted(contacts),
        "preferred_contacts": preferred,
        "risk_contacts": risk,
        "edge_contacts": edge,
        "radius_of_gyration_A": rg,
        "max_CA_span_A": span,
        "severe_clashes_lt_2A": clashes,
        "checks": checks,
        "scientific_pass": all(checks.values()),
    }


def main():
    key = (
        os.getenv("NVIDIA_API_KEY")
        or os.getenv("NVIDIA_API_KEY1")
        or os.getenv("NVIDIA_API_TOKEN")
        or os.getenv("NGC_API_KEY")
    )
    if not key:
        raise RuntimeError("Missing NVIDIA API credential")

    target = crop(INPUT.read_text(encoding="utf-8"))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "target_11QC_330_383.pdb").write_text(target, encoding="utf-8")

    attempts = []
    for family, hotspots in FAMILIES.items():
        fdir = OUT / family
        fdir.mkdir(parents=True, exist_ok=True)
        for seed in SEEDS[family]:
            result = post({
                "input_pdb": target,
                "contigs": CONTIGS,
                "hotspot_res": hotspots,
                "diffusion_steps": 50,
                "random_seed": seed,
            }, key)
            pdb = result.get("output_pdb")
            if not pdb:
                raise RuntimeError(f"No output_pdb for {family}/{seed}")
            path = fdir / f"{family}_seed_{seed}.pdb"
            path.write_text(pdb, encoding="utf-8")
            qc = analyze(pdb)
            rec = {
                "family": family,
                "hotspots": hotspots,
                "seed": seed,
                "elapsed_ms": result.get("elapsed_ms"),
                "pdb": str(path.relative_to(ROOT)),
                "qc": qc,
            }
            attempts.append(rec)
            print(json.dumps(rec, indent=2), flush=True)
            time.sleep(15)

    passing = [x for x in attempts if x["qc"]["scientific_pass"]]
    passing.sort(
        key=lambda x: (
            len(x["qc"]["preferred_contacts"]),
            -len(x["qc"]["edge_contacts"]),
            -x["qc"]["radius_of_gyration_A"],
            -x["qc"]["max_CA_span_A"],
        ),
        reverse=True,
    )

    selected_dir = OUT / "selected"
    if selected_dir.exists():
        shutil.rmtree(selected_dir)
    selected_dir.mkdir(parents=True, exist_ok=True)

    for i, x in enumerate(passing[:4], start=1):
        dst = selected_dir / f"selected_{i:02d}_{x['family']}_seed_{x['seed']}.pdb"
        shutil.copy2(ROOT / x["pdb"], dst)
        x["selected_copy"] = str(dst.relative_to(ROOT))

    summary = {
        "campaign": "local_343_window_v1",
        "fixed_target_window": [WINDOW_START, WINDOW_END],
        "contigs": CONTIGS,
        "families": FAMILIES,
        "seeds": SEEDS,
        "gate": {
            "preferred_ge_2": True,
            "risk_eq_0": True,
            "edge_contacts_le_1": True,
            "rg_le_16_5": True,
            "span_le_52": True,
        },
        "attempt_count": len(attempts),
        "passing_count": len(passing),
        "selected_count": min(4, len(passing)),
        "selected": passing[:4],
        "attempts": attempts,
    }
    (OUT / "selection_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "passing_count": len(passing),
        "selected": passing[:4],
    }, indent=2), flush=True)

    if not passing:
        raise SystemExit("No local-window backbone passed the strict gate")


if __name__ == "__main__":
    main()
