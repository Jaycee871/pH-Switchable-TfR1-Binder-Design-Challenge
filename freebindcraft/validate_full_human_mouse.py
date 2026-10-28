#!/usr/bin/env python3
"""Strict full human/mouse pose validation for FreeBindCraft TfR1 candidates.

FreeBindCraft acceptance is a generation-stage filter only. This script takes
accepted designs from a synced FreeBindCraft run, rebuilds human and mouse
full-ectodomain structural templates from each accepted apical complex, and
runs the same strict Boltz-2 cross-species pose gate used elsewhere in this repo.

No pH-switch claim is made here. pH engineering starts only after a design
passes this cross-species pose gate.
"""

from __future__ import annotations

import csv
import json
import os
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFIO, MMCIFParser, NeighborSearch, PDBParser

ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT / "freebindcraft/kaggle_smoke_results/designs"
FINAL_CSV = RUN_DIR / "final_design_stats.csv"
HUMAN_FASTA = ROOT / "data/human_TfR1/human_TfR1_P02786_121-760.fasta"
MOUSE_FASTA = ROOT / "data/mouse_TfR1/mouse_TfR1_Q62351_123-763.fasta"
MASK_CSV = ROOT / "analysis/apical_design_mask.csv"
OUT = ROOT / "predictions/freebindcraft_full_pose_screen"
URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"

HUMAN_START = 121
MOUSE_START = 123
DOMAIN_START = 189
MAX_DESIGNS = 2

PREFERRED_H = {206, 207, 209, 214, 231, 343, 344, 347, 350, 374}
RISK_H = {204, 208, 210, 211, 212, 348, 371}

OFFICIAL_OVERRIDES = {
    204: 207,
    208: 211,
    210: 213,
    211: 214,
    212: None,
    348: 350,
    371: 373,
}

AA3 = {
    "A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","E":"GLU","Q":"GLN",
    "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
    "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL",
}


def fasta(path: Path) -> str:
    return "".join(
        x.strip() for x in path.read_text(encoding="utf-8").splitlines()
        if x.strip() and not x.startswith(">")
    )


def human_mouse_map():
    hs, ms = fasta(HUMAN_FASTA), fasta(MOUSE_FASTA)
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    aln = a.align(hs, ms)[0]
    c = aln.coordinates

    mapping = {}
    for k in range(c.shape[1]-1):
        h0, h1 = int(c[0,k]), int(c[0,k+1])
        m0, m1 = int(c[1,k]), int(c[1,k+1])
        dh, dm = h1-h0, m1-m0
        if dh and dm:
            for j in range(min(dh, dm)):
                mapping[HUMAN_START+h0+j] = MOUSE_START+m0+j
        elif dh:
            for j in range(dh):
                mapping[HUMAN_START+h0+j] = None

    mapping.update(OFFICIAL_OVERRIDES)

    if MASK_CSV.exists():
        with MASK_CSV.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                try:
                    hp = int(row["human_residue"])
                except (ValueError, TypeError):
                    continue
                mp = (row.get("mouse_residue") or "").strip()
                mapping[hp] = int(mp) if mp else None
    return mapping


def rewrite(line: str, chain: str, resid: int, resname: str | None = None) -> str:
    chars = list(line.rstrip("\n"))
    if len(chars) < 80:
        chars.extend(" " * (80-len(chars)))
    chars[21] = chain
    chars[22:26] = list(f"{resid:4d}")
    if resname is not None:
        chars[17:20] = list(f"{resname:>3}")
    return "".join(chars).rstrip() + "\n"


def chain_residues(path: Path, chain: str):
    vals, seen = [], set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM") or line[21].strip() != chain:
            continue
        try:
            r = int(line[22:26])
        except ValueError:
            continue
        if r not in seen:
            vals.append(r)
            seen.add(r)
    return vals


def locate_accepted_pdb(design: str) -> Path:
    accepted = RUN_DIR / "Accepted"
    candidates = [accepted / f"{design}.pdb"]
    if accepted.exists():
        candidates.extend(accepted.glob(f"**/{design}.pdb"))
    candidates.extend(RUN_DIR.glob(f"**/{design}.pdb"))
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(f"Accepted PDB not found for {design}")


def infer_target_binder(path: Path, binder_length: int):
    chains = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        ch = line[21].strip() or "_"
        try:
            r = int(line[22:26])
        except ValueError:
            continue
        chains.setdefault(ch, set()).add(r)

    if len(chains) < 2:
        raise RuntimeError(f"Expected target+binder chains in {path}: {chains}")

    # Target settings explicitly use chain A, but retain a length-based fallback
    # for resilience across FreeBindCraft output variants.
    target = "A" if "A" in chains and len(chains["A"]) >= 150 else max(
        chains, key=lambda c: len(chains[c])
    )
    others = [c for c in chains if c != target]
    binder = min(others, key=lambda c: abs(len(chains[c]) - binder_length))
    if abs(len(chains[binder]) - binder_length) > 2:
        raise RuntimeError(
            f"Binder-chain inference mismatch for {path}: "
            f"sequence={binder_length}, chains={ {k:len(v) for k,v in chains.items()} }"
        )
    return target, binder


def make_template(
    source_pdb: Path,
    source_target: str,
    source_binder: str,
    binder_seq: str,
    species: str,
    hm_map: dict,
    prefix: Path,
):
    target_res = chain_residues(source_pdb, source_target)
    binder_res = chain_residues(source_pdb, source_binder)
    if not target_res or not binder_res:
        raise RuntimeError(f"Missing chains in {source_pdb}")

    # FreeBindCraft starts from a cropped 189-383 target that normally preserves
    # original residue numbering. If it was renumbered, remap by domain offset.
    preserve_target_numbering = min(target_res) >= 180 and max(target_res) <= 400
    binder_index = {r:i+1 for i,r in enumerate(binder_res)}

    target_lines, binder_lines = [], []
    for line in source_pdb.read_text(encoding="utf-8", errors="replace").splitlines(True):
        if not line.startswith("ATOM"):
            continue
        ch = line[21].strip()
        try:
            old = int(line[22:26])
        except ValueError:
            continue

        if ch == source_target:
            human_uniprot = (
                old if preserve_target_numbering
                else old - min(target_res) + DOMAIN_START
            )
            if not 189 <= human_uniprot <= 383:
                continue

            if species == "human":
                query_pos = human_uniprot - HUMAN_START + 1
            else:
                mouse_uniprot = hm_map.get(human_uniprot)
                if mouse_uniprot is None:
                    continue
                query_pos = mouse_uniprot - MOUSE_START + 1
            if query_pos < 1:
                continue
            target_lines.append(rewrite(line, "A", query_pos))

        elif ch == source_binder:
            pos = binder_index.get(old)
            if pos is None or pos > len(binder_seq):
                continue
            binder_lines.append(rewrite(line, "B", pos, AA3[binder_seq[pos-1]]))

    if not target_lines or not binder_lines:
        raise RuntimeError("Template construction produced an empty chain")

    prefix.parent.mkdir(parents=True, exist_ok=True)
    pdb_path = prefix.with_suffix(".pdb")
    cif_path = prefix.with_suffix(".cif")
    pdb_path.write_text(
        "".join(target_lines) + "TER\n" + "".join(binder_lines) + "TER\nEND\n",
        encoding="utf-8",
    )
    structure = PDBParser(QUIET=True).get_structure(prefix.stem, str(pdb_path))
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(cif_path))
    return cif_path


def post_json(payload: dict, api_key: str, max_attempts: int = 5):
    for attempt in range(1, max_attempts+1):
        req = Request(
            URL,
            data=json.dumps(payload).encode("utf-8"),
            headers={
                "Content-Type":"application/json",
                "Authorization":f"Bearer {api_key}",
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
                sleep_s = int(ra) if ra and ra.isdigit() else 25*attempt
                print(f"Boltz2 429 attempt {attempt}; sleep {sleep_s}s", flush=True)
                time.sleep(sleep_s)
                continue
            transient = exc.code in {500,502,503,504} or (
                exc.code == 422
                and any(x in detail for x in (
                    "InfrastructureError","AcceleratorError","CUDA error"
                ))
            )
            if transient and attempt < max_attempts:
                sleep_s = 25*attempt
                print(f"Boltz2 transient {exc.code}; sleep {sleep_s}s", flush=True)
                time.sleep(sleep_s)
                continue
            raise RuntimeError(f"Boltz2 HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("Boltz2 retries exhausted")


def contacts(cif_path: Path, offset: int, preferred: set, risk: set):
    structure = MMCIFParser(QUIET=True).get_structure("pred", str(cif_path))
    model = next(structure.get_models())
    target, binder = model["A"], model["B"]
    target_atoms = [
        a for r in target if r.id[0] == " "
        for a in r if a.element != "H"
    ]
    search = NeighborSearch(target_atoms)
    tc, bc = set(), set()

    for r in binder:
        if r.id[0] != " ":
            continue
        for a in r:
            if a.element == "H":
                continue
            close = search.search(a.coord, 5.0, level="A")
            if close:
                bc.add(r.id[1])
            for x in close:
                tc.add(x.get_parent().id[1] + offset)

    tc = sorted(tc)
    return {
        "target_contacts_uniprot": tc,
        "binder_contact_positions": sorted(bc),
        "preferred_contacts": sorted(set(tc) & preferred),
        "risk_contacts": sorted(set(tc) & risk),
        "preferred_core_retained": bool(set(tc) & preferred),
    }


def fnum(row: dict, key: str, default: float):
    try:
        return float((row.get(key) or "").strip())
    except (ValueError, TypeError):
        return default


def select_candidates():
    if not FINAL_CSV.exists():
        raise SystemExit(f"FreeBindCraft accepted-design table missing: {FINAL_CSV}")
    with FINAL_CSV.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit("FreeBindCraft produced no accepted designs")

    # Prefer ipSAE when available, then iPTM. The competition-specific structural
    # gate below still decides whether these candidates survive.
    rows.sort(
        key=lambda r: (
            -fnum(r, "Average_ipSAE", -999.0),
            -fnum(r, "Average_i_pTM", -999.0),
            fnum(r, "Average_i_pAE", 999.0),
        )
    )
    return rows[:MAX_DESIGNS]


def main():
    api_key = (
        os.getenv("NVIDIA_API_KEY")
        or os.getenv("NVIDIA_API_KEY1")
        or os.getenv("NVIDIA_API_TOKEN")
        or os.getenv("NGC_API_KEY")
    )
    if not api_key:
        raise RuntimeError("Missing NVIDIA API credential")

    rows = select_candidates()
    hm = human_mouse_map()
    mh = {m:h for h,m in hm.items() if m is not None}
    preferred_m = {hm[h] for h in PREFERRED_H if hm.get(h) is not None}
    risk_m = {hm[h] for h in RISK_H if hm.get(h) is not None}

    targets = {
        "human": {
            "sequence":fasta(HUMAN_FASTA),
            "offset":HUMAN_START-1,
            "preferred":PREFERRED_H,
            "risk":RISK_H,
        },
        "mouse": {
            "sequence":fasta(MOUSE_FASTA),
            "offset":MOUSE_START-1,
            "preferred":preferred_m,
            "risk":risk_m,
        },
    }

    OUT.mkdir(parents=True, exist_ok=True)
    batch = {}

    for row in rows:
        cid = row["Design"].strip()
        binder = row["Sequence"].strip()
        source = locate_accepted_pdb(cid)
        starget, sbinder = infer_target_binder(source, len(binder))
        cdir = OUT / cid

        cresult = {
            "candidate_id":cid,
            "binder_sequence":binder,
            "source_freebindcraft_pdb":str(source.relative_to(ROOT)),
            "freebindcraft_metrics":{
                "rank":row.get("Rank"),
                "average_i_pTM":row.get("Average_i_pTM"),
                "average_ipSAE":row.get("Average_ipSAE"),
                "average_i_pAE":row.get("Average_i_pAE"),
                "average_pLDDT":row.get("Average_pLDDT"),
                "hotspot_RMSD_A":row.get("Average_Hotspot_RMSD"),
            },
            "species":{},
        }

        for species, spec in targets.items():
            template = make_template(
                source, starget, sbinder, binder, species, hm,
                cdir / f"{cid}_{species}_apical_template",
            )
            template_text = template.read_text(encoding="utf-8")
            payload = {
                "polymers":[
                    {
                        "id":"A","molecule_type":"protein","sequence":spec["sequence"],
                        "structural_templates":[{
                            "structure":template_text,
                            "format":"cif",
                            "name":f"{cid}_{species}_target",
                            "chain_id":"A",
                        }],
                    },
                    {
                        "id":"B","molecule_type":"protein","sequence":binder,
                        "structural_templates":[{
                            "structure":template_text,
                            "format":"cif",
                            "name":f"{cid}_{species}_binder",
                            "chain_id":"B",
                        }],
                    },
                ],
                "recycling_steps":3,
                "sampling_steps":50,
                "diffusion_samples":1,
                "step_scale":1.638,
                "output_format":"mmcif",
                "write_full_pae":False,
            }
            result = post_json(payload, api_key)
            structs = result.get("structures", [])
            if not structs:
                raise RuntimeError(f"No Boltz-2 structure for {cid}/{species}")

            spdir = cdir / species
            spdir.mkdir(parents=True, exist_ok=True)
            cif = spdir / "structure_1.cif"
            cif.write_text(structs[0]["structure"], encoding="utf-8")

            cm = contacts(cif, spec["offset"], spec["preferred"], spec["risk"])
            iptm = result.get("protein_iptm_scores", [None])[0]
            apical = any(189 <= r <= 383 for r in cm["target_contacts_uniprot"])
            pose_pass = bool(
                iptm is not None
                and float(iptm) >= 0.50
                and cm["preferred_core_retained"]
                and apical
            )
            summary = {
                "confidence":result.get("confidence_scores", [None])[0],
                "ptm":result.get("ptm_scores", [None])[0],
                "protein_iptm":iptm,
                "complex_plddt":result.get("complex_plddt_scores", [None])[0],
                "complex_iplddt":result.get("complex_iplddt_scores", [None])[0],
                "complex_pde":result.get("complex_pde_scores", [None])[0],
                "complex_ipde":result.get("complex_ipde_scores", [None])[0],
                "contacts":cm,
                "apical_contact":apical,
                "pose_gate":{
                    "min_protein_iptm":0.50,
                    "requires_preferred_contact":True,
                    "pass":pose_pass,
                },
            }
            (spdir / "result_summary.json").write_text(
                json.dumps(summary, indent=2) + "\n", encoding="utf-8"
            )
            cresult["species"][species] = summary
            print(cid, species, json.dumps(summary, indent=2), flush=True)
            time.sleep(20)

        human = cresult["species"]["human"]
        mouse = cresult["species"]["mouse"]
        h_contacts = set(human["contacts"]["target_contacts_uniprot"])
        m_as_h = {
            mh[x] for x in mouse["contacts"]["target_contacts_uniprot"] if x in mh
        }
        overlap = sorted(h_contacts & m_as_h)
        union = h_contacts | m_as_h
        jaccard = len(overlap)/len(union) if union else 0.0
        same_epitope = len(overlap) >= 3 and jaccard >= 0.15

        cresult["cross_species_epitope"] = {
            "human_mouse_homologous_overlap_human_numbering":overlap,
            "overlap_count":len(overlap),
            "jaccard":round(jaccard, 4),
            "same_epitope_pass":same_epitope,
        }
        cresult["cross_species_pose_pass"] = bool(
            human["pose_gate"]["pass"]
            and mouse["pose_gate"]["pass"]
            and same_epitope
        )
        batch[cid] = cresult

    (OUT / "batch_summary.json").write_text(
        json.dumps(batch, indent=2) + "\n", encoding="utf-8"
    )
    winners = [k for k,v in batch.items() if v["cross_species_pose_pass"]]
    (OUT / "pose_pass_ids.json").write_text(
        json.dumps({"pose_pass_ids":winners}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "validated":len(batch),
        "cross_species_pose_pass":winners,
    }, indent=2))


if __name__ == "__main__":
    main()
