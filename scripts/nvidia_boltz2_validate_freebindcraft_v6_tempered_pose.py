#!/usr/bin/env python3
"""Phase 4M: full human/mouse pose validation for the FreeBindCraft v6 tempered-binary trajectory.

Identical gate to Phase 4F/4G so the result is directly comparable to the seven
earlier validated candidates. Nothing here is loosened:

- per species: protein ipTM >= 0.50 AND a preferred apical contact AND a contact
  inside the 189-383 apical domain;
- cross species: >= 3 homologous overlapping contacts AND Jaccard >= 0.15;
- cross_species_pose_pass requires all of the above.

The only preprocessing difference versus 4F/4G is that the FreeBindCraft source
PDB is an OpenMM-relaxed model carrying explicit hydrogens, while the RFdiffusion
sources were heavy-atom only. Hydrogens are stripped during template construction
so the Boltz-2 template is built the same way as in the earlier phases.

Provenance caveat recorded deliberately: this candidate is a FreeBindCraft
trajectory output (Trajectory/Relaxed). It did NOT clear FreeBindCraft own
final design filters - final_design_stats.csv is header-only - so its reported
i_pTM 0.69 is the AF2 hallucination loop scoring its own optimisation target and
must not be treated as independent evidence. This script is that independent check.
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
CANDIDATES = ROOT / "designs/family_D_freebindcraft/v6_tempered_binary_pose_screen/pose_screen_candidates.csv"
HUMAN_FASTA = ROOT / "data/human_TfR1/human_TfR1_P02786_121-760.fasta"
MOUSE_FASTA = ROOT / "data/mouse_TfR1/mouse_TfR1_Q62351_123-763.fasta"
MASK_CSV = ROOT / "analysis/apical_design_mask.csv"
OUT = ROOT / "predictions/freebindcraft_v6_tempered_binary_pose_screen"
URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"

HUMAN_START = 121
MOUSE_START = 123
DOMAIN_START = 189

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


def fasta(path):
    return "".join(
        x.strip() for x in path.read_text(encoding="utf-8").splitlines()
        if x.strip() and not x.startswith(">")
    )


def human_mouse_map():
    hs = fasta(HUMAN_FASTA)
    ms = fasta(MOUSE_FASTA)
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    aln = a.align(hs, ms)[0]
    c = aln.coordinates

    mapping = {}
    for k in range(c.shape[1] - 1):
        h0, h1 = int(c[0,k]), int(c[0,k+1])
        m0, m1 = int(c[1,k]), int(c[1,k+1])
        dh, dm = h1-h0, m1-m0
        if dh and dm:
            n = min(dh, dm)
            for j in range(n):
                mapping[HUMAN_START+h0+j] = MOUSE_START+m0+j
        elif dh:
            for j in range(dh):
                mapping[HUMAN_START+h0+j] = None

    mapping.update(OFFICIAL_OVERRIDES)

    # The committed design mask contains additional structure-aware conserved mappings.
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


def rewrite(line, chain, resid, resname=None):
    chars = list(line.rstrip("\n"))
    if len(chars) < 80:
        chars.extend(" " * (80-len(chars)))
    chars[21] = chain
    chars[22:26] = list(f"{resid:4d}")
    if resname is not None:
        chars[17:20] = list(f"{resname:>3}")
    return "".join(chars).rstrip() + "\n"


def source_chain_residues(pdb, chain):
    vals = []
    seen = set()
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM") or line[21].strip() != chain:
            continue
        r = int(line[22:26])
        if r not in seen:
            vals.append(r)
            seen.add(r)
    return vals


def make_template(source_pdb, source_target, source_binder, binder_seq, species, hm_map, prefix):
    target_residues = source_chain_residues(source_pdb, source_target)
    binder_residues = source_chain_residues(source_pdb, source_binder)
    if not target_residues or not binder_residues:
        raise RuntimeError(f"Missing target/binder residues in {source_pdb}")
    if len(binder_residues) != len(binder_seq):
        raise RuntimeError(
            f"Binder/template length mismatch: {len(binder_residues)} vs {len(binder_seq)}"
        )

    target_first = min(target_residues)
    binder_index = {r:i+1 for i,r in enumerate(binder_residues)}
    target_lines, binder_lines = [], []

    for line in source_pdb.read_text(encoding="utf-8", errors="replace").splitlines(True):
        if not line.startswith("ATOM"):
            continue
        if line[76:78].strip() == "H" or line[12:16].strip().startswith("H"):
            continue
        ch = line[21].strip()
        old = int(line[22:26])

        if ch == source_target:
            human_uniprot = old - target_first + DOMAIN_START
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
            pos = binder_index[old]
            binder_lines.append(rewrite(line, "B", pos, AA3[binder_seq[pos-1]]))

    if not target_lines or not binder_lines:
        raise RuntimeError("Template construction produced an empty chain")

    pdb_path = prefix.with_suffix(".pdb")
    cif_path = prefix.with_suffix(".cif")
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(
        "".join(target_lines) + "TER\n" + "".join(binder_lines) + "TER\nEND\n",
        encoding="utf-8",
    )
    st = PDBParser(QUIET=True).get_structure(prefix.stem, str(pdb_path))
    io = MMCIFIO()
    io.set_structure(st)
    io.save(str(cif_path))
    return cif_path


def post_json(payload, api_key, max_attempts=5):
    for attempt in range(1, max_attempts+1):
        req = Request(
            URL,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type":"application/json","Authorization":f"Bearer {api_key}"},
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
            transient = (
                exc.code in {500,502,503,504}
                or (
                    exc.code == 422
                    and any(x in detail for x in ("InfrastructureError","AcceleratorError","CUDA error"))
                )
            )
            if transient and attempt < max_attempts:
                sleep_s = 25*attempt
                print(f"Boltz2 transient {exc.code}; sleep {sleep_s}s", flush=True)
                time.sleep(sleep_s)
                continue
            raise RuntimeError(f"Boltz2 HTTP {exc.code}: {detail}") from exc
    raise RuntimeError("Boltz2 retries exhausted")


def contacts(cif_path, offset, preferred, risk):
    st = MMCIFParser(QUIET=True).get_structure("pred", str(cif_path))
    model = next(st.get_models())
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
        "preferred_contacts": sorted(set(tc)&preferred),
        "risk_contacts": sorted(set(tc)&risk),
        "preferred_core_retained": bool(set(tc)&preferred),
    }


def select_one_per_backbone(rows):
    groups = {}
    for row in rows:
        groups.setdefault(row["backbone_id"], []).append(row)
    selected = []
    for _, rr in sorted(groups.items()):
        rr.sort(key=lambda x: float(x["mpnn_score"]) if x["mpnn_score"] else 999)
        selected.append(rr[0])
    return selected[:4]


def main():
    api_key = (
        os.getenv("NVIDIA_API_KEY")
        or os.getenv("NVIDIA_API_KEY1")
        or os.getenv("NVIDIA_API_TOKEN")
        or os.getenv("NGC_API_KEY")
    )
    if not api_key:
        raise RuntimeError("Missing NVIDIA API credential")
    if not CANDIDATES.exists():
        raise SystemExit(f"Candidate table not present: {CANDIDATES}")

    with CANDIDATES.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    rows = [r for r in rows if str(r.get("prefilter_pass","")).lower() in {"true","1","yes"}]
    # Single hotspot-engaging FreeBindCraft trajectory candidate; validate it directly.
    selected = rows[:4]
    if not selected:
        raise SystemExit("No prefiltered pose-screen candidate available")

    hm = human_mouse_map()
    mh = {m:h for h,m in hm.items() if m is not None}
    preferred_m = {hm[h] for h in PREFERRED_H if hm.get(h) is not None}
    risk_m = {hm[h] for h in RISK_H if hm.get(h) is not None}

    targets = {
        "human": {
            "sequence": fasta(HUMAN_FASTA),
            "offset": HUMAN_START-1,
            "preferred": PREFERRED_H,
            "risk": RISK_H,
        },
        "mouse": {
            "sequence": fasta(MOUSE_FASTA),
            "offset": MOUSE_START-1,
            "preferred": preferred_m,
            "risk": risk_m,
        },
    }

    OUT.mkdir(parents=True, exist_ok=True)
    batch = {}

    for row in selected:
        cid = row["candidate_id"]
        binder = row["sequence"]
        source = ROOT / row["source_backbone"]
        sb = row["binder_chain"]
        starget = row["target_chain"]
        cdir = OUT / cid
        cresult = {
            "candidate_id": cid,
            "backbone_id": row["backbone_id"],
            "source_backbone": row["source_backbone"],
            "binder_sequence": binder,
            "mpnn_score": row.get("mpnn_score") or None,
            "species": {},
        }

        for species, spec in targets.items():
            template = make_template(
                source, starget, sb, binder, species, hm,
                cdir / f"{cid}_{species}_apical_template",
            )
            ttext = template.read_text(encoding="utf-8")
            payload = {
                "polymers": [
                    {
                        "id":"A","molecule_type":"protein","sequence":spec["sequence"],
                        "structural_templates":[{
                            "structure":ttext,"format":"cif","name":f"{cid}_{species}_target","chain_id":"A"
                        }],
                    },
                    {
                        "id":"B","molecule_type":"protein","sequence":binder,
                        "structural_templates":[{
                            "structure":ttext,"format":"cif","name":f"{cid}_binder","chain_id":"B"
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
                raise RuntimeError(f"No structure for {cid}/{species}")
            spdir = cdir/species
            spdir.mkdir(parents=True, exist_ok=True)
            cif = spdir/"structure_1.cif"
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
                "confidence": result.get("confidence_scores", [None])[0],
                "ptm": result.get("ptm_scores", [None])[0],
                "protein_iptm": iptm,
                "complex_plddt": result.get("complex_plddt_scores", [None])[0],
                "complex_iplddt": result.get("complex_iplddt_scores", [None])[0],
                "complex_pde": result.get("complex_pde_scores", [None])[0],
                "complex_ipde": result.get("complex_ipde_scores", [None])[0],
                "contacts": cm,
                "apical_contact": apical,
                "pose_gate": {
                    "min_protein_iptm":0.50,
                    "requires_preferred_contact":True,
                    "pass":pose_pass,
                },
            }
            (spdir/"result_summary.json").write_text(
                json.dumps(summary, indent=2)+"\n", encoding="utf-8"
            )
            cresult["species"][species] = summary
            print(cid, species, json.dumps(summary, indent=2), flush=True)
            time.sleep(20)

        h = cresult["species"]["human"]
        m = cresult["species"]["mouse"]
        h_contacts = set(h["contacts"]["target_contacts_uniprot"])
        m_as_h = {mh[x] for x in m["contacts"]["target_contacts_uniprot"] if x in mh}
        overlap = sorted(h_contacts & m_as_h)
        union = h_contacts | m_as_h
        jaccard = len(overlap)/len(union) if union else 0.0
        same_epitope = len(overlap) >= 3 and jaccard >= 0.15
        cresult["cross_species_epitope"] = {
            "human_mouse_homologous_overlap_human_numbering": overlap,
            "overlap_count": len(overlap),
            "jaccard": round(jaccard, 4),
            "same_epitope_pass": same_epitope,
        }
        cresult["cross_species_pose_pass"] = bool(
            h["pose_gate"]["pass"]
            and m["pose_gate"]["pass"]
            and same_epitope
        )
        batch[cid] = cresult

    (OUT/"batch_summary.json").write_text(
        json.dumps(batch, indent=2)+"\n", encoding="utf-8"
    )
    winners = [k for k,v in batch.items() if v["cross_species_pose_pass"]]
    (OUT/"pose_pass_ids.json").write_text(
        json.dumps({"pose_pass_ids":winners}, indent=2)+"\n", encoding="utf-8"
    )
    print(json.dumps({"validated":len(batch),"cross_species_pose_pass":winners}, indent=2))


if __name__ == "__main__":
    main()
