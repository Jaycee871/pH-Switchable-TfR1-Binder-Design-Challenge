#!/usr/bin/env python3
"""Analyze an RFD4 TfR1 joint-design mmCIF without assuming output chain IDs.

The fixed target chain is identified by sequence similarity to human TfR1
residues 190-400. All target contacts are then mapped by target-chain order back
to UniProt numbering, so analysis does not depend on generated residue IDs.
"""

from __future__ import annotations

import argparse
import json
import math
from difflib import SequenceMatcher
from pathlib import Path

from Bio.PDB import MMCIFParser, NeighborSearch

AA3 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}
PREFERRED = {343,344,347,350}
CROSS_SPECIES_RISK = {204,208,210,211,212,348,371}


def parse_fasta(path: Path) -> str:
    return "".join(
        line.strip() for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    )


def protein_residues(chain):
    return [
        r for r in chain
        if r.id[0] == " " and r.resname in AA3
    ]


def sequence_of(chain) -> str:
    return "".join(AA3[r.resname] for r in protein_residues(chain))


def distance(a,b):
    return math.sqrt(sum((x-y)**2 for x,y in zip(a,b)))


def compactness(chain):
    cas=[r["CA"].coord for r in protein_residues(chain) if "CA" in r]
    if not cas:
        return {}
    centroid=[sum(float(p[i]) for p in cas)/len(cas) for i in range(3)]
    rg=math.sqrt(
        sum(sum((float(p[i])-centroid[i])**2 for i in range(3)) for p in cas)
        / len(cas)
    )
    max_span=0.0
    for i in range(len(cas)):
        for j in range(i+1,len(cas)):
            max_span=max(max_span,distance(cas[i],cas[j]))
    end=distance(cas[0],cas[-1]) if len(cas)>1 else 0.0
    return {
        "radius_of_gyration_A":round(rg,3),
        "max_CA_span_A":round(max_span,3),
        "end_to_end_A":round(end,3),
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--cif",required=True)
    ap.add_argument("--human-fasta",required=True)
    ap.add_argument("--output",required=True)
    args=ap.parse_args()

    cif=Path(args.cif)
    human_full=parse_fasta(Path(args.human_fasta))
    target_ref=human_full[190-121:400-121+1]  # FASTA itself is residues 121-760
    assert len(target_ref)==211

    structure=MMCIFParser(QUIET=True).get_structure("rfd4",str(cif))
    model=next(structure.get_models())
    chains=[c for c in model if protein_residues(c)]
    if len(chains)<2:
        raise RuntimeError(f"Expected >=2 protein chains, found {len(chains)}")

    chain_info=[]
    for chain in chains:
        seq=sequence_of(chain)
        similarity=SequenceMatcher(None,seq,target_ref).ratio()
        chain_info.append((chain,seq,similarity))

    chain_info.sort(key=lambda x:x[2],reverse=True)
    target,target_seq,target_similarity=chain_info[0]
    binder_candidates=[x for x in chain_info[1:] if 10 <= len(x[1]) <= 250]
    if not binder_candidates:
        raise RuntimeError("Could not identify a plausible designed binder chain")
    binder,binder_seq,_=binder_candidates[0]

    target_res=protein_residues(target)
    binder_res=protein_residues(binder)

    target_atoms=[
        atom for r in target_res for atom in r
        if getattr(atom,"element","") != "H"
    ]
    search=NeighborSearch(target_atoms)
    target_index={id(r):i for i,r in enumerate(target_res)}

    contact_target=set()
    contact_binder=set()
    residue_min_distance={u:999.0 for u in PREFERRED | CROSS_SPECIES_RISK}

    # General 5-A contact set.
    for bi,br in enumerate(binder_res, start=1):
        for atom in br:
            if getattr(atom,"element","")=="H":
                continue
            near=search.search(atom.coord,5.0,level="A")
            if near:
                contact_binder.add(bi)
            for ta in near:
                tr=ta.get_parent()
                idx=target_index[id(tr)]
                contact_target.add(190+idx)

    # Exact minimum heavy-atom distances for preferred/risk target residues.
    binder_atoms=[
        a for r in binder_res for a in r if getattr(a,"element","")!="H"
    ]
    for idx,tr in enumerate(target_res):
        uniprot=190+idx
        if uniprot not in residue_min_distance:
            continue
        best=999.0
        for ta in tr:
            if getattr(ta,"element","")=="H":
                continue
            for ba in binder_atoms:
                best=min(best,distance(ta.coord,ba.coord))
        residue_min_distance[uniprot]=round(best,3)

    preferred_contacts=sorted(set(contact_target)&PREFERRED)
    risk_contacts=sorted(set(contact_target)&CROSS_SPECIES_RISK)

    result={
        "cif":str(cif),
        "protein_chains":[
            {
                "chain_id":c.id,
                "length":len(seq),
                "target_similarity":round(sim,4),
            }
            for c,seq,sim in chain_info
        ],
        "target":{
            "chain_id":target.id,
            "length":len(target_seq),
            "reference_length":len(target_ref),
            "sequence_similarity_to_human_190_400":round(target_similarity,4),
        },
        "binder":{
            "chain_id":binder.id,
            "length":len(binder_seq),
            "sequence":binder_seq,
            "histidine_count":binder_seq.count("H"),
            "aromatic_count":sum(binder_seq.count(x) for x in "FWY"),
            "hydrophobic_count":sum(binder_seq.count(x) for x in "AILMVFWY"),
            "acidic_count":binder_seq.count("D")+binder_seq.count("E"),
            "basic_count":binder_seq.count("K")+binder_seq.count("R"),
            "compactness":compactness(binder),
        },
        "interface":{
            "target_contacts_uniprot":sorted(contact_target),
            "binder_contact_positions":sorted(contact_binder),
            "preferred_residues":sorted(PREFERRED),
            "preferred_contacts":preferred_contacts,
            "preferred_contact_count":len(preferred_contacts),
            "cross_species_risk_contacts":risk_contacts,
            "min_distance_A":{
                str(k):v for k,v in sorted(residue_min_distance.items())
            },
        },
    }

    out=Path(args.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
