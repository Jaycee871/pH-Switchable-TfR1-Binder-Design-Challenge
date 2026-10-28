#!/usr/bin/env python3
"""Templated Boltz2 cross-species test for TfR1 Family-A candidate A29_03."""

from __future__ import annotations

import json
import os
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
OUTDIR = ROOT / "predictions/family_A_A29_03_templated"

CANDIDATE_ID = "TFR1_A29_MPNN_03"
BINDER = (
    "ALDEEVEKLLKEVQALEKEVAALKKKKKELEEKAKKEKNEEAEKIKEEAKKLGKEKEEKEEKL"
    "KELEKKLSEKAKKYGDKLEEEKKK"
)

TARGETS = {
    "human": ROOT / "data/human_TfR1/human_TfR1_P02786_121-760.fasta",
    "mouse": ROOT / "data/mouse_TfR1/mouse_TfR1_Q62351_123-763.fasta",
}

DEFAULT_TEMPLATE = ROOT / "predictions/templates/A29_03_rf_seed20261029_template.pdb"


def parse_fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    )


def post_json(payload: dict, api_key: str) -> dict:
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
        raise RuntimeError(f"Boltz2 HTTP {exc.code}: {detail}") from exc


def main() -> None:
    api_key = os.getenv("NVIDIA_API_KEY") or os.getenv("NGC_API_KEY")
    if not api_key:
        raise RuntimeError("Missing NVIDIA_API_KEY / NGC_API_KEY")

    template_path = Path(os.getenv("TEMPLATE_PDB", str(DEFAULT_TEMPLATE)))
    if not template_path.exists():
        raise RuntimeError(f"Template not found: {template_path}")
    template = template_path.read_text(encoding="utf-8")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    overall = {
        "candidate_id": CANDIDATE_ID,
        "binder_sequence": BINDER,
        "binder_length": len(BINDER),
        "template": str(template_path),
        "species": {},
    }

    for species, fasta_path in TARGETS.items():
        target = parse_fasta(fasta_path)
        payload = {
            "polymers": [
                {
                    "id": "A",
                    "molecule_type": "protein",
                    "sequence": target,
                    "structural_templates": [
                        {
                            "structure": template,
                            "format": "pdb",
                            "name": "rfseed20261029_target",
                            "chain_id": "A",
                        }
                    ],
                },
                {
                    "id": "B",
                    "molecule_type": "protein",
                    "sequence": BINDER,
                    "structural_templates": [
                        {
                            "structure": template,
                            "format": "pdb",
                            "name": "rfseed20261029_binder",
                            "chain_id": "B",
                        }
                    ],
                },
            ],
            "recycling_steps": 3,
            "sampling_steps": 50,
            "diffusion_samples": 1,
            "step_scale": 1.638,
            "output_format": "mmcif",
            "write_full_pae": False,
        }

        result = post_json(payload, api_key)
        structures = result.get("structures", [])
        if not structures:
            raise RuntimeError(
                f"No templated Boltz2 structure returned for {species}; "
                f"keys={sorted(result.keys())}"
            )

        spdir = OUTDIR / species
        spdir.mkdir(parents=True, exist_ok=True)
        for i, item in enumerate(structures, start=1):
            (spdir / f"structure_{i}.cif").write_text(
                item["structure"], encoding="utf-8"
            )

        summary = {
            "target_length": len(target),
            "binder_length": len(BINDER),
            "confidence_scores": result.get("confidence_scores", []),
            "structure_count": len(structures),
            "metrics": result.get("metrics", {}),
        }
        (spdir / "result_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        overall["species"][species] = summary
        print(species, json.dumps(summary, indent=2))

    (OUTDIR / "cross_species_summary.json").write_text(
        json.dumps(overall, indent=2) + "\n", encoding="utf-8"
    )
    print("Templated Boltz2 cross-species test complete")


if __name__ == "__main__":
    main()
