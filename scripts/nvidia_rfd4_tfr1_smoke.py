#!/usr/bin/env python3
"""One-replicate RFD4-Proteina joint sequence/structure smoke test for TfR1.

The target is held by internal backbone-coordinate conditioning (C_CRDGRP) and
full target sequence conditioning (C_SEQ). The conserved apical subcluster is
marked with atom-level C_HOT selections.

This script performs one POST only. It does not retry ambiguous 202/429/5xx
responses, preventing accidental duplicate generation.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "structures/rfd4/human_TfR1_11QC_A_190-400.cif"
OUTDIR = ROOT / "designs/rfd4_joint/smoke_001"
URL = "https://health.api.nvidia.com/v1/biology/ipd/rfd4-proteina/generative"

EXAMPLE_ID = "tfr1_apical_joint_001"
BINDER_LENGTH = "80-100"

# Atom-level hotspot definitions. G347 has no side chain, so CA marks its
# conserved backbone location while E/K hotspots use terminal side-chain atoms.
HOTSPOT_SELECTIONS = [
    "sel('A/*/343/[CD,OE1,OE2]')",
    "sel('A/*/344/[CE,NZ]')",
    "sel('A/*/347/[CA]')",
    "sel('A/*/350/[CD,OE1,OE2]')",
]


def build_payload(target_text: str) -> dict:
    # Hosted NIM condition_spec expects the design specification itself here;
    # unlike a local multi-example YAML file, do NOT wrap it in an example-name key.
    spec = {
        "input": TARGET.name,
        "contig": f"{BINDER_LENGTH}, /0, A190-400",
        "conditions": {
            "C_CRDGRP": {
                "1": [
                    {
                        "select": (
                            "sel('A/*/*/[N,CA,C,O]') "
                            "and not has_nan_coord()"
                        )
                    }
                ]
            },
            "C_SEQ": {
                "true": [{"select": "sel('A')"}]
            },
            "C_HOT": {
                "true": [
                    {"select": selection}
                    for selection in HOTSPOT_SELECTIONS
                ]
            },
        },
    }

    return {
        "task": "generation",
        "input": {
            "format": "condition_spec",
            "example_id": EXAMPLE_ID,
            "spec": spec,
            "files": [
                {
                    "filename": TARGET.name,
                    "content": target_text,
                }
            ],
        },
        "sampling": {
            "recycle_strategy": "thinking_medium",
            "num_replicates": 1,
            "batch_size": 1,
            "seed": 20261103,
        },
        "metrics": {
            "tags": [
                "design",
                "hotspot_conditioning",
                "interchain_interactions",
            ]
        },
        "output": {
            "structure_format": "mmcif",
            "include_sequence": True,
        },
    }


def main() -> None:
    api_key = os.getenv("NVIDIA_API_KEY") or os.getenv("NGC_API_KEY")
    if not api_key:
        raise RuntimeError("Missing NVIDIA_API_KEY / NGC_API_KEY")
    if not TARGET.exists():
        raise RuntimeError(
            f"Missing RFD4 target crop: {TARGET}. Run rfd4_prepare_tfr1_target.py"
        )

    OUTDIR.mkdir(parents=True, exist_ok=True)
    payload = build_payload(TARGET.read_text(encoding="utf-8"))

    # Persist the exact non-secret request for reproducibility.
    (OUTDIR / "request.json").write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )

    req = Request(
        URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {api_key}",
            "NVCF-POLL-SECONDS": "300",
        },
        method="POST",
    )

    try:
        with urlopen(req, timeout=900) as resp:
            status = resp.status
            body = resp.read().decode("utf-8")
            headers = dict(resp.headers.items())
    except HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        (OUTDIR / "transport_error.txt").write_text(
            body + "\n", encoding="utf-8"
        )
        details = {
            "http_status": exc.code,
            "retry_after": exc.headers.get("Retry-After"),
            "nvcf_reqid": exc.headers.get("NVCF-REQID"),
            "x_request_id": exc.headers.get("X-Request-ID"),
        }
        (OUTDIR / "transport_error.json").write_text(
            json.dumps(details, indent=2) + "\n", encoding="utf-8"
        )
        raise RuntimeError(
            f"RFD4-Proteina HTTP {exc.code}; diagnostics saved"
        ) from exc

    (OUTDIR / "response.raw").write_text(body, encoding="utf-8")
    transport = {
        "http_status": status,
        "nvcf_reqid": headers.get("NVCF-REQID"),
        "x_request_id": headers.get("X-Request-ID"),
        "retry_after": headers.get("Retry-After"),
    }
    (OUTDIR / "transport.json").write_text(
        json.dumps(transport, indent=2) + "\n", encoding="utf-8"
    )

    if status == 202:
        raise SystemExit(
            "RFD4 request is pending (HTTP 202); request ID was preserved. "
            "Do not re-POST automatically."
        )
    if status != 200:
        raise RuntimeError(f"Unexpected RFD4 HTTP status {status}")

    response = json.loads(body)
    (OUTDIR / "response.json").write_text(
        json.dumps(response, indent=2) + "\n", encoding="utf-8"
    )

    designs = response.get("designs", [])
    if not designs:
        raise RuntimeError("RFD4 returned HTTP 200 but no designs")

    summaries = []
    fasta_rows = []
    for i, design in enumerate(designs):
        diag = design.get("diagnostics") or {}
        if diag.get("all_coordinates_finite") is False:
            raise RuntimeError(f"Design {i} has non-finite coordinates")

        structure = design.get("structure") or {}
        content = structure.get("content")
        if content:
            (OUTDIR / f"design_{i:04d}.cif").write_text(
                content, encoding="utf-8"
            )

        sequence = design.get("sequence")
        if sequence:
            fasta_rows.append(f">RFD4_TFR1_{i:04d}\n{sequence}\n")

        summaries.append(
            {
                "index": i,
                "id": design.get("id"),
                "replicate_index": design.get("replicate_index"),
                "origin": design.get("origin"),
                "seed": design.get("seed"),
                "sequence": sequence,
                "sequence_length": len(sequence) if sequence else None,
                "chain_of_interest": design.get("chain_of_interest"),
                "metrics": design.get("metrics"),
                "diagnostics": diag,
            }
        )

    if fasta_rows:
        (OUTDIR / "sequences.fasta").write_text(
            "".join(fasta_rows), encoding="utf-8"
        )

    summary = {
        "status": "completed",
        "request_id": response.get("request_id"),
        "task": response.get("task"),
        "duration_ms": response.get("duration_ms"),
        "model": response.get("model"),
        "runtime": response.get("runtime"),
        "schedule": response.get("schedule"),
        "warnings": response.get("warnings"),
        "aggregated_metrics": response.get("aggregated_metrics"),
        "design_count": len(designs),
        "designs": summaries,
    }
    (OUTDIR / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print("RFD4-Proteina TfR1 joint-design smoke test completed")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
