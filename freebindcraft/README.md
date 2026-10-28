# FreeBindCraft TfR1 fallback lane

This lane is the method-switch response to repeated RFdiffusion hotspot-control failures.

## Why switch generators?

RFdiffusion campaigns in this repository established three recurring failure modes:

1. full/apical targets often produce extended backbones or drift to unintended surfaces;
2. local epitope windows improve compactness but do not reliably enforce the requested hotspot patch;
3. sequences built on apparently acceptable RFdiffusion backbones frequently lose the designed pose when predicted against the full TfR1 ectodomain.

FreeBindCraft optimizes the complex through AlphaFold2 hallucination and then performs ProteinMPNN refinement plus AF2 revalidation, so binding geometry is part of the optimization loop rather than a downstream check only.

## Target

- Human TfR1 apical domain: residues 189–383 from PDB 11QC chain A.
- Hotspots: **343, 344, 347, 350**.
- Binder length smoke test: **55–75 aa**.
- Requested accepted designs: **1**.
- Cysteine excluded for the initial smoke run.

The target is deliberately the complete apical domain, not the 330–383 micro-window. BindCraft documentation recommends retaining enough surrounding structural context to avoid artifacts from over-trimming.

## Compute routes

### Preferred hosted route
Proto tool: `freebindcraft-design`.

The user's Proto account currently has no Modal workspace connected, so the tool is discoverable but not runnable until Modal is linked and the user explicitly approves deployment.

### Existing-resource fallback
Kaggle T4 may be viable for this small target+binder system. A smoke run should be attempted before any production campaign; the FreeBindCraft installation and AF2 weights are substantially heavier than the RFdiffusion Kaggle kernel.

## Non-negotiable downstream validation

A FreeBindCraft-accepted design is **not** automatically a competition candidate.

Every accepted sequence must pass:

1. full human TfR1 structural prediction;
2. full mouse TfR1 structural prediction;
3. protein ipTM >= 0.50 in both;
4. retained conserved apical contacts;
5. homologous human/mouse epitope overlap;
6. only then: pH-switch engineering and pH 6.0 negative-state evaluation.

This avoids repeating the earlier failure mode where confident monomer folding was mistaken for correct docking.
