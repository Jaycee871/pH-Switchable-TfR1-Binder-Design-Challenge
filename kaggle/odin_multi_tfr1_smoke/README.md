# Odin-Multi TfR1 cross-reactivity smoke

True multi-context design after the FreeBindCraft v6 backbone and all subsequent
single-state / fixed-backbone sequence rescues failed to preserve the homologous
mouse apical pose.

## Contexts

- Human TfR1 apical domain, residues 189-383
  - hotspots: 343, 344, 347, 350
- Mouse TfR1 homologous apical domain, residues 191-385
  - hotspots: 345, 346, 349, 352

Both contexts use Odin-Multi's **target** role, so their gradients jointly update
one shared binder sequence.

## Smoke budget

- binder: fixed 60 aa
- 1 trajectory
- 30 soft / 20 temporary / 10 hard / 5 greedy iterations
- contact cutoff: 8 A
- inter-contact weight: 1.25
- one required contact per hotspot
- selection: `best_i_ptm` (maximize the weakest on-target iPTM)

This run tests whether true multi-context optimisation can solve the cross-species
pose failure that single-target FreeBindCraft and post-hoc ProteinMPNN could not.

Odin-Multi selection is not the final gate. Selected sequences must still pass
the repository's independent full human/mouse Boltz-2 validation before any
pH-switch engineering.

## Input packaging (fixed 2026-10-29)

`kernel_type` is `script`, so `kaggle kernels push -p <dir>` uploads **only** the file
named by `code_file`. The first smoke run (2026-10-29 20:04 Taipei) therefore died with

```
FileNotFoundError: [Errno 2] No such file or directory: '/kaggle/src/human_tfr1_apical_189_383.pdb'
```

after ~12 minutes of conda installation, because the two target PDBs in this directory
were never uploaded.

The two PDBs are now **embedded inside `tfr1_odin_multi_smoke.py`** as gzip+base64
payloads, regenerated with:

```
python scripts/embed_odin_multi_targets.py
```

The `.pdb` files in this directory remain the canonical source of truth; the embedded
copies are verified against them by SHA-256 inside the kernel.

Alternatives considered and rejected:

- **Raw GitHub URLs** — this repository is private, so the kernel would need a GitHub
  token shipped into a third-party notebook. Not acceptable.
- **RCSB download at kernel start** (the pattern used by
  `kaggle/freebindcraft_v6_tempered_binary` and `kaggle/rfdiffusion_pilot_v1`) — works for
  the human crop of 11QC chain A 189-383, but `mouse_tfr1_apical_homolog.pdb` is a curated
  homology model that exists only in this repository and cannot be regenerated from RCSB.
- **A Kaggle dataset in `dataset_sources`** — viable, but adds an external asset plus
  create-vs-version branching, and the Kaggle credential exists only as a GitHub Actions
  secret, so the dataset could not be prepared or verified locally.

`main()` now runs `materialize_targets()` as **Stage 0, before the environment install**.
It decodes both PDBs, checks SHA-256, checks CA residue ranges (human 189-383 n=195,
mouse 191-385 n=194) and asserts all four human/mouse hotspot pairs
(343/345, 344/346, 347/349, 350/352) are present. Any packaging regression now fails in
seconds instead of after a 12-minute install.
