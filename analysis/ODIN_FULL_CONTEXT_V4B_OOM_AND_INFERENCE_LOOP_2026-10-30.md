# Odin full-context v4b — monomer/PTM memory probe

Date: 2026-10-30.

## Question

V4 established that full human+mouse TfR1 ectodomain gradient optimisation with
the AF2-multimer design family could not fit on a Kaggle T4. The first design
step requested approximately 36.4 GB of device memory.

V4B tested whether changing only the AF2 design family to the lower-memory
monomer/PTM binder model would make the same full-context gradient objective
feasible.

Scientific target definitions were unchanged:

- complete human TfR1 ectodomain;
- complete mouse TfR1 ectodomain;
- human hotspots 347/350;
- mouse homologues 349/352;
- v6 60-aa warm-start;
- shared sequence across human and mouse contexts.

The probe budget was deliberately tiny:

- use_multimer_design = false;
- sample_models = false;
- 3 logits / 2 soft / 1 STE iterations;
- one trajectory.

## Result

The environment and GPU checks succeeded. JAX detected the T4:

`JAX devices: [CudaDevice(id=0)]`

The design nevertheless failed before producing a usable trajectory with:

`RESOURCE_EXHAUSTED: Out of memory while trying to allocate 35642827880 bytes`

Requested allocation:

**35,642,827,880 bytes (~35.6 GB decimal)**

The original multimer V4 request was ~36.4 GB. Switching to the monomer/PTM
design family therefore reduced the peak allocation only marginally and remains
far above the ~16 GB T4 device-memory budget.

No `selection.csv` or valid trajectory was produced.

## Decision

**Close full-ectodomain AF2 gradient backpropagation on Kaggle T4.**

Do not continue reducing iteration counts or opening V4C on this hardware. The
memory bottleneck occurs in a single full-context gradient step and is therefore
not solved by shortening the trajectory.

This is a hardware feasibility result, not a scientific rejection of full-context
multi-state design.

## Next architecture: full-context inference-in-the-loop

The project now uses:

1. tractable sequence/backbone proposal mechanisms;
2. full human and mouse ectodomain inference as the actual selection function;
3. predicted full-context poses fed back into the next sequence-design round;
4. strict homologous-apical and cross-species gates after every round.

Inference Loop R1 is seeded from two independently predicted full-context states:

- **human state:** FreeBindCraft v6 full-human prediction, protein iPTM 0.690,
  preferred contacts 231/347/350;
- **mouse state:** Odin V3 early-logit F39 full-mouse prediction, protein iPTM
  0.223, but physically contacting all four desired mouse apical residues
  345/346/349/352.

The two binder copies are tied in ProteinMPNN and only the union of their
full-context observed interface positions is redesignable. Proposed sequences
are then re-predicted from sequence only against the complete human and mouse
ectodomains.

This establishes the first actual **full-context inference feedback loop** in the
TfR1 campaign.
