# Odin-Multi full-context v4b memory probe

V4 failed before producing a trajectory because AF2-multimer full-context
backpropagation attempted a ~36.4 GB allocation on a Kaggle T4.

This probe changes **only the AF2 design family and trajectory budget**:

- complete human TfR1 ectodomain retained;
- complete mouse TfR1 ectodomain retained;
- same 347/350 ↔ 349/352 hotspots;
- same v6 warm-start;
- same human/mouse gradient weights;
- `use_multimer_design: false` to use the lower-memory monomer/PTM binder family;
- `sample_models: false`;
- 3 logits / 2 soft / 1 STE iterations;
- one trajectory.

This is a hardware-feasibility probe, not a production candidate campaign.

If it still OOMs, full-ectodomain AF2 gradient design is considered infeasible
on the available T4 and the project should switch to full-context
**inference-in-the-loop** rather than continue reducing the same gradient model.

If it succeeds, the measured runtime and memory feasibility justify a bounded
production v4c using the same full structural contexts.
