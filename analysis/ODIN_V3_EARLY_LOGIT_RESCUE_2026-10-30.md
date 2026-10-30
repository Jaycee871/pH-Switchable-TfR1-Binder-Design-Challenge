# Odin anti-escape v3 early-logit rescue

The official `best_clipped_i_pae_ratio` selector chose STE frame 73 because
the earlier logits frames are not yet one-hot discrete sequences.

However, the complete trajectory artifact shows substantially better joint
target/off-target geometry earlier in the logits stage.

The four best discrete argmax rescues are:

| ID | frame | worst target iPTM | clipped specificity ratio |
|---|---:|---:|---:|
| ODIN_V3_LOGIT_F39 | 39 | 0.7636 | 1.8995 |
| ODIN_V3_LOGIT_F37 | 37 | 0.7612 | 1.8061 |
| ODIN_V3_LOGIT_F31 | 31 | 0.7150 | 1.5565 |
| ODIN_V3_LOGIT_F36 | 36 | 0.7131 | 1.6249 |

Each sequence is produced by argmax over the saved 20-AA logits using
ColabDesign/AlphaFold's internal residue order:

`ARNDCQEGHILKMFPSTWYV`

These sequences are hypotheses only. They are sent directly to independent
sequence-only full human/full mouse Boltz-2 validation. No Odin pose is reused.
