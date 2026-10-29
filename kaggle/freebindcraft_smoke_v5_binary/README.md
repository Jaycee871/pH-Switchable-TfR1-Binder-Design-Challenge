# FreeBindCraft v5 — binary hard-hotspot fallback

Prepared only as a fallback if v4 fails the independent 5 A hotspot gate.

Same as v4:
- TfR1 11QC chain A residues 189–383
- hotspot residues 343/344/347/350
- fixed 60-aa binder
- one attempt
- inter-contact cutoff 8 A
- inter-contact loss weight 2.0
- two binder contacts per hotspot position
- independent final 5 A hotspot geometry gate

Single method change:
- FreeBindCraft/ColabDesign `i_con` is patched from `binary=False` to
  `binary=True`.

In ColabDesign this changes the contact objective from categorical distogram
cross-entropy to `-log P(distance < cutoff)`. It is a deliberate source-level
experiment, not a standard FreeBindCraft preset.

Do not launch this campaign if v4 already achieves a true 5 A hotspot gate pass.
