# Native baseline fixture

This is the focused end-to-end reference for Native Slice N0. It captures the
committed pre-extraction CLI at revision
`f46a1e1dce1606d0ea644f4f8f1cf682597ba65c`.

The fixture uses a transcript-native prior, lets the CLI resolve
`n_cells_init`, disables optional molecule clustering and NCV colour generation,
and writes the legacy bundle so `segmentation.csv` retains `transcript_id`.

The reference was captured twice in fresh processes with the pre-extraction
CLI's implicit segmentation seed of `1`, `OMP_NUM_THREADS=1`, and
`OMP_DYNAMIC=FALSE`. The molecule table, cell statistics, count matrix, and
polygon output were byte-identical across both runs. The native regression test
compares their parsed scientific content and ignores log text and the
environment-specific CLI comment in the generated parameter dump.

The locked invocation is:

```text
baysor run molecules.csv --config config.toml --output <output-directory> \
  --output-style legacy --polygon-format FeatureCollection \
  --count-matrix-format tsv --skip-ncv-color
```
