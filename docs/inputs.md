# Input Data

Baysor consumes a molecule table plus an optional prior segmentation input.

## Supported Molecule Inputs

The `coordinates` positional argument may be:

- CSV molecule table
- Parquet molecule table
- `experiment.xenium`

When `experiment.xenium` is used, Baysor resolves the adjacent Xenium
transcript table automatically and keeps enough source context to make the
legacy outputs compatible with `xeniumranger import-segmentation`.

## Required Molecule Columns

At minimum, Baysor needs:

- `x`
- `y`
- `gene`

If the data are 3D, it can also use:

- `z`

Column names are configurable through CLI flags or the config file.

## Common Column Mapping Flags

- `--x-column`
- `--y-column`
- `--z-column`
- `--gene-column`
- `--qv-column`

For Xenium, [configs/xenium.toml](../configs/xenium.toml) already maps:

- `x_location`
- `y_location`
- `z_location`
- `feature_name`
- `qv`

## Optional Source Columns

Depending on the workflow, the input table may also carry:

- a transcript-native prior column such as `cell_id`
- confidence / quality columns
- Xenium `transcript_id`

For Xenium-origin inputs, Baysor preserves `transcript_id` so it can be written
back to `legacy` output when appropriate.

## Cropping During Input Load

The C++ branch supports spatial filtering during input loading:

- `--x-min`, `--x-max`
- `--y-min`, `--y-max`
- `--z-min`, `--z-max`

This is useful for:

- quick development runs
- protocol debugging
- visualization crops
- testing large datasets without loading the full field of view

## Protocol Notes

### Xenium

Preferred input:

```bash
./build/baysor run -c configs/xenium.toml data/experiment.xenium :cell_id
```

### ISS / osm-FISH / STARmap

These usually use CSV molecule tables plus either:

- no prior, with explicit `--scale`
- an image mask prior
- a boundary file prior

See [Examples](examples.md) for runnable commands.
