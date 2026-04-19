# Xenium Workflow

This page describes the recommended workflow for 10x Xenium data in the C++
branch.

## Recommended Input

Use the Xenium manifest as the main input:

```bash
./build/baysor run -c configs/xenium.toml data/experiment.xenium :cell_id
```

Why:

- Baysor resolves the transcript table automatically
- Xenium-specific column mapping is already captured in `configs/xenium.toml`
- Xenium `transcript_id` can be preserved in `legacy` output

Passing `transcripts.parquet` directly still works, but `experiment.xenium` is
the preferred Xenium-aware entrypoint.

## Common Prior Modes

### Transcript-Native Prior

```bash
./build/baysor run -c configs/xenium.toml data/experiment.xenium :cell_id
```

This is usually the best default for Xenium.

### Cell Boundary Prior

```bash
./build/baysor run -c configs/xenium.toml \
  data/experiment.xenium data/cell_boundaries.parquet
```

### Nucleus Boundary Prior

```bash
./build/baysor run -c configs/xenium.toml \
  data/experiment.xenium data/nucleus_boundaries.parquet
```

## Recommended Output Style

For Xenium runs that may be handed off to Xenium Explorer, use `legacy`
output:

```bash
./build/baysor run -c configs/xenium.toml \
  --output-style legacy \
  -o out_dir \
  data/experiment.xenium :cell_id
```

For Xenium-origin inputs, this automatically produces Ranger-friendly:

- `segmentation.csv`
- `segmentation_polygons_2d.json`

## Xenium Explorer Handoff

The recommended Explorer path is:

1. run Baysor on the original Xenium bundle in `legacy` mode
2. run `xeniumranger import-segmentation`

Example:

```bash
xeniumranger import-segmentation \
  --id baysor_xenium \
  --xenium-bundle data \
  --transcript-assignment out_dir/segmentation.csv \
  --viz-polygons out_dir/segmentation_polygons_2d.json \
  --units microns
```

This is the preferred route instead of direct Baysor-side Xenium bundle
generation.

## Full Runs vs Crops

Use full runs for Xenium Ranger handoff.

Cropped runs are still useful for:

- development
- debugging
- visualization
- performance profiling

but they are not the right input to `xeniumranger import-segmentation`.

## Example

For a full runnable Xenium example, see:

- [examples/Xenium_pancreas_membrane_377/README.md](../examples/Xenium_pancreas_membrane_377/README.md)
