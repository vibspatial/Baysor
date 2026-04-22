# Example Run On 10x Xenium Pancreas Data

This example uses a Xenium v1 dataset with the cell-segmentation-kit prior. The
dataset is based on the 10x Genomics example
[FFPE Human Pancreas with Xenium Multimodal Cell Segmentation](https://www.10xgenomics.com/datasets/ffpe-human-pancreas-with-xenium-multimodal-cell-segmentation-1-standard).

The panel has 377 genes. The dataset page reports 140,702 detected cells and
7,170,423 high-quality decoded transcripts.

## Get The Data

Use `data/` for the raw Xenium bundle and `tests/` for Baysor outputs.

```bash
mkdir -p data tests
```

Download the Xenium output bundle from 10x:

```bash
curl -L https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Pancreas_FFPE/Xenium_V1_human_Pancreas_FFPE_outs.zip \
  -o data/Xenium_V1_human_Pancreas_FFPE_outs.zip
unzip data/Xenium_V1_human_Pancreas_FFPE_outs.zip -d data
```

The main files used below are:

- `data/experiment.xenium`
- `data/transcripts.parquet`
- `data/cell_boundaries.parquet`
- `data/nucleus_boundaries.parquet`

## Recommended Input Form

For Xenium runs, prefer passing `experiment.xenium` as the main input:

```bash
../../build/baysor run ... ./data/experiment.xenium ...
```

That gives Baysor enough source context to:

- resolve the bundled transcript table automatically
- preserve `transcript_id` in `legacy` output
- make the legacy outputs directly usable by
  `xeniumranger import-segmentation`

Passing `transcripts.parquet` directly still works, but `experiment.xenium` is
the cleaner Xenium-aware path.

## CLI Runs

These commands assume you run them from this directory. Use `../../build/baysor`
if running the binary from the repo build tree, or replace it with `baysor` if
the executable is already on your `PATH`.

### Full Run With Transcript-Native Prior Labels

This uses the `cell_id` column already present in the Xenium transcript table.

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  -o ./tests/full_cellid \
  ./data/experiment.xenium \
  :cell_id
```

This is the recommended full Xenium run if you intend to hand the result off to
Xenium Ranger / Xenium Explorer later.

### Alternative Clustering Priors

By default, Baysor uses the legacy `mrf` clustering prior with `4` clusters.
The recommended graph-based alternative is:

- `louvain`

Graph-based clustering uses NCV basis anchors and then transfers the final
labels back to all molecules. For `louvain`, `--n-clusters` is the target final
number of coarse clusters after anchor communities are merged.

Examples:

```bash
../../build/baysor run \
  -p \
  -c ../../configs/xenium.toml \
  -o ./tests/full_cellid_mrf \
  ./data/experiment.xenium \
  :cell_id
```

```bash
../../build/baysor run \
  -p \
  --cluster-method louvain \
  --n-clusters 10 \
  -c ../../configs/xenium.toml \
  -o ./tests/full_cellid_louvain_k10 \
  ./data/experiment.xenium \
  :cell_id
```

`-p` adds the HTML run report, which is useful for checking the NCV manifold,
segmentation summary plots, and output bundle before exporting to Xenium Ranger.

### Full Run With Cell Boundary Priors

This assigns molecules to prior segments by point-in-polygon against the Xenium
cell boundary file.

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  -o ./tests/full_cell_boundaries \
  ./data/experiment.xenium \
  ./data/cell_boundaries.parquet
```

### Full Run With Nucleus Boundary Priors

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  -o ./tests/full_nucleus_boundaries \
  ./data/experiment.xenium \
  ./data/nucleus_boundaries.parquet
```

### Full Run Without Priors

Without any prior segmentation input, Baysor needs an explicit scale.

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  --scale 4.5 \
  -o ./tests/full_no_prior \
  ./data/experiment.xenium
```

## Output Modes

The C++ branch currently supports two output styles:

- `legacy`
- `parquet`

`legacy` is the default. For Xenium-origin inputs, it automatically adds the
extra fields needed by `xeniumranger import-segmentation`:

- `segmentation.csv`
  - includes `transcript_id`
  - writes `is_noise` as `true` / `false`
- `segmentation_polygons_2d.json`
  - stays a GeoJSON `FeatureCollection`
  - includes `properties.cell` on each feature

### Legacy Output

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  --output-style legacy \
  -o ./tests/full_cellid_legacy \
  ./data/experiment.xenium \
  :cell_id
```

### Parquet Output

This emits:

- `molecules.parquet`
- `cells.parquet`
- `cell_boundaries.parquet`
- `cell_boundaries_3d.parquet`
- `feature_matrix.h5`

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  --output-style parquet \
  -o ./tests/full_cellid_parquet \
  ./data/experiment.xenium \
  :cell_id
```

## Xenium Explorer Handoff

The current recommended path to Xenium Explorer is:

1. run Baysor in `legacy` mode on the full Xenium bundle
2. pass Baysor outputs into `xeniumranger import-segmentation`

Use these two Baysor outputs for the handoff:

- `./tests/full_cellid/segmentation.csv`
- `./tests/full_cellid/segmentation_polygons_2d.json`

Example conversion command:

```bash
xeniumranger import-segmentation \
  --id baysor_pancreas \
  --xenium-bundle ./data \
  --transcript-assignment ./tests/full_cellid/segmentation.csv \
  --viz-polygons ./tests/full_cellid/segmentation_polygons_2d.json \
  --units microns
```

Additional Xenium Ranger flags can be added as needed, for example:

- `--localcores 16`
- `--localmem 64`

Notes:

- use a full run, not a crop, for Xenium Ranger import
- the recommended Baysor input for this workflow is `experiment.xenium`
- the legacy outputs now contain the extra Ranger-compatible fields automatically

## Cropped Runs

The C++ branch supports on-the-fly spatial filtering during input loading. This
is useful for testing and for very large datasets.

### Crop With Transcript-Native Prior Labels

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  --x-min 0 --x-max 2000 \
  --y-min 0 --y-max 2000 \
  -o ./tests/crop_cellid \
  ./data/experiment.xenium \
  :cell_id
```

### Crop With Cell Boundary Priors

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  --x-min 0 --x-max 2000 \
  --y-min 0 --y-max 2000 \
  -o ./tests/crop_cell_boundaries \
  ./data/experiment.xenium \
  ./data/cell_boundaries.parquet
```

### Crop Without Priors

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  --x-min 0 --x-max 2000 \
  --y-min 0 --y-max 2000 \
  --scale 4.5 \
  -o ./tests/crop_no_prior \
  ./data/experiment.xenium
```

Do not use cropped runs as the input to `xeniumranger import-segmentation`.

## Visual QC Helper

The repo also includes a Xenium visualization helper for rendering fixed example
regions from a Baysor run:

```bash
.venv-vis/bin/python ../../scripts/visualize_xenium_examples.py \
  -n 10 \
  ./tests/full_cellid/segmentation.csv
```

This renders side-by-side Baysor / Xenium comparison panels using the Xenium
morphology images as background when the required Python packages are available.

## Notes

- `:cell_id` is usually the fastest Xenium prior mode, because the prior labels
  are already attached to each transcript row.
- `cell_boundaries.parquet` and `nucleus_boundaries.parquet` are direct
  geometric priors and can be useful when you want the segmentation to follow
  the published Xenium polygons more explicitly.
- `configs/xenium.toml` contains the Xenium column mapping and default
  transcript filters used in these commands.
