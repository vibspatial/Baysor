# Example run on 10x Xenium pancreas data

This is an example of a Xenium v1 dataset with cell segmentation kit prior. The example is based on the 10x Genomics dataset [FFPE Human Pancreas with Xenium Multimodal Cell Segmentation](https://www.10xgenomics.com/datasets/ffpe-human-pancreas-with-xenium-multimodal-cell-segmentation-1-standard).

The panel has 377 genes. The dataset page reports 140,702 detected cells and 7,170,423 high-quality decoded transcripts. 

## Get the data

We'll use `data/` for the raw 10x bundle and `tests/` for Baysor outputs.

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

- `data/transcripts.parquet`
- `data/cell_boundaries.parquet`
- `data/nucleus_boundaries.parquet`
- `data/experiment.xenium`

## CLI run

These commands assume you run them from this directory. Use `../../build/baysor` if running the compiled binary from the repo build tree, or replace it with `baysor` if the executable is already on your `PATH`.

### With transcript-native prior labels

This uses the `cell_id` column already present in `transcripts.parquet`.

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  -o ./tests/full_cellid \
  ./data/transcripts.parquet \
  :cell_id
```

### With cell boundary priors

This assigns molecules to prior segments by point-in-polygon against the Xenium cell boundary file.

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  -o ./tests/full_cell_boundaries \
  ./data/transcripts.parquet \
  ./data/cell_boundaries.parquet
```

### With nucleus boundary priors

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  -o ./tests/full_nucleus_boundaries \
  ./data/transcripts.parquet \
  ./data/nucleus_boundaries.parquet
```

### Without priors

Without any prior segmentation input, Baysor needs an explicit scale.

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  --scale 4.5 \
  -o ./tests/full_no_prior \
  ./data/transcripts.parquet
```

## Cropped runs

baysor-cpp supports on-the-fly spatial filtering during input loading. This is useful for testing and for very large datasets. Here are examples of how to limit process to a particular region:

### Crop with transcript-native prior labels

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  --x-min 0 --x-max 2000 \
  --y-min 0 --y-max 2000 \
  -o ./tests/crop_cellid \
  ./data/transcripts.parquet \
  :cell_id
```

### Crop with cell boundary priors

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  --x-min 0 --x-max 2000 \
  --y-min 0 --y-max 2000 \
  -o ./tests/crop_cell_boundaries \
  ./data/transcripts.parquet \
  ./data/cell_boundaries.parquet
```

### Crop without priors

```bash
../../build/baysor run \
  -c ../../configs/xenium.toml \
  --x-min 0 --x-max 2000 \
  --y-min 0 --y-max 2000 \
  --scale 4.5 \
  -o ./tests/crop_no_prior \
  ./data/transcripts.parquet
```

## Notes

- `:cell_id` is usually the fastest Xenium prior mode, because the prior labels are already attached to each transcript row.
- `cell_boundaries.parquet` and `nucleus_boundaries.parquet` are direct geometric priors and can be useful when you want the segmentation to follow the published Xenium polygons more explicitly.
- `configs/xenium.toml` contains the Xenium column mapping and default transcript filters used in these commands.
