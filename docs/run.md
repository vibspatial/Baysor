# Running Baysor

The main CLI entrypoint in the C++ branch is:

```bash
./build/baysor [SUBCOMMAND] ...
```

For segmentation, use:

```bash
./build/baysor run [OPTIONS] coordinates [prior_segmentation]
```

## Minimal Forms

With an explicit scale and no prior:

```bash
./build/baysor run --scale 8 molecules.csv
```

With a prior segmentation:

```bash
./build/baysor run molecules.csv segmentation_mask.tif
```

With a config file:

```bash
./build/baysor run -c configs/xenium.toml data/transcripts.parquet :cell_id
```

## Positionals

`coordinates`
- required
- molecule table in CSV or Parquet form
- may also be `experiment.xenium` for Xenium-aware runs

`prior_segmentation`
- optional
- one of:
  - image mask
  - boundary CSV / Parquet
  - `:column_name` for transcript-native prior labels in the molecule table

## Common Options

Geometry / column mapping:

- `-x,--x-column`
- `-y,--y-column`
- `-z,--z-column`
- `-g,--gene-column`
- `--qv-column`
- `--force-2d`

Segmentation scale / convergence:

- `-s,--scale`
- `--scale-std`
- `--iters`
- `--tol`
- `--n-cells-init`
- `--cluster-method`
- `--n-clusters`
- `--cluster-resolution`
- `--cluster-graph-k`
- `--cluster-n-dims`
- `--cluster-basis-sample-size`

Filtering:

- `--min-molecules-per-cell`
- `--min-molecules-per-gene`
- `--exclude-genes`
- `--min-qv`

Cropping:

- `--x-min`, `--x-max`
- `--y-min`, `--y-max`
- `--z-min`, `--z-max`

Outputs:

- `-o,--output`
- `--output-style legacy|parquet`
- `--polygon-format`
- `--count-matrix-format`
- `-p,--plot`
- `--skip-ncv-color`

For the files written by each bundle, see [Outputs](outputs.md) and
[Output Files](output_files.md).

## Clustering Methods

`run` supports four clustering modes:

- `mrf`
  - the legacy Baysor molecule clustering path
  - default mode
  - `--n-clusters` is the exact number of molecule clusters to fit
  - default `n_clusters`: `4`
- `louvain`
  - graph clustering on NCV basis anchors, followed by label transfer to all
    molecules
  - `--n-clusters` is the target final coarse cluster count after anchor
    communities are merged
  - default `n_clusters`: `10`
- `leiden`
  - same anchor-based NCV workflow as `louvain`, with Leiden refinement
  - `--n-clusters` is again the target final coarse cluster count
  - default `n_clusters`: `10`
- `none`
  - disables the clustering prior entirely

Examples:

```bash
./build/baysor run --cluster-method mrf --n-clusters 4 ...
./build/baysor run --cluster-method louvain --n-clusters 10 ...
./build/baysor run --cluster-method leiden --n-clusters 10 ...
```

`--cluster-resolution` and `--cluster-graph-k` are advanced controls for the
graph methods. In most datasets, the defaults are a reasonable starting point.

## Config Files

Most common protocol setups are captured in:

- [configs/xenium.toml](../configs/xenium.toml)
- [configs/iss.toml](../configs/iss.toml)
- [configs/starmap.toml](../configs/starmap.toml)
- [configs/osm_fish.toml](../configs/osm_fish.toml)

The example config is:

- [configs/example_config.toml](../configs/example_config.toml)

CLI flags override the config file.

## Other Subcommands

The current C++ CLI also provides:

- [Dataset Preview](preview.md)
- [Segmentation-Free NCVs](segfree.md)

## Full Help

For the up-to-date option list, use:

```bash
./build/baysor run --help
```
