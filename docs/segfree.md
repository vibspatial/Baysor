# Segmentation-Free NCVs

`segfree` extracts neighborhood-composition vectors (NCVs) without running the
cell segmentation algorithm.

CLI shape:

```bash
./build/baysor segfree [OPTIONS] coordinates
```

## Typical Use

```bash
./build/baysor segfree \
  -c configs/xenium.toml \
  -k 100 \
  -o ncvs.loom \
  data/transcripts.parquet
```

## What It Computes

The current pipeline:

- loads molecules
- builds neighborhood-composition counts
- log-transforms the neighborhood matrix
- estimates molecule confidences
- computes NCV colors
- writes a Loom file with per-neighborhood vectors and attributes

## Common Options

- `-c,--config`
- `-x,--x-column`
- `-y,--y-column`
- `-z,--z-column`
- `-g,--gene-column`
- `--qv-column`
- `--min-qv`
- `--x-min`, `--x-max`
- `--y-min`, `--y-max`
- `--z-min`, `--z-max`
- `-k,--k-neighbors`
- `-o,--output`
- `--force-2d`

## Output

The default output is:

```text
ncvs.loom
```

The Loom file includes:

- the NCV matrix
- `ncv_color`
- per-molecule confidence

## Notes

- `segfree` accepts a transcript table directly.
- for Xenium datasets, it also accepts `experiment.xenium` and resolves the
  underlying transcript table automatically.
