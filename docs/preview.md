# Dataset Preview

`preview` generates an HTML overview of a dataset without running full
segmentation.

CLI shape:

```bash
./build/baysor preview [OPTIONS] coordinates
```

## Typical Use

```bash
./build/baysor preview \
  -c configs/xenium.toml \
  -o preview.html \
  data/transcripts.parquet
```

or for Xenium-style columns:

```bash
./build/baysor preview \
  -c configs/xenium.toml \
  data/transcripts.parquet
```

## What It Computes

The preview pipeline currently does the following:

- loads molecules
- estimates a molecule-confidence / noise model
- computes neighborhood-composition colors
- estimates a gene-structure embedding
- writes an HTML summary report

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
- `-o,--output`
- `--force-2d`

## Output

The default output is:

```text
preview.html
```

## Notes

- `preview` accepts a transcript table directly.
- for Xenium datasets, it also accepts `experiment.xenium` and resolves the
  underlying transcript table automatically.
- this is useful for quick sanity checks before a full segmentation run.
